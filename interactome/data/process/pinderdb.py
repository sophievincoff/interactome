"""Code for processing the PINDER database, extracting both positives and negatives."""

from pinder.core import get_metadata
import numpy as np
import pandas as pd
import os
from protparser.rcsb import RCSBStructure
import time
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path
import logging
from functools import lru_cache
from tqdm.contrib.concurrent import thread_map
import time, threading, requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from typing import List, Optional, TextIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import math

import requests_cache
requests_cache.install_cache("rcsb_cache", expire_after=24*3600)  # 24h cache

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

# Simple global rate limiter (e.g., 2 req/sec)
MIN_INTERVAL = float(os.getenv("RCSB_MIN_INTERVAL", "0.5"))   # seconds between GETs
MAX_WORKERS  = int(os.getenv("RCSB_MAX_WORKERS", "8"))
def _rate_wait():
    _last = [0.0]
    _lock = threading.Lock()
    with _lock:
        now = time.perf_counter()
        wait = max(0.0, MIN_INTERVAL - (now - _last[0]))
        if wait > 0: time.sleep(wait)
        _last[0] = time.perf_counter()

@lru_cache(maxsize=8192)
def chain_info_df(pdb_id: str) -> pd.DataFrame:
    """
    Fetch the RCSBStructure and related info. Thread-safe operation.
    """
    # There may be an error downloading data for a certain sequence 
    try:
        r = RCSBStructure(pdb_id, download_struct=False)
        df = r.get_chain_info_df().copy()
        df["Chain"] = df["Chain"].apply(lambda x: x.split(","))
        df["Auth"]  = df["Auth"].apply(lambda x: x.split(","))
        return df.explode(["Chain", "Auth"], ignore_index=True)
    except Exception as e:
        return pd.DataFrame(columns=["Chain","Auth"])

@lru_cache(maxsize=8192)
def row_to_record(row) -> dict:
    """
    Convert row to record using cached per-PDB table
    """
    df = chain_info_df(row.entry_id)
    chain1 = df.loc[df["Chain"] == row.asym_id_1]
    chain2 = df.loc[df["Chain"] == row.asym_id_2]
    if chain1.empty or chain2.empty:
        return {"entry_id": row.entry_id, "asym_id_1": row.asym_id_1,
                "asym_id_2": row.asym_id_2, "uniprot_1": np.nan,
                "uniprot_2": np.nan, "pdb_seq_1": np.nan, "pdb_seq_2": np.nan}
    chain1, chain2 = chain1.iloc[0], chain2.iloc[0]
    return {
        "entry_id": row.entry_id,
        "asym_id_1": row.asym_id_1,
        "asym_id_2": row.asym_id_2,
        "uniprot_1": chain1.get("UniProt IDs", np.nan),
        "uniprot_2": chain2.get("UniProt IDs", np.nan),
        "pdb_seq_1": chain1.get("PDB Seq", np.nan),
        "pdb_seq_2": chain2.get("PDB Seq", np.nan),
    }
    
def setup_process_logger(out_dir: str, chunk_idx: int) -> logging.Logger:
    """
    Give each process its own file handler so .info lines don't collide.
    """
    lg = logging.getLogger(f"pinder.process.{chunk_idx:04d}")
    lg.propagate = False
    lg.setLevel(logging.INFO)
    # Clear any old handlers
    for h in list(lg.handlers):
        lg.removeHandler(h)
        
    logger_path = Path(HydraConfig.get().runtime.output_dir) / f"worker_{chunk_idx:04d}.log"
    fh = logging.FileHandler(logger_path)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    lg.addHandler(fh)
    return lg, logger_path
    
def process_pairs(df_pairs: pd.DataFrame, max_workers: int = 8, tqdm_file: Optional[TextIO] = None,
    tqdm_desc: Optional[str] = "Fetching chain info") -> pd.DataFrame:
    """
    Multithreaded processing of data, instead of using df.apply
    Keeping max_workers modest (e.g. 4-8) is ideal to avoid rate-limiting. 
    Threads will share requests-cache 
    """
    rows = list(df_pairs[["entry_id","asym_id_1","asym_id_2"]]
                .itertuples(index=False, name="R"))
    out = thread_map(
        row_to_record, rows,
        max_workers=max_workers,
        chunksize=32,
        desc=tqdm_desc,
        mininterval=1.0,
        smoothing=0.1,
        # >>> all below are forwarded to tqdm() <<<
        file=tqdm_file,            # <<< send bar to per-process log file
        dynamic_ncols=False,       # stable width in files
        ascii=True,                # cleaner in logs
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
        leave=True
    )
    return pd.DataFrame(out)

def setup_requests_for_process():
    """Ensure each process has its own session + patched requests.get."""
    global sess, _orig_get
    
    # unique cache name per PID
    requests_cache.install_cache(
        f"rcsb_cache_{os.getpid()}",
        expire_after=24*3600, fast_save=True, stale_if_error=True
    )
    
    # Rebuild session per process to avoid cross-process state issues
    retry = Retry(
        total=6, connect=6, read=6, status=6,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset({"GET","HEAD","OPTIONS"})
    )
    sess_local = requests.Session()
    adapter = HTTPAdapter(max_retries=retry, pool_connections=100, pool_maxsize=100)
    sess_local.mount("https://", adapter); sess_local.mount("http://", adapter)
    sess_local.headers.update({"User-Agent": "pinder-batch/1.0 (contact: you@upenn.edu)"})

    # Use the same rate limiter globals (_lock/_last) defined above
    def robust_get_local(url, **kwargs):
        timeout = kwargs.pop("timeout", (5, 30))
        _rate_wait()
        return sess_local.get(url, timeout=timeout, **kwargs)
    requests.get = robust_get_local  # monkey-patch for this process
    return sess_local

def _split_disjoint(
    df: pd.DataFrame,
    per_chunk: int,
    n_chunks: int,
    seed: int = 42,
    cover_all: bool = True,
) -> List[pd.DataFrame]:
    """
    Shuffle once, then split into up to `n_chunks` pieces.

    If cover_all=True:
        - Assign *all* rows.
        - Create up to `n_chunks` chunks; the last chunk takes the remainder.
        - If total rows > n_chunks*per_chunk, the last chunk may be larger than `per_chunk`.

    If cover_all=False (original behavior):
        - Create only full chunks of exactly `per_chunk` rows (up to `n_chunks`).
        - Remainder rows are ignored.
        - Raises ValueError if there aren't enough rows to make even one full chunk.

    Returns a list of copies of DataFrame slices.
    """
    # Basic guards
    if len(df) == 0 or n_chunks <= 0 or per_chunk <= 0:
        return []

    # Deterministic shuffle, then slice
    df_shuf = df.sample(frac=1.0, random_state=seed, replace=False).reset_index(drop=True)
    n_rows = len(df_shuf)

    if cover_all:
        needed = math.ceil(n_rows / per_chunk) if per_chunk > 0 else 0
        use = min(n_chunks, needed) if n_chunks > 0 else 0
        chunks: List[pd.DataFrame] = []
        for i in range(use):
            start = i * per_chunk
            # last chunk takes the remainder (could be < per_chunk or > per_chunk if rows > n_chunks*per_chunk)
            end = (start + per_chunk) if i < use - 1 else n_rows
            if start >= n_rows:
                break
            chunks.append(df_shuf.iloc[start:end].copy())
        return chunks

    # --- original strict behavior: only full chunks, ignore remainder ---
    n_full = min(n_chunks, n_rows // per_chunk)
    if n_full == 0:
        raise ValueError(f"Not enough rows ({n_rows}) to form any {per_chunk}-row chunk.")
    chunks: List[pd.DataFrame] = []
    for i in range(n_full):
        start = i * per_chunk
        end = start + per_chunk
        chunks.append(df_shuf.iloc[start:end].copy())
    return chunks

def _process_chunk(chunk_df: pd.DataFrame, chunk_idx: int, out_dir: str, max_workers_inner: int = 1, chunktype="neg") -> str:
    """Worker run in a separate process: fetch chain info and write a per-chunk CSV."""
    setup_requests_for_process()  # per-process session + patched requests.get
    plogger, plogger_path = setup_process_logger(out_dir, chunk_idx)
    
    # Open a per-process PROGRESS file for tqdm
    with open(plogger_path, "a", buffering=1) as tqdm_fh:  # line-buffered
        t0 = time.time()
        df = process_pairs(
            chunk_df,
            max_workers=max_workers_inner,
            tqdm_file=tqdm_fh,
            tqdm_desc=f"chunk {chunk_idx:04d}"
        )
        path = Path(out_dir) / f"rcsb_info_pinder_{chunktype}_chunk_{chunk_idx:04d}.csv"
        df.to_csv(path, index=False)
        plogger.info(f"[chunk {chunk_idx:04d}] wrote {len(df)} rows to {path.name} in {time.time()-t0:.1f}s")
    return str(path)

def run_multiprocess_outer(chunks, per_core, out_dir, inner_threads, fname="rcsb_info_pinder_no_prodigy_label.csv", chunktype="neg", keep_chunks=False):
    n_procs = len(chunks)
    logger.info(f"Launching {n_procs} process(es), {per_core} negatives each (total {n_procs * per_core}).")
    
    # Fan out across processes
    t_start = time.time()
    out_paths = []
    with ProcessPoolExecutor(max_workers=n_procs) as ex:
        futures = {
            ex.submit(_process_chunk, chunk, i, str(out_dir), inner_threads, chunktype): i
            for i, chunk in enumerate(chunks)
        }
        for fut in as_completed(futures):
            out_paths.append(fut.result())
    
    # Concatenate all chunk CSVs → one final CSV (optional but convenient)
    final_path = out_dir / fname
    all_df = pd.concat([pd.read_csv(p) for p in sorted(out_paths)], ignore_index=True)
    all_df.to_csv(final_path, index=False)

    logger.info(f"Wrote {len(all_df)} rows total to {final_path} in {time.time() - t_start:.1f}s")
    
    # Cleanup
    if not keep_chunks:
        deleted = 0
        for p in sorted(set(out_paths)):
            pp = Path(p)
            try:
                # Safety: only delete the files we just created and that live in out_dir
                if pp.parent == out_dir and (pp.name.startswith("rcsb_info_pinder_neg_chunk_") or pp.name.startswith("rcsb_info_pinder_pos_chunk_")) and pp.suffix == ".csv":
                    pp.unlink()
                    deleted += 1
            except FileNotFoundError:
                pass
            except Exception as e:
                logger.warning(f"Could not delete chunk file {pp}: {e}")
        logger.info(f"Deleted {deleted} temporary chunk CSV(s).")
    else:
        logger.info("Keeping temporary chunk CSVs per config (process.keep_chunk_csvs=True)")

def main(cfg: DictConfig):
    """
    Main method. Process PINDER data, either positives or negatives or both.
    """
    mode = cfg.process.mode
    logger.info(f"Running IntAct processor in mode: {mode}")

    # If we haven't already downloaded PINDER metadata, download it 
    if not(os.path.exists(cfg.process.input_metadata)):
        metadata = get_metadata()
    else:
        metadata = pd.read_csv(Path(root) / cfg.process.input_metadata)
    metadata["resolution"] = metadata["resolution"].astype(float)
    metadata = metadata.replace("", np.nan)

    # Find negatives and drop duplicates
    metadata["temp_id"] = metadata["entry_id"].astype(str) + "_" + metadata["chain1_id"].astype(str) + "_" + metadata["chain2_id"].astype(str) + "_" + metadata["asym_id_1"].astype(str) + "_" + metadata["asym_id_2"].astype(str)
    negatives = (
        metadata.loc[(metadata["label"].isna()) & (metadata["resolution"] <= 3.5)]
        .drop_duplicates(
            subset=["temp_id"]
        )
        .reset_index(drop=True)
    )
    logger.info(f"total negatives: {len(negatives):,}")
    total_negs_no_req = len(
        metadata.loc[metadata["label"].isna()]
        .drop_duplicates(
            subset=["temp_id"]
        )
        .reset_index(drop=True)
    )
    logger.info(f"total negatives without resolution req: {total_negs_no_req:,}")
    assert (negatives["probability"] == 0).all()

    #### positives
    assert len(metadata) == len(metadata["id"].unique())
    positives = (
        metadata.loc[(~metadata["label"].isna()) & (metadata["resolution"] <= 3.5)]
        .drop_duplicates(
            subset=["entry_id", "chain1_id", "chain2_id", "asym_id_1", "asym_id_2"]
        )
        .reset_index(drop=True)
    )
    logger.info(f"\ntotal positives: {len(positives):,}")
    total_pos_no_req = len(
        metadata.loc[~metadata["label"].isna()]
        .drop_duplicates(
            subset=["entry_id", "chain1_id", "chain2_id", "asym_id_1", "asym_id_2"]
        )
        .reset_index(drop=True)
    )
    logger.info(f"total positives without resolution req: {total_pos_no_req:,}")
    
    ### make sure that negatives and positives don't have overlap of same entry and same chains
    negatives_temp_ids = negatives["temp_id"].unique().tolist()
    positives_temp_ids = positives["temp_id"].unique().tolist()
    negatives = negatives.loc[
        (~negatives["temp_id"].isin(positives_temp_ids))
    ].reset_index(drop=True)
    positives = positives.loc[
        (~positives["temp_id"].isin(negatives_temp_ids))
    ].reset_index(drop=True)
    logger.info(f"After eliminating overlap of entry_id+chain1_id+chain2_id+asym_id_1+asym_id_2:")
    logger.info(f"Negatives size: {len(negatives)}")
    logger.info(f"Positives size: {len(positives)}")
    
    # make the folder for saving final files
    os.makedirs(cfg.process.output_dir, exist_ok=True)
    
    #### negatives
    per_core = getattr(cfg.process, "per_core_samples", 1000)
    req_cores = getattr(cfg.process, "n_cores", mp.cpu_count() or 1)
    inner_threads = getattr(cfg.process, "max_workers_inner", 1)  # keep low to avoid API rate limits
    keep_chunks = getattr(cfg.process, "keep_chunks", False)
    
    # Define out dir
    out_dir = Path(root) / cfg.process.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Run the multiprocess for all negatives or all positives depending on what the user asked for 
    if mode in (["all","negative"]):
        # Split negatives into chunks
        chunks = _split_disjoint(negatives, per_chunk=per_core, n_chunks=req_cores, seed=42, cover_all=cfg.process.cover_all)
        run_multiprocess_outer(chunks, per_core, out_dir, inner_threads, fname="rcsb_info_pinder_no_prodigy_label.csv", chunktype="neg", keep_chunks=keep_chunks)
    if mode in (["all","positive"]):
        # Split positives into chunks
        chunks = _split_disjoint(positives, per_chunk=per_core, n_chunks=req_cores, seed=42, cover_all=cfg.process.cover_all)
        run_multiprocess_outer(chunks, per_core, out_dir, inner_threads, fname="rcsb_info_pinder_bioxtal_prodigy_label.csv", chunktype="pos", keep_chunks=keep_chunks)