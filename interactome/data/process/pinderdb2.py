"""
Code for processing the PINDER database, extracting both positives and negatives.
"""
from pinder.core import get_metadata
import numpy as np
import pandas as pd
from protparser.rcsb import RCSBStructure
import time

import os, time, random, logging, json
import numpy as np
import pandas as pd
import requests, urllib3
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from functools import lru_cache
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
log = logging.getLogger("pinder_full")

# -----------------------------
# Make requests.get robust even for libraries that call it directly
# -----------------------------
def _make_retrying_session(total=5, backoff=1.0):
    retry = Retry(
        total=total,
        connect=total,
        read=total,
        status=total,
        backoff_factor=backoff,                   # 1, 2, 4, 8...
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET","POST","PUT","DELETE","HEAD","OPTIONS","TRACE"])
    )
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=retry, pool_connections=100, pool_maxsize=100)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    return s

_retry_session = _make_retrying_session()

# monkey-patch requests.get so protparser benefits from retry+timeout
_orig_get = requests.get
def _robust_get(url, **kwargs):
    # connection timeout: how long we'll wait to establish the TCP connection 
    # read timeout: how long we'll wait for the server to send data fater we've connected
    timeout = kwargs.pop("timeout", (5, 30))      # (connect, read) seconds
    return _retry_session.get(url, timeout=timeout, **kwargs)
requests.get = _robust_get  # <-- affects protparser.rcsb internally

# -----------------------------
# Backoff helpers for outer retries around RCSBStructure
# -----------------------------
RETRY_EXCS = (
    requests.exceptions.RequestException,
    urllib3.exceptions.HTTPError,
    TimeoutError,
    ConnectionError,
    ValueError,  # e.g., .item() failures if row missing; we will handle gracefully
)

def _backoff_delay(attempt, base=1.0, cap=20.0, jitter=0.25):
    d = min(cap, base * (2 ** (attempt - 1)))
    return d * (1.0 + jitter * random.random())

@lru_cache(maxsize=8192)
def fetch_chain_info_df(pdb_id: str, max_tries: int = 6):
    """
    Fetch & cache per-PDB chain info
    """
    last = None
    for attempt in range(1, max_tries + 1):
        try:
            rcsb = RCSBStructure(pdb_id, download_struct=False)
            df = rcsb.get_chain_info_df().copy()
            # Normalize multi-valued fields to rows
            for col in ("Chain", "Auth"):
                if col in df.columns:
                    df[col] = df[col].apply(lambda x: x.split(",") if isinstance(x, str) else x)
            df = df.explode(["Chain","Auth"], ignore_index=True)
            return df
        except RETRY_EXCS as e:
            last = e
            wait = _backoff_delay(attempt)
            log.warning(f"[{pdb_id}] attempt {attempt} failed: {e}. Retrying in {wait:.1f}s")
            time.sleep(wait)
    raise last

def get_chain_sequences(pdb_id: str, chain1: str, chain2: str, max_tries: int = 6):
    df = fetch_chain_info_df(pdb_id, max_tries=max_tries)
    # robust selection; sometimes there are duplicates or missing rows
    c1 = df.loc[df["Chain"] == chain1]
    c2 = df.loc[df["Chain"] == chain2]
    if c1.empty or c2.empty:
        raise KeyError(f"Missing chain(s) in {pdb_id}: wanted ({chain1},{chain2}), "
                       f"got c1={len(c1)} c2={len(c2)}")
    # be tolerant to duplicates: take first row
    c1 = c1.iloc[0]
    c2 = c2.iloc[0]
    return {
        "entry_id": pdb_id,
        "asym_id_1": chain1,
        "asym_id_2": chain2,
        "uniprot_1": c1.get("UniProt IDs", np.nan),
        "uniprot_2": c2.get("UniProt IDs", np.nan),
        "pdb_seq_1": c1.get("PDB Seq", np.nan),
        "pdb_seq_2": c2.get("PDB Seq", np.nan),
    }

def process_pairs(df: pd.DataFrame, out_csv: str, fail_csv: str,
                  max_tries: int = 6, checkpoint_every: int = 2000):
    """
    Processing helpers (full DB, checkpointing)
    """
    # Make directory to save results 
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    
    # Initialize results
    results = []
    failures = []
    n = len(df)
    header_written = os.path.exists(out_csv) and os.path.getsize(out_csv) > 0

    def flush():
        nonlocal results, failures, header_written
        if results:
            out_df = pd.DataFrame(results)
            out_df.to_csv(out_csv, mode="a", header=not header_written, index=False)
            header_written = True
            results = []
        if failures:
            pd.DataFrame(failures).to_csv(fail_csv, mode="a",
                                          header=not os.path.exists(fail_csv), index=False)
            failures = []

    for i, row in tqdm(df.iterrows(), total=n, desc=f"processing → {out_csv}"):
        try:
            rec = get_chain_sequences(row.entry_id, row.asym_id_1, row.asym_id_2, max_tries=max_tries)
            results.append(rec)
        except Exception as e:
            failures.append({
                "entry_id": row.entry_id,
                "asym_id_1": row.asym_id_1,
                "asym_id_2": row.asym_id_2,
                "error": repr(e),
            })
        if (i + 1) % checkpoint_every == 0:
            flush()
    flush()

# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    # Load metadata once
    md = get_metadata()
    md["resolution"] = md["resolution"].astype(float)
    md = md.replace("", np.nan)

    key_subset = ["entry_id","chain1_id","chain2_id","asym_id_1","asym_id_2"]
    # NEGATIVES
    negatives = (
        md.loc[md["label"].isna() & (md["resolution"] <= 3.5)]
          .drop_duplicates(subset=key_subset)
          .reset_index(drop=True)
    )
    log.info(f"total negatives: {len(negatives):,} "
             f"(no-res-filter={len(md.loc[md['label'].isna()].drop_duplicates(subset=key_subset)):,})")
    # POSITIVES
    positives = (
        md.loc[~md["label"].isna() & (md["resolution"] <= 3.5)]
          .drop_duplicates(subset=key_subset)
          .reset_index(drop=True)
    )
    log.info(f"total positives: {len(positives):,} "
             f"(no-res-filter={len(md.loc[~md['label'].isna()].drop_duplicates(subset=key_subset)):,})")

    # Process full datasets with retries
    # (use .itertuples for speed; we access by attribute above)
    negatives = negatives[key_subset].itertuples(index=False, name="Row")
    #positives = positives[key_subset].itertuples(index=False, name="Row")

    # Convert to DataFrames again for our iterrows loop (tqdm needs length)
    negatives_df = pd.DataFrame(list(negatives), columns=key_subset)
    #positives_df = pd.DataFrame(list(positives), columns=key_subset)

    process_pairs(
        negatives_df,
        out_csv="rcsb_info_pinder_no_prodigy_label.csv",
        fail_csv="rcsb_info_pinder_neg_failures.csv",
        max_tries=6,
        checkpoint_every=2000,
    )
    if False:
        process_pairs(
            positives_df,
            out_csv="rcsb_info_pinder_bioxtal_prodigy_label.csv",
            fail_csv="rcsb_info_pinder_pos_failures.csv",
            max_tries=6,
            checkpoint_every=2000,
        )