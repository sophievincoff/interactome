from __future__ import annotations
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path
import os
import pandas as pd
import logging
import re
import numpy as np
from typing import Optional, Dict, Any, Iterable, Tuple, List
import concurrent.futures as cf
import math
import traceback
import multiprocessing as mp
import concurrent.futures as cf
from pathlib import Path
from pandas.api.types import is_scalar
from typing import Iterable, Dict, Any, List, Optional

from protparser.rcsb import download_rcsb

# Chain-aware mmCIF parsing
import gemmi

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

# Natural 20 + selenocysteine U
VALID_AAS = set("ACDEFGHIKLMNPQRSTVWYU")
DISALLOWED_ALWAYS = set(["X", "O"])  # X unknown, O pyrrolysine (you said don't keep)

def harmonize_nulls_to_nan(df: pd.DataFrame, *, also_blank_strings=True, keep_datetime=False) -> pd.DataFrame:
    out = df.copy()

    # 1) Convert common sentinels to real missing
    if also_blank_strings:
        out = out.replace({"": pd.NA, "None": pd.NA, "nan": pd.NA})

    # 2) Normalize to pandas NA first (unifies None/NaN/<NA>)
    out = out.convert_dtypes()

    # 3) Cast extension dtypes -> object so np.nan can live there.
    for c in out.columns:
        dt = out[c].dtype
        is_ext = isinstance(dt, pd.api.extensions.ExtensionDtype)
        if keep_datetime and pd.api.types.is_datetime64_any_dtype(dt):
            # keep datetimes as datetime64 with NaT
            continue
        if is_ext:
            out[c] = out[c].astype(object)

    # 4) Finally: make ALL missings = np.nan
    out = out.where(~out.isna(), np.nan)

    return out

def get_unique_id(row, colA="ID(s) interactor A", colB="ID(s) interactor B"):
    """
    Create a unique ID for the pair of interactors in the row, so that order does not matter
    """
    intA = row[colA]
    intB = row[colB]
    
    if intA is None or (type(intA)==float and np.isnan(intA)):
        intA=""
    if intB is None or (type(intB)==float and np.isnan(intB)):
        intB=""
    
    if intA <= intB:
        return f"{intA}_{intB}"
    return f"{intB}_{intA}"

# before we save merged, must correct invalid aas
def find_invalid_chars(seq: str, valid_chars: set) -> set:
    """
    Find and return a set of invalid characters in a sequence.

    Args:
        seq (str): The sequence you wish to search for invalid characters.
        valid_chars (set): A set of valid characters.

    Returns:
        set: A set of characters in the sequence that are not in the set of valid characters.
    """
    if seq is None or (isinstance(seq, float) and pd.isna(seq)):
        return np.nan
    s = re.sub(r"\s+", "", str(seq)).upper()
    if s == "":
        return np.nan
    unique_chars = set(s)
    if unique_chars.issubset(valid_chars):
        return np.nan
    bad = sorted(unique_chars.difference(valid_chars))
    return ",".join(bad)

def get_unique_id(row, colA="ID(s) interactor A", colB="ID(s) interactor B"):
    """
    Create a unique ID for the pair of interactors in the row, so that order does not matter
    """
    intA = row[colA]
    intB = row[colB]
    
    if intA is None or (type(intA)==float and np.isnan(intA)):
        intA=""
    if intB is None or (type(intB)==float and np.isnan(intB)):
        intB=""
    
    if intA <= intB:
        return f"{intA}_{intB}"
    return f"{intB}_{intA}"

def initial_filtering(
    complexes: pd.DataFrame,
    pepnn_test: pd.DataFrame,
    pepnn_train: pd.DataFrame,
    propedia_test: pd.DataFrame 
    ):
    complexes["Resolution"] = complexes["Resolution"].astype(float)
    
    # Find invalids
    pepnn_test["invalid_binder_chars"] = pepnn_test["Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))
    pepnn_test["invalid_receptor_chars"] = pepnn_test["Receptor Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))

    pepnn_train["invalid_binder_chars"] = pepnn_train["Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))
    pepnn_train["invalid_receptor_chars"] = pepnn_train["Receptor Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))

    propedia_test["invalid_binder_chars"] = propedia_test["Binder"].apply(lambda x: find_invalid_chars(x,VALID_AAS))
    propedia_test["invalid_receptor_chars"] = propedia_test["Receptor Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))

    complexes["invalid_binder_chars"] = complexes["Peptide Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))
    complexes["invalid_receptor_chars"] = complexes["Receptor Sequence"].apply(lambda x: find_invalid_chars(x,VALID_AAS))

    # Get seq_sort
    pepnn_test["seq_sort"] = pepnn_test.apply(lambda row: get_unique_id(row,colA="Sequence",colB="Receptor Sequence"),axis=1)
    pepnn_train["seq_sort"] = pepnn_train.apply(lambda row: get_unique_id(row,colA="Sequence",colB="Receptor Sequence"),axis=1)
    propedia_test["seq_sort"] = propedia_test.apply(lambda row: get_unique_id(row,colA="Binder",colB="Receptor Sequence"),axis=1)
    complexes["seq_sort"] = complexes.apply(lambda row: get_unique_id(row,colA="Peptide Sequence",colB="Receptor Sequence"),axis=1)

    # get unique sequence pairs in each, with and without X's 
    test1 = len(pepnn_test["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in PepNN-test: {test1}")
    test1 = len(pepnn_test.loc[
        (pepnn_test["invalid_binder_chars"].isna()) &
        (pepnn_test["invalid_receptor_chars"].isna())
    ]["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in PepNN-test with no invalid chars: {test1}")

    test1 = len(pepnn_train["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in PepNN-train: {test1}")
    test1 = len(pepnn_train.loc[
        (pepnn_train["invalid_binder_chars"].isna()) &
        (pepnn_train["invalid_receptor_chars"].isna())
        
    ]["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in PepNN-train with no invalid chars: {test1}")

    test1 = len(propedia_test["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in Propedia-test: {test1}")
    test1 = len(propedia_test.loc[
        (propedia_test["invalid_binder_chars"].isna()) &
        (propedia_test["invalid_receptor_chars"].isna())
        
    ]["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in Propedia-test with no invalid chars: {test1}")

    test1 = len(complexes["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in Propedia-test: {test1}")
    test1 = len(complexes.loc[
        (complexes["invalid_binder_chars"].isna()) &
        (complexes["invalid_receptor_chars"].isna())
        
    ]["seq_sort"].unique())
    logger.info(f"Unique sequence pairs in Propedia-test with no invalid chars: {test1}")
    
    # get total unique from all of these
    unique_seq_pairs = pd.concat(
        [
            pepnn_test.loc[
                (pepnn_test["invalid_binder_chars"].isna()) &
                (pepnn_test["invalid_receptor_chars"].isna())
            ].rename(columns={"Sequence":"Binder Sequence",
                            "Receptor": "Receptor Chain",
                            "Peptide": "Peptide Chain"})[["PDB","Binder Sequence","Receptor Sequence","Receptor Chain","Peptide Chain","seq_sort"]],
            pepnn_train.loc[
                (pepnn_train["invalid_binder_chars"].isna()) &
                (pepnn_train["invalid_receptor_chars"].isna())
            ].rename(columns={"Sequence":"Binder Sequence",
                            "Receptor": "Receptor Chain",
                            "Peptide": "Peptide Chain"})[["PDB","Binder Sequence","Receptor Sequence","Receptor Chain","Peptide Chain","seq_sort"]],
            propedia_test.loc[
                (propedia_test["invalid_binder_chars"].isna()) &
                (propedia_test["invalid_receptor_chars"].isna())
            ].rename(columns={"Binder":"Binder Sequence"})[["PDB","Binder Sequence","Receptor Sequence","seq_sort"]],
            complexes.loc[
                (complexes["invalid_binder_chars"].isna()) &
                (complexes["invalid_receptor_chars"].isna())
            ].rename(columns={"Peptide Sequence":"Binder Sequence",})[["PDB","Binder Sequence","Receptor Sequence","Receptor Chain","Peptide Chain","seq_sort"]],
        ]
    ).drop_duplicates().reset_index(drop=True)
    unique_seq_pairs["Binder Length"] = unique_seq_pairs["Binder Sequence"].apply(lambda x: len(x))
    unique_seq_pairs["Receptor Length"] = unique_seq_pairs["Receptor Sequence"].apply(lambda x: len(x))
    print(f"Total unique sequence pairs that are completely valid across all PepMLM constituent databases: {len(unique_seq_pairs)}")
    
    return unique_seq_pairs

def _normalize_chain(chain: Optional[str]) -> Optional[str]:
    if chain is None:
        return None
    s = str(chain).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "-"}:
        return None
    return s

def _seq_is_clean(seq: str) -> Tuple[bool, Optional[str]]:
    """
    Returns (ok, reason_if_not_ok)
    """
    seq = re.sub(r"\s+", "", seq).upper()
    if not seq:
        return False, "no_sequence"

    bad = sorted(set(seq) - VALID_AAS)
    if any(ch in DISALLOWED_ALWAYS for ch in seq):
        # prefer explicit reason for X/O
        if "X" in seq:
            return False, "contains_X"
        if "O" in seq:
            return False, "contains_O"

    if bad:
        # includes any letters besides allowed set (including B, Z, J, etc.)
        return False, f"contains_non_allowed_letters:{''.join(bad)}"

    return True, None

def _find_downloaded_path(output_dir: Path, pdb_id: str) -> Optional[Path]:
    """
    Return an existing downloaded file path if present (.cif preferred), case-insensitive.
    """
    pdb_id = pdb_id.lower()

    for p in output_dir.rglob("*"):
        if not p.is_file():
            continue
        name = p.stem.lower()
        ext = p.suffix.lower()
        if name == pdb_id and ext in {".cif", ".pdb"}:
            return p

    return None

def _download_one(
    pdb_id: str,
    output_dir: Path,
    primary_format: str = "cif",
) -> Dict[str, Any]:
    """
    Download one structure via protparser.download_rcsb, with format fallback.
    Returns a log row dict.
    """
    pdb_id_norm = str(pdb_id).strip().lower()
    other_format = None if primary_format == "cif" else "cif"

    # If already present, mark as skipped
    existing = _find_downloaded_path(output_dir, pdb_id_norm)
    if existing is not None:
        return {
            "pdb_id": pdb_id_norm,
            "downloaded_path": str(existing),
            "cleaned_path": None,
            "clean": None,
            "reason": "already_downloaded",
        }

    # Try primary, then fallback
    try:
        download_rcsb(
            pdb_id_norm,
            struct_format=primary_format,
            convert_if_fail=False,
            output_dir=str(output_dir),
        )
        existing = _find_downloaded_path(output_dir, pdb_id_norm)
        if existing is None:
            return {
                "pdb_id": pdb_id_norm,
                "downloaded_path": None,
                "cleaned_path": None,
                "clean": None,
                "reason": f"download_succeeded_but_file_missing:{primary_format}",
            }
        return {
            "pdb_id": pdb_id_norm,
            "downloaded_path": str(existing),
            "cleaned_path": None,
            "clean": None,
            "reason": f"downloaded:{primary_format}",
        }
    except Exception as e1:
        try:
            download_rcsb(
                pdb_id_norm,
                struct_format=other_format,
                convert_if_fail=False,
                output_dir=str(output_dir),
            )
            existing = _find_downloaded_path(output_dir, pdb_id_norm)
            if existing is None:
                return {
                    "pdb_id": pdb_id_norm,
                    "downloaded_path": None,
                    "cleaned_path": None,
                    "clean": None,
                    "reason": f"download_succeeded_but_file_missing:{other_format}",
                }
            return {
                "pdb_id": pdb_id_norm,
                "downloaded_path": str(existing),
                "cleaned_path": None,
                "clean": None,
                "reason": f"downloaded_fallback:{other_format}",
            }
        except Exception as e2:
            return {
                "pdb_id": pdb_id_norm,
                "downloaded_path": None,
                "cleaned_path": None,
                "clean": None,
                "reason": f"download_failed:{type(e1).__name__}->{type(e2).__name__}",
            }

def _existing_cif_stems(output_dir: Path) -> set[str]:
    # case-insensitive stems (handles 1A1M.cif vs 1a1m.cif)
    return {p.stem.strip().lower() for p in output_dir.rglob("*.cif")}

def filter_pdbs_to_download(
    pdb_ids: Iterable[str],
    output_dir: Path,
    redownload: bool = False,
) -> List[str]:
    pdb_ids_norm = sorted({str(x).strip().lower() for x in pdb_ids if str(x).strip()})
    if redownload:
        return pdb_ids_norm

    have = _existing_cif_stems(output_dir)
    missing = [p for p in pdb_ids_norm if p not in have]
    logger.info(f"Download skip enabled: have {len(have)} CIFs, need {len(pdb_ids_norm)}, missing {len(missing)}")
    return missing

def _chunkify(items: List[str], n_chunks: int) -> List[List[str]]:
    """Split items into n_chunks as evenly as possible."""
    n_chunks = max(1, int(n_chunks))
    n = len(items)
    if n == 0:
        return [[] for _ in range(n_chunks)]
    # ceil(n / n_chunks) sized blocks
    k = math.ceil(n / n_chunks)
    return [items[i:i+k] for i in range(0, n, k)]

def _download_many_worker(
    pdb_ids_chunk: List[str],
    output_dir: str,
    primary_format: str,
    worker_id: int,
    worker_log_dir: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Worker: download a chunk sequentially and return list of log rows.
    Writes optional per-worker log file (highly recommended on HPC).
    """
    out = []
    log_fp = None
    if worker_log_dir is not None:
        Path(worker_log_dir).mkdir(parents=True, exist_ok=True)
        log_fp = Path(worker_log_dir) / f"download_worker_{worker_id}.log"

    def _wlog(msg: str):
        if log_fp is not None:
            with open(log_fp, "a") as f:
                f.write(msg.rstrip() + "\n")

    _wlog(f"[worker {worker_id}] starting chunk size={len(pdb_ids_chunk)} output_dir={output_dir}")

    # IMPORTANT: do NOT let any exception kill the worker silently
    for j, pdb_id in enumerate(pdb_ids_chunk, 1):
        try:
            row = _download_one(pdb_id, Path(output_dir), primary_format=primary_format)
        except BaseException as e:
            # Catch BaseException in case something calls SystemExit, etc.
            row = {
                "pdb_id": str(pdb_id).strip().lower(),
                "downloaded_path": None,
                "cleaned_path": None,
                "clean": None,
                "reason": f"worker_crash:{type(e).__name__}",
            }
            _wlog(f"[worker {worker_id}] EXC pdb_id={pdb_id} err={type(e).__name__}: {e}")
            _wlog(traceback.format_exc())

        out.append(row)

        if j % 50 == 0 or j == len(pdb_ids_chunk):
            _wlog(f"[worker {worker_id}] progress {j}/{len(pdb_ids_chunk)} last={pdb_id} reason={row.get('reason')}")

    _wlog(f"[worker {worker_id}] done")
    return out

def _inventory_many_worker(
    pdb_ids_chunk: List[str],
    paths_map: Dict[str, str],   # {pdb_id -> cif_path_str}
    worker_id: int,
    worker_log_dir: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Worker: build inventory for a chunk sequentially and return list of result dicts:
      - ok rows: {"ok": True, "pdb_id": ..., "df": <pd.DataFrame>}
      - fail rows: {"ok": False, "pdb_id": ..., "error": "..."}
    """
    out: List[Dict[str, Any]] = []

    log_fp = None
    if worker_log_dir is not None:
        Path(worker_log_dir).mkdir(parents=True, exist_ok=True)
        log_fp = Path(worker_log_dir) / f"inventory_worker_{worker_id}.log"

    def _wlog(msg: str):
        if log_fp is not None:
            with open(log_fp, "a") as f:
                f.write(msg.rstrip() + "\n")

    _wlog(f"[worker {worker_id}] starting chunk size={len(pdb_ids_chunk)}")

    for j, pid in enumerate(pdb_ids_chunk, 1):
        cif_path = Path(paths_map[pid])
        try:
            df = build_chain_inventory_for_cif(cif_path, pid)
            out.append({"ok": True, "pdb_id": pid, "df": df})
        except BaseException as e:
            out.append({"ok": False, "pdb_id": pid, "error": f"{type(e).__name__}: {e}"})
            _wlog(f"[worker {worker_id}] EXC pdb_id={pid} path={cif_path} err={type(e).__name__}: {e}")
            _wlog(traceback.format_exc())

        if j % 50 == 0 or j == len(pdb_ids_chunk):
            last = out[-1]
            status = "ok" if last.get("ok") else "fail"
            _wlog(f"[worker {worker_id}] progress {j}/{len(pdb_ids_chunk)} last={pid} status={status}")

    _wlog(f"[worker {worker_id}] done")
    return out


def download_structures_parallel_chunked(
    pdb_ids: Iterable[str],
    output_dir: str,
    primary_format: str = "cif",
    max_workers: int = 16,
    log_csv: Optional[str] = None,
    worker_log_dir: Optional[str] = None,   # e.g. "raw_data/download_worker_logs"
) -> pd.DataFrame:
    """
    Process-parallel, chunked: split PDB IDs into ~max_workers chunks, each worker loops sequentially.
    More stable than submitting one future per PDB.
    """
    output_dir_p = Path(output_dir)
    output_dir_p.mkdir(parents=True, exist_ok=True)

    # Normalize + deduplicate
    pdb_ids_norm = sorted({str(x).strip().lower() for x in pdb_ids if str(x).strip()})
    if not pdb_ids_norm:
        return pd.DataFrame(columns=["pdb_id", "downloaded_path", "cleaned_path", "clean", "reason"])

    # Choose worker count
    cpu_count = os.cpu_count() or 1
    n_workers = min(int(max_workers), cpu_count, len(pdb_ids_norm))
    chunks = _chunkify(pdb_ids_norm, n_workers)

    rows: List[Dict[str, Any]] = []

    # On some clusters, forcing "spawn" can avoid weird fork issues with native libs.
    # If you *know* fork is fine, you can remove this.
    ctx = mp.get_context("spawn")

    with cf.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        futs = []
        for wid, chunk in enumerate(chunks):
            if not chunk:
                continue
            futs.append(ex.submit(
                _download_many_worker,
                chunk,
                str(output_dir_p),
                primary_format,
                wid,
                worker_log_dir,
            ))

        for fut in cf.as_completed(futs):
            rows.extend(fut.result())

    df = pd.DataFrame(rows).sort_values("pdb_id").reset_index(drop=True)

    if log_csv is not None:
        Path(log_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(log_csv, index=False)

    return df

def _col(cat: gemmi.cif.Table, name: str):
    c = cat.find_column(name)
    if c is None:
        raise KeyError(f"Missing mmCIF column '{name}' in category '{cat.name}'")
    return c

import re
from typing import Optional

def clean_mmcif_seq_text(x: Optional[str]) -> str:
    """
    Clean an mmCIF polymer sequence field (e.g. _entity_poly.pdbx_seq_one_letter_code_can).

    Handles:
      - semicolon-delimited multi-line text blocks
      - embedded newlines/spaces
      - stray leading/trailing semicolons produced by some parsers
    """
    if x is None:
        return ""
    s = str(x)

    # Trim outer whitespace first
    s = s.strip()

    # If the *value* includes mmCIF semicolon block markers, remove them.
    # In practice (e.g., via gemmi), the returned string often literally begins/ends with ';'.
    if s.startswith(";"):
        s = s[1:]  # drop the opening ';'
    # drop a trailing ';' if present after stripping
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]

    # Remove all whitespace/newlines inside the sequence
    s = re.sub(r"\s+", "", s)

    return s


def build_chain_inventory_for_cif(cif_path: Path, pdb_id: str) -> pd.DataFrame:
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()

    asym = block.find_mmcif_category("_struct_asym")
    asym_id_col = _col(asym, "id")
    ent_id_col  = _col(asym, "entity_id")

    ent_poly = block.find_mmcif_category("_entity_poly")
    ent_poly_ent = _col(ent_poly, "entity_id")
    seq_can_col  = _col(ent_poly, "pdbx_seq_one_letter_code_can")
    type_col     = ent_poly.find_column("type")  # optional; may be None

    entity_to_seq: Dict[str, str] = {}
    entity_to_type: Dict[str, str] = {}

    for i in range(len(ent_poly_ent)):
        eid = str(ent_poly_ent[i])
        #seq1 = _clean_mmcif_seq_text(seq_can_col[i])
        #seq = str(seq_can_col[i]).replace("\n", "").replace(" ", "")
        seq = clean_mmcif_seq_text(seq_can_col[i])
        entity_to_seq[eid] = seq
        if type_col is not None:
            entity_to_type[eid] = str(type_col[i])

    rows = []
    for i in range(len(asym_id_col)):
        chain_id = str(asym_id_col[i])
        entity_id = str(ent_id_col[i])

        seq = entity_to_seq.get(entity_id, "")
        ok, why = _seq_is_clean(seq) if seq else (False, "no_sequence")

        invalid = ""
        if not ok:
            s = re.sub(r"\s+", "", seq).upper() if seq else ""
            invalid_set = sorted(set(s) - VALID_AAS)
            invalid = "".join(invalid_set) if invalid_set else (why or "")
            
        resolved_pos = get_resolved_positions_from_atom_site(cif_path, chain_id)
        seq_res_masked, seq_res_only, n_res, frac_res = build_resolved_sequences(seq, resolved_pos)

        rows.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "entity_id": entity_id,
            "seq_can": seq,
            "seq_len": len(seq) if seq else 0,
            "is_valid_aa": bool(ok),
            "invalid_letters": invalid,

            # resolved part
            "seq_resolved_masked": seq_res_masked,
            "seq_resolved_only": seq_res_only,
            "n_resolved": n_res,
            "frac_resolved": frac_res,
        })

    return pd.DataFrame(rows)

def list_cifs_in_dir(output_dir: Path) -> Dict[str, Path]:
    """
    Return {pdb_id_lower: path} for all *.cif found under output_dir.
    If duplicates exist, keeps the first one encountered.
    """
    out = {}
    for p in output_dir.rglob("*.cif"):
        pid = p.stem.strip().lower()
        out.setdefault(pid, p)
    return out

def _parse_chain_set(x) -> set[str]:
    """
    Turns 'A', 'A,B', 'A B', 'A;B' into {'A','B'}.
    Returns empty set if missing/invalid.
    """
    x = _normalize_chain(x)
    if x is None:
        return set()
    parts = re.split(r"[,\s;/|]+", x.strip())
    return {p for p in (s.strip() for s in parts) if p}

def _norm_chain_list(chains) -> set[str]:
    if not isinstance(chains, list):
        return set()
    out = set()
    for c in chains:
        c2 = _normalize_chain(c)
        if c2 is not None:
            out.add(c2)
    return out


def build_full_chain_inventory(
    download_log: pd.DataFrame,
    inventory_out: Path,
    output_dir: Path,
    max_workers: int = 16,
    append: bool = True,
) -> pd.DataFrame:
    """
    Build chain inventory for all CIFs we have access to.
    - Uses any downloaded_path values from download_log (if present)
    - ALSO scans output_dir for existing CIFs (critical when skipping downloads)
    - Optionally appends to an existing inventory CSV instead of overwriting
    """

    # 1) Load existing inventory if desired
    existing = None
    if append and inventory_out.exists():
        existing = pd.read_csv(inventory_out)
        # normalize key
        if "pdb_id" in existing.columns:
            existing["pdb_id"] = existing["pdb_id"].astype(str).str.lower()

    # 2) Collect CIF paths from BOTH sources
    paths: Dict[str, Path] = {}

    # From download_log (if it has usable paths)
    if download_log is not None and len(download_log) > 0 and "downloaded_path" in download_log.columns:
        for _, r in download_log.iterrows():
            p = r.get("downloaded_path", None)
            if p is None or (isinstance(p, float) and pd.isna(p)):
                continue
            pp = Path(p)
            if pp.exists():
                paths[str(r["pdb_id"]).strip().lower()] = pp

    # From output_dir scan (covers already-present CIFs)
    paths.update(list_cifs_in_dir(output_dir))

    if not paths:
        raise RuntimeError(
            "No CIF files found to build inventory. "
            "Check output_dir and download_log. (Maybe output_dir is wrong.)"
        )

    # 3) Decide which pdb_ids actually need to be (re)built
    to_build = sorted(paths.keys())
    if existing is not None and "pdb_id" in existing.columns:
        already_have = set(existing["pdb_id"].astype(str).str.lower().unique())
        # Only build new ones; if you want to rebuild everything, remove this filter
        to_build = [pid for pid in to_build if pid not in already_have]

    logger.info(f"Building chain inventory for {len(to_build)} CIFs (append={append}).")

    # 4) Build new inventory rows (parallel chunked, like downloads)
    if not to_build:
        logger.info("[inventory] nothing new to build")
        inv = existing if existing is not None else pd.DataFrame()
        inv.to_csv(inventory_out, index=False)
        return inv

    # Make a simple serializable dict for workers (no Paths)
    paths_map = {pid: str(paths[pid]) for pid in to_build}

    cpu_count = os.cpu_count() or 1
    n_workers = min(max_workers, cpu_count, len(to_build))
    chunks = _chunkify(to_build, n_workers)

    logger.info(f"[inventory] building {len(to_build)} CIFs with n_workers={n_workers}")

    # where to put logs (mirror downloads)
    hydra_dir = HydraConfig.get().run.dir
    worker_log_dir = str(Path(hydra_dir) / "process_worker_logs")

    # if gemmi/native libs are weird under fork, keep spawn like downloads
    ctx = mp.get_context("spawn")

    results: List[Dict[str, Any]] = []
    with cf.ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        futs = []
        for wid, chunk in enumerate(chunks):
            if not chunk:
                continue
            futs.append(ex.submit(
                _inventory_many_worker,
                chunk,
                paths_map,
                wid,
                worker_log_dir,
            ))

        done = 0
        total_futs = len(futs)
        for fut in cf.as_completed(futs):
            results.extend(fut.result())
            done += 1
            logger.info(f"[inventory] finished worker {done}/{total_futs}")

    # Collate
    new_parts = [r["df"] for r in results if r.get("ok")]
    failures = [r for r in results if not r.get("ok")]

    if failures:
        fail_out = inventory_out.with_suffix(".inventory_failures.csv")
        pd.DataFrame(failures).to_csv(fail_out, index=False)
        logger.warning(f"[inventory] failures={len(failures)} wrote {fail_out}")

    new_inv = pd.concat(new_parts, ignore_index=True) if new_parts else pd.DataFrame()

    # 5) Combine + deduplicate + save
    if existing is not None and not existing.empty:
        inv = pd.concat([existing, new_inv], ignore_index=True)
    else:
        inv = new_inv

    # Dedup at the “chain” level (safe key)
    if not inv.empty:
        inv["pdb_id"] = inv["pdb_id"].astype(str).str.lower()
        inv["chain_id"] = inv["chain_id"].astype(str)
        inv = inv.drop_duplicates(subset=["pdb_id", "chain_id"], keep="last")

    inventory_out.parent.mkdir(parents=True, exist_ok=True)
    inv.to_csv(inventory_out, index=False)
    return inv
 
def match_interactions_to_inventory(
    combined_df: pd.DataFrame,
    inventory_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Match interaction sequences against BOTH:
      - full canonical chain sequence (seq_can)
      - resolved-only chain sequence (seq_resolved_only)

    Preference order per side: full -> resolved.

    Adds:
      receptor_match_type: full|resolved|none
      binder_match_type:   full|resolved|none
      interaction_match_type: "<receptor>/<binder>"
    """

    df = combined_df.copy()

    # ---- normalize interaction-side ----
    df["pdb_id"] = df["PDB"].astype(str).str.strip().str.lower()
    df["Receptor Sequence"] = (
        df["Receptor Sequence"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    )
    df["Binder Sequence"] = (
        df["Binder Sequence"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    )

    # ---- validate inventory early (before any column access) ----
    if inventory_df is None or inventory_df.empty:
        raise RuntimeError(
            "inventory_df is empty (likely inventory build failed). "
            "Check *.inventory_failures.csv and worker logs."
        )
    if "is_valid_aa" not in inventory_df.columns:
        raise RuntimeError(
            f"inventory_df missing 'is_valid_aa'. Columns={list(inventory_df.columns)}"
        )
    if "seq_resolved_only" not in inventory_df.columns:
        raise RuntimeError(
            f"inventory_df missing 'seq_resolved_only'. Columns={list(inventory_df.columns)}"
        )
    if "seq_can" not in inventory_df.columns or "pdb_id" not in inventory_df.columns or "chain_id" not in inventory_df.columns:
        raise RuntimeError(
            f"inventory_df missing required columns. Columns={list(inventory_df.columns)}"
        )

    # ---- normalize + restrict inventory ----
    inv = inventory_df[inventory_df["is_valid_aa"]].copy()
    inv["pdb_id"] = inv["pdb_id"].astype(str).str.strip().str.lower()
    inv["chain_id"] = inv["chain_id"].astype(str).str.strip()

    # drop junk chain ids
    inv = inv.dropna(subset=["pdb_id", "chain_id"])
    inv = inv[inv["chain_id"].ne("") & inv["chain_id"].str.lower().ne("nan")]

    # normalize sequences
    inv["seq_can"] = inv["seq_can"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    inv["seq_resolved_only"] = (
        inv["seq_resolved_only"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    )

    # drop rows missing the full sequence
    inv = inv.dropna(subset=["seq_can"])
    inv = inv[inv["seq_can"].ne("") & inv["seq_can"].str.lower().ne("nan")]

    # resolved-only can be empty; ignore empties for resolved matching
    inv_res = inv.dropna(subset=["seq_resolved_only"]).copy()
    inv_res = inv_res[inv_res["seq_resolved_only"].ne("") & inv_res["seq_resolved_only"].str.lower().ne("nan")]

    # ---- build lookups ----
    inv_full_grouped = (
        inv.groupby(["pdb_id", "seq_can"])["chain_id"]
        .agg(lambda s: sorted(set(s.tolist())))
        .reset_index()
        .rename(columns={"seq_can": "query_seq", "chain_id": "matching_chains_full"})
    )

    inv_res_grouped = (
        inv_res.groupby(["pdb_id", "seq_resolved_only"])["chain_id"]
        .agg(lambda s: sorted(set(s.tolist())))
        .reset_index()
        .rename(columns={"seq_resolved_only": "query_seq", "chain_id": "matching_chains_resolved"})
    )

    def _pick_match(full_list, res_list):
        if isinstance(full_list, list) and len(full_list) > 0:
            return full_list, "full"
        if isinstance(res_list, list) and len(res_list) > 0:
            return res_list, "resolved"
        return np.nan, "none"

    # ---- receptor matching: full then resolved ----
    df = df.merge(
        inv_full_grouped.rename(columns={"query_seq": "Receptor Sequence"}),
        on=["pdb_id", "Receptor Sequence"],
        how="left",
    )
    df = df.merge(
        inv_res_grouped.rename(columns={"query_seq": "Receptor Sequence"}),
        on=["pdb_id", "Receptor Sequence"],
        how="left",
    )

    rec_pick = df.apply(
        lambda r: _pick_match(r.get("matching_chains_full", np.nan),
                              r.get("matching_chains_resolved", np.nan)),
        axis=1,
        result_type="expand",
    )
    df["receptor_matching_chains"] = rec_pick[0]
    df["receptor_match_type"] = rec_pick[1]

    # ---- binder matching: full then resolved ----
    df = df.merge(
        inv_full_grouped.rename(columns={
            "query_seq": "Binder Sequence",
            "matching_chains_full": "binder_matching_chains_full"
        }),
        on=["pdb_id", "Binder Sequence"],
        how="left",
    )
    df = df.merge(
        inv_res_grouped.rename(columns={
            "query_seq": "Binder Sequence",
            "matching_chains_resolved": "binder_matching_chains_resolved"
        }),
        on=["pdb_id", "Binder Sequence"],
        how="left",
    )

    bind_pick = df.apply(
        lambda r: _pick_match(r.get("binder_matching_chains_full", np.nan),
                              r.get("binder_matching_chains_resolved", np.nan)),
        axis=1,
        result_type="expand",
    )
    df["binder_matching_chains"] = bind_pick[0]
    df["binder_match_type"] = bind_pick[1]

    # ---- found matches? ----
    df["has_receptor_match"] = df["receptor_match_type"].ne("none")
    df["has_binder_match"]   = df["binder_match_type"].ne("none")
    df["interaction_match_type"] = df["receptor_match_type"] + "/" + df["binder_match_type"]

    # ---- chain reporting (optional) ----
    df["Receptor Chain Correct"] = df["receptor_matching_chains"].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else np.nan
    )
    df["Peptide Chain Correct"] = df["binder_matching_chains"].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else np.nan
    )

    df["Receptor Chain Provided"] = df.get("Receptor Chain", np.nan)
    df["Peptide Chain Provided"]  = df.get("Peptide Chain", np.nan)

    # ---- mismatch logic: membership against all matching chains ----
    prov_rec = df["Receptor Chain Provided"].apply(_parse_chain_set)
    prov_pep = df["Peptide Chain Provided"].apply(_parse_chain_set)

    match_rec = df["receptor_matching_chains"].apply(_norm_chain_list)
    match_pep = df["binder_matching_chains"].apply(_norm_chain_list)

    df["receptor_chain_mismatch"] = (
        df["has_receptor_match"]
        & prov_rec.apply(len).gt(0)
        & match_rec.apply(len).gt(0)
        & (~(prov_rec & match_rec).apply(bool))
    )

    df["peptide_chain_mismatch"] = (
        df["has_binder_match"]
        & prov_pep.apply(len).gt(0)
        & match_pep.apply(len).gt(0)
        & (~(prov_pep & match_pep).apply(bool))
    )

    # ---- rejection reasons ----
    df["reject_reason"] = "accepted"
    df.loc[~df["has_receptor_match"], "reject_reason"] = "no_receptor_match_full_or_resolved"
    df.loc[df["has_receptor_match"] & ~df["has_binder_match"], "reject_reason"] = "no_binder_match_full_or_resolved"

    reject_on_chain_mismatch = False
    if reject_on_chain_mismatch:
        df.loc[df["has_receptor_match"] & df["has_binder_match"] & df["receptor_chain_mismatch"],
               "reject_reason"] = "receptor_chain_mismatch"
        df.loc[df["has_receptor_match"] & df["has_binder_match"] & df["peptide_chain_mismatch"],
               "reject_reason"] = "peptide_chain_mismatch"

        both_mis = (
            df["has_receptor_match"] & df["has_binder_match"]
            & df["receptor_chain_mismatch"] & df["peptide_chain_mismatch"]
        )
        df.loc[both_mis, "reject_reason"] = "both_chain_mismatch"
        
    list_cols = ["matching_chains_full", "matching_chains_resolved",
       "receptor_matching_chains",
       "binder_matching_chains_full", "binder_matching_chains_resolved",
       "binder_matching_chains"]
    
    for x in list_cols:
        df[x] = df[x].apply(lambda x: ",".join(x) if (type(x)!=float and type(x)!=str) else x)

    kept_df = df[df["reject_reason"].eq("accepted")].copy()
    rejected_df = df[~df["reject_reason"].eq("accepted")].copy()

    reason_counts = (
        rejected_df["reject_reason"]
        .value_counts(dropna=False)
        .rename_axis("reject_reason")
        .reset_index(name="n")
    )

    return kept_df, rejected_df, reason_counts

def _one_letter_from_comp_id(comp_id: str) -> str:
    """3-letter (or modified) residue name -> 1-letter; fall back to X."""
    try:
        r = gemmi.find_tabulated_residue(str(comp_id).upper())
        ol = r.one_letter_code
        return ol if ol and ol != "?" else "X"
    except Exception:
        return "X"

def get_resolved_positions_from_atom_site(cif_path: Path, chain_id: str) -> set[int]:
    doc = gemmi.cif.read_file(str(cif_path))
    block = doc.sole_block()
    atom = block.find_mmcif_category("_atom_site")
    if atom is None:
        return set()

    col_group = atom.find_column("group_PDB")
    col_asym  = atom.find_column("label_asym_id")
    col_seq   = atom.find_column("label_seq_id")

    if col_asym is None or col_seq is None:
        return set()

    resolved = set()
    n = len(atom) if hasattr(atom, "__len__") else atom.get_row_count()
    for i in range(n):   # <-- FIX: was atom.length()
        if col_group is not None:
            g = str(col_group[i])
            if g != "ATOM":
                continue

        if str(col_asym[i]).strip() != str(chain_id).strip():
            continue

        s = str(col_seq[i]).strip()
        if not s or s in {".", "?"}:
            continue
        try:
            resolved.add(int(s))
        except ValueError:
            pass

    return resolved

def build_resolved_sequences(full_seq: str, resolved_positions: set[int]) -> tuple[str, str, int, float]:
    """
    full_seq is the canonical polymer sequence (length L).
    resolved_positions are 1..L indices that have coordinates.
    """
    L = len(full_seq)
    masked = []
    only = []
    for i in range(1, L + 1):
        aa = full_seq[i - 1]
        if i in resolved_positions:
            masked.append(aa)
            only.append(aa)
        else:
            masked.append("-")
    n_res = len(resolved_positions)
    frac = (n_res / L) if L else 0.0
    return "".join(masked), "".join(only), n_res, frac

def convert_to_full_sequences(matched: pd.DataFrame, inv: pd.DataFrame):
    """
    In the original dataset, sometimes the sequence kept was just the resolved one.
    Sometimes it was the full sequence. 
    Here we harmonize the data so that we only keep the *full sequences* for both binder and receptor.
    """
    matched["PDB"] = matched["PDB"].str.upper()
    inv["PDB"] = inv["pdb_id"].str.upper()

    # for now, convert list columns (like ['A','B',F']) to ","-concatenated like "A,B,F"
    list_cols = ['matching_chains_full', 'matching_chains_resolved',
       'receptor_matching_chains',
       'binder_matching_chains_full', 'binder_matching_chains_resolved',
       'binder_matching_chains']
    for x in list_cols:
        matched[x] = matched[x].apply(lambda x: ast.literal_eval(x) if (type(x)!=float and x.startswith("[")) else x)
        matched[x] = matched[x].apply(lambda x: ",".join(x) if (type(x)!=float and type(x)!=str) else x)
        
    ## BEFORE harmonizing to full
    # Build up to a group on seq_sorts where every row is a different seq_sort
    agg_spec = {
        "Binder Sequence": ("Binder Sequence", "first"),
        "Receptor Sequence": ("Receptor Sequence", "first"),
        "Binder Length": ("Binder Length", "first"),
        "Receptor Length": ("Receptor Length", "first"),
        "Receptor Chain": ("Receptor Chain", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain": ("Peptide Chain", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "matching_chains_full": ("matching_chains_full", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "matching_chains_resolved": ("matching_chains_resolved", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "receptor_matching_chains": ("receptor_matching_chains", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "receptor_match_type": ("receptor_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains_full" : ("binder_matching_chains_full", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains_resolved": ("binder_matching_chains_resolved", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains": ("binder_matching_chains", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_match_type": ("binder_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "has_receptor_match": ("has_receptor_match", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "has_binder_match": ("has_binder_match", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "interaction_match_type": ("interaction_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Receptor Chain Correct": ("Receptor Chain Correct", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain Correct": ("Peptide Chain Correct", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Receptor Chain Provided": ("Receptor Chain Provided", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain Provided": ("Peptide Chain Provided", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "reject_reason": ("reject_reason", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
    }
    gb = matched.groupby(["seq_sort", "PDB"]).agg(**agg_spec).reset_index()
    gb = harmonize_nulls_to_nan(gb)
    gb["has_receptor_match"] = gb["has_receptor_match"].astype(bool)
    gb["has_binder_match"] = gb["has_binder_match"].astype(bool)
    
    agg_spec = {
        "Binder Sequence": ("Binder Sequence", "first"),
        "Receptor Sequence": ("Receptor Sequence", "first"),
        "Binder Length": ("Binder Length", "first"),
        "Receptor Length": ("Receptor Length", "first"),
        "PDB": ("PDB", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain": ("Receptor Chain", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain": ("Binder Chain", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "matching_chains_full": ("matching_chains_full", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "matching_chains_resolved": ("matching_chains_resolved", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "receptor_matching_chains": ("receptor_matching_chains", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "receptor_match_type": ("receptor_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains_full" : ("binder_matching_chains_full", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains_resolved": ("binder_matching_chains_resolved", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains": ("binder_matching_chains", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_match_type": ("binder_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "has_receptor_match": ("has_receptor_match", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "has_binder_match": ("has_binder_match", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "interaction_match_type": ("interaction_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain Correct": ("Receptor Chain Correct", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain Correct": ("Binder Chain Correct", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain Provided": ("Receptor Chain Provided", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain Provided": ("Binder Chain Provided", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "reject_reason": ("reject_reason", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
    }
    gb_unique_pairs = gb.groupby(["seq_sort"]).agg(**agg_spec).reset_index()
    gb_unique_pairs = harmonize_nulls_to_nan(gb_unique_pairs)
    gb_unique_pairs["has_receptor_match"] = gb_unique_pairs["has_receptor_match"].astype(bool)
    gb_unique_pairs["has_binder_match"] = gb_unique_pairs["has_binder_match"].astype(bool)

    test1 = len(gb_unique_pairs.loc[
        (gb_unique_pairs["receptor_matching_chains"].isna()) | 
        (gb_unique_pairs["binder_matching_chains"].isna())
        ])==0
    logger.info(f"All sequence pairs match at least one receptor chain and binder chain in one PDB file: {test1}")
    test1 = len(gb.loc[
        (gb["receptor_matching_chains"].isna()) | 
        (gb["binder_matching_chains"].isna())
        ])==0
    logger.info(f"All sequence+PDB pairs match at least one receptor chain and binder chain in one PDB file: {test1}")

    ### MERGE with inv to get full sequences only
    gb["receptor_matching_chains"] = gb["receptor_matching_chains"].apply(lambda x: x.split(",") if type(x)==str else x)
    gb["binder_matching_chains"] = gb["binder_matching_chains"].apply(lambda x: x.split(",") if type(x)==str else x)
    gb_expl = gb.explode("receptor_matching_chains").reset_index(drop=True)
    gb_expl["PDB"] = gb_expl["PDB"].str.upper()
    gb_expl = gb_expl.explode("binder_matching_chains").reset_index(drop=True)
    gb_expl = pd.merge(
        gb_expl,
        inv[["PDB","chain_id","seq_can"]].rename(columns={
            "chain_id":"receptor_matching_chains",
            "seq_can": "receptor_seq_full"}),
        on=["PDB","receptor_matching_chains"],
        how="left"
    )
    gb_expl = pd.merge(
        gb_expl,
        inv[["PDB","chain_id","seq_can"]].rename(columns={
            "chain_id":"binder_matching_chains",
            "seq_can": "binder_seq_full"}),
        on=["PDB","binder_matching_chains"],
        how="left"
    ).reset_index(drop=True)
    gb_expl["receptor_seq_full_length"] = gb_expl["receptor_seq_full"].apply(lambda x: len(x))
    gb_expl["binder_seq_full_length"] = gb_expl["binder_seq_full"].apply(lambda x: len(x))
    gb_expl["seq_sort_full"] = gb_expl.apply(lambda row: get_unique_id(row, colA="receptor_seq_full",colB="binder_seq_full"),axis=1)

    agg_spec = {
        "Binder Sequence PepMLM": ("Binder Sequence", "first"),
        "Receptor Sequence PepMLM": ("Receptor Sequence", "first"),
        "Binder Sequence Full": ("binder_seq_full", "first"),
        "Receptor Sequence Full": ("receptor_seq_full","first"),
        "Binder Length PepMLM": ("Binder Length", "first"),
        "Receptor Length PepMLM": ("Receptor Length", "first"),
        "Binder Length Full": ("binder_seq_full_length", "first"),
        "Receptor Length Full": ("receptor_seq_full_length", "first"),
        "Receptor Chain": ("Receptor Chain", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain": ("Binder Chain", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "matching_chains_full": ("matching_chains_full", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "matching_chains_resolved": ("matching_chains_resolved", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "receptor_matching_chains": ("receptor_matching_chains", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "receptor_match_type": ("receptor_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains_full" : ("binder_matching_chains_full", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains_resolved": ("binder_matching_chains_resolved", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_matching_chains": ("binder_matching_chains", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "binder_match_type": ("binder_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "has_receptor_match": ("has_receptor_match", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "has_binder_match": ("has_binder_match", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "interaction_match_type": ("interaction_match_type", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Receptor Chain Correct": ("Receptor Chain Correct", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain Correct": ("Binder Chain Correct", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Receptor Chain Provided": ("Receptor Chain Provided", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "Binder Chain Provided": ("Binder Chain Provided", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
        "reject_reason": ("reject_reason", lambda x: ",".join(sorted([str(y) for y in set(x) if type(x)!=float]))),
    }
    gb_full = gb_expl.groupby(["seq_sort_full", "PDB"]).agg(**agg_spec).reset_index()
    gb_full = harmonize_nulls_to_nan(gb_full)
    gb_full["has_receptor_match"] = gb_full["has_receptor_match"].astype(bool)
    gb_full["has_binder_match"] = gb_full["has_binder_match"].astype(bool)
    
    agg_spec = {
        "Binder Sequence PepMLM": ("Binder Sequence PepMLM", "first"),
        "Receptor Sequence PepMLM": ("Receptor Sequence PepMLM", "first"),
        "Binder Sequence Full": ("Binder Sequence Full", "first"),
        "Receptor Sequence Full": ("Receptor Sequence Full","first"),
        "Binder Length PepMLM": ("Binder Length PepMLM", "first"),
        "Receptor Length PepMLM": ("Receptor Length PepMLM", "first"),
        "Binder Length Full": ("Binder Length Full", "first"),
        "Receptor Length Full": ("Receptor Length Full", "first"),
        "PDB": ("PDB", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain": ("Receptor Chain", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain": ("Binder Chain", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "matching_chains_full": ("matching_chains_full", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "matching_chains_resolved": ("matching_chains_resolved", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "receptor_matching_chains": ("receptor_matching_chains", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "receptor_match_type": ("receptor_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains_full" : ("binder_matching_chains_full", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains_resolved": ("binder_matching_chains_resolved", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_matching_chains": ("binder_matching_chains", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "binder_match_type": ("binder_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "has_receptor_match": ("has_receptor_match", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "has_binder_match": ("has_binder_match", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "interaction_match_type": ("interaction_match_type", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain Correct": ("Receptor Chain Correct", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain Correct": ("Binder Chain Correct", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Receptor Chain Provided": ("Receptor Chain Provided", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "Binder Chain Provided": ("Binder Chain Provided", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
        "reject_reason": ("reject_reason", lambda x: "|".join([str(y) if type(y)==str else "" for y in x])),
    }
    gb_full_unique_pairs = gb_full.groupby(["seq_sort_full"]).agg(**agg_spec).reset_index()
    gb_full_unique_pairs = harmonize_nulls_to_nan(gb_full_unique_pairs)
    gb_full_unique_pairs["has_receptor_match"] = gb_full_unique_pairs["has_receptor_match"].astype(bool)
    gb_full_unique_pairs["has_binder_match"] = gb_full_unique_pairs["has_binder_match"].astype(bool)

    return gb_unique_pairs, gb_full_unique_pairs

def main(cfg: DictConfig):
    # --- Load inputs ---
    complex_path = Path(root) / cfg.process.complex_file
    pepnn_test_path = Path(root) / cfg.process.pepnn_test
    pepnn_train_path = Path(root) / cfg.process.pepnn_train
    propedia_test_path = Path(root) / cfg.process.propedia_test

    complexes = pd.read_csv(complex_path)
    pepnn_test = pd.read_csv(pepnn_test_path)
    pepnn_train = pd.read_csv(pepnn_train_path)
    propedia_test = pd.read_csv(propedia_test_path)

    combined_df = initial_filtering(
        complexes=complexes,
        pepnn_test=pepnn_test,
        pepnn_train=pepnn_train,
        propedia_test=propedia_test,
    )

    # Normalize PDB ids
    combined_df["pdb_id"] = combined_df["PDB"].astype(str).str.strip().str.lower()

    # --- Parallel download (one row per PDB id) ---
    out_dir = Path(root) / cfg.process.struct_out_dir          # e.g. "raw_data/rcsb_cif"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    redownload = bool(getattr(cfg.process, "redownload", False))

    unique_pdb_ids = sorted(combined_df["pdb_id"].unique().tolist())
    logger.info(f"Total unique PDB IDs: {len(unique_pdb_ids)}")
    to_download = filter_pdbs_to_download(unique_pdb_ids, out_dir, redownload=redownload)
    logger.info(f"Total to download (have not already been downloaded): {len(to_download)}")
    
    # worker directory is going to be the Hydra directory
    hydra_dir = HydraConfig.get().run.dir

    download_log = download_structures_parallel_chunked(
        to_download,
        output_dir=str(out_dir),
        primary_format="cif",
        max_workers=int(cfg.process.max_workers),
        log_csv=str(Path(root) / cfg.process.download_log_csv),  # optional
        worker_log_dir=str(Path(hydra_dir) / "download_worker_logs"), 
    )
    
    # 1) download_log from download_structures_parallel(...)
    # 2) build chain inventory
    inventory_path = Path(root) / cfg.process.chain_inventory_csv
    inventory_df = build_full_chain_inventory(
        download_log=download_log,
        inventory_out=inventory_path,
        output_dir=out_dir,
        max_workers=int(cfg.process.max_workers),
        append=True,  # this is what you want
    )
    logger.info(f"Saved inventory_df to: {inventory_path} ")

    # 3) match your interaction DB back to inventory
    kept_df, rejected_df, reason_counts = match_interactions_to_inventory(combined_df, inventory_df)

    matched_out = Path(root) / cfg.process.matched_interactions_csv
    matched_out.parent.mkdir(parents=True, exist_ok=True)
    kept_df.to_csv(matched_out, index=False)
    logger.info(f"Saved kept_df to: {matched_out} (n={len(kept_df)})")

    # put rejects next to matched file
    reject_out = matched_out.with_suffix(".rejected.csv")
    rejected_df.to_csv(reject_out, index=False)
    logger.info(f"Saved rejected_df to: {reject_out} (n={len(rejected_df)})")

    counts_out = matched_out.with_suffix(".rejected_reason_counts.csv")
    reason_counts.to_csv(counts_out, index=False)
    logger.info(f"Saved reject reason counts to: {counts_out}")
    logger.info(f"Reject reason summary:\n{reason_counts.to_string(index=False)}")
    
    # do some post-processing to recombine everything with matched sequences
    unique_pairs, full_unique_pairs = convert_to_full_sequences(matched=kept_df, inv=inventory_df)
    unique_pairs_out = Path(root) / cfg.process.pepmlm_pairs_csv
    unique_pairs_out.parent.mkdir(parents=True, exist_ok=True)
    unique_pairs.to_csv(unique_pairs_out, index=False)
    logger.info(f"Saved unique_pairs to: {unique_pairs_out} (n={len(unique_pairs)})")
    
    full_unique_pairs_out = Path(root) / cfg.process.full_pairs_csv
    full_unique_pairs_out.parent.mkdir(parents=True, exist_ok=True)
    full_unique_pairs.to_csv(full_unique_pairs_out, index=False)
    logger.info(f"Saved full_unique_pairs to: {full_unique_pairs_out} (n={len(full_unique_pairs)})")

    
if __name__ == "__main__":
    main()