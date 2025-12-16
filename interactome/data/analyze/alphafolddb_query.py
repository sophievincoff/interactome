#!/usr/bin/env python3
import concurrent.futures
import json
import time
from pathlib import Path
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
import os
import requests
import rootutils
import logging
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

BASE_URL = "https://alphafold.ebi.ac.uk/api/prediction"

def fetch_single_uniprot(uniprot_id: str, timeout: float = 20.0):
    """
    Fetch AlphaFoldDB prediction for a single UniProt ID.

    Returns a dict with the fields of interest (or None if not found / error).
    """
    uniprot_id = uniprot_id.strip()
    if not uniprot_id:
        return None

    url = f"{BASE_URL}/{uniprot_id}"
    try:
        resp = requests.get(url, headers={"accept": "application/json"}, timeout=timeout)
        if resp.status_code == 404:
            # No prediction available
            return [{
                "uniprot_id": uniprot_id,
                "found": False,
                "globalMetricValue": None,
                "fractionPlddtVeryLow": None,
                "fractionPlddtLow": None,
                "fractionPlddtConfident": None,
                "fractionPlddtVeryHigh": None,
                "sequence": None,
                "sequenceStart": None,
                "sequenceEnd": None,
                "uniprotAccession": None,
                "modelEntityId": None,
            }]
        resp.raise_for_status()
        data = resp.json()

        # API returns a list of dicts: typically take the first prediction
        if not isinstance(data, list) or len(data) == 0:
            return [{
                "uniprot_id": uniprot_id,
                "found": False,
                "globalMetricValue": None,
                "fractionPlddtVeryLow": None,
                "fractionPlddtLow": None,
                "fractionPlddtConfident": None,
                "fractionPlddtVeryHigh": None,
                "sequence": None,
                "sequenceStart": None,
                "sequenceEnd": None,
                "uniprotAccession": None,
                "modelEntityId": None,
            }]

        all_predictions = []
        for pred in data:
            # Safely extract fields (use .get so missing keys don't crash)
            all_predictions.append({
                "uniprot_id": uniprot_id,
                "found": True,
                "globalMetricValue": pred.get("globalMetricValue"),
                "fractionPlddtVeryLow": pred.get("fractionPlddtVeryLow"),
                "fractionPlddtLow": pred.get("fractionPlddtLow"),
                "fractionPlddtConfident": pred.get("fractionPlddtConfident"),
                "fractionPlddtVeryHigh": pred.get("fractionPlddtVeryHigh"),
                "sequence": pred.get("sequence"),
                "sequenceStart": pred.get("sequenceStart"),
                "sequenceEnd": pred.get("sequenceEnd"),
                "uniprotAccession": pred.get("uniprotAccession"),
                "modelEntityId": pred.get("modelEntityId"),
            })
        return all_predictions

    except requests.RequestException as e:
        # Network / HTTP error -> log as not found with error message if you like
        return [{
            "uniprot_id": uniprot_id,
            "found": False,
            "globalMetricValue": None,
            "fractionPlddtVeryLow": None,
            "fractionPlddtLow": None,
            "fractionPlddtConfident": None,
            "fractionPlddtVeryHigh": None,
            "sequence": None,
            "sequenceStart": None,
            "sequenceEnd": None,
            "uniprotAccession": None,
            "modelEntityId": None,
            "error": str(e),
        }]

def read_uniprot_ids(path: Path):
    if str(path).endswith(".txt"):
        with path.open() as f:
            ids = [line.strip() for line in f if line.strip()]
        # de-dup while preserving order
        seen = set()
        uniq_ids = []
        for uid in ids:
            if uid not in seen:
                seen.add(uid)
                uniq_ids.append(uid)
        return uniq_ids
    else:
        raise Exception("Reading this filetype not implemented.")
    return []

def write_tsv(results, out_path: Path):
    import csv

    all_keys = set()
    for r in results:
        if isinstance(r, dict):
            all_keys.update(r.keys())

    base_order = [
        "uniprot_id",
        "found",
        "globalMetricValue",
        "fractionPlddtVeryLow",
        "fractionPlddtLow",
        "fractionPlddtConfident",
        "fractionPlddtVeryHigh",
        "sequence",
        "sequenceStart",
        "sequenceEnd",
        "uniprotAccession",
        "modelEntityId",
        "error",
    ]
    fieldnames = [k for k in base_order if k in all_keys] + [
        k for k in sorted(all_keys) if k not in base_order
    ]

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            if isinstance(row, dict):
                writer.writerow(row)

def main(cfg: DictConfig):
    output_path = Path(root) / cfg.analyze.output_path
    os.makedirs(output_path.parent, exist_ok=True)
    if cfg.analyze.debug:
        debug_path = Path(root) / cfg.analyze.debug_file
        uniprot_ids = read_uniprot_ids(debug_path)
        logger.info(f"Loaded {len(uniprot_ids)} unique UniProt IDs from debug_file: {debug_path}")
        use_n = min(len(uniprot_ids), 16)
        uniprot_ids = uniprot_ids[0:use_n]
        logger.info(f"DEBUG MODE: {cfg.analyze.debug} -> trimmed to {use_n} sequences")
    else:
        input_path = Path(root) / cfg.analyze.input_path
        uniprot_ids = read_uniprot_ids(input_path)
        logger.info(f"Loaded {len(uniprot_ids)} unique UniProt IDs")

    results = []
    # Use ThreadPoolExecutor because this is network-bound I/O.
    with concurrent.futures.ThreadPoolExecutor(max_workers=cfg.analyze.workers) as executor:
        future_to_id = {}
        for uid in uniprot_ids:
            fut = executor.submit(fetch_single_uniprot, uid)
            future_to_id[fut] = uid
            if cfg.analyze.sleep > 0:
                time.sleep(cfg.analyze.sleep)

        for i, fut in enumerate(concurrent.futures.as_completed(future_to_id), 1):
            uid = future_to_id[fut]
            try:
                res = fut.result()
            except Exception as e:
                res = [{
                    "uniprot_id": uid,
                    "found": False,
                    "globalMetricValue": None,
                    "fractionPlddtVeryLow": None,
                    "fractionPlddtLow": None,
                    "fractionPlddtConfident": None,
                    "fractionPlddtVeryHigh": None,
                    "sequence": None,
                    "sequenceStart": None,
                    "sequenceEnd": None,
                    "uniprotAccession": None,
                    "modelEntityId": None,
                    "error": f"unhandled: {e}",
                }]
            if res is None:
                continue
            elif isinstance(res, dict):
                results.append(res)
            elif isinstance(res, list):
                results.extend(res)
            else:
                logger.warning(f"Unexpected result type for {uid}: {type(res)}")
            if i % 50 == 0:
                logger.info(f"Completed {i}/{len(uniprot_ids)} IDs")

    write_tsv(results, output_path)
    logger.info(f"Wrote {len(results)} rows to {output_path}")


if __name__ == "__main__":
    main()
