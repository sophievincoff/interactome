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
from concurrent.futures import ThreadPoolExecutor

import requests_cache
requests_cache.install_cache("rcsb_cache", expire_after=24*3600)  # 24h cache

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

@lru_cache(maxsize=8192)
def get_chain_sequences(pdb_id, chain1, chain2):
    rcsbstructure = RCSBStructure(pdb_id, download_struct=False)
    chaindf = rcsbstructure.get_chain_info_df()
    chaindf["Chain"] = chaindf["Chain"].apply(lambda x: x.split(","))
    chaindf["Auth"] = chaindf["Auth"].apply(lambda x: x.split(","))
    chaindf = chaindf.explode(["Chain", "Auth"])
    d = {}
    d["entry_id"] = pdb_id
    d["asym_id_1"] = chain1
    d["asym_id_2"] = chain2
    chain1_info = chaindf.loc[chaindf["Chain"] == chain1]
    chain2_info = chaindf.loc[chaindf["Chain"] == chain2]
    d["uniprot_1"] = chain1_info["UniProt IDs"].item()
    d["uniprot_2"] = chain2_info["UniProt IDs"].item()
    d["pdb_seq_1"] = chain1_info["PDB Seq"].item()
    d["pdb_seq_2"] = chain2_info["PDB Seq"].item()
    return d

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
    neg_subset = negatives.sample(1000, random_state=42)
    start = time.time()
    d = neg_subset.apply(
        lambda row: get_chain_sequences(
            row["entry_id"], row["asym_id_1"], row["asym_id_2"]
        ),
        axis=1,
    )
    end = time.time()
    logger.info(f"Total time to gather RCSB PDB data for all negatives: {end-start:.2f}s")

    # save df
    df = pd.DataFrame(d.tolist())
    logger.info(f"Length of negative DF from unique entry-chain1-chain2 combinations: {len(df)}")
    logger.info(df)
    df.to_csv(Path(cfg.process.output_dir) / "rcsb_info_pinder_no_prodigy_label.csv", index=False)
    
    pos_subset = positives.sample(2, random_state=42)
    start = time.time()
    d = pos_subset.apply(
        lambda row: get_chain_sequences(
            row["entry_id"], row["asym_id_1"], row["asym_id_2"]
        ),
        axis=1,
    )
    end = time.time()
    logger.info(f"Total time to gather RCSB PDB data for all positives: {end-start:.2f}s")
    # save df
    df = pd.DataFrame(d.tolist())
    logger.info(f"Length of positive DF from unique entry-chain1-chain2 combinations: {len(df)}")
    logger.info(df)
    df.to_csv(Path(cfg.process.output_dir) / "rcsb_info_pinder_bioxtal_prodigy_label.csv", index=False)
