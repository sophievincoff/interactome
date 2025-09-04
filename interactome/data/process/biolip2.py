from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path
import os
import pandas as pd
import logging

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

def main(cfg: DictConfig):
    # Set up needed paths
    database_path = Path(root) / cfg.process.database.input_file # raw biogrid download
    protein_seqs_path = Path(root) / cfg.process.protein_seqs.input_file
    peptide_seqs_path = Path(root) / cfg.process.peptide_seqs.input_file
    example_output_dir = Path(root) / cfg.process.example_output_dir    # output directory for example
    os.makedirs(example_output_dir, exist_ok=True)
    
    data = pd.read_csv(database_path, sep='\t', header=0, low_memory=False)
    columns = [
        "rcsb", 
        "receptor_chain",
        "resolution", # -1.00 means lack of resolution information
        "binding_site_code",
        "ligand_id",
        "ligand_chain",
        "ligand_serial",
        "binding_site_residues_pdbnum",
        "binding_site_residues_1index",
        "catalytic_site_residues_pdbnum", # different sites separated by ;
        "catalytic_site_residues_1index", # different sites separated by ;
        "ec",
        "go",
        "affinity_litsurvey",
        "affinity_BindingMOAD",
        "affinity_PDBbindCN",
        "affinity_BindingDB",
        "uniprotkb",
        "pubmed",
        "ligand_seq_no",
        "receptor_sequence" 
    ]
    data.columns = columns
    example_subset = data.sample(n=1000, random_state=42).reset_index(drop=True) # Sample 1000 rows for the example
    example_subset.to_csv(example_output_dir / "biolip2_example.tsv", sep='\t', index=False)
    
    logging.info(data["ligand_id"].value_counts())
    return

if __name__ == "__main__":
    main()