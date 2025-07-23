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
    biogrid_fpath = Path(root) / cfg.process.input_file # raw biogrid download
    example_output_dir = Path(root) / cfg.process.example_output_dir    # output directory for example
    os.makedirs(example_output_dir, exist_ok=True)
    
    example_subset = pd.read_csv(biogrid_fpath, sep='\t', header=0, low_memory=False)
    example_subset = example_subset.sample(n=1000, random_state=42).reset_index(drop=True) # Sample 1000 rows for the example
    example_subset.to_csv(example_output_dir / "biogrid_example.tsv", sep='\t', index=False)
    return

if __name__ == "__main__":
    main()