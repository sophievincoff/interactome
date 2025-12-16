import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

import hydra
from omegaconf import DictConfig
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# import your processing entry points here
from interactome.data.analyze.alphafolddb_query import main as alphafolddb_query_main

@hydra.main(config_path=str(root / "configs"), config_name="analyze_main", version_base="1.3")
def main(cfg: DictConfig):
    name = cfg.analyze.name.lower()
    logger.info(f"Running processor for database: {name}")
    
    if name == "alphafolddb_query":
        alphafolddb_query_main(cfg)
    else:
        raise ValueError(f"No processor defined for: {name}")

if __name__ == "__main__":
    main()
