import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

import hydra
from omegaconf import DictConfig
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# import your processing entry points here
from interactome.data.process_intact import main as process_intact_main
from interactome.data.process_biogrid import main as process_biogrid_main

@hydra.main(config_path=str(root / "configs"), config_name="process", version_base="1.3")
def main(cfg: DictConfig):
    db_name = cfg.process.name.lower()
    logger.info(f"Running processor for database: {db_name}")
    
    if db_name == "intact":
        process_intact_main(cfg)
    elif db_name == "biogrid":
        process_biogrid_main(cfg)
    else:
        raise ValueError(f"No processor defined for database: {db_name}")

if __name__ == "__main__":
    main()
