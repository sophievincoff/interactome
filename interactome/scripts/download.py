import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

import hydra
from omegaconf import DictConfig
import logging

logger = logging.getLogger(__name__)

# import your processing entry points here
from interactome.data.download.intact import main as download_intact_main
from interactome.data.download.biogrid import main as download_biogrid_main
from interactome.data.download.hippie import main as download_hippie_main
from interactome.data.download.biolip2 import main as download_biolip2_main

@hydra.main(config_path=str(root / "configs"), config_name="download_main", version_base="1.3")
def main(cfg: DictConfig):
    db_name = cfg.download.name.lower()
    logger.info(f"Running download for database: {db_name}")
    
    if db_name == "intact":
        download_intact_main(cfg)
    elif db_name == "biogrid":
        download_biogrid_main(cfg)
    elif db_name == "hippie":
        download_hippie_main(cfg)
    elif db_name == "biolip2":
        download_biolip2_main(cfg)
    else:
        raise ValueError(f"No processor defined for database: {db_name}")

if __name__ == "__main__":
    main()
