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
from interactome.data.download.mentha import main as download_mentha_main
from interactome.data.download.negatome2 import main as download_negatome2_main
from interactome.data.download.uniref50 import main as download_uniref50_main
from interactome.data.download.uniref90 import main as download_uniref90_main
from interactome.data.download.gene2pubmed import main as download_gene2pubmed_main
from interactome.data.download.dibs import main as download_dibs_main
from interactome.data.download.mfib2 import main as download_mfib2_main
from interactome.data.download.pepmlm import main as download_pepmlm_main

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
    elif db_name == "mentha":
        download_mentha_main(cfg)
    elif db_name == "negatome2":
        download_negatome2_main(cfg)
    elif db_name == "uniref50":
        download_uniref50_main(cfg)
    elif db_name == "uniref90":
        download_uniref90_main(cfg)
    elif db_name == "gene2pubmed":
        download_gene2pubmed_main(cfg)
    elif db_name == "dibs":
        download_dibs_main(cfg)
    elif db_name == "mfib2":
        download_mfib2_main(cfg)
    elif db_name == "pepmlm":
        download_pepmlm_main(cfg)
    else:
        raise ValueError(f"No processor defined for database: {db_name}")

if __name__ == "__main__":
    main()
