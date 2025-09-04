from omegaconf import DictConfig
import logging

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
logger = logging.getLogger(__name__)

from .utils import DownloadObject

def main(cfg: DictConfig):
    """
    Download IntAct file, which is one zip file
    """
    
    # Instantiate downloader object
    downloader = DownloadObject(
        logger=logger,
        name=cfg.download.name,
        script_path=cfg.script_path,
        delete_zip=cfg.download.delete_zip
    )
    
if __name__ == "__main__":
    main()