from omegaconf import DictConfig
from pathlib import Path
import logging
import subprocess

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
        script_path=cfg.script_path
    )
    
    # Download main databse
    if False:
        downloader.download(
            url=cfg.download.url,
            output_dir=cfg.download.output_dir,
            filename=cfg.download.filename,
            unzip=cfg.download.unzip,
            delete_zip=cfg.download.delete_zip
        )
        
    # Download XML and don't unzip it
    downloader.download(
        url=cfg.download.url_xml,
        output_dir=cfg.download.output_dir,
        filename=cfg.download.filename_xml,
        unzip=cfg.download.unzip_xml,
        delete_zip=cfg.download.delete_zip_xml
    )
    
if __name__ == "__main__":
    main()