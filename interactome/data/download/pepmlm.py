from omegaconf import DictConfig
from pathlib import Path
import logging
import subprocess
import os

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
        
    # Download the raw data in txt format
    downloader.download(
        url=cfg.download.url,
        output_dir=cfg.download.output_dir,
        filename=cfg.download.filename,
        unzip=cfg.download.unzip,
        delete_zip=cfg.download.delete_zip
    )
    
    created_dir = str(Path(root)/cfg.download.output_dir/"Train and Test Data")
    print(f"Created dir: {created_dir}")
    new_dir = created_dir.replace("Train and Test Data","train_test_data")
    os.rename(created_dir, new_dir)
    print(f"Moved to: {new_dir}")
    
if __name__ == "__main__":
    main()