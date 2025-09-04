from omegaconf import DictConfig
from pathlib import Path
import logging
import subprocess
import pandas as pd

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
    
    # Download main BioGRID database
    downloader.download(
        url=cfg.download.url,
        output_dir=cfg.download.output_dir,
        filename=cfg.download.filename
    )
    
    # Download MI terms map
    downloader.download(
        url=cfg.download.mi_map_url,
        output_dir=cfg.download.output_dir,
        filename=cfg.download.mi_map_filename.replace(".txt", ".xls")
    )
    
    path_to_mimap = Path(root) / cfg.download.mi_map_filename.replace(".txt", ".xls")
    save_path = Path(root) / cfg.download.mi_map_filename
    df = pd.read_excel(path_to_mimap, sheet_name="Mapping")  # You can also specify the sheet name
    df.to_csv(save_path, sep="\t", index=False)

    
if __name__ == "__main__":
    main()