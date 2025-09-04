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
    Download BioLiP2 files
    """
    
    # Instantiate downloader object
    downloader = DownloadObject(
        logger=logger,
        name=cfg.download.name,
        script_path=cfg.script_path,
        delete_zip=cfg.download.delete_zip
    )
    
    # Download main BioLiP2 database
    downloader.download(
        url=cfg.download.url,
        output_dir=cfg.download.output_dir,
        filename=cfg.download.filename
    )
    
    # Download protein sequences
    downloader.download(
        url=cfg.download.prot_seqs_url,
        output_dir=cfg.download.prot_seqs_output_dir,
        filename=cfg.download.prot_seqs_filename
    )
    
    # Download peptide sequences
    downloader.download(
        url=cfg.download.pep_seqs_url,
        output_dir=cfg.download.pep_seqs_output_dir,
        filename=cfg.download.pep_seqs_filename
    )
    
if __name__ == "__main__":
    main()