from omegaconf import DictConfig
from pathlib import Path
import logging
import subprocess

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

def main(cfg: DictConfig):
    """
    Download IntAct file, which is one zip file
    """
    logger = logging.getLogger(__name__)

    script_path = root / "data/download/download_unzip.sh"
   
    url = cfg.download.url
    mi_map_url = cfg.download.mi_map_url
    
    filename = cfg.download.filename
    mi_map_filename = cfg.download.mi_map_filename
    
    output_dir = root / cfg.download.output_dir

    ##### 
    #### Download full database
    ####
    logger.info(f"Running download for {cfg.download.name} full database")
    logger.info(f"Script: {script_path}")
    logger.info(f"URL: {url}")
    logger.info(f"Output: {output_dir / filename}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the download_unzip script as a subproces 
    result = subprocess.run(
        ["bash", str(script_path), url, str(output_dir), filename],
        capture_output=True,
        text=True
    )

    logger.info("STDOUT:\n" + result.stdout)
    if result.stderr:
        logger.warning("STDERR:\n" + result.stderr)

    if result.returncode != 0:
        logger.error(f"Download script exited with code {result.returncode}")
    else:
        logger.info("Download completed successfully.")
    
    ##### 
    #### PSI-MI 2.5 Mapping file
    ####
    logger.info(f"Running download for {cfg.download.name} MI Map file")
    logger.info(f"Script: {script_path}")
    logger.info(f"URL: {mi_map_url}")
    logger.info(f"Output: {output_dir / mi_map_filename}")
    
    # Run the download.sh script as a subproces 
    result = subprocess.run(
        ["bash", str(script_path), mi_map_url, str(output_dir), mi_map_filename],
        capture_output=True,
        text=True
    )

    logger.info("STDOUT:\n" + result.stdout)
    if result.stderr:
        logger.warning("STDERR:\n" + result.stderr)

    if result.returncode != 0:
        logger.error(f"Download script exited with code {result.returncode}")
    else:
        logger.info("Download completed successfully.")
    
if __name__ == "__main__":
    main()