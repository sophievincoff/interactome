import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

import hydra
from omegaconf import DictConfig
from pathlib import Path
import subprocess
import logging

@hydra.main(config_path=str(root / "configs"), config_name="download", version_base="1.3")
def main(cfg: DictConfig):
    """
    Download a database of your choice to a location of your choice.
    """
    log = logging.getLogger(__name__)

    script_path = root / cfg.database.script_path
    url = cfg.database.url
    output_dir = root / cfg.database.output_dir
    filename = cfg.database.filename

    log.info(f"Running download for {cfg.database.name}")
    log.info(f"Script: {script_path}")
    log.info(f"URL: {url}")
    log.info(f"Output: {output_dir / filename}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the download.sh script as a subproces 
    result = subprocess.run(
        ["bash", str(script_path), url, str(output_dir), filename],
        capture_output=True,
        text=True
    )

    log.info("STDOUT:\n" + result.stdout)
    if result.stderr:
        log.warning("STDERR:\n" + result.stderr)

    if result.returncode != 0:
        log.error(f"Download script exited with code {result.returncode}")
    else:
        log.info("Download completed successfully.")


if __name__ == "__main__":
    main()
