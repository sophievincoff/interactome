from pathlib import Path
import subprocess
import os
import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

class DownloadObject:
    def __init__(self, logger, name: str, script_path: Path, delete_zip: bool=True):
        self.logger = logger
        self.name = name
        self.script_path = Path(root) / script_path
        self.delete_zip = str(delete_zip).lower()
        assert self.delete_zip in ("true", "false")

    def download(self, url: str, output_dir: str, filename: str):
        """
        Download file to the specified output dir and save it with filename. 
        """
        output_dir = Path(root) / output_dir
        output_path = output_dir / filename

        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info(f"Running download for {self.name}")
        self.logger.info(f"Script: {self.script_path}")
        self.logger.info(f"URL: {url}")
        self.logger.info(f"Output: {output_path}")

        result = subprocess.run(
            ["bash", str(self.script_path), url, str(output_dir), filename, self.delete_zip],
            capture_output=True,
            text=True
        )

        self.logger.info("STDOUT:\n" + result.stdout)
        if result.stderr:
            self.logger.warning("STDERR:\n" + result.stderr)

        if result.returncode != 0:
            self.logger.error(f"Download script exited with code {result.returncode}")
        else:
            self.logger.info("Download completed successfully.")
