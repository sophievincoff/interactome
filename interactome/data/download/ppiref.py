# Python script for collecting raw data needed for processing
"""
This is the original script I used in fuson_flow to download the files I would need to process PPIRef
It has not been edited at all other than having fuson_flow imports fixed and removing log_update
TODO: integrate this with the hydra framework of the rest of the repository. 
"""
from ppiref.utils.misc import download_from_zenodo
from ppiref.split import read_split
from ppiref.utils.ppi import PPI
from protparser.rcsb import download_rcsb
from pathlib import Path
from multiprocessing import Manager, Pool
from functools import partial
import os
import shutil
import json
import requests

def flatten_directory(parent_dir: str):
    """
    Move all files from subdirectories into the parent directory 
    and delete the now-empty subdirectories, only if subdirectories exist.

    Args:
        parent_dir (str): Path to the parent directory.
    """
    # Check if there are any subdirectories
    subdirs = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    
    if not subdirs:
        print("No subdirectories found. Nothing to flatten.")
        return
    
    print(f"Found {len(subdirs)} subdirectories. Flattening directory: {parent_dir}")
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(parent_dir, topdown=False):
        for file in files:
            # Construct the full path to the file
            src = os.path.join(root, file)
            dest = os.path.join(parent_dir, file)
            
            # Ensure no filename collision
            if os.path.exists(dest):
                base, extension = os.path.splitext(file)
                dest = os.path.join(parent_dir, f"{base}_copy{extension}")
            
            # Move the file
            shutil.move(src, dest)
            print(f"Moved: {src} → {dest}")
        
        # Remove the subdirectory if it's now empty
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            try:
                os.rmdir(dir_path)
                print(f"Deleted empty directory: {dir_path}")
            except OSError:
                print(f"Failed to delete non-empty directory: {dir_path}")

def download_single_pdb(pdb_id, struct_format, output_dir, log_list):
    """
    Downloads a single PDB file and logs the progress.
    """
    other_struct_format="cif"
    if struct_format=="cif": other_struct_format="pdb"
    try:
        download_rcsb(pdb_id, struct_format=struct_format, convert_if_fail=False, output_dir=output_dir)
    except Exception as e:
        try:
            print(f"Failed to download {pdb_id} with {struct_format}. Trying with {other_struct_format}.")
            download_rcsb(pdb_id, struct_format=other_struct_format, convert_if_fail=False, output_dir=output_dir)
        except Exception as e:
            print(f"Failed to run download_rcsb for {pdb_id}: {e}")

def download_original_structures(pdb_ids, struct_format="pdb", output_dir="raw_data/ppi_6A_RCSB_full", num_processes=64):
    """
    Parallelized downloading of PDB files using multiprocessing.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nDownloading {len(pdb_ids)} PDB ids from the RCSB website.")
    
    # Check which ones are already there
    base_path = Path(output_dir)
    pdbs_already_downloaded = [str(x.stem) for x in (base_path.rglob(f"*.pdb"))]
    cifs_already_downloaded = [str(x.stem) for x in (base_path.rglob(f"*.cif"))]
    already_downloaded = set(pdbs_already_downloaded) | set(cifs_already_downloaded)
    remaining_pdbs = list(set(pdb_ids) - already_downloaded) 
    
    print(f"\tAlready downloaded: {len(already_downloaded)} ({len(pdbs_already_downloaded)} PDBs, {len(cifs_already_downloaded)} CIFs). Remaining: {len(remaining_pdbs)}")
    
    if not remaining_pdbs:
        print("All PDBs are already downloaded.")
        return
    
    # Use Manager for shared logging
    with Manager() as manager:
        log_list = manager.list()
        
        # Parallel downloading using Pool
        with Pool(processes=num_processes) as pool:
            pool.map(
                partial(download_single_pdb, struct_format=struct_format, output_dir=output_dir, log_list=log_list),
                remaining_pdbs
            )
        
        # print the logs after multiprocessing
        for log in log_list:
            print(log)
            
def download_file(url, save_path):
    """
    Download a file from a given URL and save it locally.

    Args:
        url (str): URL to the file.
        save_path (str): Local path to save the downloaded file.
    
    Returns:
        None
    """
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, 'wb') as file:
            file.write(response.content)
        print(f"File downloaded successfully and saved to {save_path}")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")
            
def download_ppiref_split_json(splitname='ppiref_6A_filtered'):
    """
    Download splits. We will be interested in ppiref_6A_filtered (PPIRef300K) and ppiref_6A_filtered_clustered_04 (PPIRef50K)
    """
    file_url = f"https://raw.githubusercontent.com/anton-bushuiev/PPIRef/main/ppiref/data/splits/{splitname}.json"
    save_location = f"raw_data/{splitname}.json"

    # Download the file (if we haven't yet)
    if not(os.path.exists(save_location)):
        download_file(file_url, save_location)
    
    # Read the download 
    with open(save_location, 'r') as f:
        data = json.load(f)['folds']['whole']

    return data

def download_ppiref_data(remove_noppi=True):
    """
    Download the PDBs in ppi_6A from PPIRef, using the PPIRef package from GitHub.
    Delete .noppi files
    """
    os.makedirs("raw_data",exist_ok=True)
    move_path = "raw_data/ppi_6A"
    
    if not(os.path.exists(move_path)):
        # download 6A dataset from PPIRef
        download_folder = "/usr/local/lib/python3.10/dist-packages/ppiref/data/ppiref"
        download_path = f"{download_folder}/ppi_6A" # the zenodo download unzips it
        
        if not(os.path.exists(download_path)):
            print(f"Downloading ppi_6A.zip to {download_folder}...")
            download_from_zenodo('ppi_6A.zip')
            print(f"\tDone.")
        else:
            print(f"Already downloaded ppi_6A.zip to {download_folder}. Continuing")
        
        # move all the data into raw data
        shutil.copytree(download_path,move_path)
        print(f"Copied PPIRef data from {download_folder} to {move_path}")
        
        # now remove the downloaded data to save space
        shutil.rmtree(download_path)
        print(f"Deleted original downloads at {download_path}")
        
        # remove .noppi files
        # also remove empty folders
        if remove_noppi:
            print("\nRemoving .noppi files and any folders with ONLY .noppi files...")
            for root, dirs, files in os.walk(move_path):
                for file in files:
                    if file.endswith(".noppi"):
                        file_path = os.path.join(root, file)
                        os.remove(file_path)
                        #print(f"Removed {file_path}")
            for root, dirs, files in os.walk(move_path):
                for dir in dirs:
                    if len(os.listdir(os.path.join(root,dir))) == 0:
                        os.rmdir(os.path.join(root,dir))
                        #print(f"Removed empty folder {os.path.join(root,dir)}")
    
    else:
        print(f"Already downloaded and deposited data in {move_path}. Delete this folder and rerun the script to redownload.")
        
def filter_to_splits(keep_ids):
    """
    Goes through the raw_data/ppi_6A folder and keeps only the PDBs that are in the keep_ids set
    Should be based on splits downloaded with download_ppiref_split_json
    """
    move_path = "raw_data/ppi_6A"
    
    # Make a Path object so we can easily glob pdb files, count how many we have
    base_path = Path(move_path)
    total_pdb_files = len(list(base_path.rglob("*.pdb")))
    print(f"\nTotal .pdb files in raw_data/ppi_6A: {total_pdb_files}")
    if total_pdb_files<700000:
        print(f"\tSmall size of raw_data/ppi_6A indicates this filtering has already been done. Repeating:")
    
    # Filter to keep only PDBs we want
    print(f"Filtering to only {len(keep_ids)} PDBs...")
    # Delete PDB files
    for root, dirs, files in os.walk(move_path):
        for file in files:
            if file.split(".pdb")[0] not in keep_ids:
                file_path = os.path.join(root, file)
                os.remove(file_path)
                #print(f"Removed {file_path}")
    # Delete empty directories
    for root, dirs, files in os.walk(move_path):
        for dir in dirs:
            if len(os.listdir(os.path.join(root,dir))) == 0:
                os.rmdir(os.path.join(root,dir))
                #print(f"Removed empty folder {os.path.join(root,dir)}")
    
    filtered_pdb_files = list(base_path.rglob("*.pdb"))
    filtered_pdb_files = set([str(pdb_file).split("/")[-1].split(".pdb")[0] for pdb_file in filtered_pdb_files])
    total_pdb_files = len(filtered_pdb_files)
    keep_ids_not_found = keep_ids - filtered_pdb_files
    print(f"Total .pdb files remaining in raw_data/ppi_6A: {total_pdb_files}")
    print(f"Total .pdb files not found: {len(keep_ids_not_found)}. Writing to raw_data/ppi_6A_filtered_pdbs_not_found.txt")
    
    with open("raw_data/ppi_6A_filtered_pdbs_not_found.txt", "w") as f:
        for pdb_id in keep_ids_not_found:
            f.write(f"{pdb_id}\n")

def main():
    #with open_logfile("collection_log.txt"):
    if True:
        # Download ppiref data
        download_ppiref_data()
        
        # Get PPIRef300K and PPIRef50K
        ppiref300k = download_ppiref_split_json(splitname='ppiref_6A_filtered')
        ppiref50k = download_ppiref_split_json(splitname='ppiref_6A_filtered_clustered_04')
        print("\nDownloaded PPIRef300K and PPIRef50K")
        print(f"\tPPIRef300K: {len(ppiref300k)} PPIs")
        print(f"\tPPIRef50K: {len(ppiref50k)} PPIs")
        
        # Find the union of these two, and filter down to only those files 
        keep_ids = set(ppiref300k) | set(ppiref50k)
        print(f"\tTotal unique IDs in either PPIRef300K OR PPIRef50K: {len(keep_ids)}. Keeping these IDs only.")
        filter_to_splits(keep_ids)
        
        # Download the keep_ids, but take off the chain identifiers from their names (e.g. 10gs_A_B --> 10gs)
        original_ids = list(set([x.split('_')[0].upper() for x in keep_ids]))
        original_pdbs_output_dir = "raw_data/ppi_6A_RCSB_full"
        download_original_structures(original_ids, struct_format="pdb", output_dir=original_pdbs_output_dir)
        
        # Get rid of the subfolders 
        flatten_directory("raw_data/ppi_6A")
        
if __name__ == "__main__":
    main()