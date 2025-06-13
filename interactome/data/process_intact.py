from os import listdir
import os
from os.path import isfile, join
import pandas as pd
from lxml import etree as ET
import logging
import pymex.mif
import xml.etree.ElementTree as ET
import builtins
from interactome.utils import get_git_root
import os
import pandas as pd
import multiprocessing
from functools import partial
from math import ceil
from datetime import datetime


def organize_cv(
    path="../data_files/raw/intact/cv/intact.obo",
    starting_id="MI:0190",
    out_dir="../data_files/processed/intact/cv",
):
    """
    Organize controlled vocabulary file for IntAct into a DataFrame, save as csv
    """

    terms = {}  # ID → {"name": ..., "parent_ids": [...]}

    with open(path, "r") as f:
        current = {}
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current = {}
            elif line == "":
                if "id" in current and "name" in current:
                    terms[current["id"]] = {
                        "name": current["name"],
                        "parent_ids": current.get("is_a", []),
                    }
            elif line.startswith("id:"):
                current["id"] = line[3:].strip()
            elif line.startswith("name:"):
                current["name"] = line[5:].strip()
            elif line.startswith("is_a:"):
                parent_id = line[5:].split("!")[0].strip()
                current.setdefault("is_a", []).append(parent_id)

    # Build reverse index: parent_id → list of child_ids
    children_map = {}
    for term_id, term in terms.items():
        for parent_id in term["parent_ids"]:
            children_map.setdefault(parent_id, []).append(term_id)

    # Recursively collect descendants starting from MI:2232
    results = []

    def traverse(term_id, parent_id=None):
        if term_id not in terms:
            return
        label = terms[term_id]["name"]
        results.append({"label": label, "id": term_id, "parent_id": parent_id})
        for child_id in children_map.get(term_id, []):
            traverse(child_id, parent_id=term_id)

    traverse(starting_id)

    df = pd.DataFrame(results)
    print(df.head())

    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(
        f"{out_dir}/{starting_id.lower().replace(':','_')}_subtree.csv", index=False
    )


def see_all_keys(file, output_dir="../data_files/processed/intact/examples"):
    rec = pymex.mif.Record()
    rec.parseMif(file, "mif300")
    with open(f"{output_dir}/available_keys.txt", "w") as f:
        for interaction in rec.interactions:
            f.write("--- INTERACTION ---\n")
            f.write("\n\t".join((dir(interaction))) + "\n")
            # f.write(str(interaction.__dict__) + "\n")

            f.write("--- EXPERIMENT ---\n")
            f.write("\n\t".join((dir(interaction.experiment))) + "\n")

            f.write("--- PARTICIPANT 1 ---\n")
            f.write("\n\t".join((dir(interaction.participants[0]))) + "\n")
            # f.write(str(interaction.participants[0].__dict__) + "\n")

            f.write("--- INTERACTOR 1 ---\n")
            f.write("\n\t".join((dir(interaction.participants[0].interactor))) + "\n")
            # f.write(str(interaction.participants[0].interactor.__dict__) + "\n")


def get_all_paths_and_values(obj, prefix="root", visited=None, depth=0, max_depth=20):
    if visited is None:
        visited = set()

    if depth > max_depth or id(obj) in visited:
        return []

    visited.add(id(obj))
    results = []

    if isinstance(obj, dict):
        for k, v in obj.items():
            new_prefix = f"{prefix}[{repr(k)}]"
            results.append((new_prefix, v))
            results.extend(
                get_all_paths_and_values(v, new_prefix, visited, depth + 1, max_depth)
            )

    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            new_prefix = f"{prefix}[{i}]"
            results.append((new_prefix, item))
            results.extend(
                get_all_paths_and_values(
                    item, new_prefix, visited, depth + 1, max_depth
                )
            )

    elif hasattr(obj, "__dict__"):
        for attr in dir(obj):
            if attr.startswith("_"):
                continue
            try:
                val = getattr(obj, attr)
            except Exception:
                continue
            new_prefix = f"{prefix}.{attr}"
            results.append((new_prefix, val))
            results.extend(
                get_all_paths_and_values(val, new_prefix, visited, depth + 1, max_depth)
            )

    return results


def record_all_paths_and_values(
    file, output_dir="../data_files/processed/intact/examples", interaction_no=0
):
    """
    Records all the paths and values for the first interaction as an example.
    """
    rec = pymex.mif.Record()
    filename = file.split("/")[-1].replace(".xml", "")
    rec.parseMif(file, "mif300")
    print(file)

    # Just get the first interaction
    interaction = rec.interactions[interaction_no]

    if len(interaction.participants) == 2:
        if type(interaction.participants[0].interactor) != type(None) and type(
            interaction.participants[1].interactor
        ) != type(None):
            try:
                # if both are proteins
                if (
                    interaction.participants[0].interactor.type.label == "protein"
                    and interaction.participants[1].interactor.type.label == "protein"
                ):
                    paths_and_vals = get_all_paths_and_values(interaction)

                    os.makedirs(output_dir, exist_ok=True)
                    with open(
                        f"{output_dir}/xml_{filename}_interaction{interaction_no}_access_paths_with_values.txt",
                        "w",
                        encoding="utf-8",
                    ) as f:
                        f.write("\n--PATHS ONLY--\n")
                        for path, val in paths_and_vals:
                            f.write(f"{path}\n")

                        f.write("\n--PATHS WITH VALUES--\n")

                        for path, val in paths_and_vals:
                            try:
                                val_str = str(val)
                                if len(val_str) > 200:
                                    val_str = val_str[:200] + "..."
                            except Exception as e:
                                val_str = f"[Error converting value: {e}]"
                            f.write(f"{path}: {val_str}\n")
                else:
                    logging.warning(
                        f"Interaction in {file} does not contain two protein interactors."
                    )
            except Exception as e:
                logging.error(f"Error processing interaction: {e}")
        else:
            logging.warning(
                f"Interaction is not between two interactors of type != None."
            )
    else:
        logging.warning(
            f"Interaction in {file} does not have exactly two participants. Has {len(interaction.participants)}"
        )


def parse_psi30_interactor(interactor):
    """
    Parse a single interactor from PSI-MI 3.0 format
    """
    gene_name = None
    gene_symbol = None
    length = None
    protein = None
    primaryref_db = None
    primaryref_id = None
    host_taxid = None
    host_cell_type = None
    host_compartment = None
    host_tissue = None
    go = []
    uniprotkb = []
    ensp = []
    ensg = []
    enst = []
    interpro = []
    rscbpdb = []

    # Sequence, name, label, symbol
    ### YES: SEQUENCE
    try:
        gene_name = interactor.label
    except:
        gene_name = None
    
    try:
        protein = interactor.sequence
        length = len(interactor.sequence)
    except:
        protein = None
        length = None

    try:
        gene_symbol = interactor.gene
    except:
        gene_symbol = None

    # Try adding host info
    try:
        host_taxid = interactor.host.taxid
    except:
        host_taxid = None
    
    try:
        host_cell_type = interactor.host.cellType
    except:
        host_cell_type = None
        
    try:
        host_compartment = interactor.host.compartment
    except:
        host_compartment = None
    
    try:
        host_tissue = interactor.host.tissue.label
    except:
        host_tissue = None

    # Identifiers: uniprotkb, ensembl (g, t, p), go, interpro, rscbpdb
    if type(interactor.primaryRef) != type(None):
        try:
            database = interactor.primaryRef.db
            primaryref_db = database
            primaryref_id = interactor.primaryRef.ac
            if database == "uniprotkb":
                uniprotkb.append(interactor.primaryRef.ac)
            elif database == "ensembl":
                if interactor.primaryRef.ac.startswith("ENSP"):
                    ensp.append(interactor.primaryRef.ac)
                elif interactor.primaryRef.ac.startswith("ENSG"):
                    ensg.append(interactor.primaryRef.ac)
                elif interactor.primaryRef.ac.startswith("ENST"):
                    enst.append(interactor.primaryRef.ac)
            elif database == "interpro":
                interpro.append(interactor.primaryRef.ac)
            elif database == "rscb pdb":
                rscbpdb.append(interactor.primaryRef.ac)
            elif database == "go":
                go.append(interactor.primaryRef.ac)
        except Exception as e:
            logger.warning(f"Error processing interactor primaryRef: {e}")

    if len(interactor.secondaryRef) > 0:
        for sec in interactor.secondaryRef:
            try:
                database = sec.db
                if database == "uniprotkb":
                    uniprotkb.append(sec.ac)
                elif database == "ensembl":
                    if sec.ac.startswith("ENSP"):
                        ensp.append(sec.ac)
                    elif sec.ac.startswith("ENSG"):
                        ensg.append(sec.ac)
                    elif sec.ac.startswith("ENST"):
                        enst.append(sec.ac)
                elif database == "interpro":
                    interpro.append(sec.ac)
                elif database == "rscb pdb":
                    rscbpdb.append(sec.ac)
                elif database == "go":
                    go.append(sec.ac)
            except Exception as e:
                logger.warning(f"Error processing interactor secondaryRef: {e}")

    ret_dict = {
        "gene_name": gene_name,
        "gene_symbol": gene_symbol,
        "length": length,
        "protein": protein,
        "uniprotkb": ",".join(uniprotkb) if len(uniprotkb)>0 else None,
        "ensp": ",".join(ensp) if len(ensp)>0 else None,
        "ensg": ",".join(ensg) if len(ensg)>0 else None,
        "enst": ",".join(enst) if len(enst)>0 else None,
        "interpro": ",".join(interpro) if len(interpro)>0 else None,
        "rscbpdb": ",".join(rscbpdb) if len(rscbpdb)>0 else None,
        "primaryref_db": primaryref_db,
        "primaryref_id": primaryref_id,
        "go": ",".join(go) if len(go)>0 else None,
        "host_taxid": host_taxid,
        "host_cell_type": host_cell_type,
        "host_compartment": host_compartment,
        "host_tissue": host_tissue
    }
    return ret_dict

def parse_psi30_experiment(experiment):
    """
    Parse a single PSI-30 experiment
    """
    new_experiment = {}
    if type(experiment.bibref.primaryRef) is not None:
        try:
            # PubMed ID
            if experiment.bibref.primaryRef.db == "pubmed":
                new_experiment["pubmed"] = experiment.bibref.primaryRef.ac
            else:
                pubmed.append(None)
                
            # Experiment info
            if experiment.method.name is not None:
                new_experiment["method"] = experiment.method.name
            
            try:
                if experiment.partmethod.name is not None:
                    new_experiment["partmethod"] = experiment.partmethod.name
            except:
                logger.warning(f"Error processing experiment partmethod: {e}")
            
            if len(experiment.hosts) > 0:
                new_experiment["hosts"] = []
                for host in experiment.hosts:
                    new_host = {}
                    # Host organism tax ID  
                    if host.taxid is not None:
                        new_host["taxid"] = host.taxid
                    
                        if not(host.taxid=="-1"):
                            # Host organism label
                            if host.label is not None:
                                new_host["label"] = host.label

                            # Host organism cell type
                            if host.cellType is not None:
                                new_host["cell_type"] = host.cellType
                            
                            # Host organism cell compartment
                            if host.compartment is not None:
                                new_host["compartment"] = host.compartment
                            
                            # Host organism tissue
                            if host.tissue.label is not None:
                                new_host["tissue"] = host.tissue.label 
                            
                            new_experiment["hosts"].append(new_host)      
                        else:
                            logger.warning(f"Host taxid is -1, not full host") 
                                         
                        
        except Exception as e:
            logging.error(f"Error processing experiment: {e}")

    return new_experiment

def parse_psi30(files, acceptable_labels):
    """
    Parse PSI-MI 3.0 formatted files to get interaction information.
    """
    rec = pymex.mif.Record()

    counter_of_short = 0

    interaction_label = []
    interaction_mi = []
    # There can be multiple experiments per interaction 
    # Each experiment entry will be a list of dictionaries
    experiments = []
    
    gene_name_1 = []
    gene_symbol_1 = []
    length_1 = []
    protein_1 = []
    uniprotkb_1 = []
    ensp_1 = []
    ensg_1 = []
    enst_1 = []
    interpro_1 = []
    rscbpdb_1 = []
    primaryref_db_1 = []
    primaryref_id_1 = []
    go_1 = []
    host_taxid_1 = []
    host_cell_type_1 = []
    host_compartment_1 = []
    host_tissue_1 = []

    gene_name_2 = []
    gene_symbol_2 = []
    length_2 = []
    protein_2 = []
    uniprotkb_2 = []
    ensp_2 = []
    ensg_2 = []
    enst_2 = []
    interpro_2 = []
    rscbpdb_2 = []
    primaryref_db_2 = []
    primaryref_id_2 = []
    go_2 = []
    host_taxid_2 = []
    host_cell_type_2 = []
    host_compartment_2 = []
    host_tissue_2 = []

    for file in files:
        rec.parseMif(file, "mif300")
        # try:
        for i, interaction in enumerate(rec.interactions):
            logger.info(f"Processing interaction {i} in {file}")
            # Binary interaction
            if len(interaction.participants) == 2:
                if type(interaction.participants[0].interactor) != type(None) and type(
                    interaction.participants[1].interactor
                ) != type(None):
                    try:
                        # if both are proteins
                        if (
                            interaction.participants[0].interactor.type.label
                            == "protein"
                            and interaction.participants[1].interactor.type.label
                            == "protein"
                        ):
                            label = interaction.type.label

                            # just make sure it's truly a molecular interaction (physical associatin, direct, etc)
                            if label in acceptable_labels:
                                interaction_label.append(label)
                                
                                # Add MI
                                if type(interaction.type.primaryRef) != type(None):
                                    try:
                                        if interaction.type.primaryRef.db == "psi-mi":
                                            interaction_mi.append(interaction.type.primaryRef.ac)
                                    except Exception as e:
                                        logging.error(f"Error getting MI: {e}")
                                        interaction_mi.append(None)
                                        
                                # Get info for interactors 1 and 2
                                info_1 = parse_psi30_interactor(interaction.participants[0].interactor)
                                info_2 = parse_psi30_interactor(interaction.participants[1].interactor)
                                
                                # Add information
                                gene_name_1.append(info_1["gene_name"])
                                gene_symbol_1.append(info_1["gene_symbol"])
                                length_1.append(info_1["length"])
                                protein_1.append(info_1["protein"])
                                uniprotkb_1.append(info_1["uniprotkb"])
                                ensp_1.append(info_1["ensp"])
                                ensg_1.append(info_1["ensg"])
                                enst_1.append(info_1["enst"])
                                interpro_1.append(info_1["interpro"])
                                rscbpdb_1.append(info_1["rscbpdb"])
                                primaryref_db_1.append(info_1["primaryref_db"])
                                primaryref_id_1.append(info_1["primaryref_id"])
                                go_1.append(info_1["go"])
                                host_taxid_1.append(info_1["host_taxid"])
                                host_cell_type_1.append(info_1["host_cell_type"])
                                host_compartment_1.append(info_1["host_compartment"])
                                host_tissue_1.append(info_1["host_tissue"])
                                
                                gene_name_2.append(info_2["gene_name"])
                                gene_symbol_2.append(info_2["gene_symbol"])
                                length_2.append(info_2["length"])
                                protein_2.append(info_2["protein"])
                                uniprotkb_2.append(info_2["uniprotkb"])
                                ensp_2.append(info_2["ensp"])
                                ensg_2.append(info_2["ensg"])
                                enst_2.append(info_2["enst"])
                                interpro_2.append(info_2["interpro"])
                                rscbpdb_2.append(info_2["rscbpdb"])
                                primaryref_db_2.append(info_2["primaryref_db"])
                                primaryref_id_2.append(info_2["primaryref_id"])
                                go_2.append(info_2["go"])
                                host_taxid_2.append(info_2["host_taxid"])
                                host_cell_type_2.append(info_2["host_cell_type"])
                                host_compartment_2.append(info_2["host_compartment"])
                                host_tissue_2.append(info_2["host_tissue"])
                                
                                # Get experiment-related info 
                                if len(interaction.experiments) > 0:
                                    cur_experiments = []
                                    for experiment in interaction.experiments:
                                        new_experiment = parse_psi30_experiment(experiment)
                                        cur_experiments.append(new_experiment)
                                    
                                    experiments.append(cur_experiments)
                                else:
                                    logging.info(f"No experiments found for interaction {i} in {file}")
                                    experiments.append(None)
                                    
                    except Exception as e:
                        logging.error(f"Error processing interaction: {e}")
                        continue
            # If not a binary interaction, see how many there are! And what the deal is
            else:
                # get the types of all the participants
                part_types = []
                for participant in interaction.participants:
                    if participant.interactor and participant.interactor.type:
                        part_types.append(participant.interactor.type.label)
                part_types = ",".join(part_types)
                logging.warning(
                    f"Interaction {i} at {file} has {len(interaction.participants)} participants.\n\tTypes: {part_types}"
                )
                continue
            
    ret_dict = {
        "interaction_label": interaction_label,
        "interaction_mi": interaction_mi,
        "experiments": experiments,
        "gene_name_1": gene_name_1,
        "gene_symbol_1": gene_symbol_1,
        "length_1": length_1,
        "protein_1": protein_1,
        "uniprotkb_1": uniprotkb_1,
        "ensp_1": ensp_1,
        "ensg_1": ensg_1,
        "enst_1": enst_1,
        "interpro_1": interpro_1,
        "rscbpdb_1": rscbpdb_1,
        "primaryref_db_1": primaryref_db_1,
        "primaryref_id_1": primaryref_id_1,
        "go_1": go_1,
        "host_taxid_1": host_taxid_1,
        "host_cell_type_1": host_cell_type_1,
        "host_compartment_1": host_compartment_1,
        "host_tissue_1": host_tissue_1,
        "gene_name_2": gene_name_2,
        "gene_symbol_2": gene_symbol_2,
        "length_2": length_2,
        "protein_2": protein_2,
        "uniprotkb_2": uniprotkb_2,
        "ensp_2": ensp_2,
        "ensg_2": ensg_2,
        "enst_2": enst_2,
        "interpro_2": interpro_2,
        "rscbpdb_2": rscbpdb_2,
        "primaryref_db_2": primaryref_db_2,
        "primaryref_id_2": primaryref_id_2,
        "go_2": go_2,
        "host_taxid_2": host_taxid_2,
        "host_cell_type_2": host_cell_type_2,
        "host_compartment_2": host_compartment_2,
        "host_tissue_2": host_tissue_2,
    }
    
    return ret_dict


def get_pos_neg_files(folder):
    """
    Get positive and negative files from the given path.

    Args:
        path (str): Path to the directory containing files.

    Returns:
        tuple: Two lists, one for positive files and one for negative files.
    """
    posfiles = []
    negfiles = []

    for path, subdirs, files in os.walk(folder):
        for name in files:
            file_name = os.path.join(path, name)
            if "negative" not in file_name:
                posfiles.append(file_name)
            else:
                negfiles.append(file_name)

    return posfiles, negfiles

def parse_and_save(file_idx, file_batch, acceptable_labels, output_dir):
    """
    Parses a batch of PSI-MI 3.0 files and saves output to a CSV file.
    """
    psi30_info = parse_psi30(file_batch, acceptable_labels)

    df = pd.DataFrame.from_dict(psi30_info)
    df = df.rename(columns={
        "gene_name_1": "protein_1",
        "gene_name_2": "protein_2",
        "protein_1": "aa_1",
        "protein_2": "aa_2",
    })

    out_path = os.path.join(output_dir, f"interactome_part_{file_idx}.csv")
    df.to_csv(out_path, index=False)
    return out_path  # Return path for tracking

def split_into_chunks(lst, n_chunks):
    """
    Splits a list into approximately equal-sized chunks.
    """
    chunk_size = ceil(len(lst) / n_chunks)
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

def run_parallel_parsing(posfiles, acceptable_labels, output_dir="../data_files/processed/intact/output_chunks"):
    os.makedirs(output_dir, exist_ok=True)

    n_cores = multiprocessing.cpu_count()
    file_chunks = split_into_chunks(posfiles, n_cores)

    with multiprocessing.Pool(processes=n_cores) as pool:
        tasks = [
            (i, chunk)
            for i, chunk in enumerate(file_chunks)
        ]
        results = pool.starmap(
            partial(parse_and_save, acceptable_labels=acceptable_labels, output_dir=output_dir),
            tasks
        )
    return results

def merge_outputs(intermediate_output_dir="../data_files/processed/intact/output_chunks", interaction_type="positive", final_output_dir="../data_files/processed/intact"):
    dfs = []
    for file in sorted(os.listdir(intermediate_output_dir)):
        if file.endswith(".csv"):
            df = pd.read_csv(os.path.join(intermediate_output_dir, file))
            dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)
    
    # get today's date 
    today = datetime.now().strftime("%Y-%m-%d")
    merged_filename = os.path.join(final_output_dir, f"intact_processed_{interaction_type}PPIs_{today}.csv")
    os.makedirs(final_output_dir, exist_ok=True)
    # Save the merged DataFrame to a CSV file
    merged.to_csv(merged_filename, index=False)
    logger.info(f"Merged output written to {merged_filename}")


def main():
    """
    Main method. Process interactome data.
    """
    # create dataset_50, species_50, and pmid_50 one by one and then merge them to the full set
    git_root = get_git_root()
    intact_folder = f"{git_root}/interactome/data_files/raw/intact"
    mypath_datasets = f"{intact_folder}/psi30/datasets/"
    mypath_species = f"{intact_folder}/psi30/species/"
    mypath_pmid = f"{intact_folder}/psi30/pmid/"

    posfiles, negfiles = get_pos_neg_files(mypath_pmid)
    logger.info(
        f"Total positive files: {len(posfiles)}\nTotal negative files: {len(negfiles)}"
    )

    # read in the CSV file with the acceptable labels and make a label list
    mypath_terms = (
        f"{git_root}/interactome/data_files/processed/intact/cv/mi_0190_subtree.csv"
    )
    if not (os.path.exists(mypath_terms)):
        logger.info(f"Organizing controlled vocabulary file")
        organize_cv()

    acceptable_labels = pd.read_csv(mypath_terms)
    acceptable_labels = list(acceptable_labels["label"])

    logger.info(
        f"Getting all available keys from the first positive file: {posfiles[0]}"
    )
    see_all_keys(posfiles[0])
    interaction_no = 1
    logger.info(
        f"Recording all paths from first positive file, interaction {interaction_no}: {posfiles[0]}"
    )
    #record_all_paths_and_values(posfiles[0], interaction_no=interaction_no)
    record_all_paths_and_values(negfiles[0], interaction_no=0)

    # Do a sample run
    example_fname = posfiles[0].split("/")[-1].replace(".xml", "")
    test_psi30_info = parse_psi30([posfiles[0]], acceptable_labels)
    test_interactome = pd.DataFrame.from_dict(test_psi30_info)
    test_interactome = test_interactome.rename(
        columns={
            "gene_name_1": "protein_1",
            "gene_name_2": "protein_2",
            "protein_1": "aa_1",
            "protein_2": "aa_2",
        }
    )
    test_interactome.to_csv(f"../data_files/processed/intact/examples/positive_{example_fname}xml_interactome.csv", index=False)

    negfile_intermediate_path = "../data_files/processed/intact/negative_output_chunks"
    posfile_intermediate_path = "../data_files/processed/intact/positive_output_chunks"
    
    # Do negative parallel processing
    #utput_files = run_parallel_parsing(negfiles, acceptable_labels, output_dir=negfile_intermediate_path)
    #merge_outputs(negfile_intermediate_path, interaction_type="negative")


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename="process_intact.log",
        encoding="utf-8",
        level=logging.DEBUG,
        filemode="w",
    )

    main()
