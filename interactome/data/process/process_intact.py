import os
import pandas as pd
from lxml import etree as ET
import logging
import pymex.mif
import xml.etree.ElementTree as ET
import os
import pandas as pd
import multiprocessing
from functools import partial
from math import ceil
from datetime import datetime
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

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
    logger.info(df.head())

    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(
        f"{out_dir}/{starting_id.lower().replace(':','_')}_subtree.csv", index=False
    )


def see_all_keys(file, output_dir:Path):
    output_dir = Path(output_dir)
    
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
    file, output_dir:Path, interaction_no=0
):
    """
    Records all the paths and values for the first interaction as an example.
    """
    output_dir = Path(output_dir)
    rec = pymex.mif.Record()
    filename = file.split("/")[-1].replace(".xml", "")
    rec.parseMif(file, "mif300")
    logger.info(file)

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
                    logger.warning(
                        f"Interaction in {file} does not contain two protein interactors."
                    )
            except Exception as e:
                logger.error(f"Error processing interaction: {e}")
        else:
            logger.warning(
                f"Interaction is not between two interactors of type != None."
            )
    else:
        logger.warning(
            f"Interaction in {file} does not have exactly two participants. Has {len(interaction.participants)}"
        )

def query_all_interaction_keys(interaction):
    """
    Helper method for debug mode. Run this to see all the keys in an interaction. 
    """
    # parse the interaction - every key 
    try:
      logger.info(f"\tInteraction alias: {interaction.alias}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction alias: {e}")

    try:  
      logger.info(f"\tInteraction attribs: {interaction.attribs}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction attribs: {e}")

    try:
      logger.info(f"\tInteraction availability: {interaction.availability}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction availability: {e}")
    
    try:
      logger.info(f"\tInteraction confidence: {interaction.confidence}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction confidence: {e}")

    try:
      logger.info(f"\tInteraction direct: {interaction.direct}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction direct: {e}")

    try:
      logger.info(f"\tInteraction expand: {interaction.expand}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction expand: {e}")

    try:
      logger.info(f"\tInteraction experiment: {interaction.experiment}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction experiment: {e}")

    try:
      logger.info(f"\tInteraction experimentCount: {interaction.experimentCount}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction experimentCount: {e}")

    try:
      logger.info(f"\tInteraction experiments: {interaction.experiments}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction experiments: {e}")

    try:
      logger.info(f"\tInteraction gene: {interaction.gene}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction gene: {e}")

    try:
      logger.info(f"\tInteraction getExperiment: {interaction.getExperiment}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction getExperiment: {e}")

    try:
      logger.info(f"\tInteraction getParticipant: {interaction.getParticipant}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction getParticipant: {e}")
    
    try:
      logger.info(f"\tInteraction imexid: {interaction.imexid}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction imexid: {e}")

    try:
      logger.info(f"\tInteraction intramolecular: {interaction.intramolecular}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction intramolecular: {e}")

    try:
      logger.info(f"\tInteraction label: {interaction.label}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction label: {e}")

    try:  
      logger.info(f"\tInteraction modelled: {interaction.modelled}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction modelled: {e}")

    try:
      logger.info(f"\tInteraction name: {interaction.name}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction name: {e}")

    try:
      logger.info(f"\tInteraction params: {interaction.params}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction params: {e}")
    
    try:
      logger.info(f"\tInteraction participantCount: {interaction.participantCount}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction participantCount: {e}")
    
    try:
      logger.info(f"\tInteraction participants: {interaction.participants}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction participants: {e}")
    
    try:
      logger.info(f"\tInteraction physical: {interaction.physical}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction physical: {e}")
    
    try:
      logger.info(f"\tInteraction primaryRef: {interaction.primaryRef}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction primaryRef: {e}")
    
    try:
      logger.info(f"\tInteraction secondaryRef: {interaction.secondaryRef}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction secondaryRef: {e}")
    
    try:
      logger.info(f"\tInteraction source: {interaction.source}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction source: {e}")
    
    try:
      logger.info(f"\tInteraction type: {interaction.type}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction type: {e}")
    
    try:
      logger.info(f"\tInteraction types: {interaction.types}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction types: {e}")
    
    try:
      logger.info(f"\tInteraction xrefCount: {interaction.xrefCount}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction xrefCount: {e}")
    
    try:
      logger.info(f"\tInteraction xrefs: {interaction.xrefs}")
    except Exception as e:
      logger.info(f"\tERROR processing Interaction xrefs: {e}")
        
def run_debug_mode(debug_file: str, acceptable_labels: list, output_dir: Path):
    """
    Debug mode: pick one xml file and parse with lots of print statements to understand what is going wrong.
    """
    logger.info(f"--- Running DEBUG mode on: {debug_file} ---")
    example_fname = Path(debug_file).stem
    rec = pymex.mif.Record()
    rec.parseMif(debug_file, "mif300")

    logger.info(f"File has {len(rec.interactions)} interactions.")
    for i, interaction in enumerate(rec.interactions):
        logger.info(f"\nInteraction {i}")
        logger.info(f"  Interaction imexid: {interaction.imexid}")
        logger.info(f"  Interaction label: {interaction.label}")
        logger.info(f"  Number of participants: {len(interaction.participants)}")
        try:
            label = interaction.type.label
            logger.info(f"  Interaction type label: {label}")
        except:
            logger.info("  No interaction type label found.")

        if len(interaction.participants) == 2:
            p1 = interaction.participants[0]
            p2 = interaction.participants[1]
            try:
                logger.info(f"\tInteractor 1 type: {p1.interactor.type.label}")
                logger.info(f"\tInteractor 2 type: {p2.interactor.type.label}")
            except:
                logger.info("\tCould not access interactor type.")
        else:
            logger.info("\tInteraction does not have exactly two participants.")
            for part in interaction.participants:
                try:
                    logger.info(f"\t\tParticipant type: {part.interactor.type.label}")
                except:
                    logger.info("\t\tCould not access participant interactor type.")

                try:
                    part_info = parse_psi30_interactor(part.interactor)
                    logger.info(f"\t\tParticipant gene symbol and uniprot: {part_info['gene_symbol']}, {part_info['uniprotkb']}")
                except Exception as e:
                    logger.error(f"\t\tError parsing participant interactor: {e}")

        query_all_interaction_keys(interaction)
    
    
    logger.info("\n--- Running parse_psi30 ---")
    test_info = parse_psi30([debug_file], acceptable_labels)

    df = pd.DataFrame.from_dict(test_info).rename(columns={
        "gene_name_1": "protein_1",
        "gene_name_2": "protein_2",
        "protein_1": "aa_1",
        "protein_2": "aa_2",
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    debug_path = output_dir / f"DEBUG_{example_fname}_parsed.csv"
    df.to_csv(debug_path, index=False)
    logger.info(f"Saved debug output to {debug_path}")

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
                new_experiment["pubmed"] = None
                
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
                            if host.tissue is not None:
                                if host.tissue.label is not None:
                                    new_host["tissue"] = host.tissue.label 
                            
                            new_experiment["hosts"].append(new_host)      
                        else:
                            logger.warning(f"Host taxid is -1, not full host")             
                        
        except Exception as e:
            logger.error(f"Error processing experiment: {e}")

    return new_experiment

def parse_psi30(files, acceptable_labels):
    """
    Parse PSI-MI 3.0 formatted files to get interaction information.
    """
    rec = pymex.mif.Record()

    interaction_label = []
    interaction_mi = []
    # There can be multiple experiments per interaction 
    # Each experiment entry will be a list of dictionaries
    experiments = []
    year = []
    
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
        # get the year. Example file path ends in psi30/pmid/2003/14609943.xml
        fyear = int(file.split("psi30/pmid/")[-1].split("/")[0])
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
                                year.append(fyear)  # we need the year for each interaction
                                
                                # Add MI
                                if type(interaction.type.primaryRef) != type(None):
                                    try:
                                        if interaction.type.primaryRef.db == "psi-mi":
                                            interaction_mi.append(interaction.type.primaryRef.ac)
                                    except Exception as e:
                                        logger.error(f"Error getting MI: {e}")
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
                                    logger.info(f"No experiments found for interaction {i} in {file}")
                                    experiments.append(None)
                    
                    except Exception as e:
                        logger.error(f"Error processing interaction: {e}")
                        continue
            # If not a binary interaction, see how many there are! And what the deal is
            else:
                # get the types of all the participants
                part_types = []
                for participant in interaction.participants:
                    if participant.interactor and participant.interactor.type:
                        part_types.append(participant.interactor.type.label)
                part_types = ",".join(part_types)
                logger.warning(
                    f"Interaction {i} at {file} has {len(interaction.participants)} participants.\n\tTypes: {part_types}"
                )
                continue
            
    ret_dict = {
        "interaction_label": interaction_label,
        "interaction_mi": interaction_mi,
        "experiments": experiments,
        "year": year,
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
    # for safety, leave one out
    n_cores = n_cores-1
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

def merge_outputs(intermediate_output_dir: Path, interaction_type: str, final_output_dir: Path):
    dfs = []

    intermediate_output_dir = Path(intermediate_output_dir)
    final_output_dir = Path(final_output_dir)

    for file in sorted(intermediate_output_dir.glob("*.csv")):
        df = pd.read_csv(file)
        dfs.append(df)

    if not dfs:
        logger.warning(f"No CSV files found in {intermediate_output_dir}. Skipping merge.")
        return

    merged = pd.concat(dfs, ignore_index=True)

    today = datetime.now().strftime("%Y-%m-%d")
    merged_filename = final_output_dir / f"intact_processed_{interaction_type}PPIs_{today}.csv"
    final_output_dir.mkdir(parents=True, exist_ok=True)

    merged.to_csv(merged_filename, index=False)
    logger.info(f"Merged output written to {merged_filename}")


def main(cfg: DictConfig):
    """
    Main method. Process interactome data.
    """
    mode = cfg.process.mode
    logger.info(f"Running IntAct processor in mode: {mode}")

    intact_folder = Path(root) / cfg.process.input_dir
    mypath_pmid = intact_folder / "psi30/pmid"
    mypath_terms = Path(root) / cfg.process.label_terms_file
    pos_output_dir = Path(root) / cfg.process.pos_output_dir
    neg_output_dir = Path(root) / cfg.process.neg_output_dir
    example_output_dir = Path(root) / cfg.process.example_output_dir

    posfiles, negfiles = get_pos_neg_files(mypath_pmid)
    logger.info(
        f"Total positive files: {len(posfiles)}\nTotal negative files: {len(negfiles)}"
    )

    # read in the CSV file with the acceptable labels and make a label list
    if not mypath_terms.exists():
        logger.info("Generating controlled vocabulary terms...")
        organize_cv()

    acceptable_labels = pd.read_csv(mypath_terms)["label"].tolist()
    
    # If running in debug mode, just do the debug 
    if cfg.process.debug:
        run_debug_mode(
            debug_file=str(Path(root) / cfg.process.debug_file),
            acceptable_labels=pd.read_csv(mypath_terms)["label"].tolist(),
            output_dir=Path(HydraConfig.get().runtime.output_dir)
        )
        return

    if mode in ["positive", "all"]:
        logger.info(f"Seeing keys in first positive file: {posfiles[0]}")
        see_all_keys(posfiles[0], output_dir=example_output_dir)
        record_all_paths_and_values(posfiles[0], output_dir=example_output_dir, interaction_no=1)
        
        logger.info("Running sample positive file")
        example_fname = Path(posfiles[0]).stem
        test_info = parse_psi30([posfiles[0]], acceptable_labels)
        df = pd.DataFrame.from_dict(test_info).rename(
            columns={
                "gene_name_1": "protein_1",
                "gene_name_2": "protein_2",
                "protein_1": "aa_1",
                "protein_2": "aa_2",
            }
        )
        example_output_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(example_output_dir / f"positive_{example_fname}_interactome.csv", index=False)

        logger.info("Running parallel processing for positives")
        run_parallel_parsing(posfiles, acceptable_labels, output_dir=pos_output_dir)
        merge_outputs(
            intermediate_output_dir=pos_output_dir,
            interaction_type="positive",
            final_output_dir=Path(root) / cfg.process.final_output_dir
        )

    if mode in ["negative", "all"]:
        logger.info("Recording values in first negative file")
        record_all_paths_and_values(negfiles[0], output_dir=example_output_dir, interaction_no=0)

        logger.info("Running parallel processing for negatives")
        run_parallel_parsing(negfiles, acceptable_labels, output_dir=neg_output_dir)
        merge_outputs(
            intermediate_output_dir=neg_output_dir,
            interaction_type="negative",
            final_output_dir=Path(root) / cfg.process.final_output_dir
        )

if __name__ == "__main__":
    main()
