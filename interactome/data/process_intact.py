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

def organize_cv(path= "../data_files/raw/intact/cv/intact.obo", starting_id = "MI:0190", out_dir="../data_files/processed/intact/cv"):
    """
    Organize controlled vocabulary file for IntAct into a DataFrame, save as csv 
    """
    import pandas as pd

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
                        "parent_ids": current.get("is_a", [])
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
    df.to_csv(f"{out_dir}/{starting_id.lower().replace(':','_')}_subtree.csv", index=False)

def see_all_keys(file):
    rec = pymex.mif.Record()
    rec.parseMif(file, 'mif300' )
    with open("available_keys.log", "w") as f:
        for interaction in rec.interactions:
            f.write("--- INTERACTION ---\n")
            f.write("\n\t".join((dir(interaction))) + "\n")
            #f.write(str(interaction.__dict__) + "\n")
            f.write("EXPERIMENTS:")
            f.write(str(interaction.experiments))
            
            f.write("--- EXPERIMENT ---\n")
            f.write("\n\t".join((dir(interaction.experiment))) + "\n")

            f.write("--- PARTICIPANT 1 ---\n")
            f.write("\n\t".join((dir(interaction.participants[0]))) + "\n")
            #f.write(str(interaction.participants[0].__dict__) + "\n")

            f.write("--- INTERACTOR 1 ---\n")
            f.write("\n\t".join((dir(interaction.participants[0].interactor))) + "\n")
            #f.write(str(interaction.participants[0].interactor.__dict__) + "\n")

def get_all_paths_and_values(obj, prefix="root", visited=None, depth=0, max_depth=10):
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
            results.extend(get_all_paths_and_values(v, new_prefix, visited, depth+1, max_depth))

    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            new_prefix = f"{prefix}[{i}]"
            results.append((new_prefix, item))
            results.extend(get_all_paths_and_values(item, new_prefix, visited, depth+1, max_depth))

    elif hasattr(obj, "__dict__"):
        for attr in dir(obj):
            if attr.startswith("_") or attr in dir(builtins):
                continue
            try:
                val = getattr(obj, attr)
            except Exception:
                continue
            new_prefix = f"{prefix}.{attr}"
            results.append((new_prefix, val))
            results.extend(get_all_paths_and_values(val, new_prefix, visited, depth+1, max_depth))

    return results

def record_all_paths_and_values(file, output_dir="../data_files/processed/intact/examples"):
    """
    Records all the paths and values for the first interaction as an example. 
    """
    rec = pymex.mif.Record()
    filename = file.split("/")[-1]
    rec.parseMif(file, 'mif300')

    # Just get the first interaction
    interaction = rec.interactions[0]

    paths_and_vals = get_all_paths_and_values(interaction)

    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/xml_{filename}_access_paths_with_values.txt", "w", encoding="utf-8") as f:
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
        
def parse_psi30(files, acceptable_labels):
    """
    Parse PSI-MI 3.0 formatted files to get interaction information. 
    """
    rec = pymex.mif.Record()

    counter_of_short = 0

    interaction_label = []
    gene_name_1 = []
    gene_name_2 = []
    gene_symbol_1 = []
    gene_symbol_2 = [] 
    length_1 = []
    length_2 = []
    protein_1 = []
    protein_2 = []

    for file in files:
        rec.parseMif(file, 'mif300' )
        # try:
        for interaction in rec.interactions:
            if len(interaction.participants) == 2:
                if type(interaction.participants[0].interactor) != type(None) and type(interaction.participants[1].interactor) != type(None):
                    try:
                        # if both are proteins
                        if interaction.participants[0].interactor.type.label == 'protein' and interaction.participants[1].interactor.type.label == 'protein':
                            label = interaction.type.label
                            
                            # if it's actually a direct interaction
                            if label in acceptable_labels:
                                interaction_label.append(label)
                                # if there's a sequence
                                if type(interaction.participants[0].interactor.sequence)!= type(None):
                                    # if len(interaction.participants[0].interactor.sequence) < 1000:

                                    gene_name_1.append(interaction.participants[0].interactor.label)
                                    protein_1.append(interaction.participants[0].interactor.sequence)
                                    length_1.append(len(interaction.participants[0].interactor.sequence))

                                    try:
                                        representation = interaction.participants[0].interactor.alias.__repr__()
                                        gene_symbol_1.append(ast.literal_eval(representation)[0]['value'])
                                    except:
                                        gene_symbol_1.append(None)
                                else:
                                    try:
                                        gene_name_1.append(interaction.participants[0].interactor.label)
                                    except:
                                        gene_name_1.append(None)

                                    try:
                                        protein_1.append(interaction.participants[0].interactor.sequence)
                                    except:
                                        protein_1.append(None)

                                    try:
                                        length_1.append(len(interaction.participants[0].interactor.sequence))
                                    except:
                                        length_1.append(None)

                                    try:
                                        representation = interaction.participants[0].interactor.alias.__repr__()
                                        gene_symbol_1.append(ast.literal_eval(representation)[0]['value'])
                                    except:
                                        gene_symbol_1.append(None)

                                if type(interaction.participants[1].interactor.sequence)!= type(None):
                                    # if len(interaction.participants[1].interactor.sequence) < 1000:
                                    gene_name_2.append(interaction.participants[1].interactor.label)
                                    protein_2.append(interaction.participants[1].interactor.sequence)
                                    length_2.append(len(interaction.participants[1].interactor.sequence))

                                    try:
                                        representation = interaction.participants[1].interactor.alias.__repr__()
                                        gene_symbol_2.append(ast.literal_eval(representation)[0]['value'])
                                    except:
                                        gene_symbol_2.append(None)
                                else:
                                    try:
                                        gene_name_2.append(interaction.participants[1].interactor.label)
                                    except:
                                        gene_name_2.append(None)

                                    try:
                                        protein_2.append(interaction.participants[1].interactor.sequence)
                                    except:
                                        protein_2.append(None)

                                    try:
                                        length_2.append(len(interaction.participants[1].interactor.sequence))
                                    except:
                                        length_2.append(None)

                                    try:
                                        representation = interaction.participants[1].interactor.alias.__repr__()
                                        gene_symbol_2.append(ast.literal_eval(representation)[0]['value'])
                                    except:
                                        gene_symbol_2.append(None)

                    except Exception as e:
                        logging.error(f"Error processing interaction: {e}")
                        continue
    
    return gene_name_1, gene_name_2, gene_symbol_1, gene_symbol_2, length_1, length_2, protein_1, protein_2, interaction_label

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

    import os as os
    for path, subdirs, files in os.walk(folder):
        for name in files:
            file_name = os.path.join(path, name)
            if 'negative' not in file_name:
                posfiles.append(file_name)
            else:
                negfiles.append(file_name)
    
    return posfiles, negfiles

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
    logger.info(f"Total positive files: {len(posfiles)}\nTotal negative files: {len(negfiles)}")

    # read in the CSV file with the acceptable labels and make a label list
    mypath_terms = f"{git_root}/interactome/data_files/processed/intact/cv/mi_0190_subtree.csv"
    if not(os.path.exists(mypath_terms)):
        logger.info(f"Organizing controlled vocabulary file")
        organize_cv()
        
    acceptable_labels = pd.read_csv(mypath_terms)
    acceptable_labels = list(acceptable_labels['label'])
    
    logger.info(f"Getting all available keys from the first positive file: {posfiles[0]}")
    see_all_keys(posfiles[0])
    record_all_paths_and_values(posfiles[0])
    
    #gene_name_1, gene_name_2, gene_symbol_1, gene_symbol_2, length_1, length_2, protein_1, protein_2, interaction_label = parse_psi30(posfiles)     
    gene_name_1, gene_name_2, gene_symbol_1, gene_symbol_2, length_1, length_2, protein_1, protein_2, interaction_label= parse_psi30([posfiles[0]], acceptable_labels)                 
    test_interactome = pd.DataFrame({'protein_1': gene_name_1,'protein_2': gene_name_2, 'length_1': length_1,'length_2': length_2, 'aa_1':protein_1, 'aa_2':protein_2, 'label':interaction_label})
    test_interactome.to_csv("test_interactome.csv",index=False)
    #interactome = pd.DataFrame({'protein_1': gene_name_1,'protein_2': gene_name_2, 'length_1': length_1,'length_2': length_2, 'aa_1':protein_1, 'aa_2':protein_2, 'label':interaction_label})
    #interactome.to_csv('pmid_interactomes.csv')

    #interactome = pd.DataFrame({'gene_name_1': gene_symbol_1, 'gene_name_2': gene_symbol_2, 'protein_1': gene_name_1,'protein_2': gene_name_2, 'length_1': length_1,'length_2': length_2, 'aa_1':protein_1, 'aa_2':protein_2, 'label':interaction_label})
    #interactome.to_csv('pmid_interactomes2.csv')

if __name__ == '__main__':
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename='process_intact.log', encoding='utf-8', level=logging.DEBUG, filemode='w')

    main()