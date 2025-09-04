from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig
from pathlib import Path
import os
import pandas as pd
import logging

import rootutils
root = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)

logger = logging.getLogger(__name__)

def clean_biogrid(df: pd.DataFrame, output_dir) -> pd.DataFrame:
    """
    Clean the database and rename its columns to match IntAct format where possible
    """
    # Define paths to make things easier later
    output_dir = Path(root) / output_dir
    
    col_map = {
        "#BioGRID Interaction ID": "biogrid_interaction", 
        "Entrez Gene Interactor A": "entrez_1", 
        "Entrez Gene Interactor B": "entrez_2", 
        "BioGRID ID Interactor A": "biogrid_1", 
        "BioGRID ID Interactor B": "biogrid_2", 
        "Systematic Name Interactor A": "systematic_name_1", 
        "Systematic Name Interactor B": "systematic_name_2", 
        "Official Symbol Interactor A": "gene_symbol_1", 
        "Official Symbol Interactor B": "gene_symbol_2", 
        "Synonyms Interactor A": "synonyms_1", 
        "Synonyms Interactor B": "synonyms_2", 
        "Experimental System": "experimental_system", 
        "Experimental System Type": "experimental_system_type", 
        "Author": "author", 
        "Publication Source": "pubmed", 
        "Organism ID Interactor A": "organism_id_1", 
        "Organism ID Interactor B": "organism_id_2", 
        "Throughput": "throughput", 
        "Score": "score", 
        "Modification": "modifications", 
        "Qualifications": "qualifications", 
        "Tags": "tags", 
        "Source Database": "source_db", 
        "SWISS-PROT Accessions Interactor A": "uniprotkb_1", 
        "TREMBL Accessions Interactor A": "trembl_1", 
        "REFSEQ Accessions Interactor A": "refseq_1" ,
        "SWISS-PROT Accessions Interactor B": "uniprotkb_2", 
        "TREMBL Accessions Interactor B": "trembl_2", 
        "REFSEQ Accessions Interactor B": "refseq_2", 
        "Ontology Term IDs": "ontology_term_ids", 
        "Ontology Term Names": "ontology_term_names", 
        "Ontology Term Categories": "ontology_term_categories",
        "Ontology Term Qualifier IDs": "ontology_term_qualifier_ids", 
        "Ontology Term Qualifier Names": "ontology_term_qualifier_names", 
        "Ontology Term Types": "ontology_term_types", 
        "Organism Name Interactor A": "organism_name_1", 
        "Organism Name Interactor B": "organism_id_2"
    }
    
    df = df.rename(columns=col_map)
    
    # Make sure everything is either a valid entry or np.nan
    df = df.replace({"": pd.NA, " ": pd.NA, "-": pd.NA})
    
    # How many things don't have UniProt IDs? 
    logger.info(f"Total interactions: {len(df)}")
    logger.info(f"Interactions without UniProt ID for interactor 1: {df['uniprotkb_1'].isna().sum()}")
    logger.info(f"Interactions without UniProt ID for interactor 2: {df['uniprotkb_2'].isna().sum()}")
    logger.info(f"Interactions without a RefSeq for interactor 1: {df['refseq_1'].isna().sum()}")
    logger.info(f"Interactions without a RefSeq for interactor 2: {df['refseq_2'].isna().sum()}")
    logger.info(f"Interactions with multiple RefSeqs for interactor 1: {df['refseq_1'].str.contains('|').sum()}")
    logger.info(f"Interactions with multiple RefSeqs for interactor 2: {df['refseq_2'].str.contains('|').sum()}")

    # make directory for ID mapping
    idmap_dir = output_dir / "id_mapping"
    idmap_input_dir = idmap_dir / "inputs"
    idmap_output_dir = idmap_dir / "outputs"
    os.makedirs(idmap_dir, exist_ok=True)
    os.makedirs(idmap_input_dir, exist_ok=True)
    os.makedirs(idmap_output_dir, exist_ok=True)
    
    # Make input files
    idmap_inputs_trembl = pd.concat([df["trembl_1"],df["trembl_2"]]).dropna().unique().tolist()
    idmap_inputs_trembl = list(set("|".join(idmap_inputs_trembl).split("|")))
    idmap_inputs_entrez = pd.concat([df["entrez_1"],df["entrez_2"]]).dropna().unique().tolist()
    idmap_inputs_refseq = pd.concat([df["refseq_1"],df["refseq_2"]]).dropna().unique().tolist()  
    idmap_inputs_refseq = list(set("|".join(idmap_inputs_refseq).split("|")))
    idmap_inputs_refseq_batch1 = idmap_inputs_refseq[0:75000]
    idmap_inputs_refseq_batch2 = idmap_inputs_refseq[75000::]
     
    with open(idmap_input_dir / "biogrid_trembl_ids.txt", "w") as f:
        for item in idmap_inputs_trembl:
            f.write(f"{item}\n")

    with open(idmap_input_dir / "biogrid_entrez_ids.txt", "w") as f:
        for item in idmap_inputs_entrez:
            f.write(f"{item}\n")
            
    with open(idmap_input_dir / "biogrid_refseq_ids_batch1.txt", "w") as f:
        for item in idmap_inputs_refseq_batch1:
            f.write(f"{item}\n")    
    
    with open(idmap_input_dir / "biogrid_refseq_ids_batch2.txt", "w") as f:
        for item in idmap_inputs_refseq_batch2:
            f.write(f"{item}\n")

def main(cfg: DictConfig):
    # Set up needed paths
    biogrid_fpath = Path(root) / cfg.process.input_file # raw biogrid download
    example_output_dir = Path(root) / cfg.process.example_output_dir    # output directory for example
    os.makedirs(example_output_dir, exist_ok=True)
    
    biogrid_db = pd.read_csv(biogrid_fpath, sep="\t", header=0, low_memory=False)
    example_subset = biogrid_db.sample(n=1000, random_state=42).reset_index(drop=True) # Sample 1000 rows for the example
    example_subset.to_csv(example_output_dir / "biogrid_example.tsv", sep="\t", index=False)

    clean_biogrid(df=biogrid_db, output_dir=cfg.process.output_dir)
    
    return

if __name__ == "__main__":
    main()