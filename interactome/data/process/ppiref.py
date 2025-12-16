# Script for cleaning the parsed PPIRef data and keeping only the clean rows
"""
This is the original file from fusonflow directory. Have not changed it other than removing print and open_logfile
"""
import pandas as pd
import numpy as np
import os
import re
import json

def remove_batch_files(batch_file_dir="processed_data/ppi_6A"):
    """
    Remove all batch files (.pkl and .json) in batch_file_dir, assuming they have already been aggregated 

    Args:
        batch_file_dir (str, optional): _description_. Defaults to "processed_data/ppi_6A".
    """
    
    full_results_json = f"{batch_file_dir}/processed_6A_results.json"
    full_results_pickle = f"{batch_file_dir}/processed_6A_results.pkl"

    # If the aggregated files exist, we are free to remove the batch files
    if os.path.exists(full_results_json) and os.path.exists(full_results_pickle):
        files = os.listdir(batch_file_dir)
        pattern = r'^processed_6A_results_batch_\d+\.(json|pkl)$'
        matched_files = [f for f in files if re.match(pattern, f)]
        
        if len(matched_files)>0:
            print(f"Removing all batch files in {batch_file_dir}...")
            # Remove any file that ends in _batch_N.pkl or _batch_N.json where N can be anything
            for filename in matched_files:
                os.remove(os.path.join(batch_file_dir, filename))
                print(f"\tRemoved {filename}")

def get_alphabetical_id(row, target_chain_col="Target", binder_chain_col="Binder"):
    pdb_id = row["PDB_ID"]
    target = row[target_chain_col]
    binder = row[binder_chain_col]
    if target < binder:
        return pdb_id + "_" + target + "_" + binder
    else:
        return pdb_id + "_" + binder + "_" + target

def get_correct_chain(row, provided_name_col="Target_chain_provided", error_col="Target_chain_errors"):
    if pd.isna(row[error_col]):
        return row[provided_name_col]
    else:
        error_str = row[error_col]
        return error_str.split('true_')[1].strip()
    
def combine_motifs(motif_list):
    # Step 1: Flatten the list and extract unique motifs
    motif_set = set(motif for sublist in motif_list for motif in sublist.split(','))

    # Step 2: Parse motifs into residue and position
    motifs_by_position = {}
    for motif in motif_set:
        residue, position = motif.split('_')
        position = int(position)
        if position in motifs_by_position:
            motifs_by_position[position].add(residue)
        else:
            motifs_by_position[position] = {residue}
    
    # Step 3: Check for conflicting motifs
    bad_positions = [pos for pos, residues in motifs_by_position.items() if len(residues) > 1]
    if bad_positions:
        bad_motifs = [
            f"{','.join(residues)}_{pos}"
            for pos, residues in motifs_by_position.items() if pos in bad_positions
        ]
        print(f"WARNING: {len(bad_positions)} positions have multiple conflicting motifs: \n\t{bad_motifs}")
    
    # Step 4: Sort and combine unique motifs
    sorted_motifs = sorted(
        (f"{list(residues)[0]}_{pos}" for pos, residues in motifs_by_position.items()),
        key=lambda x: int(x.split('_')[1])
    )
    return ",".join(sorted_motifs)

def make_flipped_target_binder_df(df):
    """
    Flip the order of the chains in the dfframe, effectively doubling it.
    Change Chain1 --> Target and Chain 2 --> Binder

    Args:
        df (pandas.DataFrame): some version of processed_6A
    """
    # Make a copy and log initial stats
    flipped_df = df.copy(deep=True)  # make a copy: this one will get flipped 
    print(f"Rows in df: {len(df)}\tRows in flipped_df: {len(flipped_df)}\tEqual length (should be): {len(df) == len(flipped_df)}")

    # Change the columns to Target and Binder, in opposite orders for df and flipped_df
    flipped_df.columns = flipped_df.columns.str.replace("Chain1", "Binder").str.replace("Chain2", "Target")
    df.columns = df.columns.str.replace("Chain1", "Target").str.replace("Chain2","Binder")
    df = pd.concat([df, flipped_df], ignore_index=True)
    print(f"Rows in df+flipped_df: {len(df)}")
    
    # Reset the full_name row, and also makea  full_name_unique column to group together things that came from the same PPIRef file
    df["full_name"] = df["PDB_ID"] + "_" + df["Target"] + "_" + df["Binder"]
    df["full_name_unique"] = df.apply(lambda row: get_alphabetical_id(row), axis=1) # unique string where the alphabetically smaller chain is first 
    df['full_sequence'] = df[f'Target_og_sequence'] + '|' + df[f'Binder_og_sequence']             # e.g. MLKJF|MEAPL
    df['full_offset'] = df[f'Target_offset'].astype(str) + '|' + df[f'Binder_offset'].astype(str) # e.g. A|B
    df["full_seqs_and_motifs"] = df["full_sequence"] + "|" + df[f"Target_motifs"] + "|" + df[f"Binder_motifs"]       # e.g. MLKJF|MEAPL|M_1,K_3|L_2,L_3
    
    columns_to_not_consider = ['PDB_ID', 'Target', 'Binder', 'full_name','full_name_unique']
    columns_to_consider = df.columns.drop(columns_to_not_consider).to_list()
    df = df.drop_duplicates(subset=columns_to_consider, keep='first').reset_index(drop=True)
    print(f"Kept only the first of multiple rows with the same of ALL of the following:\n\t{','.join(columns_to_consider)}")
    print(f"Rows in df: {len(df)}")
    return df
    
def check_for_duplicate_ppis(df, chain1_name="Chain1", chain2_name="Chain2", motif_suffix="motifs", offset_suffix="offset", drop_analysis_columns=False):
    """
    Find (and print out information about) any rows that appear to be duplicates: same chain 1 and 2 sequence, and possibly same motifs

    Args:
        df (pandas.DataFrame): a version of the processed 6A dataframe
        chain1_name (str): the name of chain 1 in columns referring to it (usually Chain1, but later will be 'target')
        chain2_name (str): the name of chain 2 in columns referring to it (usually Chain2, but later will be 'binder')
        motif_suffix (str): the suffix of the columns that contain the motifs (e.g. 'motifs', 'offset_motifs)
    """
    print(f"DUPLICATE PPI CHECK{'-'*100}")
    # Add columns needed for the analysis
    analysis_columns = ['full_sequence','full_seqs_and_motifs','full_offset']
    df['full_sequence'] = df[f'{chain1_name}_og_sequence'] + '|' + df[f'{chain2_name}_og_sequence']             # e.g. MLKJF|MEAPL
    if offset_suffix is not None: 
        df['full_offset'] = df[f'{chain1_name}_{offset_suffix}'].astype(str) + '|' + df[f'{chain2_name}_{offset_suffix}'].astype(str) # e.g. A|B
    else: 
        df['full_offset'] = [np.nan]*len(df)
    df["full_seqs_and_motifs"] = df["full_sequence"] + "|" + df[f"{chain1_name}_{motif_suffix}"] + "|" + df[f"{chain2_name}_{motif_suffix}"]       # e.g. MLKJF|MEAPL|M_1,K_3|L_2,L_3
    
    # dups: rows with the same full_sequence and potentially different other information 
    dup_sequences = df.loc[df.duplicated(subset=['full_sequence'])]['full_sequence'].to_list()
    dups = df.loc[df['full_sequence'].isin(dup_sequences)]
    dups = dups.groupby("full_sequence").agg(
        full_name=("full_name", ",".join),
        PDB_ID=("PDB_ID", set),
        full_offset=("full_offset", set),
    ).reset_index()
    dups[[f'{chain1_name}_og_sequence',f'{chain2_name}_og_sequence']] = dups['full_sequence'].str.split("|", expand=True)
    dups["homomer"] = dups[f"{chain1_name}_og_sequence"] == dups[f"{chain2_name}_og_sequence"]
    dups["PDB_ID"] = dups["PDB_ID"].apply(lambda x: ",".join(list(x)))
    dups["n_PDB_IDs"] = dups["PDB_ID"].str.count(",") + 1
    dups["n_full_names"] = dups["full_name"].str.count(",") + 1
    print("\nInteractions with same (1) Chain 1 sequence, (2) Chain 2 sequence")
    print(f"Total unique interactions affected: {sum(dups['n_full_names'])}")
    print(f"Total unique PDBs affected: {len(set(','.join(dups['PDB_ID'].tolist()).split(',')))}")
    if offset_suffix is not None:
        dups["full_offset"] = dups["full_offset"].apply(lambda x: ",".join(list(x)))
        dups["n_unique_offsets"] = dups["full_offset"].str.count(",") + 1
        print(f"Total cases where the sequences are the same but the offsets are different: {len(dups[dups['n_unique_offsets'] > 1])}")

    # superdups: rows with the same chain 1 and chain 2 sequence, and same chain 1 and chain 2 motifs 
    dup_ppis = df.loc[df.duplicated(subset=['full_seqs_and_motifs'])]['full_seqs_and_motifs'].to_list()
    superdups = df.loc[df['full_seqs_and_motifs'].isin(dup_ppis)]
    superdups = superdups.groupby("full_seqs_and_motifs").agg(
        full_name=("full_name", ",".join),
        PDB_ID=("PDB_ID", set),
    ).reset_index()
    superdups[[f'{chain1_name}_og_sequence',f'{chain2_name}_og_sequence', f'{chain1_name}_motifs', f'{chain2_name}_motifs']] = superdups['full_seqs_and_motifs'].str.split("|", expand=True)
    superdups["homomer"] = superdups[f"{chain1_name}_og_sequence"] == superdups[f"{chain2_name}_og_sequence"]
    superdups["PDB_ID"] = superdups["PDB_ID"].apply(lambda x: ",".join(list(x)))
    superdups["n_PDB_IDs"] = superdups["PDB_ID"].str.count(",") + 1
    superdups["n_full_names"] = superdups["full_name"].str.count(",") + 1
    print("\nInteractions with same (1) Chain 1 sequence, (2) Chain 2 sequence, (3) Chain 1 motifs, (4) Chain 2 motifs")
    print(f"Total unique interactions affected: {sum(superdups['n_full_names'])}")
    print(f"Total unique PDBs affected: {len(set(','.join(superdups['PDB_ID'].tolist()).split(',')))}")

    ## What is different about these??? Who's in one and not the other?
    dups_full_ids = set(",".join(dups["full_name"].tolist()).split(","))
    superdups_full_ids = set(",".join(superdups["full_name"].tolist()).split(","))
    dups_not_superdups = list(dups_full_ids.difference(superdups_full_ids))
    print(f"\nUnique PDBs in dups: {len(dups_full_ids)}")
    print(f"Unique PDBs in superdups: {len(superdups_full_ids)}")
    print(f"Unique PDBs in both: {len(dups_full_ids.intersection(superdups_full_ids))}")
    print(f"Unique PDBs in dups but not in superdups: {len(dups_not_superdups)} total. First 10: {dups_not_superdups[0:10]}")
    print(f"{'-'*100}")
    
    if drop_analysis_columns:
        df = df.drop(columns=analysis_columns)

def check_data_coverage(df, fatal_errors):
    """
    Checks whether the combination of all IDS (e.g. 6xc4_H_L) in df and in fatal_errors is all the IDs in ppiref_6A_filtered.json
    This should be true, because all of these IDs were submitted to the parse.py script. If false, something was erroneously lost along the way. 
    
    Args:
        df (pd.DataFrame): holds results from parse.py
        fatal_errors (list): holds fatal errors from parse.py
    """
    fatal_error_ids = [line.split(":")[0] for line in fatal_errors]
    success_ids = df["full_name"].unique().tolist()
    
    with open("raw_data/ppiref_6A_filtered.json", "r") as file:
        original_ids = json.load(file)["folds"]["whole"]
        
    # Check that the union of success_ids and fatal_error_ids is equal to original_ids
    covered_all_ids = set(success_ids).union(set(fatal_error_ids)) == set(original_ids)
    no_successes_and_fatals = len( set(success_ids).intersection(set(fatal_error_ids)) ) == 0
    if covered_all_ids and no_successes_and_fatals:
        print("Successful parse!\n\t1. All IDs in ppiref_6A_filtered.json were successfully parsed\n\t2. No IDs in parsing_fatal_errors.txt appear in the processed results DataFrame.")
    else:
        print("WARNING")
        print(f"covered_all_ids: {covered_all_ids}\tno_successes_and_fatals: {no_successes_and_fatals}")

def check_error_columns(df, chain1_name="Chain1", chain2_name="Chain2"):
    # Add any columns that are no longer there for unnecessary analyses
    cols_to_delete = []
    for col in [f'{chain1_name}_og_missing', f'{chain2_name}_og_missing', 
                f'{chain1_name}_residue_errors', f'{chain2_name}_residue_errors']:
        if col not in list(df.columns):
            print(f"{col} not in columns... adding [np.nan]*len(df) temporarily")
            df[col] = [np.nan]*len(df)
            cols_to_delete.append(col)
    
    ### Check out the X's and missing sequences situation
    print(f"\nInvestigating sequences containing X's")
    n_seqs_with_X = len(df.loc[
    (df[f'{chain1_name}_og_sequence'].str.contains('X')) | 
    (df[f'{chain2_name}_og_sequence'].str.contains('X'))
    ])
    n_seqs_with_missing = len(df.loc[
        (df[f'{chain1_name}_og_missing'].notna()) | 
        (df[f'{chain2_name}_og_missing'].notna())
    ])
    # who has X's but does not have any missing positions? 
    n_chain1_X_not_missing = len(df.loc[
        (
            (df[f'{chain1_name}_og_sequence'].str.contains('X')) &  # sequence has an X
            (df[f'{chain1_name}_og_missing'].isna())    # sequence has no missing spots
        )
    ])
    n_chain2_X_not_missing = len(df.loc[
        (
            (df[f'{chain2_name}_og_sequence'].str.contains('X')) &  # sequence has an X
            (df[f'{chain2_name}_og_missing'].isna())    # sequence has no missing spots
        )
    ])
    print(f"\tTotal rows with X's: {n_seqs_with_X}\tTotal rows with X's because of a missing positions: {n_seqs_with_missing}")
    print(f"\tTotal chain 1 sequences with X's but none are from a missing position: {n_chain1_X_not_missing}")
    print(f"\tTotal chain 2 sequences with X's but none are from a missing position: {n_chain2_X_not_missing}")
    
    # How many sequences have insertions? 
    n_chain1_with_insertions = len(df.loc[
        df[f'{chain1_name}_og_insertions'].notna()
    ])
    n_chain2_with_insertions = len(df.loc[
        df[f'{chain2_name}_og_insertions'].notna()
    ])
    print(f"\tTotal rows where Chain 1 has insertions: {n_chain1_with_insertions}")
    print(f"\tTotal rows where Chain 2 has insertions: {n_chain2_with_insertions}")

    # How many sequences have residue errors? 
    n_chain1_with_residue_errors = len(df.loc[
        df[f'{chain1_name}_residue_errors'].notna()
    ])
    n_chain2_with_residue_errors = len(df.loc[
        df[f'{chain2_name}_residue_errors'].notna()
    ])
    print(f"\tTotal rows where Chain 1 has residue errors: {n_chain1_with_residue_errors}")
    print(f"\tTotal rows where Chain 2 has residue errors: {n_chain2_with_residue_errors}")

    # How many sequences have chain errors? 
    n_chain1_with_chain_errors = len(df.loc[
        df[f'{chain1_name}_chain_errors'].notna()
    ])
    n_chain2_with_chain_errors = len(df.loc[
        df[f'{chain2_name}_chain_errors'].notna()
    ])
    print(f"\tTotal rows where Chain 1 has chain errors: {n_chain1_with_chain_errors}")
    print(f"\tTotal rows where Chain 2 has chain errors: {n_chain2_with_chain_errors}")
    
    # How many sequences have both sequences, have no residues missing, and no residue errors in either chain? 
    print(f"\nFinding clean rows:")
    n_clean = len(df.loc[
        (df[f'{chain1_name}_og_sequence'].str.len()>0) & (df[f'{chain2_name}_og_sequence'].str.len()>0) &   # sequence
        (df[f'{chain1_name}_og_missing'].isna()) & (df[f'{chain2_name}_og_missing'].isna()) &   # nothing missing
        (df[f'{chain1_name}_residue_errors'].isna()) & (df[f'{chain2_name}_residue_errors'].isna())    # no residue errors
    ])
    print(f"\tClean: both Chains 1 and 2 have sequences, no missing residues, and no residue errors = {n_clean}")
    # How many sequences have both sequences, have no residues missing, no insertions, and no residue errors in either chain?
    n_super_clean = len(df.loc[
        (df[f'{chain1_name}_og_sequence'].str.len()>0) & (df[f'{chain2_name}_og_sequence'].str.len()>0) &   # sequence
        (df[f'{chain1_name}_og_missing'].isna()) & (df[f'{chain2_name}_og_missing'].isna()) &   # nothing missing
        (df[f'{chain1_name}_og_insertions'].isna()) & (df[f'{chain2_name}_og_insertions'].isna()) &     # no insertions
        (df[f'{chain1_name}_residue_errors'].isna()) & (df[f'{chain2_name}_residue_errors'].isna())    # no residue errors
    ])
    print(f"\tSuper clean: both Chains 1 and 2 have sequences, no missing residues, no insertions, and no residue errors = {n_super_clean}")
    
    # How many sequences have all the above ^ criteria AND no chain errors?
    n_perfect = len(df.loc[
        (df[f'{chain1_name}_og_sequence'].str.len()>0) & (df[f'{chain2_name}_og_sequence'].str.len()>0) &   # sequence
        (df[f'{chain1_name}_og_missing'].isna()) & (df[f'{chain2_name}_og_missing'].isna()) &   # nothing missing
        (df[f'{chain1_name}_og_insertions'].isna()) & (df[f'{chain2_name}_og_insertions'].isna()) &     # no insertions
        (df[f'{chain1_name}_residue_errors'].isna()) & (df[f'{chain2_name}_residue_errors'].isna()) &   # no residue errors
        (df[f'{chain1_name}_chain_errors'].isna()) & (df[f'{chain2_name}_chain_errors'].isna())
    ])
    print(f"\tPerfect: both Chains 1 and 2 have sequences, no missing residues, no insertions, and no residue errors = {n_perfect}")

    # remove any columns that we added
    if len(cols_to_delete)>0:
        print(f"Removing the temporarily added columns {cols_to_delete}.")
        df = df.drop(columns=cols_to_delete)
    
def check_homomers(df, chain1_name="Chain1", chain2_name="Chain2"):
    """
    Check how many interactions are between two identical chains

    Args:
        df (pd.DataFrame): holds results from parse.py
    """
    
    df['Homomer'] = df[f'{chain1_name}_og_sequence'] == df[f'{chain2_name}_og_sequence']
    n_homomers = len(df.loc[
        df['Homomer']
    ])
    print(f"\tTotal rows where both chains are the same (homomers): {n_homomers}")
    
def get_data_overview(df, chain1_name="Chain1", chain2_name="Chain2"):
    """
    Check out what's inside of processed_6A_results
    """
    print(f"DATA OVERVIEW {'-'*150}")
    joined_columns = "\n\t".join(list(df.columns))
    print(f"Columns: {joined_columns}")
    print(f"Dataset size: {len(df)}")
    print(f"Unique PDB IDs: {len(df['PDB_ID'].unique())}")
    
    check_error_columns(df, chain1_name=chain1_name, chain2_name=chain2_name)
    
    check_homomers(df, chain1_name=chain1_name, chain2_name=chain2_name)
    
    print('-'*164)
    
def combine_duplicates(df):
    """
    Combine rows that have the same Target and Binder sequence. 

    Args:
        df (pandas.DataFrame): processed version of PPIRef ppi_6A results in Target/Binder format, rather than Chain1/Chain2 format
    """
    
    # We shouldn't have any missing sequences OR residue errors here. 
    assert df['Target_og_missing'].isna().all()
    assert df['Binder_og_missing'].isna().all()
    assert df['Target_residue_errors'].isna().all()
    assert df['Binder_residue_errors'].isna().all()

    df = df.drop(columns=['Target_og_missing','Binder_og_missing','Target_residue_errors','Binder_residue_errors',
                        'Target_motifs','Binder_motifs','Target_offset', 'Binder_offset', 'Target_start_id',
                        'Binder_start_id','n_Target_motifs', 'n_Binder_motifs', 'full_offset',
                        'full_seqs_and_motifs','Homomer'])

    # Change names for clarity
    # _provided means it's the chain in the name of the original PPIRef file
    # _correct means it's the chain we were actually able to find and process in the PDB file 
    df = df.rename(columns={
        'Target': 'Target_chain_provided',
        'Binder': 'Binder_chain_provided',
        'full_name': 'full_name_provided',
        'full_name_unique': 'full_name_provided_unique',
    })

    df['Target_chain_correct'] = df.apply(lambda row: get_correct_chain(row, provided_name_col="Target_chain_provided", error_col="Target_chain_errors"), axis=1)
    df['Binder_chain_correct'] = df.apply(lambda row: get_correct_chain(row, provided_name_col="Binder_chain_provided", error_col="Binder_chain_errors"), axis=1)
    df['full_name_correct'] = df['PDB_ID'] + '_' + df['Target_chain_correct'] + '_' + df['Binder_chain_correct']
    df['full_name_correct_unique'] = df.apply(lambda row: get_alphabetical_id(row, target_chain_col="Target_chain_correct", binder_chain_col="Binder_chain_correct"), axis=1)

    # Combine duplicate rows into one groped dataframe 
    for col in ['Target_og_insertions','Binder_og_insertions','Target_chain_errors','Binder_chain_errors']:
        df[col] = df[col].fillna("none")

    grouped_df = df.groupby(['full_sequence']).agg({
        'full_name_correct': lambda x: '|'.join(x), # this should be the ID used to identify the entry
        'PDB_ID': lambda x: '|'.join(x.unique()), 
        'Target_offset_motifs': lambda x: combine_motifs(x),
        'Binder_offset_motifs': lambda x: combine_motifs(x),
        #'Target_offset_motifs': lambda x: ','.join(x),
        #'Binder_offset_motifs': lambda x: ','.join(x),
        'Target_og_sequence': lambda x: '|'.join(x),
        'Binder_og_sequence': lambda x: '|'.join(x),
        'Target_og_insertions': lambda x: '|'.join(x),
        'Binder_og_insertions': lambda x: '|'.join(x),
        'Target_chain_errors': lambda x: '|'.join(x),
        'Binder_chain_errors': lambda x: '|'.join(x),
        #'Target_chain_provided': lambda x: '|'.join(x),
        #'Binder_chain_provided': lambda x: '|'.join(x),
        'full_name_provided': lambda x: '|'.join(x),
        'full_name_provided_unique': lambda x: '|'.join(x),
        'full_name_correct_unique': lambda x: '|'.join(x),
    }).reset_index(drop=True) # do not want the full_sequence column in the end

    # if there are no insertions, we don't need a bunch of none's. Only keep the none's if they denote the ratio of yes error/no error entries
    for col in ['Target_og_insertions','Binder_og_insertions','Target_chain_errors','Binder_chain_errors']:
        grouped_df[col] = grouped_df[col].apply(lambda x: np.nan if set(x.split("|"))=={'none'} else x)
        
    return grouped_df
   
def main():
    #with open_logfile("cleaning_log.txt"):
    if True:
        # Start by removing unnecessary files from processed data, if needed
        remove_batch_files()
        
        # Read in the processed data 
        processed_data_path = "processed_data/ppi_6A/processed_6A_results.json"
        with open(processed_data_path, "r") as f:
            data = json.load(f)
        
        # Convert the data to a dataframe
        df = pd.DataFrame(data)
        
        list_columns=['Chain1_motifs', 'Chain2_motifs', 'Chain1_offset_motifs', 'Chain2_offset_motifs', 'Chain1_og_missing', 'Chain2_og_missing', 'Chain1_og_insertions', 'Chain2_og_insertions']
        list_comment_columns = ['Chain1_residue_errors', 'Chain2_residue_errors', 'Chain1_chain_errors', 'Chain2_chain_errors']
        for col in list_columns:
            df[col] = df[col].apply(lambda x: ','.join([str(y) for y in x]))
        for col in list_comment_columns:
            df[col] = df[col].apply(lambda x: '|'.join(x))
        
        # Now replace all the empty lists and None's with np.nan
        df = df.replace('',np.nan)
        df = df.fillna(np.nan)  # extra safeguard
        print(f"\nPreview of PPIRef parsing results:")
        print(df.head())
        # Add some interesting columns
        df['n_Chain1_motifs'] = df['Chain1_motifs'].apply(lambda x: x.count(',')+1)
        df['n_Chain2_motifs'] = df['Chain2_motifs'].apply(lambda x: x.count(',')+1) 
        df['full_name'] = df['PDB_ID'] + '_' + df['Chain1'] + '_' + df['Chain2']
        df.head(1000).to_csv("temp_processed_6A.csv",index=False)

        get_data_overview(df)
        
        # Read in the fatal errors from parsing
        with open("parsing_fatal_errors.txt","r") as f:
            fatal_errors = f.readlines()
        
        # Check whether all of the sequences either failed or succeeded, and neither both failed AND suceeded(hopefully yes!)
        check_data_coverage(df, fatal_errors)

        # Clean up the dataset
        clean_round1 = df.loc[
        (df['Chain1_og_sequence'].str.len()>0) & (df['Chain2_og_sequence'].str.len()>0) &   # sequence
        (df['Chain1_og_missing'].isna()) & (df['Chain2_og_missing'].isna()) &   # nothing missing
        (df['Chain1_residue_errors'].isna()) & (df['Chain2_residue_errors'].isna())    # no residue errors
        ]
        
        # check again
        print("\nRound 1 cleaning: selected rows with sequences for both chains, no missing residues, no residue errors.")
        print("Post-cleaning overview:")
        get_data_overview(clean_round1)
        check_for_duplicate_ppis(clean_round1)
        
        clean_round1.to_csv("processed_data/clean_6A_results.csv",index=False)
        
        # flip and double the dataset
        target_binder_df = make_flipped_target_binder_df(clean_round1)
        check_for_duplicate_ppis(target_binder_df, chain1_name="Target", chain2_name="Binder", motif_suffix="motifs")
        
        # now, do a second round of cleaning
        # combine all rows that have the same Target and Binder sequences and combine all possible motifs 
        target_binder_df.iloc[0:1000,:].to_csv("temp_target_binder_df.csv",index=False)
        clean_round2 = combine_duplicates(target_binder_df)
        # Check for duplicates, and drop_analysis_columns=True because we don't need full_sequence, etc. 
        print(f"{'-'*100}")
        print(f"Round 2 cleaning: flipped and doubled database, then combined rows where Target and Binder were the same")
        print("Post-cleaning overview:")  
        get_data_overview(clean_round2, chain1_name="Target", chain2_name="Binder") 
        
        # quickly make sure there are no duplicates
        clean_round2['full_sequence'] = clean_round2['Target_og_sequence'] + '|' + clean_round2['Binder_og_sequence']
        print(f"Total duplicates: {len(clean_round2.loc[clean_round2.duplicated(subset=['full_sequence'])])}")
        clean_round2 = clean_round2.drop(columns=['full_sequence'])
        
        clean_round2.to_csv("processed_data/all_clean_unique_6A_targets_binders.csv",index=False)   
        
if __name__ == "__main__":
    main()