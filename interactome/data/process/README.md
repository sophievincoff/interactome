Processing each data to have shared columns so they can be combined easily 

# IntAct
```
interaction_label,
interaction_mi
experiments
year
protein_1
gene_symbol_1
length_1
aa_1
uniprotkb_1
ensp_1
ensg_1
enst_1
interpro_1
rscbpdb_1
primaryref_db_1
primaryref_id_1
go_1
host_taxid_1
host_cell_type_1
host_compartment_1
host_tissue_1
protein_2
gene_symbol_2
length_2
aa_2
uniprotkb_2
ensp_2
ensg_2
enst_2
interpro_2
rscbpdb_2
primaryref_db_2
primaryref_id_2
go_2
host_taxid_2
host_cell_type_2
host_compartment_2
host_tissue_2
```

UniProt ID mapping for IntAct was done on October 15th, 2025. 
87,256 IDs mapped to 87,417 results
* Downloaded: FASTA canonical and isoform
* Downloaded: TSV with default columns, plus Chain and peptide-related columns from the PTM/Processing section (Peptide, Propeptide, Signal Peptide, Transit Peptide)


# PepMLM data
Here is what PepNN said: 
"To generate the fragment dataset, we scanned all protein–protein complex interfaces in the protein databank (PDB) that were deposited before 2018/04/30 using the Peptiderive Rosetta protocol25 to identify protein fragments of length 5-25 amino acids that contribute to a large portion of the complex interface energy (Supplementary Fig. 1). These fragment-protein complexes were filtered based on their estimated interface energy as well as the buried surface area to ensure that they had binding properties that were reasonably close to that of peptide-protein complexes."