# SNOOPPI: Sequence-Normalized Database of On- and Off-Target Protein Protein Interactions
This repository holds code for curating the SNOOPPI PPI database. 

## Environment setup
Run the following command to build a conda environment with all required dependencies.
```
conda env create -f environment.yaml
```
Make sure to activate the environment before running any code: 
```
conda activate interact
```

```
configs/
│
├── download.yaml              # top-level config
├── download/
│   └── intact.yaml
├── process.yaml               # top-level config
└── process/
    └── intact.yaml

scripts/
├── download.py                # runs appropriate .sh script based on config
└── process.py                 # runs appropriate .py or logic
└── clean_intact.ipynb         # final filtering of IntAct database

data/
└── download_intact.sh
```

To replicate the SNOOPPI curation pipeline, please run the following code.
```
cd interactome

# Download intact
python scripts/download.py database=intact

# Run initial processing on IntAct - extract data from XML files
python scripts/process.py database=intact
```

To complete the processing pipeline, run all cells in `scripts/clean_intact.ipynb`