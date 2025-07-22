# interactome
Curating PPI positive and negative interactome

## Environment setup
```
conda env create -f environment.yaml
```

```
configs/
│
├── download.yaml              # top-level config
├── download/
│   ├── biogrid.yaml
│   └── intact.yaml
├── process.yaml               # top-level config
└── process/
    ├── biogrid.yaml
    └── intact.yaml

scripts/
├── download.py                # runs appropriate .sh script based on config
└── process.py                 # runs appropriate .py or logic

data/
├── download_biogrid.sh
└── download_intact.sh
```

Full pipeline:
```
cd interactome

# Download biogrid
python scripts/download.py database=biogrid

# Download intact
python scripts/download.py database=intact
```