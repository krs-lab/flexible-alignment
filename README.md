## Flexible alignment 

Computaitonally screening for PIP-binding domains across *Legionella pneumophila* effector space through flexible domain alignment.  

### Steps 

- Get the latest **AlphaFold models** with their corresponding PAE (predicted alignemnt error) files (.pdb + .json) with `fetch_pdb_pae.py`. 
    - If needed convert the pdb/pae files' IDs into UniProt or Locus IDs with `convert_ids_uniprot_locus.py`

### Commands

```bash
# get pdb/pae files 
python3 fetch_pdb_pae.py -i legionella_effectors_368.csv

# convert IDs to corresponding uniprot/locus IDs
python3 convert_ids_uniprot_locus.py -i pdb_pae_dir -k key.csv -d u2l
```

