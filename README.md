## Flexible alignment 

Computaitonally screening for PIP-binding domains across *Legionella pneumophila* effector space through flexible domain alignment.  

### Steps 

- Get the latest **AlphaFold models** with their corresponding PAE (predicted alignemnt error) files (.pdb + .json) with `fetch_pdb_pae.py`. 
    - If needed convert the pdb/pae files' IDs into UniProt or Locus IDs with `convert_ids_uniprot_locus.py`
- Slice domains out with given start-end amino acid residues with `domain_slicer.py`
    - Domain slicer script requires biopython. Make sure it is installed before running (for any script requiring biopython)
    - Alternatively, have biopython installed in a virtual environment and run the script from that env. 
- Prepare all curated query domains in a query_dir. 
- Prepare the test set 
    - Domain-vs-domain (dvd): 
        - Use `merizo` to segment the domains of the test protein set. 
    - Domain-vs-200aa (dv200): 
        - Slice the pdbs into 200aa long segments with `fixed-length_slicer.py`

### Commands

```bash
# get pdb/pae files 
python3 fetch_pdb_pae.py -i legionella_effectors_368.csv

# convert IDs to corresponding uniprot/locus IDs
python3 convert_ids_uniprot_locus.py -i pdb_pae_dir -k key.csv -d u2l

# ----------
### query domain set 
# slice domains w coordinates
source /home/visnu/venvPyAnalysis/bin/activate # activates virtual env
python3 domain_slicer.py -i pdb_container_dir -c residue_coordinates.csv

# ----------
### dvd
# segment domains with merizo 
PDB_DIR="/dir/with/pdb/files"
TIMESTAMP=$(date +%Y_%m_%d_%H%M%S)
TSV_FILE="${PDB_DIR}/merizo_results_${TIMESTAMP}.tsv"

nohup python /fast_disk/visnu/merizo_work/Merizo/predict.py \
  -d cpu \
  -i "$PDB_DIR"/*.pdb \
  --iterate \
  --save_domains \
  > "$TSV_FILE" 2>&1 &

# check job status 
ps -p <job ID>

# manual check in the dir should show domain files being generated

# check merizo segmentation outputs 
for d in *.domains; do echo "=== $d ==="; cat "$d"; echo ""; done

# segment domains according to merizo output; min domain taken 83aa 
python3 slice_merizo_domains_v3_discontDom_minLen.py -m 83

# ----------
### dv200
python3 fixed-length_slicer.py -i dir_w_FL_pdbs/ -s 200
```

