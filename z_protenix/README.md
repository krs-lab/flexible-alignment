# Protenix batch: domain vs. PIP moiety screen

Protenix predictions to rank which protein domains show promising ipTM scores against PIP moieties.

- **Author:** Visnu
- **Lab / PI:** Schmitz Lab/ Karl Schmitz
- **Status:** testing
- **Last updated:** YYYY-MM-DD

---

## 0. To do

- [x] Reboot server and confirm GPU works
- [x] Install on lab server (WS2) 
- [ ] Test run on WS2 
- [x] Install on HPC cluster (BIOMIX)
- [ ] Test run on BIOMIX
- [ ] Install on HPC cluster + test run 
- [ ] Decide: lab server or HPC for the full batch
- [ ] Run batch
- [ ] Summarise ipTM and pick domains for bench validation

---

## 1. Log

```bash
## install on WS2 (lab server)
# 1. Navigate to your target directory
cd /fast_disk/visnu

# 2. Download and install Miniconda silently into /fast_disk/visnu/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /fast_disk/visnu/miniconda3
rm miniconda.sh

# 3. Create a clean Python 3.11 environment named protenix_env
/fast_disk/visnu/miniconda3/bin/conda create -n protenix_env python=3.11 -y

# 4. Clone the Protenix repository
git clone https://github.com/bytedance/Protenix.git
cd Protenix

# 5. Perform an editable installation using the environment's explicit pip binary
/fast_disk/visnu/miniconda3/envs/protenix_env/bin/pip install -e .

# verifying 
/fast_disk/visnu/miniconda3/envs/protenix_env/bin/python -c "
import torch
import protenix
print('CUDA Available:', torch.cuda.is_available())
print('Protenix Path:', protenix.__file__)
"

# ---

## install on BIOMIX 
# Set up local directory structure
mkdir -p /home/pritom/schmitzLab/tools
cd /home/pritom/schmitzLab/tools

# Download and install Miniconda (if not already present)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /home/pritom/schmitzLab/tools/miniconda3
rm miniconda.sh

# Source conda and create a clean Python 3.11 environment
source /home/pritom/schmitzLab/tools/miniconda3/etc/profile.d/conda.sh
conda create -n protenix_env python=3.11 -y
conda activate protenix_env

# Clone the repository (or navigate to your existing local clone)
cd /home/pritom/schmitzLab/tools
git clone https://github.com/bytedance/Protenix.git
cd Protenix

# Perform an editable installation
pip install -e .

# Save this script as /home/pritom/schmitzLab/tools/run_protenix.sbatch:
#!/bin/bash
#SBATCH --job-name=protenix_job
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=protenix_%j.out
#SBATCH --error=protenix_%j.err

# Source local conda profile and activate environment
source /home/pritom/schmitzLab/tools/miniconda3/etc/profile.d/conda.sh
conda activate protenix_env

# Run your Protenix inference or command here
protenix --help

sbatch /home/pritom/schmitzLab/tools/run_protenix.sbatch
squeue -u pritom


# ===

### Protenix: Test run 
# Prepare JSON file: 
    # containing ligand chemical strucgture and domain AA residues 

# Domain AA seq: 
    # from domain PDBs through ChmeraX
    # Open PDB in UCSF ChimeraX; in terminal select full seq 
    sequence chain #1/A 
    # save as fasta 
    save <lpg2529_l4h.fasta> alignment 1/A

# loading AA seq from FASTA containing dir 
python3 make_input_v0.py dir_having_domain_fasta/ 
    
# get ligand SMILES (Simplified Molecular Input Line Entry System) from PubChem 
    # manually entered in `make_input_v0.py`
    # PI3P:
        # query: phosphatidylinositol 3-phosphate; Compound CID:446005
    # PI4P: 
        # query: phosphatidylinositol 4-phosphate; Compound CID:643965
    # PI(4,5)P2: 
        # query: phosphatidylinositol 4,5-phosphate; Compound CID:178202108
    # PI (-ve control): 
        # query: phosphatidylinositol; Compound CID:15488999

# running on WS2 (lab server)
    # initiate protenix env 
    source /fast_disk/visnu/miniconda3/bin/activate protenix_env
    # checking GPU status 
    nvidia-smi
    # check for JSON batch file in the dir 
    ls batch.json

# running protenix engine 
PYTHONNOUSERSITE=1 CUDA_VISIBLE_DEVICES=0 nohup \
  /fast_disk/visnu/miniconda3/envs/protenix_env/bin/protenix pred \
  -i batch.json \
  -o ./protenix_out \
  --seeds 101 \
  > protenix.log 2>&1 &
    # PYTHONNOUSERSITE=1 :to make Python ignore user site dir 
    # /absolute/path/of/protenix :skips path lookup and to use the working version 
# check progress
tail -f protenix.log

# ===

### Protenix: Interpretation 
# tabulate the data 
cd protenix_out
python - <<'EOF' > scores.tsv
import json, glob, os
print("complex\tsample\trank\tiptm\tptm\tclash\tchain_pair_iptm")
for f in sorted(glob.glob("*/seed_101/predictions/*summary_confidence*.json")):
    d = json.load(open(f))
    name = f.split("/")[0]
    s = f.split("sample_")[-1].replace(".json","")
    print(name, s, round(d["ranking_score"],3), round(d["iptm"],3), round(d["ptm"],3),
          d.get("has_clash"), d.get("chain_pair_iptm"), sep="\t")
EOF
# GFP is showing an overall high ipTM with all ligands
# but for the other PIP binders there is a decreasing trend of ipTMs relative to the negative control 


```

What question is this answering? One or two sentences.

> Example: Which protein domains show high ipTM with PIP moieties, so we can prioritise them for bench validation?

## 2. Where it runs

| Machine | Path | Notes |
|---|---|---|
| Lab server (WS2) | `/fast_disk/visnu/protenix` | RTX 4070 Ti SUPER (16 GB VRAM), no sudo |
| UD HPC | `<path on cluster>` | fallback for large batches (SLURM) |

## 3. Environment

```bash
conda create -p /fast_disk/visnu/protenix_env python=3.11 -y
conda activate /fast_disk/visnu/protenix_env
git clone https://github.com/bytedance/protenix.git
cd protenix
pip install -e .
```
- Protenix version / commit: `85767b8`
- Python: `Python 3.11.11`
- PyTorch / CUDA: `2.7.1+cu126 12.6`
- NVIDIA driver: `580.178.04, NVIDIA GeForce RTX 4070 Ti SUPER`
- Model checkpoint used: `Protenix v2.0.0 official pretrained weights`
- Weights / data location: `/fast_disk/visnu/Protenix/`

Quick check that the GPU is visible:

```bash
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

## 4. Inputs

- Protein domains: `<where the sequences/domain definitions come from, e.g. domain segmentation output>`
- PIP moieties (ligands): `<how they are specified, e.g. CCD code or SMILES>`
- Input format: Protenix JSON, one entry per domain/ligand pair
- Number of jobs: `<N>`
- Size range (residues): `<min> to <max>`

## 5. Directory layout

```
.
├── README.md
├── inputs/          # JSON files and domain FASTA files
├── msa/             # precomputed MSAs (if used)
├── scripts/         # numbered scripts, run in order
│   ├── 01_make_inputs.py
│   ├── 02_run_protenix.sh
│   └── 03_collect_scores.py
├── outputs/         # predictions (do not commit large files)
├── results/         # summary tables and plots
└── logs/
```

## 6. How to run

```bash
conda activate /fast_disk/visnu/protenix_env

# 1. build inputs
python scripts/01_make_inputs.py

# 2. run predictions
bash scripts/02_run_protenix.sh

# 3. collect ipTM and other confidence scores
python scripts/03_collect_scores.py
```

Settings used (fill in): seeds `<>`, dtype `<bf16/fp32>`, diffusion samples `<>`, other flags `<>`.

## 7. Outputs

- Predictions: `outputs/<job_name>/seed_<n>/predictions/`
- Summary table: `results/iptm_summary.csv` (columns: domain, ligand, ipTM, pTM, ranking score, ...)

## 8. Results so far

| Domain | Ligand | ipTM | Notes |
|---|---|---|---|
| | | | |

Short interpretation: `<what stands out, what to validate next>`

## 9. Log of what I did

Newest first. One line per commit-worthy step.

- 2026-11-10: 


- YYYY-MM-DD: GPU not usable (driver/library mismatch), reboot requested
- YYYY-MM-DD: ...

## 10. Problems and fixes

| Problem | Fix |
|---|---|
| `nvidia-smi` driver/library mismatch | server needs a reboot (admin) |
| | |

## 12. References

- Chen X. et al. Protenix - Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction. bioRxiv (2025). https://doi.org/10.1101/2025.01.08.631967
- Protenix: https://github.com/bytedance/protenix
