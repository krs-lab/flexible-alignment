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

Record what you actually ended up with (this is what makes it reproducible):

- Protenix version / commit: `<git rev-parse --short HEAD>`
- Python: `<python --version>`
- PyTorch / CUDA: `<python -c "import torch; print(torch.__version__, torch.version.cuda)">`
- NVIDIA driver: `<from nvidia-smi once it works>`
- Model checkpoint used: `<name and version>`
- Weights / data location: `<path, and the env variable if one was set>`

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
