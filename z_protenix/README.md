# Protenix batch: domain vs. PIP moiety screen

Protenix predictions to rank which protein domains show promising ipTM scores against PIP moieties.

- **Author:** Visnu
- **Lab / PI:** Schmitz Lab/ Karl Schmitz
- **Status:** testing
- **Last updated:** YYYY-MM-DD

---

## 0. To do

- [x] Reboot server and confirm GPU works
- [ ] Test run on small batch on WS2 
- [ ] If successful run on full/larger set 
- [ ] Else, install on BIOMIX 
- [ ] Test run 
- [ ] Run on full/actual data set 
- [ ] Decide: lab server or HPC for the full batch
- [ ] Run batch
- [ ] Summarise ipTM and pick domains for bench validation

---

## 1. Goal

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

- YYYY-MM-DD: installed Protenix in `/fast_disk/visnu/protenix_env`
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
