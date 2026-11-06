# US-align: Flexible Domain Alignment

Structural alignment of query domains against target domains / segments with [US-align](https://github.com/pylelab/USalign). This is step **04 (flexible alignment)** of the pipeline, run alongside Foldseek and FATCAT.

> **Status:** template (v0.0)  
> **Input from:** step 03 (sliced domains)  
> **Output to:** step 05 (find consensus)

---

## Contents

1. [Installation](#1-installation)
2. [Inputs](#2-inputs)
3. [Running US-align](#3-running-us-align)
4. [Getting the data](#4-getting-the-data)
5. [Analysis](#5-analysis)
6. [Plotting](#6-plotting)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Installation

### US-align

```bash
# TODO: install / compile instructions
```

### Python dependencies (analysis and plotting)

```bash
# TODO
```

| Tool     | Version | Notes |
|----------|---------|-------|
| US-align | TODO    | output of `USalign -h` or release tag |
| Python   | TODO    | |
| pandas   | TODO    | |

---

## 2. Inputs

| Input              | Location          | Description |
|--------------------|-------------------|-------------|
| Query domains      | `query_domains`   | curated domains (`.pdb`), from `domain_slicer.py` |
| Target set (dvd)   | `target_domains`  | Merizo-segmented domains |
| Target set (dv200) | `target_segments` | 200-aa segments from `fixed-length_slicer.py` |

> TODO: note how the file lists are built (e.g. `query_list.txt`, `target_list.txt`).

---

## 3. Running US-align

### 3.1 Domain vs domain (dvd)

```bash
# TODO: command
```

| Flag | Why |
|------|-----|
| TODO | |

> TODO: record the exact parameters used for the final run.

### 3.2 Run in the background (HPC)

```bash
# TODO: nohup / SLURM example
```

### 3.3 dv200

```bash
# TODO: same as above with the 200-aa segments
```

---

## 4. Getting the data

> TODO: describe the output format and column order.

Quick checks:

```bash
# TODO: number of hits, head of output, queries with >= 1 hit
```

Convert to CSV:

```bash
# TODO: parser script or snippet
```

---

## 5. Analysis

```python
python3 usalign_batch_v3.0.py -q dir_query_set -t dir_target_set
# runs us-align through the script, with query and target directories as input 
# on KRS-lab server 
```

Cutoff notes:
- TM-score > 0.5 generally indicates the same fold.
- TODO: state the thresholds used for this project and why.

---

## 6. Plotting

```python
import pandas as pd
import matplotlib.pyplot as plt

# TODO
```

Planned plots (TODO):
- [ ] TM-score histogram: dvd vs dv200
- [ ] Best-hit TM-score per query domain (bar)
- [ ] Heatmap of query x target TM-scores

---

## 7. Troubleshooting

| Problem | Fix |
|---------|-----|
| TODO    | TODO |

---

## Output files

| File | Description |
|------|-------------|
| `results/usalign_dvd.*` | raw domain-vs-domain alignments |
| `results/usalign_dv200.*` | raw domain-vs-200aa alignments |
| `results/*_hits.csv` | filtered hits |
| `results/*_best_per_query.csv` | top hit per query |
| `results/*.png` | plots |

## References

- Zhang C, Shine M, Pyle AM, Zhang Y. US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes. *Nat Methods* 19, 1109–1115 (2022). https://doi.org/10.1038/s41592-022-01585-1
- US-align GitHub: https://github.com/pylelab/USalign
