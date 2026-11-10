# Foldseek: Domain Alignment

Structural alignment of query domains against target domains / segments with [Foldseek](https://github.com/steineggerlab/foldseek). This is step **04 (flexible alignment)** of the pipeline, run alongside FATCAT and US-align.

- **Author:** Visnu
- **Lab / PI:** Schmitz Lab/ Karl Schmitz
- **Status:** testing
- **Last updated:** YYYY-MM-DD

---

## 0. To do

- [x] Install on lab server (WS2)
- [x] File transfer 
- [x] Test run on WS2 - TM-align mode
- [x] Test run on WS2 - 3Di mode
- [ ] Run batch on DvD 
- [ ] Run batch on Dv200aa
- [ ] Plot with histograms 
- [ ] Summarise to present/bench validation


---


## 1. Log

```bash
### file transfer 
# remote to local 
rsync -avP krslab@server.ip:/remote/location/address/ .

# local to remote 
rsync -avP local_dir krslab@server.ip:/remote/location/address/

# ---

### Foldseek 
## TM-align mode 
# get alignemnt scores 
rm -rf tmp
foldseek easy-search query_domains2_curated target_domains5_merizo aln.tsv tmp \
  --alignment-type 1 --exhaustive-search 1 -e inf \
  --exact-tmscore 1 --tmalign-fast 0 \
  --tmalign-hit-order 1 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,rmsd,fident,alnlen,qstart,qend,tstart,tend,qaln,taln"

# add header to the .tsv 
printf 'query\ttarget\talntmscore\tqtmscore\tttmscore\tlddt\trmsd\tfident\talnlen\tqstart\tqend\ttstart\ttend\tqaln\ttaln\n' | cat - aln.tsv > aln_header.tsv

# get superimposition
mkdir -p superimposed_pdbs
rm -rf tmp2
foldseek easy-search query_domains2_curated target_domains5_merizo superimposed_pdbs_ tmp2 \
  --alignment-type 1 --exhaustive-search 1 -e inf \
  --exact-tmscore 1 --tmalign-fast 0 \
  --tmalign-hit-order 1 \
  --tmscore-threshold 0.1 --tmscore-threshold-mode 1 \
  --format-mode 5

## 3Di mode 


```


## 1. Installation

### In progress

```bash
in progress
```
> in progress

### Python dependencies (analysis and plotting)

```bash
in progress
```

| Tool       | Version | Notes |
|------------|---------|-------|
| Foldseek   | TODO    | output of `foldseek version` |
| Python     | TODO    | |
| pandas     | TODO    | |

---

## 2. Inputs

| Input                 | Location              | Description |
|-------                |----------             |-------------|
| Query domains         | `query_domains`       | curated domains (`.pdb`), from `domain_slicer.py` |
| Target set (dvd)      | `target_domains`      | Merizo-segmented domains |
| Target set (dv200)    | `target_segments`     | 200-aa segments from `fixed-length_slicer.py` |

Foldseek reads `.pdb` files directly, so no conversion is needed.

---

## 3. Running Foldseek

### 3.2 Search (`easy-search`)

```bash
mkdir -p results tmp

foldseek easy-search \
    query_domains/ \
    targetDB \
    results/foldseek_dvd.m8 \
    tmp/ \
    --alignment-type 1 \
    --exhaustive-search 1 \
    -e 10 \
    --tmscore-threshold 0.0 \
    --threads 8 \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,fident,alnlen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits,rmsd"
```

| Flag                      | Why |
|------                     |-----|
| `--alignment-type 1`      | TM-align based alignment (default `2` is 3Di+AA) |
| `--exhaustive-search 1`   | skips the prefilter and aligns all query-target pairs, which matters for distant orthologs |
| `-e 10`                   | loose E-value so weak hits are kept |
| `--tmscore-threshold 0.0` | no TM-score cutoff, so filtering happens in analysis |
| `--format-output` | columns for downstream analysis |

> TODO: record the exact parameters used for the final run.

### 3.3 Run in the background (HPC)

```bash
nohup foldseek easy-search query_domains/ targetDB results/foldseek_dvd.m8 tmp/ \
    --alignment-type 1 --exhaustive-search 1 -e 10 --threads 8 \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,fident,alnlen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits,rmsd" \
    > results/foldseek_dvd.log 2>&1 &

ps -p <job ID>
tail -f results/foldseek_dvd.log
```

### 3.4 dv200

Same as above with the 200-aa segments:

```bash
foldseek createdb target_segments/ targetDB_200
foldseek easy-search query_domains/ targetDB_200 results/foldseek_dv200.m8 tmp/ \
    --alignment-type 1 --exhaustive-search 1 -e 10 --threads 8 \
    --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,fident,alnlen,qstart,qend,tstart,tend,qcov,tcov,evalue,bits,rmsd"
```

---

## 4. Getting the data

Output is a tab-separated `.m8` file with no header. Column order follows `--format-output`:

```
query  target  alntmscore  qtmscore  ttmscore  lddt  fident  alnlen  qstart  qend  tstart  tend  qcov  tcov  evalue  bits  rmsd
```

Quick checks:

```bash
wc -l results/foldseek_dvd.m8                      # number of hits
head -5 results/foldseek_dvd.m8
cut -f1 results/foldseek_dvd.m8 | sort -u | wc -l  # queries with >= 1 hit
```

Add a header and convert to CSV:

```bash
python3 foldseek_m8_to_csv.py -i results/foldseek_dvd.m8 -o results/foldseek_dvd.csv
```

> TODO: add `foldseek_m8_to_csv.py`, or use the snippet in section 5.

---

## 5. Analysis

```python
import pandas as pd

cols = ["query", "target", "alntmscore", "qtmscore", "ttmscore", "lddt",
        "fident", "alnlen", "qstart", "qend", "tstart", "tend",
        "qcov", "tcov", "evalue", "bits", "rmsd"]

df = pd.read_csv("results/foldseek_dvd.m8", sep="\t", names=cols)

# Use the max of query- and target-normalised TM-score per pair
df["tm_max"] = df[["qtmscore", "ttmscore"]].max(axis=1)

# Best hit per query
best = df.sort_values("tm_max", ascending=False).groupby("query").head(1)

# Candidate hits (TODO: confirm cutoffs)
hits = df[(df["tm_max"] >= 0.5) & (df["evalue"] <= 1e-3)]

hits.to_csv("results/foldseek_dvd_hits.csv", index=False)
best.to_csv("results/foldseek_dvd_best_per_query.csv", index=False)
```

Cutoff notes:
- TM-score > 0.5 generally indicates the same fold.
- TODO: state the thresholds used for this project and why.

---

## 6. Plotting

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("results/foldseek_dvd_hits.csv")

# Histogram of TM-scores
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(df["tm_max"], bins=30, edgecolor="black")
ax.set_xlabel("TM-score (max of query/target normalised)")
ax.set_ylabel("Number of hits")
ax.set_title("Foldseek: domain vs domain")
fig.tight_layout()
fig.savefig("results/foldseek_dvd_tm_hist.png", dpi=300)

# TM-score vs RMSD
fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(df["rmsd"], df["tm_max"], s=8, alpha=0.5)
ax.set_xlabel("RMSD (Å)")
ax.set_ylabel("TM-score")
fig.tight_layout()
fig.savefig("results/foldseek_dvd_tm_vs_rmsd.png", dpi=300)
```

Planned plots (TODO):
- [ ] TM-score histogram: dvd vs dv200
- [ ] Best-hit TM-score per query domain (bar)
- [ ] Heatmap of query x target TM-scores

---

## 7. Troubleshooting

| Problem | Fix |
|---------|-----|
| `Illegal instruction` | CPU lacks AVX2, so use the SSE2 binary |
| Empty output | loosen `-e`, add `--exhaustive-search 1`, check that the target DB was built |
| Slow runs | increase `--threads`, or split the query set |
| Missing hits for short domains | short domains give low-confidence scores, so check `qcov`/`tcov` |

---

## Output files

| File | Description |
|------|-------------|
| `results/foldseek_dvd.m8` | raw domain-vs-domain hits |
| `results/foldseek_dv200.m8` | raw domain-vs-200aa hits |
| `results/*_hits.csv` | filtered hits |
| `results/*_best_per_query.csv` | top hit per query |
| `results/*.png` | plots |

## References :white_check_mark:

- van Kempen M. et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* (2024). https://doi.org/10.1038/s41587-023-01773-0
- Foldseek GitHub: https://github.com/steineggerlab/foldseek
