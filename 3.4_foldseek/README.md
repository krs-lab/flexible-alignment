# Foldseek: Domain Alignment :white_check_mark:

Structural alignment of query domains against target domains / segments with [Foldseek](https://github.com/steineggerlab/foldseek). 

- **Author:** Visnu
- **Lab / PI:** Schmitz Lab/ Karl Schmitz
- **Status:** testing
- **Last updated:** YYYY-MM-DD

---

## 0. To do :white_check_mark:

- [x] Install on lab server (WS2)
- [x] File transfer 
- [x] Test run on WS2 (TMalign)
- [x] Test run on WS2 (3Di)
- [x] Run batch on DvD (TMaling + 3Di)
- [x] Run batch on Dv200aa (TMalign + 3Di)
- [x] Plot with histograms (TMalign + 3Di)
- [x] Summarise to present/bench validation


---

## 1. Log :white_check_mark:

```bash
# ===

### file transfer 
# remote to local 
rsync -avP krslab@server.ip:/remote/location/address/ .

# local to remote 
rsync -avP local_dir krslab@server.ip:/remote/location/address/

# ===

### Foldseek 
# enable foldseek env; as installed inside conda env 
source /fast_disk/visnu/foldseek_env.sh

## Scores: TMalign mode
rm -rf tmp
foldseek easy-search query_domain_dir/ target_domain_dir/ aln_tmalign-mode.tsv tmp \
  --alignment-type 1 --exhaustive-search 1 -e inf \
  --exact-tmscore 1 --tmalign-fast 0 \
  --tmalign-hit-order 1 \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,rmsd,fident,alnlen,qstart,qend,tstart,tend,qaln,taln"
    # --alignment-type 1 :uses TM-align mode 
    # --exhaustive-search 1 :enables all-vs-all alignment 
    # -e inf :no matches are filtered out 
    # --exact-tmscore 1 :exact TM-score enabled 
    # --tmalign-fast 0 :disables TMalign's fast mode 
    # --tmalign-hit-order 1 :orders hits relative to query structure 
    # --format-output "..." :outputs specified cols written in output 

# add header to the .tsv; clean up 
printf 'query\ttarget\talntmscore\tqtmscore\tttmscore\tlddt\trmsd\tfident\talnlen\tqstart\tqend\ttstart\ttend\tqaln\ttaln\n' | cat - aln_tmalign-mode.tsv > aln_tmalign-mode_header.tsv
rm -rf tmp

# ---

# Superimposition: TMalign mode
rm -rf tmp2
foldseek easy-search query_domain_dir/ target_domain_dir/ sup_ tmp2 \
  --alignment-type 1 --exhaustive-search 1 -e inf \
  --exact-tmscore 1 --tmalign-fast 0 \
  --tmalign-hit-order 1 \
  --tmscore-threshold 0.1 --tmscore-threshold-mode 1 \
  --format-mode 5
    # sup_ :output prefix string 
    # tmp2  :temporary working dir for foldseek 
    # --tmscore-threshold 0.1 :minimum TM-score cutoff
    # --tmscore-threshold-mode 1 :mode 1 normalizes TM-score against query len
    # --format-mode 5 :mode 5 outputs structural superimposition pdb files

# create target directory, move all generated pdbs inside; clean up tmps
mkdir -p superimposed_pdbs
mv sup_*.pdb superimposed_pdbs/
rm -rf tmp2

# ===

### analysis 
# split the foldseek output for each query 
python3 foldseek_output_splitter_v0.py -i aln_fs_tma_header.tsv -c query_col -o individual_fs_output_dir

# plot histograms for each query 
python3 histogram_foldseek_v1.6.py -i individual_fs_output/ -tmq 0.4 -v hits_to_highlight.csv -dpi 600

# ===

## Alignment scores: 3Di mode 
rm -rf tmp
foldseek easy-search query_domain_dir/ target_domain_dir/ aln_3di_mode.tsv tmp \
--alignment-type 2 \
--exhaustive-search 1 -e 1000 -c 0 \
--format-output "query,target,bits,evalue,fident,alnlen,qstart,qend,tstart,tend,qaln,taln,lddt"
    # --alignment-type 2 :to run 3Di mode
    # --exhaustive-search 1 :enables all-vs-all
    # -e 1000 :E-value threshold 
    # -c 0: min coverage threshold (0% indicates none discarded)

# adding header to csv 
printf 'query\ttarget\tbits\tevalue\tfident\talnlen\tqstart\tqend\ttstart\ttend\tqaln\ttaln\tlddt\n' | cat - aln_3di.tsv > aln_3di_header.tsv

# ---

# superimposition: 3Di mode
rm -rf tmp2
foldseek easy-search query_domain_dir/ target_domain_dir/ sup_3di_ tmp2 \
  --alignment-type 2 \
  --exhaustive-search 1 -e 1000 -c 0 \
  --format-mode 5
    # --format-mode 5 :mode 5 outputs structural superimposition pdb files

# create target directory, move all generated pdbs inside; clean up tmps
mkdir -p superimposed_pdbs_3di
mv sup_*.pdb superimposed_pdbs_3di/
rm -rf tmp2

# ===

### analysis 
# split the .tsv same as before 
# plot histogram with raw bit score >=13 and modified zscore of raw bit >=2
python3 histogram_foldseek-3di_v1.8.py -i individual_fs_output_3di/ -bit-raw 13 -z-bit 2 -v hits_to_highlight.csv -dpi 600

# ===

## 2D plot: tmalign vs 3di 
python3 condense_foldseek_v0.1.py -i-tma foldseek_results_for_query1_tma.tsv -qtm 0.4 -i-3di foldseek_results_for_query1_3di.tsv -bit-raw 13 -z-bit 2 -v-list hits_to_highlight.csv -dpi 600

```

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
