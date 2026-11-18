# Foldseek: Domain Alignment :white_check_mark:

Structural alignment of query domains against target domains / segments with [Foldseek](https://github.com/steineggerlab/foldseek). 

- **Author:** Visnu
- **Lab / PI:** Schmitz Lab/ Karl Schmitz
- **Status:** testing
- **Last updated:** 2026-11-18

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
- [x] Simplified script to run Foldseek 
- [x] 2D plot to show common 3Di and TMalign hits


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
# to enable foldseek env; as installed inside conda env 
source /fast_disk/visnu/foldseek_env.sh

## Run Foldseek in TM-align and 3Di mode to get alignment scores + superimposed files at once 
chmod +x foldseek_dual_mode.sh
./foldseek_dual_mode.sh query_domain_dir/ target_domain_dir/ <label, dvd/dv200> <threads, default 8>
    # flags/arguments used: 
    # TMP :temporary working dir and subdirs for foldseek
    # 
    # in TMalign mode: 
    # --alignment-type 1 :uses TM-align mode
    # --exhaustive-search 1 :enables all-vs-all alignment
    # -e inf :no matches are filtered out
    # --exact-tmscore 1 :exact TM-score enabled
    # --tmalign-fast 0 :disables TMalign's fast mode
    # --tmalign-hit-order 1 :orders hits relative to query structure
    # --format-output "..." :outputs specified cols written in output
    # --threads :CPU threads intend to use 
    # --tmscore-threshold 0.1 :minimum TM-score cutoff
    # --tmscore-threshold-mode 1 :mode 1 normalizes TM-score against query len
    # --format-mode 5 :mode 5 outputs structural superimposition pdb files
    # 
    # in 3Di mode:
    # --alignment-type 2 :to run 3Di mode
    # -e 1000 :E-value threshold
    # -c 0: min coverage threshold (0% coverage threshold so none discarded)
    # --format-mode 5 :mode 5 outputs structural superimposition pdb files

## Analysis    
# transfer results to the local machine 
# individualize the TSV for each domains 
python3 foldseek_output_splitter_v0.py -i <aln_file_from_each_mode.tsv> -c <query_col> -o <output_dir>

# make a dir contain info required to make 2D plot
python3 condense_foldseek_v0.1.py -i-tma <tma-mode_aln_for_individual_query.tsv> -qtm <0.4> -i-3di <3di-mode_aln_for_individual_query.tsv> -bit-raw <13> -z-bit <2> -v-list <hits_to_highlight.csv> -dpi <600>

# ===

## Individual analyses
# to plot histograms for each query (TMalign mode) 
python3 histogram_foldseek_v1.6.py -i individual_fs_output/ -tmq 0.4 -v hits_to_highlight.csv -dpi 600

# plot histogram with raw bit score >=13 and modified zscore of raw bit >=2 (3Di mode)
python3 histogram_foldseek-3di_v1.8.py -i individual_fs_output_3di/ -bit-raw 13 -z-bit 2 -v hits_to_highlight.csv -dpi 600
```

---

## References :white_check_mark:

- van Kempen M. et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* (2024). https://doi.org/10.1038/s41587-023-01773-0
- Foldseek GitHub: https://github.com/steineggerlab/foldseek
