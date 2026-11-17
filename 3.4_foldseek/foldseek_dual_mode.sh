#!/usr/bin/env bash
# foldseek_dual_mode.sh  <query_dir> <target_dir> <label> [threads]
# e.g. ./foldseek_dual_mode.sh query_domain_dir/ target_domain_dir/ domainwise 16
set -euo pipefail

usage() {
  echo "Usage: $0 <query_dir> <target_dir> <label> [threads=8]"
  echo "Example: $0 query_domain_dir/ target_domain_dir/ domainwise 16"
}
if [[ $# -lt 3 || "$1" == "-h" || "$1" == "--help" ]]; then
  usage; exit 0
fi

Q=$1; T=$2; LABEL=$3; THREADS=${4:-8}
OUT=results_${LABEL}
mkdir -p "$OUT"/{superimposed_pdbs_tma,superimposed_pdbs_3di}
TMP=$(mktemp -d tmp_XXXX); trap 'rm -rf "$TMP"' EXIT

TMA_HDR='query\ttarget\talntmscore\tqtmscore\tttmscore\tlddt\trmsd\tfident\talnlen\tqstart\tqend\ttstart\ttend\tqaln\ttaln'
DI_HDR='query\ttarget\tbits\tevalue\tfident\talnlen\tqstart\tqend\ttstart\ttend\tqaln\ttaln\tlddt'

TMA_OPTS="--alignment-type 1 --exhaustive-search 1 -e inf --exact-tmscore 1 --tmalign-fast 0 --tmalign-hit-order 1 --threads $THREADS"
DI_OPTS="--alignment-type 2 --exhaustive-search 1 -e 1000 -c 0 --threads $THREADS"

# ---- TMalign mode: scores
foldseek easy-search "$Q" "$T" "$TMP/tma.tsv" "$TMP/w" $TMA_OPTS \
  --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,rmsd,fident,alnlen,qstart,qend,tstart,tend,qaln,taln"
{ printf "$TMA_HDR\n"; cat "$TMP/tma.tsv"; } > "$OUT/aln_tmalign-mode_header.tsv"

# ---- TMalign mode: superimposition
foldseek easy-search "$Q" "$T" "$OUT/superimposed_pdbs_tma/sup_" "$TMP/w2" $TMA_OPTS \
  --tmscore-threshold 0.1 --tmscore-threshold-mode 1 --format-mode 5

# ---- 3Di mode: scores
foldseek easy-search "$Q" "$T" "$TMP/3di.tsv" "$TMP/w3" $DI_OPTS \
  --format-output "query,target,bits,evalue,fident,alnlen,qstart,qend,tstart,tend,qaln,taln,lddt"
{ printf "$DI_HDR\n"; cat "$TMP/3di.tsv"; } > "$OUT/aln_3di_header.tsv"

# ---- 3Di mode: superimposition
foldseek easy-search "$Q" "$T" "$OUT/superimposed_pdbs_3di/sup_3di_" "$TMP/w4" $DI_OPTS \
  --format-mode 5

echo "Done -> $OUT"
