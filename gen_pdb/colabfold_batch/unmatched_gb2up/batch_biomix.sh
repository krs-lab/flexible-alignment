#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_DIR="$SCRIPT_DIR"
OUTPUT_DIR="$SCRIPT_DIR/output"
TEMP_FASTA_DIR="$SCRIPT_DIR/.tmp_single_fastas"
LOG_FILE="$SCRIPT_DIR/colabfold_biomix_run.log"

export PATH="/home/pritom/schmitzLab/main/pipBinders/colabfold/localcolabfold/.pixi/envs/default/bin:$PATH"

# Clean and recreate temporary directory safely
rm -rf "$TEMP_FASTA_DIR"
mkdir -p "$OUTPUT_DIR" "$TEMP_FASTA_DIR"

echo "==========================================" | tee -a "$LOG_FILE"
echo "Starting ColabFold batch on BIOMIX at $(date)" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

echo "Splitting multi-FASTA sequences into single FASTA files..." | tee -a "$LOG_FILE"

# Parse multi-FASTA files and split into individual FASTA files
for ffile in "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fas; do
    [ -e "$ffile" ] || continue
    awk -v outdir="$TEMP_FASTA_DIR" '
        /^>/ {
            split($1, a, ">");
            seq_id = a[2];
            gsub(/[[:space:]\r\n]/, "", seq_id);
            outfile = outdir "/" seq_id ".fasta";
            print $0 > outfile;
            next;
        }
        {
            if (outfile != "") {
                print $0 >> outfile;
            }
        }
    ' "$ffile"
done

# Check if any fasta files were actually generated before running ColabFold
count=$(ls -1 "$TEMP_FASTA_DIR"/*.fasta 2>/dev/null | wc -l)
if [ "$count" -eq 0 ]; then
    echo "ERROR: No single FASTA files were generated in $TEMP_FASTA_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Generated $count single FASTA files. Running colabfold_batch..." | tee -a "$LOG_FILE"

# Run colabfold_batch on split sequence files
cd "$TEMP_FASTA_DIR" || exit 1
colabfold_batch . "$OUTPUT_DIR" --max-msa 256:512 --save-all >> "$LOG_FILE" 2>&1

echo "Processing output files in $OUTPUT_DIR..." | tee -a "$LOG_FILE"

# Safely rename top-ranked PDB files to <GenBank_ID>.pdb inside output/
# 	1st iteration: renames unrelaxed to gb_id.pdb
# 	2nd iteration: renames relaxed to gb_id.pdb
for unrelaxed_pdb in "$OUTPUT_DIR"/*_unrelaxed_rank_001_*.pdb; do
    [ -f "$unrelaxed_pdb" ] || continue

    filename=$(basename "$unrelaxed_pdb")
    genbank_id="${filename%%_unrelaxed*}"

    # Check if a corresponding relaxed file exists
    relaxed_pdb=$(ls "$OUTPUT_DIR"/${genbank_id}_relaxed_rank_001_*.pdb 2>/dev/null | head -n 1)

    if [ -f "$relaxed_pdb" ]; then
        target_pdb="$relaxed_pdb"
    else
        target_pdb="$unrelaxed_pdb"
    fi

    mv "$target_pdb" "$OUTPUT_DIR/${genbank_id}.pdb"
    echo "Renamed $(basename "$target_pdb") -> ${genbank_id}.pdb" >> "$LOG_FILE"
done

# Clean up rank 2-5 models and temporary split fastas
rm -f "$OUTPUT_DIR"/*_rank_00[2-5]_*.pdb
rm -rf "$TEMP_FASTA_DIR"

echo "Finished processing at $(date)" | tee -a "$LOG_FILE"
