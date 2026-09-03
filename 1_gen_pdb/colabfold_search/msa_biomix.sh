#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_DIR="$SCRIPT_DIR"
OUTPUT_DIR="$SCRIPT_DIR/msa_output"
TEMP_FASTA_DIR="$SCRIPT_DIR/.tmp_single_fastas"
LOG_FILE="$SCRIPT_DIR/colabfold_search_biomix.log"

export PATH="/home/pritom/schmitzLab/main/pipBinders/colabfold/localcolabfold/.pixi/envs/default/bin:$PATH"
export MMSEQS_IGNORE_INDEX=1

# Clean and recreate temporary directories
rm -rf "$TEMP_FASTA_DIR"
mkdir -p "$OUTPUT_DIR" "$TEMP_FASTA_DIR"

echo "==========================================" | tee -a "$LOG_FILE"
echo "Starting ColabFold MSA search on BIOMIX at $(date)" | tee -a "$LOG_FILE"
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

count=$(ls -1 "$TEMP_FASTA_DIR"/*.fasta 2>/dev/null | wc -l)
if [ "$count" -eq 0 ]; then
    echo "ERROR: No single FASTA files were generated in $TEMP_FASTA_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Generated $count single FASTA files. Running colabfold_search..." | tee -a "$LOG_FILE"

# Run colabfold_search across all split fastas against the system database
colabfold_search \
    "$TEMP_FASTA_DIR" \
    /mnt/dbases/colabfold_db \
    "$OUTPUT_DIR" >> "$LOG_FILE" 2>&1

# Clean up temporary split fastas
rm -rf "$TEMP_FASTA_DIR"

echo "MSA generation complete. Outputs saved in $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Finished processing at $(date)" | tee -a "$LOG_FILE"
