#!/usr/bin/env python3
"""
Downloads PDB models and PAE JSON files from AlphaFold DB using UniProt IDs from a CSV file.
Saves them in `pdb_pae/` named using the locus tag/name (e.g., lpg0140.pdb) or UniProt ID if locus is missing.

Usage:
    python fetch_pdb_pae.py -i input.csv
    python fetch_pdb_pae.py --input input.csv
"""

import argparse
import csv
import sys
import time
from pathlib import Path
import requests

ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
OUT_DIR = Path("pdb_pae")
REQUEST_DELAY = 0.15  # Politeness delay between API calls


def fetch_alphafold_files(uniprot_id, output_name):
    """Query AlphaFold API and download PDB + PAE JSON files, naming them with output_name."""
    try:
        resp = requests.get(ALPHAFOLD_API_URL.format(uniprot_id=uniprot_id), timeout=30)
        if resp.status_code != 200 or not resp.json():
            print(f"  [MISSING] No AlphaFold model found for: {uniprot_id} ({output_name})")
            return False

        # Get latest primary prediction model entry returned by API
        data = resp.json()[0]
        pdb_url = data.get("pdbUrl")
        pae_url = data.get("paeDocUrl")

        if not pdb_url:
            print(f"  [WARN] Missing PDB URL for: {uniprot_id} ({output_name})")
            return False

        pdb_file = OUT_DIR / f"{output_name}.pdb"
        pae_file = OUT_DIR / f"{output_name}_pae.json"

        # Download PDB
        if not pdb_file.exists():
            pdb_resp = requests.get(pdb_url, timeout=30)
            if pdb_resp.status_code == 200:
                pdb_file.write_bytes(pdb_resp.content)

        # Download PAE JSON
        if pae_url and not pae_file.exists():
            pae_resp = requests.get(pae_url, timeout=30)
            if pae_resp.status_code == 200:
                pae_file.write_bytes(pae_resp.content)

        print(f"  [SUCCESS] Saved: {output_name}.pdb & {output_name}_pae.json")
        return True

    except Exception as e:
        print(f"  [ERROR] Failed to download {uniprot_id} ({output_name}): {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download PDB models and PAE JSON files from AlphaFold DB using UniProt IDs from a CSV file."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="CSV_FILE",
        help="Path to the input CSV file containing UniProt IDs and optional Locus names.",
    )

    args = parser.parse_args()
    csv_path = args.input

    if not csv_path.exists():
        print(f"Error: File {csv_path} does not exist.")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    valid_uniprot_cols = {"uniprot_ids", "uniprot_id", "uniprot"}
    valid_locus_cols = {"locus", "locus name", "locus_name", "locus tag", "locus_tag"}

    entries = []  # List of tuples: (uniprot_id, output_name)

    with open(csv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        uniprot_col = None
        locus_col = None

        if reader.fieldnames:
            for field in reader.fieldnames:
                clean_field = field.strip().lower()
                if clean_field in valid_uniprot_cols and not uniprot_col:
                    uniprot_col = field
                elif clean_field in valid_locus_cols and not locus_col:
                    locus_col = field

        if not uniprot_col:
            print("Error: Input CSV must contain a UniProt ID column (e.g., 'uniprot_ids', 'UNIPROT_IDs', 'uniprot', etc.).")
            sys.exit(1)

        seen_names = set()
        for row in reader:
            uid = row[uniprot_col].strip()
            if not uid:
                continue

            # Determine file label (locus name preferred, fallback to UniProt ID)
            out_name = uid
            if locus_col and row.get(locus_col):
                locus_val = row[locus_col].strip()
                if locus_val:
                    out_name = locus_val

            if out_name not in seen_names:
                seen_names.add(out_name)
                entries.append((uid, out_name))

    print(f"Found {len(entries)} unique target entries. Starting downloads into '{OUT_DIR}/'...\n")

    success_count = 0
    for uid, out_name in entries:
        if fetch_alphafold_files(uid, out_name):
            success_count += 1
        time.sleep(REQUEST_DELAY)

    print(f"\nFinished! Downloaded {success_count}/{len(entries)} entries.")


if __name__ == "__main__":
    main()
