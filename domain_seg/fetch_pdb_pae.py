#!/usr/bin/env python3
"""
Downloads PDB models and PAE JSON files from AlphaFold DB using UniProt IDs from a CSV file.
Saves them in `pdb_pae/` renamed as <uniprot_id>.pdb and <uniprot_id>_pae.json.

Usage:
    python fetch_pdb_pae.py input.csv
"""

import csv
import sys
import time
from pathlib import Path
import requests

ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
OUT_DIR = Path("pdb_pae")
REQUEST_DELAY = 0.15  # Politeness delay between API calls


def fetch_alphafold_files(uniprot_id):
    """Query AlphaFold API and download PDB + PAE JSON files."""
    try:
        resp = requests.get(ALPHAFOLD_API_URL.format(uniprot_id=uniprot_id), timeout=30)
        if resp.status_code != 200 or not resp.json():
            print(f"  [MISSING] No AlphaFold model found for: {uniprot_id}")
            return False

        data = resp.json()[0]  # Get primary prediction model
        pdb_url = data.get("pdbUrl")
        pae_url = data.get("paeDocUrl")

        if not pdb_url:
            print(f"  [WARN] Missing PDB URL for: {uniprot_id}")
            return False

        pdb_file = OUT_DIR / f"{uniprot_id}.pdb"
        pae_file = OUT_DIR / f"{uniprot_id}_pae.json"

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

        print(f"  [SUCCESS] Saved: {uniprot_id}.pdb & {uniprot_id}_pae.json")
        return True

    except Exception as e:
        print(f"  [ERROR] Failed to download {uniprot_id}: {e}")
        return False


def main():
    if len(sys.argv) < 2:
        print("Usage: python fetch_pdb_pae.py <input_csv_file>")
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    if not csv_path.exists():
        print(f"Error: File {csv_path} does not exist.")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read unique UniProt IDs from CSV
    uniprot_ids = []
    valid_column_names = {"uniprot_ids", "uniprot_id", "uniprot"}

    with open(csv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        
        # Dynamically match column name regardless of capitalization
        uniprot_col = None
        if reader.fieldnames:
            for field in reader.fieldnames:
                if field.strip().lower() in valid_column_names:
                    uniprot_col = field
                    break

        if not uniprot_col:
            print("Error: Input CSV must contain a UniProt ID column (e.g., 'uniprot_ids', 'UNIPROT_IDs', 'uniprot', etc.).")
            sys.exit(1)
            
        for row in reader:
            uid = row[uniprot_col].strip()
            if uid and uid not in uniprot_ids:
                uniprot_ids.append(uid)

    print(f"Found {len(uniprot_ids)} unique UniProt IDs. Starting downloads into '{OUT_DIR}/'...\n")

    success_count = 0
    for uid in uniprot_ids:
        if fetch_alphafold_files(uid):
            success_count += 1
        time.sleep(REQUEST_DELAY)

    print(f"\nFinished! Downloaded {success_count}/{len(uniprot_ids)} entries.")


if __name__ == "__main__":
    main()
