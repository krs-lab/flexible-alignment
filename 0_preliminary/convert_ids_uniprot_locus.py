#!/usr/bin/env python3
"""
Copies and renames .pdb, .pse, and PAE .json files from an input directory into a new subdirectory
(uniprot_ids or locus_names) based on a CSV mapping key file.

Usage:
    python3 convert_uniprot_locusname.py -i input_dir -k key.csv --direction u2l
    python3 convert_uniprot_locusname.py --input input_dir --key-file key.csv --direction locusname-to-uniprot
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path


def load_mapping_dict(csv_path, direction):
    """
    Reads CSV key file and builds a mapping lookup dictionary.
    Handles field matching dynamically for UniProt and Locus columns.
    """
    valid_uniprot_cols = {"uniprot_ids", "uniprot_id", "uniprot", "uniprotid"}
    valid_locus_cols = {"locus", "locus name", "locus_name", "locus tag", "locus_tag", "locusname", "locustag"}

    mapping = {}

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

        if not uniprot_col or not locus_col:
            print(
                "Error: Key CSV file must contain both a UniProt column (e.g., 'uniprot_id') "
                "and a Locus column (e.g., 'locus_name')."
            )
            sys.exit(1)

        for row in reader:
            uid = row[uniprot_col].strip()
            locus = row[locus_col].strip()

            if not uid or not locus:
                continue

            if direction in ["u2l", "uniprot-to-locusname", "unitprot-to-locusname"]:
                mapping[uid] = locus
            else:  # l2u / locusname-to-uniprot
                mapping[locus] = uid

    return mapping


def copy_and_convert_filenames(input_dir, output_dir, mapping):
    """
    Scans target directory for .pdb, .pse, and .json files and copies them to output_dir
    with new converted prefixes.
    """
    supported_extensions = {".pdb", ".pse", ".json"}
    matched_count = 0
    skipped_count = 0

    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort mapping keys by length descending to prevent partial prefix mis-matches
    sorted_keys = sorted(mapping.keys(), key=len, reverse=True)

    for filepath in sorted(input_dir.iterdir()):
        if filepath.is_file() and filepath.suffix.lower() in supported_extensions:
            filename = filepath.name
            matched_key = None

            # Find matching base identifier
            for key in sorted_keys:
                if filename == f"{key}{filepath.suffix}" or filename.startswith(f"{key}_"):
                    matched_key = key
                    break

            if matched_key:
                new_prefix = mapping[matched_key]
                new_filename = filename.replace(matched_key, new_prefix, 1)
                new_filepath = output_dir / new_filename

                if new_filepath.exists():
                    print(f"  [SKIP] Target file already exists: {new_filename}")
                    skipped_count += 1
                else:
                    shutil.copy2(filepath, new_filepath)
                    print(f"  [COPIED & RENAMED] {filename} -> {new_filename}")
                    matched_count += 1
            else:
                skipped_count += 1

    print(f"\nFinished! Copied {matched_count} file(s) to '{output_dir}'. Skipped/Unmatched: {skipped_count}.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert filenames (.pdb, .pse, .json) between UniProt IDs and Locus Names using a key CSV file.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="DIR",
        help="Path to the directory containing .pdb, .pse, or PAE .json files.",
    )

    parser.add_argument(
        "-k",
        "--key-file",
        required=True,
        type=Path,
        metavar="CSV_FILE",
        help="Path to the CSV file mapping UniProt IDs and Locus names.",
    )

    parser.add_argument(
        "-d",
        "--direction",
        required=True,
        choices=[
            "u2l",
            "l2u",
            "uniprot-to-locusname",
            "unitprot-to-locusname",
            "locusname-to-uniprot",
        ],
        help=(
            "Conversion direction:\n"
            "  u2l / uniprot-to-locusname  : Rename UniProt ID -> Locus Name (output: <input>/locus_names)\n"
            "  l2u / locusname-to-uniprot  : Rename Locus Name -> UniProt ID (output: <input>/uniprot_ids)"
        ),
    )

    args = parser.parse_args()

    if not args.input.exists() or not args.input.is_dir():
        print(f"Error: Input directory '{args.input}' does not exist or is not a directory.")
        sys.exit(1)

    if not args.key_file.exists():
        print(f"Error: Key CSV file '{args.key_file}' does not exist.")
        sys.exit(1)

    # Determine output subdirectory name based on target goal
    if args.direction in ["u2l", "uniprot-to-locusname", "unitprot-to-locusname"]:
        output_dir = args.input / "locus_names"
    else:
        output_dir = args.input / "uniprot_ids"

    mapping = load_mapping_dict(args.key_file, args.direction)
    print(f"Loaded {len(mapping)} identifier mappings from key file.")
    print(f"Processing files in '{args.input}' -> Output directory: '{output_dir}'...\n")

    copy_and_convert_filenames(args.input, output_dir, mapping)


if __name__ == "__main__":
    main()
