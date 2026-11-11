#!/usr/bin/env python3
import argparse
import csv
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Split a Foldseek TSV output file into separate files for each unique query domain, keeping headers intact."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input TSV file."
    )
    parser.add_argument(
        "-c", "--column",
        default="query",
        help="Column name or 0-based index representing the query domain. Default is 'query'."
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=".",
        help="Directory to save the split TSV files. Default is the current directory."
    )

    args = parser.parse_args()

    input_file = args.input
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        with open(input_file, mode="r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            
            if not header:
                print("Error: Input file is empty or missing a header.", file=sys.stderr)
                sys.exit(1)

            # Determine column index
            col_arg = args.column
            col_idx = None
            
            try:
                col_idx = int(col_arg)
            except ValueError:
                for idx, h in enumerate(header):
                    if h.lower() == col_arg.lower():
                        col_idx = idx
                        break
                if col_idx is None:
                    for idx, h in enumerate(header):
                        if col_arg.lower() in h.lower():
                            col_idx = idx
                            break

            if col_idx is None or col_idx >= len(header):
                print(f"Error: Column '{col_arg}' not found in headers: {header}", file=sys.stderr)
                sys.exit(1)

            print(f"Using column '{header[col_idx]}' (index {col_idx}) for splitting.")

            writers = {}
            file_handles = {}
            row_count = 0

            for row in reader:
                if not row or len(row) <= col_idx:
                    continue
                domain_val = row[col_idx].strip()
                safe_domain = "".join(c if c.isalnum() or c in ('_', '-', '.') else '_' for c in domain_val)
                if not safe_domain:
                    safe_domain = "unknown"

                if safe_domain not in file_handles:
                    out_path = os.path.join(args.output_dir, f"{safe_domain}.tsv")
                    out_f = open(out_path, mode="w", encoding="utf-8", newline="")
                    writer = csv.writer(out_f, delimiter="\t")
                    writer.writerow(header)
                    file_handles[safe_domain] = out_f
                    writers[safe_domain] = writer

                writers[safe_domain].writerow(row)
                row_count += 1

            for f_handle in file_handles.values():
                f_handle.close()

            print(f"Successfully processed {row_count} rows into {len(file_handles)} unique files in '{args.output_dir}'.")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()