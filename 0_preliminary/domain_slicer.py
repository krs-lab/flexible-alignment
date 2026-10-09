import os
import csv
import glob
import re
import argparse
from Bio.PDB import PDBParser, PDBIO, Select

# ==============================================================================
# CLASS: DomainSelect
# Purpose: Filters PDB residues based on a list of residue range tuples.
# ==============================================================================
class DomainSelect(Select):
    def __init__(self, ranges):
        # Store residue range tuples, e.g., [(1, 20), (244, 318)]
        self.ranges = ranges

    def accept_residue(self, residue):
        # Filter out heteroatoms, water molecules, and ligands
        if residue.id[0] != " ":
            return 0
            
        # Extract residue sequence number
        res_num = residue.id[1]
        
        # Keep residue if it falls into ANY of the specified segments
        for start, end in self.ranges:
            if start <= res_num <= end:
                return 1
                
        return 0


# ==============================================================================
# FUNCTION: parse_range_string
# Purpose: Extracts numerical start-end pairs from continuous or discontinuous 
# Merizo range strings (e.g., "21-107", "1-20_244-318", or "10-50,100-150").
# ==============================================================================
def parse_range_string(range_str):
    ranges = []
    # Find all pairs of numbers separated by a hyphen
    # Match pattern: one or more digits, a hyphen, one or more digits
    matches = re.findall(r'(\d+)-(\d+)', range_str)
    
    for start, end in matches:
        ranges.append((int(start), int(end)))
        
    return ranges


# ==============================================================================
# FUNCTION: process_merizo_domains
# Purpose: Parses Merizo domains file and extracts individual PDB domain files.
# ==============================================================================
def process_merizo_domains(domains_file, pdb_dir, out_dir, min_length):
    parser = PDBParser(QUIET=True)
    
    filename = os.path.basename(domains_file)
    accession_id = filename.replace("_merizo_v2.domains", "").replace(".domains", "")
    
    possible_names = [f"{accession_id}.pdb", f"{accession_id}_merizo_v2.pdb"]
    pdb_path = None
    
    for search_dir in [pdb_dir, os.getcwd()]:
        for name in possible_names:
            candidate = os.path.join(search_dir, name)
            if os.path.exists(candidate):
                pdb_path = candidate
                break
        if pdb_path:
            break
            
    if not pdb_path:
        print(f"[SKIP] Source PDB for {filename} not found.")
        return

    structure = parser.get_structure(accession_id, pdb_path)
    
    with open(domains_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('='):
                continue
            
            parts = line.split()
            if len(parts) >= 7:
                dom_num = parts[1]      # e.g., "1"
                length = int(parts[2])  # e.g., 95
                conf = parts[3]         # e.g., "0.806"
                range_str = parts[6]    # e.g., "1-20_244-318"
                
                # Skip domains that fall below the minimum residue length threshold
                if length < min_length:
                    print(f"[SKIP] Domain {dom_num} of {accession_id} length ({length} aa) < threshold ({min_length} aa)")
                    continue
                
                ranges = parse_range_string(range_str)
                if not ranges:
                    continue
                
                io = PDBIO()
                io.set_structure(structure)
                
                # Output file name schema: e.g., "Q5ZSI8_dom1_1-20_244-318.pdb"
                out_name = f"{accession_id}_dom{dom_num}_{range_str}.pdb"
                out_path = os.path.join(out_dir, out_name)
                
                io.save(out_path, DomainSelect(ranges))
                print(f"[SAVED] {out_name} | Ranges: {ranges} | Length: {length} aa | Conf: {conf}")


# ==============================================================================
# FUNCTION: find_column
# Purpose: Finds a CSV column by trying several accepted names (case-insensitive).
# ==============================================================================
def find_column(fieldnames, candidates):
    lookup = {name.strip().lower(): name for name in fieldnames if name}
    for cand in candidates:
        if cand.lower() in lookup:
            return lookup[cand.lower()]
    return None


# ==============================================================================
# FUNCTION: process_curated_domains
# Purpose: Reads a CSV of curated domain start/end coordinates, matches each row
# to a PDB file in the input directory via the locus ID (case-insensitive), and
# saves the sliced domain, e.g., "gfp_14-227.pdb".
# ==============================================================================
def process_curated_domains(csv_file, pdb_dir, out_dir, min_length=None):
    parser = PDBParser(QUIET=True)

    # Index PDB files in the input directory by lowercase stem (GFP.pdb -> "gfp")
    pdb_index = {}
    for path in glob.glob(os.path.join(pdb_dir, "*.pdb")):
        stem = os.path.splitext(os.path.basename(path))[0]
        pdb_index[stem.lower()] = path

    matched_loci = set()

    with open(csv_file, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []

        # Accepted spellings for the locus column (matching is case-insensitive,
        # so locus / Locus / LOCUS all work, as do the _id variants)
        locus_col = find_column(fields, ["locus_id", "locus", "locusid", "locus id"])
        start_col = find_column(fields, ["start_res", "start", "start_residue"])
        end_col = find_column(fields, ["end_res", "end", "end_residue"])

        if not locus_col or not start_col or not end_col:
            print(f"[ERROR] CSV must contain locus (locus/locus_id), start (start_res) and end (end_res) columns.")
            print(f"        Found columns: {fields}")
            return

        for row in reader:
            locus = (row.get(locus_col) or "").strip()
            if not locus:
                continue

            try:
                start = int(str(row[start_col]).strip())
                end = int(str(row[end_col]).strip())
            except (ValueError, TypeError):
                print(f"[SKIP] {locus}: invalid start/end coordinates ({row.get(start_col)}, {row.get(end_col)})")
                continue

            if start > end:
                start, end = end, start

            length = end - start + 1
            if min_length is not None and length < min_length:
                print(f"[SKIP] {locus} {start}-{end} length ({length} aa) < threshold ({min_length} aa)")
                continue

            pdb_path = pdb_index.get(locus.lower())
            if not pdb_path:
                print(f"[SKIP] Source PDB for locus '{locus}' not found in {pdb_dir}")
                continue

            matched_loci.add(locus.lower())
            structure = parser.get_structure(locus, pdb_path)

            # Check requested range against the actual residues in the PDB
            res_nums = [r.id[1] for r in structure.get_residues() if r.id[0] == " "]
            if not res_nums:
                print(f"[SKIP] {locus}: no standard residues found in {os.path.basename(pdb_path)}")
                continue
            prot_min, prot_max = min(res_nums), max(res_nums)
            if start < prot_min or end > prot_max:
                print(f"[SKIP] {locus}: requested range {start}-{end} is out of range "
                      f"(protein residues {prot_min}-{prot_max}, {len(res_nums)} aa)")
                continue

            io = PDBIO()
            io.set_structure(structure)

            # Output file name schema: e.g., "gfp_14-227.pdb"
            out_name = f"{locus}_{start}-{end}.pdb"
            out_path = os.path.join(out_dir, out_name)

            io.save(out_path, DomainSelect([(start, end)]))
            print(f"[SAVED] {out_name} | Range: {start}-{end} | Length: {length} aa")

    # Report PDBs in the input dir that had no matching CSV entry
    for stem, path in sorted(pdb_index.items()):
        if stem not in matched_loci:
            print(f"[NOTE] {os.path.basename(path)} had no sliced domain (no matching/valid CSV entry).")


if __name__ == "__main__":
    # Parse CLI arguments
    cli_parser = argparse.ArgumentParser(
        description="Slice PDB files by Merizo domain predictions or curated coordinates (CSV).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 slicer.py -i test_slicing_dir -c pipBinder_domain_info_Bolt.csv\n"
            "  python3 slicer.py --input test_slicing_dir --coordinate pipBinder_domain_info_Bolt.csv\n"
            "  python3 slicer.py -i test_slicing_dir -c pipBinder_domain_info_Bolt.csv -m 50\n"
        ),
    )
    cli_parser.add_argument("-i", "--input", type=str, default=None, help="Directory containing the source PDB files.")
    cli_parser.add_argument("-c", "--coordinate", type=str, default=None, help="CSV file with curated domain start/end coordinates (curated mode).")
    cli_parser.add_argument("-m", "--min-length", type=int, default=None, help="Minimum amino acid length threshold for domain slicing.")
    args = cli_parser.parse_args()

    CURRENT_DIR = os.getcwd()

    # Choose segmentation mode: providing a CSV selects curated mode, otherwise ask
    if args.coordinate:
        mode = "curated"
    else:
        while True:
            choice = input("Segment domains by (1) Merizo or (2) curated coordinates CSV? [1/2]: ").strip().lower()
            if choice in ("1", "merizo", "m"):
                mode = "merizo"
                break
            elif choice in ("2", "curated", "c"):
                mode = "curated"
                break
            print("Invalid input. Please enter 1 (Merizo) or 2 (curated).")

    # ----------------------------------------------------------------------
    # CURATED MODE: slice using start/end coordinates from a CSV file
    # ----------------------------------------------------------------------
    if mode == "curated":
        csv_path = args.coordinate
        if not csv_path:
            csv_path = input("Enter path to coordinates CSV file: ").strip()
        if not os.path.isfile(csv_path):
            raise SystemExit(f"[ERROR] CSV file not found: {csv_path}")

        pdb_dir = args.input
        if not pdb_dir:
            pdb_dir = input("Enter path to directory containing PDB files: ").strip()
        if not os.path.isdir(pdb_dir):
            raise SystemExit(f"[ERROR] Input directory not found: {pdb_dir}")

        OUTPUT_DIR = os.path.join(CURRENT_DIR, "domains_sliced")
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        print(f"Mode: curated coordinates")
        print(f"Input PDB Directory: {pdb_dir}")
        print(f"Coordinates CSV: {csv_path}")
        if args.min_length is not None:
            print(f"Minimum Length Threshold: {args.min_length} aa")
        print()

        process_curated_domains(csv_path, pdb_dir, OUTPUT_DIR, args.min_length)

        print(f"\nDone! Sliced PDBs saved to: {OUTPUT_DIR}")

    # ----------------------------------------------------------------------
    # MERIZO MODE: slice using .domains files in the current directory
    # ----------------------------------------------------------------------
    else:
        # Prompt interactively if argument was not provided via command line
        min_length = args.min_length
        if min_length is None:
            while True:
                try:
                    user_input = input("Enter minimum amino acid length threshold (e.g., 87): ").strip()
                    min_length = int(user_input)
                    break
                except ValueError:
                    print("Invalid input. Please enter a valid integer.")

        # Update this path if needed (can be overridden with -i/--input)
        PDB_SOURCE_DIR = args.input or "/home/krslab/projects/visnu/main/exp/8_3D_homology_search/pdb_pae_query-pr-set"
        
        OUTPUT_DIR = os.path.join(CURRENT_DIR, "domain_pdbs")
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        
        domain_files = sorted(glob.glob(os.path.join(CURRENT_DIR, "*.domains")))
        
        print(f"Working Directory: {CURRENT_DIR}")
        print(f"Minimum Length Threshold: {min_length} aa")
        print(f"Found {len(domain_files)} .domains file(s). Processing...\n")
        
        for d_file in domain_files:
            process_merizo_domains(d_file, PDB_SOURCE_DIR, OUTPUT_DIR, min_length)
            
        print(f"\nDone! Sliced PDBs saved to: {OUTPUT_DIR}")
