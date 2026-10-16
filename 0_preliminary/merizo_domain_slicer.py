import os
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


if __name__ == "__main__":
    # Parse CLI argument for minimum amino acid length threshold
    cli_parser = argparse.ArgumentParser(description="Slice PDB files based on Merizo domain predictions.")
    cli_parser.add_argument("-m", "--min-length", type=int, default=None, help="Minimum amino acid length threshold for domain slicing.")
    args = cli_parser.parse_args()

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

    CURRENT_DIR = os.getcwd()
    
    # Update this path if needed
    PDB_SOURCE_DIR = "/home/krslab/projects/visnu/main/exp/8_3D_homology_search/pdb_pae_query-pr-set"
    
    OUTPUT_DIR = os.path.join(CURRENT_DIR, "domain_pdbs")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    domain_files = sorted(glob.glob(os.path.join(CURRENT_DIR, "*.domains")))
    
    print(f"Working Directory: {CURRENT_DIR}")
    print(f"Minimum Length Threshold: {min_length} aa")
    print(f"Found {len(domain_files)} .domains file(s). Processing...\n")
    
    for d_file in domain_files:
        process_merizo_domains(d_file, PDB_SOURCE_DIR, OUTPUT_DIR, min_length)
        
    print(f"\nDone! Sliced PDBs saved to: {OUTPUT_DIR}")
