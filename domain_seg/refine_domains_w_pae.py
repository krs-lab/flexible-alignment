"""
To refine merizo domain segmentation outputs with PAE 

Usage: 
        python3 refine_domains_6.py
        python3 refine_domains_6.py /path/to/pdb/pae/input_dir
"""

#!/usr/bin/env python3
import os
import sys
import json
import glob
import numpy as np
import pandas as pd

# ---
# loads PAE matrix and converts into NumPy array 
# ---
def load_pae_matrix(json_path):
    """Loads and normalizes PAE matrix from AlphaFold/ColabFold JSON files."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    if isinstance(data, list):
        data = data[0]
    if "predicted_aligned_error" in data:
        return np.array(data["predicted_aligned_error"])
    elif "pae" in data:
        return np.array(data["pae"])
    else:
        raise KeyError(f"No PAE matrix found in {json_path}")

# ---
# reads Merizo output files (.domains)
# ---
def parse_merizo_domains(domain_file_path):
    """Parses Merizo .domains files into 1-based residue index lists."""
    domains = []
    with open(domain_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('name'):
                continue
            parts = line.split()
            for part in parts:
                if '-' in part and not part.startswith('v'):
                    dom_blocks = part.split(',')
                    for block in dom_blocks:
                        residues = []
                        sub_blocks = block.split('_')
                        for sb in sub_blocks:
                            if '-' in sb:
                                try:
                                    s, e = map(int, sb.split('-'))
                                    residues.extend(range(s, e + 1))
                                except ValueError:
                                    continue
                        if residues:
                            domains.append(sorted(residues))
    return domains

# ---
# calculate avg PAE between 2 protein domains 
# ---
def calculate_inter_domain_pae(pae_matrix, dom_a, dom_b):
    """Calculates symmetric mean inter-domain PAE between two sets of residue numbers."""
    idx_a = [r - 1 for r in dom_a if 0 <= r - 1 < pae_matrix.shape[0]]
    idx_b = [r - 1 for r in dom_b if 0 <= r - 1 < pae_matrix.shape[1]]
    if not idx_a or not idx_b:
        return 0.0
    pae_a2b = pae_matrix[np.ix_(idx_a, idx_b)].mean()
    pae_b2a = pae_matrix[np.ix_(idx_b, idx_a)].mean()
    return (pae_a2b + pae_b2a) / 2.0

def main():
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = os.path.dirname(os.path.abspath(__file__))

    domain_files = glob.glob(os.path.join(target_dir, "*.domains"))
    
    MERGE_PAE_CUTOFF = 2.7 # calibrated based on only 4 control proteins
    
    summary_results = []
    json_export = {}
    
    for dom_path in sorted(domain_files):
        base_name = os.path.basename(dom_path).split('_merizo')[0]
        json_path = os.path.join(target_dir, f"{base_name}_pae.json")
        
        if not os.path.exists(json_path):
            continue
            
        pae_matrix = load_pae_matrix(json_path)
        merizo_doms = parse_merizo_domains(dom_path)
        initial_count = len(merizo_doms)
        
        action = "Unchanged"
        final_doms = merizo_doms
        
        # Apply PAE Merge Filter if Merizo predicted > 1 domain
        if initial_count > 1:
            i = 0
            merged_doms = [list(d) for d in merizo_doms]
            while i < len(merged_doms):
                j = i + 1
                while j < len(merged_doms):
                    mean_pae = calculate_inter_domain_pae(pae_matrix, merged_doms[i], merged_doms[j])
                    if mean_pae < MERGE_PAE_CUTOFF:
                        merged_doms[i] = sorted(merged_doms[i] + merged_doms[j])
                        merged_doms.pop(j)
                        action = "Merged (PAE)"
                    else:
                        j += 1
                i += 1
            final_doms = merged_doms
        
        # Save boundary assignments to dictionary
        json_export[base_name] = final_doms
        
        summary_results.append({
            "protein_id": base_name,
            "merizo_raw": initial_count,
            "pae_consensus": len(final_doms),
            "action_taken": action
        })
        
    # Write full boundary mappings to JSON
    output_json_path = "refined_domains.json"
    with open(output_json_path, 'w') as f:
        json.dump(json_export, f, indent=2)
        
    # Console Summary Table
    df = pd.DataFrame(summary_results)
    print("\n" + "="*65)
    print(df.to_string(index=False))
    print("="*65)
    print(f"\n[SUCCESS] Wrote refined domain residue boundaries to {output_json_path}")

if __name__ == "__main__":
    main()
