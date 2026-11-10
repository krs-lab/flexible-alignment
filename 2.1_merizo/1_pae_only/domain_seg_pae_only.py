#!/usr/bin/env python3
"""
Batch PAE-based domain segmentation for AlphaFold/ColabFold rank-1 models.

For every rank-1 PAE JSON found under INPUT_DIR (matching PAE_GLOB), this
script:
  1. Loads the PAE matrix (handles both ColabFold and AlphaFold-DB schemas)
  2. Detects domain boundaries from the cross-block PAE profile using
     PEAK detection (robust to single noisy threshold-crossings)
  3. Merges any domain fragment shorter than MIN_DOMAIN_LEN into its
     neighbor, so you never get spurious 2-residue "domains"
  4. Saves an annotated PAE heatmap PNG per protein
  5. Saves a per-protein domain table as .txt
  6. Appends one row per protein to a combined all-proteins .tsv summary

Requires: numpy, matplotlib, scipy
"""

import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# ==============================================================================
# CONFIG
# ==============================================================================
INPUT_DIR = "."                                # where the rank_001 JSONs live
PAE_GLOB = "*_scores_rank_001_*.json"           # pattern to auto-discover proteins
OUTPUT_DIR = "domain_segmentation_output"

WINDOW = 15                          # flanking residue window for cross-block PAE
PAE_CUTOFF = 12.0                    # minimum cross-block PAE to count as a boundary
MIN_DOMAIN_LEN = 40                  # minimum residues per domain; shorter -> merged
PEAK_MIN_DISTANCE = MIN_DOMAIN_LEN   # minimum spacing between boundary peaks

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ==============================================================================
# 1. DISCOVER INPUT FILES
# ==============================================================================
def discover_pae_files(input_dir, pattern):
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not files:
        raise FileNotFoundError(
            f"No PAE JSON files matched '{pattern}' in '{input_dir}'."
        )
    return files


def protein_id_from_filename(path):
    """Pull the leading protein/UniProt ID off a ColabFold-style filename,
    e.g. 'Q189K5_scores_rank_001_alphafold2_ptm_model_2_seed_000.json' -> 'Q189K5'."""
    base = os.path.basename(path)
    match = re.match(r"^([A-Za-z0-9_.\-]+?)_scores_rank", base)
    return match.group(1) if match else os.path.splitext(base)[0]


# ==============================================================================
# 2. LOAD PAE MATRIX (ColabFold or AlphaFold-DB schema)
# ==============================================================================
def load_pae(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)

    pae_data = data[0] if isinstance(data, list) else data

    if "predicted_aligned_error" in pae_data:
        pae = np.array(pae_data["predicted_aligned_error"], dtype=float)
    elif "pae" in pae_data:
        pae = np.array(pae_data["pae"], dtype=float)
    else:
        raise KeyError(f"No PAE matrix found in {json_path}")

    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE matrix in {json_path} is not square: {pae.shape}")

    return pae


# ==============================================================================
# 3. DOMAIN BOUNDARY DETECTION
# ==============================================================================
def cross_block_profile(pae, window):
    """Mean PAE between the `window` residues just before i and the `window`
    residues just after i. High values mean i sits between two poorly
    cross-aligned (likely different) domains."""
    num_res = pae.shape[0]
    profile = np.zeros(num_res)
    for i in range(window, num_res - window):
        profile[i] = np.mean(pae[i - window:i, i:i + window])
    return profile


def detect_boundaries(profile, cutoff, min_distance):
    """Peak-based boundary detection instead of first-threshold-crossing:
    finds LOCAL MAXIMA of the cross-block PAE profile above `cutoff`, at
    least `min_distance` residues apart. Much less sensitive to noisy,
    single-residue crossings than a greedy left-to-right scan, and it
    places the cut at the true center of the high-PAE band rather than
    its leading edge."""
    peaks, _ = find_peaks(profile, height=cutoff, distance=min_distance)
    return sorted(peaks.tolist())


def build_domains(cuts, num_res, min_domain_len):
    """Turn boundary positions into 1-indexed inclusive (start, end) domain
    spans, then merge any span shorter than min_domain_len into a neighbor
    so noise never produces a tiny spurious domain."""
    bounds = [0] + cuts + [num_res]
    spans = [(bounds[i] + 1, bounds[i + 1]) for i in range(len(bounds) - 1)]

    merged = []
    for span in spans:
        length = span[1] - span[0] + 1
        if merged and length < min_domain_len:
            prev_start, _ = merged[-1]
            merged[-1] = (prev_start, span[1])
        else:
            merged.append(span)

    # a short FIRST domain has no left neighbor to merge into -> fold it right
    if len(merged) > 1 and (merged[0][1] - merged[0][0] + 1) < min_domain_len:
        merged[1] = (merged[0][0], merged[1][1])
        merged.pop(0)

    return merged


# ==============================================================================
# 4. PLOTTING
# ==============================================================================
def plot_pae_with_domains(pae, domains, protein_id, out_png):
    num_res = pae.shape[0]
    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    cax = ax.imshow(pae, cmap="bwr_r", vmin=0, vmax=30, origin="upper")

    major_step = max(100, (num_res // 10 // 100) * 100 or 100)
    minor_step = max(20, num_res // 40)
    ax.set_xticks(np.arange(0, num_res, major_step))
    ax.set_yticks(np.arange(0, num_res, major_step))
    ax.set_xticks(np.arange(0, num_res, minor_step), minor=True)
    ax.set_yticks(np.arange(0, num_res, minor_step), minor=True)
    ax.grid(which="minor", color="black", linestyle=":", linewidth=0.5, alpha=0.3)
    ax.grid(which="major", color="black", linestyle="--", linewidth=0.8, alpha=0.6)

    for start, end in domains[:-1]:
        ax.axvline(x=end, color="yellow", linestyle="-.", linewidth=1.2, alpha=0.85)
        ax.axhline(y=end, color="yellow", linestyle="-.", linewidth=1.2, alpha=0.85)

    ax.set_title(f"{protein_id} — Rank 1 PAE Domain Segmentation ({num_res} aa)",
                 fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Residue Position", fontsize=12)
    ax.set_ylabel("Residue Position", fontsize=12)

    cbar = fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Predicted Aligned Error (Å)", fontsize=11)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close(fig)


# ==============================================================================
# 5. TEXT OUTPUT (per protein + combined summary)
# ==============================================================================
def write_domain_table(protein_id, domains, num_res, out_txt):
    with open(out_txt, "w") as f:
        f.write(f"Protein: {protein_id}\n")
        f.write(f"Length: {num_res} aa\n")
        f.write(f"Domains detected: {len(domains)}\n")
        f.write("=" * 50 + "\n")
        for idx, (start, end) in enumerate(domains, 1):
            length = end - start + 1
            f.write(f"Domain {idx:02d}: Residues {start:4d}-{end:4d}  (Length: {length:3d} aa)\n")
        f.write("=" * 50 + "\n")


def append_summary_row(summary_path, protein_id, num_res, domains, write_header):
    with open(summary_path, "a") as f:
        if write_header:
            f.write("protein_id\tlength_aa\tn_domains\tdomain_boundaries\n")
        boundary_str = ";".join(f"{s}-{e}" for s, e in domains)
        f.write(f"{protein_id}\t{num_res}\t{len(domains)}\t{boundary_str}\n")


# ==============================================================================
# 6. MAIN
# ==============================================================================
def process_one(json_path, protein_id):
    print(f"[*] {protein_id}: loading {os.path.basename(json_path)}")
    pae = load_pae(json_path)
    num_res = pae.shape[0]

    profile = cross_block_profile(pae, WINDOW)
    cuts = detect_boundaries(profile, PAE_CUTOFF, PEAK_MIN_DISTANCE)
    domains = build_domains(cuts, num_res, MIN_DOMAIN_LEN)

    out_png = os.path.join(OUTPUT_DIR, f"{protein_id}_Rank1_PAE_Domains.png")
    out_txt = os.path.join(OUTPUT_DIR, f"{protein_id}_domains.txt")

    plot_pae_with_domains(pae, domains, protein_id, out_png)
    write_domain_table(protein_id, domains, num_res, out_txt)

    print(f"    -> {len(domains)} domain(s); PNG: {out_png}; table: {out_txt}")
    return num_res, domains


def main():
    files = discover_pae_files(INPUT_DIR, PAE_GLOB)
    summary_path = os.path.join(OUTPUT_DIR, "all_proteins_domain_summary.tsv")
    if os.path.exists(summary_path):
        os.remove(summary_path)

    write_header = True
    for path in files:
        protein_id = protein_id_from_filename(path)
        try:
            num_res, domains = process_one(path, protein_id)
            append_summary_row(summary_path, protein_id, num_res, domains, write_header)
            write_header = False
        except Exception as e:
            print(f"[!] Skipping {protein_id}: {e}", file=sys.stderr)

    print(f"\nDone. Combined summary: {summary_path}")


if __name__ == "__main__":
    main()
