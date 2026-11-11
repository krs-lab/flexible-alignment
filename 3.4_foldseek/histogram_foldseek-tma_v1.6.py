import argparse
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Set global matplotlib style for publication-ready figures
plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["font.family"] = "sans-serif"

# Maximum number of hits listed in the summary panel on each plot
MAX_SUMMARY_HITS = 30


def read_table(path):
    """Read a delimited file: tab-separated for .tsv, comma-separated otherwise."""
    sep = "\t" if os.path.splitext(path)[1].lower() == ".tsv" else ","
    return pd.read_csv(path, sep=sep)


def find_column(df, candidates):
    """Return the actual column name matching any candidate (case-insensitive), else None."""
    lookup = {str(c).strip().lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lookup:
            return lookup[cand.lower()]
    return None


def draw_bars(ax, ranks, scores, is_hit, is_validated, n_hits, label=False):
    """Draw background (grey), significant hit (red), and validated/highlighted (bright blue) bars on a given axes."""
    bg_mask = ~is_hit & ~is_validated
    hit_mask = is_hit & ~is_validated
    val_mask = is_validated

    # 1. Background / Insignificant data points: Grey Bars
    if np.any(bg_mask):
        ax.bar(
            ranks[bg_mask],
            scores[bg_mask],
            width=0.8 * 0.6,  # 60% of hit bar width
            color="#9fa5ac",
            edgecolor="none",
            alpha=0.85,
            zorder=2,
            label="Background Targets" if label else None,
        )

    # 2. Significant Hit data points: Red Bars
    if np.any(hit_mask):
        ax.bar(
            ranks[hit_mask],
            scores[hit_mask],
            width=0.8,
            color="#d34447",
            edgecolor="none",
            zorder=3,
            label=f"Significant Hits (Rank 1-{n_hits})" if label else None,
        )

    # 3. Validated / Highlighted Targets from CSV: Bright Blue (#2563EB)
    if np.any(val_mask):
        ax.bar(
            ranks[val_mask],
            scores[val_mask],
            width=0.8,
            color="#2563EB",
            edgecolor="none",
            zorder=4,
            label="Validated Locus IDs" if label else None,
        )


# 1. Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Generate target domain plots from US-align CSV files (one CSV per query) using the query TM-score."
)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Path to input directory containing US-align CSV files (one per query)",
)
parser.add_argument(
    "-tmq",
    "--tm-score-query",
    "-tm-score-query",
    dest="tm_score_query",
    type=float,
    default=0.5,
    help="Minimum query TM-score required to consider a hit valid (default: 0.5)",
)
parser.add_argument(
    "-v",
    "--validated-list",
    type=str,
    default="pipBinder_domains_validated_by_Bolt-2026.csv",
    help="Path to validation CSV file containing highlight column",
)
parser.add_argument(
    "-dpi",
    "--dpi",
    type=int,
    default=1200,
    help="Output figure resolution in DPI (default: 1200)",
)
args = parser.parse_args()

input_dir = args.input
tmq_thr = args.tm_score_query
output_dir = os.path.join(input_dir, "plots")
os.makedirs(output_dir, exist_ok=True)

# Load validation locus IDs if file exists
validated_locus_set = set()
validation_csv = args.validated_list
if os.path.exists(validation_csv):
    try:
        val_df = read_table(validation_csv)
        if "highlight" in val_df.columns:
            validated_locus_set = set(
                val_df["highlight"].dropna().astype(str).str.replace(".pdb", "", regex=False).str.lower().str.strip()
            )
            print(f"Loaded {len(validated_locus_set)} highlight IDs from '{validation_csv}'.")
        else:
            print(f"Warning: 'highlight' column not found in '{validation_csv}'.")
    except Exception as e:
        print(f"Error reading validation CSV: {e}")
else:
    print(f"Notice: Validation CSV '{validation_csv}' not found. Highlighting skipped.")

csv_files = sorted(
    glob.glob(os.path.join(input_dir, "*.csv")) + glob.glob(os.path.join(input_dir, "*.tsv"))
)

if not csv_files:
    print(f"Error: No CSV files found in '{input_dir}'")
    exit(1)

# 2. Process each CSV file (each CSV holds the alignments of one query)
for csv_file in csv_files:
    usalign_df = read_table(csv_file)

    target_col = find_column(usalign_df, ["Target"])
    tmq_col = find_column(usalign_df, ["TM_Score_Query", "qtmscore"])
    query_col = find_column(usalign_df, ["Query"])

    if target_col is None or tmq_col is None:
        print(f"Skipping '{csv_file}': missing 'Target' or 'TM_Score_Query'/'qtmscore' column")
        continue

    # Query domain ID from the CSV (fallback: file name)
    if query_col is not None and len(usalign_df) > 0:
        query_domain = usalign_df[query_col].iloc[0]
    else:
        query_domain = os.path.splitext(os.path.basename(csv_file))[0]

    tmq_values = pd.to_numeric(usalign_df[tmq_col], errors="coerce").fillna(0.0).values
    row_series = pd.Series(tmq_values, index=usalign_df[target_col].values)

    # Clean query domain ID (strip .pdb if present)
    clean_query = str(query_domain).replace(".pdb", "")

    # Rank targets by query TM-score (highest first)
    sorted_series = row_series.sort_values(ascending=False, kind="stable")
    # Clean target IDs (strip .pdb if present)
    sorted_targets = [str(col).replace(".pdb", "") for col in sorted_series.index]
    sorted_scores = sorted_series.values
    m = len(sorted_scores)

    if m == 0:
        print(f"Skipping '{csv_file}': no targets")
        continue

    ranks = np.arange(1, m + 1)

    is_hit = sorted_scores >= tmq_thr
    n_hits = int(np.sum(is_hit))

    # Determine which targets match the validated locus IDs from CSV
    is_validated = np.array([
        t.lower() in validated_locus_set for t in sorted_targets
    ])

    # Figure with main plot (left) and summary panel (right)
    fig = plt.figure(figsize=(16, 6), dpi=args.dpi)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.03)
    ax_main = fig.add_subplot(gs[0, 0])
    ax_txt = fig.add_subplot(gs[0, 1])
    ax_txt.axis("off")
    plt.sca(ax_main)

    draw_bars(ax_main, ranks, sorted_scores, is_hit, is_validated, n_hits, label=True)

    # Tiny locus_id labels below each validated (blue) bar; font height is scaled to the
    # bar width so the labels are only legible when zoomed in
    if np.any(is_validated):
        bar_width_pt = (
            ax_main.get_position().width * fig.get_figwidth() * 72.0 / (m + 1) * 0.8
        )
        locus_fontsize = bar_width_pt * 0.9
        for r, tgt in zip(ranks[is_validated], np.array(sorted_targets)[is_validated]):
            ax_main.annotate(
                tgt.split("_")[0],
                xy=(r, 0),
                xycoords=ax_main.get_xaxis_transform(),
                xytext=(0, -1),
                textcoords="offset points",
                rotation=90,
                fontsize=locus_fontsize,
                color="#2563EB",
                ha="center",
                va="top",
                annotation_clip=False,
                zorder=5,
            )

    # Build summary entries using maximum target ID length for exact column alignment
    top_hits_summary = []
    validated_line_flags = []
    n_listed = min(n_hits, MAX_SUMMARY_HITS)  # summary panel lists only the top hits
    max_id_len = max([len(t) for t in sorted_targets[:n_listed]], default=0)

    for i in range(n_listed):
        top_hits_summary.append(
            f"[{i+1:>2}] {sorted_targets[i]:<{max_id_len}} | TM_Q={sorted_scores[i]:.3f}"
        )
        validated_line_flags.append(bool(is_validated[i]))

    # Reference Threshold Line
    plt.axhline(
        y=tmq_thr,
        color="#991B1B",
        linestyle="--",
        linewidth=0.9,
        label=f"TM-score Cutoff (TM_Q={tmq_thr})",
    )

    # Display All Significant Hits in Summary Panel (right of the plot)
    if top_hits_summary:
        header = f"Significant Hits Summary (min TM_Q = {tmq_thr}):"
        full_text = header + "\n" + "\n".join(top_hits_summary)

        # Box: white (invisible) text sized to the full summary, carrying the original bbox
        base_text = ax_txt.text(
            0.0,
            1.0,
            full_text,
            transform=ax_txt.transAxes,
            fontsize=8,
            fontfamily="monospace",
            color="#FFFFFF",
            va="top",
            ha="left",
            bbox=dict(
                boxstyle="round,pad=0.5",
                facecolor="#FFFFFF",
                edgecolor="#991B1B",
                alpha=0.92,
                lw=0.8,
            ),
            zorder=6,
        )

        # Measure the exact line step from the real summary text (in points)
        renderer = fig.canvas.get_renderer()
        probe_1 = ax_txt.text(0, 1, header, fontsize=8, fontfamily="monospace", va="top")
        line_step_px = (
            base_text.get_window_extent(renderer).height
            - probe_1.get_window_extent(renderer).height
        ) / len(top_hits_summary)
        probe_1.remove()
        line_step_pt = line_step_px * 72.0 / fig.dpi

        # Draw header and each hit line separately on a fixed vertical grid;
        # validated hits are drawn in bright blue (#2563EB)
        summary_lines = [header] + top_hits_summary
        line_colors = ["black"] + [
            "#2563EB" if flag else "black" for flag in validated_line_flags
        ]
        for k, (line, line_color) in enumerate(zip(summary_lines, line_colors)):
            ax_txt.annotate(
                line,
                xy=(0.0, 1.0),
                xycoords="axes fraction",
                xytext=(0, -k * line_step_pt),
                textcoords="offset points",
                fontsize=8,
                fontfamily="monospace",
                color=line_color,
                va="top",
                ha="left",
                zorder=7,
            )

    # Axis scaling and padding (TM-score is bounded, so a single linear axis is enough)
    y_max = max(max(sorted_scores), tmq_thr)
    pad = y_max * 0.15
    plt.ylim(0, y_max + pad)
    plt.xlim(0, m + 1)

    # Axis labels and title
    plt.xlabel("Target Domain Rank Order", fontsize=10, fontweight="bold")
    plt.ylabel("Query TM-score", fontsize=10, fontweight="bold")
    ax_main.set_title(
        f"Target Domain Alignment TM-scores\nQuery Domain: {clean_query}",
        fontsize=11,
        pad=12,
    )

    # Legend
    handles, labels = ax_main.get_legend_handles_labels()
    ax_main.legend(handles, labels, loc="upper right", fontsize=8, framealpha=0.9)

    # Save plot
    clean_query_file = clean_query.replace("/", "_").replace("\\", "_")
    csv_basename = os.path.splitext(os.path.basename(csv_file))[0]
    save_path = os.path.join(output_dir, f"{csv_basename}_{clean_query_file}.png")

    plt.savefig(save_path, dpi=args.dpi, bbox_inches="tight")
    plt.close()

print(f"Finished. Generated plots inside '{output_dir}'.")
