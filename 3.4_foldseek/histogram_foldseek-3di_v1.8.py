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
    description="Generate target domain plots from Foldseek CSV/TSV files (one file per query) using the raw bit score and its robust z-score."
)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Path to input directory containing US-align CSV files (one per query)",
)
parser.add_argument(
    "-bit-raw",
    "--bit-raw",
    dest="bit_raw",
    type=float,
    default=15.0,
    help="Raw bit score must be at or above this value to be a significant hit (default: 15)",
)
parser.add_argument(
    "-z-bit",
    "--z-bit",
    dest="z_bit",
    type=float,
    default=2.0,
    help="Median-based (robust) z-score of the bit score must be above this value (default: 2)",
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
bit_thr = args.bit_raw
z_thr = args.z_bit
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
    bit_col = find_column(usalign_df, ["bit", "bits", "bitscore", "bit_score"])
    query_col = find_column(usalign_df, ["Query"])

    if target_col is None or bit_col is None:
        print(f"Skipping '{csv_file}': missing 'Target' or 'bit'/'bits' column")
        continue

    # Query domain ID from the CSV (fallback: file name)
    if query_col is not None and len(usalign_df) > 0:
        query_domain = usalign_df[query_col].iloc[0]
    else:
        query_domain = os.path.splitext(os.path.basename(csv_file))[0]

    bit_values = pd.to_numeric(usalign_df[bit_col], errors="coerce").fillna(0.0).values
    row_series = pd.Series(bit_values, index=usalign_df[target_col].values)

    # Clean query domain ID (strip .pdb if present)
    clean_query = str(query_domain).replace(".pdb", "")

    # Rank targets by raw bit score (highest first)
    sorted_series = row_series.sort_values(ascending=False, kind="stable")
    # Clean target IDs (strip .pdb if present)
    sorted_targets = [str(col).replace(".pdb", "") for col in sorted_series.index]
    sorted_bits = sorted_series.values
    m = len(sorted_bits)

    if m == 0:
        print(f"Skipping '{csv_file}': no targets")
        continue

    ranks = np.arange(1, m + 1)

    # Median-based (robust) z-score of the bit score across all targets of this query:
    # z = (bits - median) / (1.4826 * MAD); falls back to the standard deviation if MAD is 0
    med = np.median(sorted_bits)
    scale = 1.4826 * np.median(np.abs(sorted_bits - med))
    if scale == 0:
        scale = np.std(sorted_bits)
    z_scores = (sorted_bits - med) / scale if scale > 0 else np.zeros(m)
    sorted_scores = z_scores  # plotted on the y-axis

    # Significant hit: raw bit score >= threshold AND median z-score > threshold
    is_hit = (sorted_bits >= bit_thr) & (sorted_scores > z_thr)
    n_hits = int(np.sum(is_hit))

    # Determine which targets match the validated locus IDs from CSV
    is_validated = np.array([
        t.lower() in validated_locus_set for t in sorted_targets
    ])

    # Break the y-axis only when the top score (e.g. self-alignment) dwarfs the rest
    broken = m > 1 and sorted_scores[0] > 1.5 * max(sorted_scores[1], z_thr)

    # Figure with main plot (left) and summary panel (right)
    fig = plt.figure(figsize=(16, 6), dpi=args.dpi)
    if broken:
        gs = fig.add_gridspec(
            2,
            2,
            width_ratios=[2.2, 1],
            height_ratios=[1, 3],
            wspace=0.03,
            hspace=0.06,
        )
        ax_top = fig.add_subplot(gs[0, 0])
        ax_main = fig.add_subplot(gs[1, 0], sharex=ax_top)
        ax_txt = fig.add_subplot(gs[:, 1])
    else:
        gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.03)
        ax_top = None
        ax_main = fig.add_subplot(gs[0, 0])
        ax_txt = fig.add_subplot(gs[0, 1])
    ax_txt.axis("off")
    plt.sca(ax_main)

    # Bars: lower (main) axes carries the legend labels; top axes shows the outlier
    draw_bars(ax_main, ranks, sorted_scores, is_hit, is_validated, n_hits, label=True)
    if broken:
        draw_bars(ax_top, ranks, sorted_scores, is_hit, is_validated, n_hits)

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
            f"[{i+1:>2}] {sorted_targets[i]:<{max_id_len}} | Bit={sorted_bits[i]:.0f} | Z={sorted_scores[i]:.2f}"
        )
        validated_line_flags.append(bool(is_validated[i]))

    # Reference Threshold Line
    plt.axhline(
        y=z_thr,
        color="#991B1B",
        linestyle="--",
        linewidth=0.9,
        label=f"Z-score Cutoff (z={z_thr:g})",
    )

    # Display All Significant Hits in Summary Panel (right of the plot)
    if top_hits_summary:
        header = f"Significant Hits Summary (Bit >= {bit_thr:g}, Z > {z_thr:g}):"
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

    # Axis scaling and padding (lower axes excludes the outlier when the axis is broken)
    if broken:
        y_max = max(max(sorted_scores[1:]), z_thr)
    else:
        y_max = max(max(sorted_scores), z_thr) if m > 0 else 1
    y_min = min(min(sorted_scores), 0) if m > 0 else 0
    pad = (y_max - y_min) * 0.15
    plt.ylim(y_min - (pad if y_min < 0 else 0), y_max + pad)
    plt.xlim(0, m + 1)

    # Broken-axis styling: top segment shows the outlier, diagonal marks at the break
    if broken:
        ax_top.set_ylim(
            max(sorted_scores[0] * 0.85, y_max + pad), sorted_scores[0] * 1.08
        )
        ax_top.spines["bottom"].set_visible(False)
        ax_main.spines["top"].set_visible(False)
        ax_top.tick_params(bottom=False, labelbottom=False)

        dx, dy_top = 0.006, 0.04
        dy_main = dy_top / 3  # same physical size given the 1:3 height ratio
        brk = dict(color="black", clip_on=False, linewidth=1.0, zorder=10)
        for xc in (0, 1):
            ax_top.plot(
                (xc - dx, xc + dx),
                (-dy_top, dy_top),
                transform=ax_top.transAxes,
                **brk,
            )
            ax_main.plot(
                (xc - dx, xc + dx),
                (1 - dy_main, 1 + dy_main),
                transform=ax_main.transAxes,
                **brk,
            )

    # Axis labels and title
    plt.xlabel("Target Domain Rank Order", fontsize=10, fontweight="bold")
    if broken:
        # Single y-label centered across both axis segments
        pos_top, pos_main = ax_top.get_position(), ax_main.get_position()
        fig.text(
            pos_main.x0 - 0.04,
            (pos_main.y0 + pos_top.y1) / 2,
            "Modified z-score (bit score)",
            rotation=90,
            va="center",
            ha="center",
            fontsize=10,
            fontweight="bold",
        )
    else:
        plt.ylabel("Modified z-score (bit score)", fontsize=10, fontweight="bold")
    title_ax = ax_top if broken else ax_main
    title_ax.set_title(
        f"Target Domain Alignment Z-scores\nQuery Domain: {clean_query}",
        fontsize=11,
        pad=12,
    )

    # Legend (no background grid); sits in the top segment when the axis is broken
    legend_ax = ax_top if broken else ax_main
    handles, labels = ax_main.get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="upper right", fontsize=8, framealpha=0.9)

    # Save plot
    clean_query_file = clean_query.replace("/", "_").replace("\\", "_")
    csv_basename = os.path.splitext(os.path.basename(csv_file))[0]
    save_path = os.path.join(output_dir, f"{csv_basename}_{clean_query_file}.png")

    plt.savefig(save_path, dpi=args.dpi, bbox_inches="tight")
    plt.close()

print(f"Finished. Generated plots inside '{output_dir}'.")
