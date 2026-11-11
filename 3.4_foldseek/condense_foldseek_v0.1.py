#!/usr/bin/env python3
"""
condense_foldseek.py

Pool Foldseek hits from two search modes (TM-align and 3Di) into one 2D scatter
plot per query domain, highlighting hits that pass BOTH filters and (optionally)
hits previously validated in a CSV.

Significance filters (applied here, per query):
  TM-align file : qtmscore >= --tm-cutoff                      (default 0.4)
  3Di file      : bits >= --min-bits  AND
                  modified z-score of bits >= --z-cutoff       (defaults 13, 2)
                  modified z = 0.6745 * (bits - median) / MAD, computed over ALL
                  targets of that query in the 3Di file (median/MAD based).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def parse_args():
    p = argparse.ArgumentParser(
        prog="condense_foldseek.py",
        description=("Pool Foldseek TM-align and 3Di results into a 2D scatter plot "
                     "(x = TM-align qTM-score, y = 3Di bit score) per query domain. "
                     "Common hits (significant in both modes) and previously "
                     "validated hits are highlighted."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i-tma", "--input-tma", dest="input_tma", required=True,
                   metavar="TSV", help="Foldseek TM-align TSV with header "
                   "(needs columns: query, target, qtmscore)")
    p.add_argument("-i-3di", "--input-3di", dest="input_3di", required=True,
                   metavar="TSV", help="Foldseek 3Di TSV with header "
                   "(needs columns: query, target, bits)")
    p.add_argument("-dpi", "--dpi", type=int, default=300,
                   help="Resolution of the output PNG")
    p.add_argument("-v-list", "--validated-list", "-c", "--highlight-csv",
                   dest="highlight_csv", default=None, metavar="CSV",
                   help="Optional CSV of previously identified hits to highlight; "
                   "targets are taken from its 'highlight' column "
                   "(needs columns: Query, highlight)")
    p.add_argument("-q", "--query", nargs="+", default=None,
                   help="Only process these query names (default: all queries "
                   "found in both files)")
    p.add_argument("-o", "--outdir", default=".", help="Output directory")
    p.add_argument("-qtm", "--qtm-cutoff", "--tm-cutoff", dest="tm_cutoff",
                   type=float, default=0.4, metavar="FLOAT",
                   help="TM-align filter: keep hits with query TM-score "
                   "(qtmscore) >= this value")
    p.add_argument("-bit-raw", "--bit-raw", "--min-bits", dest="min_bits",
                   type=float, default=13, metavar="FLOAT",
                   help="3Di filter: keep hits with raw bit score >= this value")
    p.add_argument("-z-bit", "--z-bit", "--z-cutoff", dest="z_cutoff",
                   type=float, default=2.0, metavar="FLOAT",
                   help="3Di filter: keep hits with median-based modified "
                   "z-score of the raw bit score >= this value (applied "
                   "together with -bit-raw)")
    p.add_argument("--linear-y", action="store_true",
                   help="Use a linear instead of log y-axis")
    p.add_argument("--no-labels", action="store_true",
                   help="Do not label points")
    return p.parse_args()


def read_tsv(path, needed):
    # skip the huge alignment-string columns for speed/memory
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in set(needed))
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {path} is missing column(s): {', '.join(missing)}")
    return df


def modified_z(x):
    """Median/MAD based modified z-score (Iglewicz & Hoaglin)."""
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:  # fall back to mean absolute deviation
        mad = np.mean(np.abs(x - med)) / 1.253314
    if mad == 0:
        return np.zeros_like(x, dtype=float)
    return 0.6745 * (x - med) / mad


def safe_name(s):
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in s)


def process_query(q, tma, d3, validated, args):
    t = tma[tma["query"] == q][["target", "qtmscore"]].copy()
    d = d3[d3["query"] == q][["target", "bits"]].copy()
    # keep the best score if a target appears more than once
    t = t.groupby("target", as_index=False)["qtmscore"].max()
    d = d.groupby("target", as_index=False)["bits"].max()

    d["zscore"] = modified_z(d["bits"].to_numpy(float))
    t["sig_tma"] = t["qtmscore"] >= args.tm_cutoff
    d["sig_3di"] = (d["bits"] >= args.min_bits) & (d["zscore"] >= args.z_cutoff)

    m = t.merge(d, on="target", how="outer")
    m["sig_tma"] = m["sig_tma"].fillna(False).astype(bool)
    m["sig_3di"] = m["sig_3di"].fillna(False).astype(bool)
    m["common"] = m["sig_tma"] & m["sig_3di"]
    m["validated"] = m["target"].isin(validated)
    m.insert(0, "query", q)

    n_missing = int(m[["qtmscore", "bits"]].isna().any(axis=1).sum())
    plot_df = m.dropna(subset=["qtmscore", "bits"])

    # ---------------- plot ----------------
    fig, ax = plt.subplots(figsize=(9, 7.5))
    bg = plot_df[~plot_df.sig_tma & ~plot_df.sig_3di]
    only_t = plot_df[plot_df.sig_tma & ~plot_df.sig_3di]
    only_d = plot_df[~plot_df.sig_tma & plot_df.sig_3di]
    both = plot_df[plot_df.common]

    ax.scatter(bg.qtmscore, bg.bits, s=14, c="#b5b9c0", alpha=.6, lw=0, zorder=2)
    ax.scatter(only_t.qtmscore, only_t.bits, s=28, c="#e9a23b", alpha=.85,
               edgecolor="k", lw=.3, zorder=3)
    ax.scatter(only_d.qtmscore, only_d.bits, s=28, c="#3aa6a0", alpha=.85,
               edgecolor="k", lw=.3, zorder=3)
    ax.scatter(both.qtmscore, both.bits, s=60, c="#d62f2f", edgecolor="k",
               lw=.6, zorder=4)

    val = plot_df[plot_df.validated]
    ax.scatter(val.qtmscore, val.bits, s=190, facecolors="none",
               edgecolor="#2563eb", lw=1.8, zorder=5)

    ax.axvline(args.tm_cutoff, color="#8b1a1a", ls="--", lw=1, zorder=1)
    ax.axhline(args.min_bits, color="#8b1a1a", ls=":", lw=1, zorder=1)
    # bit score equivalent of the z-cutoff
    zb = d.loc[d.zscore >= args.z_cutoff, "bits"]
    if len(zb):
        ax.axhline(zb.min(), color="#8b1a1a", ls="--", lw=1, zorder=1)

    if not args.linear_y:
        ax.set_yscale("log")
    ax.set_xlabel(f"TM-align hit strength: query TM-score (cutoff {args.tm_cutoff})",
                  fontweight="bold")
    ax.set_ylabel("3Di hit strength: raw bit score "
                  f"(bits>={args.min_bits:g}, mod. z>={args.z_cutoff:g})",
                  fontweight="bold")
    ax.set_title(f"Foldseek TM-align vs 3Di pooled hits\nQuery: {q}",
                 fontweight="bold")
    ax.grid(alpha=.2)

    if not args.no_labels:
        try:
            from adjustText import adjust_text
        except ImportError:
            adjust_text = None
        # without adjustText, labels of crowded common hits would overlap,
        # so only previously identified hits are labelled
        lab = plot_df[plot_df.common | plot_df.validated] if adjust_text \
            else plot_df[plot_df.validated]
        texts = []
        for _, r in lab.iterrows():
            col = "#2563eb" if r.validated else "#7a1010"
            texts.append(ax.text(r.qtmscore, r.bits, r.target, fontsize=6.5,
                                 color=col, zorder=6))
        if adjust_text and texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="gray", lw=.4))
        else:
            for i, tx in enumerate(texts):
                tx.set_position((tx.get_position()[0] + 0.006, tx.get_position()[1]))
                tx.set_va("bottom" if i % 2 else "top")

    handles = [
        Line2D([], [], marker="o", ls="", c="#b5b9c0", label=f"Not significant ({len(bg)})"),
        Line2D([], [], marker="o", ls="", c="#e9a23b", mec="k", label=f"TM-align only ({len(only_t)})"),
        Line2D([], [], marker="o", ls="", c="#3aa6a0", mec="k", label=f"3Di only ({len(only_d)})"),
        Line2D([], [], marker="o", ls="", c="#d62f2f", mec="k", ms=8, label=f"Common to both ({len(both)})"),
        Line2D([], [], marker="o", ls="", mfc="none", mec="#2563eb", mew=1.8, ms=12,
               label=f"Previously identified ({len(val)})"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, framealpha=.95)
    fig.tight_layout()

    base = os.path.join(args.outdir, f"{safe_name(q)}_foldseek_condensed")
    fig.savefig(base + ".png", dpi=args.dpi)
    plt.close(fig)

    out = m[m.sig_tma | m.sig_3di | m.validated].sort_values(
        ["common", "qtmscore"], ascending=[False, False])
    out.to_csv(base + ".csv", index=False)

    print(f"[{q}] TM-align hits: {int(m.sig_tma.sum())} | 3Di hits: {int(m.sig_3di.sum())} "
          f"| common: {len(both)} | validated (in plot): {len(val)}")
    if n_missing:
        print(f"    note: {n_missing} target(s) absent from one file were not plotted")
    miss_val = set(validated) - set(plot_df.target)
    if miss_val:
        print(f"    note: validated target(s) not found in the TSVs: {sorted(miss_val)}")
    print(f"    -> {base}.png / {base}.csv")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tma = read_tsv(args.input_tma, ["query", "target", "qtmscore"])
    d3 = read_tsv(args.input_3di, ["query", "target", "bits"])

    hl = None
    if args.highlight_csv:
        hl = pd.read_csv(args.highlight_csv)
        for c in ("Query", "highlight"):
            if c not in hl.columns:
                sys.exit(f"ERROR: {args.highlight_csv} is missing column '{c}'")

    queries = [q for q in tma["query"].unique() if q in set(d3["query"])]
    if args.query:
        queries = [q for q in queries if q in args.query]
    if not queries:
        sys.exit("ERROR: no matching query names found in both files")

    for q in queries:
        validated = set()
        if hl is not None:
            validated = set(hl.loc[hl["Query"] == q, "highlight"].dropna())
        process_query(q, tma, d3, validated, args)


if __name__ == "__main__":
    main()
