#!/usr/bin/env python
"""Observed vs predicted log-count scatter for one model (count-head QC).

Every point is a test peak: x = observed log-counts, y = the model's predicted
log-counts. Density-shaded and annotated with the Pearson r that count_pearson
summarises, so a weak count head is visible at a glance.
"""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--counts", required=True,
        help="Per-region counts TSV (group, obs_logcount, pred_logcount).")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("--sample", default="", help="Sample name for the title.")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy as np
    import pandas as pd
    from _style import apply_style, save_figure
    apply_style()
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.counts, sep="\t")
    x = df["obs_logcount"].to_numpy(float)
    y = df["pred_logcount"].to_numpy(float)
    r = np.corrcoef(x, y)[0, 1]

    fig, ax = plt.subplots(figsize=(4.2, 4))
    hb = ax.hexbin(x, y, gridsize=45, mincnt=1, cmap="viridis", bins="log")
    fig.colorbar(hb, ax=ax, label="peaks")

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], color="#d1495b", lw=1, ls="--", label="y = x")

    ax.set_xlabel("observed log-counts")
    ax.set_ylabel("predicted log-counts")
    # Wrap the (long) sample name so the title fits the figure width.
    import textwrap
    header = "\n".join(textwrap.wrap(args.sample, 30)) + "\n" if args.sample else ""
    ax.set_title(f"{header}count Pearson r = {r:.3f}  (n={len(x)})", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
