#!/usr/bin/env python
"""Scatter a performance metric against a per-model covariate (depth / peaks).

Answers "does prediction quality track data depth?" -- e.g. count_pearson vs
n_fragments or vs n_peaks. Annotates the Spearman correlation and a binned
median trend so the relationship is legible under many overplotted points.
"""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--metrics", required=True, help="metrics.tsv.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-x", "--covariate", default="n_fragments",
        help="Column for the x-axis (default n_fragments).")
    parser.add_argument("-y", "--metric", default="count_pearson",
        help="Performance column for the y-axis (default count_pearson).")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy as np
    import pandas as pd
    from _style import apply_style, save_figure
    apply_style()
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.metrics, sep="\t").dropna(subset=[args.covariate, args.metric])
    x, y = df[args.covariate].to_numpy(float), df[args.metric].to_numpy(float)

    # Spearman = Pearson on ranks (avoids a scipy dependency).
    rho = np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1]

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(x, y, s=6, alpha=0.25, edgecolor="none", color="#3b6ea5")

    # Binned median trend on a log-x grid (depth/peaks span orders of magnitude).
    lx = np.log10(x)
    edges = np.linspace(lx.min(), lx.max(), 13)
    idx = np.clip(np.digitize(lx, edges) - 1, 0, len(edges) - 2)
    centers, meds = [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() >= 5:
            centers.append(10 ** ((edges[b] + edges[b + 1]) / 2))
            meds.append(np.median(y[m]))
    ax.plot(centers, meds, color="#d1495b", lw=2, marker="o", ms=3, label="binned median")

    ax.set_xscale("log")
    ax.set_xlabel(args.covariate)
    ax.set_ylabel(args.metric)
    ax.set_title(f"{args.metric} vs {args.covariate}  (Spearman ρ={rho:.2f})")
    ax.legend(frameon=False, fontsize=8)
    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
