#!/usr/bin/env python
"""Rank plot of models by a metric, highlighting the flagged outliers.

Reads metrics.tsv (which already carries the boolean `outlier` column from
gather_metrics) and draws every model ranked by `metric`, with the Tukey fence and
the flagged (lower-tail) models highlighted. No table is written here; the flagged
models live in the `outlier` column of metrics.tsv.
"""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--metrics", required=True, help="metrics.tsv.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-m", "--metric", default="count_pearson",
        help="Metric to rank on (default count_pearson).")
    parser.add_argument("-k", "--iqr_mult", type=float, default=1.5,
        help="Tukey fence multiplier for the drawn fence line (default 1.5).")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy as np
    import pandas as pd
    from _style import apply_style, despine, save_figure
    apply_style()
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.metrics, sep="\t").dropna(subset=[args.metric])
    df = df.sort_values(args.metric, ascending=True).reset_index(drop=True)

    m = df[args.metric]
    q1, q3 = m.quantile(0.25), m.quantile(0.75)
    fence = q1 - args.iqr_mult * (q3 - q1)
    is_out = df["outlier"].to_numpy()
    y = m.to_numpy()

    fig, ax = plt.subplots(figsize=(6, 3.6))
    rank = np.arange(len(df))
    ax.scatter(rank[~is_out], y[~is_out], s=5, color="#3b6ea5",
        alpha=0.5, edgecolor="none", label="pass")
    ax.scatter(rank[is_out], y[is_out], s=14, color="#d1495b",
        edgecolor="none", label=f"outlier (< {fence:.2f})")
    ax.axhline(fence, color="#d1495b", lw=1, ls="--")
    ax.set_xlabel(f"models ranked by {args.metric}")
    ax.set_ylabel(args.metric)
    ax.set_title(f"{args.metric}: {int(is_out.sum())} of {len(df)} below Tukey fence")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    despine(ax)
    save_figure(fig, args.output)

    print(f"{int(is_out.sum())} outliers (fence={fence:.3f})")


if __name__ == "__main__":
    main()
