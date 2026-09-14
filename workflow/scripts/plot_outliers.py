#!/usr/bin/env python
"""Flag under-performing models and write the full per-model metrics table.

Outliers are the lower tail of `metric` by the Tukey rule (below Q1 - k*IQR).
Writes a rank plot with the fence + flagged models highlighted, and a TSV of every
model (all metrics + covariates) carrying an `outlier` boolean column.
"""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--metrics", required=True, help="metrics.tsv.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-t", "--table", required=True,
        help="Destination TSV: every model, all metrics, with an `outlier` bool.")
    parser.add_argument("-m", "--metric", default="count_pearson",
        help="Metric to screen on (default count_pearson).")
    parser.add_argument("-k", "--iqr_mult", type=float, default=1.5,
        help="Tukey fence multiplier (default 1.5).")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy as np
    import pandas as pd
    from _style import apply_style
    apply_style()
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.metrics, sep="\t").dropna(subset=[args.metric])
    v = df[args.metric]
    q1, q3 = v.quantile(0.25), v.quantile(0.75)
    fence = q1 - args.iqr_mult * (q3 - q1)
    df = df.sort_values(args.metric, ascending=True).reset_index(drop=True)
    df["outlier"] = df[args.metric] < fence

    # Every model, all metrics/covariates, flagged by `outlier`.
    df.to_csv(args.table, sep="\t", index=False)

    is_out = df["outlier"].to_numpy()
    y = df[args.metric].to_numpy()
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
    fig.savefig(args.output)

    print(f"{int(is_out.sum())} outliers (fence={fence:.3f}); wrote {args.table}")


if __name__ == "__main__":
    main()
