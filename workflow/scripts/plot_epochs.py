#!/usr/bin/env python
"""Plot a model's training/validation metrics over epochs."""

import argparse

# Per-epoch columns written by cherimoya's training logger (<sample>.log).
METRICS = [
    "Training MNLL",
    "Training Count MSE",
    "Validation MNLL",
    "Validation Profile Pearson",
    "Validation Count Pearson",
    "Validation Count MSE",
]


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--history", required=True,
        help="Training log TSV (<sample>.log) from the fit step.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("--sample", default="", help="Sample label for the title.")
    return parser


def main():
    args = build_parser().parse_args()

    from _style import apply_style
    apply_style()

    import math

    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.history, sep="\t")
    metrics = [m for m in METRICS if m in df.columns]

    # One subplot per metric (own y-scale); constrained_layout spaces the
    # suptitle so it does not collide with the axis titles.
    ncols = min(3, len(metrics))
    nrows = math.ceil(len(metrics) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows),
        squeeze=False)
    axes = axes.ravel()
    for ax, metric in zip(axes, metrics):
        ax.plot(df["Epoch"], df[metric], marker="o")
        ax.set_title(metric)
        ax.set_xlabel("Epoch")
    for ax in axes[len(metrics):]:
        ax.set_visible(False)

    if args.sample:
        fig.suptitle(args.sample)
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
