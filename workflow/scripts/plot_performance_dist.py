#!/usr/bin/env python
"""Plot the distribution of model performance, faceted across metrics."""

import argparse
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--performance", nargs="+", required=True,
        help="Per-model performance TSVs (<sample>.performance.tsv).")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    return parser


def main():
    args = build_parser().parse_args()

    from _style import apply_style, despine, save_figure
    apply_style()

    import math

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Each performance TSV holds one row per signal group; tag with the model id
    # (<sample>/fold_<k> from the path) so every model contributes to the distribution.
    frames = []
    for path in args.performance:
        d = pd.read_csv(path, sep="\t")
        d["model"] = f"{Path(path).parents[1].name}/{Path(path).parent.name}"
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    metrics = [c for c in df.columns if c != "model"]

    # One subplot per metric, each with its own x-scale -- the metrics span very
    # different ranges (e.g. profile_mnll ~130 vs spearman ~0.08), so a shared
    # axis crushes them. KDE only when a metric has more than one distinct value.
    ncols = min(4, len(metrics))
    nrows = math.ceil(len(metrics) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.8 * nrows),
        squeeze=False)
    axes = axes.ravel()
    for ax, metric in zip(axes, metrics):
        sns.histplot(data=df, x=metric, kde=df[metric].nunique() > 1, ax=ax)
        ax.set_title(metric)
        ax.set_xlabel("")
        despine(ax)
    for ax in axes[len(metrics):]:
        ax.set_visible(False)

    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
