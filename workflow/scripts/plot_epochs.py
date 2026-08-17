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

    import pandas as pd
    import seaborn as sns

    df = pd.read_csv(args.history, sep="\t")
    metrics = [m for m in METRICS if m in df.columns]
    long = df.melt(
        id_vars="Epoch", value_vars=metrics, var_name="metric", value_name="value"
    )

    g = sns.relplot(
        data=long, x="Epoch", y="value", col="metric", col_wrap=3,
        kind="line", marker="o", facet_kws={"sharey": False},
    )
    g.set_titles("{col_name}")
    if args.sample:
        g.figure.suptitle(args.sample)
    g.figure.savefig(args.output)


if __name__ == "__main__":
    main()
