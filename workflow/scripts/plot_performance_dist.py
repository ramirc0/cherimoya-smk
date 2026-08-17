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

    from _style import apply_style
    apply_style()

    import pandas as pd
    import seaborn as sns

    # Each performance TSV holds one row per signal group; tag with the sample
    # (its parent dir name) so every model contributes to the distribution.
    frames = []
    for path in args.performance:
        d = pd.read_csv(path, sep="\t")
        d["sample"] = Path(path).parent.name
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)

    metrics = [c for c in df.columns if c != "sample"]
    long = df.melt(
        id_vars="sample", value_vars=metrics, var_name="metric", value_name="value"
    )

    g = sns.displot(
        data=long, x="value", col="metric", col_wrap=4, kind="hist", kde=True,
        facet_kws={"sharex": False, "sharey": False},
    )
    g.set_titles("{col_name}")
    g.figure.savefig(args.output)


if __name__ == "__main__":
    main()
