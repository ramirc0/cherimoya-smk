#!/usr/bin/env python
"""Plot the motifs that best match the most seqlets."""

import argparse


def build_parser():
    """Build the plot_motif_counts.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_motif_counts.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--counts", required=True,
        help="Seqlets per motif from annotate.py (motif_seqlet_count.tsv).")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-n", "--n_motifs", type=int, default=20,
        help="Top motifs to show.")
    parser.add_argument("--sample", default="", help="Sample label for the title.")
    return parser


def main():
    """Plot one model's top motifs by seqlet count."""
    args = build_parser().parse_args()

    from _style import apply_style, despine, save_figure
    apply_style()

    import polars as pl
    import matplotlib.pyplot as plt

    counts = pl.read_csv(args.counts, separator="\t")
    total = counts["count"].sum()
    top = (counts.with_columns(pl.col("motifs").str.strip_chars())
        .sort("count", descending=True).head(args.n_motifs).reverse())

    fig, ax = plt.subplots(figsize=(4.5, 0.6 + 0.18 * len(top)))
    ax.barh(top["motifs"], top["count"] / total)
    ax.set_xlabel("Fraction of seqlets")
    ax.set_title(f"{args.sample} (n={total:,} seqlets)".strip())
    despine(ax, categorical_y=True)

    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
