#!/usr/bin/env python
"""Plot the seqlet length distribution."""

import argparse


def build_parser():
    """Build the plot_seqlet_lengths.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_seqlet_lengths.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--seqlets", required=True,
        help="Seqlet BED from seqlets.py.")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("--sample", default="", help="Sample label for the title.")
    return parser


def main():
    """Plot one model's seqlet length distribution."""
    args = build_parser().parse_args()

    from _style import apply_style, despine, save_figure
    apply_style()

    import numpy
    import polars as pl
    import matplotlib.pyplot as plt

    bed = pl.read_csv(args.seqlets, separator="\t", has_header=False,
        new_columns=["chrom", "start", "end", "attribution", "p"])
    lengths = (bed["end"] - bed["start"]).to_numpy()

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(lengths, bins=numpy.arange(lengths.min(), lengths.max() + 2))
    ax.set_xlabel("Seqlet length (bp)")
    ax.set_ylabel("Seqlets")
    ax.set_title(f"{args.sample} (n={len(lengths):,} seqlets)".strip())
    despine(ax)

    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
