#!/usr/bin/env python
"""Plot the mean attribution by position across all attributed peaks."""

import argparse


def build_parser():
    """Build the plot_attribution_profile.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_attribution_profile.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ohe_filename", required=True,
        help="One-hot sequences from attribute.py (.npz).")
    parser.add_argument("--attr_filename", required=True,
        help="Hypothetical attributions from attribute.py (.npz).")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("--sample", default="", help="Sample label for the title.")
    return parser


def main():
    """Plot one model's mean attribution by position."""
    args = build_parser().parse_args()

    from _style import apply_style, despine, save_figure
    apply_style()

    import numpy
    import matplotlib.pyplot as plt

    X = numpy.load(args.ohe_filename)["arr_0"]
    X_attr = numpy.load(args.attr_filename)["arr_0"]
    # Attribution of the observed base, as seqlets.py scores it.
    attr = (X_attr * X).sum(axis=1)
    pos = numpy.arange(attr.shape[1]) - attr.shape[1] // 2

    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.plot(pos, numpy.abs(attr).mean(axis=0), lw=1, label="Mean absolute")
    ax.plot(pos, attr.mean(axis=0), lw=1, label="Mean")
    ax.axhline(0, color="0.6", lw=0.5, zorder=0)
    ax.set_xlabel("Position from peak center (bp)")
    ax.set_ylabel("Attribution")
    ax.set_title(f"{args.sample} (n={attr.shape[0]:,} peaks)".strip())
    ax.legend(frameon=False)
    despine(ax)

    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
