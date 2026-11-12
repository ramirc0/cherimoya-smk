"""Observed vs predicted log-count scatter for one model (count-head QC).

Every point is a test peak: x = observed log-counts, y = the model's predicted
log-counts. Points are colored by density as the lab's R heatscatters are, and
the title has the Pearson r that count_pearson summarises, so a weak count head
is visible at a glance.
"""

import argparse


def build_parser():
    """Build the plot_count_scatter.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_count_scatter.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--counts",
        required=True,
        help="Per-region counts TSV (group, obs_logcount, pred_logcount).",
    )
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("--sample", default="", help="Sample name for the title.")
    return parser


def main():
    """Plot one model's observed vs predicted log-counts."""
    args = build_parser().parse_args()

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import poikilos as pk
    from heatscatter import heatscatter

    pk.use("zed")
    zed = pk.PALETTES["zed"].colors

    df = pd.read_csv(args.counts, sep="\t")
    x = df["obs_logcount"].to_numpy(float)
    y = df["pred_logcount"].to_numpy(float)
    r = np.corrcoef(x, y)[0, 1]

    fig, ax = plt.subplots(figsize=(4, 4))
    heatscatter(ax, x, y, size=2)

    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax.plot([lo, hi], [lo, hi], color=zed["black"], label="y = x")

    ax.set_xlabel("observed log-counts")
    ax.set_ylabel("predicted log-counts")
    # Wrap the (long) sample name so the title fits the figure width.
    import textwrap

    header = "\n".join(textwrap.wrap(args.sample, 30)) + "\n" if args.sample else ""
    ax.set_title(f"{header}count Pearson r = {r:.3f}  (n={len(x)})")
    ax.legend(loc="upper left")
    pk.despine(ax, offset=7.5)
    pk.save_figure(fig, args.output)


if __name__ == "__main__":
    main()
