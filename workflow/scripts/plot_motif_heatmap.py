#!/usr/bin/env python
"""Heatmap of the top motifs across samples, as each sample's fraction of seqlets.

Motifs are the `--n_motifs` with the highest fraction in any sample, so a motif
strong in a single sample still shows. Rows and columns are ordered by
hierarchical clustering. Row labels are hidden past `--max_labels` samples.
"""

import argparse


def build_parser():
    """Build the plot_motif_heatmap.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_motif_heatmap.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--counts", nargs="+", required=True,
        help="Seqlets per motif (<results>/<sample>/fold_<k>/<sample>.motif_seqlet_count.tsv).")
    parser.add_argument("-o", "--output", required=True, help="Destination SVG.")
    parser.add_argument("-n", "--n_motifs", type=int, default=30)
    parser.add_argument("--max_labels", type=int, default=40,
        help="Most samples to label by name.")
    return parser


def _order(M):
    """Order the rows of `M` by average-linkage clustering.

    Parameters
    ----------
    M : numpy.ndarray
        Matrix whose rows are ordered.

    Returns
    -------
    numpy.ndarray
        Row indices in leaf order. The identity order below 3 rows.
    """
    import numpy
    from scipy.cluster.hierarchy import leaves_list, linkage

    if len(M) < 3:
        return numpy.arange(len(M))
    return leaves_list(linkage(M, method="average", optimal_ordering=True))


def main():
    """Plot the top motifs across models as a heatmap."""
    args = build_parser().parse_args()

    from _style import apply_style, despine, save_figure
    apply_style()

    from pathlib import Path

    import numpy
    import polars as pl
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm

    frames = []
    for path in args.counts:
        path = Path(path)
        frames.append(pl.read_csv(path, separator="\t").with_columns(
            pl.col("motifs").str.strip_chars(),
            (pl.col("count") / pl.col("count").sum()).alias("fraction"),
            sample=pl.lit(path.parents[1].name),
            fold=pl.lit(path.parent.name)))
    df = pl.concat(frames)
    if df.select("sample").n_unique() < len(args.counts):
        df = df.with_columns(sample=pl.format("{}/{}", "sample", "fold"))

    top = (df.group_by("motifs").agg(pl.col("fraction").max())
        .sort("fraction", descending=True).head(args.n_motifs)["motifs"])
    wide = (df.filter(pl.col("motifs").is_in(top.implode()))
        .pivot(on="motifs", index="sample", values="fraction")
        .fill_null(0))
    samples = wide["sample"].to_numpy()
    motifs = numpy.array(wide.columns[1:])
    M = wide.drop("sample").to_numpy()
    # sqrt keeps a dominant motif from flattening the distances between the rest.
    rows, cols = _order(numpy.sqrt(M)), _order(numpy.sqrt(M.T))
    M, samples, motifs = M[rows][:, cols], samples[rows], motifs[cols]

    n = len(samples)
    # Rows are 0.22 in tall until the heatmap reaches 9 in; the rest is labels.
    # The colorbar gets its own cell so its size does not follow the row count.
    height = min(0.22 * n, 9)
    fig, axes = plt.subplot_mosaic([[".", "cbar", "."], ["heat"] * 3],
        figsize=(2.5 + 0.2 * len(motifs), 3 + height), height_ratios=[0.12, height])
    ax = axes["heat"]
    mesh = ax.pcolormesh(M, cmap="viridis", norm=PowerNorm(0.5, vmin=0))
    ax.invert_yaxis()
    ax.set_xticks(numpy.arange(len(motifs)) + 0.5, motifs, rotation=90)
    if n <= args.max_labels:
        ax.set_yticks(numpy.arange(n) + 0.5, samples)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"Samples (n={n})")
    cbar = fig.colorbar(mesh, cax=axes["cbar"], orientation="horizontal")
    cbar.ax.set_title("Fraction of seqlets")
    # Evenly spaced under the sqrt norm; linear ticks crowd its top end.
    ticks = M.max() * numpy.linspace(0, 1, 5) ** 2
    cbar.set_ticks(ticks, labels=[f"{t:.2g}" for t in ticks])
    despine(ax, categorical_x=True, categorical_y=True)

    save_figure(fig, args.output)


if __name__ == "__main__":
    main()
