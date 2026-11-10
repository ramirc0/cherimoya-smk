"""Plot the distribution of model performance, one PDF page per metric family."""

import argparse
from pathlib import Path

from performance_metrics import METRICS

# Page slug -> (title, rows of metrics). Counts come first, peaks above all loci
# so each column compares the same statistic.
PAGES = {
    "counts": (
        "Counts",
        [
            ["count_pearson", "count_spearman", "count_mse"],
            ["all_count_pearson", "all_count_spearman", "all_count_mse"],
        ],
    ),
    "profile": (
        "Profile shape",
        [["profile_mnll", "profile_jsd", "profile_pearson", "profile_spearman"]],
    ),
    "peak_vs_negative": ("Peak vs negative", [["auroc", "auprc"]]),
}


def build_parser():
    """Build the plot_performance_dist.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_performance_dist.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--performance",
        nargs="+",
        required=True,
        help="Per-model performance TSVs (<sample>.performance.tsv).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Destination PDF. Each page is also saved as <stem>.<page>.{svg,png}.",
    )
    return parser


def main():
    """Plot the distribution of model performance, one page per metric family."""
    args = build_parser().parse_args()

    import matplotlib.pyplot as plt
    import pandas as pd
    import poikilos as pk
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import to_rgb
    from matplotlib.ticker import MaxNLocator
    from matplotlib.transforms import offset_copy

    pk.use("zed")
    zed = pk.PALETTES["zed"].colors
    center_styles = {
        "median": {"color": zed["red"]},
        "mean": {"color": zed["blue"], "linestyle": "--"},
    }
    # Value pills sit side by side under each panel title, so each names its line
    # and no legend is needed.
    pill_sides = {"median": (-6, "right"), "mean": (6, "left")}

    df = pd.concat(
        [pd.read_csv(path, sep="\t") for path in args.performance], ignore_index=True
    )
    stem = Path(args.output).with_suffix("")

    # A peak with no reads on a strand makes its JSD NaN, and the per-model mean
    # with it. Runs from before 2026-10 lack the all-loci and peak-vs-negative
    # columns.
    notes = {}
    for metric in METRICS:
        if metric not in df:
            notes[metric] = "not in this run's evaluate output"
        elif df[metric].isna().all():
            notes[metric] = "NaN in every model"

    with PdfPages(args.output) as pdf:
        for slug, (title, rows) in PAGES.items():
            ncols = max(len(row) for row in rows)
            # An empty panel would take its column's shared tick labels with it.
            shared = len(rows) > 1 and not any(m in notes for row in rows for m in row)
            fig, axes = plt.subplots(
                len(rows),
                ncols,
                figsize=(3.2 * ncols, 2.8 * len(rows)),
                sharex="col" if shared else False,
                squeeze=False,
            )
            fig.suptitle(title)
            for r, row in enumerate(rows):
                for ax, metric in zip(axes[r], row):
                    group, name, unit = METRICS[metric]
                    ax.set_title(name, pad=24)
                    if metric in notes:
                        ax.text(
                            0.5,
                            0.5,
                            notes[metric],
                            transform=ax.transAxes,
                            ha="center",
                            va="center",
                            color=zed["k50"],
                        )
                        ax.set_axis_off()
                        continue
                    sns.histplot(
                        data=df,
                        x=metric,
                        kde=df[metric].nunique() > 1,
                        color=zed["k20"],
                        linewidth=0,
                        ax=ax,
                    )
                    # histplot draws the KDE in the bar color and ignores line_kws color.
                    for line in ax.lines:
                        line.set_color(zed["black"])
                    for stat, style in center_styles.items():
                        value = getattr(df[metric], stat)()
                        ax.axvline(value, **style)
                        dx, ha = pill_sides[stat]
                        color = style["color"]
                        ax.text(
                            0.5,
                            1.0,
                            f"{stat} {value:.3g}",
                            transform=offset_copy(
                                ax.transAxes, fig=fig, x=dx, y=6, units="points"
                            ),
                            ha=ha,
                            va="bottom",
                            fontsize="small",
                            bbox={
                                "boxstyle": "round,pad=0.3",
                                "facecolor": [
                                    1 - 0.18 * (1 - c) for c in to_rgb(color)
                                ],
                                "edgecolor": color,
                                "linewidth": 0.75,
                            },
                        )
                    ax.xaxis.set_major_locator(MaxNLocator(5, steps=[1, 2, 5, 10]))
                    ax.yaxis.set_major_locator(
                        MaxNLocator(integer=True, steps=[1, 2, 5, 10])
                    )
                    ax.set_xlabel(unit if not shared or r == len(rows) - 1 else "")
                    ax.set_ylabel("")
                    pk.despine(ax, offset=7.5)
                # Several rows on a page differ by group, so the row label names it.
                label = f"{group}\nmodels" if len(rows) > 1 else "models"
                axes[r][0].set_ylabel(label)
                for ax in axes[r][len(row) :]:
                    ax.set_visible(False)
            pk.save_figure(fig, f"{stem}.{slug}")
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
