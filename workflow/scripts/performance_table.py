"""Summarize model performance per metric as an HTML table and a PNG of it."""

import argparse

from performance_metrics import METRICS

STATS = ["mean", "std", "min", "25%", "50%", "75%", "max"]


def build_parser():
    """Build the performance_table.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the performance_table.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--performance",
        nargs="+",
        required=True,
        help="Per-model performance TSVs (<sample>.performance.tsv).",
    )
    parser.add_argument("-o", "--output", required=True, help="Destination HTML.")
    parser.add_argument(
        "--png", required=True, help="Destination PNG, rendered in headless Firefox."
    )
    return parser


def main():
    """Write one row per metric with its spread across models."""
    args = build_parser().parse_args()

    import poikilos as pk
    import polars as pl
    from great_tables import GT

    # Match the zed figures' fonts and sizes.
    zed = pk.rc_params("zed")

    # evaluate.py writes "nan" when a metric is undefined. Read as text, since a
    # column of only "nan" would infer as str, then count NaN as missing.
    performance = (
        pl.concat(
            [
                pl.read_csv(path, separator="\t", infer_schema=False)
                for path in args.performance
            ]
        )
        .cast(pl.Float64)
        .with_columns(pl.all().fill_nan(None))
    )
    summary = (
        performance.describe()
        .transpose(include_header=True, header_name="metric", column_names="statistic")
        .sort(
            pl.col("metric").replace_strict(
                {m: i for i, m in enumerate(METRICS)}, default=len(METRICS)
            )
        )
        .with_columns(
            pl.col("metric")
            .replace_strict({m: g for m, (g, _, _) in METRICS.items()}, default="Other")
            .alias("group"),
            pl.col("metric")
            .replace_strict({m: n for m, (_, n, _) in METRICS.items()}, default=None)
            .fill_null(pl.col("metric"))
            .alias("name"),
        )
    )

    table = (
        GT(summary, rowname_col="name", groupname_col="group")
        .tab_header(
            title="Model performance",
            subtitle=f"{len(args.performance):,} models, on each fold's test chromosomes",
        )
        .cols_hide("metric")
        .cols_move_to_start(["count", "null_count"])
        .tab_spanner("Across models", STATS)
        .cols_label(
            count="Models",
            null_count="NaN",
            mean="Mean",
            std="SD",
            min="Min",
            **{"25%": "Q1", "50%": "Median", "75%": "Q3"},
            max="Max",
        )
        .fmt_integer(["count", "null_count"])
        .fmt_number(STATS, decimals=3)
        .fmt_number(STATS, rows=pl.col("metric") == "profile_mnll", decimals=1)
        .sub_missing(missing_text="")
        .tab_source_note(
            "Lower is better for NLL, Jensen-Shannon and MSE. Higher is better "
            "for the rest. NaN counts models where the metric is undefined."
        )
        .opt_table_font(font=[*zed["font.sans-serif"], "sans-serif"])
        .tab_options(
            table_font_size=f"{zed['font.size']}pt",
            heading_title_font_size=f"{zed['figure.titlesize']}pt",
            heading_subtitle_font_size=f"{zed['font.size']}pt",
        )
    )
    table.write_raw_html(args.output, make_page=True)
    table.save(args.png, web_driver="firefox", scale=2)


if __name__ == "__main__":
    main()
