#!/usr/bin/env python
"""Collect per-model performance into one tidy table for run-level QC.

One row per (sample, CV fold): joins each model's `<sample>.performance.tsv`
(fold taken from its `fold_<k>/` parent dir) with the sample sheet's `genome` and
the requested QC covariates. `dataset` is the sample_id prefix before the first
`__`. Covariate columns are emitted only when named in `--covariates`: `n_peaks`
(lines in the narrowPeak) and `n_fragments` (from the count_fragments output);
each is NaN when its source file is absent (e.g. a bigWig signal has no depth).
"""

import argparse
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--performance", nargs="+", required=True,
        help="Per-model performance TSVs (<sample>.performance.tsv).")
    parser.add_argument("-s", "--samples", required=True,
        help="Sample sheet (sample_id, genome, ...).")
    parser.add_argument("-r", "--results_dir", required=True,
        help="Run output dir (<results>/<run_id>) holding per-sample QC files.")
    parser.add_argument("--covariates", nargs="*", default=[],
        choices=["n_peaks", "n_fragments"],
        help="QC covariate columns to emit.")
    parser.add_argument("-o", "--output", required=True, help="Destination TSV.")
    return parser


def _peak_count(results_dir, sample):
    """Lines in the sample's peak file (macs3 or provided), else NaN."""
    for name in (f"{sample}_peaks.narrowPeak", f"{sample}.peaks.narrowPeak"):
        path = Path(results_dir) / sample / name
        if path.exists():
            with open(path) as fh:
                return sum(1 for _ in fh)
    return float("nan")


def _fragment_count(results_dir, sample):
    """Fragment count from count_fragments (<sample>.n_fragments.txt), else NaN."""
    path = Path(results_dir) / sample / f"{sample}.n_fragments.txt"
    if path.exists():
        return int(path.read_text().strip())
    return float("nan")


def main():
    args = build_parser().parse_args()

    import pandas as pd

    sheet = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample_id")
    has_genome = "genome" in sheet.columns

    rows = []
    for path in args.performance:
        sample = Path(path).name[: -len(".performance.tsv")]
        # Fold lives in the parent dir (<sample>/fold_<k>/), not the filename.
        row_fold = Path(path).parent.name.removeprefix("fold_")
        in_sheet = sample in sheet.index
        metrics = pd.read_csv(path, sep="\t")
        # One row per signal group; these models are single-group, keep the mean.
        row = metrics.mean(numeric_only=True).to_dict()
        row["sample"] = sample
        row["fold"] = row_fold
        row["dataset"] = sample.split("__")[0]
        row["genome"] = sheet.at[sample, "genome"] if has_genome and in_sheet else None
        if "n_peaks" in args.covariates:
            row["n_peaks"] = _peak_count(args.results_dir, sample)
        if "n_fragments" in args.covariates:
            row["n_fragments"] = _fragment_count(args.results_dir, sample)
        rows.append(row)

    df = pd.DataFrame(rows)
    lead = ["sample", "fold", "dataset", "genome", *args.covariates]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(df)} rows -> {args.output}")


if __name__ == "__main__":
    main()
