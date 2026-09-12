#!/usr/bin/env python
"""Collect per-model performance into one tidy table for run-level QC.

Joins every model's `<sample>.performance.tsv` with the sample sheet (genome,
n_fragments) and the peak count (lines in the sample's narrowPeak), keyed by
sample. `dataset` is the sample_id prefix before the first `__`. The `genome` and
`n_fragments` sheet columns are optional (NaN/None when absent).
"""

import argparse
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--performance", nargs="+", required=True,
        help="Per-model performance TSVs (<sample>.performance.tsv).")
    parser.add_argument("-s", "--samples", required=True,
        help="Sample sheet (sample_id, genome, [n_fragments], ...).")
    parser.add_argument("-r", "--results_dir", required=True,
        help="Run output dir (<results>/<run_id>) holding <sample>/ peak files.")
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


def main():
    args = build_parser().parse_args()

    import pandas as pd

    sheet = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample_id")
    has_genome = "genome" in sheet.columns
    has_nfrag = "n_fragments" in sheet.columns

    rows = []
    for path in args.performance:
        sample = Path(path).name[: -len(".performance.tsv")]
        in_sheet = sample in sheet.index
        metrics = pd.read_csv(path, sep="\t")
        # One row per signal group; these models are single-group, keep the mean.
        row = metrics.mean(numeric_only=True).to_dict()
        row["sample"] = sample
        row["dataset"] = sample.split("__")[0]
        row["genome"] = sheet.at[sample, "genome"] if has_genome and in_sheet else None
        nfrag = sheet.at[sample, "n_fragments"] if has_nfrag and in_sheet else None
        row["n_fragments"] = float(nfrag) if nfrag not in (None, "") else float("nan")
        row["n_peaks"] = _peak_count(args.results_dir, sample)
        rows.append(row)

    df = pd.DataFrame(rows)
    lead = ["sample", "dataset", "genome", "n_fragments", "n_peaks"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(df)} rows -> {args.output}")


if __name__ == "__main__":
    main()
