#!/usr/bin/env python
"""Identify GC-matched negative regions for a set of peaks."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--peaks", required=True,
        help="A peak file (narrowPeak/BED) to draw GC-matched negatives against.")
    parser.add_argument("-f", "--fasta", required=True,
        help="A FASTA file of the genome the peaks are defined on.")
    parser.add_argument("-b", "--bigwig", default=None,
        help="Optional bigWig used to also match on signal.")
    parser.add_argument("-o", "--output", required=True,
        help="Destination BED file for the matched negative loci.")
    parser.add_argument("-l", "--bin_width", type=float, default=0.02,
        help="GC bin width for matching. Default 0.02.")
    parser.add_argument("-n", "--max_n_perc", type=float, default=0.1,
        help="Maximum fraction of Ns allowed in a locus. Default 0.1.")
    parser.add_argument("-a", "--beta", type=float, default=0.5,
        help="Signal-matching strength when a bigWig is given. Default 0.5.")
    parser.add_argument("-w", "--in_window", type=int, default=2114,
        help="Input window size. Default 2114.")
    parser.add_argument("-x", "--out_window", type=int, default=1000,
        help="Output window size. Default 1000.")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--on_missing_contig", choices=("filter", "fail"),
        default="filter",
        help="What to do with peaks on contigs absent from the FASTA. "
             "filter: drop and warn (default); fail: error out.")
    return parser


def resolve_present_chroms(peaks_path, fasta_path, policy):
    """Return the sorted peak contigs that are present in the FASTA.

    Peaks on contigs missing from the FASTA would crash extract_matching_loci
    (a pyfaidx KeyError when it builds chrom sizes from the peak contigs). This
    filters them out under ``policy="filter"`` (warning to stderr) or errors
    under ``policy="fail"``. It always errors if no peaks remain, since that
    signals the wrong reference genome rather than a stray contig.
    """
    import sys

    import pandas as pd
    import pyfaidx

    peak_chroms = pd.read_csv(peaks_path, sep="\t", usecols=[0], header=None,
        names=["chrom"], dtype=str)["chrom"]
    fasta_contigs = set(map(str, pyfaidx.Fasta(fasta_path).keys()))

    present = sorted(set(peak_chroms) & fasta_contigs)
    missing = sorted(set(peak_chroms) - fasta_contigs)

    if missing:
        n_dropped = int(peak_chroms.isin(missing).sum())
        detail = (f"{n_dropped} peak(s) on {len(missing)} contig(s) absent from the "
                  f"FASTA: {', '.join(missing)}")
        if policy == "fail":
            raise SystemExit(
                f"[negatives] ERROR: {detail}. Fix the peaks/reference or set "
                "peaks.on_missing_contig=filter."
            )
        print(f"[negatives] WARNING: dropping {detail}.", file=sys.stderr)

    if not present:
        raise SystemExit(
            "[negatives] ERROR: no peaks remain on FASTA contigs; wrong reference genome?"
        )

    return present


def main():
    args = build_parser().parse_args()

    from tangermeme.match import extract_matching_loci

    present = resolve_present_chroms(args.peaks, args.fasta, args.on_missing_contig)

    matched_loci = extract_matching_loci(
        loci=args.peaks,
        fasta=args.fasta,
        gc_bin_width=args.bin_width,
        max_n_perc=args.max_n_perc,
        bigwig=args.bigwig,
        signal_beta=args.beta,
        in_window=args.in_window,
        out_window=args.out_window,
        chroms=present,
        verbose=args.verbose,
        n_jobs=1,
    )

    matched_loci.to_csv(args.output, header=False, sep="\t", index=False)


if __name__ == "__main__":
    main()
