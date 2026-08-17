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
    return parser


def main():
    args = build_parser().parse_args()

    from tangermeme.match import extract_matching_loci

    matched_loci = extract_matching_loci(
        loci=args.peaks,
        fasta=args.fasta,
        gc_bin_width=args.bin_width,
        max_n_perc=args.max_n_perc,
        bigwig=args.bigwig,
        signal_beta=args.beta,
        in_window=args.in_window,
        out_window=args.out_window,
        chroms=None,
        verbose=args.verbose,
        n_jobs=1,
    )

    matched_loci.to_csv(args.output, header=False, sep="\t", index=False)


if __name__ == "__main__":
    main()
