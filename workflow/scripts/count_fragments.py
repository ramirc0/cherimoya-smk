#!/usr/bin/env python
"""Count the fragments in a signal file (depth QC covariate).

Writes a single integer. For a fragments file, that is the number of non-comment
lines; for a BAM, the read-pair count (read1 records) under `--paired_end`, else
the read count. bigWig signals have no recoverable count and never reach here.
"""

import argparse


def build_parser():
    """Build the count_fragments.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the count_fragments.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--signal", required=True,
        help="Signal file: a BAM or a (bgzipped) fragments file.")
    parser.add_argument("-o", "--output", required=True,
        help="Destination text file (one integer).")
    parser.add_argument("--fragments", action="store_true", default=False,
        help="Signal is a fragments file (count lines), not a BAM.")
    parser.add_argument("--paired_end", action="store_true", default=False,
        help="BAM is paired-end: count read pairs (read1) as fragments.")
    return parser


def count(args):
    """Count the fragments of the signal, per its type.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed flags `signal`, `fragments` and `paired_end`.

    Returns
    -------
    int
        Non-comment lines of a fragments file. For a BAM, the read1 records
        under `paired_end`, else all reads.
    """
    if args.fragments:
        import gzip

        with gzip.open(args.signal, "rt") as fh:
            return sum(1 for line in fh if not line.startswith("#"))

    import pysam

    flags = ["-c", "-f", "64"] if args.paired_end else ["-c"]
    return int(pysam.view(*flags, args.signal))


def main():
    """Count the signal's fragments and write the count."""
    args = build_parser().parse_args()
    n = count(args)
    with open(args.output, "w") as fh:
        fh.write(f"{n}\n")
    print(f"{args.signal}: {n} fragments")


if __name__ == "__main__":
    main()
