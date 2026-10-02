#!/usr/bin/env python
"""Annotate seqlets with their nearest motifs by TomTom and count seqlets per motif."""

import argparse


def build_parser():
    """Build the annotate.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the annotate.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True, help="Genome FASTA.")
    parser.add_argument("-b", "--seqlet_filename", required=True,
        help="Seqlet BED from seqlets.py.")
    parser.add_argument("-t", "--motifs", required=True, help="Motif database (MEME).")
    parser.add_argument("--output_filename", default="seqlets_annotated.bed")
    parser.add_argument("--count_filename", default="motif_seqlet_count.tsv",
        help="Seqlets per best-matching motif.")
    parser.add_argument("--n_score_bins", type=int, default=100)
    parser.add_argument("--n_median_bins", type=int, default=1000)
    parser.add_argument("--n_target_bins", type=int, default=100,
        help="0 disables approximate target hashing (exact p-values).")
    parser.add_argument("--n_cache", type=int, default=250)
    parser.add_argument("--n_nearest", type=int, default=1,
        help="Motifs reported per seqlet, as (motif, -log p) column pairs.")
    parser.add_argument("--no_reverse_complement", dest="reverse_complement",
        action="store_false", default=True)
    parser.add_argument("--n_jobs", type=int, default=-1)
    return parser


def main():
    """Annotate the seqlets and count them per motif."""
    args = build_parser().parse_args()

    import math

    import pandas
    import pyfaidx

    from memelite.io import read_meme
    from memelite.tomtom import tomtom
    from memelite.utils import one_hot_encode

    df = pandas.read_csv(args.seqlet_filename, sep="\t", usecols=(0, 1, 2),
        header=None, names=["chrom", "start", "end"])
    fa = pyfaidx.Fasta(args.sequences)
    targets = read_meme(args.motifs)
    names = list(targets)

    seqs = [one_hot_encode(fa[c][s:e].seq.upper())
            for c, s, e in zip(df["chrom"], df["start"], df["end"])]
    p, _, _, _, _, idxs = tomtom(seqs, list(targets.values()),
        n_nearest=args.n_nearest, n_score_bins=args.n_score_bins,
        n_median_bins=args.n_median_bins,
        n_target_bins=args.n_target_bins or None, n_cache=args.n_cache,
        reverse_complement=args.reverse_complement, n_jobs=args.n_jobs)

    # Fixed-width columns, as `ttl -f -b` prints them.
    motif = "\t{:" + str(max(map(len, names))) + "}\t{:6.6}"
    fmt = ("{:" + str(max(map(len, fa.keys()))) + "}"
        + "\t{:" + str(len(str(int(df["start"].max())))) + "}"
        + "\t{:" + str(len(str(int(df["end"].max())))) + "}"
        + motif * args.n_nearest + "\n")

    best = []
    with open(args.output_filename, "w") as out:
        for i, row in enumerate(zip(df["chrom"], df["start"], df["end"])):
            hits = []
            for k in range(args.n_nearest):
                hits += [names[int(idxs[i, k])],
                         -math.log(p[i, k]) if p[i, k] > 0 else float("inf")]
            line = fmt.format(*row, *hits)
            best.append(line.split("\t")[3])
            out.write(line)

    counts = pandas.DataFrame({"motifs": best}).value_counts()
    counts.to_csv(args.count_filename, sep="\t")


if __name__ == "__main__":
    main()
