#!/usr/bin/env python
"""Call seqlets from attributions and write them as genome coordinates (BED)."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l", "--loci", required=True,
        help="Peak file (narrowPeak/BED) that was attributed.")
    parser.add_argument("--ohe_filename", required=True,
        help="One-hot sequences from attribute.py (.npz).")
    parser.add_argument("--attr_filename", required=True,
        help="Hypothetical attributions from attribute.py (.npz).")
    parser.add_argument("--idx_filename", required=True,
        help="Mask of attributed loci from attribute.py (.npy).")
    parser.add_argument("--output_filename", default="seqlets.bed")
    parser.add_argument("--chroms", nargs="+", default=[
        "chr2", "chr4", "chr5", "chr7", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr21", "chr22", "chrX", "chrY", "chr8", "chr20"],
        help="Chromosomes attribute.py was run on.")
    parser.add_argument("--threshold", type=float, default=0.01)
    parser.add_argument("--min_seqlet_len", type=int, default=4)
    parser.add_argument("--max_seqlet_len", type=int, default=25)
    parser.add_argument("--additional_flanks", type=int, default=3)
    parser.add_argument("--n_bins", type=int, default=1000,
        help="Histogram bins for the attribution null distribution.")
    return parser


def main():
    args = build_parser().parse_args()

    import numpy
    import torch

    from tangermeme.io import _interleave_loci
    from tangermeme.seqlet import recursive_seqlets
    from tangermeme.utils import example_to_fasta_coords

    loci = _interleave_loci([args.loci], args.chroms)
    loci = loci.iloc[numpy.load(args.idx_filename)]

    X = torch.from_numpy(numpy.load(args.ohe_filename)["arr_0"])
    X_attr = torch.from_numpy(numpy.load(args.attr_filename)["arr_0"])
    X_attr = (X_attr * X).sum(dim=1)

    seqlets = recursive_seqlets(
        X_attr,
        threshold=args.threshold,
        min_seqlet_len=args.min_seqlet_len,
        max_seqlet_len=args.max_seqlet_len,
        additional_flanks=args.additional_flanks,
        n_bins=args.n_bins,
    ).sort_values("attribution", ascending=False)

    # An empty result has object-dtype columns that example_to_fasta_coords
    # cannot index loci with; the empty BED is still written.
    if len(seqlets) > 0:
        seqlets = example_to_fasta_coords(seqlets, loci, X.shape[-1])

    seqlets.to_csv(args.output_filename, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
