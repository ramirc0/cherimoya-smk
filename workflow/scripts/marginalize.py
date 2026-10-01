#!/usr/bin/env python
"""Report the marginal effect of inserting each motif into background loci."""

import argparse
import os


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True, help="Genome FASTA.")
    parser.add_argument("-l", "--loci", required=True,
        help="Loci (narrowPeak/BED) to insert motifs into.")
    parser.add_argument("-m", "--model", required=True,
        help="Trained model checkpoint (<name>.torch).")
    parser.add_argument("-t", "--motifs", required=True, help="MEME motif file.")
    parser.add_argument("-e", "--exclusion_lists", nargs="+", default=None,
        help="Optional BED files of regions to exclude.")
    parser.add_argument("-o", "--output_filename", default="marginalize/",
        help="Report directory.")
    parser.add_argument("--chroms", nargs="+", default=[
        "chr2", "chr4", "chr5", "chr7", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr21", "chr22", "chrX", "chrY"])
    parser.add_argument("--n_loci", type=int, default=100,
        help="Loci used, in file order unless --shuffle.")
    parser.add_argument("--shuffle", action="store_true", default=False,
        help="Draw the loci at random from all loci on --chroms.")
    parser.add_argument("--random_state", type=int, default=0,
        help="Seed of --shuffle.")
    parser.add_argument("--attributions", action="store_true", default=False,
        help="Also report DeepLIFT/SHAP attributions before and after insertion.")
    parser.add_argument("--no_minimal", dest="minimal", action="store_false", default=True,
        help="Also plot the outputs before and after insertion, not only the difference.")
    parser.add_argument("--batch_size", type=int, default=512)
    parser.add_argument("--in_window", type=int, default=2114)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--no_compile", dest="compile", action="store_false", default=True)
    parser.add_argument("--compile_mode", default="max-autotune")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    return parser


def main():
    args = build_parser().parse_args()

    import numpy

    from bpnetlite.marginalize import marginalization_report
    from tangermeme.io import extract_loci

    from cherimoya import Cherimoya, ControlWrapper

    model = Cherimoya.load(args.model, device=args.device, compile=args.compile,
        compile_mode=args.compile_mode)
    if model.n_control_tracks > 0:
        model = ControlWrapper(model)

    # `extract_loci` keeps the first `n_loci` rows, so a shuffle must see them all.
    X = extract_loci(
        sequences=args.sequences,
        loci=args.loci,
        chroms=args.chroms,
        in_window=args.in_window,
        max_jitter=0,
        exclusion_lists=args.exclusion_lists,
        ignore=list("QWERYUIOPSDFHJKLZXVBNM"),
        n_loci=None if args.shuffle else args.n_loci,
        verbose=args.verbose,
    ).float()

    if args.shuffle:
        idxs = numpy.arange(X.shape[0])
        numpy.random.RandomState(args.random_state).shuffle(idxs)
        X = X[idxs]
    X = X[:args.n_loci]

    # The report joins `output_dir + name`, so it needs the trailing separator.
    marginalization_report(
        model,
        args.motifs,
        X,
        os.path.join(args.output_filename, ""),
        attributions=args.attributions,
        batch_size=args.batch_size,
        minimal=args.minimal,
        device=args.device,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
