#!/usr/bin/env python
"""Attribute a trained Cherimoya model over peaks with `cherimoya attribute`."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True, help="Genome FASTA.")
    parser.add_argument("-l", "--loci", required=True,
        help="Peak file (narrowPeak/BED) to attribute.")
    parser.add_argument("-m", "--model", required=True,
        help="Trained model checkpoint (<name>.torch).")
    parser.add_argument("--ohe_filename", default="attributions.ohe.npz")
    parser.add_argument("--attr_filename", default="attributions.attr.npz")
    parser.add_argument("--idx_filename", default="attributions.idx.npy")
    parser.add_argument("--chroms", nargs="+", default=[
        "chr2", "chr4", "chr5", "chr7", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr21", "chr22", "chrX", "chrY", "chr8", "chr20"])
    parser.add_argument("--algorithm", default="deep_lift_shap",
        choices=["deep_lift_shap", "saturation_mutagenesis"])
    parser.add_argument("--output", default="counts", choices=["counts", "profile"])
    parser.add_argument("--group", type=int, default=0,
        help="Signal group of the output head to attribute.")
    parser.add_argument("--attr_window", type=int, default=400,
        help="Width of the saved central slice.")
    parser.add_argument("--n_shuffles", type=int, default=20)
    parser.add_argument("--warning_threshold", type=float, default=1e-3)
    parser.add_argument("--random_state", type=int, default=0)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--in_window", type=int, default=2114)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dtype", default="float32")
    parser.add_argument("--compile", action="store_true", default=False,
        help="Compile the model (saturation_mutagenesis only).")
    parser.add_argument("--compile_mode", default="max-autotune")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    return parser


def main():
    args = build_parser().parse_args()

    from cherimoya_cli.commands import attribute

    attribute.run(argparse.Namespace(parameters=vars(args)))


if __name__ == "__main__":
    main()
