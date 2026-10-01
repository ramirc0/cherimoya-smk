#!/usr/bin/env python
"""Compute hypothetical attributions of a trained Cherimoya model over peaks."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True, help="Genome FASTA.")
    parser.add_argument("-l", "--loci", required=True,
        help="Peak file (narrowPeak/BED) to attribute.")
    parser.add_argument("-m", "--model", required=True,
        help="Trained model checkpoint (<name>.torch).")
    parser.add_argument("-e", "--exclusion_lists", nargs="+", default=None,
        help="Optional BED files of regions to exclude.")
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

    import numpy

    from tangermeme.deep_lift_shap import deep_lift_shap
    from tangermeme.io import extract_loci
    from tangermeme.saturation_mutagenesis import saturation_mutagenesis

    from cherimoya import Cherimoya, ControlWrapper, LogCountWrapper, ProfileWrapper
    from cherimoya.deep_lift_shap import attribution_ops

    # DeepLIFT is never compiled: its backward hooks cause graph breaks.
    model = Cherimoya.load(args.model, device=args.device,
        compile=args.compile and args.algorithm == "saturation_mutagenesis",
        compile_mode=args.compile_mode)

    X, idxs = extract_loci(
        sequences=args.sequences,
        loci=args.loci,
        chroms=args.chroms,
        in_window=args.in_window,
        max_jitter=0,
        exclusion_lists=args.exclusion_lists,
        ignore=list("QWERYUIOPSDFHJKLZXVBNM"),
        return_mask=True,
        verbose=args.verbose,
    )

    # Drop loci containing N; `idxs` marks the loci kept.
    no_n = X.sum(dim=(1, 2)) == X.shape[-1]
    X = X[no_n]
    idxs[idxs.clone()] = no_n

    head = LogCountWrapper if args.output == "counts" else ProfileWrapper
    wrapper = head(ControlWrapper(model), group=args.group)

    if args.attr_window > X.shape[-1]:
        raise ValueError(f"attr_window ({args.attr_window}) is wider than "
            f"in_window ({X.shape[-1]})")
    start = X.shape[-1] // 2 - args.attr_window // 2
    end = start + args.attr_window

    if args.algorithm == "deep_lift_shap":
        X_attr = deep_lift_shap(
            wrapper,
            X,
            hypothetical=True,
            n_shuffles=args.n_shuffles,
            batch_size=args.batch_size,
            warning_threshold=args.warning_threshold,
            additional_nonlinear_ops=attribution_ops(),
            dtype=args.dtype,
            device=args.device,
            random_state=args.random_state,
            verbose=args.verbose,
        )[:, :, start:end].float()
    else:
        X_attr = saturation_mutagenesis(
            wrapper,
            X,
            dtype=args.dtype,
            device=args.device,
            batch_size=args.batch_size,
            verbose=args.verbose,
            hypothetical=True,
            start=start,
            end=end,
        ).float()

    numpy.savez_compressed(args.ohe_filename, X[:, :, start:end])
    numpy.savez_compressed(args.attr_filename, X_attr)
    numpy.save(args.idx_filename, idxs)


if __name__ == "__main__":
    main()
