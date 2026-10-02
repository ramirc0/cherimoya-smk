#!/usr/bin/env python
"""Discover motifs with TF-MoDISco from one-hot sequences and hypothetical attributions."""

import argparse


def build_parser():
    """Build the modisco_motifs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the modisco_motifs.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True,
        help="One-hot sequences (.npz) from attribute.py.")
    parser.add_argument("-a", "--attributions", required=True,
        help="Hypothetical attributions (.npz) from attribute.py.")
    parser.add_argument("-o", "--output_filename", default="modisco_results.h5")
    parser.add_argument("-n", "--n_seqlets", type=int, default=100000,
        help="Maximum seqlets per metacluster.")
    parser.add_argument("-l", "--n_leiden", type=int, default=2,
        help="Leiden clusterings, each with its own seed.")
    parser.add_argument("-w", "--window", type=int, default=400,
        help="Central window of each locus used for discovery.")
    parser.add_argument("-z", "--size", type=int, default=20,
        help="Seqlet core size (sliding_window_size).")
    parser.add_argument("-t", "--trim_size", type=int, default=30,
        help="Pattern trim size (trim_to_window_size).")
    parser.add_argument("-f", "--seqlet_flank_size", type=int, default=5)
    parser.add_argument("-g", "--initial_flank_to_add", type=int, default=10)
    parser.add_argument("-j", "--final_flank_to_add", type=int, default=0)
    parser.add_argument("--target_seqlet_fdr", type=float, default=0.05,
        help="FDR for the seqlet attribution threshold.")
    parser.add_argument("--min_metacluster_size", type=int, default=100,
        help="Minimum seqlets for a metacluster to be clustered.")
    parser.add_argument("--n_leiden_iterations", type=int, default=-1,
        help="Leiden iterations per run; -1 runs until stable.")
    parser.add_argument("--final_min_cluster_size", type=int, default=20,
        help="Minimum seqlets for a final pattern.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    return parser


def main():
    """Run TF-MoDISco and write its results to h5."""
    args = build_parser().parse_args()

    import numpy
    import modiscolite

    sequences = numpy.load(args.sequences)["arr_0"]
    attributions = numpy.load(args.attributions)["arr_0"]
    if sequences.shape[2] < args.window:
        raise ValueError(f"window {args.window} exceeds sequence length {sequences.shape[2]}")

    start, end = modiscolite.util.calculate_window_offsets(sequences.shape[2] // 2, args.window)
    sequences = sequences[:, :, start:end].transpose(0, 2, 1).astype("float32")
    attributions = attributions[:, :, start:end].transpose(0, 2, 1).astype("float32")

    pos_patterns, neg_patterns = modiscolite.tfmodisco.TFMoDISco(
        hypothetical_contribs=attributions,
        one_hot=sequences,
        max_seqlets_per_metacluster=args.n_seqlets,
        sliding_window_size=args.size,
        flank_size=args.seqlet_flank_size,
        trim_to_window_size=args.trim_size,
        initial_flank_to_add=args.initial_flank_to_add,
        final_flank_to_add=args.final_flank_to_add,
        target_seqlet_fdr=args.target_seqlet_fdr,
        min_metacluster_size=args.min_metacluster_size,
        n_leiden_runs=args.n_leiden,
        n_leiden_iterations=args.n_leiden_iterations,
        final_min_cluster_size=args.final_min_cluster_size,
        verbose=args.verbose)

    modiscolite.io.save_hdf5(args.output_filename, pos_patterns, neg_patterns, args.window)


if __name__ == "__main__":
    main()
