#!/usr/bin/env python
"""Train a Cherimoya model on one sample's signal."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)

    # Data
    parser.add_argument("-s", "--sequences", required=True,
        help="Genome FASTA.")
    parser.add_argument("-l", "--loci", required=True,
        help="Peak file (narrowPeak/BED) of positive loci.")
    parser.add_argument("-neg", "--negatives", required=True,
        help="BED file of GC-matched negative loci.")
    parser.add_argument("-sig", "--signals", nargs="+", required=True,
        help="One or more signal bigWigs. A flat list of N files is N "
             "independent unstranded groups; with --stranded the files are "
             "one stranded (+, -) group.")
    parser.add_argument("-c", "--controls", nargs="+", default=None,
        help="Optional control bigWig(s), grouped like --signals.")
    parser.add_argument("--stranded", action="store_true", default=False,
        help="Treat the signal (and control) files as a single stranded "
             "(+, -) group instead of independent unstranded tracks.")
    parser.add_argument("-e", "--exclusion_lists", nargs="+", default=None,
        help="Optional BED files of regions to exclude.")
    parser.add_argument("-o", "--name", required=True,
        help="Output prefix; the best checkpoint is written to <name>.torch.")

    # Model
    parser.add_argument("--n_filters", type=int, default=128)
    parser.add_argument("--n_layers", type=int, default=9)
    parser.add_argument("--expansion", type=int, default=2)
    parser.add_argument("--residual_scale", type=float, default=0.15)

    # Windows / sampling
    parser.add_argument("--in_window", type=int, default=2114)
    parser.add_argument("--out_window", type=int, default=1000)
    parser.add_argument("--max_jitter", type=int, default=500)
    parser.add_argument("--negative_ratio", type=float, default=0.25)
    parser.add_argument("--reverse_complement", action="store_true", default=True)
    parser.add_argument("--no_reverse_complement", dest="reverse_complement",
        action="store_false")
    parser.add_argument("--summits", action="store_true", default=False)

    # Optimization
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--max_epochs", type=int, default=20)
    parser.add_argument("--min_total_steps", type=int, default=20000,
        help="Raise max_epochs until training takes this many steps; 0 "
             "disables the floor.")
    parser.add_argument("--loss_weights", nargs=2, type=float, default=None,
        metavar=("PROFILE", "COUNT"),
        help="Fixed loss weights replacing the learned Kendall weights.")
    parser.add_argument("--n_warmup_epochs", type=int, default=2)
    parser.add_argument("--early_stopping", type=int, default=None)
    parser.add_argument("--muon_lr", type=float, default=0.025)
    parser.add_argument("--muon_wd", type=float, default=0.03)
    parser.add_argument("--adam_lr", type=float, default=0.001)
    parser.add_argument("--adam_wd", type=float, default=0.0)
    parser.add_argument("--lw_lr", type=float, default=0.001)
    parser.add_argument("--lw_wd", type=float, default=0.0)
    parser.add_argument("--lw_momentum", type=float, default=0.9)

    # Runtime
    parser.add_argument("--num_workers", type=int, default=1)
    parser.add_argument("--dtype", default="float32")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--compile", action="store_true", default=True)
    parser.add_argument("--no_compile", dest="compile", action="store_false")
    parser.add_argument("--compile_mode", default="max-autotune")
    parser.add_argument("--random_state", type=int, default=0)
    parser.add_argument("--training_chroms", nargs="+", default=[
        "chr2", "chr4", "chr5", "chr7", "chr9", "chr10", "chr11", "chr12",
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr21",
        "chr22", "chrX", "chrY"])
    parser.add_argument("--validation_chroms", nargs="+",
        default=["chr8", "chr20"])
    parser.add_argument("-v", "--verbose", action="store_true", default=False)

    return parser


def main():
    args = build_parser().parse_args()

    import os

    os.environ["TORCH_CUDNN_V8_API_ENABLED"] = "1"

    import lightning
    import torch

    from cherimoya import Cherimoya
    from cherimoya.io import PeakGenerator, normalize_signal_groups
    from cherimoya.training import fit

    from tangermeme.io import _interleave_loci, extract_loci

    lightning.seed_everything(args.random_state, verbose=False)

    # --stranded wraps the flat file lists into one (+, -) group each; else
    # each file is its own unstranded group. The structured spec preserves
    # grouping for the model; the flat file list feeds extract_loci.
    signals = [args.signals] if args.stranded else args.signals
    controls = None if args.controls is None else (
        [args.controls] if args.stranded else args.controls)
    signal_files, signal_groups = normalize_signal_groups(signals)
    control_files, control_groups = normalize_signal_groups(controls)

    training_data = PeakGenerator(
        peaks=args.loci,
        negatives=args.negatives,
        sequences=args.sequences,
        signals=signals,
        controls=controls,
        chroms=args.training_chroms,
        in_window=args.in_window,
        out_window=args.out_window,
        max_jitter=args.max_jitter,
        negative_ratio=args.negative_ratio,
        reverse_complement=args.reverse_complement,
        summits=args.summits,
        exclusion_lists=args.exclusion_lists,
        random_state=args.random_state,
        verbose=args.verbose,
        signal_groups=signal_groups,
        control_groups=control_groups,
    ).dataset

    def _extract_valid(loci, summits):
        return extract_loci(
            sequences=args.sequences,
            signals=signal_files,
            in_signals=control_files,
            loci=loci,
            chroms=args.validation_chroms,
            in_window=args.in_window,
            out_window=args.out_window,
            max_jitter=0,
            summits=summits,
            exclusion_lists=args.exclusion_lists,
            ignore=list("QWERYUIOPSDFHJKLZXVBNM"),
            verbose=args.verbose,
        )

    # Peaks are extracted as in training; negatives have no summit column.
    valid_data = _extract_valid(args.loci, args.summits)

    # The validation negatives, labeled 0, feed the measures that separate
    # peaks from negatives. The other measures use the peaks alone.
    # `extract_loci` raises when no negative falls on the validation chroms.
    valid_labels = None
    if len(_interleave_loci(args.negatives, args.validation_chroms)) > 0:
        negative_data = _extract_valid(args.negatives, False)
        valid_labels = torch.cat([torch.ones(len(valid_data[0])),
            torch.zeros(len(negative_data[0]))])
        valid_data = [torch.cat(pair) for pair in zip(valid_data, negative_data)]

    if control_files is not None:
        valid_sequences, valid_signals, valid_controls = valid_data
        n_control_tracks = len(control_files)
    else:
        valid_sequences, valid_signals = valid_data
        valid_controls = None
        n_control_tracks = 0

    trimming = (args.in_window - args.out_window) // 2

    model = Cherimoya(
        n_filters=args.n_filters,
        n_layers=args.n_layers,
        signal_groups=signal_groups,
        n_control_tracks=n_control_tracks,
        expansion=args.expansion,
        residual_scale=args.residual_scale,
        trimming=trimming,
        name=args.name,
        verbose=args.verbose,
        compile=args.compile,
        compile_mode=args.compile_mode,
        random_state=args.random_state,
    )

    # Raise max_epochs to reach min_total_steps; the LR schedules stretch with it.
    # A trailing partial batch counts as a step, as upstream counts it.
    steps_per_epoch = -(-len(training_data) // args.batch_size)
    max_epochs = args.max_epochs
    if steps_per_epoch > 0 and steps_per_epoch * max_epochs < args.min_total_steps:
        max_epochs = -(-args.min_total_steps // steps_per_epoch)

    fit(
        model,
        training_data,
        valid_sequences,
        valid_signals,
        X_ctl_valid=valid_controls,
        labels_valid=valid_labels,
        max_epochs=max_epochs,
        early_stopping=args.early_stopping,
        dtype=args.dtype,
        accelerator={"cuda": "gpu"}.get(args.device, args.device),
        verbose=args.verbose,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        n_warmup_steps=steps_per_epoch * args.n_warmup_epochs,
        n_decay_steps=steps_per_epoch * max(1, max_epochs - args.n_warmup_epochs),
        muon_lr=args.muon_lr,
        muon_wd=args.muon_wd,
        adam_lr=args.adam_lr,
        adam_wd=args.adam_wd,
        lw_lr=args.lw_lr,
        lw_wd=args.lw_wd,
        lw_momentum=args.lw_momentum,
        loss_weights=args.loss_weights,
    )


if __name__ == "__main__":
    main()
