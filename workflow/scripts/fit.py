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
             "independent unstranded groups.")
    parser.add_argument("-c", "--controls", nargs="+", default=None,
        help="Optional control bigWig(s).")
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
    parser.add_argument("--n_warmup_epochs", type=int, default=2)
    parser.add_argument("--early_stopping", type=int, default=5)
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
    parser.add_argument("--random_state", type=int, default=None)
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

    import torch

    torch.backends.cudnn.benchmark = True
    torch.set_float32_matmul_precision("high")

    from torch.optim import Muon
    from torch.optim.lr_scheduler import (
        LinearLR,
        CosineAnnealingLR,
        ConstantLR,
        SequentialLR,
    )

    from cherimoya import Cherimoya
    from cherimoya.io import PeakGenerator, normalize_signal_groups

    from tangermeme.io import extract_loci

    # Flat file list feeds extract_loci; group sizes drive RC permutation.
    signal_files, signal_groups = normalize_signal_groups(args.signals)
    control_files, control_groups = normalize_signal_groups(args.controls)

    training_data = PeakGenerator(
        peaks=args.loci,
        negatives=args.negatives,
        sequences=args.sequences,
        signals=args.signals,
        controls=args.controls,
        chroms=args.training_chroms,
        in_window=args.in_window,
        out_window=args.out_window,
        max_jitter=args.max_jitter,
        negative_ratio=args.negative_ratio,
        reverse_complement=args.reverse_complement,
        summits=args.summits,
        exclusion_lists=args.exclusion_lists,
        random_state=args.random_state,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        verbose=args.verbose,
        signal_groups=signal_groups,
        control_groups=control_groups,
    )

    valid_data = extract_loci(
        sequences=args.sequences,
        signals=signal_files,
        in_signals=control_files,
        loci=args.loci,
        chroms=args.validation_chroms,
        in_window=args.in_window,
        out_window=args.out_window,
        max_jitter=0,
        exclusion_lists=args.exclusion_lists,
        ignore=list("QWERYUIOPSDFHJKLZXVBNM"),
        verbose=args.verbose,
    )

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
    ).to(args.device)

    num_warmup_iters = len(training_data) * args.n_warmup_epochs
    num_decay_iters = len(training_data) * max(1, args.max_epochs - args.n_warmup_epochs)

    # 2D projection weights -> Muon; lw0/lw1 -> SGD; everything else -> AdamW.
    muon_params, adam_params, lw_params = [], [], []
    for name, p in model.named_parameters():
        if name in ("lw0", "lw1"):
            lw_params.append(p)
        elif (
            p.ndim == 2
            and "weight" in name
            and name != "linear.weight"
            and "conv_weight" not in name
        ):
            muon_params.append(p)
        else:
            adam_params.append(p)

    muon_optimizer = Muon(muon_params, lr=args.muon_lr, weight_decay=args.muon_wd)
    muon_scheduler = SequentialLR(
        muon_optimizer,
        schedulers=[
            LinearLR(muon_optimizer, start_factor=0.01, total_iters=num_warmup_iters),
            CosineAnnealingLR(muon_optimizer, T_max=num_decay_iters, eta_min=1e-5),
        ],
        milestones=[num_warmup_iters],
    )

    adam_optimizer = torch.optim.AdamW(adam_params, lr=args.adam_lr, weight_decay=args.adam_wd)
    adam_scheduler = SequentialLR(
        adam_optimizer,
        schedulers=[
            LinearLR(adam_optimizer, start_factor=0.01, total_iters=num_warmup_iters),
            CosineAnnealingLR(adam_optimizer, T_max=num_decay_iters, eta_min=1e-5),
        ],
        milestones=[num_warmup_iters],
    )

    # lw weights hold a flat rate after warmup (no cosine decay).
    lw_optimizer = torch.optim.SGD(
        lw_params, lr=args.lw_lr, weight_decay=args.lw_wd, momentum=args.lw_momentum
    )
    lw_scheduler = SequentialLR(
        lw_optimizer,
        schedulers=[
            LinearLR(lw_optimizer, start_factor=0.01, total_iters=num_warmup_iters),
            ConstantLR(lw_optimizer, factor=1.0, total_iters=1),
        ],
        milestones=[num_warmup_iters],
    )

    model.fit(
        training_data,
        muon_optimizer,
        adam_optimizer,
        lw_optimizer,
        muon_scheduler,
        adam_scheduler,
        lw_scheduler,
        X_valid=valid_sequences,
        X_ctl_valid=valid_controls,
        y_valid=valid_signals,
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        early_stopping=args.early_stopping,
        dtype=args.dtype,
        device=args.device,
    )


if __name__ == "__main__":
    main()
