#!/usr/bin/env python
"""Evaluate a trained Cherimoya model and write per-group performance metrics."""

import argparse

MEASURE_NAMES = [
    "profile_mnll",
    "profile_jsd",
    "profile_pearson",
    "profile_spearman",
    "count_pearson",
    "count_spearman",
    "count_mse",
]


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sequences", required=True, help="Genome FASTA.")
    parser.add_argument("-l", "--loci", required=True,
        help="Peak file (narrowPeak/BED) to evaluate on.")
    parser.add_argument("-sig", "--signals", nargs="+", required=True,
        help="Signal bigWig(s) matching the model's signal groups.")
    parser.add_argument("-c", "--controls", nargs="+", default=None,
        help="Optional control bigWig(s).")
    parser.add_argument("-e", "--exclusion_lists", nargs="+", default=None)
    parser.add_argument("-m", "--model", required=True,
        help="Trained model checkpoint (<name>.torch).")
    parser.add_argument("-o", "--performance_filename", required=True,
        help="Destination TSV for the performance table.")
    parser.add_argument("--counts_filename", default=None,
        help="Optional TSV of per-region observed vs predicted log-counts.")
    parser.add_argument("--chroms", nargs="+", default=["chr8", "chr20"])
    parser.add_argument("--batch_size", type=int, default=512)
    parser.add_argument("--in_window", type=int, default=2114)
    parser.add_argument("--out_window", type=int, default=1000)
    parser.add_argument("--reverse_complement_average", action="store_true",
        default=False)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dtype", default="float32")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    return parser


def main():
    args = build_parser().parse_args()

    import torch

    from tangermeme.io import extract_loci
    from tangermeme.predict import predict

    from cherimoya import Cherimoya, ControlWrapper
    from cherimoya.io import normalize_signal_groups
    from cherimoya.performance import calculate_performance_measures

    # Flatten structured signals/controls for extract_loci.
    signal_files, signal_groups = normalize_signal_groups(args.signals)
    control_files, _ = normalize_signal_groups(args.controls)

    model = Cherimoya.load(args.model, device=args.device)

    examples = extract_loci(
        sequences=args.sequences,
        signals=signal_files,
        in_signals=control_files,
        loci=args.loci,
        chroms=args.chroms,
        in_window=args.in_window,
        out_window=args.out_window,
        exclusion_lists=args.exclusion_lists,
        max_jitter=0,
        ignore=list("QWERYUIOPSDFHJKLZXVBNM"),
        verbose=args.verbose,
    )

    if control_files is None:
        X, y = examples
        X_ctl = None
        if model.n_control_tracks > 0:
            model = ControlWrapper(model)
    else:
        X, y, X_ctl = examples
        X_ctl = (X_ctl,)

    y_hat_logits, y_hat_logcounts = predict(
        model, X, args=X_ctl, batch_size=args.batch_size,
        device=args.device, dtype=args.dtype, verbose=args.verbose,
    )

    if args.reverse_complement_average:
        X_rc = torch.flip(X, dims=(-1, -2))
        X_ctl_rc = None if X_ctl is None else (torch.flip(X_ctl[0], dims=(-1, -2)),)

        y_hat_logits_rc, y_hat_logcounts_rc = predict(
            model, X_rc, args=X_ctl_rc, batch_size=args.batch_size,
            device=args.device, dtype=args.dtype, verbose=args.verbose,
        )

        y_hat_logits_rc = torch.flip(y_hat_logits_rc, dims=(-1, -2))
        y_hat_logits = (y_hat_logits + y_hat_logits_rc) / 2
        y_hat_logcounts = (y_hat_logcounts + y_hat_logcounts_rc) / 2

    # Prefer the checkpoint's own count-head layout.
    model_signal_groups = getattr(model, "signal_groups", None)
    if model_signal_groups is None:
        model_signal_groups = signal_groups

    measures = calculate_performance_measures(
        y_hat_logits, y, y_hat_logcounts, signal_groups=model_signal_groups
    )

    # One row per signal group: average profile metrics over the group's
    # channel slice; count metrics are already per-group.
    groups = model_signal_groups or [y_hat_logits.shape[1]]
    rows = []
    offset = 0
    for i, g in enumerate(groups):
        row = []
        for name in MEASURE_NAMES:
            value = measures[name]
            if name.startswith("profile_"):
                row.append(value[:, offset:offset + g].mean().item())
            else:
                row.append(value[i].item() if value.ndim >= 1 else value.item())
        rows.append(row)
        offset += g

    def _format_rows():
        yield "\t".join(MEASURE_NAMES)
        for row in rows:
            yield "\t".join(str(v) for v in row)

    with open(args.performance_filename, "w") as outfile:
        outfile.write("\n".join(_format_rows()))

    # Per-region observed vs predicted log-counts, the points behind count_pearson.
    if args.counts_filename:
        counts = y_hat_logcounts if y_hat_logcounts.ndim > 1 else y_hat_logcounts[:, None]
        lines = ["group\tobs_logcount\tpred_logcount"]
        offset = 0
        for i, g in enumerate(groups):
            obs = torch.log1p(y[:, offset:offset + g, :].sum(dim=(1, 2)))
            lines += [f"{i}\t{o}\t{p}"
                      for o, p in zip(obs.tolist(), counts[:, i].tolist())]
            offset += g
        with open(args.counts_filename, "w") as outfile:
            outfile.write("\n".join(lines))

    if args.verbose:
        for line in _format_rows():
            print(line)


if __name__ == "__main__":
    main()
