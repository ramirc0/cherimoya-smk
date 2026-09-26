#!/usr/bin/env python
"""Write a torchinfo layer summary and parameter count for a trained model."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m", "--model", required=True,
        help="Checkpoint written by fit (<name>.torch).")
    parser.add_argument("-o", "--output", required=True,
        help="Text file for the summary table.")
    parser.add_argument("--in_window", type=int, default=2114,
        help="Input sequence length of the dummy forward pass.")
    parser.add_argument("--depth", type=int, default=2,
        help="Module nesting depth shown in the table.")
    return parser


def main():
    args = build_parser().parse_args()

    import torch
    import torchinfo

    from cherimoya import Cherimoya

    # The checkpoint stores its constructor config, so this is the trained
    # architecture. Eager forward: torch.compile hides the hooks torchinfo uses.
    model = Cherimoya.load(args.model, device="cpu", compile=False)

    # torchinfo rejects None entries, so the control input is only passed
    # when the model has control tracks.
    inputs = [torch.zeros(1, 4, args.in_window)]
    if model.n_control_tracks:
        inputs.append(torch.zeros(1, model.n_control_tracks, args.in_window))

    stats = torchinfo.summary(model, input_data=inputs, depth=args.depth,
        verbose=0)
    with open(args.output, "w") as fh:
        fh.write(f"{stats}\n")


if __name__ == "__main__":
    main()
