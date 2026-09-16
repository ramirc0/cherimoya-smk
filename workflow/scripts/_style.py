"""Shared matplotlib style for cherimoya-smk plots."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.style  # noqa: F401  (registers mpl.style)


def apply_style():
    """Set the project plotting style (call before importing pyplot)."""
    mpl.use("Agg")
    mpl.style.use("default")
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["savefig.format"] = "svg"
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["figure.constrained_layout.use"] = True


def save_figure(fig, output):
    """Write `output` (SVG) plus a 300-DPI PNG sibling (same stem, .png)."""
    fig.savefig(output)
    fig.savefig(Path(output).with_suffix(".png"), dpi=300)
