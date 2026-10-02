"""Shared matplotlib style for cherimoya-smk plots (pipeline + report).

Local copy of the matplotlib-style skill: the single source of plot style. Change
it here, never per script. `apply_style()` MUST run before `import matplotlib.pyplot`.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.style  # noqa: F401  (registers mpl.style)

FONT_STACK = [
    "Anthropic Sans Text",
    "Google Sans Flex",
    "Arimo",
    "Arial",
    "DejaVu Sans",
]

FORMATS = ("svg", "png")


def _register_fonts():
    """Register Anthropic Sans Text OTFs that matplotlib's cache does not see."""
    from matplotlib import font_manager

    for d in ("~/.local/share/fonts", "~/.fonts"):
        for otf in Path(d).expanduser().glob("**/AnthropicSans-Text-*.otf"):
            font_manager.fontManager.addfont(str(otf))


def apply_style():
    """Set the project plotting style (call before importing pyplot)."""
    mpl.use("Agg")
    mpl.style.use("default")
    _register_fonts()
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = FONT_STACK
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["figure.constrained_layout.use"] = True
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"


def despine(ax, categorical_x=False, categorical_y=False):
    """Offset left and bottom spines by 10 pt and trim them to the end ticks.

    Each continuous axis is fitted to its data without margins, then widened
    to the nearest ticks enclosing the data, so the trimmed spine never ends
    short of the data. When the locator has no enclosing tick (e.g. dates),
    the spine runs to the data edge without adding a tick. A categorical x or
    y axis (bar, heatmap) has no spine or tick marks; its labels carry the
    categories.

    Data at the limits sits on the axes edge, where matplotlib would clip half
    of each marker and line width. When both axes are fitted to the data, the
    plotted artists are unclipped so they draw whole into the spine offset.
    Explicit limits set by the caller keep clipping on.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to despine.
    categorical_x : bool, default False
        Treat the x axis as categorical.
    categorical_y : bool, default False
        Treat the y axis as categorical.
    """
    import seaborn as sns

    fitted = ax.get_autoscalex_on() and ax.get_autoscaley_on()
    ax.margins(0)
    ax.autoscale_view()
    sns.despine(ax=ax, bottom=categorical_x, left=categorical_y, offset=10)
    axes = []
    if not categorical_y:
        axes.append((ax.yaxis, ax.get_ylim, ax.set_ylim, "left"))
    if not categorical_x:
        axes.append((ax.xaxis, ax.get_xlim, ax.set_xlim, "bottom"))
    for axis, get_lim, set_lim, spine in axes:
        lo, hi = sorted(get_lim())
        # Calling the locator reads the view limits itself; date locators reject raw floats.
        ticks = axis.get_major_locator()()
        # lo and hi become the enclosing ticks, or stay at the data edges when
        # the locator has none (e.g. dates).
        lo = max((t for t in ticks if t <= lo), default=lo)
        hi = min((t for t in ticks if t >= hi), default=hi)
        axis.set_ticks([t for t in ticks if lo <= t <= hi])
        set_lim(lo, hi)
        ax.spines[spine].set_bounds(lo, hi)
    if fitted:
        for artist in [*ax.lines, *ax.collections, *ax.patches]:
            artist.set_clip_on(False)
    if categorical_x:
        ax.tick_params(axis="x", length=0)
    if categorical_y:
        ax.tick_params(axis="y", length=0)


def save_figure(fig, path, **kwargs):
    """Write `fig` to `path` as both SVG and PNG.

    `path` may end in `.svg`/`.png` (that extension is dropped) or carry none;
    the format extension is appended, so a dotted stem like `x.epochs` keeps its
    dots (`x.epochs.svg`, `x.epochs.png`).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    path : str or pathlib.Path
        Output path, with or without a `.svg` or `.png` extension.
    **kwargs
        Passed to `Figure.savefig`.

    Returns
    -------
    list of pathlib.Path
        The SVG and PNG paths written.
    """
    p = Path(path)
    base = p.with_suffix("") if p.suffix.lower() in (".svg", ".png") else p
    base.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in FORMATS:
        out = Path(f"{base}.{fmt}")
        fig.savefig(out, format=fmt, **kwargs)
        written.append(out)
    return written
