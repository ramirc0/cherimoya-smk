"""Density-colored scatters ("heatscatters") as the Weng lab draws them in R.

A port of the lab's ggplot2 recipe (Pratt et al. 2024, Sci Adv): each point
takes the ``MASS::kde2d(x, y, n = 256)`` density at the grid node at or below
it, and points are drawn in input order on a reversed Spectral ramp.
``heatscatter-vignette`` checks the port against R to 4e-7 of the peak density.
"""

import numpy as np


def bandwidth_nrd(values):
    """Return R's ``MASS::bandwidth.nrd``, the normal-reference bandwidth.

    Parameters
    ----------
    values : numpy.ndarray
        One coordinate of the points.

    Returns
    -------
    float
        ``4 * 1.06 * min(sd, IQR / 1.34) * n ** (-1/5)``.
    """
    q25, q75 = np.quantile(values, [0.25, 0.75])
    spread = min(np.std(values, ddof=1), (q75 - q25) / 1.34)
    return 4 * 1.06 * spread * len(values) ** (-1 / 5)


def kde2d(x, y, n=256, chunk=200_000):
    """Return ``MASS::kde2d(x, y, n = n)`` with its default bandwidths.

    The Gaussian kernel matrices are built in chunks of points, so millions of
    points fit in memory.

    Parameters
    ----------
    x : numpy.ndarray
        Point x coordinates.
    y : numpy.ndarray
        Point y coordinates.
    n : int
        Grid points per axis. The grid spans the range of each coordinate.
    chunk : int
        Points per kernel block.

    Returns
    -------
    tuple of numpy.ndarray
        Grid ``gx`` (n,), grid ``gy`` (n,) and density ``z`` (n, n), indexed
        ``z[ix, iy]`` as in R.
    """
    gx = np.linspace(x.min(), x.max(), n)
    gy = np.linspace(y.min(), y.max(), n)
    hx, hy = bandwidth_nrd(x) / 4, bandwidth_nrd(y) / 4
    z = np.zeros((n, n))
    for start in range(0, len(x), chunk):
        kx = np.exp(-0.5 * ((gx[:, None] - x[start : start + chunk]) / hx) ** 2)
        ky = np.exp(-0.5 * ((gy[:, None] - y[start : start + chunk]) / hy) ** 2)
        z += kx @ ky.T
    return gx, gy, z / (2 * np.pi * len(x) * hx * hy)


def point_density(x, y, n=256):
    """Return each point's density, as the lab's R ``get_density`` does.

    Each point takes the value of the grid node at or below it on both axes,
    as R's ``findInterval`` picks it. There is no interpolation.

    Parameters
    ----------
    x : numpy.ndarray
        Point x coordinates.
    y : numpy.ndarray
        Point y coordinates.
    n : int
        Grid points per axis passed to `kde2d`.

    Returns
    -------
    numpy.ndarray
        One density per point.
    """
    gx, gy, z = kde2d(x, y, n)
    ix = np.searchsorted(gx, x, side="right") - 1
    iy = np.searchsorted(gy, y, side="right") - 1
    return z[ix, iy]


def heatscatter(ax, x, y, size=1.0):
    """Draw points colored by `point_density` on a reversed Spectral ramp.

    Points keep their input order, as ggplot2 draws them. Color scales
    linearly from the lowest to the highest density, as
    ``scale_color_gradientn`` does.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : numpy.ndarray
        Point x coordinates.
    y : numpy.ndarray
        Point y coordinates.
    size : float
        Marker area in pt².

    Returns
    -------
    matplotlib.collections.PathCollection
        The drawn points.
    """
    # A constant axis has zero bandwidth, where R's kde2d fails. Color every
    # point alike there instead of drawing NaN colors, which are invisible.
    with np.errstate(divide="ignore", invalid="ignore"):
        density = point_density(x, y)
    if not np.isfinite(density).all():
        density = np.zeros_like(x)
    return ax.scatter(
        x, y, c=density, cmap="Spectral_r", s=size, linewidths=0, rasterized=True
    )
