"""Keep PyTorch's CPU thread count through tangermeme's numba kernels."""

import os
import sys


def limit_cpu_threads(n_threads=None):
    """Cap PyTorch and numba at one CPU thread count.

    numba sizes its thread pool to the visible CPUs and ignores
    ``OMP_NUM_THREADS``. On its OpenMP layer, which shares PyTorch's libgomp,
    the first parallel kernel also resizes the calling thread's OpenMP team
    to that pool, so PyTorch's own setting is lost. Capping the pool at
    PyTorch's count keeps both at that count on every threading layer.

    Parameters
    ----------
    n_threads : int, optional
        Thread count for both. Default: PyTorch's current count, which
        follows ``OMP_NUM_THREADS``.

    Raises
    ------
    RuntimeError
        If numba is already imported, since it reads its pool size once, at
        import.
    """
    if "numba" in sys.modules:
        raise RuntimeError("limit_cpu_threads must run before numba is imported")

    import torch

    if n_threads is not None:
        torch.set_num_threads(n_threads)
    os.environ["NUMBA_NUM_THREADS"] = str(torch.get_num_threads())
