"""tangermeme's numba kernels must not resize PyTorch's CPU thread pool."""

import os
import re
import subprocess
import sys

import pytest
from conftest import SCRIPTS

TORCH_AND_NUMBA_SCRIPTS = [
    "fit.py",
    "evaluate.py",
    "attribute.py",
    "marginalize.py",
    "seqlets.py",
]

needs_spare_cpus = pytest.mark.skipif(
    len(os.sched_getaffinity(0)) <= 3, reason="needs over 3 CPUs to show a resize"
)

PROBE = """
import sys
sys.path.insert(0, {scripts!r})
{setup}
import torch
from tangermeme.utils import one_hot_encode
one_hot_encode("ACGT" * 2**18, n_jobs=-1)
print(torch.get_num_threads())
"""


def torch_threads_after_numba(setup):
    """PyTorch's thread count after one parallel tangermeme kernel."""
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in ("OMP_NUM_THREADS", "NUMBA_NUM_THREADS")
    }
    # numba prefers TBB where the host has libtbb, and TBB leaves OpenMP alone.
    env["NUMBA_THREADING_LAYER"] = "omp"
    result = subprocess.run(
        [sys.executable, "-c", PROBE.format(scripts=str(SCRIPTS), setup=setup)],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return int(result.stdout.split()[-1])


@needs_spare_cpus
def test_numba_resizes_torch_threads_without_a_cap():
    assert torch_threads_after_numba("import torch; torch.set_num_threads(2)") > 2


@needs_spare_cpus
def test_limit_cpu_threads_holds_through_numba():
    setup = "from cpu_threads import limit_cpu_threads; limit_cpu_threads(2)"
    assert torch_threads_after_numba(setup) == 2


@needs_spare_cpus
def test_limit_cpu_threads_defaults_to_torch_count():
    setup = (
        "import torch; torch.set_num_threads(3)\n"
        "from cpu_threads import limit_cpu_threads; limit_cpu_threads()"
    )
    assert torch_threads_after_numba(setup) == 3


def test_limit_cpu_threads_refuses_after_numba_import():
    setup = (
        "import numba\nfrom cpu_threads import limit_cpu_threads; limit_cpu_threads(2)"
    )
    with pytest.raises(subprocess.CalledProcessError):
        torch_threads_after_numba(setup)


@pytest.mark.parametrize("script", TORCH_AND_NUMBA_SCRIPTS)
def test_script_limits_threads_before_importing_numba(script):
    source = (SCRIPTS / script).read_text()
    limit = source.index("limit_cpu_threads(")
    numba_import = re.search(
        r"^\s*(?:import|from) (?:lightning|cherimoya|tangermeme|bpnetlite)\b",
        source,
        re.MULTILINE,
    )
    assert limit < numba_import.start()
