"""End-to-end fit -> evaluate smoke test.

Skipped unless CHERIMOYA_SMK_SMOKE points at a dir with genome.fa(.fai),
peaks.narrowPeak, negatives.bed, signal.bw. Needs a GPU for device=cuda.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_SMOKE")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_SMOKE to run"),
]

DEVICE = os.environ.get("CHERIMOYA_SMK_SMOKE_DEVICE", "cuda")


def test_fit_then_evaluate(tmp_path):
    fx = Path(FIXTURES)
    prefix = tmp_path / "smoke"

    subprocess.run([
        sys.executable, str(SCRIPTS / "fit.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "peaks.narrowPeak"),
        "-neg", str(fx / "negatives.bed"),
        "-sig", str(fx / "signal.bw"),
        "-o", str(prefix),
        "--max_epochs", "1", "--n_warmup_epochs", "0", "--batch_size", "8",
        "--device", DEVICE,
    ], check=True)

    assert prefix.with_suffix(".torch").exists()

    perf = tmp_path / "smoke.performance.tsv"
    subprocess.run([
        sys.executable, str(SCRIPTS / "evaluate.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "peaks.narrowPeak"),
        "-sig", str(fx / "signal.bw"),
        "-m", str(prefix.with_suffix(".torch")),
        "-o", str(perf),
        "--device", DEVICE,
    ], check=True)

    header = perf.read_text().splitlines()[0].split("\t")
    assert len(header) == 7
