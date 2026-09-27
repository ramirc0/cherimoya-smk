"""Parity of attribute.py with the attribute step of `cherimoya pipeline`.

Skipped unless CHERIMOYA_SMK_ATTR points at a dir with genome.fa(.fai),
peaks.narrowPeak, model.torch (a trained checkpoint), and fold.json (its CV
fold). Needs a GPU for device=cuda.
"""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_ATTR")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_ATTR to run"),
]

DEVICE = os.environ.get("CHERIMOYA_SMK_SMOKE_DEVICE", "cuda")


def _official(fx, out):
    # The JSON `cherimoya pipeline` writes for step 2, then its CLI.
    from cherimoya_cli.defaults import (
        default_attribute_parameters, default_pipeline_parameters)
    from cherimoya_cli.utils import _check_set, _extract_set

    parameters = copy.deepcopy(default_pipeline_parameters)
    parameters.update(sequences=str(fx / "genome.fa"),
        loci=[str(fx / "peaks.narrowPeak")], model=str(fx / "model.torch"),
        device=DEVICE)
    attribute_parameters = _extract_set(
        parameters, default_attribute_parameters, "attribute_parameters")
    _check_set(attribute_parameters, "ohe_filename", str(out / "o.ohe.npz"))
    _check_set(attribute_parameters, "attr_filename", str(out / "o.attr.npz"))
    _check_set(attribute_parameters, "idx_filename", str(out / "o.idxs.npy"))

    name = out / "attribute.json"
    name.write_text(json.dumps(attribute_parameters, sort_keys=True, indent=4))
    cli = Path(sys.executable).parent / "cherimoya"
    subprocess.run([str(cli), "attribute", "-p", str(name)], check=True)


def _ours(fx, out):
    # attribute.py as the rule calls it: chroms are the fold's train + valid.
    fold = json.loads((fx / "fold.json").read_text())
    subprocess.run([
        sys.executable, str(SCRIPTS / "attribute.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "peaks.narrowPeak"),
        "-m", str(fx / "model.torch"),
        "--chroms", *fold["train"], *fold["valid"],
        "--device", DEVICE,
        "--ohe_filename", str(out / "s.ohe.npz"),
        "--attr_filename", str(out / "s.attr.npz"),
        "--idx_filename", str(out / "s.idxs.npy"),
    ], check=True)


def test_attribute_matches_pipeline(tmp_path):
    fx = Path(FIXTURES)
    _official(fx, tmp_path)
    _ours(fx, tmp_path)

    idxs = numpy.load(tmp_path / "o.idxs.npy")
    numpy.testing.assert_array_equal(numpy.load(tmp_path / "s.idxs.npy"), idxs)

    ohe = numpy.load(tmp_path / "o.ohe.npz")["arr_0"]
    numpy.testing.assert_array_equal(numpy.load(tmp_path / "s.ohe.npz")["arr_0"], ohe)
    assert ohe.shape[0] == idxs.sum() > 0

    attr = numpy.load(tmp_path / "o.attr.npz")["arr_0"]
    ours = numpy.load(tmp_path / "s.attr.npz")["arr_0"]
    print("max |attr diff|:", numpy.abs(ours - attr).max())
    numpy.testing.assert_allclose(ours, attr, rtol=1e-4, atol=1e-6)
