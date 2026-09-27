"""Parity of seqlets.py with the seqlets step of `cherimoya pipeline`.

Skipped unless CHERIMOYA_SMK_SEQLETS points at a dir with peaks.narrowPeak,
fold.json (the attributed CV fold), and attribute.py outputs
attributions.{ohe.npz,attr.npz,idxs.npy}. CPU only.
"""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_SEQLETS")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_SEQLETS to run"),
]


def _chroms(fx):
    fold = json.loads((fx / "fold.json").read_text())
    return fold["train"] + fold["valid"]


def _official(fx, out):
    # The JSON `cherimoya pipeline` writes for step 3.1, then its CLI.
    from cherimoya_cli.defaults import (
        default_pipeline_parameters, default_seqlet_parameters)
    from cherimoya_cli.utils import _check_set, _extract_set

    parameters = copy.deepcopy(default_pipeline_parameters)
    parameters.update(loci=[str(fx / "peaks.narrowPeak")])
    seqlet_parameters = _extract_set(
        parameters, default_seqlet_parameters, "seqlet_parameters")
    _check_set(seqlet_parameters, "ohe_filename", str(fx / "attributions.ohe.npz"))
    _check_set(seqlet_parameters, "attr_filename", str(fx / "attributions.attr.npz"))
    _check_set(seqlet_parameters, "idx_filename", str(fx / "attributions.idxs.npy"))
    _check_set(seqlet_parameters, "output_filename", str(out / "o.bed"))
    _check_set(seqlet_parameters, "chroms", _chroms(fx))

    name = out / "seqlets.json"
    name.write_text(json.dumps(seqlet_parameters, sort_keys=True, indent=4))
    cli = Path(sys.executable).parent / "cherimoya"
    subprocess.run([str(cli), "seqlets", "-p", str(name)], check=True)


def _ours(fx, out):
    subprocess.run([
        sys.executable, str(SCRIPTS / "seqlets.py"),
        "-l", str(fx / "peaks.narrowPeak"),
        "--ohe_filename", str(fx / "attributions.ohe.npz"),
        "--attr_filename", str(fx / "attributions.attr.npz"),
        "--idx_filename", str(fx / "attributions.idxs.npy"),
        "--output_filename", str(out / "s.bed"),
        "--chroms", *_chroms(fx),
    ], check=True)


def test_seqlets_match_pipeline(tmp_path):
    fx = Path(FIXTURES)
    _official(fx, tmp_path)
    _ours(fx, tmp_path)

    official = (tmp_path / "o.bed").read_bytes()
    print("seqlets:", official.count(b"\n"))
    assert official.count(b"\n") > 0
    assert (tmp_path / "s.bed").read_bytes() == official
