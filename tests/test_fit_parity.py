"""Parity of fit.py and evaluate.py with `cherimoya fit`, bitwise on CPU.

Skipped unless CHERIMOYA_SMK_FIT points at a dir with genome.fa(.fai),
peaks.narrowPeak, negatives.bed, blacklist.bed, fold.json (a CV fold), and
either signal.bw (unstranded) or signal.+.bw, signal.-.bw, control.+.bw,
control.-.bw (stranded). Every subprocess runs eager (TORCH_COMPILE_DISABLE=1),
and both sides must run on the same node: Muon's bf16 math depends on the CPU.
"""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import evaluate

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_FIT")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_FIT to run"),
]

ENV = {**os.environ, "TORCH_COMPILE_DISABLE": "1"}
LOSS_WEIGHTS = [1.333, 0.274]
TIME_COLUMNS = {"Training Time", "Validation Time"}


def _tracks(fx):
    """Return the flat signal and control files and whether they are stranded."""
    if (fx / "signal.+.bw").exists():
        signals = [str(fx / "signal.+.bw"), str(fx / "signal.-.bw")]
        controls = [str(fx / "control.+.bw"), str(fx / "control.-.bw")]
        return signals, controls, True
    return [str(fx / "signal.bw")], None, False


def _official(fx, out, loss_weights):
    """Run `cherimoya fit` on its JSON and return the output prefix."""
    from cherimoya_cli.defaults import default_fit_parameters

    fold = json.loads((fx / "fold.json").read_text())
    signals, controls, stranded = _tracks(fx)
    parameters = copy.deepcopy(default_fit_parameters)
    parameters.update(sequences=str(fx / "genome.fa"),
        loci=str(fx / "peaks.narrowPeak"),
        negatives=str(fx / "negatives.bed"),
        signals=[signals] if stranded else signals,
        controls=[controls] if stranded else controls,
        exclusion_lists=[str(fx / "blacklist.bed")],
        training_chroms=fold["train"], validation_chroms=fold["valid"],
        device="cpu", max_epochs=2, min_total_steps=0, n_warmup_epochs=1,
        batch_size=16, random_state=0, loss_weights=loss_weights,
        name=str(out / "o"))

    name = out / "fit.json"
    name.write_text(json.dumps(parameters, sort_keys=True, indent=4))
    cli = Path(sys.executable).parent / "cherimoya"
    subprocess.run([str(cli), "fit", "-p", str(name)], check=True, env=ENV)
    return out / "o"


def _ours(fx, out, loss_weights, random_state=0):
    """Run fit.py as the rule calls it and return the output prefix."""
    fold = json.loads((fx / "fold.json").read_text())
    signals, controls, stranded = _tracks(fx)
    cmd = [
        sys.executable, str(SCRIPTS / "fit.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "peaks.narrowPeak"),
        "-neg", str(fx / "negatives.bed"),
        "-sig", *signals,
        "--exclusion_lists", str(fx / "blacklist.bed"),
        "--training_chroms", *fold["train"],
        "--validation_chroms", *fold["valid"],
        "--device", "cpu",
        "--max_epochs", "2",
        "--min_total_steps", "0",
        "--n_warmup_epochs", "1",
        "--batch_size", "16",
        "--random_state", str(random_state),
        "-o", str(out / "s"),
    ]
    if controls:
        cmd += ["-c", *controls]
    if stranded:
        cmd.append("--stranded")
    if loss_weights:
        cmd += ["--loss_weights", *map(str, loss_weights)]
    subprocess.run(cmd, check=True, env=ENV)
    return out / "s"


@pytest.fixture(scope="module")
def fx():
    """The fixture dir."""
    return Path(FIXTURES)


@pytest.fixture(scope="module")
def official(fx, tmp_path_factory):
    """Return a function giving the official prefix, one fit per loss_weights."""
    runs = {}

    def run(loss_weights):
        key = str(loss_weights)
        if key not in runs:
            runs[key] = _official(fx, tmp_path_factory.mktemp("official"),
                loss_weights)
        return runs[key]

    return run


def _load(path):
    """Load a checkpoint payload on CPU."""
    import torch

    return torch.load(path, map_location="cpu", weights_only=True)


def _log_rows(path):
    """Return a log's tab-split rows without the timing columns."""
    rows = [line.split("\t") for line in Path(path).read_text().splitlines()]
    keep = [i for i, name in enumerate(rows[0]) if name not in TIME_COLUMNS]
    assert len(keep) == len(rows[0]) - len(TIME_COLUMNS)
    return [[row[i] for i in keep] for row in rows]


@pytest.mark.parametrize("loss_weights", [None, LOSS_WEIGHTS])
def test_fit_matches_cli(fx, official, tmp_path, loss_weights):
    import torch

    theirs = official(loss_weights)
    ours = _ours(fx, tmp_path, loss_weights)

    for suffix in (".torch", ".final.torch"):
        o = _load(f"{theirs}{suffix}")
        s = _load(f"{ours}{suffix}")
        assert {**s["config"], "name": None} == {**o["config"], "name": None}
        assert s["state_dict"].keys() == o["state_dict"].keys()
        for key, tensor in o["state_dict"].items():
            assert torch.equal(s["state_dict"][key], tensor), (suffix, key)

    for suffix in (".log", ".detailed.log"):
        assert _log_rows(f"{ours}{suffix}") == _log_rows(f"{theirs}{suffix}")


def test_fit_seed_matters(fx, official, tmp_path):
    import torch

    o = _load(f"{official(None)}.torch")["state_dict"]
    s = _load(f"{_ours(fx, tmp_path, None, random_state=1)}.torch")["state_dict"]
    assert any(not torch.equal(s[key], o[key]) for key in o)


def test_evaluate_matches_cli(fx, official, tmp_path):
    theirs = official(None)
    fold = json.loads((fx / "fold.json").read_text())
    signals, controls, stranded = _tracks(fx)
    cmd = [
        sys.executable, str(SCRIPTS / "evaluate.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "peaks.narrowPeak"),
        "-sig", *signals,
        "--exclusion_lists", str(fx / "blacklist.bed"),
        "-m", f"{theirs}.torch",
        "-o", str(tmp_path / "s.performance.tsv"),
        "--chroms", *fold["valid"],
        "--batch_size", "16",
        "--no_compile",
        "--device", "cpu",
    ]
    if controls:
        cmd += ["-c", *controls]
    if stranded:
        cmd.append("--stranded")
    subprocess.run(cmd, check=True, env=ENV)

    def columns(path):
        header, *rows = [line.split("\t")
            for line in Path(path).read_text().splitlines()]
        return {name: [row[header.index(name)] for row in rows]
            for name in evaluate.MEASURE_NAMES}

    o = columns(f"{theirs}.performance.tsv")
    assert all(o.values())
    assert columns(tmp_path / "s.performance.tsv") == o
