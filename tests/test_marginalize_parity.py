"""Parity of marginalize.py with the marginalize step of `cherimoya pipeline`.

Skipped unless CHERIMOYA_SMK_MARGINALIZE points at a dir with genome.fa(.fai),
peaks.narrowPeak, negatives.bed, blacklist.bed, model.torch (a trained
checkpoint), fold.json (its CV fold) and motifs.meme. The report test needs
a GPU for device=cuda.
"""

import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import torch

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_MARGINALIZE")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_MARGINALIZE to run"),
]

DEVICE = os.environ.get("CHERIMOYA_SMK_SMOKE_DEVICE", "cuda")


def _official_json(fx, out):
    # The JSON `cherimoya pipeline` writes for step 5. Motifs go into the
    # negatives, with the blacklist excluded.
    from cherimoya_cli.defaults import (
        default_marginalize_parameters, default_pipeline_parameters)
    from cherimoya_cli.utils import _check_set, _extract_set

    parameters = copy.deepcopy(default_pipeline_parameters)
    parameters.update(sequences=str(fx / "genome.fa"),
        loci=[str(fx / "peaks.narrowPeak")], negatives=[str(fx / "negatives.bed")],
        exclusion_lists=[str(fx / "blacklist.bed")], model=str(fx / "model.torch"),
        motifs=str(fx / "motifs.meme"), device=DEVICE)
    marginalize_parameters = _extract_set(
        parameters, default_marginalize_parameters, "marginalize_parameters")
    marginalize_parameters["loci"] = (
        parameters["marginalize_parameters"]["loci"] or parameters["negatives"])
    _check_set(marginalize_parameters, "output_filename", str(out) + "/")
    _check_set(marginalize_parameters, "motifs", parameters["motifs"])

    name = out.parent / "marginalize.json"
    name.write_text(json.dumps(marginalize_parameters, sort_keys=True, indent=4))
    return name


def _ours_argv(fx, out):
    # marginalize.py as the rule calls it: chroms are the fold's train set.
    fold = json.loads((fx / "fold.json").read_text())
    return [str(SCRIPTS / "marginalize.py"),
        "-s", str(fx / "genome.fa"),
        "-l", str(fx / "negatives.bed"),
        "-e", str(fx / "blacklist.bed"),
        "-m", str(fx / "model.torch"),
        "-t", str(fx / "motifs.meme"),
        "-o", str(out),
        "--chroms", *fold["train"],
        "--device", DEVICE,
        "--verbose"]


def test_report_inputs_match_pipeline(tmp_path, monkeypatch):
    # Both steps end in one marginalization_report call; record its arguments.
    import argparse

    import bpnetlite.marginalize
    import marginalize
    from cherimoya_cli.commands import marginalize as official

    calls = []
    monkeypatch.setattr(bpnetlite.marginalize, "marginalization_report",
        lambda *args, **kwargs: calls.append((args, kwargs)))

    fx = Path(FIXTURES)
    official.run(argparse.Namespace(parameters=str(_official_json(fx, tmp_path / "o"))))
    monkeypatch.setattr(sys, "argv", _ours_argv(fx, tmp_path / "s"))
    marginalize.main()

    (o_model, o_motifs, o_X, o_dir), o_kwargs = calls[0]
    (s_model, s_motifs, s_X, s_dir), s_kwargs = calls[1]
    assert o_X.shape == (100, 4, 2114)
    assert torch.equal(s_X, o_X)
    assert s_motifs == o_motifs
    assert (o_dir, s_dir) == (str(tmp_path / "o") + "/", str(tmp_path / "s") + "/")
    assert s_kwargs == o_kwargs
    assert type(s_model) is type(o_model)
    o_state, s_state = o_model.state_dict(), s_model.state_dict()
    assert s_state.keys() == o_state.keys()
    assert all(torch.equal(s_state[k], o_state[k]) for k in o_state)


def test_report_matches_pipeline(tmp_path):
    fx = Path(FIXTURES)
    official, ours = tmp_path / "o", tmp_path / "s"
    # Inductor autotunes some kernels by timing, so fresh compiles can pick
    # different block sizes and move counts PNGs by a few pixels. One shared
    # cache makes both runs use the kernels the official run compiled.
    env = {**os.environ, "TORCHINDUCTOR_CACHE_DIR": str(tmp_path / "inductor")}
    cli = Path(sys.executable).parent / "cherimoya"
    subprocess.run([str(cli), "marginalize", "-p",
        str(_official_json(fx, official))], check=True, env=env)
    subprocess.run([sys.executable, *_ours_argv(fx, ours)], check=True, env=env)

    names = sorted(p.name for p in official.iterdir())
    assert sorted(p.name for p in ours.iterdir()) == names
    for name in names:
        if name.endswith(".png"):
            assert (ours / name).read_bytes() == (official / name).read_bytes(), name

    # Motifs sharing a consensus insert the same sequence, so they tie on the
    # counts ranking and GPU noise can swap their rows (1 run in 17).
    o_rows = (official / "marginalization.html").read_text().split("<tr")
    s_rows = (ours / "marginalization.html").read_text().split("<tr")
    assert len(o_rows) > 2
    assert sorted(s_rows) == sorted(o_rows)
    sequence = lambda row: re.findall(r"<td>(.*?)</td>", row)[1:2]
    assert [sequence(r) for r in s_rows] == [sequence(r) for r in o_rows]
