"""make_folds.py on synthetic chrom.sizes files."""

import json
import os

import make_folds
import pytest


@pytest.fixture
def refs(tmp_path, monkeypatch):
    """Point make_folds at a tmp resources dir, and return a chrom.sizes writer."""
    monkeypatch.setattr(make_folds, "REFS", tmp_path / "refs")
    monkeypatch.setattr(make_folds, "OUT", tmp_path / "folds")
    (tmp_path / "refs").mkdir()

    def write(genome, chroms):
        lines = [f"{c}\t{1_000_000 + 10_000 * i}" for i, c in enumerate(chroms)]
        (tmp_path / "refs" / f"{genome}.chrom.sizes").write_text(
            "\n".join(lines) + "\n"
        )

    return write


def test_rerun_keeps_unchanged_fold_files(refs):
    refs("toy", [f"chr{n}" for n in range(1, 8)])
    make_folds.main()
    path = make_folds.OUT / "toy" / "fold_0.json"
    os.utime(path, (0, 0))
    make_folds.main()
    assert path.stat().st_mtime == 0
    assert json.loads(path.read_text())["test"]


def test_mm10_uses_the_lab_folds(refs):
    refs("mm10", [*(f"chr{n}" for n in range(1, 20)), "chrX", "chrY", "chrM"])
    fold = make_folds.folds_for("mm10")[0]
    assert (fold["test"], fold["valid"]) == (
        ["chr1", "chr3", "chr6"],
        ["chr8", "chr16"],
    )
