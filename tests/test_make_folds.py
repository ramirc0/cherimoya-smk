"""make_folds.py on synthetic chrom.sizes files."""

import json
import os

import make_folds
import pytest

FLY = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"]
YEAST = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII"]
YEAST += ["chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]


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


@pytest.mark.parametrize("main", [FLY, YEAST], ids=["fly", "yeast"])
def test_every_main_chromosome_is_tested_once(refs, main):
    refs("toy", [*main, "chrM", "chrUn_DS483562v1"])
    folds = make_folds.folds_for("toy")
    assert len(folds) == make_folds.N_FOLDS
    assert sorted(c for fold in folds for c in fold["test"]) == sorted(main)
    for fold in folds:
        assert sorted(fold["test"] + fold["valid"] + fold["train"]) == sorted(main)


def test_genome_with_fewer_chromosomes_than_folds_is_skipped(refs, capsys):
    refs("pombe", ["chr1", "chr2", "chr3"])
    make_folds.main()
    assert not (make_folds.OUT / "pombe").exists()
    assert "pombe" in capsys.readouterr().out
