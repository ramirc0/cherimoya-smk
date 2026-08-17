"""Contig-filtering guard in negatives.py: drop/fail on peaks off the FASTA."""

import pytest

import negatives


def _write_fixtures(tmp_path, peak_chroms):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 10 + "\n")  # only chr1 exists
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("".join(f"{c}\t10\t20\n" for c in peak_chroms))
    return str(peaks), str(fasta)


def test_filter_drops_absent_contigs(tmp_path, capsys):
    peaks, fasta = _write_fixtures(tmp_path, ["chr1", "chrZ", "chr1"])
    assert negatives.resolve_present_chroms(peaks, fasta, "filter") == ["chr1"]
    err = capsys.readouterr().err
    assert "chrZ" in err and "WARNING" in err


def test_fail_policy_raises_on_absent_contig(tmp_path):
    peaks, fasta = _write_fixtures(tmp_path, ["chr1", "chrZ"])
    with pytest.raises(SystemExit) as exc:
        negatives.resolve_present_chroms(peaks, fasta, "fail")
    assert "chrZ" in str(exc.value)


def test_all_absent_always_raises(tmp_path):
    peaks, fasta = _write_fixtures(tmp_path, ["chrZ", "chrY"])
    with pytest.raises(SystemExit):
        negatives.resolve_present_chroms(peaks, fasta, "filter")


def test_no_missing_returns_all(tmp_path, capsys):
    peaks, fasta = _write_fixtures(tmp_path, ["chr1", "chr1"])
    assert negatives.resolve_present_chroms(peaks, fasta, "filter") == ["chr1"]
    assert capsys.readouterr().err == ""  # no warning when nothing dropped
