"""Dry runs on a toy workdir: which sheet rows a run takes, and how their tracks enter it."""

import re
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = ROOT / ".pixi" / "envs" / "default" / "bin" / "snakemake"
FIT = "results/default/S/fold_0/S.torch"

pytestmark = pytest.mark.skipif(
    not SNAKEMAKE.exists(),
    reason="needs the launcher env (`pixi install` in the repo root)",
)


def _dry_run(tmp_path, *, stranded, **row):
    """Dry-run fit for one sheet row `row` on a toy genome, and return the output."""
    config = yaml.safe_load(
        (ROOT / "config" / "config.chipseq-tf.yaml.template").read_text()
    )
    config["samples"] = "samples.tsv"
    config["genomes"] = {
        "toy": {
            "fasta": "toy.fa",
            "fai": "toy.fa.fai",
            "chrom_sizes": "toy.chrom.sizes",
            "gsize": 1000,
        }
    }
    config["references"]["blacklist"] = None
    config["preprocess"]["unstranded"] = not stranded
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "config.yaml").write_text(yaml.safe_dump(config))
    folds = tmp_path / "resources" / "folds" / "toy"
    folds.mkdir(parents=True)
    (folds / "fold_0.json").write_text(
        '{"train": ["chr1"], "valid": ["chr2"], "test": ["chr3"]}'
    )

    row = {"sample_id": "S", "genome": "toy", **row}
    cols = ["sample_id", "signal", "genome", "control", "peaks"]
    sheet = ["\t".join(cols), "\t".join(row.get(c, "") for c in cols)]
    (tmp_path / "samples.tsv").write_text("\n".join(sheet) + "\n")
    cells = ",".join(row.get(c, "") for c in ("signal", "control", "peaks"))
    for name in ["toy.fa", "toy.fa.fai", "toy.chrom.sizes", *cells.split(",")]:
        if name.strip():
            (tmp_path / name.strip()).touch()

    result = subprocess.run(
        [
            SNAKEMAKE,
            "-n",
            "-p",
            "--snakefile",
            ROOT / "workflow" / "Snakefile",
            "--directory",
            tmp_path,
            FIT,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout + result.stderr


def _fit_command(out):
    """Return the fit command that a dry run printed."""
    start = out.index("scripts/fit.py")
    return out[start : out.index("\n\n", start)]


def test_stranded_run_takes_a_prebuilt_pair_in_order(tmp_path):
    code, out = _dry_run(
        tmp_path,
        stranded=True,
        signal="pos.bw,neg.bw",
        control="ctl_pos.bw, ctl_neg.bw",
        peaks="peaks.narrowPeak",
    )
    assert code == 0, out
    assert "rule bam2bw" not in out
    fit = _fit_command(out)
    assert "--stranded" in fit
    assert "-sig pos.bw neg.bw " in fit, fit
    assert "-c ctl_pos.bw ctl_neg.bw " in fit, fit


def test_stranded_run_takes_a_prebuilt_pair_without_control(tmp_path):
    code, out = _dry_run(
        tmp_path, stranded=True, signal="pos.bw,neg.bw", peaks="peaks.narrowPeak"
    )
    assert code == 0, out
    fit = _fit_command(out)
    assert "-sig pos.bw neg.bw " in fit, fit
    assert not re.search(r"\s-c\s", fit), fit


@pytest.mark.parametrize(
    "stranded,signal",
    [
        (True, "x.bw"),  # one bigWig cannot be split into strands
        (True, "x.bam,neg.bw"),  # a pair must be two bigWigs
        (False, "pos.bw,neg.bw"),  # an unstranded run takes one track
    ],
)
def test_rejects_a_signal_the_run_cannot_take(tmp_path, stranded, signal):
    code, out = _dry_run(
        tmp_path, stranded=stranded, signal=signal, peaks="peaks.narrowPeak"
    )
    assert code != 0
    assert "cannot take" in out, out


def test_bigwig_control_needs_peaks(tmp_path):
    # macs3 reads BAM/fragments only, so it cannot use a bigWig control.
    code, out = _dry_run(tmp_path, stranded=False, signal="x.bam", control="ctl.bw")
    assert code != 0
    assert "no peaks" in out, out
