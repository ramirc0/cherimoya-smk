"""Dry runs on a toy workdir: which sheet rows a run takes, and how their tracks enter it."""

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
    cells = [row.get(c, "") for c in ("signal", "control", "peaks")]
    for name in [
        "toy.fa",
        "toy.fa.fai",
        "toy.chrom.sizes",
        *",".join(cells).split(","),
    ]:
        if name:
            (tmp_path / name).touch()

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


def test_bigwig_control_needs_peaks(tmp_path):
    # macs3 reads BAM/fragments only, so it cannot use a bigWig control.
    code, out = _dry_run(tmp_path, stranded=False, signal="x.bam", control="ctl.bw")
    assert code != 0
    assert "no peaks" in out, out
