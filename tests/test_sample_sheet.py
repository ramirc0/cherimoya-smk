"""Dry runs on a toy workdir: which sheet rows a run takes, and how their tracks enter it."""

import re
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = ROOT / ".pixi" / "envs" / "default" / "bin" / "snakemake"

pytestmark = pytest.mark.skipif(
    not SNAKEMAKE.exists(),
    reason="needs the launcher env (`pixi install` in the repo root)",
)


def _dry_run(tmp_path, *rows, stranded=True, genomes=None, extra=None):
    """Dry-run fit and evaluate for sheet `rows`, and return the exit code and output.

    A row is sample `S` on genome `toy` unless it says otherwise. `genomes` maps
    each genome to extra keys of its config entry. `extra` adds top-level config.
    """
    rows = [{"sample_id": "S", "genome": "toy", **row} for row in rows]
    genomes = genomes or {"toy": {}}
    config = yaml.safe_load(
        (ROOT / "config" / "config.chipseq-tf.yaml.template").read_text()
    )
    config["samples"] = "samples.tsv"
    config["genomes"] = {
        g: {
            "fasta": f"{g}.fa",
            "fai": f"{g}.fa.fai",
            "chrom_sizes": f"{g}.chrom.sizes",
            "gsize": 1000,
            **keys,
        }
        for g, keys in genomes.items()
    }
    config["preprocess"]["unstranded"] = not stranded
    config.update(extra or {})
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "config.yaml").write_text(yaml.safe_dump(config))

    files = [p for entry in config["genomes"].values() for p in entry.values()]
    for g in genomes:
        folds = tmp_path / "resources" / "folds" / g
        folds.mkdir(parents=True)
        (folds / "fold_0.json").write_text(
            '{"train": ["chr1"], "valid": ["chr2"], "test": ["chr3"]}'
        )
    cols = ["sample_id", "signal", "genome", "control", "peaks"]
    sheet = ["\t".join(cols)]
    for row in rows:
        sheet.append("\t".join(row.get(c, "") for c in cols))
        for c in ("signal", "control", "peaks"):
            files += row.get(c, "").split(",")
    (tmp_path / "samples.tsv").write_text("\n".join(sheet) + "\n")
    for name in files:
        if isinstance(name, str) and name.strip():
            (tmp_path / name.strip()).touch()

    targets = [
        f"results/default/{r['sample_id']}/fold_0/{r['sample_id']}.performance.tsv"
        for r in rows
    ]
    result = subprocess.run(
        [
            SNAKEMAKE,
            "-n",
            "-p",
            "--snakefile",
            ROOT / "workflow" / "Snakefile",
            "--directory",
            tmp_path,
            *targets,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.returncode, result.stdout + result.stderr


def _command(out, script, sample="S"):
    """Return the `script` command that a dry run printed for `sample`."""
    for block in out.split("\n\n"):
        if f"scripts/{script}" in block and f"/{sample}/fold_0/{sample}" in block:
            return block
    raise AssertionError(f"no {script} command for {sample}:\n{out}")


def test_stranded_run_takes_a_prebuilt_pair_in_order(tmp_path):
    code, out = _dry_run(
        tmp_path,
        {
            "signal": "pos.bw,neg.bw",
            "control": "ctl_pos.bw, ctl_neg.bw",
            "peaks": "peaks.narrowPeak",
        },
    )
    assert code == 0, out
    assert "rule bam2bw" not in out
    fit = _command(out, "fit.py")
    assert "--stranded" in fit
    assert "-sig pos.bw neg.bw " in fit, fit
    assert "-c ctl_pos.bw ctl_neg.bw " in fit, fit


def test_stranded_run_takes_a_prebuilt_pair_without_control(tmp_path):
    code, out = _dry_run(
        tmp_path, {"signal": "pos.bw,neg.bw", "peaks": "peaks.narrowPeak"}
    )
    assert code == 0, out
    fit = _command(out, "fit.py")
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
        tmp_path, {"signal": signal, "peaks": "peaks.narrowPeak"}, stranded=stranded
    )
    assert code != 0
    assert "cannot take" in out, out


def test_bigwig_control_needs_peaks(tmp_path):
    # macs3 reads BAM/fragments only, so it cannot use a bigWig control.
    code, out = _dry_run(
        tmp_path, {"signal": "x.bam", "control": "ctl.bw"}, stranded=False
    )
    assert code != 0
    assert "no peaks" in out, out


def test_each_genome_excludes_its_own_blacklist(tmp_path):
    code, out = _dry_run(
        tmp_path,
        *(
            {
                "sample_id": g.upper(),
                "genome": g,
                "signal": f"{g}.bam",
                "peaks": f"{g}.bed",
            }
            for g in ("a", "b", "c")
        ),
        genomes={
            "a": {"blacklist": "a.bl.bed"},
            "b": {"blacklist": "b.bl.bed"},
            "c": {},
        },
    )
    assert code == 0, out
    for script in ("fit.py", "evaluate.py"):
        assert "--exclusion_lists a.bl.bed" in _command(out, script, "A")
        assert "--exclusion_lists b.bl.bed" in _command(out, script, "B")
        assert "--exclusion_lists" not in _command(out, script, "C")


def test_top_level_blacklist_is_rejected(tmp_path):
    code, out = _dry_run(
        tmp_path,
        {"signal": "x.bam", "peaks": "peaks.narrowPeak"},
        extra={"references": {"blacklist": "bl.bed"}},
    )
    assert code != 0
    assert "genomes.<genome>.blacklist" in out, out
