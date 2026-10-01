"""negatives.py draws the same loci on every run unless the seed changes."""

import subprocess
import sys
from pathlib import Path

import negatives
import numpy

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "negatives.py"


def _run(tmp_path, name, *extra):
    rng = numpy.random.default_rng(1)
    fasta = tmp_path / "genome.fa"
    if not fasta.exists():
        seqs = {c: "".join(rng.choice(list("ACGT"), 300_000)) for c in ("chr1", "chr2")}
        fasta.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
        starts = rng.integers(5_000, 290_000, 40)
        (tmp_path / "peaks.bed").write_text(
            "".join(f"chr{i % 2 + 1}\t{s}\t{s + 500}\n" for i, s in enumerate(starts))
        )
    out = tmp_path / name
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-i",
            str(tmp_path / "peaks.bed"),
            "-f",
            str(fasta),
            "-o",
            str(out),
            *extra,
        ],
        check=True,
    )
    return out.read_text()


def test_default_seed():
    assert negatives.build_parser().get_default("random_state") == 243746692


def test_reruns_draw_the_same_loci(tmp_path):
    first = _run(tmp_path, "a.bed")
    assert first and first == _run(tmp_path, "b.bed")


def test_seed_changes_the_draw(tmp_path):
    assert _run(tmp_path, "a.bed") != _run(tmp_path, "b.bed", "--random_state", "1")
