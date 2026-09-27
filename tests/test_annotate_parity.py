"""Parity of annotate.py with the annotation step of `cherimoya pipeline`.

Skipped unless CHERIMOYA_SMK_ANNOTATE points at a dir with seqlets.bed,
genome.fa (+ .fai), and motifs.meme. CPU only.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
FIXTURES = os.environ.get("CHERIMOYA_SMK_ANNOTATE")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_ANNOTATE to run"),
]


def _official(fx, out):
    # Step 3.2 of `cherimoya pipeline` at the pin: ttl, then a per-motif count.
    import pandas
    from cherimoya_cli.defaults import default_annotation_parameters as a

    cmd = [str(Path(sys.executable).parent / "ttl"),
        "-f", str(fx / "genome.fa"),
        "-b", str(fx / "seqlets.bed"),
        "-s", str(a["n_score_bins"]),
        "-m", str(a["n_median_bins"]),
        "-a", str(a["n_target_bins"]),
        "-c", str(a["n_cache"]),
        "-j", str(a["n_jobs"]),
        "-t", str(fx / "motifs.meme")]
    assert a["reverse_complement"]
    with open(out / "o.bed", "w") as f:
        subprocess.run(cmd, check=True, stdout=f)

    annotated = pandas.read_csv(out / "o.bed", sep="\t", header=None,
        usecols=(3,), names=["motifs"])
    annotated.value_counts().to_csv(out / "o.tsv", sep="\t")


def _ours(fx, out):
    subprocess.run([
        sys.executable, str(SCRIPTS / "annotate.py"),
        "-s", str(fx / "genome.fa"),
        "-b", str(fx / "seqlets.bed"),
        "-t", str(fx / "motifs.meme"),
        "--output_filename", str(out / "s.bed"),
        "--count_filename", str(out / "s.tsv"),
    ], check=True)


def test_annotate_matches_pipeline(tmp_path):
    fx = Path(FIXTURES)
    _official(fx, tmp_path)
    _ours(fx, tmp_path)

    official = (tmp_path / "o.bed").read_bytes()
    print("seqlets:", official.count(b"\n"))
    assert official.count(b"\n") > 0
    assert (tmp_path / "s.bed").read_bytes() == official
    assert (tmp_path / "s.tsv").read_bytes() == (tmp_path / "o.tsv").read_bytes()
