"""Parity of the modisco rules with steps 4.1 and 4.2 of `cherimoya pipeline`.

Skipped unless CHERIMOYA_SMK_MODISCO points at a dir with attributions.ohe.npz,
attributions.attr.npz, and motifs.meme from the pipeline inputs, and
modisco_results.h5 and modisco/ from the workflow. CPU only. The motifs test
runs the official `modisco motifs` (about 1 h on 16 cores for atac). The report
test needs the official MEME `tomtom` on PATH and an AVX512 CPU, as the rule's
cascadelake preset: its seqlet example picks differ on broadwell.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy
import pytest

FIXTURES = os.environ.get("CHERIMOYA_SMK_MODISCO")
MODISCO = str(Path(sys.executable).parent / "modisco")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not FIXTURES, reason="set CHERIMOYA_SMK_MODISCO to run"),
]


def _h5_items(path):
    """Map every dataset and attribute path in an HDF5 file to its value."""
    import h5py

    items = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            for key, value in obj.attrs.items():
                items[f"{name}@{key}"] = numpy.asarray(value)
            if isinstance(obj, h5py.Dataset):
                items[name] = obj[()]
        visit("/", f)
        f.visititems(visit)
    return items


def _strip_images(html):
    return re.sub(r"data:image/png;base64,[A-Za-z0-9+/=]*", "", html)


def test_motifs_match_pipeline(tmp_path):
    # Step 4.1 at the pin: only -n is set; everything else is the CLI default.
    from cherimoya_cli.defaults import default_pipeline_parameters as p

    fx = Path(FIXTURES)
    subprocess.run([MODISCO, "motifs",
        "-s", str(fx / "attributions.ohe.npz"),
        "-a", str(fx / "attributions.attr.npz"),
        "-n", str(p["modisco_motifs_parameters"]["n_seqlets"]),
        "-o", str(tmp_path / "o.h5")], check=True)

    official = _h5_items(tmp_path / "o.h5")
    ours = _h5_items(fx / "modisco_results.h5")
    print("patterns:", sum(k.endswith("/sequence") for k in official))
    assert official.keys() == ours.keys()
    for key in official:
        assert numpy.array_equal(official[key], ours[key]), key


def test_report_matches_pipeline(tmp_path):
    # Step 4.2 at the pin: `modisco report -i -o -s ./ -m`, MEME tomtom.
    if not shutil.which("tomtom"):
        pytest.skip("MEME tomtom not on PATH")
    fx = Path(FIXTURES)
    subprocess.run([MODISCO, "report",
        "-i", str(fx / "modisco_results.h5"),
        "-o", str(tmp_path / "o"),
        "-s", "./",
        "-m", str(fx / "motifs.meme")], check=True)

    official = _strip_images((tmp_path / "o" / "report.html").read_text())
    ours = _strip_images((fx / "modisco" / "report.html").read_text())
    assert "<code>" in official
    assert ours == official
