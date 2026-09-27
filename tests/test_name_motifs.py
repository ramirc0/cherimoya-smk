"""name_motifs.py rewrites JASPAR IDs and leaves everything else alone."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "name_motifs.py"

MEME = """MEME version 4

MOTIF MA0074.1 RXRA::VDR
letter-probability matrix: alength= 4 w= 1 nsites= 1 E= 0
 0.25 0.25 0.25 0.25
MOTIF MA1854.2 Etv1/4/5
MOTIF TCF-BHLH-SNAI_0
"""


def test_ids_become_name_accession(tmp_path):
    src, out = tmp_path / "in.meme", tmp_path / "out.meme"
    src.write_text(MEME)
    subprocess.run([sys.executable, str(SCRIPT), str(src), str(out)], check=True)
    assert out.read_text() == (MEME
        .replace("MA0074.1 RXRA::VDR", "RXRA+VDR_MA0074.1")
        .replace("MA1854.2 Etv1/4/5", "Etv1-4-5_MA1854.2"))
