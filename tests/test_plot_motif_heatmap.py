"""plot_motif_heatmap.py on synthetic runs, to check it scales to 100+ samples."""

import subprocess
import sys
from pathlib import Path

import numpy
import polars as pl

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "plot_motif_heatmap.py"


def synthetic_counts(root, n_samples, names, n_programs=6, seed=0):
    """Write motif_seqlet_count.tsv per sample under root/<sample>/fold_0/.

    Each sample mixes a background with one of `n_programs` motif programs,
    so clustering has groups to recover. Returns the paths written.
    """
    rng = numpy.random.default_rng(seed)
    background = rng.dirichlet(numpy.full(len(names), 0.3))
    programs = [rng.choice(len(names), 5, replace=False) for _ in range(n_programs)]
    paths = []
    for i in range(n_samples):
        p = background.copy()
        p[programs[i % n_programs]] += rng.uniform(0.02, 0.15, 5)
        counts = rng.multinomial(rng.integers(20_000, 200_000), p / p.sum())
        keep = counts > 0
        df = pl.DataFrame({"motifs": numpy.array(names)[keep], "count": counts[keep]})
        path = root / f"S{i:03d}" / "fold_0" / f"S{i:03d}.motif_seqlet_count.tsv"
        path.parent.mkdir(parents=True)
        df.sort("count", descending=True).write_csv(path, separator="\t")
        paths.append(path)
    return paths


def _run(paths, out):
    subprocess.run([sys.executable, str(SCRIPT), "-i", *map(str, paths), "-o", str(out)],
        check=True)


def test_heatmap_scales_to_many_samples(tmp_path):
    names = [f"MA{i:04d}.1 TF{i}" for i in range(500)]
    _run(synthetic_counts(tmp_path, 120, names), tmp_path / "h.svg")
    assert (tmp_path / "h.svg").exists() and (tmp_path / "h.png").exists()


def test_heatmap_single_sample(tmp_path):
    names = [f"MA{i:04d}.1 TF{i}" for i in range(50)]
    _run(synthetic_counts(tmp_path, 1, names), tmp_path / "h.svg")
    assert (tmp_path / "h.png").exists()
