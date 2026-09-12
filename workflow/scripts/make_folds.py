"""Generate chrombpnet-style 5-fold CV chromosome splits, tailored per genome.

chrombpnet's convention (https://zenodo.org/records/7445373): each fold is
{"test", "valid", "train"}; across the 5 folds every chromosome is in `test`
exactly once, and the test sets are balanced by *length* (~20% of the genome
each), not by count. `valid` is a small (~7%) early-stopping holdout; `train` is
the rest. There is no official split for non-human genomes, so we reproduce the
same structure from each genome's chrom.sizes:

- hg38: chrombpnet's published folds verbatim (the established standard).
- others: greedy length-balanced 5-way `test` partition (each chrom once), plus a
  rotating `valid` (smallest chroms of the next fold's test bin, up to ~8% of the
  genome). Only main chromosomes (autosomes incl. chr2A/2B/chr1b/2b + X/Y) are
  split; scaffolds and chrM stay in `train`.

Writes resources/folds/<genome>/fold_{0..4}.json for every genome whose
chrom.sizes is present under resources/refs/.
"""

import json
import re
from pathlib import Path

REFS = Path("resources/refs")
OUT = Path("resources/folds")
N_FOLDS = 5
VALID_FRACTION = 0.08  # ~8% of the genome held for validation per fold
MAIN = re.compile(r"chr([0-9]+[AB]?|[0-9]+b|X|Y)")

# chrombpnet's published hg38 folds (zenodo 7445373), used verbatim.
HG38 = [
    {"test": ["chr1", "chr3", "chr6"], "valid": ["chr8", "chr20"]},
    {"test": ["chr2", "chr8", "chr9", "chr16"], "valid": ["chr12", "chr17"]},
    {"test": ["chr4", "chr11", "chr12", "chr15", "chrY"], "valid": ["chr22", "chr7"]},
    {"test": ["chr5", "chr10", "chr14", "chr18", "chr20", "chr22"], "valid": ["chr6", "chr21"]},
    {"test": ["chr7", "chr13", "chr17", "chr19", "chr21", "chrX"], "valid": ["chr10", "chr18"]},
]


def main_chroms(genome):
    """[(chrom, size)] for a genome's main chromosomes, from chrom.sizes."""
    out = []
    with (REFS / f"{genome}.chrom.sizes").open() as fh:
        for line in fh:
            name, size = line.split("\t")[:2]
            if MAIN.fullmatch(name):
                out.append((name, int(size)))
    return out


def balanced_test_bins(sizes):
    """Greedy length-balanced partition of chroms into N_FOLDS test sets."""
    bins, totals = [[] for _ in range(N_FOLDS)], [0] * N_FOLDS
    for chrom, size in sorted(sizes.items(), key=lambda kv: -kv[1]):
        j = min(range(N_FOLDS), key=lambda k: totals[k])
        bins[j].append(chrom)
        totals[j] += size
    return bins


def folds_for(genome):
    """The 5 {test, valid, train} folds for a genome."""
    sizes = dict(main_chroms(genome))
    all_chroms = list(sizes)
    if genome == "hg38":
        specs = HG38
    else:
        bins = balanced_test_bins(sizes)
        target = VALID_FRACTION * sum(sizes.values())
        specs = []
        for i in range(N_FOLDS):
            nxt = sorted(bins[(i + 1) % N_FOLDS], key=lambda c: sizes[c])
            valid, acc = [], 0
            for c in nxt:  # smallest chroms of the next bin until ~8%
                valid.append(c)
                acc += sizes[c]
                if acc >= target:
                    break
            specs.append({"test": bins[i], "valid": valid})

    folds = []
    for spec in specs:
        held = set(spec["test"]) | set(spec["valid"])
        folds.append({
            "test": spec["test"],
            "valid": spec["valid"],
            "train": [c for c in all_chroms if c not in held],
        })
    return folds


def main():
    """Write fold JSONs for every genome under resources/refs."""
    genomes = sorted(p.name[: -len(".chrom.sizes")]
                     for p in REFS.glob("*.chrom.sizes"))
    for genome in genomes:
        folds = folds_for(genome)
        d = OUT / genome
        d.mkdir(parents=True, exist_ok=True)
        for i, fold in enumerate(folds):
            (d / f"fold_{i}.json").write_text(json.dumps(fold, indent=2) + "\n")
        # sanity: every main chrom in test exactly once
        counts = {}
        for fold in folds:
            for c in fold["test"]:
                counts[c] = counts.get(c, 0) + 1
        bad = {c: n for c, n in counts.items() if n != 1}
        print(f"{genome}: {len(folds)} folds, {len(counts)} chroms"
              + (f"  WARN test-count!=1: {bad}" if bad else "  (each chrom test x1)"))


if __name__ == "__main__":
    main()
