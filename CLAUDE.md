# cherimoya-smk

A Snakemake workflow wrapping the `cherimoya` PyTorch library. It replaces the
`cherimoya pipeline` command with a file-based DAG: per-sample fan-out from a
sample sheet, per-fold cross-validation, incremental re-execution, SLURM support.
It does not reimplement the model; each stage is a small flag-driven script that
imports the library.

## DAG

Preprocessing is fold-agnostic (one bigWig/peaks/negatives per sample); `fit` and
`evaluate` fan out per CV fold. `config["folds"]` (e.g. `[0]` or `[0,1,2,3,4]`)
sets which folds run.

```
macs3 (only if a sample has no peaks) ┐
signal(bam) ── bam2bw ────────────────┼─ negatives ─ fit* ─ evaluate* ─ performance.tsv + counts.tsv
peaks (provided) ── prep_peaks ───────┘                        │
                                                               ├─ epochs.svg, count_scatter.svg  (per model)
                                                               └─ gather_metrics ─ metrics.tsv ─┬─ performance_distribution.svg
                                                                                                ├─ count_pearson_vs_<cov>.svg
                                                                                                └─ outliers.svg + worst_offenders.{tsv,html}
```

`*` fan out per fold. Out of scope: attribution, seqlets, tomtom, modisco,
marginalize.

## Layout

```
workflow/
  Snakefile              includes + rule all (per-model QC + run-level QC)
  rules/
    common.smk           config, sample sheet, genomes/folds, path constants, flag builders
    preprocess.smk       prep_peaks, macs3, bam2bw
    negatives.smk        GC-matched negatives
    train.smk            fit, evaluate (per fold; evaluate also emits counts + count_scatter)
    report.smk           gather_metrics, perf_vs_covariate, outliers, performance_distribution
  scripts/               fit.py, evaluate.py, negatives.py, make_folds.py,
                         gather_metrics.py, plot_*.py, _style.py
  envs/cherimoya.yaml    single conda env (pins cherimoya from a git commit)
config/
  config.yaml.template   tracked; copy to config.yaml (gitignored), then edit
  samples.*.tsv          sample sheet (gitignored)
profiles/{local,slurm}/  execution profiles (slurm-aicr is a reference only)
resources/               gitignored data: refs/<genome>.{fa,fa.fai,chrom.sizes},
                         folds/<genome>/fold_{0..4}.json, blacklist BED
tests/                   parser-drift + argparse smoke + negatives + slow e2e
scratch/                 HANDOFF*.md, RUN.md (run notes, session history)
```

## Multispecies + CV folds

Each sample picks its assembly via the sample sheet's `genome` column. Assemblies
are defined in `config["genomes"]` (`fasta`, `fai`, `chrom_sizes`, macs3 `gsize`).
Train/valid/test chrom splits live in per-genome fold JSONs
(`resources/folds/<genome>/fold_<k>.json`, schema `{train, valid, test}`),
generated once by `python workflow/scripts/make_folds.py` (hg38 = chrombpnet's
published folds; others = length-balanced greedy partition). `fit` reads
train+valid; `evaluate` reads test. There are no chrom lists in config.

Adding a genome: add a `genomes:` block, provide
`resources/refs/<g>.{fa,fa.fai,chrom.sizes}`, run `make_folds.py`, set the sheet
column.

## Running

```bash
cp config/config.yaml.template config/config.yaml    # one-time
python workflow/scripts/make_folds.py                # one-time per genome (needs chrom.sizes)
snakemake --profile profiles/local                   # local (device: cpu for CPU box)
snakemake --profile profiles/slurm                    # SLURM; GPU rules -> gpuh200
snakemake -n -p --profile profiles/local              # dry run
snakemake --report report.html                        # collect report() outputs
```

Override per run: `--config samples=config/samples.mohd.tsv run_id=mytag`.
Outputs group under `results/<run_id>/<sample>/fold_<fold>/`; logs and benchmarks
mirror it.

## Conventions

Rules follow the [Nextstrain Snakemake style guide][sg]. When editing rules:

- Raw triple-quoted `shell:` blocks (`r"""`), one option per line.
- Log via `exec &> >(tee {log:q})`; every interpolation quoted with `:q`.
- Multi-value args are lists (see `fit_flags`/`eval_flags` in common.smk), so
  `{params.flags:q}` quotes each token. Never build space-joined flag strings.
- Config mapped into `shell` through `params:` (lambdas), never interpolated
  directly. `config[key]` for required keys; never bare `config.get(key)`.
- Every rule has `log:`, `benchmark:`, `conda:`. No `run:` blocks, no `message:`.

The `snakemake-template` repo (`~/Projects/snakemake-template`) is the canonical
demonstrator of these idioms; keep the two in sync.

[sg]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html

## Environment

Single conda env at `workflow/envs/cherimoya.yaml` (cherimoya pip-installed from
a pinned commit; pulls torch/CUDA, tangermeme, macs3, bam2bw). Snakemake builds
it under `.conda/` when run with a profile. Locked pins in
`cherimoya.{conda,pip}-lock.txt`. The SLURM executor plugin must live in the env
that launches snakemake, not the per-rule env. The launcher env is the
`snakemake` conda env (`conda activate snakemake`); pytest with the library deps
runs in the built per-rule env under `.conda/`.

## Tests

```bash
pytest              # parser-drift guards + argparse smoke tests (44 pass, 1 skip)
pytest -m slow      # end-to-end fit->evaluate (needs CHERIMOYA_SMK_SMOKE fixtures + GPU)
```

`test_parser_drift.py` asserts script defaults still match `cherimoya_cli.defaults`.
`evaluate.py`'s `--counts_filename` is a local output path (not a library
hyperparameter), so it is deliberately absent from the drift `EVAL_KEYS`.

## Gotchas

- **No GPU on the dev/login node.** GPU is only via SLURM `gpuh200`. `fit`/
  `evaluate` set `device: cuda`; set to `cpu` in config for a CPU-only run.
- **Fold JSONs must exist before running.** `common.smk` reads them in `params`
  lambdas at DAG-build time (even for `-n`). Run `make_folds.py` first.
- **chrom.sizes are provided, not derived.** There is no `chrom_sizes` rule; each
  genome needs `resources/refs/<g>.chrom.sizes` (e.g. `cut -f1,2 <g>.fa.fai`).
- **Profile precedence.** Use `--profile profiles/{local,slurm}`. Never a
  `profiles/default` (Snakemake auto-loads it as a workflow profile that
  outranks `--profile` and would strip GPU routing).
- **`tasks_per_gpu: 0`** on fit/evaluate suppresses `--ntasks-per-gpu`, which
  conflicts with `--cpus-per-task` on this SLURM.
- **SLURM resources are flat presets, not scaled.** Run with defaults; to catch
  stragglers, bump the preset in `profiles/slurm` and rerun (only failed jobs
  re-run). No `--stats` in Snakemake 9.
- See `scratch/HANDOFF-2026-09-12.md` for the fullest session history.
```
