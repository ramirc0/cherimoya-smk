# cherimoya-smk

Snakemake workflow around the `cherimoya` PyTorch library. It replaces
`cherimoya pipeline` with a file-based DAG: per-sample fan-out from a sample
sheet, per-fold CV, SLURM support. It never reimplements the model. Each stage is
a flag-driven script in `workflow/scripts/` that imports the library. User-facing
docs are in `README.md`; this file covers what an agent needs to edit safely.

## Commands

```bash
conda activate snakemake                              # launcher env (has the SLURM plugin)
snakemake -n -p --profile profiles/local              # dry run; builds the full DAG
snakemake --profile profiles/local                    # local
snakemake --profile profiles/slurm                    # SLURM; fit/evaluate -> gpuh200
snakemake <target> --profile profiles/local --config samples=... run_id=...
python workflow/scripts/make_folds.py                 # fold JSONs, once per genome
.conda/<hash>_/bin/python -m pytest                   # 53 pass, 1 skip; the env built from cherimoya.yaml
.conda/<hash>_/bin/python -m pytest -m slow           # e2e; needs CHERIMOYA_SMK_SMOKE fixtures + GPU
```

Put the target **before** `--config`; otherwise Snakemake parses it as a config
entry. `config/config.yaml` points `samples:` at a local sheet that may not
exist. Pass `--config samples=config/samples.atac.tsv` for a dry run.

## DAG

```
macs3 (no peaks) ┐
signal ─ bam2bw ─┼─ negatives ─ fit* ─ evaluate* ─ performance.tsv + counts.tsv + per-model plots
peaks ─ prep_peaks ┘                        └─ gather_metrics ─ metrics.tsv ─ run-level plots
```

`*` fans out per `config["folds"]`. `bam2bw_control` builds the control track
for fit/evaluate. `count_fragments` feeds `gather_metrics`. `config_snapshot`
writes the resolved config to `results/<run_id>/config.snapshot.json`.
`model_summary` (CPU) loads each checkpoint and writes a torchinfo table with
the param count. Preprocessing is fold-agnostic. Outputs go under `results/<run_id>/<sample>/[fold_<k>/]`; logs and benchmarks mirror that
layout. Out of scope: attribution, seqlets, tomtom, modisco, marginalize.

## Where things live

- `workflow/rules/common.smk`: config, sample sheet validation, path helpers,
  `fit_flags`/`eval_flags`. Single gates: `STRANDED` (from
  `preprocess.unstranded`), `COVARIATES` (from `config["qc"]`), `_is_bigwig`
  (entry-point detection by extension).
- `preprocess.smk` (prep_peaks, macs3, bam2bw, bam2bw_control, count_fragments),
  `negatives.smk`, `train.smk` (fit, evaluate), `report.smk` (config_snapshot,
  gather_metrics, plots).
- `workflow/scripts/_style.py`: shared figure style. `save_figure` writes SVG
  and PNG together. Every `plot_*.py` MUST use it and follow the
  `matplotlib-style` skill.
- `config/config.<assay>.yaml.template`: tracked (atac, dnase, chipseq-tf).
  `config.yaml` and `samples.*.tsv` are gitignored.
- `resources/` (gitignored): `refs/<g>.{fa,fa.fai,chrom.sizes}`,
  `folds/<g>/fold_<k>.json` (`{train, valid, test}`), blacklist BED.
- `docs/run-paths.dot`: entry-path diagram. Rerender the SVG after editing it.
- `workflow/envs/cherimoya.yaml`: the one per-rule env, fully pinned (conda
  `name=version=build`, pip `==`, transitive deps included). A new dependency
  MUST go in with its exact version and any new transitive pins. There are no
  lock files.

## Behavior to preserve

- A sample enters where its inputs allow. A bigWig `signal` skips bam2bw. An
  empty `peaks` cell triggers macs3. bigWig without peaks is a `WorkflowError`.
  So is any bigWig signal or control under a stranded run (only bam2bw can
  build the (+, -) pair).
- Stranded runs train with `signal_groups=[2]` and `n_control_tracks=2`.
  The `control` column feeds both `macs3 -c` and the model.
- `n_peaks` is counted from the narrowPeak. `n_fragments` is a full signal
  scan, NaN for bigWig input. QC covariates are derived, never sheet columns.
- Fold JSONs are read in `params` lambdas at DAG-build time, so even `-n` fails
  without them.
- Script defaults MUST match `cherimoya_cli.defaults` at the **pinned** commit
  (in `.conda/`), not a dev checkout. `test_parser_drift.py` enforces this.
  `evaluate.py --counts_filename` is a local output path and is deliberately
  absent from the drift `EVAL_KEYS`.

## Rule conventions

Follow the [Nextstrain Snakemake style guide][sg]. Keep
`~/Projects/snakemake-template` in sync; it is the reference for these idioms.

- Raw `shell: r"""` blocks, one option per line.
- Log with `exec &> >(tee {log:q})`. Quote every interpolation with `:q`.
- Multi-value args are lists (`{params.flags:q}` quotes each token). Never build
  space-joined flag strings.
- Config reaches `shell` only through `params:` lambdas. Use `config[key]` for
  required keys, never bare `config.get(key)`.
- Every rule has `log:`, `benchmark:`, `conda:`. No `run:` blocks, no `message:`.

[sg]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html

## Gotchas

- No GPU on the login node. GPU is only via SLURM `gpuh200`. Set
  `fit.device`/`evaluate.device` to `cpu` for a CPU run.
- Never add `profiles/default`. Snakemake auto-loads it as a workflow profile
  that outranks `--profile` and strips GPU routing.
- `tasks_per_gpu: 0` on fit/evaluate suppresses `--ntasks-per-gpu`, which
  conflicts with `--cpus-per-task` on this cluster.
- SLURM resources are flat presets. To fix stragglers, raise the preset and
  rerun. Snakemake 9 has no `--stats`.
- The SLURM executor plugin lives in the launcher env, not the per-rule env.
- `scratch/` (gitignored) may hold `HANDOFF*.md` notes. `_archive/` holds
  shelved worktrees. Neither exists in a fresh clone.
