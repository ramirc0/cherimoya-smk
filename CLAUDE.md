# cherimoya-smk

Snakemake workflow around the `cherimoya` PyTorch library. It is a full
replacement for `cherimoya_cli` and all its subcommands, with a file-based DAG:
per-sample fan-out from a sample sheet, per-fold CV, SLURM support. It never
reimplements the model. Each stage is a flag-driven script in `workflow/scripts/`
that calls the library (`cherimoya`, `tangermeme`) directly.

Workflow code MUST NOT import from `cherimoya_cli`, not even private helpers;
port the logic into the script instead, so any stage can take extra params.
`cherimoya_cli` is only the reference: tests use it for default drift and parity. User-facing
docs are in `README.md`; this file covers what an agent needs to edit safely.

## Commands

```bash
conda activate snakemake                              # launcher env (has the SLURM plugin)
snakemake -n -p --profile profiles/local              # dry run; builds the full DAG
snakemake --profile profiles/local                    # local
snakemake --profile profiles/slurm                    # SLURM; fit/evaluate -> gpuh200
snakemake <target> --profile profiles/local --config samples=... run_id=...
python workflow/scripts/make_folds.py                 # fold JSONs, once per genome
python workflow/scripts/name_motifs.py <in.meme> <out.meme>  # NAME_ACCESSION IDs for modisco report
.conda/<hash>_/bin/python -m pytest                   # 186 pass, 13 skip; the env built from cherimoya.yaml
.conda/<hash>_/bin/python -m pytest -m slow           # e2e; needs CHERIMOYA_SMK_SMOKE fixtures + GPU
CHERIMOYA_SMK_FIT=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_fit_parity.py  # vs cherimoya fit + evaluate; CPU
CHERIMOYA_SMK_ATTR=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_attribute_parity.py  # vs cherimoya attribute; GPU; atol 1e-4 counts, 2e-4 profile (official noise)
CHERIMOYA_SMK_SEQLETS=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_seqlets_parity.py  # vs cherimoya seqlets; CPU
CHERIMOYA_SMK_ANNOTATE=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_annotate_parity.py  # vs pipeline ttl step; CPU
CHERIMOYA_SMK_MODISCO=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_modisco_parity.py  # vs pipeline modisco; cascadelake, MEME tomtom on PATH
CHERIMOYA_SMK_MARGINALIZE=<dir> .conda/<hash>_/bin/python -m pytest -m slow tests/test_marginalize_parity.py  # vs pipeline marginalize; GPU
```

Put the target **before** `--config`; otherwise Snakemake parses it as a config
entry. `config/config.yaml` points `samples:` at a local sheet that may not
exist. Pass `--config samples=config/samples.atac.tsv` for a dry run.

## DAG

```
macs3 (no peaks) ┐
signal ─ bam2bw ─┼─ negatives ─ fit* ─ evaluate* ─ performance.tsv + counts.tsv + per-model plots
peaks ─ prep_peaks ┘                 │      └─ gather_metrics ─ metrics.tsv ─ run-level plots
                                     └─ attribute* ─ attributions.{ohe.npz,attr.npz,idxs.npy} ─┬─ seqlets* ─ seqlets.bed ─ annotate* ─ seqlets_annotated.bed + motif_seqlet_count.tsv
                                     │                                                        └─ modisco_motifs* ─ modisco_results.h5 ─ modisco_report* ─ modisco/report.html
                                     └─ marginalize* ─ marginalize/marginalization.html
```

`*` fans out per `config["folds"]`. `bam2bw_control` builds the control track
for fit/evaluate. `count_fragments` feeds `gather_metrics`. `config_snapshot`
writes the resolved config to `results/<run_id>/config.snapshot.json`.
`model_summary` (CPU) loads each checkpoint and writes a torchinfo table with
the param count.
`attribution_profile`, `seqlet_lengths`, `motif_counts` plot each
model; `motif_heatmap` plots the run. They are separate rules so a plot change
never reruns attribute. Preprocessing is fold-agnostic. Outputs go under `results/<run_id>/<sample>/[fold_<k>/]`; logs and benchmarks mirror that
layout.

## Where things live

- `workflow/rules/common.smk`: config, sample sheet validation, path helpers,
  `fit_flags`/`eval_flags`/`attr_flags`/`seqlet_flags`/`annot_flags`/`modisco_flags`/`modisco_report_flags`/`marginalize_flags`. Single gates: `STRANDED` (from
  `preprocess.unstranded`), `COVARIATES` (from `config["qc"]`), `_is_bigwig`
  (entry-point detection by extension).
- `preprocess.smk` (prep_peaks, macs3, bam2bw, bam2bw_control, count_fragments),
  `negatives.smk`, `train.smk` (fit, evaluate), `interpret.smk` (attribute, seqlets, annotate, modisco_motifs, modisco_report, marginalize),
  `report.smk` (config_snapshot, gather_metrics, plots).
- `workflow/scripts/_style.py`: shared figure style. `save_figure` writes SVG
  and PNG together. Every `plot_*.py` MUST use it and follow the
  `matplotlib-style` skill. `despine(ax, categorical_y=True)` is a local addition
  for heatmaps and horizontal bars.
- `config/config.<assay>.yaml.template`: tracked (atac, dnase, chipseq-tf).
  `config.yaml` and `samples.*.tsv` are gitignored.
- `resources/` (gitignored): `refs/<g>.{fa,fa.fai,chrom.sizes}`,
  `folds/<g>/fold_<k>.json` (`{train, valid, test}`), blacklist BED,
  `motifs/`: the JASPAR2026 MEME file (symlink) and its `_named` copy from
  `name_motifs.py`.
- `docs/run-paths.dot`: entry-path diagram. Rerender the SVG after editing it.
- `workflow/envs/cherimoya.yaml`: the main per-rule env, fully pinned (conda
  `name=version=build`, pip `==`, transitive deps included). A new dependency
  MUST go in with its exact version and any new transitive pins. There are no
  lock files. Editing it rebuilds the env and reruns every rule (software-env
  trigger).
- `workflow/envs/cherimoya-sm70.yaml`: the main env for V100 (`sm_70`) GPUs,
  with `profiles/slurm-v100`. torch comes from the CUDA 12.6 index. It is an
  exact `mamba env export` of a built env, not hand-written: rebuild it the
  same way. Selected with `--config conda_env=...`, since a CLI `--config`
  replaces a profile's `config:`.
- `workflow/envs/modisco_report.yaml`: `modisco_report` only. It adds MEME
  `tomtom` 5.5.9, whose `icu<76` conflicts with the main env. Its pip pins
  MUST match `cherimoya.yaml`.

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
  `evaluate.py --counts_filename` and `annotate.py --count_filename` are local
  output paths and are deliberately absent from the drift keys.
- fit.py calls `cherimoya.training.fit` (Lightning). Its validation set is the
  peaks plus every negative on `validation_chroms`. evaluate adds the fold's
  test-chrom negatives. Checkpoint selection stays peaks-only count Pearson.
  `evaluate.py --negatives`, like `--counts_filename`, is absent from the
  drift keys.
- The templates set `compile_mode: max-autotune-no-cudagraphs` for evaluate
  and marginalize. CUDA graph capture under `max-autotune` crashed
  intermittently, and the outputs are bitwise identical without it. The
  scripts keep the CLI default.
- `annotate` passes `--n_jobs {threads}` (profile `set-threads`), not the
  official `-1`. TomTom output is byte-identical across thread counts; memory
  grows about 115 MB per thread.
- `modisco_motifs.py` ports `modisco motifs` onto `modiscolite` (modisco 2.5.2
  owns the installed files, not modisco-lite). Its h5 is byte-identical to the
  official one and deterministic across runs and CPUs. `modisco_report` runs
  `modisco report` without `-l`, so it needs MEME `tomtom`, as the official run
  does. It goes through `modisco_report.py`, which turns off pandas 3 string
  inference: otherwise a pattern with fewer TomTom matches than `n_matches`
  crashes modisco 2.5.2 (the official pipeline too). Output is otherwise
  identical to `modisco report`. Its seqlet example picks depend on AVX512, so
  the SLURM preset pins it to cascadelake.
- `marginalize` inserts motifs into the **negatives**, minus the blacklist, as
  the pipeline does at the pin (its `loci` fall back to `negatives`). Chroms
  are the fold's train set (`training_chroms` in the pipeline). The report
  comes from `bpnetlite.marginalize.marginalization_report`. Motifs sharing a
  consensus tie on the ranking, so GPU noise can swap their HTML rows.
  The report is byte-identical when both runs share compiled kernels. Fresh
  compiles can differ, because inductor autotunes by timing.

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
- Concurrent runs in one workdir need distinct run_ids and staggered starts:
  simultaneous launches race on `.snakemake/iocache/latest.pkl`.
- Never add `profiles/default`. Snakemake auto-loads it as a workflow profile
  that outranks `--profile` and strips GPU routing.
- `tasks_per_gpu: 0` on fit/evaluate suppresses `--ntasks-per-gpu`, which
  conflicts with `--cpus-per-task` on this cluster.
- SLURM resources are flat presets. To fix stragglers, raise the preset and
  rerun. Snakemake 9 has no `--stats`.
- Never add `--signal` to fit's SLURM preset. On SIGUSR1, Lightning writes
  `hpc_ckpt_*` into the fold dir, and the next fit resumes from it and crashes.
- The SLURM executor plugin lives in the launcher env, not the per-rule env.
- `scratch/` (gitignored) may hold `HANDOFF*.md` notes. `_archive/` holds
  shelved worktrees. Neither exists in a fresh clone.
