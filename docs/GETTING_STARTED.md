# Getting started

This guide covers a run from setup to outputs. The [README](../README.md) has the overview.

## Install

Install [pixi](https://pixi.sh) and put it on `PATH`. Then install the launcher env from the repository root:

```bash
pixi install
```

The launcher (`pixi.toml` at the root) has Snakemake, the pixi software-deployment plugin and the SLURM executor plugins. The rule envs live in a second pixi workspace, `workflow/envs/`. Each job activates its env from `workflow/envs/.pixi/envs/`. An `onstart` hook installs the selected env and `modisco-report` before the first job, because parallel SLURM jobs would otherwise race to install them. To install them ahead of time:

```bash
pixi install --manifest-path workflow/envs/pixi.toml --locked -e cherimoya -e modisco-report
```

## Set up a run

Each run keeps its inputs in `config/<run_id>/`. A run_id is the start date, then the assay and what varies, such as `2026-11-11-k562-atac`.

```bash
mkdir -p config/<run_id>
cp config/config.atac.yaml.template config/<run_id>/config.yaml   # or .dnase / .chipseq-tf
cp config/samples.tsv.template config/<run_id>/samples.tsv
```

In `config.yaml`, set `run_id: <run_id>` and `samples: config/<run_id>/samples.tsv`. Only the templates in `config/` are tracked.

### Sample sheet

The sheet is a TSV with one row per sample:

| column      | required | meaning |
|-------------|----------|---------|
| `sample_id` | yes      | output name |
| `signal`    | yes      | BAM, fragments file, or a prebuilt bigWig (`.bw`/`.bigwig`). A stranded run takes prebuilt tracks as a (+, -) pair, `plus.bw,minus.bw` |
| `genome`    | yes      | a key under `genomes:` in config |
| `control`   | no       | control BAM, or prebuilt bigWigs as for `signal`. Fed to the model. A BAM also goes to `macs3 -c` |
| `peaks`     | no       | narrowPeak. If empty, macs3 calls peaks from the BAM |

A bigWig signal or control skips `bam2bw` but can't be used for peak calling, so its row needs a `peaks` file.

### Config

Start from one of the three templates:

| template | tracks | shift | control |
|----------|--------|-------|---------|
| `config.atac.yaml.template` | unstranded, paired-end | Tn5 +4/-4 | none |
| `config.dnase.yaml.template` | unstranded, single-end | none | none |
| `config.chipseq-tf.yaml.template` | stranded (+, -) | none | input |

Model and preprocessing defaults match `cherimoya_cli/defaults.py` at the pinned cherimoya commit. `folds:` picks which CV folds to run, such as `[0]` or `[0, 1, 2, 3, 4]`. `random_state` seeds the negatives, fit, attribute and marginalize.

### References

Each genome needs these files under `resources/`:

- `refs/<g>.fa`, `refs/<g>.fa.fai`, `refs/<g>.chrom.sizes` (`cut -f1,2 <g>.fa.fai`)
- `folds/<g>/fold_{0..4}.json` with `{train, valid, test}` chromosome lists. Write them once per genome with `python workflow/scripts/make_folds.py`. hg38 uses the published chrombpnet folds and mm10 the lab's folds. Other genomes get a length-balanced split. A genome with fewer than five main chromosomes, such as S. pombe, gets no folds.

To add a genome, add a `genomes:` entry, provide the three ref files, rerun `make_folds.py`, and use the name in the sheet. An entry may also set `blacklist`, a BED of regions that fit, evaluate, attribute and marginalize exclude. Each genome carries its own, so one run can mix genomes.

## Run

Work in the launcher env (`pixi shell`, or prefix each command with `pixi run`). Dry-run first, then run locally or on SLURM:

```bash
snakemake -n -p --profile profiles/local --configfile config/<run_id>/config.yaml
snakemake --profile profiles/local --configfile config/<run_id>/config.yaml
snakemake --profile profiles/slurm --configfile config/<run_id>/config.yaml
```

Every run names its config. There is no `config/config.yaml` on purpose, so a run without `--configfile` fails instead of falling back to a stale default.

### Run part of the pipeline

Name a stage to stop there. Each stage also builds every stage above it.

| target | adds |
|---|---|
| `preprocess` | bigWigs, peaks, negatives and fragment counts per sample; `config.snapshot.json` |
| `train` | model, evaluation, plots and `summary.txt` per fold; `metrics.tsv` and the run-level performance plots |
| `attributions` | attributions and the attribution profile plot per fold |
| `motifs` | seqlets, TomTom annotation and their plots per fold; the motif heatmap |
| `modisco` | TF-MoDISco patterns and report per fold |
| `all` | marginalization per fold. This is the default |

Put the target before `--config` or `--configfile`. Otherwise, Snakemake parses it as a config entry:

```bash
snakemake train --profile profiles/slurm --configfile config/<run_id>/config.yaml
```

Later, run `all` with the same config. Only the missing stages run. Swap sample sheets with `--config samples=config/<run_id>/samples.other.tsv`.

## Outputs

Everything goes under `results/<run_id>/`. Logs and benchmarks use the same layout under `logs/` and `benchmarks/`.

- `<sample>/fold_<k>/`: model (`.torch`), `performance.tsv`, `counts.tsv`, training curve and count scatter plots, and `summary.txt` (torchinfo layer table with the total parameter count). `performance.tsv` also has `all_count_*` over peaks plus negatives, and `auroc` and `auprc` of peaks vs negatives. The negatives are those on the fold's test chroms.
- `<sample>/fold_<k>/` also holds the interpretation outputs:
  - DeepLIFT/SHAP attributions: `attributions.{ohe.npz,attr.npz,idxs.npy}`.
  - `seqlets.bed`: chrom, start, end, attribution and p-value, sorted by attribution.
  - `seqlets_annotated.bed`: each seqlet's nearest TomTom motif in `annotate.motifs`, with its -log p-value. `motif_seqlet_count.tsv` counts seqlets per motif.
  - Plots of the mean attribution by position, the seqlet lengths and the top motifs.
  - TF-MoDISco patterns (`modisco_results.h5`) and their HTML report (`modisco/report.html`). The report lists each pattern's top TomTom matches in `annotate.motifs`.
  - The marginalization report (`marginalize/marginalization.html`). It inserts each `annotate.motifs` consensus at the center of the first 100 training-chrom negatives and ranks the motifs by the mean change in predicted counts.
- `report/`: `metrics.tsv` with an outlier flag per model, the performance distribution and summary table, count Pearson vs `n_peaks`/`n_fragments`, the outlier plot, and a heatmap of the top motifs by sample. The distribution, covariate and outlier plots need `qc.min_models` models (samples x folds, default 10).
- `config.snapshot.json`: the fully resolved config, including `--config` overrides. Rerun it with `--configfile results/<run_id>/config.snapshot.json`.

Plots are written as SVG and PNG.

## Profiles

`profiles/local` runs every job on the current host. Without a GPU, set `fit.device`, `evaluate.device`, `attribute.device` and `marginalize.device` to `cpu`.

`profiles/slurm` sends `fit`, `evaluate`, `attribute` and `marginalize` to the `gpuh200` partition with one GPU, and everything else to CPU partitions. `modisco_motifs` takes 16 cores on cascadelake nodes (minutes per model). `modisco_report` also runs on cascadelake: without AVX512, the seqlet examples it picks differ. Resources are fixed per rule. If a job runs out of memory or time, raise its value in the profile and rerun. Only failed jobs rerun. Set your own `slurm_account` before using it.

`profiles/slurm-v100` is the same profile on the V100 `gpu` partition. V100s are `sm_70`, which the main env's CUDA 13 torch cannot run, so pair it with the CUDA 12.6 env: `--config pixi_env=cherimoya-sm70`. A cold compile makes the first jobs slower than on an H200.

## Environments

`workflow/envs/pixi.toml` defines three rule envs, and `pixi.lock` pins every conda and PyPI package. Releases newer than 7 days are excluded when the lock is solved.

- `cherimoya` is the main env. cherimoya, macs3 and modisco come from pinned git commits. torch (CUDA 13), tangermeme, figwig and the other Python packages are PyPI releases.
- `cherimoya-sm70` is the same env with torch from the CUDA 12.6 index, for V100s.
- `modisco-report` serves `modisco_report` alone. The MEME-suite `tomtom` it calls needs an `icu` that conflicts with the main env.

Rules run with `locked=True`, so a manifest edit needs `pixi lock` in `workflow/envs/`. Any change to the manifest or lock reruns every rule (software-env trigger).

## Tests

Run the tests from the `cherimoya` env (`pixi shell -e cherimoya` in `workflow/envs/`):

```bash
pytest            # drift guards against cherimoya defaults, CLI smoke tests
pytest -m slow    # end-to-end fit + evaluate, needs CHERIMOYA_SMK_SMOKE fixtures and a GPU
CHERIMOYA_SMK_FIT=<dir> pytest -m slow tests/test_fit_parity.py  # vs cherimoya fit + evaluate, CPU
CHERIMOYA_SMK_ATTR=<dir> pytest -m slow tests/test_attribute_parity.py    # vs cherimoya attribute, GPU
CHERIMOYA_SMK_SEQLETS=<dir> pytest -m slow tests/test_seqlets_parity.py   # vs cherimoya seqlets, CPU
CHERIMOYA_SMK_ANNOTATE=<dir> pytest -m slow tests/test_annotate_parity.py # vs cherimoya pipeline annotation, CPU
CHERIMOYA_SMK_MODISCO=<dir> pytest -m slow tests/test_modisco_parity.py   # vs cherimoya pipeline modisco, CPU (AVX512), MEME tomtom on PATH
CHERIMOYA_SMK_MARGINALIZE=<dir> pytest -m slow tests/test_marginalize_parity.py  # vs cherimoya pipeline marginalize, GPU
```

The parity tests run the official `cherimoya` command and the workflow script on the same inputs. Their docstrings list the fixture files.

`pixi run lint` (root workspace) runs `ruff check`, `ruff format --check` and `numpydoc lint` on the scripts. `pyproject.toml` holds their config. Neither tool parses `.smk` files, so keep the rule helpers numpydoc by hand.
