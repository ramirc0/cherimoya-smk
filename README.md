# cherimoya-smk

A Snakemake workflow that trains and evaluates
[cherimoya](https://cherimoya.readthedocs.io) models across many samples. It
replaces the orchestration in `cherimoya pipeline` with a file-based DAG. Only
stale outputs re-run, samples fan out from a sheet, and each CV fold is its
own job. A SLURM profile is included. The model itself stays in the library.
Each stage is a small script in `workflow/scripts/` that imports it.

```
macs3 (no peaks provided) ┐
signal ── bam2bw ─────────┼─ negatives ─ fit* ─ evaluate* ─ performance + counts
peaks  ── prep_peaks ─────┘                 │      └─ per-model and run-level QC plots
                                            ├─ attribute* ─┬─ seqlets* ─ annotate*
                                            │              └─ modisco_motifs* ─ modisco_report*
                                            └─ marginalize*
```

`*` runs once per CV fold in `folds`.

## Quick start

```bash
cp config/config.atac.yaml.template config/config.yaml   # or .dnase / .chipseq-tf
pixi shell                                               # launcher env
python workflow/scripts/make_folds.py                    # once per genome
snakemake -n -p --profile profiles/local                 # dry run
snakemake --profile profiles/local                       # run locally
snakemake --profile profiles/slurm                       # run on SLURM
```

The launcher env (`pixi.toml` at the root) has Snakemake, the pixi
software-deployment plugin and the SLURM executor plugins. Rule envs live in a
second pixi workspace, `workflow/envs/`. pixi must be on `PATH`. Each job
activates its env from `workflow/envs/.pixi/envs/` and installs it on first use.

## Inputs

**Sample sheet** (TSV, path set by `samples:` in config). One row per sample:

| column      | required | meaning |
|-------------|----------|---------|
| `sample_id` | yes      | output name |
| `signal`    | yes      | BAM, fragments file, or a pre-built bigWig (`.bw`/`.bigwig`) |
| `genome`    | yes      | a key under `genomes:` in config |
| `control`   | no       | control BAM. Passed to `macs3 -c` and fed to the model |
| `peaks`     | no       | narrowPeak. If empty, macs3 calls peaks from the BAM |

A bigWig signal skips `bam2bw` but can't be used for peak calling, so it needs
a `peaks` file. See `docs/run-paths.svg` for every entry path.

**Config** (`config/config.yaml`, gitignored). Start from one of the three
templates:

| template | tracks | shift | control |
|----------|--------|-------|---------|
| `config.atac.yaml.template` | unstranded, paired-end | Tn5 +4/-4 | none |
| `config.dnase.yaml.template` | unstranded, single-end | none | none |
| `config.chipseq-tf.yaml.template` | stranded (+, -) | none | input |

Model and preprocessing defaults match `cherimoya_cli/defaults.py` at the pinned
commit. `folds:` picks which CV folds to run, e.g. `[0]` or `[0, 1, 2, 3, 4]`.

**References**, per genome under `resources/`:

- `refs/<g>.fa`, `refs/<g>.fa.fai`, `refs/<g>.chrom.sizes` (`cut -f1,2 <g>.fa.fai`)
- `folds/<g>/fold_{0..4}.json` with `{train, valid, test}` chromosome lists,
  written by `make_folds.py`. hg38 uses the published chrombpnet folds. Other
  genomes get a length-balanced split.

To add a genome, add a `genomes:` entry, provide the three ref files, rerun
`make_folds.py`, and use the name in the sheet.

## Outputs

Everything goes under `results/<run_id>/` (`run_id` defaults to `default`).
Logs and benchmarks use the same layout under `logs/` and `benchmarks/`.

- `<sample>/fold_<k>/`: model (`.torch`), `performance.tsv`, `counts.tsv`,
  training curve and count scatter plots, and `summary.txt` (torchinfo layer
  table with the total parameter count). `performance.tsv` also has
  `all_count_*` over peaks plus negatives, and `auroc` and `auprc` of peaks vs
  negatives. The negatives are those on the fold's test chroms.
- `<sample>/fold_<k>/` also: DeepLIFT/SHAP attributions
  (`attributions.{ohe.npz,attr.npz,idxs.npy}`) and `seqlets.bed` (chrom,
  start, end, attribution, p-value; sorted by attribution),
  `seqlets_annotated.bed` (each seqlet's nearest TomTom motif in
  `annotate.motifs` and its -log p-value), `motif_seqlet_count.tsv`
  (seqlets per motif), and plots of the mean attribution by position, the
  seqlet lengths, and the top motifs. TF-MoDISco patterns
  (`modisco_results.h5`) and their HTML report (`modisco/report.html`), with
  each pattern's top TomTom matches in `annotate.motifs`. The marginalization
  report (`marginalize/marginalization.html`): each `annotate.motifs`
  consensus inserted at the center of the first 100 training-chrom negatives,
  ranked by the mean change in predicted counts
- `report/`: `metrics.tsv` with an outlier flag per model, the performance
  distribution, count Pearson vs `n_peaks`/`n_fragments`, the outlier plot,
  and a heatmap of the top motifs by sample
- `config.snapshot.json`: the fully resolved config, including `--config`
  overrides. Rerun it with `--configfile results/<run_id>/config.snapshot.json`.

Plots are written as SVG and PNG. `snakemake --report report.html` bundles them.

Keep runs apart with `--config run_id=mytag`. Swap sample sheets with
`--config samples=config/samples.other.tsv`.

## Run part of the pipeline

Name a stage to stop there. Each stage also builds every stage above it.

| target | adds |
|---|---|
| `preprocess` | bigWigs, peaks and negatives per sample; `config.snapshot.json` |
| `train` | model, evaluation, plots and `summary.txt` per fold; `report/` metrics and plots |
| `attributions` | attributions and the attribution profile plot per fold |
| `motifs` | seqlets, TomTom annotation and their plots per fold; the motif heatmap |
| `modisco` | TF-MoDISco patterns and report per fold |
| `all` | marginalization per fold. This is the default |

Put the target before `--config`:

```bash
snakemake train --profile profiles/slurm --config run_id=mytag
```

Run `all` later on the same `run_id` and only the remaining stages run.

## SLURM

`profiles/slurm` sends `fit`, `evaluate`, `attribute` and `marginalize` to the `gpuh200` partition with one
GPU and everything else to CPU partitions. `modisco_motifs` takes 16 cores on
cascadelake nodes (1 to 5 h). `modisco_report` also runs on cascadelake: without
AVX512, the seqlet examples it picks differ. Resources are fixed per rule. If a
job runs out of memory or time, raise its value in the profile and rerun. Only
failed jobs re-run. Set your own `slurm_account` before using it.

`profiles/slurm-v100` is the same profile on the V100 `gpu` partition. V100s
are `sm_70`, which the main env's CUDA 13 torch cannot run, so pair it with the
CUDA 12.6 env: `--config pixi_env=cherimoya-sm70`. A cold
compile makes the first jobs slower than on an H200.

With no GPU, set `fit.device`, `evaluate.device`, `attribute.device` and
`marginalize.device` to `cpu`.

## Environment

`workflow/envs/pixi.toml` defines three rule envs, and `pixi.lock` pins
every conda and PyPI package. `cherimoya` is the main env. cherimoya comes from
a pinned git commit. torch (CUDA), tangermeme, macs3, bam2bw and pytest are
PyPI packages in the same env. `cherimoya-sm70` is the same env with torch from
the CUDA 12.6 index. `modisco_report` uses `modisco-report`: the MEME-suite
`tomtom` it calls needs an `icu` that conflicts with the main env. Rules run
with `locked=True`, so a manifest edit needs `pixi lock` in `workflow/envs/`.
Any change to the manifest or lock reruns every rule (software-env trigger).

## Tests

Run from the `cherimoya` env (`pixi shell -e cherimoya` in `workflow/envs/`):

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

The parity tests run the official `cherimoya` command and the workflow script
on the same inputs. Their docstrings list the fixture files.
