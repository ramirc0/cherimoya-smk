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
peaks  ── prep_peaks ─────┘                        └─ per-model and run-level QC plots
```

`*` runs once per CV fold. Attribution, seqlets, TomTom, modisco and
marginalize are not implemented.

## Quick start

```bash
cp config/config.atac.yaml.template config/config.yaml   # or .dnase / .chipseq-tf
python workflow/scripts/make_folds.py                    # once per genome
snakemake -n -p --profile profiles/local                 # dry run
snakemake --profile profiles/local                       # run locally
snakemake --profile profiles/slurm                       # run on SLURM
```

Launch from an env that has Snakemake and, for SLURM,
`snakemake-executor-plugin-slurm`. Snakemake builds the per-rule conda env
under `.conda/` on first run.

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
  table with the total parameter count)
- `report/`: `metrics.tsv` with an outlier flag per model, the performance
  distribution, count Pearson vs `n_peaks`/`n_fragments`, and the outlier plot
- `config.snapshot.json`: the fully resolved config, including `--config`
  overrides. Rerun it with `--configfile results/<run_id>/config.snapshot.json`.

Plots are written as SVG and PNG. `snakemake --report report.html` bundles them.

Keep runs apart with `--config run_id=mytag`. Swap sample sheets with
`--config samples=config/samples.other.tsv`.

## SLURM

`profiles/slurm` sends `fit` and `evaluate` to the `gpuh200` partition with one
GPU and everything else to CPU partitions. Resources are fixed per rule. If a
job runs out of memory or time, raise its value in the profile and rerun. Only
failed jobs re-run. Set your own `slurm_account` before using it.

With no GPU, set `fit.device` and `evaluate.device` to `cpu`.

## Environment

`workflow/envs/cherimoya.yaml` is the only env definition, and it is fully
pinned. Conda packages carry version and build, pip packages an exact `==`,
transitive dependencies included. cherimoya comes from a pinned git commit.
torch (CUDA), tangermeme, macs3, bam2bw and pytest are pip packages in the same
env. Any change to the file makes Snakemake build a fresh env under `.conda/`.

Newer cherimoya adds `min_total_steps`, fixed `loss_weights` and fit-time
seeding. The pinned commit has none of them. Wire them into `fit.py` when the
pin moves.

## Tests

Run from the built env under `.conda/`:

```bash
pytest            # drift guards against cherimoya defaults, CLI smoke tests
pytest -m slow    # end-to-end fit + evaluate, needs CHERIMOYA_SMK_SMOKE fixtures and a GPU
```
