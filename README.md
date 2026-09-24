# cherimoya-smk

A standalone Snakemake workflow that replaces the `cherimoya pipeline` command's
orchestration with a file-based DAG: incremental re-execution, multi-sample
fan-out from a sample sheet, and first-class SLURM support. It does **not**
replace the `cherimoya` library (the PyTorch model) -- each stage is a small
flag-driven script in `workflow/scripts/` that imports the library directly.

## Scope (this phase)

Preprocessing is per sample; `fit`/`evaluate` fan out per CV fold. See
`docs/run-paths.svg` for how a sample enters given its inputs.

```
macs3 (no peaks provided) ┐
signal ── bam2bw ─────────┼─ negatives ─ fit* ─ evaluate* ─ performance.tsv + counts.tsv
peaks  ── prep_peaks ─────┘                        └─ per-model + run-level QC (svgs, metrics)
```

`*` per fold. Attribution, seqlets, tomtom annotation, modisco, and marginalize
are out of scope for now.

## cherimoya pipeline parity

Each stage mirrors the matching step of `cherimoya pipeline` for the pinned
cherimoya commit, so results track the library rather than a reimplementation:

- Peak calling (`macs3`), `bam2bw`, and GC-matched negatives take the same flags
  and defaults; `fit`/`evaluate` mirror the library `fit`/`evaluate` commands.
- Stranded runs (`preprocess.unstranded: false`, the TF ChIP-seq template) build
  the `(+, -)` bigWig pair for both signal and control and train with
  `signal_groups=[2]` and `n_control_tracks=2`; unstranded runs use one track and
  no control track. The `control` column feeds both `macs3 -c` and the model.

Newer cherimoya adds a `min_total_steps` epoch floor, fixed `loss_weights`, and
fit-time RNG seeding with a model `random_state`. These are absent from the
pinned commit; wire them into `fit.py` when the env pin is bumped.

## Inputs

- `config/samples.tsv` — one row per sample: `sample_id`, `signal`, `genome`,
  optional `control`, optional `peaks`. `signal` is a BAM/fragments file
  (`bam2bw` + `macs3`) or an already-built bigWig (`.bw`/`.bigwig`, used as-is;
  requires provided peaks). An empty `peaks` cell triggers macs3 peak-calling.
- `config/config.yaml` — per-assembly `genomes` (fasta/fai/chrom_sizes/gsize),
  the CV `folds` to run, and every model/preprocess parameter, passed to the
  scripts as explicit CLI flags. Defaults mirror `cherimoya_cli/defaults.py`.
  Train/valid/test chroms come from per-genome fold JSONs, not config. Three
  assay-specific starting points are tracked: `config.atac.yaml.template`
  (unstranded, paired-end, Tn5 shift), `config.dnase.yaml.template` (unstranded,
  single-end, no shift), and `config.chipseq-tf.yaml.template` (stranded, no
  shift, input control).

## Environment

A single conda env (`workflow/envs/cherimoya.yaml`) pip-installs cherimoya from a
pinned git commit, which pulls torch (CUDA), tangermeme, macs3, and bam2bw.
Snakemake builds and caches it under `.conda/` when run with the profiles below.
Exact pins from a solved build are in `workflow/envs/cherimoya.{conda,pip}-lock.txt`
(recreate with `conda create --file cherimoya.conda-lock.txt` + `pip install -r
cherimoya.pip-lock.txt`).

## Running

```bash
# one-time: copy the template for your assay, then edit
cp config/config.atac.yaml.template config/config.yaml        # ATAC
# cp config/config.dnase.yaml.template config/config.yaml       # DNase-seq
# cp config/config.chipseq-tf.yaml.template config/config.yaml  # TF ChIP-seq

# one-time per genome: build fold JSONs from resources/refs/<g>.chrom.sizes
python workflow/scripts/make_folds.py

# local
snakemake --profile profiles/local

# SLURM (GPU rules -> gpuh200 with --gres=gpu:1)
snakemake --profile profiles/slurm
```

`snakemake -n -p --profile profiles/local` shows the plan without executing.
Because dependencies are file-based, editing one sample's bigWig re-runs only
that sample's `fit`/`evaluate`; everything else stays cached.

The SLURM profile needs `snakemake-executor-plugin-slurm` in the environment
that launches snakemake (not in the per-rule cherimoya env). Set an account/QOS
in `profiles/slurm/config.yaml` if your cluster requires them.

## Tests

```bash
pytest              # parser-drift guards + argparse smoke tests
pytest -m slow      # end-to-end fit->evaluate (needs CHERIMOYA_SMK_SMOKE fixtures + GPU)
```
