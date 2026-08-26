# cherimoya-smk

A standalone Snakemake workflow that replaces the `cherimoya pipeline` command's
orchestration with a file-based DAG: incremental re-execution, multi-sample
fan-out from a sample sheet, and first-class SLURM support. It does **not**
replace the `cherimoya` library (the PyTorch model) -- each stage is a small
flag-driven script in `workflow/scripts/` that imports the library directly.

## Scope (this phase)

```
macs3 (only if a sample has no peaks) ┐
signal(bam) ── bam2bw ────────────────┼─ negatives ─ fit ─ evaluate ─ performance.tsv
peaks (provided) ── prep_peaks ───────┘
```

Attribution, seqlets, tomtom annotation, modisco, and marginalize are out of
scope for now.

## Inputs

- `config/samples.tsv` — one row per sample: `sample_id`, `signal` (BAM /
  fragments), optional `control`, optional `peaks`. An empty `peaks` cell makes
  the sample go through macs3 peak-calling.
- `config/config.yaml` — references (FASTA + `.fai`) and every model/ATAC
  parameter, passed to the scripts as explicit CLI flags. Defaults mirror
  `cherimoya_cli/defaults.py`.

## Environment

A single conda env (`workflow/envs/cherimoya.yaml`) pip-installs cherimoya from a
pinned git commit, which pulls torch (CUDA), tangermeme, macs3, and bam2bw.
Snakemake builds and caches it under `.conda/` when run with the profiles below.
Exact pins from a solved build are in `workflow/envs/cherimoya.{conda,pip}-lock.txt`
(recreate with `conda create --file cherimoya.conda-lock.txt` + `pip install -r
cherimoya.pip-lock.txt`).

## Running

```bash
# one-time: create your config from the tracked template, then edit
cp config/config.yaml.template config/config.yaml

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
