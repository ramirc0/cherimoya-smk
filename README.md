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
activates its env from `workflow/envs/.pixi/envs/`. An `onstart` hook installs
the selected env and `modisco-report` before the first job. Parallel SLURM jobs
would otherwise race to install them.

## Inputs

**Sample sheet** (TSV, path set by `samples:` in config). One row per sample:

| column      | required | meaning |
|-------------|----------|---------|
| `sample_id` | yes      | output name |
| `signal`    | yes      | BAM, fragments file, or a pre-built bigWig (`.bw`/`.bigwig`). A stranded run takes pre-built tracks as a (+, -) pair, `plus.bw,minus.bw` |
| `genome`    | yes      | a key under `genomes:` in config |
| `control`   | no       | control BAM, or pre-built bigWigs as for `signal`. Fed to the model. A BAM also goes to `macs3 -c` |
| `peaks`     | no       | narrowPeak. If empty, macs3 calls peaks from the BAM |

A bigWig signal or control skips `bam2bw` but can't be used for peak calling,
so its row needs a `peaks` file. See `docs/run-paths.svg` for every entry path.

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
  written by `make_folds.py`. hg38 uses the published chrombpnet folds and mm10
  the lab's folds. Other genomes get a length-balanced split. A genome with
  fewer than five main chromosomes, such as S. pombe, gets no folds.

To add a genome, add a `genomes:` entry, provide the three ref files, rerun
`make_folds.py`, and use the name in the sheet. An entry may also set
`blacklist`, a BED of regions that fit, evaluate, attribute and marginalize
exclude. Each genome carries its own, so one run can mix genomes.

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
  and a heatmap of the top motifs by sample. The distribution, covariate and
  outlier plots need `qc.min_models` models (samples x folds, default 10)
- `config.snapshot.json`: the fully resolved config, including `--config`
  overrides. Rerun it with `--configfile results/<run_id>/config.snapshot.json`.

Plots are written as SVG and PNG.

Keep runs apart with `--config run_id=mytag`. Swap sample sheets with
`--config samples=config/samples.other.tsv`.

## Run part of the pipeline

Name a stage to stop there. Each stage also builds every stage above it.

| target | adds |
|---|---|
| `preprocess` | bigWigs, peaks, negatives and fragment counts per sample; `config.snapshot.json` |
| `train` | model, evaluation, plots and `summary.txt` per fold; `metrics.tsv` and the run-level performance plots |
| `attributions` | attributions and the attribution profile plot per fold |
| `motifs` | seqlets, TomTom annotation and their plots per fold; the motif heatmap |
| `modisco` | TF-MoDISco patterns and report per fold |
| `all` | marginalization per fold. This is the default |

Put the target before `--config`:

```bash
snakemake train --profile profiles/slurm --config run_id=mytag
```

Later, run `all` with the same `run_id`. Only the missing stages run.

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
a pinned git commit. torch (CUDA), tangermeme, macs3, figwig and pytest are
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

`pixi run lint` (root workspace) runs `ruff check`, `ruff format --check` and
`numpydoc lint` on the scripts. `pyproject.toml` holds their config. Neither
tool parses `.smk` files, so keep the rule helpers numpydoc by hand.

## References

Amemiya, Haley M., Anshul Kundaje, and Alan P. Boyle. 2019. “The ENCODE Blacklist: Identification of Problematic Regions of the Genome.” Scientific Reports 9 (1): 9354. https://doi.org/10.1038/s41598-019-45839-z.

Avsec, Žiga, Melanie Weilert, Avanti Shrikumar, et al. 2021. “Base-Resolution Models of Transcription-Factor Binding Reveal Soft Motif Syntax.” Nature Genetics 53 (3): 354–66. https://doi.org/10.1038/s41588-021-00782-6.

Gupta, Shobhit, John A. Stamatoyannopoulos, Timothy L. Bailey, and William Stafford Noble. 2007. “Quantifying Similarity between Motifs.” Genome Biology 8 (2): R24. https://doi.org/10.1186/gb-2007-8-2-r24.

Lundberg, Scott M., and Su-In Lee. 2017. “A Unified Approach to Interpreting Model Predictions.” In Advances in Neural Information Processing Systems 30. https://proceedings.neurips.cc/paper/2017/hash/8a20a8621978632d76c43dfd28b67767-Abstract.html.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (January): 33. https://doi.org/10.12688/f1000research.29032.3.

Ovek Baydar, Damla, Ieva Rauluseviciute, Dina R. Aronsen, et al. 2025. “JASPAR 2026: Expansion of Transcription Factor Binding Profiles and Integration of Deep Learning Models.” Nucleic Acids Research 54 (D1): D184–93. https://doi.org/10.1093/nar/gkaf1209.

Pampari, Anusri, Anna Shcherbina, Evgeny Z. Kvon, et al. 2024. “ChromBPNet: Bias Factorized, Base-Resolution Deep Learning Models of Chromatin Accessibility Reveal Cis-Regulatory Sequence Syntax, Transcription Factor Footprints and Regulatory Variants.” Preprint, bioRxiv. https://doi.org/10.1101/2024.12.25.630221.

Schreiber, Jacob. 2025. “Tomtom-Lite: Accelerating Tomtom Enables Large-Scale and Real-Time Motif Similarity Scoring.” Bioinformatics 41 (11): btaf577. https://doi.org/10.1093/bioinformatics/btaf577.

Schreiber, Jacob. 2026. “Tangermeme: A Toolkit for Understanding Cis-Regulatory Logic Using Deep Learning Models.” Nature Methods, ahead of print, October 8. https://doi.org/10.1038/s41592-026-03254-z.

Schreiber, Jacob. n.d. Bpnet-Lite. Github. Accessed October 14, 2026. https://github.com/jmschrei/bpnet-lite.

Schreiber, Jacob. n.d. Cherimoya. Github. Accessed October 14, 2026. https://github.com/jmschrei/cherimoya.

Schreiber, Jacob. n.d. Figwig. Github. Accessed October 14, 2026. https://github.com/jmschrei/figwig.

Shrikumar, Avanti, Peyton Greenside, and Anshul Kundaje. 2017. “Learning Important Features through Propagating Activation Differences.” In Proceedings of the 34th International Conference on Machine Learning, 70:3145–53. PMLR. https://proceedings.mlr.press/v70/shrikumar17a.html.

Shrikumar, Avanti, Katherine Tian, Žiga Avsec, et al. 2018. “Technical Note on Transcription Factor Motif Discovery from Importance Scores (TF-MoDISco) Version 0.5.6.5.” Preprint, arXiv. https://doi.org/10.48550/arXiv.1811.00416.

Zhang, Yong, Tao Liu, Clifford A. Meyer, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9 (9): R137. https://doi.org/10.1186/gb-2008-9-9-r137.
