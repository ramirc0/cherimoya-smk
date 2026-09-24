# Load config + sample sheet, derive per-sample lookups and the rule helpers.

import json
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"        # per-sample outputs, grouped per run
# assembly -> {fasta, fai, chrom_sizes, gsize}. Each sample picks its genome via
# the sample sheet's `genome` column (multi-species data).
GENOMES = config["genomes"]
BLACKLIST = config["references"].get("blacklist") or None
FOLDS = [str(f) for f in config.get("folds", [0])]   # CV folds to train/evaluate
# Stranded runs emit a (+, -) bigWig pair per signal/control (one group of 2);
# unstranded runs emit a single track. The single gate for the whole workflow.
STRANDED = not config["preprocess"]["unstranded"]
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# Absolute so the conda: directive resolves the same from any rule file.
CONDA_ENV = str((Path(workflow.basedir).parent / config["conda_env"]).resolve())


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Create it (columns: sample_id, signal, genome, [control], [peaks]); "
        "signal is a BAM/fragments file (bam2bw + macs3) or an already-built "
        "bigWig (requires provided peaks). See config/samples.tsv for the layout."
    )

_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0)


def _column(name):
    """sample_id -> value for `name`, treating empty strings as absent (None)."""
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


def _is_bigwig(path):
    """True when `path` is already a bigWig (bam2bw can be skipped)."""
    return path.lower().endswith((".bw", ".bigwig"))


SAMPLES = _manifest.get_column("sample_id").to_list()
SIGNAL_OF = _column("signal")
CONTROL_OF = _column("control")
PEAKS_OF = _column("peaks")
GENOME_OF = _column("genome")

# QC covariates for perf_vs_covariate, gated by config["qc"] and computability:
# n_peaks is always derivable; n_fragments needs a scannable (non-bigWig) signal.
_qc = config["qc"]
_has_depth = any(not _is_bigwig(SIGNAL_OF[s]) for s in SAMPLES)
COVARIATES = [c for c, on in (("n_peaks", _qc["n_peaks"]),
                              ("n_fragments", _qc["n_fragments"] and _has_depth))
              if on]

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")

_no_genome = [s for s in SAMPLES if not GENOME_OF.get(s)]
if _no_genome:
    raise WorkflowError(
        f"{len(_no_genome)} sample(s) have no `genome` (e.g. {_no_genome[:3]}); "
        "every row needs a genome present in config['genomes']."
    )
_unknown_genomes = sorted({g for g in GENOME_OF.values() if g not in GENOMES})
if _unknown_genomes:
    raise WorkflowError(
        f"genome(s) not defined in config['genomes']: {_unknown_genomes}"
    )
# A bigWig signal skips bam2bw; peaks cannot be called from it, so they must be
# provided (macs3 needs a BAM/fragments file).
_bw_no_peaks = [s for s in SAMPLES
                if _is_bigwig(SIGNAL_OF[s]) and not PEAKS_OF.get(s)]
if _bw_no_peaks:
    raise WorkflowError(
        f"{len(_bw_no_peaks)} sample(s) give a bigWig signal but no peaks "
        f"(e.g. {_bw_no_peaks[:3]}); provide a `peaks` file or a BAM/fragments "
        "signal so macs3 can call peaks."
    )
# A stranded run needs bam2bw to build the (+, -) pair; a provided bigWig is a
# single track that cannot be split into strands.
if STRANDED:
    _bw_stranded = [s for s in SAMPLES if _is_bigwig(SIGNAL_OF[s])
                    or (CONTROL_OF.get(s) and _is_bigwig(CONTROL_OF[s]))]
    if _bw_stranded:
        raise WorkflowError(
            f"{len(_bw_stranded)} sample(s) give a pre-built bigWig under a "
            f"stranded run (e.g. {_bw_stranded[:3]}); a stranded (+, -) pair "
            "can only come from bam2bw. Provide a BAM signal/control or set "
            "preprocess.unstranded: true."
        )


wildcard_constraints:
    sample=r"[A-Za-z0-9][A-Za-z0-9_.-]*",
    fold=r"[0-4]",


def prefix(sample):
    """Canonical output prefix for a sample: results/<sample>/<sample>."""
    return f"{OUTDIR}/{sample}/{sample}"


def signal_bw(sample):
    """The signal bigWig(s) for `sample` as a fit/evaluate input list.

    A provided bigWig is used as-is (bam2bw skipped); otherwise the bam2bw
    output(s): a stranded run gives the (+, -) pair, an unstranded run a
    single track.
    """
    signal = SIGNAL_OF[sample]
    if _is_bigwig(signal):
        return [signal]
    p = prefix(sample)
    return [f"{p}.+.bw", f"{p}.-.bw"] if STRANDED else [f"{p}.bw"]


def control_bw(sample):
    """The control bigWig(s) for `sample` as a fit/evaluate model input, or []
    when the sample has no control.

    Mirrors `signal_bw` strandedness: a stranded run gives the control (+, -)
    pair, an unstranded run a single track. A provided bigWig is used as-is.
    """
    ctl = CONTROL_OF.get(sample)
    if not ctl:
        return []
    if _is_bigwig(ctl):
        return [ctl]
    p = f"{prefix(sample)}.control"
    return [f"{p}.+.bw", f"{p}.-.bw"] if STRANDED else [f"{p}.bw"]


def n_fragments_file(sample):
    """The count_fragments output (fragment count) for `sample`."""
    return f"{prefix(sample)}.n_fragments.txt"


# Per-sample genome lookups (the sample's assembly from the sheet).
def fasta_of(sample):
    """Genome FASTA for a sample (negatives/fit/evaluate)."""
    return GENOMES[GENOME_OF[sample]]["fasta"]


def gsize_of(sample):
    """macs3 effective genome size for a sample (preset like hs/mm, or a number)."""
    return GENOMES[GENOME_OF[sample]]["gsize"]


def chrom_sizes_of(sample):
    """Pre-generated chrom.sizes for the sample's genome (bam2bw -s)."""
    return GENOMES[GENOME_OF[sample]]["chrom_sizes"]


# CV fold chromosome splits (chrombpnet-style), per genome, from resources/folds.
def fold_json(sample, fold):
    """Fold-split JSON for a sample's genome (also a fit/evaluate input)."""
    return f"resources/folds/{GENOME_OF[sample]}/fold_{fold}.json"


def _fold(sample, fold):
    with open(fold_json(sample, fold)) as fh:
        return json.load(fh)


def fold_prefix(sample, fold):
    """Per-(sample, fold) output prefix: results/<sample>/fold_<fold>/<sample>."""
    return f"{OUTDIR}/{sample}/fold_{fold}/{sample}"


# Provided peaks, else macs3 output.
def peaks_for(wildcards):
    """Peak file for a sample: a provided (normalized) file, else macs3 output."""
    if PEAKS_OF.get(wildcards.sample):
        return f"{prefix(wildcards.sample)}.peaks.narrowPeak"
    return f"{prefix(wildcards.sample)}_peaks.narrowPeak"


def macs3_control_input(wildcards):
    """Control file(s) fed to macs3 as -c, or [] when the sample has none."""
    ctl = CONTROL_OF.get(wildcards.sample)
    return [ctl] if ctl else []


def macs3_format(wildcards):
    """The macs3 -f format string derived from the signal file + config."""
    pp = config["preprocess"]
    if pp["fragments"]:
        return "FRAG"
    fname = SIGNAL_OF[wildcards.sample]
    ext = fname.split(".")[-2] if fname.endswith(".gz") else fname.split(".")[-1]
    ext = ext.upper()
    if pp["paired_end"]:
        ext += "PE"
    return ext


def _list_flag(name, values):
    return [f"--{name}", *(str(v) for v in values)]


# The blacklist as a rule input (for dependency tracking) or [] when unset.
def blacklist_input(wildcards):
    return [BLACKLIST] if BLACKLIST else []


def fit_flags(sample, fold):
    """All fit.py hyperparameter flag tokens from config['fit'] (use with :q).

    Training/validation chroms come from the sample's genome CV `fold`.
    """
    f = config["fit"]
    flags = [
        "--n_filters", f["n_filters"],
        "--n_layers", f["n_layers"],
        "--expansion", f["expansion"],
        "--residual_scale", f["residual_scale"],
        "--in_window", f["in_window"],
        "--out_window", f["out_window"],
        "--max_jitter", f["max_jitter"],
        "--negative_ratio", f["negative_ratio"],
        "--batch_size", f["batch_size"],
        "--max_epochs", f["max_epochs"],
        "--n_warmup_epochs", f["n_warmup_epochs"],
        "--early_stopping", f["early_stopping"],
        "--muon_lr", f["muon_lr"],
        "--muon_wd", f["muon_wd"],
        "--adam_lr", f["adam_lr"],
        "--adam_wd", f["adam_wd"],
        "--lw_lr", f["lw_lr"],
        "--lw_wd", f["lw_wd"],
        "--lw_momentum", f["lw_momentum"],
        "--num_workers", f["num_workers"],
        "--dtype", f["dtype"],
        "--device", f["device"],
        *_list_flag("training_chroms", _fold(sample, fold)["train"]),
        *_list_flag("validation_chroms", _fold(sample, fold)["valid"]),
    ]
    if STRANDED:
        flags.append("--stranded")
    if not f["reverse_complement"]:
        flags.append("--no_reverse_complement")
    if f["summits"]:
        flags.append("--summits")
    if f["random_state"] is not None:
        flags += ["--random_state", f["random_state"]]
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]


def eval_flags(sample, fold):
    """All evaluate.py flag tokens from config['evaluate'] (use with :q).

    Eval chroms are the sample's genome CV `fold` test set.
    """
    e = config["evaluate"]
    flags = [
        "--batch_size", e["batch_size"],
        "--in_window", e["in_window"],
        "--out_window", e["out_window"],
        "--dtype", e["dtype"],
        "--device", e["device"],
        *_list_flag("chroms", _fold(sample, fold)["test"]),
    ]
    if STRANDED:
        flags.append("--stranded")
    if e["reverse_complement_average"]:
        flags.append("--reverse_complement_average")
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]
