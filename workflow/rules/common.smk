# Load config + sample sheet, derive per-sample lookups and the rule helpers.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"        # per-sample outputs, grouped per run
FASTA = config["references"]["fasta"]
FAI = config["references"]["fai"]
BLACKLIST = config["references"].get("blacklist") or None
CHROM_SIZES = f"{RESULTS}/refs/chrom.sizes"   # run-independent, shared across runs
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# Absolute so the conda: directive resolves the same from any rule file.
CONDA_ENV = str((Path(workflow.basedir).parent / config["conda_env"]).resolve())


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Create it (columns: sample_id, signal, [control], [peaks]); see "
        "config/samples.tsv for the expected layout."
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


SAMPLES = _manifest.get_column("sample_id").to_list()
SIGNAL_OF = _column("signal")
CONTROL_OF = _column("control")
PEAKS_OF = _column("peaks")

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")


wildcard_constraints:
    sample=r"[A-Za-z0-9][A-Za-z0-9_.-]*",


def prefix(sample):
    """Canonical output prefix for a sample: results/<sample>/<sample>."""
    return f"{OUTDIR}/{sample}/{sample}"


def signal_bw(sample):
    """The unstranded bigWig produced by the bam2bw rule for `sample`."""
    return f"{prefix(sample)}.bw"


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
    return f"--{name} " + " ".join(str(v) for v in values)


# The blacklist as a rule input (for dependency tracking) or [] when unset.
def blacklist_input(wildcards):
    return [BLACKLIST] if BLACKLIST else []


def fit_flags():
    """All fit.py hyperparameter flags built from config['fit']."""
    f = config["fit"]
    flags = [
        f"--n_filters {f['n_filters']}",
        f"--n_layers {f['n_layers']}",
        f"--expansion {f['expansion']}",
        f"--residual_scale {f['residual_scale']}",
        f"--in_window {f['in_window']}",
        f"--out_window {f['out_window']}",
        f"--max_jitter {f['max_jitter']}",
        f"--negative_ratio {f['negative_ratio']}",
        f"--batch_size {f['batch_size']}",
        f"--max_epochs {f['max_epochs']}",
        f"--n_warmup_epochs {f['n_warmup_epochs']}",
        f"--early_stopping {f['early_stopping']}",
        f"--muon_lr {f['muon_lr']}",
        f"--muon_wd {f['muon_wd']}",
        f"--adam_lr {f['adam_lr']}",
        f"--adam_wd {f['adam_wd']}",
        f"--lw_lr {f['lw_lr']}",
        f"--lw_wd {f['lw_wd']}",
        f"--lw_momentum {f['lw_momentum']}",
        f"--num_workers {f['num_workers']}",
        f"--dtype {f['dtype']}",
        f"--device {f['device']}",
        _list_flag("training_chroms", f["training_chroms"]),
        _list_flag("validation_chroms", f["validation_chroms"]),
    ]
    if not f["reverse_complement"]:
        flags.append("--no_reverse_complement")
    if f["summits"]:
        flags.append("--summits")
    if f.get("random_state") is not None:
        flags.append(f"--random_state {f['random_state']}")
    if BLACKLIST:
        flags.append(f"--exclusion_lists {BLACKLIST}")
    return " ".join(flags)


def eval_flags():
    """All evaluate.py flags built from config['evaluate']."""
    e = config["evaluate"]
    flags = [
        f"--batch_size {e['batch_size']}",
        f"--in_window {e['in_window']}",
        f"--out_window {e['out_window']}",
        f"--dtype {e['dtype']}",
        f"--device {e['device']}",
        _list_flag("chroms", e["chroms"]),
    ]
    if e["reverse_complement_average"]:
        flags.append("--reverse_complement_average")
    if BLACKLIST:
        flags.append(f"--exclusion_lists {BLACKLIST}")
    return " ".join(flags)
