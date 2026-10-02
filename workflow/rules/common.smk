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
FOLDS = [str(f) for f in config["folds"]]   # CV folds; every per-fold rule fans out over them
if "folds" in config["attribute"]:
    raise WorkflowError("attribute.folds was removed; interpretation runs on every fold in `folds`.")
# Stranded runs emit a (+, -) bigWig pair per signal/control (one group of 2);
# unstranded runs emit a single track. The single gate for the whole workflow.
STRANDED = not config["preprocess"]["unstranded"]
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# The workspace path is relative to the rule file; every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


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
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


def _is_bigwig(path):
    """Check whether `path` is already a bigWig, so bam2bw can be skipped.

    Parameters
    ----------
    path : str
        Signal or control file.

    Returns
    -------
    bool
        True for a `.bw` or `.bigwig` extension, in any case.
    """
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
    """Return the output prefix of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<sample>/<sample>`.
    """
    return f"{OUTDIR}/{sample}/{sample}"


def signal_bw(sample):
    """Return the signal bigWig(s) of a sample, as fit and evaluate take them.

    A provided bigWig is used as-is and bam2bw is skipped. Otherwise these
    are the bam2bw outputs: the (+, -) pair for a stranded run, one track for
    an unstranded run.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    list of str
        One or two bigWig paths.
    """
    signal = SIGNAL_OF[sample]
    if _is_bigwig(signal):
        return [signal]
    p = prefix(sample)
    return [f"{p}.+.bw", f"{p}.-.bw"] if STRANDED else [f"{p}.bw"]


def control_bw(sample):
    """Return the control bigWig(s) of a sample, as fit and evaluate take them.

    Strandedness follows `signal_bw`: the (+, -) pair for a stranded run, one
    track for an unstranded run. A provided bigWig is used as-is.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    list of str
        Up to two bigWig paths. Empty when the sample has no control.
    """
    ctl = CONTROL_OF.get(sample)
    if not ctl:
        return []
    if _is_bigwig(ctl):
        return [ctl]
    p = f"{prefix(sample)}.control"
    return [f"{p}.+.bw", f"{p}.-.bw"] if STRANDED else [f"{p}.bw"]


def n_fragments_file(sample):
    """Return the count_fragments output of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    list of str
        The fragment-count file. Empty for a bigWig signal or when
        `qc.n_fragments` is off.
    """
    if not config["qc"]["n_fragments"] or _is_bigwig(SIGNAL_OF[sample]):
        return []
    return [f"{prefix(sample)}.n_fragments.txt"]


# Per-sample genome lookups (the sample's assembly from the sheet).
def fasta_of(sample):
    """Return the genome FASTA of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        FASTA of the sample's genome.
    """
    return GENOMES[GENOME_OF[sample]]["fasta"]


def gsize_of(sample):
    """Return the macs3 effective genome size of a sample.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str or int
        A macs3 preset such as `hs` or `mm`, or a size in bp.
    """
    return GENOMES[GENOME_OF[sample]]["gsize"]


def chrom_sizes_of(sample):
    """Return the chrom.sizes of a sample's genome, for bam2bw `-s`.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.

    Returns
    -------
    str
        Pre-generated chrom.sizes file.
    """
    return GENOMES[GENOME_OF[sample]]["chrom_sizes"]


# CV fold chromosome splits (chrombpnet-style), per genome, from resources/folds.
def fold_json(sample, fold):
    """Return the CV fold JSON of a sample's genome.

    It is also an input of the per-fold rules.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    str
        Path `resources/folds/<genome>/fold_<fold>.json`.
    """
    return f"resources/folds/{GENOME_OF[sample]}/fold_{fold}.json"


def _fold(sample, fold):
    """Load the CV fold JSON of a sample's genome.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    dict of str to list of str
        Chromosomes under `train`, `valid` and `test`.
    """
    with open(fold_json(sample, fold)) as fh:
        return json.load(fh)


def fold_prefix(sample, fold):
    """Return the output prefix of a sample's model for one CV fold.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<sample>/fold_<fold>/<sample>`.
    """
    return f"{OUTDIR}/{sample}/fold_{fold}/{sample}"


def per_fold(*suffixes):
    """Return per-model output paths for every sample and fold.

    Parameters
    ----------
    *suffixes : str
        File suffixes, each appended to `fold_prefix` after a dot.

    Returns
    -------
    list of str
        Paths `<fold_prefix>.<suffix>` for every sample, fold and suffix.
    """
    return [f"{fold_prefix(s, fold)}.{x}"
            for s in SAMPLES for fold in FOLDS for x in suffixes]


def peaks_for(wildcards):
    """Return the peak file of a sample.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Rule wildcards with `sample`.

    Returns
    -------
    str
        The prep_peaks copy of the provided peaks, or the macs3 output when
        the sheet gives none.
    """
    if PEAKS_OF.get(wildcards.sample):
        return f"{prefix(wildcards.sample)}.peaks.narrowPeak"
    return f"{prefix(wildcards.sample)}_peaks.narrowPeak"


def macs3_control_input(wildcards):
    """Return the control file(s) macs3 takes with `-c`.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Rule wildcards with `sample`.

    Returns
    -------
    list of str
        The sample's control. Empty when it has none.
    """
    ctl = CONTROL_OF.get(wildcards.sample)
    return [ctl] if ctl else []


def macs3_format(wildcards):
    """Return the macs3 `-f` format of a sample's signal.

    `FRAG` when `preprocess.fragments` is set. Otherwise the signal's file
    extension without `.gz`, upper-cased, with `PE` appended for paired-end
    data.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Rule wildcards with `sample`.

    Returns
    -------
    str
        A macs3 format such as `BAM`, `BAMPE` or `FRAG`.
    """
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
    """Return a multi-value flag as tokens.

    Parameters
    ----------
    name : str
        Flag name without the leading dashes.
    values : iterable
        Flag values.

    Returns
    -------
    list of str
        The flag `--<name>` followed by each value.
    """
    return [f"--{name}", *(str(v) for v in values)]


def blacklist_input(wildcards):
    """Return the blacklist as a rule input, so rules track it.

    Parameters
    ----------
    wildcards : snakemake.io.Wildcards
        Rule wildcards. Unused.

    Returns
    -------
    list of str
        The blacklist BED. Empty when `references.blacklist` is unset.
    """
    return [BLACKLIST] if BLACKLIST else []


def fit_flags(sample, fold):
    """Return the fit.py flag tokens from `config["fit"]`.

    Training and validation chroms come from the sample's genome CV fold.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
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
        "--min_total_steps", f["min_total_steps"],
        "--n_warmup_epochs", f["n_warmup_epochs"],
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
        "--compile_mode", f["compile_mode"],
        "--random_state", config["random_state"],
        *_list_flag("training_chroms", _fold(sample, fold)["train"]),
        *_list_flag("validation_chroms", _fold(sample, fold)["valid"]),
    ]
    if STRANDED:
        flags.append("--stranded")
    if not f["reverse_complement"]:
        flags.append("--no_reverse_complement")
    if f["summits"]:
        flags.append("--summits")
    if not f["compile"]:
        flags.append("--no_compile")
    if f["early_stopping"] is not None:
        flags += ["--early_stopping", f["early_stopping"]]
    if f["loss_weights"] is not None:
        flags += _list_flag("loss_weights", f["loss_weights"])
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]


def eval_flags(sample, fold):
    """Return the evaluate.py flag tokens from `config["evaluate"]`.

    Eval chroms are the test set of the sample's genome CV fold. Windows and
    summits are fit's, as the CLI's evaluate inherits them.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    e = config["evaluate"]
    flags = [
        "--batch_size", e["batch_size"],
        "--in_window", config["fit"]["in_window"],
        "--out_window", config["fit"]["out_window"],
        "--dtype", e["dtype"],
        "--device", e["device"],
        "--compile_mode", e["compile_mode"],
        *_list_flag("chroms", _fold(sample, fold)["test"]),
    ]
    if STRANDED:
        flags.append("--stranded")
    if e["reverse_complement_average"]:
        flags.append("--reverse_complement_average")
    if config["fit"]["summits"]:
        flags.append("--summits")
    if not e["compile"]:
        flags.append("--no_compile")
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]


def attr_flags(sample, fold):
    """Return the attribute.py flag tokens from `config["attribute"]`.

    Chroms are the train and valid sets of the sample's genome CV fold.
    `in_window` is fit's. `cherimoya pipeline` shares both.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    a = config["attribute"]
    flags = [
        "--algorithm", a["algorithm"],
        "--output", a["output"],
        "--group", a["group"],
        "--attr_window", a["attr_window"],
        "--n_shuffles", a["n_shuffles"],
        "--warning_threshold", a["warning_threshold"],
        "--random_state", config["random_state"],
        "--batch_size", a["batch_size"],
        "--in_window", config["fit"]["in_window"],
        "--dtype", a["dtype"],
        "--device", a["device"],
        "--compile_mode", a["compile_mode"],
        *_list_flag("chroms", _fold(sample, fold)["train"] + _fold(sample, fold)["valid"]),
    ]
    if a["compile"]:
        flags.append("--compile")
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]


def seqlet_flags(sample, fold):
    """Return the seqlets.py flag tokens from `config["seqlets"]`.

    Chroms are attribute's, so loci line up with its index mask.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    s = config["seqlets"]
    flags = [
        "--threshold", s["threshold"],
        "--min_seqlet_len", s["min_seqlet_len"],
        "--max_seqlet_len", s["max_seqlet_len"],
        "--additional_flanks", s["additional_flanks"],
        "--n_bins", s["n_bins"],
        *_list_flag("chroms", _fold(sample, fold)["train"] + _fold(sample, fold)["valid"]),
    ]
    return [str(x) for x in flags]


def annot_flags():
    """Return the annotate.py flag tokens from `config["annotate"]`.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    a = config["annotate"]
    flags = [
        "--n_score_bins", a["n_score_bins"],
        "--n_median_bins", a["n_median_bins"],
        "--n_target_bins", a["n_target_bins"],
        "--n_cache", a["n_cache"],
        "--n_nearest", a["n_nearest"],
    ]
    if not a["reverse_complement"]:
        flags.append("--no_reverse_complement")
    return [str(x) for x in flags]


def modisco_flags():
    """Return the modisco_motifs.py flag tokens from `config["modisco"]`.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    m = config["modisco"]
    flags = [
        "--n_seqlets", m["n_seqlets"],
        "--n_leiden", m["n_leiden"],
        "--window", m["window"],
        "--size", m["size"],
        "--trim_size", m["trim_size"],
        "--seqlet_flank_size", m["seqlet_flank_size"],
        "--initial_flank_to_add", m["initial_flank_to_add"],
        "--final_flank_to_add", m["final_flank_to_add"],
        "--target_seqlet_fdr", m["target_seqlet_fdr"],
        "--min_metacluster_size", m["min_metacluster_size"],
        "--n_leiden_iterations", m["n_leiden_iterations"],
        "--final_min_cluster_size", m["final_min_cluster_size"],
    ]
    return [str(x) for x in flags]


def modisco_report_flags():
    """Return the `modisco report` flag tokens from `config["modisco_report"]`.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    r = config["modisco_report"]
    flags = [
        "--n_matches", r["n_matches"],
        "--n_examples", r["n_examples"],
        "--trim_threshold", r["trim_threshold"],
    ]
    if r["lite"]:
        flags.append("--lite")
    return [str(x) for x in flags]


def marginalize_flags(sample, fold):
    """Return the marginalize.py flag tokens from `config["marginalize"]`.

    Chroms are the train set of the sample's genome CV fold, because
    `cherimoya pipeline` leaves marginalize at its `training_chroms` default.
    `in_window` is fit's.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    fold : str
        CV fold index.

    Returns
    -------
    list of str
        Flag tokens. `shell:` quotes each with `:q`.
    """
    m = config["marginalize"]
    flags = [
        "--n_loci", m["n_loci"],
        "--random_state", config["random_state"],
        "--batch_size", m["batch_size"],
        "--in_window", config["fit"]["in_window"],
        "--device", m["device"],
        "--compile_mode", m["compile_mode"],
        *_list_flag("chroms", _fold(sample, fold)["train"]),
    ]
    if m["shuffle"]:
        flags.append("--shuffle")
    if m["attributions"]:
        flags.append("--attributions")
    if not m["minimal"]:
        flags.append("--no_minimal")
    if not m["compile"]:
        flags.append("--no_compile")
    if BLACKLIST:
        flags += ["--exclusion_lists", BLACKLIST]
    return [str(x) for x in flags]
