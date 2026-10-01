"""Assert the script defaults still match cherimoya_cli/defaults.py."""

import inspect
import re
import subprocess
import sys
from pathlib import Path

import pytest

import fit
import evaluate
import attribute
import seqlets
import annotate
import modisco_motifs
import marginalize

defaults = pytest.importorskip("cherimoya_cli.defaults")


# Keys present in both our parser and the cherimoya_cli default dict.
FIT_KEYS = [
    "n_filters", "n_layers", "expansion", "residual_scale", "batch_size",
    "in_window", "out_window", "max_jitter", "negative_ratio",
    "reverse_complement", "summits", "max_epochs", "min_total_steps",
    "loss_weights", "n_warmup_epochs", "early_stopping", "muon_lr",
    "muon_wd", "adam_lr", "adam_wd", "lw_lr", "lw_wd", "lw_momentum",
    "num_workers", "dtype", "device",
    "random_state", "training_chroms", "validation_chroms",
]

EVAL_KEYS = [
    "batch_size", "in_window", "out_window", "reverse_complement_average",
    "summits", "device", "dtype", "compile", "compile_mode", "chroms",
]

ATTR_KEYS = [
    "batch_size", "in_window", "chroms", "compile", "compile_mode",
    "algorithm", "output", "group", "attr_window", "n_shuffles",
    "warning_threshold", "random_state", "ohe_filename", "attr_filename",
    "idx_filename", "dtype", "device", "verbose",
]

SEQLET_KEYS = [
    "threshold", "min_seqlet_len", "max_seqlet_len", "additional_flanks",
    "chroms", "output_filename",
]

ANNOT_KEYS = [
    "n_score_bins", "n_median_bins", "n_target_bins", "n_cache",
    "reverse_complement", "n_jobs", "output_filename",
]

MODISCO_KEYS = ["n_seqlets"]

# `modisco motifs` flags the pipeline leaves at the CLI default.
MODISCO_CLI_KEYS = [
    "n_leiden", "window", "size", "trim_size", "seqlet_flank_size",
    "initial_flank_to_add", "final_flank_to_add",
]

# Extras passed to TFMoDISco; the CLI hard-codes target_seqlet_fdr=0.05.
MODISCO_LIB_KEYS = ["min_metacluster_size", "n_leiden_iterations", "final_min_cluster_size"]

MARGINALIZE_KEYS = [
    "batch_size", "in_window", "chroms", "n_loci", "shuffle", "random_state",
    "attributions", "minimal", "compile", "compile_mode", "output_filename",
    "device", "verbose",
]

TEMPLATES = sorted((Path(__file__).resolve().parents[1] / "config").glob("*.yaml.template"))


def _template_block(template, name):
    """Map each key of a flat top-level block in a config template to its value string."""
    block = re.search(rf"^{name}:.*?\n((?:  .*\n)+)", template.read_text(), re.M).group(1)
    return dict(re.findall(r"^  (\w+): *([^\s#]+)", block, re.M))


def _cli_defaults(*cmd):
    """Map each `--flag` of a `modisco` subcommand to its help default string."""
    exe = Path(sys.executable).parent / "modisco"
    text = subprocess.run([str(exe), *cmd, "--help"], capture_output=True, text=True,
        check=True).stdout
    blocks = [" ".join(b.split()) for b in re.split(r"\n(?=\s{2}-)", text)]
    return {m.group(1): re.findall(r"\(default: ([^)]*)\)", b)[-1]
            for b in blocks if (m := re.search(r"--(\w+)", b))
            and "(default:" in b}


@pytest.mark.parametrize("key", FIT_KEYS)
def test_fit_defaults_match(key):
    parser = fit.build_parser()
    assert parser.get_default(key) == defaults.default_fit_parameters[key]


def test_fit_library_defaults_match():
    from cherimoya.training import fit as library_fit

    params = inspect.signature(library_fit).parameters
    for key in ("devices", "progress_bar"):
        assert params[key].default == defaults.default_fit_parameters[key], key


@pytest.mark.parametrize("key", EVAL_KEYS)
def test_evaluate_defaults_match(key):
    parser = evaluate.build_parser()
    assert parser.get_default(key) == defaults.default_evaluate_parameters[key]


@pytest.mark.parametrize("key", ATTR_KEYS)
def test_attribute_defaults_match(key):
    parser = attribute.build_parser()
    assert parser.get_default(key) == defaults.default_attribute_parameters[key]


@pytest.mark.parametrize("key", SEQLET_KEYS)
def test_seqlets_defaults_match(key):
    parser = seqlets.build_parser()
    assert parser.get_default(key) == defaults.default_seqlet_parameters[key]


@pytest.mark.parametrize("key", ANNOT_KEYS)
def test_annotate_defaults_match(key):
    parser = annotate.build_parser()
    assert parser.get_default(key) == defaults.default_annotation_parameters[key]


@pytest.mark.parametrize("key", MODISCO_KEYS)
def test_modisco_defaults_match(key):
    parser = modisco_motifs.build_parser()
    expected = defaults.default_pipeline_parameters["modisco_motifs_parameters"][key]
    assert parser.get_default(key) == expected


@pytest.mark.parametrize("key", MODISCO_CLI_KEYS)
def test_modisco_cli_defaults_match(key):
    parser = modisco_motifs.build_parser()
    assert str(parser.get_default(key)) == _cli_defaults("motifs")[key]


def test_modisco_library_defaults_match():
    from modiscolite.tfmodisco import TFMoDISco

    params = inspect.signature(TFMoDISco).parameters
    parser = modisco_motifs.build_parser()
    for key in MODISCO_LIB_KEYS:
        assert parser.get_default(key) == params[key].default, key
    assert parser.get_default("target_seqlet_fdr") == 0.05


@pytest.mark.parametrize("template", TEMPLATES, ids=[t.name for t in TEMPLATES])
def test_modisco_report_config_matches_cli(template):
    report = _template_block(template, "modisco_report")
    cli = _cli_defaults("report")
    for key in ("n_matches", "n_examples", "trim_threshold", "lite"):
        assert report[key].lower() == cli[key].lower(), key


@pytest.mark.parametrize("template", TEMPLATES, ids=[t.name for t in TEMPLATES])
def test_modisco_config_matches_script(template):
    config = _template_block(template, "modisco")
    parser = modisco_motifs.build_parser()
    assert len(config) == 12
    for key, value in config.items():
        assert str(parser.get_default(key)) == value, key


@pytest.mark.parametrize("key", MARGINALIZE_KEYS)
def test_marginalize_defaults_match(key):
    parser = marginalize.build_parser()
    assert parser.get_default(key) == defaults.default_marginalize_parameters[key]


@pytest.mark.parametrize("template", TEMPLATES, ids=[t.name for t in TEMPLATES])
def test_marginalize_config_matches_script(template):
    config = _template_block(template, "marginalize")
    parser = marginalize.build_parser()
    assert len(config) == 9
    for key, value in config.items():
        assert str(parser.get_default(key)).lower() == value, key


def test_negatives_defaults():
    import negatives

    parser = negatives.build_parser()
    assert parser.get_default("bin_width") == 0.02
    assert parser.get_default("max_n_perc") == 0.1
    assert parser.get_default("beta") == 0.5
    assert parser.get_default("in_window") == 2114
    assert parser.get_default("out_window") == 1000
