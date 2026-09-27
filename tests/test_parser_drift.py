"""Assert the script defaults still match cherimoya_cli/defaults.py."""

import pytest

import fit
import evaluate
import attribute
import seqlets
import annotate

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
    "device", "dtype", "compile", "compile_mode", "chroms",
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


@pytest.mark.parametrize("key", FIT_KEYS)
def test_fit_defaults_match(key):
    parser = fit.build_parser()
    assert parser.get_default(key) == defaults.default_fit_parameters[key]


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


def test_negatives_defaults():
    import negatives

    parser = negatives.build_parser()
    assert parser.get_default("bin_width") == 0.02
    assert parser.get_default("max_n_perc") == 0.1
    assert parser.get_default("beta") == 0.5
    assert parser.get_default("in_window") == 2114
    assert parser.get_default("out_window") == 1000
