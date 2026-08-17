"""Assert the script defaults still match cherimoya_cli/defaults.py."""

import pytest

import fit
import evaluate

defaults = pytest.importorskip("cherimoya_cli.defaults")


# Keys present in both our parser and the cherimoya_cli default dict.
FIT_KEYS = [
    "n_filters", "n_layers", "expansion", "residual_scale", "batch_size",
    "in_window", "out_window", "max_jitter", "negative_ratio",
    "reverse_complement", "summits", "max_epochs", "n_warmup_epochs",
    "early_stopping", "muon_lr", "muon_wd", "adam_lr", "adam_wd", "lw_lr",
    "lw_wd", "lw_momentum", "num_workers", "dtype", "device",
    "training_chroms", "validation_chroms",
]

EVAL_KEYS = [
    "batch_size", "in_window", "out_window", "reverse_complement_average",
    "device", "dtype", "chroms",
]


@pytest.mark.parametrize("key", FIT_KEYS)
def test_fit_defaults_match(key):
    parser = fit.build_parser()
    assert parser.get_default(key) == defaults.default_fit_parameters[key]


@pytest.mark.parametrize("key", EVAL_KEYS)
def test_evaluate_defaults_match(key):
    parser = evaluate.build_parser()
    assert parser.get_default(key) == defaults.default_evaluate_parameters[key]


def test_negatives_defaults():
    import negatives

    parser = negatives.build_parser()
    assert parser.get_default("bin_width") == 0.02
    assert parser.get_default("max_n_perc") == 0.1
    assert parser.get_default("beta") == 0.5
    assert parser.get_default("in_window") == 2114
    assert parser.get_default("out_window") == 1000
