"""Guards for the stranded (+, -) grouping wired for cherimoya pipeline parity."""

import pytest

import fit
import evaluate

FIT_MIN = ["-s", "g.fa", "-l", "p.np", "-neg", "n.bed", "-sig", "a.bw", "-o", "o"]
EVAL_MIN = ["-s", "g.fa", "-l", "p.np", "-sig", "a.bw", "-m", "m.torch", "-o", "o"]


def test_fit_stranded_default_false():
    assert fit.build_parser().parse_args(FIT_MIN).stranded is False


def test_evaluate_stranded_default_false():
    assert evaluate.build_parser().parse_args(EVAL_MIN).stranded is False


@pytest.mark.parametrize("parser,argv", [(fit, FIT_MIN), (evaluate, EVAL_MIN)])
def test_stranded_flag_parses(parser, argv):
    assert parser.build_parser().parse_args(argv + ["--stranded"]).stranded is True


def test_stranded_grouping_contract():
    # The library grouping our --stranded relies on: a wrapped pair is one
    # (+, -) group; a flat list is independent unstranded groups.
    nsg = pytest.importorskip("cherimoya.io").normalize_signal_groups
    assert nsg([["a.+.bw", "a.-.bw"]])[1] == [2]
    assert nsg(["a.+.bw", "a.-.bw"])[1] == [1, 1]
