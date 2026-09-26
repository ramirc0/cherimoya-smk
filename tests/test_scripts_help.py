"""Smoke-test each script's argparse CLI without heavy imports."""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"


@pytest.mark.parametrize("script", ["negatives.py", "fit.py", "evaluate.py", "model_summary.py"])
def test_help_exits_zero(script):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / script), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.lower()


@pytest.mark.parametrize("script", ["negatives.py", "fit.py", "evaluate.py", "model_summary.py"])
def test_missing_required_args_exits_two(script):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / script)],
        capture_output=True, text=True,
    )
    assert result.returncode == 2  # argparse's exit code for missing required args
