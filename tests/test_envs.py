"""The V100 env must lock the same packages as the main env, CUDA stack aside."""

import re
from pathlib import Path

import yaml

LOCK = Path(__file__).resolve().parents[1] / "workflow" / "envs" / "pixi.lock"
CUDA = re.compile(r"(torch|triton|nvidia_\w+|cuda_\w+)-")


def _packages(env):
    entries = yaml.safe_load(LOCK.read_text())["environments"][env]["packages"][
        "linux-64"
    ]
    return {
        url
        for entry in entries
        for url in entry.values()
        if not CUDA.match(url.rsplit("/", 1)[-1])
    }


def test_sm70_env_matches_main_env_outside_cuda():
    assert _packages("cherimoya-sm70") == _packages("cherimoya")
