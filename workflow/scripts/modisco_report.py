#!/usr/bin/env python
"""Run `modisco report` with pandas' pre-3.0 string handling; takes its arguments."""

import runpy
import sys
from pathlib import Path

import pandas

# modisco 2.5.2 fills missing TomTom matches with None and tests them for
# truthiness. pandas 3 turns them into NaN, which crashes the report whenever
# a pattern has fewer than --n_matches matches.
pandas.set_option("future.infer_string", False)
sys.argv[0] = str(Path(sys.executable).parent / "modisco")
runpy.run_path(sys.argv[0], run_name="__main__")
