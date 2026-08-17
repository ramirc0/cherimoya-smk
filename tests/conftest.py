"""Make the workflow scripts importable as top-level modules in tests."""

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
