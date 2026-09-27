"""Static checks on the rules' shell blocks."""

import re
from pathlib import Path

import pytest

RULES = Path(__file__).resolve().parents[1] / "workflow" / "rules"


def _shell_blocks():
    for smk in sorted(RULES.glob("*.smk")):
        for name, body in re.findall(r"^rule (\w+):(.*?)(?=^rule |\Z)",
                                     smk.read_text(), re.S | re.M):
            block = re.search(r'shell:\s*r"""(.*?)"""', body, re.S)
            if block:
                yield name, block.group(1)


BLOCKS = list(_shell_blocks())


@pytest.mark.parametrize("name,shell", BLOCKS, ids=[n for n, _ in BLOCKS])
def test_command_does_not_end_on_a_params_only_line(name, shell):
    # Snakemake strips trailing whitespace from the command. If the last line
    # is a params list that renders empty, the command ends in a bare `\`,
    # which bash passes to the tool as a literal argument.
    last = [line.strip() for line in shell.splitlines() if line.strip()][-1]
    assert not re.fullmatch(r"(\{params\.\w+:q\}\s*)+", last), (
        f"rule {name} ends on {last!r}; move optional flags before required ones")
