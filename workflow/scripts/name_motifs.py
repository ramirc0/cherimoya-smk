#!/usr/bin/env python
"""Rewrite MEME motif IDs as NAME_ACCESSION so reports show readable names.

`modisco report` shows only the MEME motif ID, which in JASPAR is the
accession. `MOTIF MA0074.1 RXRA::VDR` becomes `MOTIF RXRA+VDR_MA0074.1`. The
accession keeps IDs unique. `::` becomes `+` and `/` becomes `-`, since the ID
is also a logo filename and a relative link in the report HTML. `MOTIF`
lines without a name are kept as-is.
"""

import argparse


def build_parser():
    """Build the name_motifs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the name_motifs.py flags.
    """
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="MEME file with `MOTIF <accession> <name>` lines.")
    parser.add_argument("output", help="Renamed MEME file.")
    return parser


def main():
    """Rewrite the motif IDs of one MEME file."""
    args = build_parser().parse_args()

    with open(args.input) as src, open(args.output, "w") as out:
        for line in src:
            fields = line.split(maxsplit=2)
            if fields[:1] == ["MOTIF"] and len(fields) == 3:
                name = fields[2].strip().replace("::", "+").replace("/", "-")
                line = f"MOTIF {name}_{fields[1]}\n"
            out.write(line)


if __name__ == "__main__":
    main()
