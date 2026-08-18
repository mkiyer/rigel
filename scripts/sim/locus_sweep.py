#!/usr/bin/env python3
"""HOW DOES ONE LOCUS BEHAVE ACROSS A COMBINATORIAL SWEEP? — a CLI wrapper, one line of its own.

Thin entry point over :mod:`rigel.sim.locus_sweep`, which owns the sweep and its documentation:
transcripts given by strand and exon coordinates, RNA placed inside introns by genomic span, and one
abundance per group. ⚠ The engine's docstring is the long form — this file exists so the sweep has a
`scripts/sim/` command like every other simulator stage, and it must not accumulate logic of its own.
"""

from __future__ import annotations

import sys


def main():
    from rigel.sim.locus_sweep import main as run_main

    return run_main()


if __name__ == "__main__":
    sys.exit(main())