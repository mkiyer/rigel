#!/usr/bin/env python3
"""BUILD A WHOLE SYNTHETIC MINI-GENOME SUITE FROM NOTHING — a CLI wrapper, one line of its own.

Thin entry point over :mod:`rigel.sim.suite`, which owns the work and its documentation: a synthetic
genome, reads simulated across a strand × gDNA × intronic-RNA grid, and an oracle BAM plus paired-end
FASTQ per condition. ⚠ It carves its OWN genome, so it is not the gDNA ladder's builder — that is
`scripts/sim/panel.py build`, which works from a real human backbone (`docs/TESTING.md` §1).
"""

from __future__ import annotations

import sys


def main():
    from rigel.sim.suite import main as run_main

    return run_main()


if __name__ == "__main__":
    sys.exit(main())