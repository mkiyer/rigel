#!/usr/bin/env python3
"""Generate whole-genome simulated reads from an existing FASTA/GTF config."""

from __future__ import annotations

import sys


def main() -> int:
    from rigel.sim.whole_genome import main as run_main

    return run_main()


if __name__ == "__main__":
    sys.exit(main())