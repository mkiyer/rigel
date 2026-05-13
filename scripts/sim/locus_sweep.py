#!/usr/bin/env python3
"""Run small locus-level simulation diagnostics."""

from __future__ import annotations

import sys


def main():
    from rigel.sim.locus_sweep import main as run_main

    return run_main()


if __name__ == "__main__":
    sys.exit(main())