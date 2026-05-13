#!/usr/bin/env python3
"""Generate a complete synthetic benchmark suite and manifest."""

from __future__ import annotations

import sys


def main():
    from rigel.sim.suite import main as run_main

    return run_main()


if __name__ == "__main__":
    sys.exit(main())