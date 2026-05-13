#!/usr/bin/env python3
"""Run Rigel and evaluate a simulated benchmark suite."""

from __future__ import annotations


def main() -> None:
    from rigel.sim.analysis import main as run_main

    run_main()


if __name__ == "__main__":
    main()