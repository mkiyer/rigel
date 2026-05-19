#!/usr/bin/env python3
"""Run Rigel regional-exposure quantification on the synthetic-24 suite."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_OUT_NAME = "rigel_regional_v3_out"
DEFAULT_ANNOTATED_NAME = "annotated_regional_v3.bam"


def load_conditions(base: Path, selected: list[str] | None) -> list[str]:
    if selected:
        return selected
    manifest = json.loads((base / "manifest.json").read_text())
    return [str(row["name"]) for row in manifest.get("conditions", [])]


def run_condition(
    base: Path,
    index_dir: Path,
    condition: str,
    out_name: str,
    annotated_name: str,
    seed: int,
    threads: int,
    force: bool,
) -> dict[str, object]:
    cond_dir = base / condition
    out_dir = cond_dir / out_name
    annotated_bam = cond_dir / annotated_name
    log_path = cond_dir / f"{out_name}.log"

    if (out_dir / "quant.feather").exists() and not force:
        return {
            "condition": condition,
            "status": "skipped",
            "returncode": 0,
            "out_dir": str(out_dir),
            "annotated_bam": str(annotated_bam),
            "log": str(log_path),
            "elapsed_sec": 0.0,
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "rigel",
        "quant",
        "--bam",
        str(cond_dir / "sim_oracle.bam"),
        "--index",
        str(index_dir),
        "-o",
        str(out_dir),
        "--annotated-bam",
        str(annotated_bam),
        "--sj-strand-tag",
        "auto",
        "--regional-exposure",
        "auto",
        "--seed",
        str(seed),
        "--threads",
        str(threads),
        "--emit-locus-stats",
        "--tsv",
    ]

    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    elapsed = time.time() - start
    log_path.write_text(
        "COMMAND\n" + " ".join(cmd) + "\n\nSTDOUT\n" + result.stdout + "\n\nSTDERR\n" + result.stderr
    )

    return {
        "condition": condition,
        "status": "ok" if result.returncode == 0 else "failed",
        "returncode": result.returncode,
        "out_dir": str(out_dir),
        "annotated_bam": str(annotated_bam),
        "log": str(log_path),
        "elapsed_sec": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-name", default=DEFAULT_OUT_NAME)
    parser.add_argument("--annotated-name", default=DEFAULT_ANNOTATED_NAME)
    parser.add_argument("--seed", type=int, default=20260518)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--conditions", nargs="*", default=None)
    args = parser.parse_args()

    base = args.base.resolve()
    index_dir = base / "rigel_index"
    conditions = load_conditions(base, args.conditions)
    if not conditions:
        raise SystemExit(f"No conditions found in {base}")
    if not index_dir.exists():
        raise SystemExit(f"Rigel index not found: {index_dir}")

    print(f"Synthetic base: {base}")
    print(f"Output name:    {args.out_name}")
    print(f"Conditions:     {len(conditions)}")

    run_rows: list[dict[str, object]] = []
    suite_start = time.time()
    for idx, condition in enumerate(conditions, start=1):
        print(f"[{idx:02d}/{len(conditions):02d}] {condition}", flush=True)
        row = run_condition(
            base,
            index_dir,
            condition,
            args.out_name,
            args.annotated_name,
            args.seed,
            args.threads,
            args.force,
        )
        run_rows.append(row)
        print(
            f"    {row['status']} rc={row['returncode']} elapsed={float(row['elapsed_sec']):.1f}s",
            flush=True,
        )
        if row["returncode"] != 0:
            print(f"    log: {row['log']}", flush=True)
            break

    summary = {
        "base": str(base),
        "out_name": args.out_name,
        "annotated_name": args.annotated_name,
        "seed": args.seed,
        "threads": args.threads,
        "elapsed_sec": time.time() - suite_start,
        "conditions": run_rows,
    }
    summary_path = base / f"{args.out_name}_run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Run summary: {summary_path}")

    return 0 if all(row["returncode"] == 0 for row in run_rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())