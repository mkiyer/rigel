#!/usr/bin/env python3
"""Rerun Rigel on the synthetic 24-condition suite into a separate output name."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from rigel.sim.manifest import condition_manifest_map, load_manifest


DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")


def log_tail(path: Path, n_lines: int = 80) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n_lines:])


def run_condition(
    *,
    base: Path,
    condition: str,
    index_dir: Path,
    out_name: str,
    annotated_name: str,
    threads: int | None,
    force: bool,
) -> dict[str, object]:
    cond_dir = base / condition
    out_dir = cond_dir / out_name
    annotated_bam = cond_dir / annotated_name
    log_path = cond_dir / f"{out_name}.log"

    if force:
        if out_dir.exists():
            shutil.rmtree(out_dir)
        if annotated_bam.exists():
            annotated_bam.unlink()

    if (out_dir / "quant.feather").exists():
        return {
            "condition": condition,
            "status": "skipped",
            "elapsed_sec": 0.0,
            "out_dir": str(out_dir),
            "annotated_bam": str(annotated_bam),
            "log": str(log_path),
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
        "--emit-locus-stats",
        "--tsv",
    ]
    if threads is not None and threads > 0:
        cmd.extend(["--threads", str(threads)])

    started = time.time()
    with log_path.open("w") as log_handle:
        log_handle.write("COMMAND: " + " ".join(cmd) + "\n\n")
        log_handle.flush()
        proc = subprocess.run(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    elapsed = time.time() - started
    status = "ok" if proc.returncode == 0 else "failed"
    return {
        "condition": condition,
        "status": status,
        "returncode": proc.returncode,
        "elapsed_sec": elapsed,
        "out_dir": str(out_dir),
        "annotated_bam": str(annotated_bam),
        "log": str(log_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--out-name", default="rigel_efflen_out")
    parser.add_argument("--annotated-name", default="annotated_efflen.bam")
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--conditions", nargs="*", default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    base = args.base.resolve()
    manifest = load_manifest(base)
    condition_meta = condition_manifest_map(manifest)
    conditions = args.conditions or list(condition_meta)
    index_dir = base / "rigel_index"

    if not (index_dir / "transcripts.feather").exists():
        raise FileNotFoundError(f"Rigel index not found: {index_dir}")

    print(f"base: {base}")
    print(f"conditions: {len(conditions)}")
    print(f"output: <condition>/{args.out_name}")
    print(f"annotated BAM: <condition>/{args.annotated_name}")

    rows: list[dict[str, object]] = []
    all_started = time.time()
    for index, condition in enumerate(conditions, start=1):
        print(f"[{index:02d}/{len(conditions):02d}] {condition}", flush=True)
        row = run_condition(
            base=base,
            condition=condition,
            index_dir=index_dir,
            out_name=args.out_name,
            annotated_name=args.annotated_name,
            threads=args.threads,
            force=args.force,
        )
        rows.append(row)
        if row["status"] == "failed":
            print(f"  FAILED after {row['elapsed_sec']:.1f}s; log tail:")
            print(log_tail(Path(str(row["log"]))))
        else:
            print(f"  {row['status']} ({row['elapsed_sec']:.1f}s)", flush=True)

    summary = {
        "base": str(base),
        "out_name": args.out_name,
        "annotated_name": args.annotated_name,
        "threads": args.threads,
        "force": args.force,
        "elapsed_sec": time.time() - all_started,
        "conditions": rows,
    }
    summary_path = base / f"{args.out_name}_run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"wrote {summary_path}")

    failures = [row for row in rows if row["status"] == "failed"]
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
