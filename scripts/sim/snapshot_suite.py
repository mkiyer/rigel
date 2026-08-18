#!/usr/bin/env python
"""Deterministic digest of a simulator suite output — the refactor safety gate.

The `rigel.sim` simulator generates our test/benchmark ground truth, so any refactor of it
 must be **output-identical**. This tool runs a small
pinned suite and emits a normalized JSON digest (manifest condition identity + per-condition
truth-file hashes + per-condition oracle-BAM record digest), with all absolute paths / timestamps
stripped so the digest is location-independent and reproducible from a fixed seed.

⚠ **The suite it runs is its OWN pinned smoke suite, not the 16-condition ladder**, and its condition
names (`gdna_equal_ss_0.50_nrna_none`) come from that pinned profile. That is deliberate — the digest has
to be cheap and reproducible from a fixed seed — but it means a name printed here is NOT a panel rung and
must not be looked up in `sim/configs/gdna_ladder.yaml`.

⛔ **IT PARSED NO ARGUMENTS AT ALL UNTIL 2026-08-17**: every flag was silently discarded and
`--help` ran a full simulation and printed a digest, exiting 0. A tool that answers `--help` with a
measurement is indistinguishable from one that ignored what you asked it.

Usage:
    # before a refactor phase:
    python scripts/sim/snapshot_suite.py > /tmp/sim_before.json
    # after:
    python scripts/sim/snapshot_suite.py > /tmp/sim_after.json
    diff <(jq -S . /tmp/sim_before.json) <(jq -S . /tmp/sim_after.json)   # must be empty

    # digest a suite that already exists, instead of running the pinned one
    python scripts/sim/snapshot_suite.py --suite-dir ~/Downloads/rigel_runs/suite/ladder
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

# Pinned tiny suite — small + serial (workers=1) ⇒ deterministic, ~1s. Do not change these
# without regenerating both sides of any in-flight before/after comparison.
_SUITE_ARGS = [
    "--profile", "smoke",
    "--seed", "7",
    "--n-genes", "4",
    "--target-transcripts", "8",
    "--genome-length", "60000",
    "--n-rna", "1500",
]

# Manifest condition fields that define simulation identity (everything else is a path/derived).
_STABLE_CONDITION_FIELDS = (
    "name", "gdna_label", "gdna_rate",
    "gdna_strand_overdispersion", "gdna_strand_overdispersion_label",
    "strand_specificity", "nrna_label", "nrna_mode", "nrna_ratio",
    "capture_label", "capture_enabled",
    "n_mrna", "n_nrna", "n_rna", "n_gdna", "n_total", "seed",
    "n_mrna_observed", "n_nrna_observed", "n_gdna_observed",
)

# Per-condition data files whose raw bytes are path-free (pure tabular truth).
_TRUTH_FILES = ("truth_abundances.tsv", "truth_fragment_lengths.tsv")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _bam_digest(bam_path: Path) -> dict:
    """Record count + a content hash over sorted (qname, flag, ref, pos, cigar, seq) tuples.

    Order-independent (sorted) and header-independent — captures fragment identity, not byte layout
    (BAM headers carry the command line / are not part of the simulated content)."""
    rows: list[str] = []
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            rows.append(
                f"{r.query_name}\t{r.flag}\t{r.reference_name}\t{r.reference_start}\t"
                f"{r.cigarstring}\t{r.query_sequence}"
            )
    rows.sort()
    return {"n_records": len(rows), "sha256": _sha256_bytes("\n".join(rows).encode())}


def digest_suite_dir(suite_dir: Path) -> dict:
    """Normalized, location-independent digest of a suite output directory."""
    manifest = json.loads((suite_dir / "manifest.json").read_text())
    conditions = {}
    for cond in manifest.get("conditions", []):
        name = cond["name"]
        cdir = suite_dir / name
        entry = {k: cond.get(k) for k in _STABLE_CONDITION_FIELDS}
        for fname in _TRUTH_FILES:
            fpath = cdir / fname
            entry[fname] = _sha256_bytes(fpath.read_bytes()) if fpath.exists() else None
        bam = cdir / "sim_oracle.bam"
        entry["oracle_bam"] = _bam_digest(bam) if bam.exists() else None
        conditions[name] = entry
    return {
        "n_conditions": len(conditions),
        "conditions": dict(sorted(conditions.items())),
    }


def run_pinned_suite(outdir: Path) -> None:
    """⛔ The child's stderr is CAPTURED AND REPRINTED on failure, never discarded. It was sent to
    ``DEVNULL`` on both streams, so a simulator that raised surfaced here as a bare
    ``CalledProcessError`` with the reason thrown away — the digest's whole job is to notice a change in
    the simulator, and swallowing its error message is the one way to make that impossible."""
    cmd = [
        sys.executable, str(Path(__file__).parent / "simulate_suite.py"),
        "--outdir", str(outdir), *_SUITE_ARGS,
    ]
    done = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if done.returncode != 0:
        sys.stderr.write(done.stderr or "")
        raise SystemExit(f"⛔ the pinned suite failed (rc={done.returncode}): {' '.join(cmd)}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--suite-dir", type=Path, default=None,
                    help="digest THIS suite directory instead of running the pinned one")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the digest here instead of stdout")
    args = ap.parse_args()

    if args.suite_dir is not None:
        if not (args.suite_dir / "manifest.json").is_file():
            print(f"⛔ no manifest.json under {args.suite_dir} — that is not a simulator suite dir",
                  file=sys.stderr)
            return 2
        digest = digest_suite_dir(args.suite_dir)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "suite"
            run_pinned_suite(outdir)
            digest = digest_suite_dir(outdir)

    text = json.dumps(digest, indent=1, sort_keys=True) + "\n"
    if args.out is not None:
        args.out.write_text(text)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
