#!/usr/bin/env python
"""Deterministic digest of a simulator suite output — the refactor safety gate.

The `rigel.sim` simulator generates our test/benchmark ground truth, so any refactor of it
(see docs/CARRY_FORWARD.md) must be **output-identical**. This tool runs a small
pinned suite and emits a normalized JSON digest (manifest condition identity + per-condition
truth-file hashes + per-condition oracle-BAM record digest), with all absolute paths / timestamps
stripped so the digest is location-independent and reproducible from a fixed seed.

Usage:
    # before a refactor phase:
    python scripts/sim/snapshot_suite.py > /tmp/sim_before.json
    # after:
    python scripts/sim/snapshot_suite.py > /tmp/sim_after.json
    diff <(jq -S . /tmp/sim_before.json) <(jq -S . /tmp/sim_after.json)   # must be empty

`digest_suite_dir(path)` can also digest an existing suite dir for ad-hoc comparison.
"""

from __future__ import annotations

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
    cmd = [
        sys.executable, str(Path(__file__).parent / "simulate_suite.py"),
        "--outdir", str(outdir), *_SUITE_ARGS,
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(tmp) / "suite"
        run_pinned_suite(outdir)
        digest = digest_suite_dir(outdir)
    json.dump(digest, sys.stdout, indent=1, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
