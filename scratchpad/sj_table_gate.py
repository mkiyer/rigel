"""S1 GATE — the per-junction SJ strand table's marginal 2x2 must EQUAL the existing 2x2, exactly.

This is the whole correctness argument of docs/calibration/sj_strand_table_design.md:

    sum over motif-POS junctions of n_sense      == pos_pos
    sum over motif-POS junctions of n_antisense  == neg_pos
    sum over motif-NEG junctions of n_sense      == neg_neg
    sum over motif-NEG junctions of n_antisense  == pos_neg

plus the independent count check the scanner supplies: total table depth == n_strand_trained,
i.e. one qualified fragment credits exactly one junction.

Runs a production scan (no calibration) and compares.  Usage:

    OMP_NUM_THREADS=1 python scratchpad/sj_table_gate.py --bam X.bam --index IDX
    OMP_NUM_THREADS=1 python scratchpad/sj_table_gate.py --suite   # all 32 synthetic conditions
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

import rigel.pipeline as pipe
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CFRNA = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna")
STRAND_POS, STRAND_NEG = 1, 2

_INDEX_CACHE: dict[str, TranscriptIndex] = {}


def raw_scan(bam: str, index_dir: str) -> dict:
    """Production scan; return the raw ``strand_observations`` dict from C++.

    Captured by wrapping ``StrandModels.from_scan`` so the scan runs through exactly the
    production path (``scan_and_buffer``) — the gate must test what ships, not a
    reimplementation of it.
    """
    index = _INDEX_CACHE.get(index_dir) or _INDEX_CACHE.setdefault(
        index_dir, TranscriptIndex.load(index_dir)
    )
    cfg = PipelineConfig()
    sc = dc(cfg.scan, sj_strand_tag=pipe._native_detect_sj_tag(bam))
    captured: dict = {}
    from rigel.strand_model import StrandModels

    orig = StrandModels.from_scan

    def _capture(strand_dict):
        captured.update(strand_dict)
        return orig(strand_dict)

    StrandModels.from_scan = staticmethod(_capture)
    try:
        stats, sm, _flm, _buf, _pl = pipe.scan_and_buffer(bam, index, sc)
    finally:
        StrandModels.from_scan = orig
    captured["_n_strand_trained"] = stats.n_strand_trained
    captured["_model"] = sm
    return captured


def check(label: str, sd: dict) -> bool:
    motif = np.asarray(sd["sj_motif_strand"])
    n_se = np.asarray(sd["sj_n_sense"]).astype(np.int64)
    n_as = np.asarray(sd["sj_n_antisense"]).astype(np.int64)
    mp, mn = motif == STRAND_POS, motif == STRAND_NEG
    got = {
        "pos_pos": int(n_se[mp].sum()),
        "neg_pos": int(n_as[mp].sum()),
        "neg_neg": int(n_se[mn].sum()),
        "pos_neg": int(n_as[mn].sum()),
    }

    model = sd["_model"].exonic_spliced
    ref = {
        "pos_pos": model.pos_pos,
        "pos_neg": model.pos_neg,
        "neg_pos": model.neg_pos,
        "neg_neg": model.neg_neg,
    }
    ok = ref == got and model.contingency_matches_table()
    n_j = int(motif.size)
    n_obs = sum(ref.values())
    kappa = (ref["pos_pos"] + ref["neg_neg"]) / n_obs if n_obs else 0.5
    depth = n_se + n_as
    print(
        f"  {label:<52} {'PASS' if ok else 'FAIL':>4} "
        f"n_obs={n_obs:>9,} junctions={n_j:>8,} kappa={kappa:.6f} "
        f"depth med={np.median(depth) if n_j else 0:.0f} max={depth.max() if n_j else 0:.0f} "
        f"d>=100={int((depth >= 100).sum()):>6,}"
    )
    if not ok:
        for k in ("pos_pos", "pos_neg", "neg_pos", "neg_neg"):
            flag = "" if ref[k] == got[k] else "   <-- MISMATCH"
            print(f"      {k}: 2x2={ref[k]:>10,}  table={got[k]:>10,}{flag}")
    # ⭐ THE INVARIANT: the scanner counts qualified fragments independently of the table
    # it builds, so this is a genuine cross-check that one fragment credits one junction.
    assert int(depth.sum()) == sd["_n_strand_trained"], (
        f"{label}: table total {depth.sum()} != n_strand_trained {sd['_n_strand_trained']}"
    )
    assert int(depth.sum()) == n_obs, f"{label}: table total {depth.sum()} != 2x2 total {n_obs}"
    # Motif strand must be POS or NEG only (qualification guarantees it).
    assert int(np.count_nonzero(~(mp | mn))) == 0, f"{label}: motif strand not in POS/NEG"
    # Sorted, unique keys (determinism).
    key = np.stack(
        [
            np.asarray(sd["sj_ref_id"]).astype(np.int64),
            np.asarray(sd["sj_start"]),
            np.asarray(sd["sj_end"]),
            motif.astype(np.int64),
        ]
    )
    order = np.lexsort(key[::-1])
    assert np.array_equal(order, np.arange(n_j)), f"{label}: junction keys not sorted"
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam")
    ap.add_argument("--index")
    ap.add_argument("--suite", action="store_true")
    ap.add_argument("--cfrna", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    jobs: list[tuple[str, str, str]] = []
    if args.bam:
        jobs.append((Path(args.bam).name, args.bam, args.index))
    if args.suite:
        idx = str(SUITE / "rigel_index")
        conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
        jobs += [(c, str(SUITE / c / "sim_oracle.bam"), idx) for c in conds]
    if args.cfrna:
        idx = "/Users/mkiyer/Downloads/rigel_runs/rigel_index"
        for d in sorted(CFRNA.iterdir()):
            b = d / "bam" / "star.srt.rmdup.collate.bam"
            if b.exists():
                jobs.append((d.name, str(b), idx))
    if args.limit:
        jobs = jobs[: args.limit]

    print(f"S1 GATE: table marginal 2x2 == existing 2x2 ({len(jobs)} datasets)\n")
    n_ok = 0
    for label, bam, index_dir in jobs:
        n_ok += bool(check(label, raw_scan(bam, index_dir)))
    print(f"\n{n_ok}/{len(jobs)} PASS")
    return 0 if n_ok == len(jobs) else 1


if __name__ == "__main__":
    sys.exit(main())
