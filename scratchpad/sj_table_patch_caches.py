"""Refresh the pickled ``strand_model`` in the development calibration caches.

The caches (``_selfsolve_cache``, ``_calib_cache``) hold a scanned payload, oracle truth pools and
a trained ``StrandModels``.  The payload and the truth are unchanged by the SJ strand table work,
but the strand model now carries the per-junction table the RNA overdispersion is fitted from — so
only that field needs refreshing, and re-deriving the oracle truth would be pure waste.

Each patch RE-SCANS the BAM and asserts the 2×2 is bit-identical to what the cache already held.
That is Gate 1 (κ bit-identical) run over every cached dataset as a side effect of the patch.

    OMP_NUM_THREADS=1 python scratchpad/sj_table_patch_caches.py --suite --cfrna [--dry-run]
"""

from __future__ import annotations

import argparse
import os
import pickle
import shutil
import sys
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import rigel.pipeline as pipe
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CFRNA = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna")
CFRNA_INDEX = "/Users/mkiyer/Downloads/rigel_runs/rigel_index"

_INDEX_CACHE: dict[str, TranscriptIndex] = {}


def _index(path: str) -> TranscriptIndex:
    if path not in _INDEX_CACHE:
        _INDEX_CACHE[path] = TranscriptIndex.load(path)
    return _INDEX_CACHE[path]


def _cells(sm) -> tuple[int, int, int, int]:
    m = sm.exonic_spliced
    return (m.pos_pos, m.pos_neg, m.neg_pos, m.neg_neg)


def patch(pkl: Path, bam: str, index_dir: str, dry_run: bool) -> bool:
    with open(pkl, "rb") as fh:
        inp = pickle.load(fh)
    before = _cells(inp["strand_model"])

    cfg = PipelineConfig()
    sc = dc(cfg.scan, sj_strand_tag=pipe._native_detect_sj_tag(bam))
    _stats, sm, _flm, _buf, _payload = pipe.scan_and_buffer(bam, _index(index_dir), sc)
    after = _cells(sm)

    ok = before == after and sm.exonic_spliced.contingency_matches_table()
    table = sm.sj_table
    print(
        f"  {pkl.stem:<52} {'OK ' if ok else 'FAIL'} "
        f"2x2={after}  n_obs={sm.n_observations:,}  junctions={table.n_junctions:,}"
    )
    if not ok:
        print(f"      cached 2x2 = {before}")
        print(f"      rescan 2x2 = {after}")
        return False
    if not dry_run:
        inp["strand_model"] = sm
        tmp = pkl.with_suffix(".pkl.tmp")
        with open(tmp, "wb") as fh:
            pickle.dump(inp, fh, protocol=pickle.HIGHEST_PROTOCOL)
        shutil.move(str(tmp), str(pkl))
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", action="store_true")
    ap.add_argument("--cfrna", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    jobs: list[tuple[Path, str, str]] = []
    if args.suite:
        idx = str(SUITE / "rigel_index")
        for pkl in sorted((SUITE / "_selfsolve_cache").glob("*.pkl")):
            jobs.append((pkl, str(SUITE / pkl.stem / "sim_oracle.bam"), idx))
    if args.cfrna:
        stem_to_dir = {d.name.split("_")[1]: d for d in CFRNA.iterdir() if d.name.startswith("mctp_")}
        for pkl in sorted((CFRNA / "_calib_cache").glob("*.pkl")):
            with open(pkl, "rb") as fh:
                bam = pickle.load(fh)["bam"]
            jobs.append((pkl, bam, CFRNA_INDEX))

    print(f"PATCH {len(jobs)} caches{' (DRY RUN)' if args.dry_run else ''}\n")
    n_ok = sum(patch(*j, args.dry_run) for j in jobs)
    print(f"\n{n_ok}/{len(jobs)} OK")
    return 0 if n_ok == len(jobs) else 1


if __name__ == "__main__":
    sys.exit(main())
