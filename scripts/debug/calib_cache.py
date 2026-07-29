"""Scan-once / cache-everything harness for calibration iteration.

Scanning a real BAM (the vcap DNA-mix is 3.1 GB) is the expensive step; the calibration phase we are
iterating on is milliseconds. This caches EVERY input `calibrate()` needs from a single BAM scan, so the
pass-0 KDE / prior / precision work can re-run calibration hundreds of times without re-scanning.

    # scan once -> cache
    python scripts/debug/calib_cache.py build --bam X.bam --index IDX --out cache.pkl [--verify]
    # then, in any experiment:
    from calib_cache import load, calibrate_from_cache
    inp = load("cache.pkl"); res = calibrate_from_cache(inp)          # ms, no scan

The cache is a pickle of the exact `calibrate()` kwargs (payload, region_arrays, strand_model, the two FL
pmfs) + provenance (bam/index/sj_tag). `--verify` proves a cached calibrate == a fresh scan+calibrate.
"""

from __future__ import annotations

import argparse
import os
import pickle
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.splice_graph import build_boundary_flags_array
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType


def build(bam: str, index_dir: str, out: str, cfg: PipelineConfig | None = None) -> dict:
    """Scan ``bam`` against ``index_dir`` ONCE and pickle every calibration input to ``out``."""
    cfg = cfg or PipelineConfig()
    index = TranscriptIndex.load(index_dir)
    ra = RegionArrays.from_index(index)
    sj_tag = _native_detect_sj_tag(bam)
    sc = dc(cfg.scan, sj_strand_tag=sj_tag)
    _stats, sm, flm, _buf, payload = scan_and_buffer(bam, index, sc)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    inp = dict(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=np.asarray(fl.gdna_pmf),
        rna_fl_pmf=np.asarray(fl.rna_pmf),
        sj_tag=sj_tag,
        bam=str(bam),
        index_dir=str(index_dir),
        # ⚠ PROVENANCE (2026-07-29). `index_dir` was already stored and never checked, so a cache built
        # against one partition loaded silently against another — and `region_arrays` is pickled ALONGSIDE
        # the payload, so the two stayed mutually consistent while both were stale relative to the index on
        # disk. `load(..., index=...)` now verifies this hash. See `TranscriptIndex.partition_hash`.
        partition_hash=index.partition_hash,
    )
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as fh:
        pickle.dump(inp, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"cached -> {out}  ({os.path.getsize(out) / 1e6:.1f} MB)")
    return inp


def load(path: str, index: "TranscriptIndex | None" = None) -> dict:
    """Load a calibration-input cache (returns the dict of `calibrate()` kwargs + provenance).

    Pass ``index`` to VERIFY the cache was built against that index's partition; mismatch raises. A cache
    predating the ``partition_hash`` field cannot be verified and also raises when an index is supplied —
    rebuild it rather than trust it.
    """
    with open(path, "rb") as fh:
        inp = pickle.load(fh)
    if index is not None:
        want, got = index.partition_hash, inp.get("partition_hash")
        if got != want:
            raise ValueError(
                f"{path} was built against partition {got or '<unrecorded>'} but the supplied index is "
                f"{want}. The payload's region/boundary keying is not valid for this index — rebuild "
                f"the cache (calib_cache.py build --bam ... --index ...)."
            )
    return inp


def calibrate_from_cache(
    inp: dict,
    cfg: PipelineConfig | None = None,
    _debug: dict | None = None,
    index: "TranscriptIndex | None" = None,
):
    """Run production `calibrate()` from a cache dict — the fast inner loop (no scan).

    Pass ``index`` to supply the graph's structural flags (W1c). They are derived here rather than
    cached: see the note beside the call.
    """
    cfg = cfg or PipelineConfig()
    return calibrate(
        payload=inp["payload"],
        region_arrays=inp["region_arrays"],
        strand_model=inp["strand_model"],
        gdna_fl_pmf=inp["gdna_fl_pmf"],
        rna_fl_pmf=inp["rna_fl_pmf"],
        config=cfg.calibration,
        # ⚠ DERIVED, never cached. `partition_hash` covers nodes.feather only, so a cached flags
        # array could go stale against a rebuilt edges.feather and verify CLEAN — which is exactly
        # what the 2026-07-29 flag-filter fix would have done (it rewrote every edge file and left
        # every node file byte-identical). Building them is one pass over edges; the cache exists to
        # skip the BAM scan, not that.
        boundary_flags=None if index is None else build_boundary_flags_array(index),
        _debug=_debug,
    )


def _verify(bam: str, index_dir: str, inp: dict, cfg: PipelineConfig) -> None:
    """Prove the cache round-trips: a cached calibrate matches a fresh scan+calibrate."""
    index = TranscriptIndex.load(index_dir)
    ra = RegionArrays.from_index(index)
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _s, sm, flm, _b, payload = scan_and_buffer(bam, index, sc)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    fresh = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=cfg.calibration,
    )
    cached = calibrate_from_cache(inp, cfg)
    for attr in ("rna_sense_frac", "gdna_density_global"):
        a, b = getattr(fresh, attr), getattr(cached, attr)
        ok = np.allclose(a, b, rtol=1e-9, atol=1e-12)
        print(f"  verify {attr}: fresh={a:.8g} cached={b:.8g}  match={ok}")
        assert ok, f"cache MISMATCH on {attr}"
    print("  verify: cache round-trips (cached calibrate == fresh scan+calibrate)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="scan a BAM once and cache all calibration inputs")
    b.add_argument("--bam", required=True)
    b.add_argument("--index", required=True)
    b.add_argument("--out", required=True)
    b.add_argument(
        "--verify", action="store_true", help="re-scan and confirm the cache round-trips"
    )
    args = ap.parse_args()
    cfg = PipelineConfig()
    if args.cmd == "build":
        inp = build(args.bam, args.index, args.out, cfg)
        if args.verify:
            _verify(args.bam, args.index, inp, cfg)


if __name__ == "__main__":
    main()
