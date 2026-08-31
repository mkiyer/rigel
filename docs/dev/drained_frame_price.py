#!/usr/bin/env python
"""ISSUES: instruments-calibrate-undrained-cache — PRICE the drained-vs-undrained frame gap.

⚠ A SANDBOX PROTOTYPE (docs/dev), kept so the pricing is reproducible while the frame ruling is
pending. The measured results live in the ISSUES entry; if the ruling lands, the migration replaces
this with real machinery and this file is deleted.

Per condition, two arms through calibration_vs_oracle.measure_condition — the release-metric
instrument's own scorer, unmodified:

  undrained  the instrument as it ships (calibrate on the pass-one cache; oracle = cached parts)
  drained    production's frame: the whole drained via pipeline._drain_side_buffer (production
             seed), each cached origin partition drained by REPLAYING the whole's choices
             (second_pass.lift_choices), OracleTruth.from_parts re-validating sum-to-full on the
             drained frame — the drain's own end-to-end identity gate, for free.

Also reported: how far the CERTIFIED TRUTH itself moves between frames (Σ|Δ| of the six override
masses), the drained-fragment count, and lift ambiguity.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import types
from pathlib import Path

import numpy as np

DESIGN = Path("/Users/mkiyer/proj/rigel/scripts/design")
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests")


def _sib(name):
    key = name[:-3]
    if key not in sys.modules:
        sp = importlib.util.spec_from_file_location(key, DESIGN / name)
        m = importlib.util.module_from_spec(sp)
        sys.modules[key] = m
        sp.loader.exec_module(m)
    return sys.modules[key]


CVO = _sib("calibration_vs_oracle.py")

from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.pipeline import _drain_side_buffer  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402
from rigel.second_pass import drain, lift_choices  # noqa: E402


def build_frames(suite: Path, oracle_cache: Path, cond: str, index, seed: int):
    """(undrained OracleTruth, drained OracleTruth, drained whole+strand_model, stats dict)."""
    cache = read_scan_cache(suite / "scan_cache" / cond, index)
    root = oracle_cache / cond
    parts_u = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    truth_u = OracleTruth.from_parts(cache.payload, parts_u)

    # ⛔ a FRESH load for the drain: drain()/_drain_side_buffer may share or mutate arrays, and
    # truth_u above holds references into the first load.
    cache2 = read_scan_cache(suite / "scan_cache" / cond, index)
    lift: dict = {}
    drained_whole = _drain_side_buffer(
        cache2.payload, index, cache2.strand_model, seed=seed, _lift=lift
    )
    stats = {
        "held": int(cache.payload.deferred.n_fragments),
        "deposited": 0,
        "n_ambiguous": 0,
    }
    if lift:
        # ⛔ fresh partition payloads: `parts_u` went into truth_u and drain() may mutate in place
        parts_f = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
        lifted, n_amb = lift_choices(
            lift["undrained"], [parts_f[k] for k in ORIGINS], lift["choices"]
        )
        parts_d = {
            k: drain(parts_f[k], ch, region_types=lift["region_types"], sj=lift["sj"])
            for k, ch in zip(ORIGINS, lifted)
        }
        stats["n_ambiguous"] = int(n_amb)
        stats["deposited"] = int(drained_whole.drain.deposited)
        # ⭐ THE FINDING THE FIRST RUN SURFACED: the whole-library drain resolves some TRUE-gDNA held
        # fragments to a SPLICED hypothesis, so in the drained tally the "gDNA never spliced" exact
        # zero is FALSE — production's certified-RNA channel carries a small drain-injected gDNA
        # contamination. Recorded per bank; the exact-zeros gate is degraded to this record for the
        # drained arm ONLY (sum-to-full, which runs first inside _validate, stays hard).
        import calibration._oracle as ORC

        stats["gdna_in_rna_only"] = {
            b: int(np.asarray(getattr(parts_d["gdna"], b), np.int64).sum())
            for b in ORC._RNA_ONLY_BANKS
        }
        orig_validate = OracleTruth._validate

        def tolerant(self):
            try:
                orig_validate(self)
            except AssertionError as e:
                if "gdna partition has" not in str(e):
                    raise

        OracleTruth._validate = tolerant
        try:
            truth_d = OracleTruth.from_parts(drained_whole, parts_d, n_ambiguous=n_amb)
        finally:
            OracleTruth._validate = orig_validate
    else:
        truth_d = truth_u
    return truth_u, truth_d, drained_whole, cache.strand_model, stats


def truth_movement(truth_u, truth_d, region_arrays):
    """Σ|Δ| of each override-mass field between the two frames, plus the totals."""
    mu = truth_u.override_masses(region_arrays)
    md = truth_d.override_masses(region_arrays)
    out = {}
    for k in mu:
        a, b = np.asarray(mu[k], np.float64), np.asarray(md[k], np.float64)
        out[k] = (float(np.abs(b - a).sum()), float(a.sum()), float(b.sum()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--conditions", nargs="+", required=True)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    region_arrays = RegionArrays.from_index(index)
    pc = PipelineConfig()
    oracle_cache = args.panel / "oracle_cache"

    orig_read, orig_load = CVO.read_scan_cache, CVO.load_oracle

    for cond in args.conditions:
        truth_u, truth_d, drained_whole, strand_model, st = build_frames(
            args.panel, oracle_cache, cond, index, pc.second_pass_seed
        )
        print(f"\n== {cond}: held {st['held']:,} → deposited {st['deposited']:,} "
              f"(lift ambiguity {st['n_ambiguous']})")
        if any(st.get("gdna_in_rna_only", {}).values()):
            print(f"   ⚠ drain-injected gDNA in the certified-RNA banks: {st['gdna_in_rna_only']}")

        mv = truth_movement(truth_u, truth_d, region_arrays)
        print("   certified-truth movement, Σ|Δ| (undrained total → drained total):")
        for k, (d, a, b) in mv.items():
            print(f"     {k:22s} {d:14,.1f}   ({a:14,.1f} → {b:14,.1f})")

        # arm 1: the instrument as it ships
        row_u = CVO.measure_condition(index, region_arrays, pc, args.panel, oracle_cache, cond)

        # arm 2: the drained frame — same measure_condition, two functions patched
        fired = {"read": 0, "load": 0}

        def read_patched(path, ix):
            p = Path(path)
            if p.name == cond and p.parent.name == "scan_cache":
                fired["read"] += 1
                return types.SimpleNamespace(payload=drained_whole, strand_model=strand_model)
            return orig_read(path, ix)

        def load_patched(suite, ocache, condition, ix, scan_payload):
            fired["load"] += 1
            return truth_d, "drained (lift-replayed cached parts)"

        CVO.read_scan_cache, CVO.load_oracle = read_patched, load_patched
        try:
            row_d = CVO.measure_condition(index, region_arrays, pc, args.panel, oracle_cache, cond)
        finally:
            CVO.read_scan_cache, CVO.load_oracle = orig_read, orig_load
        if fired["read"] != 1 or fired["load"] != 1:
            raise RuntimeError(f"{cond}: patches fired read={fired['read']} load={fired['load']}")
        if row_u["noop_differences"] or row_d["noop_differences"]:
            raise RuntimeError(f"{cond}: noop gate failed "
                               f"{row_u['noop_differences']} / {row_d['noop_differences']}")

        print(f"   {'':24s} {'UNDRAINED (ships)':>22s} {'DRAINED (production)':>22s} {'delta':>12s}")
        for axis in row_u["axes"]:
            au = row_u["axes"][axis]["abs_err"]
            ad = row_d["axes"][axis]["abs_err"]
            print(f"   {axis + ' abs_err(P vs O)':24s} {au:22,.1f} {ad:22,.1f} {ad - au:+12,.1f}"
                  f"   ({(ad - au) / au * 100 if au else 0:+.2f} %)")
        for name in ("library_f_gdna_P", "library_f_gdna_O"):
            print(f"   {name:24s} {row_u[name]:22.6f} {row_d[name]:22.6f} "
                  f"{row_d[name] - row_u[name]:+12.6f}")
        for name in ("ruler_abs_err", "ruler_n_moved"):
            print(f"   {name:24s} {row_u[name]:22,.1f} {row_d[name]:22,.1f} "
                  f"{row_d[name] - row_u[name]:+12,.1f}")
        pu, pd = row_u["pools"], row_d["pools"]
        for k in ("P_gdna", "O_gdna", "P_rna", "O_rna"):
            print(f"   pool {k:19s} {pu[k]:22,.1f} {pd[k]:22,.1f} {pd[k] - pu[k]:+12,.1f}")
        print(f"   pool true_gdna/true_rna  {pu['true_gdna']:,} / {pu['true_rna']:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
