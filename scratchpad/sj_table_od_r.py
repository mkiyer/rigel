"""GATE 2 / GATE 3 — the fitted RNA strand overdispersion, OLD (boundary sides) vs NEW (SJ table).

Gate 2 (real cfRNA, the primary gate): ``od_r`` must land at **0.001–0.016** and fall OFF the 0.2
ceiling on all four libraries.  That range was measured independently, from deep junctions, before
this design existed (``strand_overdispersion_design.md`` §3) — if the table reproduces it, two
unrelated routes agree.

Gate 3 (synthetic): the suite has true ``od = 0`` by construction, so ``od_r`` must stay small.
It CANNOT validate a non-zero value; only real data can.

Also reports the deep-junction reference values the gate is stated against, recomputed here from
the same table but restricted to junctions deep enough to see the minority strand and re-centred on
their own mean, so a mean misfit cannot masquerade as dispersion.

    OMP_NUM_THREADS=1 python scratchpad/sj_table_od_r.py --cfrna --suite
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rigel.calibration.gdna_strand import (
    _MAX_OVERDISPERSION,
    _null_information,
    fit_rna_strand_from_sj_table,
    fit_rna_strand_overdispersion,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CFRNA = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")


def old_od_r(inp, kappa: float):
    """The RETIRED fit: boundary-side spliced counts (SPLICED_ANNOT + UNANNOT + IMPLICIT)."""
    ra = inp.get("region_arrays")
    if ra is None:
        ra = RegionArrays.from_region_df(_suite_index().region_df, _suite_index().ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(inp["payload"], ra)
    sense = np.concatenate(
        [
            np.asarray(sub.left.n_spliced_sense, dtype=np.float64),
            np.asarray(sub.right.n_spliced_sense, dtype=np.float64),
        ]
    )
    anti = np.concatenate(
        [
            np.asarray(sub.left.n_spliced_antisense, dtype=np.float64),
            np.asarray(sub.right.n_spliced_antisense, dtype=np.float64),
        ]
    )
    return fit_rna_strand_overdispersion(sense, anti + sense, kappa)


_IDX = {}


def _suite_index():
    from rigel.index import TranscriptIndex

    if "suite" not in _IDX:
        _IDX["suite"] = TranscriptIndex.load(str(SUITE / "rigel_index"))
    return _IDX["suite"]


def _raw_mom(sense, total, mu):
    """The unshrunk pooled MoM at mean ``mu`` — what the ceiling is hiding."""
    pq = mu * (1.0 - mu)
    num = float(np.sum((sense - total * mu) ** 2 - total * pq))
    den = float(np.sum(np.maximum(total * (total - 1.0), 0.0) * pq))
    return num / den if den > 0 else float("nan")


def deep_reference(table, min_depth: int):
    """The independent route: deep junctions only, re-centred on THEIR OWN mean."""
    n_se = np.asarray(table.n_sense, dtype=np.float64)
    n_tot = n_se + np.asarray(table.n_antisense, dtype=np.float64)
    m = n_tot >= min_depth
    if m.sum() == 0 or n_tot[m].sum() == 0:
        return float("nan"), 0
    mu = float(n_se[m].sum() / n_tot[m].sum())
    return _raw_mom(n_se[m], n_tot[m], mu), int(m.sum())


def report(label: str, inp: dict) -> None:
    sm = inp["strand_model"]
    table = sm.sj_table
    kappa = float(fit_strand_balance(sm).rna_sense_frac)

    old = old_od_r(inp, kappa)
    new = fit_rna_strand_from_sj_table(table, rna_sense_frac=kappa)

    n_se = np.asarray(table.n_sense, dtype=np.float64)
    n_tot = n_se + np.asarray(table.n_antisense, dtype=np.float64)
    raw_new = _raw_mom(n_se, n_tot, kappa)
    info = _null_information(n_tot, kappa * (1.0 - kappa))
    d100, k100 = deep_reference(table, 100)
    d1000, k1000 = deep_reference(table, 1000)

    ceil = _MAX_OVERDISPERSION
    flag_old = " CEILING" if abs(old.rna_strand_overdispersion - ceil) < 1e-9 else ""
    flag_new = " CEILING" if abs(new.rna_strand_overdispersion - ceil) < 1e-9 else ""
    print(
        f"  {label:<50} "
        f"kappa={kappa:.6f} | OLD od_r={old.rna_strand_overdispersion:.4f}{flag_old:<8} "
        f"(raw {old_raw(inp, kappa):>9.3f}, {old.n_seed_nodes:>7,} sides) | "
        f"NEW od_r={new.rna_strand_overdispersion:.5f}{flag_new:<8} "
        f"(raw {raw_new:>8.5f}, {new.n_seed_nodes:>7,} junctions, I={info:,.0f}) | "
        f"deep>=100 {d100:.5f} (n={k100}) >=1000 {d1000:.5f} (n={k1000})"
    )


def old_raw(inp, kappa):
    ra = inp.get("region_arrays") or RegionArrays.from_region_df(
        _suite_index().region_df, _suite_index().ref_name_to_id
    )
    sub = CalibrationSubstrate.from_payload(inp["payload"], ra)
    sense = np.concatenate(
        [
            np.asarray(sub.left.n_spliced_sense, dtype=np.float64),
            np.asarray(sub.right.n_spliced_sense, dtype=np.float64),
        ]
    )
    anti = np.concatenate(
        [
            np.asarray(sub.left.n_spliced_antisense, dtype=np.float64),
            np.asarray(sub.right.n_spliced_antisense, dtype=np.float64),
        ]
    )
    return _raw_mom(sense, anti + sense, kappa)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cfrna", action="store_true")
    ap.add_argument("--suite", action="store_true")
    args = ap.parse_args()

    if args.cfrna:
        print("\nGATE 2 — REAL cfRNA (target: 0.001–0.016, off the 0.2 ceiling)\n")
        for pkl in sorted(CFRNA.glob("*.pkl")):
            with open(pkl, "rb") as fh:
                report(pkl.stem, pickle.load(fh))

    if args.suite:
        print("\nGATE 3 — SYNTHETIC (true od = 0 by construction; must stay small)\n")
        vals = []
        for pkl in sorted((SUITE / "_selfsolve_cache").glob("*.pkl")):
            with open(pkl, "rb") as fh:
                inp = pickle.load(fh)
            report(pkl.stem, inp)
            sm = inp["strand_model"]
            k = float(fit_strand_balance(sm).rna_sense_frac)
            vals.append(
                (
                    old_od_r(inp, k).rna_strand_overdispersion,
                    fit_rna_strand_from_sj_table(
                        sm.sj_table, rna_sense_frac=k
                    ).rna_strand_overdispersion,
                )
            )
        o, n = np.array([v[0] for v in vals]), np.array([v[1] for v in vals])
        print(
            f"\n  synthetic OLD od_r: mean {o.mean():.5f} max {o.max():.5f} | "
            f"NEW od_r: mean {n.mean():.5f} max {n.max():.5f} | max|Δ| {np.abs(n - o).max():.5f}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
