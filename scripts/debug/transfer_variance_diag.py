"""What sigma^2_transfer ACTUALLY is — measured against the oracle, per pair class.

`sigma^2_transfer` is the adjacent-pair overdispersion in a BP message::

    Var(log rho_c^src as evidence about rho_c^dst) = Var(log f_c^src) + 1/n_c^src + sigma^2_transfer
                                                     \\___ strand ___/   \\_ count _/   \\__ THIS __/

It is a property of the WORLD (do adjacent nodes have the same rate?), not of our knowledge. It was zeroed in
priority #2 because the shipped estimator (`bp_solver.adjacent_disagreement_variance`, the TOTAL-density
adjacent disagreement) is a different quantity: on a capture-OFF library the true gDNA transfer variance is
~0 while the total-density estimator reports ~1.0 nats^2 of pure RNA-expression variation, which is what
gagged every message.

This tool measures, on the ORACLE, what the honest value would be:

  * TRUE  gDNA / RNA / total adjacent transfer variance, on the corrected geometry
  * the RNA-FREE ANCHOR: pairs where RNA is absent, so rho_total == rho_g IDENTICALLY (belief-free, no
    location read, no assertion) — a candidate belief-free estimator
  * the three-bullet STRATIFICATION (depleted-depleted / enriched-enriched / MIXED), the user's own spec
  * the production sigma^2_imp for reference

THE ZERO CONVENTION — read this before quoting any number from this tool.
======================================================================
A transfer variance is a variance of LOG rates, so nodes with zero mass are `log 0 = -inf` and something must
be done with them. **The choice is worth 20x** and every number is meaningless without it. Measured on
gdna300_ss0.50_none_capON, flooring the zero-mass nodes at an arbitrary `eps`::

    eps        1e-1   1e-2   1e-3   1e-4   1e-5   1e-6   1e-9   1e-12
    sigma^2    4.92   7.46  10.02  14.12  19.76   26.9   57.7   102.4

An `eps` floor is not a convention, it is a free parameter. This tool uses neither an `eps` floor nor a
"drop the zeros" filter. Both are wrong for the same reason: **they gate on COUNT when the physics gates on
OPPORTUNITY.**

  * `eff == 0` — the node is SHORTER THAN THE SHORTEST FRAGMENT. It cannot observe anything; its zero is not
    a measurement of a rate, it is the absence of an experiment. **EXCLUDED.** (Verified invariant: `eff == 0`
    => `mass == 0` on every node, both faces, all six eff-lengths — so nothing is discarded that had data.
    `tests/calibration/test_node_chain.py::test_zero_opportunity_nodes_carry_zero_mass`.)
  * `eff > 0` and `mass == 0` — a genuine Poisson zero over real exposure. It IS a measurement: "my rate is
    below the minimum I could have observed". **KEPT**, with the rate floored at the MINIMUM OBSERVABLE
    DENSITY `1/eff` — i.e. `rho = max(mass, 1) / eff`, one fragment's worth. Not a free parameter: it is the
    experiment's own resolution limit, and it is production's existing convention (`bp_solver._scan` floors
    its message at `max(rho, 1.0/egd)`).

The true zeros are REAL and are kept. Dropping them would discard the depleted tail, which on a capture
library is exactly the gDNA signal.

⚠ Every number here divides by an effective length. `region_eff_length` was OFF BY ONE
(`E[max(0,L-l)]` for a discrete `max(0,L-l+1)` start-position count), which made it EXACTLY ZERO whenever a
region equalled the shortest fragment -> _EPS-floored -> ~1e9 densities on 211/1698 regions. That is fixed;
measured, it moves the transfer variance by 0.1% (6.826 -> 6.819) — it wrecked two NPMLE grids, not this.

⚠ Poisson removal uses the oracle's FRACTIONAL mass as the count `n`. The Kish `n_eff = (sum w)^2/sum w^2 >=
mass`, so `1/mass` OVER-states the Poisson term and the removed variances are LOWER bounds. Raw values are
reported alongside so the subtraction is visible.

    OMP_NUM_THREADS=1 python scripts/debug/transfer_variance_diag.py [--suite DIR] [--conditions a,b]
                                                                     [--zero-convention opportunity|drop|eps]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from selfsolve_diag import _scan_and_truth

from rigel.calibration.bp_solver import (
    _adjacent_edges,
    _poisson_moment_var,
    adjacent_disagreement_variance,
    node_global_geometry,
)
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def _node_truth(chain, inp, eff, convention="opportunity", eps=1e-6):
    """Per-chain-node ORACLE gDNA / RNA / total mass, and the rates in the node_global_geometry frame.

    ``convention`` decides what a zero-mass node's log-rate is (module docstring — it is worth 20x):
      * ``opportunity`` (default, the derived one): rate = ``max(mass, 1)/eff`` — floored at the node's own
        MINIMUM OBSERVABLE DENSITY, one fragment's worth. Zero-``eff`` nodes are excluded by the caller.
      * ``drop``: zero-mass nodes excluded entirely (gates on COUNT — discards the real depleted tail).
      * ``eps``: rate floored at a free parameter (what produced 17.129; shown only to expose the range).
    """
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pools(p):
        g = np.asarray(p["gdna_pos"], float) + np.asarray(p["gdna_neg"], float)
        r = sum(
            np.asarray(p[k], float)
            for k in ("mat_uns_pos", "mat_uns_neg", "nas_uns_pos", "nas_uns_neg")
        )
        return g, r

    gR, rR = pools(inp["region_pools"])
    gB, rB = pools(inp["boundary_pools"])
    ri = np.clip(idx, 0, gR.shape[0] - 1)
    bi = np.clip(idx, 0, gB.shape[0] - 1)
    g = np.where(isr, gR[ri], gB[bi])
    r = np.where(isr, rR[ri], rB[bi])
    e = np.maximum(eff, _EPS)

    def rate(mass):
        if convention == "opportunity":
            return np.maximum(mass, 1.0) / e  # min-observable floor: one fragment
        if convention == "eps":
            return np.maximum(mass / e, eps)
        return mass / e  # 'drop' — the caller masks mass == 0

    return g, r, rate(g), rate(r), rate(g + r)


def _tv(lo, src, dst, mask, n_src, n_dst):
    """Transfer variance of a log-rate over a pair mask: raw, and with Poisson sampling removed."""
    if int(mask.sum()) < 2:
        return np.nan, np.nan, 0
    d = lo[dst][mask] - lo[src][mask]
    return (
        float(np.var(d)),
        float(_poisson_moment_var(d, n_src[src][mask], n_dst[dst][mask])),
        int(mask.sum()),
    )


def measure(inp, index, cfg, convention="opportunity", eps=1e-6):
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pl = inp["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    _mass, eff = node_global_geometry(chain, geom)
    g, r, rho_g, rho_r, rho_t = _node_truth(chain, inp, eff, convention=convention, eps=eps)

    def lg(x):
        return np.where(x > 0, np.log(np.maximum(x, 1e-300)), 0.0)

    Lg, Lr, Lt = lg(rho_g), lg(rho_r), lg(rho_t)
    src, dst, _ = _adjacent_edges(chain)
    # THE GATE IS ON OPPORTUNITY, NOT ON COUNT. `eff > _EPS` == "this node could observe something".
    # A zero-mass node with real eff is KEPT: its zero is a genuine measurement, floored at 1/eff.
    # ⚠ The gate MUST be against `node_geometry._EPS`, not a local epsilon. `build_node_geometry` floors
    # every eff at 1e-9 AT CONSTRUCTION (node_geometry.py:209-214), so a zero-opportunity node arrives here
    # already wearing eff=1e-9, indistinguishable from a tiny one. Gating at a SMALLER epsilon readmits it
    # with a floored rate of 1/1e-9 = 1e9 — which measured sigma^2 = 69.06 instead of 6.82 and looked like a
    # real result. Same 1e9 pathology as the region_eff_length off-by-one, through a different door: the
    # _EPS floor erases "cannot observe" vs "barely observes", and every consumer picking its own epsilon
    # re-opens it.
    # "this face could observe something" — gate ABOVE the node_geometry floor (1e-9), NOT this script's
    # local 1e-12: gating below the geometry floor re-admits floored zero-opportunity nodes at ~1e9 (the whole
    # point of the former `has_opportunity` helper, now removed).
    obs = np.asarray(eff, float) > 1e-9 * 1.001
    live = obs[src] & obs[dst]
    if convention == "drop":  # the count-gated alternative, for contrast only
        live &= (g[src] > 0) & (g[dst] > 0)
    src, dst = src[live], dst[live]
    alle = np.ones(src.shape[0], bool)

    # enrichment strata, from the TRUE gDNA rate (an oracle read — this is the CEILING a belief-free
    # stratifier could reach, not something production can do)
    med = np.median(Lg[np.concatenate([src, dst])])
    enr_s, enr_d = Lg[src] > med, Lg[dst] > med
    rows = {
        "ALL": alle,
        "dep-dep": (~enr_s) & (~enr_d),
        "enr-enr": enr_s & enr_d,
        "MIXED": enr_s ^ enr_d,
        # the belief-free candidate: RNA absent at BOTH endpoints => rho_total == rho_g identically
        "RNA-free anchor": (r[src] <= 1e-9) & (r[dst] <= 1e-9),
    }
    out = []
    for name, m in rows.items():
        vg_raw, vg_p, n = _tv(Lg, src, dst, m, g, g)
        vt_raw, vt_p, _ = _tv(Lt, src, dst, m, g + r, g + r)
        vr_raw, vr_p, _ = _tv(Lr, src, dst, m, np.maximum(r, _EPS), np.maximum(r, _EPS))
        out.append(
            dict(
                pair_class=name,
                n_pairs=n,
                gdna_raw=vg_raw,
                gdna_pois=vg_p,
                total_raw=vt_raw,
                total_pois=vt_p,
                rna_raw=vr_raw,
                rna_pois=vr_p,
                total_over_gdna=(vt_raw / vg_raw) if vg_raw > _EPS else np.inf,
            )
        )
    return pd.DataFrame(out), float(adjacent_disagreement_variance(chain, geom))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--zero-convention",
        default="opportunity",
        choices=["opportunity", "drop", "eps"],
        help="how a zero-mass node's log-rate is defined; see the module docstring (worth 20x)",
    )
    ap.add_argument("--eps", type=float, default=1e-6, help="only for --zero-convention eps")
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))
    allr = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        df, s2imp = measure(inp, index, cfg, convention=a.zero_convention, eps=a.eps)
        print(
            f"\n{'=' * 104}\n{c}\n   zero convention = {a.zero_convention}"
            f"   |   production sigma^2_imp (total-density adjacent disagreement) = {s2imp:.4f}\n{'=' * 104}"
        )
        print(
            df.to_string(
                index=False,
                float_format=lambda x: f"{x:,.4g}",
                columns=[
                    "pair_class",
                    "n_pairs",
                    "gdna_raw",
                    "gdna_pois",
                    "total_raw",
                    "total_pois",
                    "rna_raw",
                    "total_over_gdna",
                ],
            )
        )
        allr.append(df.assign(condition=c, sigma2_imp=s2imp, zero_convention=a.zero_convention))
    res = pd.concat(allr, ignore_index=True)
    dest = Path(a.out) if a.out else suite / "transfer_variance_diag.tsv"
    res.to_csv(dest, sep="\t", index=False)
    print(f"\nfull table -> {dest}")
    print(
        "\ngdna_pois/total_pois: Poisson-removed via the production _poisson_moment_var. The oracle's"
    )
    print(
        "counts are FRACTIONAL mass and Kish n_eff >= mass, so 1/mass over-subtracts: these are LOWER"
    )
    print(
        "bounds. 'total_over_gdna' > 1 means the total-density stand-in is CONSERVATIVE for gDNA."
    )


if __name__ == "__main__":
    main()
