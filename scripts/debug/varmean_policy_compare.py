"""Compare var~mean refit GATING policies on precision (zero-gDNA phantom) AND recall (gDNA recovery).

The symmetric gDNA refit on ALL nodes is circular (learns the zero-gDNA+nascent exon phantom). But seed-only
misses the real-data gold signal (UNEXPRESSED single-strand exons = pure gDNA-on-exon, where capture enrichment
is strongest). Compromise: gate to SINGLE-STRAND nodes (strand_obs = free_pos^free_neg) — definitively
strand-solvable, includes single-strand exons, excludes the hard AMBIG nodes.

Policies (gDNA gate / RNA gate):
  all   : gDNA=all nodes,      RNA=free_s            (current nested default)
  ss_g  : gDNA=strand_obs,     RNA=free_s            (single-strand gDNA only)
  ss    : gDNA=strand_obs,     RNA=free_s&strand_obs (fully symmetric single-strand)

Metrics: zero-gDNA+nascent → total phantom gDNA (lower better). gDNA-present → assigned vs oracle (recall).
NOTE the sim has NO unexpressed genes, so it CANNOT show the single-strand policy's main real-data benefit
(training gDNA on unexpressed exons). This measures the precision/recall tradeoff on EXPRESSED data only.

    OMP_NUM_THREADS=1 python scripts/debug/varmean_policy_compare.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import _edge_varmean, _fit_offset  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

ZERO = "gdna_none_ss_0.50_nrna_rnd_capture_off"   # precision: all assigned gDNA is phantom
PRES = "gdna_gdna300_ss_0.50_nrna_none_capture_on"  # recall: flagship gDNA-present (the -11% win)
MAX_OUTER = 4


def _gdna_fit(gate):
    def f(chain, densities, geometry, statics):
        n = np.asarray(chain.kind).shape[0]
        live = np.ones(n, bool) if gate == "all" else np.asarray(statics.strand_obs, bool)
        m, r, o = _edge_varmean(chain, densities.rho_g_left, densities.rho_g_right,
                                geometry.eff_gdna_left, geometry.eff_gdna_right, live)
        return _fit_offset(m, r, o)
    return f


def _rna_fit(gate):
    def f(chain, densities, geometry, statics):
        fp, fn = np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool)
        so = np.asarray(statics.strand_obs, bool)
        lp, ln = (fp & so, fn & so) if gate == "ss" else (fp, fn)
        mp, rp, op = _edge_varmean(chain, densities.rho_pos_left, densities.rho_pos_right,
                                   geometry.eff_rna_left, geometry.eff_rna_right, lp)
        mn, rn, on = _edge_varmean(chain, densities.rho_neg_left, densities.rho_neg_right,
                                   geometry.eff_rna_left, geometry.eff_rna_right, ln)
        return _fit_offset(mp + mn, rp + rn, op + on)
    return f


POLICIES = {
    "all  (gDNA=all,  RNA=free_s)": ("all", "free"),
    "ss_g (gDNA=1str, RNA=free_s)": ("ss", "free"),
    "ss   (gDNA=1str, RNA=1str)  ": ("ss", "ss"),
}


def run(blob, ra, gpol, rpol):
    cfg = dataclasses.replace(CalibrationConfig(), sweep_max_outer=MAX_OUTER)
    og, orf = bp.fit_gdna_varmean, bp.fit_rna_varmean
    bp.fit_gdna_varmean, bp.fit_rna_varmean = _gdna_fit(gpol), _rna_fit(rpol)
    try:
        cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                        gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=cfg)
    finally:
        bp.fit_gdna_varmean, bp.fit_rna_varmean = og, orf
    return cal


def main():
    idx_z, blob_z = build_or_load_cache(ZERO, False)
    ra_z = RegionArrays.from_region_df(idx_z.region_df, idx_z.ref_name_to_id)
    idx_p, blob_p = build_or_load_cache(PRES, False)
    ra_p = RegionArrays.from_region_df(idx_p.region_df, idx_p.ref_name_to_id)
    g_oracle = np.asarray(CalibrationSubstrate.from_payload(blob_p["payload_gdna"], ra_p)
                          .contained.mass_unspliced, float).sum()

    print(f"var~mean GATING policy comparison (max_outer={MAX_OUTER})")
    print(f"  ZERO+nascent = {ZERO}   (phantom = total assigned gDNA; lower better)")
    print(f"  gDNA-present = {PRES}   (oracle contained gDNA = {g_oracle:,.0f}; assigned closer = better recall)")
    print(f"\n  {'policy':>30} {'zero phantom':>14} {'pres assigned':>14} {'pres under-call':>16} {'ρg_global(z/p)':>16}")
    for name, (gpol, rpol) in POLICIES.items():
        cz = run(blob_z, ra_z, gpol, rpol)
        cp = run(blob_p, ra_p, gpol, rpol)
        ph = np.asarray(cz.mass_gdna_contained, float).sum()
        asg = np.asarray(cp.mass_gdna_contained, float).sum()
        print(f"  {name:>30} {ph:>14,.0f} {asg:>14,.0f} {asg - g_oracle:>+16,.0f} "
              f"{cz.gdna_density_global:>7.4f}/{cp.gdna_density_global:<7.4f}")


if __name__ == "__main__":
    main()
