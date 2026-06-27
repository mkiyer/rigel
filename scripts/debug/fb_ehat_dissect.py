"""Dissect the enrichment transfer ê(z) node-by-node — the dominant calibration error (32–46%).

ê(z) = E[ρ_g | z] is fit on SINGLE-STRAND EXON 'teachers' (where gDNA density is observable) and APPLIED to
AMBIG EXONS as the prior mean μ = clip(ê(z)·E/M, 0, 1). The dissection showed it over-calls AMBIG exons
(injects phantom). This goes node-by-node to find WHY:

  TEACHER side — for the single-strand exons the fit learned from: is the (z, response ρ_g) it fit biased vs the
                 oracle gDNA density? (a biased teacher ⇒ a biased ê.)
  APPLY side   — for AMBIG exons: z, ρ̂=ê(z), implied μ, vs the ORACLE gDNA density/fraction. Is z even
                 predictive of gDNA, or is it contaminated by RNA crossing (so ê maps RNA→gDNA)?

    OMP_NUM_THREADS=1 python scripts/debug/fb_ehat_dissect.py [cond ...]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.node_chain import build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = ["gdna_gdna300_ss_0.50_nrna_none_capture_on", "gdna_none_ss_0.50_nrna_rnd_capture_off"]
_EPS = 1e-9


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    finally:
        calmod.node_sweep = orig

    ehat, z, apply, sig = cap["ehat"], cap["z_enrich"], cap["enrich_apply"], cap["enrich_sig"]
    effg, massg = cap["eff_global"], cap["mass_global"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    refidx = np.asarray(ch.ref_idx, np.int64)

    # oracle per region: gDNA / RNA contained (split BAMs); densities use the gDNA contained eff-len.
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    r_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    reg_eff_g = region_eff_length(np.asarray(ra.region_size_bp, float), blob["gdna_pmf"])
    R = g_or.shape[0]
    fg_true = g_or / np.maximum(g_or + r_or, _EPS)
    rho_g_true = g_or / np.maximum(reg_eff_g, _EPS)            # oracle gDNA crossing-scale density
    rho_tot_true = (g_or + r_or) / np.maximum(reg_eff_g, _EPS)  # oracle TOTAL density (what z ~ tracks)

    print(f"\n=== {cond} : ê(z) dissection  (significant={sig}) ===")
    if not sig:
        print("  ê collapsed to flat ρ_global (no z→ρ_g signal). Nothing to dissect.")
        return

    rows = []
    for i in np.where(apply)[0]:
        rg = int(refidx[i])
        if rg >= R or (g_or[rg] + r_or[rg]) <= _EPS:
            continue
        zi = float(z[i])
        rho_hat = float(max(ehat.predict(np.array([zi]))[0], 0.0))
        mu = float(np.clip(rho_hat * effg[i] / max(massg[i], _EPS), 0.0, 1.0))
        rows.append(dict(reg=rg, z=zi, rho_hat=rho_hat, mu=mu, fg_true=fg_true[rg],
                         rho_g_true=rho_g_true[rg], rho_tot_true=rho_tot_true[rg],
                         mass=float(g_or[rg] + r_or[rg]), gdna_over=(mu - fg_true[rg]) * (g_or[rg] + r_or[rg])))
    if not rows:
        print("  (no AMBIG-exon apply nodes)")
        return

    over = sum(r["gdna_over"] for r in rows)
    print(f"  APPLY set: {len(rows)} AMBIG exons; ê(z) implied-μ gDNA-mass error vs true = {over:>+12,.0f}")
    # is z predictive of gDNA, or of TOTAL (RNA-contaminated)? correlate z with oracle gDNA vs total density.
    za = np.array([r["z"] for r in rows])
    rgt = np.array([r["rho_g_true"] for r in rows])
    rtt = np.array([r["rho_tot_true"] for r in rows])

    def corr(a, b):
        if a.std() < _EPS or b.std() < _EPS:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])
    print(f"  corr(z, oracle gDNA density) = {corr(za, rgt):+.2f}   "
          f"corr(z, oracle TOTAL density) = {corr(za, rtt):+.2f}   "
          f"(if z tracks TOTAL not gDNA, ê maps RNA→gDNA)")
    print(f"  {'reg':>5} {'mass':>8} {'z':>8} {'ρ̂=ê(z)':>8} {'μ(impl)':>7} {'fg_true':>7} "
          f"{'ρg_true':>8} {'ρtot_tru':>8}  note")
    for r in sorted(rows, key=lambda r: -abs(r["gdna_over"]))[:15]:
        note = "ê OVER-calls" if r["mu"] - r["fg_true"] > 0.1 else ("ê under" if r["mu"] - r["fg_true"] < -0.1 else "·")
        print(f"  {r['reg']:>5} {r['mass']:>8,.0f} {r['z']:>8.2f} {r['rho_hat']:>8.2f} {r['mu']:>7.2f} "
              f"{r['fg_true']:>7.2f} {r['rho_g_true']:>8.2f} {r['rho_tot_true']:>8.2f}  {note}")


def main():
    for c in (sys.argv[1:] or DEFAULT):
        run(c)


if __name__ == "__main__":
    main()
