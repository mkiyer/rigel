"""Phase-2 gate: validate the PRODUCTION `bp_solver.fit_enrichment_transfer` on the flagship (capture on/off).

Confirms the production function (z channel + pre-loop ê fit + conditional σ²_bio + significance gate) matches
the validated harness:
  - capture ON  → significant; ê recovers the withheld AMBIG exons (net within a few %, corr ~0.6); enriched
    σ²_bio gives N≈24.
  - capture OFF → NOT significant ⇒ ê ≡ ρ_global (flat) ⇒ the z-prior equals the genome-wide global (invariance).
  - ê is fit on the strand-only init f_g (static) — deterministic, no belief feedback.

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_phase2_check.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
ON = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
OFF = "gdna_gdna300_ss_0.99_nrna_none_capture_off"


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs = CalibrationSubstrate.from_payload(payload, ra)
    bs = BoundarySubstrate.from_payload(payload)
    geom = bp.build_node_geometry(ch, cs, bs, ra, blob["gdna_pmf"], blob["rna_pmf"])
    stat = bp.build_node_statics(ch, cs, bs, ra)
    cap = {}
    orig = bp.node_sweep
    calibrate.__globals__["node_sweep"] = lambda c, s, g, b, rg, bsub, **k: (cap.update(k), orig(c, s, g, b, rg, bsub, **k))[1]
    calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    calibrate.__globals__["node_sweep"] = orig
    kappa = cap["rna_sense_frac"]
    odg, odr = cap.get("gdna_strand_overdispersion", 0.0), cap.get("rna_strand_overdispersion", 0.0)

    is_reg = np.asarray(ch.kind) == REGION
    refidx = np.asarray(ch.ref_idx, np.int64)
    EGl = np.asarray(geom.eff_gdna_left)
    Ml = np.asarray(geom.mass_left)
    T = Ml / np.maximum(EGl, _EPS)
    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    csr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    gR = np.asarray(csg.contained.mass_unspliced, float)
    rR = np.asarray(csr.contained.mass_unspliced, float)
    r = np.clip(refidx, 0, gR.shape[0] - 1)
    ofg = gR[r] / np.maximum(gR[r] + rR[r], _EPS)
    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])[r]
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])[r]

    # strand-only init f_g (the basis the production _gdna_seed_estimate / fit_enrichment_transfer use)
    lat = _simplex_lattice(60)
    fgg = lat[2]
    f_init = np.zeros(ch.n_nodes)
    for i in np.where(is_reg & ((scl == 1) | (scl == 2)) & (Ml > 0))[0]:
        psi = _local_loglik(stat.u_pos[i:i + 1], stat.u_neg[i:i + 1], stat.spliced_pos[i:i + 1],
                            stat.spliced_neg[i:i + 1], stat.free_pos[i:i + 1], stat.free_neg[i:i + 1],
                            kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i + 1])
        f_init[i] = float(_fg_median(psi, fgg)[0])

    rho_global, _gvm, _vm = bp._gdna_seed_estimate(ch, stat, geom, ra, bs, f_init, kappa)
    ehat, z, sig2, significant = bp.fit_enrichment_transfer(ch, stat, geom, ra, bs, f_init, kappa, rho_global)

    # AMBIG-exon recovery via the production ê
    amb = np.where(is_reg & (scl == 3) & (tcl == 2) & (Ml > 0) & np.isfinite(z))[0]
    trueg = gR[r[amb]]
    pred = ehat.predict(z[amb])
    fg = np.clip(pred / np.maximum(T[amb], _EPS), 0, 1)
    net = (fg * Ml[amb]).sum() - trueg.sum()
    corr = float(np.corrcoef(fg, ofg[amb])[0, 1])
    zf = z[np.isfinite(z) & is_reg]
    zlo, zhi = float(np.nanmin(zf)), float(np.nanmax(zf))
    print(f"\n=== {cond} ===")
    print(f"  ρ_global={rho_global:.3f}   significant={significant}")
    print(f"  ê over fit-z range [{zlo:.3f},{zhi:.3f}]: ê(lo)={float(ehat.predict(np.array([zlo]))[0]):.3f}  "
          f"ê(hi)={float(ehat.predict(np.array([zhi]))[0]):.3f}")
    print(f"  AMBIG-exon recovery (true {trueg.sum():,.0f}): pred {(fg*Ml[amb]).sum():,.0f}  "
          f"net {net:+,.0f} ({100*net/max(trueg.sum(),1):+.1f}%)  corr {corr:.2f}")
    if significant and sig2 is not None:
        # enriched-level precision N = ê²/(σ²_bio + sampling), sampled at a high-z enriched node
        hz = float(np.nanpercentile(z[amb], 75))
        e_hi = float(ehat.predict(np.array([hz]))[0])
        s2 = float(sig2.predict(np.array([e_hi]))[0])
        ncap = e_hi ** 2 / max(s2, _EPS)
        print(f"  enriched precision @z≈{hz:.2f} (ê≈{e_hi:.1f}): σ²_bio={s2:.3f}  N≈{ncap:.1f}")
    return cond, significant, rho_global, net, trueg.sum(), corr, ehat


def main():
    print("PHASE 2 gate — production bp_solver.fit_enrichment_transfer on the flagship")
    on = run(ON)
    off = run(OFF)
    print("\n=== GATES ===")
    g_on_sig = on[1] is True
    g_on_rec = abs(on[3]) / max(on[4], 1) <= 0.08 and on[5] >= 0.55
    g_off_flat = (off[1] is False) and np.allclose(off[6].predict(np.array([0.01, 1.0, 100.0])), off[2])
    print(f"  capture-ON significant + recovers AMBIG:  {'PASS' if g_on_sig and g_on_rec else 'FAIL'} "
          f"(sig={on[1]}, net {100*on[3]/max(on[4],1):+.1f}%, corr {on[5]:.2f})")
    print(f"  capture-OFF collapses to flat ρ_global:   {'PASS' if g_off_flat else 'FAIL'} "
          f"(sig={off[1]}, ê≡ρ_global={off[2]:.3f})")
    print(f"\n  OVERALL: {'PASS — proceed to Phase 3' if g_on_sig and g_on_rec and g_off_flat else 'gate failure'}")


if __name__ == "__main__":
    main()
