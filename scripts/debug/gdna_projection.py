"""The DNA-prior PROJECTION (the endpoint): observed DNA log10-density -> anchor mu* onto the fitted landscape.

THE ELEGANT-CORE RESULT (2026-07-21): clean unified landscape (`gdna_explore_lib.recipe`) + ASYMMETRIC-UPWARD
projection. Pass-0 systematically UNDER-calls enriched gDNA, so a node's observed density is a LOWER bracket:
trust landscape mass ABOVE it (wide hup), refuse mass BELOW it (tight hdn->0), bounded to cap_up decades. This
recovers the under-called enriched nodes (enr_recovery -0.05 -> +0.25, err improves 0.226 -> ~0.20 cross-suite),
MATCHING a 15-constant workflow synthesis at just 3 projection constants / 0 special cases.

The last remaining magic is hup/cap (the de-bias magnitude) — the open follow-up is to DERIVE them per-node from
the deconvolution/belief uncertainty (an uncertain node can be under-called more -> wider hup). hdn->0 is fully
principled (the observation is a floor). Run: OMP_NUM_THREADS=1 python scripts/debug/gdna_projection.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import gdna_explore_lib as L  # noqa: E402

GRID = L.GRID


def project_asym(P, d_obs, hup=0.70, hdn=0.02, cap_up=1.0):
    """Asymmetric-upward sampling-likelihood projection. d_obs is a LOWER bracket (pass-0 under-calls enriched
    gDNA); look UP (wide hup), not DOWN (tight hdn), bounded to cap_up decades. Responsibility-mean read-out."""
    d = np.asarray(d_obs, float)
    logP = np.log(np.maximum(np.asarray(P, float), 1e-30))
    diff = GRID[None, :] - d[:, None]  # mu - obs
    h = np.where(diff >= 0, hup, hdn)
    lk = np.where(diff > cap_up, -1e9, -0.5 * (diff / h) ** 2)  # bounded upward de-bias
    lr = logP[None, :] + lk
    lr -= lr.max(1, keepdims=True)
    r = np.exp(lr)
    r /= np.maximum(r.sum(1, keepdims=True), 1e-30)
    return (r * GRID[None, :]).sum(1)


def score_enriched(suite, recipe, proj=project_asym):
    """Enriched-sensitivity of `recipe`+`proj` on a suite. Returns (enr_recovery, enr_abs_err, fabrication).
    enr_recovery = mean(mu*-obs) over truly-enriched nodes (want POSITIVE); enr_abs_err = mean|mu*-oracle|;
    fabrication = mean(mu*) on zero-DNA capOFF exons (want low)."""
    rec, err, fab = [], [], []
    for s in L.load_scenarios(suite):
        P = recipe(s)
        if P is None:
            continue
        g, E, G = s["g_hat"], s["eff"], s["G"]
        live = E > 1e-9
        obs = np.log10(np.maximum(g, 1e-9) / np.maximum(E, 1e-9))
        tru = np.log10(np.maximum(G, 1e-9) / np.maximum(E, 1e-9))
        enr = live & (G > 0) & (tru > 0)
        if enr.sum() >= 5:
            mu = proj(P, obs[enr])
            rec.append(float((mu - obs[enr]).mean()))
            err.append(float(np.abs(mu - tru[enr]).mean()))
        if s["group"][0] == "OFF" and s["group"][1] == "none":
            m = live & (s["mass"] > 1e-12) & (s["ntype"] == 2)
            if m.sum() >= 5:
                fab.append(float(proj(P, obs[m]).mean()))
    return (float(np.mean(rec)), float(np.mean(err)),
            float(np.mean(fab)) if fab else float("nan"))


if __name__ == "__main__":
    for su in ("ambig", "quick"):
        r, e, f = score_enriched(su, L.recipe)
        print(f"{su:6s}  enr_recovery={r:+.3f}  enr_abs_err={e:.3f}  fabrication={f:+.3f}")
