"""RESURRECTION MEASUREMENT — the gDNA-hyperprior (Role B) track, re-run at the rebuilt pass-0.

The track was paused 2026-07-21 (`archive/enrichment_sensitivity_worklog.md` §6) on the owner's call: the
landscape's enriched mode was under-massed because PASS-0 itself under-called enriched gDNA, so fixing the
projection was treating a symptom. Pass-0 has since been rebuilt (suite mwae 0.15 -> 0.084). This re-measures
both goals against the numbers recorded at the pause.

  GOAL 1  the LANDSCAPE  — EMD-to-oracle of P(log10 rho_g).            anchor: 0.267 cross-suite
  GOAL 2  the PROJECTION — does an enriched node's observed density get PULLED UP to the enriched level?
          enr_recovery = mean(mu* - obs) over truly-enriched nodes     anchor: -0.05 symmetric / +0.25 asym
          enr_abs_err  = mean|mu* - oracle|                            anchor: 0.226 symmetric / ~0.20 asym
          fabrication  = mu* on zero-gDNA capOFF exons (specificity)   anchor: -0.83 -> -0.4..-0.6

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_resurrect.py [--plot]
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import gdna_explore_lib as L  # noqa: E402
from gdna_projection import project_asym  # noqa: E402

SUITES = ("ambig", "quick")


def load_recipe():
    """The 2026-07-21 'unified' landscape recipe.

    Kept as a function so the several callers do not all have to change: it now simply returns
    `gdna_explore_lib.recipe`, which was promoted out of an `exec()`'d bare-function file on 2026-07-27
    (verified byte-identical over all 48 cached scenarios)."""
    return L.recipe


def goal1(recipe):
    print("=== GOAL 1 — the LANDSCAPE (EMD to the oracle landscape; anchor 0.267 cross-suite) ===")
    tot = []
    for su in SUITES:
        try:
            scen = L.load_scenarios(su)
        except FileNotFoundError:
            continue
        r = L.evaluate(recipe, scen)
        tot.append(r["mean_emd"])
        print(f"  {su:6s}  mean_emd={r['mean_emd']:.3f}  max_emd={r['max_emd']:.3f}  "
              f"mean_enr_l1={r['mean_enr']:.3f}  max_enr_l1={r['max_enr']:.3f}")
    if tot:
        print(f"  {'CROSS':6s}  mean_emd={np.mean(tot):.3f}   (was 0.267)")


def _score(recipe, scen, proj):
    """enr_recovery / enr_abs_err over truly-enriched nodes, + the zero-gDNA fabrication canary."""
    rec, err, fab = [], [], []
    for s in scen:
        P = recipe(s)
        if P is None:
            continue
        live = s["eff"] > 1e-9
        obs = np.log10(np.maximum(s["g_hat"], 1e-9) / np.maximum(s["eff"], 1e-9))
        tru = np.log10(np.maximum(s["G"], 1e-9) / np.maximum(s["eff"], 1e-9))
        enr = live & (s["G"] > 0) & (tru > 0.0)
        if enr.sum() >= 5:
            mu = proj(P, obs[enr])
            rec.append(float((mu - obs[enr]).mean()))
            err.append(float(np.abs(mu - tru[enr]).mean()))
        if s["group"][0] == "OFF" and s["group"][1] == "none":   # zero-gDNA, capture-off: must NOT fabricate
            m = live & (s["mass"] > 1e-12) & (s["ntype"] == 2)
            if m.sum() >= 5:
                fab.append(float(proj(P, obs[m]).mean()))
    return (float(np.mean(rec)) if rec else np.nan, float(np.mean(err)) if err else np.nan,
            float(np.mean(fab)) if fab else np.nan)


def goal2(recipe):
    print("\n=== GOAL 2 — the PROJECTION (the endpoint: is an enriched node PULLED UP?) ===")
    print(f"  {'projection':22s} {'suite':6s} {'enr_recovery':>13s} {'enr_abs_err':>12s} {'fabrication':>12s}")
    projs = [("symmetric h=0.15", lambda P, d: L.project(P, d, 0.15, "mean")),
             ("symmetric h=0.40", lambda P, d: L.project(P, d, 0.40, "mean")),
             ("symmetric mode-readout", lambda P, d: L.project(P, d, 0.15, "mode")),
             ("asymmetric (Band-Aid)", project_asym)]
    for name, proj in projs:
        for su in SUITES:
            try:
                scen = L.load_scenarios(su)
            except FileNotFoundError:
                continue
            r, e, f = _score(recipe, scen, proj)
            print(f"  {name:22s} {su:6s} {r:+13.3f} {e:12.3f} {f:+12.3f}")
    print("  anchors 2026-07-21: symmetric -0.050 / 0.226 · asymmetric +0.250 / ~0.200 · fab -0.4..-0.6")


def plot(recipe, out="/Users/mkiyer/proj/rigel/docs/calibration/figures/resurrect_fit_vs_oracle.png"):
    """The fit-vs-oracle overlay the owner was reading when the track was paused."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scen = L.load_scenarios("ambig")
    want = [s for s in scen if s["group"][2] == "0.50"]           # unstranded — where the prior must work
    want = sorted(want, key=lambda s: (s["group"][0], s["group"][1]))[:12]
    fig, axes = plt.subplots(3, 4, figsize=(20, 11), sharex=True)
    for ax, s in zip(axes.ravel(), want):
        orc, fit = L.oracle_landscape(s), recipe(s)
        ax.fill_between(L.GRID, 0, orc, color="0.75", label="ORACLE")
        if fit is not None:
            ax.plot(L.GRID, fit, "C3", lw=1.6, label="fit (pass-0)")
        ax.set_title(f"{s['group'][0]} {s['group'][1]} {s['group'][3]}\nEMD={L.emd(fit, orc):.3f}", fontsize=9)
        ax.set_xlim(-5, 2.0)
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle("gDNA hyperprior landscape — fit vs ORACLE at the rebuilt pass-0 (unstranded conditions)")
    fig.supxlabel("log10 rho_g (gDNA density)")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    rec = load_recipe()
    goal1(rec)
    goal2(rec)
    if "--plot" in sys.argv:
        plot(rec)
