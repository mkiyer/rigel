"""STATE OF CALIBRATION — are we projecting OBSERVED nodes onto the CORRECT MODE?

Every measurement of the hyperprior so far has been about the FIT: does the landscape reproduce the
population's two components (N1's census, N5's shape, W3's projection).  That is the TRAINING side.  The
question never asked directly is the APPLICATION side:

    when the prior is applied at the re-solve, does a node that is TRULY enriched come out enriched?

That is a classification, so score it as one.  Ground truth per live REGION node is its oracle gDNA density
`ρ_g = G/E_g`; "enriched" is above the ORACLE landscape's own split (N1: scoring at the fit's own split
hides the defect).  Assigned is `ρ̂_g = f̂_g·M/E_g` from the solved belief, at pass-0 and at the re-solve, so
the PRIOR's contribution is the difference between the two columns.

Reported mass-weighted, because a node's influence downstream is its fragment mass.

Run: OMP_NUM_THREADS=1 python scratchpad/state_mode_assignment.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402

_EPS = 1e-12
_KNN = 0.5
REFIT = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_refit_cache.pkl")

STRATA = [
    ("ALL 32", lambda s: True),
    ("gDNA-bearing, capture ON/VSTRONG", lambda s: s["group"][1] != "none"
     and s["group"][0] in ("ON", "VSTRONG")),
    ("  ... capture ON only", lambda s: s["group"][1] != "none" and s["group"][0] == "ON"),
    ("  ... capture VSTRONG only", lambda s: s["group"][1] != "none" and s["group"][0] == "VSTRONG"),
    ("gDNA-bearing, capture OFF", lambda s: s["group"][1] != "none" and s["group"][0] == "OFF"),
    ("  ... stranded", lambda s: s["group"][1] != "none" and s["group"][2] == "0.99"),
    ("  ... unstranded", lambda s: s["group"][1] != "none" and s["group"][2] == "0.50"),
    ("ZERO-gDNA (any 'enriched' is FALSE)", lambda s: s["group"][1] == "none"),
]


def joined():
    rf = pickle.loads(REFIT.read_bytes())
    for s in L.load_scenarios("ambig"):
        r = rf.get(s["cond"])
        if r is None:
            continue
        assert np.allclose(r[0]["f0"], s["f0"]), s["cond"]
        yield s, r


def main():
    rows = []
    for s, r in joined():
        live = L.live_region(s) & (s["mass"] > _EPS) & (s["eff"] > 1e-9)
        if not live.any():
            continue
        E = np.maximum(s["eff"][live], _EPS)
        m = s["mass"][live]
        G = np.maximum(s["G"][live], 0.0)
        true_rho = G / E
        # ⚠ A FIXED split at log10 rho_g = 0 (the N5 gamut convention), NOT the per-condition
        # `two_component` valley: on a condition with only ONE true mode (capture OFF, zero-gDNA) the
        # valley is an artefact of a tail and lands anywhere, which made VSTRONG read "0 % truly
        # enriched" — impossible, and it invalidated the first version of this table.
        split = 1.0
        rho0 = np.maximum(r[0]["f0"][live], 0.0) * m / E      # pass-0
        rho1 = np.maximum(r[1]["f0"][live], 0.0) * m / E      # re-solved
        # nodes with NO true gDNA carry no log-density error — they are a FABRICATION question, scored
        # separately, because |log10(rho_hat/0)| is not a number.
        has = G > 0.0
        rows.append(dict(s=s, w=m, has=has, true_hi=(true_rho > split) & has,
                         hi0=rho0 > split, hi1=rho1 > split,
                         e0=np.abs(np.log10(np.maximum(rho0, _EPS) / np.maximum(true_rho, _EPS))),
                         e1=np.abs(np.log10(np.maximum(rho1, _EPS) / np.maximum(true_rho, _EPS))),
                         true_rho=true_rho))

    print("=== M1. MODE ASSIGNMENT — is a truly-ENRICHED node called enriched? (mass-weighted) ===")
    print("    sens = P(called enriched | truly enriched).   FPR = P(called enriched | truly depleted).\n")
    print(f"{'stratum':38s} {'n cond':>6s} {'true enr':>9s} | {'sens p0':>8s} {'sens r1':>8s} | "
          f"{'FPR p0':>7s} {'FPR r1':>7s} | {'FABRIC r1':>9s}")
    for name, sel in STRATA:
        rs = [x for x in rows if sel(x["s"])]
        if not rs:
            continue
        w = np.concatenate([x["w"] for x in rs])
        t = np.concatenate([x["true_hi"] for x in rs])
        h0 = np.concatenate([x["hi0"] for x in rs])
        h1 = np.concatenate([x["hi1"] for x in rs])
        hs = np.concatenate([x["has"] for x in rs])
        wt, wf = w[t].sum(), w[~t & hs].sum()
        sens0 = w[t & h0].sum() / max(wt, _EPS)
        sens1 = w[t & h1].sum() / max(wt, _EPS)
        fpr0 = w[~t & hs & h0].sum() / max(wf, _EPS)
        fpr1 = w[~t & hs & h1].sum() / max(wf, _EPS)
        # FABRICATION: mass with NO true gDNA at all that is nevertheless called enriched
        wz = w[~hs].sum()
        fab1 = w[~hs & h1].sum() / max(wz, _EPS) if wz > 0 else float("nan")
        ok = wt / w.sum() > 0.01
        sfmt = (f"{sens0:8.3f} {sens1:8.3f}" if ok else f"{'--':>8s} {'--':>8s}")
        print(f"{name:38s} {len(rs):6d} {100 * wt / w.sum():8.1f}% | {sfmt} | "
              f"{fpr0:7.4f} {fpr1:7.4f} | {fab1:9.4f}")

    print("\n\n=== M2. DENSITY ACCURACY on the two true populations (mass-wt mean |log10 error|, decades) ===")
    print("    G > 0 nodes only — a node with no true gDNA has no log-density error to score.\n")
    print(f"{'stratum':38s} {'ENRICHED p0':>12s} {'r1':>7s} | {'DEPLETED p0':>12s} {'r1':>7s}")
    for name, sel in STRATA:
        rs = [x for x in rows if sel(x["s"])]
        if not rs:
            continue
        w = np.concatenate([x["w"] for x in rs])
        t = np.concatenate([x["true_hi"] for x in rs])
        e0 = np.concatenate([x["e0"] for x in rs])
        e1 = np.concatenate([x["e1"] for x in rs])
        hs = np.concatenate([x["has"] for x in rs])
        t = t & hs

        def mw(v, m_):
            return np.average(v[m_], weights=w[m_]) if w[m_].sum() > 0 else np.nan
        dep = (~t) & hs
        print(f"{name:38s} {mw(e0, t):12.3f} {mw(e1, t):7.3f} | {mw(e0, dep):12.3f} {mw(e1, dep):7.3f}")

    print("\n\n=== M3. WHERE THE MISASSIGNED ENRICHED MASS GOES (gDNA-bearing capON/VSTRONG) ===")
    rs = [x for x in rows if x["s"]["group"][1] != "none" and x["s"]["group"][0] in ("ON", "VSTRONG")]
    w = np.concatenate([x["w"] for x in rs])
    t = np.concatenate([x["true_hi"] for x in rs])
    h1 = np.concatenate([x["hi1"] for x in rs])
    e1 = np.concatenate([x["e1"] for x in rs])
    miss = t & ~h1
    print(f"   truly-enriched mass          : {w[t].sum():14,.0f}")
    print(f"   ... called DEPLETED at r1    : {w[miss].sum():14,.0f}  ({100 * w[miss].sum() / max(w[t].sum(), _EPS):.1f} %)")
    if miss.any():
        print(f"   their density error (dec)    : {np.average(e1[miss], weights=w[miss]):14.3f}")
        print(f"   the correctly-called ones    : {np.average(e1[t & h1], weights=w[t & h1]):14.3f}")

    print("\n\n=== M4. PER CONDITION — sensitivity, the enriched census the prior is applied to ===")
    print(f"{'condition':50s} {'true enr%':>9s} {'sens p0':>8s} {'sens r1':>8s} {'FPR r1':>8s}")
    for x in sorted(rows, key=lambda z: -np.average(z["true_hi"], weights=z["w"])):
        w, t, h0, h1 = x["w"], x["true_hi"], x["hi0"], x["hi1"]
        if w[t].sum() <= 0:
            continue
        print(f"{x['s']['cond'][5:]:50s} {100 * w[t].sum() / w.sum():8.1f}% "
              f"{w[t & h0].sum() / w[t].sum():8.3f} {w[t & h1].sum() / w[t].sum():8.3f} "
              f"{w[~t & h1].sum() / max(w[~t].sum(), _EPS):8.4f}")


if __name__ == "__main__":
    main()
