"""N1 — WHERE DOES THE ENRICHED MASS GO?  An EXACT accounting, not an inference.

The fitted landscape is a weighted sum of per-node kernels:  P = sum_i wt_i * pn_i  (sum wt_i = 1).
So the mass it places above a threshold `t` is EXACTLY

    M(t) = sum_i  wt_i * q_i(t),        q_i(t) = fraction of node i's own kernel above t

which decomposes over any partition of the training nodes.  The oracle landscape is the same object with
uniform weights, true counts, and its own substrate.  Every number below is that identity evaluated — no
model of the deficit, just its ledger.

Three multiplicative factors separate:
    SUBSTRATE   — the share of training WEIGHT sitting on truly-enriched nodes (who is in the set)
    PLACEMENT   — of that weight, the share whose kernel lands above the split (are they put in the right place)
    LEAKAGE     — mass contributed above the split by truly-DEPLETED nodes (false enrichment)

Convention reproduced exactly from gdna_hyperprior_production_plan.md §W3: gDNA-bearing capture-ON/VSTRONG
conditions (n=13), scored at the ORACLE landscape's own split.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_massflow.py
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import gdna_explore_lib as L  # noqa: E402

GRID = L.GRID
_EPS = 1e-12
KNN = 0.5


def conds(suite="ambig"):
    """The 13 gDNA-bearing capture-ON/VSTRONG conditions the enriched-mass headline is measured on."""
    return [s for s in L.load_scenarios(suite)
            if s["group"][1] != "none" and s["group"][0] in ("ON", "VSTRONG")]


def kernels(g, E):
    """Per-node zero-native Poisson posterior on GRID, row-normalised — the recipe's own kernel."""
    lam = np.exp(GRID * L.LN10)[None, :] * E[:, None]
    ll = g[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    return pn / np.maximum(pn.sum(1, keepdims=True), _EPS)


def truth(s):
    """True log10 density per node, at the one-count resolution wall (the metrics' own convention)."""
    return np.log10(np.maximum(s["G"], 1.0)) - np.log10(np.maximum(s["eff"], _EPS))


def observed(s):
    """pass-0's log10 density per node, same convention."""
    return np.log10(np.maximum(s["g_hat"], 1.0)) - np.log10(np.maximum(s["eff"], _EPS))


def ledger(s, sel, w, g, split):
    """The exact M(t) decomposition over {truly-enriched, truly-depleted} training nodes."""
    E = s["eff"][sel]
    pn = kernels(g[sel], E)
    wt = w / max(float(w.sum()), _EPS)
    hi = GRID > split
    q = pn[:, hi].sum(1)                     # per-node share of kernel above the split
    enr = truth(s)[sel] > split
    return dict(
        mass=float((wt * q).sum()),
        w_enr=float(wt[enr].sum()),          # SUBSTRATE: weight share on truly-enriched nodes
        place=float((wt[enr] * q[enr]).sum() / max(wt[enr].sum(), _EPS)),   # PLACEMENT
        from_enr=float((wt[enr] * q[enr]).sum()),
        leak=float((wt[~enr] * q[~enr]).sum()),                             # LEAKAGE
        n_enr=int(enr.sum()), n=int(sel.sum()),
    )


def main():
    ss = conds()
    print("=== N1: the exact mass ledger, 13 gDNA-bearing capture-ON/VSTRONG conditions ===\n")
    print(f"{'condition':44s} {'split':>6s} | {'ORACLE':>27s} | {'FIT (production recipe)':>38s}")
    print(f"{'':44s} {'':>6s} | {'w_enr':>7s} {'place':>6s} {'mass':>6s}   | "
          f"{'w_enr':>7s} {'place':>6s} {'leak':>6s} {'mass':>6s} {'recov':>6s}")
    rows = []
    for s in ss:
        orc = L.oracle_landscape(s)
        split = L.two_component(orc)["split"]
        osel = L.live_region(s)
        o = ledger(s, osel, np.ones(int(osel.sum())), s["G"], split)
        mk = L.masks(s)
        fsel = L.recipe_substrate(s, mk)
        f = ledger(s, fsel, L.recipe_weights(s, fsel, mk), s["g_hat"], split)
        rows.append((o, f))
        print(f"{s['cond'][5:]:44s} {split:+6.2f} | {o['w_enr']:7.4f} {o['place']:6.3f} {o['mass']:6.4f}   | "
              f"{f['w_enr']:7.4f} {f['place']:6.3f} {f['leak']:6.4f} {f['mass']:6.4f} "
              f"{f['mass'] / max(o['mass'], 1e-9):6.2f}")

    def m(rs, k):
        return float(np.mean([r[k] for r in rs]))

    orc_rows = [r[0] for r in rows]
    fit_rows = [r[1] for r in rows]
    print(f"\n{'MEAN':44s} {'':6s} | {m(orc_rows, 'w_enr'):7.4f} {m(orc_rows, 'place'):6.3f} {m(orc_rows, 'mass'):6.4f}   | "
          f"{m(fit_rows, 'w_enr'):7.4f} {m(fit_rows, 'place'):6.3f} {m(fit_rows, 'leak'):6.4f} {m(fit_rows, 'mass'):6.4f}")
    print(f"\n  ORACLE enriched mass {m(orc_rows, 'mass'):.4f} = w_enr {m(orc_rows, 'w_enr'):.4f} x place {m(orc_rows, 'place'):.3f}"
          f" + leak {m(orc_rows, 'leak'):.4f}")
    print(f"  FIT    enriched mass {m(fit_rows, 'mass'):.4f} = w_enr {m(fit_rows, 'w_enr'):.4f} x place {m(fit_rows, 'place'):.3f}"
          f" + leak {m(fit_rows, 'leak'):.4f}")
    print("\n  => the two factors, in ratio to the oracle:")
    print(f"     SUBSTRATE (weight share on truly-enriched nodes): {m(fit_rows, 'w_enr') / m(orc_rows, 'w_enr'):.3f}")
    print(f"     PLACEMENT (that weight's kernel above the split): {m(fit_rows, 'place') / m(orc_rows, 'place'):.3f}")


if __name__ == "__main__":
    main()
