"""N1 — the ARMS implied by the ledger, one change each.

The ledger (`gdna_n1_massflow` / `_ladder` / `_factorial`) closes the enriched-mass deficit exactly:

    recovery 0.399 = SUBSTRATE 0.574 x COUNTS 1.014 x WEIGHT 0.687
    SUBSTRATE 0.574 = AMBIG-excluded 0.802 x BOUNDARIES-included 0.682

so there are exactly two live levers (counts are inert, AMBIG is a standing structural decision).  This
measures each on its own, plus the one the mechanism actually suggests:

  * **BND-OUT** — train on regions only, which is what PRODUCTION `_fit_gdna_hyperprior` already does;
    boundaries were added to the exploration recipe provisionally on 2026-07-27.
  * **FLAT** — the weight ablation, for reference (not a candidate: it discards real precision).
  * **WIDTH** — precision moved from the MASS to the KERNEL WIDTH.  A node's declared variance currently
    *deletes* it from the population (`w = ref/(v+ref)`), and never widens its kernel.  For a MIXTURE, that
    is the wrong place for it: precision-weighting is correct for estimating a common mean, but the mixing
    proportions of a population are a *census*, and down-weighting a node removes it from the census.
    Measured mechanism: on unstranded data `w` does not separate enriched from depleted *within* a class
    (exon 0.466 vs 0.459) — it separates CLASSES (exon 0.46 vs intron 0.95 vs intergenic 1.00), and every
    enriched node is an exon.  So the weight halves exactly the class that carries the enriched mode.
    The honest statement of an imprecise node is a BROAD kernel of unit mass: Poisson (x) N(0, tau^2).

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_arms.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import kernels  # noqa: E402

_EPS = 1e-12
KNN = 0.5


def tau_dec(s):
    """The node's declared width OF ITS LOG DENSITY, in decades. The density is ``f_g * M / E`` with `M`, `E`
    known, so its log-width is the belief's own ``sqrt(Var(log f_g))/ln10`` — already computed by the solver,
    no constant introduced.

    ⚠ Exactly on the zero-mass nodes (`M = 0`) that width is **0, not infinite**: ``f_g * 0 = 0`` for every
    `f_g`, so the density is known exactly even though the composition is not identified at all. Measured:
    all 928 non-finite `var_gdna` values in the substrate are precisely those nodes. This is the value of the
    quantity, not an admission gate — and it is the same structural predicate `recipe_weights` already uses
    to give the zero anchor `w = 1`."""
    v = np.nan_to_num(s["var"], nan=0.0, posinf=0.0)
    return np.where(s["mass"] > _EPS, np.sqrt(np.maximum(v, 0.0)) / L.LN10, 0.0)


def landscape(s, sel, w, *, tau=None, knn_scale=KNN, n_bins=12):
    """One landscape. `w` is per-node MASS; `tau` (decades, optional) is per-node extra kernel WIDTH.
    The kNN population resolution is applied exactly as `gdna_explore_lib.recipe` does, binning on the
    node's TOTAL width so the two width sources compose."""
    if not sel.any():
        return None
    g, E = s["g_hat"][sel], s["eff"][sel]
    pn = kernels(g, E)
    a = np.log10(np.maximum(g, 1.0)) - np.log10(np.maximum(E, _EPS))
    n = a.size
    if knn_scale > 0:
        kk = max(int(round(np.sqrt(n))), 2)
        srt = np.sort(a)
        idx = np.searchsorted(srt, a)
        lo, hi = np.maximum(idx - kk, 0), np.minimum(idx + kk, n - 1)
        h_i = knn_scale * np.maximum(np.maximum(a - srt[lo], srt[hi] - a), L.GRID_H)
    else:
        h_i = np.full(n, L.GRID_H)
    if tau is not None:
        h_i = np.sqrt(h_i**2 + np.maximum(tau[sel], 0.0) ** 2)
    edges = np.quantile(h_i, np.linspace(0.0, 1.0, n_bins + 1))
    out = np.zeros_like(L.GRID)
    for b in range(n_bins):
        m = (h_i >= edges[b]) & (h_i <= edges[b + 1] if b == n_bins - 1 else h_i < edges[b + 1])
        if m.any():
            out += L._convolve((w[m][:, None] * pn[m]).sum(0), float(np.mean(h_i[m])))
    t = float(out.sum())
    return out / t if t > 0 and np.isfinite(t) else None


def arms(s):
    """(name, selection, mass-weights, extra-width) — one change per arm from the production recipe."""
    mk = L.masks(s)
    live = s["eff"] > 1e-9
    reg = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
    bnd = mk["boundary"] & live & mk["havemass"] & (mk["single"] | mk["gonly"])
    prod = reg | bnd
    t = tau_dec(s)
    yield "BASE (production recipe)", prod, L.recipe_weights(s, prod, mk), None
    yield "BND-OUT (regions only)", reg, L.recipe_weights(s, reg, mk), None
    yield "FLAT (no precision)", prod, np.ones(int(prod.sum())), None
    yield "WIDTH (precision -> kernel)", prod, np.ones(int(prod.sum())), t
    yield "WIDTH + BND-OUT", reg, np.ones(int(reg.sum())), t
    # CONTROL: is the per-node structure of tau doing the work, or is this just global smoothing under
    # another name (which W2 refuted)? Same arm with tau replaced by its own median over the substrate.
    med = np.full_like(t, float(np.median(t[reg])) if reg.any() else 0.0)
    yield "  ctrl: CONST-tau + BND-OUT", reg, np.ones(int(reg.sum())), med


def run(ss, title):
    acc = {}
    for s in ss:
        orc = L.oracle_landscape(s)                       # REGION-only oracle: the production target
        split = L.two_component(orc)["split"]
        hi = L.GRID > split
        for name, sel, w, tau in arms(s):
            P = landscape(s, sel, w, tau=tau)
            if P is None:
                continue
            tc = L.two_component(P, split=split)
            own = L.two_component(P)
            matched = L.oracle_landscape(s, sel=sel)
            acc.setdefault(name, []).append(dict(
                mass=float(P[hi].sum()), emd=L.emd(P, orc), emd_m=L.emd(P, matched),
                dep_loc=tc["dep_loc"], dep_w=tc["dep_width"], enr_loc=tc["enr_loc"],
                enr_w=tc["enr_width"], own_split=own["split"],
            ))
    o = [L.two_component(L.oracle_landscape(s)) for s in ss]
    om = float(np.mean([x["enr_mass"] for x in o]))
    print(f"=== {title} (n={len(ss)}), kNN 0.5 kernel throughout ===")
    print(f"    ORACLE (region-only): enriched mass {om:.4f}, dep {np.mean([x['dep_loc'] for x in o]):+.3f}"
          f"/{np.mean([x['dep_width'] for x in o]):.3f}, enr "
          f"{np.mean([x['enr_loc'] for x in o]):+.3f}/{np.mean([x['enr_width'] for x in o]):.3f}, "
          f"split {np.mean([x['split'] for x in o]):+.2f}\n")
    print(f"{'arm':28s} {'enr mass':>9s} {'recov':>6s} {'EMD/reg':>8s} {'EMD/match':>10s} "
          f"{'dep loc':>8s} {'dep wid':>8s} {'enr loc':>8s} {'enr wid':>8s} {'own split':>10s}")
    for name, rs in acc.items():
        def m(k):
            return float(np.nanmean([r[k] for r in rs]))
        print(f"{name:28s} {m('mass'):9.4f} {m('mass') / om:6.2f} {m('emd'):8.4f} {m('emd_m'):10.4f} "
              f"{m('dep_loc'):+8.3f} {m('dep_w'):8.3f} {m('enr_loc'):+8.3f} {m('enr_w'):8.3f} "
              f"{m('own_split'):+10.2f}")
    print()


def main():
    for su in ("ambig", "quick"):
        allsc = L.load_scenarios(su)
        run([s for s in allsc if s["group"][1] != "none" and s["group"][0] in ("ON", "VSTRONG")],
            f"{su}: gDNA-bearing capture-ON/VSTRONG — the enriched-mass headline set")
        run(allsc, f"{su}: ALL conditions — the generalization / no-regression check")
        run([s for s in allsc if s["group"][1] == "none"],
            f"{su}: ZERO-gDNA conditions — the fabrication guard (enr mass must stay ~0)")
    print("  EMD/reg   = against the region-only oracle (production's target — production trains on regions)")
    print("  EMD/match = against an oracle on the arm's OWN nodes (removes the substrate from the score)")


if __name__ == "__main__":
    main()
