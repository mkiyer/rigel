"""N1 — WHERE the mass lands, and WHO the split-crossers are.

Part A — DESTINATION.  Mass is conserved, so every unit of training weight lands somewhere on the grid.
This partitions the training set into populations and reports, for each, its share of the total weight and
the share of its own kernel that lands in each density band.  The band edges are the two ground-truth levels
(the depleted point mass and the enriched mode) with the oracle's valley between them, so "the valley" is a
real place, not an artefact of binning.

Part B — THE SPLIT-CROSSERS.  Truly-enriched nodes inside the training substrate that pass-0 places on the
depleted side of the split (HANDOFF_17 §0: 217/875 unstranded, 9/752 stranded).  Profiled by node class,
eff-length, mass, count, declared variance and reliability weight.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_where.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds, kernels, observed, truth  # noqa: E402

_EPS = 1e-12
NTYPE = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def populations(s, mk, sel, tr, split):
    """The training set partitioned by what a node IS (truth + structure), not by what pass-0 said."""
    enr = tr > split
    return [
        ("region x truly-ENRICHED", sel & mk["region"] & enr),
        ("region x truly-depleted", sel & mk["region"] & ~enr),
        ("boundary x truly-ENRICHED", sel & mk["boundary"] & enr),
        ("boundary x truly-depleted", sel & mk["boundary"] & ~enr),
    ]


def part_a():
    ss = conds()
    # bands: [-inf, dep+0.5] depleted level | (.., split] the VALLEY | (split, inf) enriched
    tot = {}
    for s in ss:
        split = L.two_component(L.oracle_landscape(s))["split"]
        mk = L.masks(s)
        sel = L.recipe_substrate(s, mk)
        tr = truth(s)
        w_all = np.zeros(sel.shape)
        w_all[sel] = L.recipe_weights(s, sel, mk)
        w_all /= max(w_all.sum(), _EPS)
        dep_lvl = float(np.median(tr[mk["region"] & (s["ntype"] < 2) & (s["eff"] > 1e-9)]))
        edges = (dep_lvl + 0.5, split)
        for name, m in populations(s, mk, sel, tr, split):
            if not m.any():
                continue
            pn = kernels(s["g_hat"][m], s["eff"][m])
            w = w_all[m]
            b_low = float((w * pn[:, L.GRID <= edges[0]].sum(1)).sum())
            b_val = float((w * pn[:, (L.GRID > edges[0]) & (L.GRID <= edges[1])].sum(1)).sum())
            b_hi = float((w * pn[:, L.GRID > edges[1]].sum(1)).sum())
            a = tot.setdefault(name, np.zeros(5))
            a += np.array([w.sum(), b_low, b_val, b_hi, m.sum()])
    n = len(ss)
    print("=== A. DESTINATION of the training weight (13 conditions, production substrate + weights) ===")
    print("    bands: DEPLETED = within 0.5 dec of the true gDNA level | VALLEY = between it and the "
          "oracle's split | ENRICHED = above the split\n")
    print(f"{'population':28s} {'n':>7s} {'weight':>8s} | {'-> DEPLETED':>12s} {'-> VALLEY':>10s} "
          f"{'-> ENRICHED':>12s}")
    grand = np.sum([v for v in tot.values()], 0)
    for name, a in tot.items():
        print(f"{name:28s} {a[4] / n:7.0f} {a[0] / n:8.4f} | {a[1] / n:12.4f} {a[2] / n:10.4f} "
              f"{a[3] / n:12.4f}")
    print(f"{'TOTAL':28s} {grand[4] / n:7.0f} {grand[0] / n:8.4f} | {grand[1] / n:12.4f} "
          f"{grand[2] / n:10.4f} {grand[3] / n:12.4f}")
    print(f"\n  fitted enriched mass = {grand[3] / n:.4f}   (oracle, region-only: 0.1405)")


def part_b():
    ss = conds()
    print("\n=== B. THE SPLIT-CROSSERS: truly-enriched training nodes that pass-0 puts below the split ===")
    for lbl, keep in (("unstranded (ss 0.50)", lambda s: s["group"][2] == "0.50"),
                      ("stranded  (ss 0.99)", lambda s: s["group"][2] == "0.99")):
        rec = {}
        for s in ss:
            if not keep(s):
                continue
            split = L.two_component(L.oracle_landscape(s))["split"]
            mk = L.masks(s)
            sel = L.recipe_substrate(s, mk)
            tr, ob = truth(s), observed(s)
            w = np.zeros(sel.shape)
            w[sel] = L.recipe_weights(s, sel, mk)
            enr = sel & (tr > split)
            cross = enr & (ob <= split)
            for k, m in (("all truly-enriched", enr), ("CROSSERS", cross)):
                d = rec.setdefault(k, {c: [] for c in
                                       ("n", "shift", "eff", "mass", "G", "g_hat", "var", "w", "ntype")})
                d["n"].append(int(m.sum()))
                d["shift"].extend((ob - tr)[m])
                d["eff"].extend(s["eff"][m])
                d["mass"].extend(s["mass"][m])
                d["G"].extend(s["G"][m])
                d["g_hat"].extend(s["g_hat"][m])
                d["var"].extend(np.nan_to_num(s["var"][m], nan=1e9, posinf=1e9))
                d["w"].extend(w[m])
                d["ntype"].extend(s["ntype"][m])
        n_enr = sum(rec["all truly-enriched"]["n"])
        n_cr = sum(rec["CROSSERS"]["n"])
        print(f"\n  {lbl}: {n_cr}/{n_enr} cross ({100 * n_cr / max(n_enr, 1):.1f} %)")
        print(f"    {'group':20s} {'n':>6s} {'med shift':>10s} {'med eff':>8s} {'med mass':>9s} "
              f"{'med G':>8s} {'med w':>7s} {'med var':>9s}")
        for k in ("all truly-enriched", "CROSSERS"):
            d = rec[k]
            print(f"    {k:20s} {sum(d['n']):6d} {np.median(d['shift']):+10.3f} "
                  f"{np.median(d['eff']):8.1f} {np.median(d['mass']):9.1f} {np.median(d['G']):8.1f} "
                  f"{np.median(d['w']):7.3f} {np.median(d['var']):9.4f}")
        for k in ("all truly-enriched", "CROSSERS"):
            nt = np.array(rec[k]["ntype"])
            comp = "  ".join(f"{NTYPE[t]}={int((nt == t).sum())}" for t in sorted(set(nt.tolist())))
            print(f"    {k:20s} class mix: {comp}")
        # does the crossing actually cost the fit? their weight share of the enriched population
        we, wc = np.array(rec["all truly-enriched"]["w"]), np.array(rec["CROSSERS"]["w"])
        print(f"    crossers hold {wc.sum() / max(we.sum(), _EPS):.1%} of the truly-enriched weight "
              f"({n_cr / max(n_enr, 1):.1%} of the nodes)")


if __name__ == "__main__":
    part_a()
    part_b()
