"""DISSECTION 6 — THE SOURCE: the mature-RNA GRAFT under-claims, and the shortfall becomes gDNA.

`gdna_d5_zero_source.py` localized it: on a ZERO-gDNA library, 84 % of the false-positive gDNA sits on
exons that DID receive an RNA imputation message, and that message claims **f_rna = 0.646** where the truth
is exactly **1.000** — flat across every precision quintile (0.635 / 0.536 / 0.658 / 0.661 / 0.662), so it
is a systematic MODE bias, not scatter. Whatever the message does not claim is left to gDNA by construction,
and the solve lands at f_g = 0.271 just under the message's own 0.354 floor.

The candidate mechanism is already on the books. The graft's claim is

        rho_R(exon)  >=  rho_nu(B) + rho_mu(B)          (a LOWER BOUND)

asserted as an equality (`phi == 1` to 2.2e-16) — `ROADMAP.md`, the `omega_graft` debt. The junction flux at
a boundary measures the RNA that CROSSES it; an exon also carries mature fragments that never cross, and
those are invisible to the graft. So the graft must under-claim, and the shortfall is exactly the mass the
solver has no RNA explanation for.

This tests that directly, and tests the structural split the debt predicts: the miss should be worst where
RNA does not flow through the seam — a transcript TERMINUS, which the region map cannot currently see (P1g).
A terminus is detectable here as a flanking boundary carrying NO spliced flux.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d6_graft.py [cond]
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from gdna_d5_zero_source import load  # noqa: E402

from rigel.calibration.node_geometry import _node_region_type  # noqa: E402

_EPS = 1e-9


def main():
    index, ra, inp, dbg, cc = load()
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    fg = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    eff = np.asarray(cap["eff_global"])
    spl = np.asarray(st.mass_spliced, float)
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    pp, pn = np.asarray(cap["prec_p"], float), np.asarray(cap["prec_n"], float)
    mp, mn = np.asarray(cap["mode_p"], float), np.asarray(cap["mode_n"], float)
    has_r = (pp > 0) | (pn > 0)
    live = mass > _EPS
    fp = fg * mass
    claim = np.exp(np.clip(np.where(pp > 0, mp, np.where(pn > 0, mn, np.nan)), -50.0, 0.0))

    def flank(idx, arr):
        """The value of `arr` at each node's left / right neighbour (0 at a reference terminal)."""
        lo = np.where(idx >= 0, arr[np.clip(idx, 0, len(arr) - 1)], 0.0)
        return lo

    spl_l, spl_r = flank(left, spl), flank(right, spl)
    # A boundary with NO spliced flux carries no junction: on that side the exon is a transcript TERMINUS
    # (or abuts an unexpressed neighbour). This is exactly the structural bit the region map lacks — P1g.
    n_open = (spl_l <= _EPS).astype(int) + (spl_r <= _EPS).astype(int)

    sel = live & (cls == 2) & has_r
    print(f"# {len(sel)} nodes; exons with an RNA message: {int(sel.sum())}\n")
    print("=== A. The graft's claim vs the truth, by how many flanks carry NO junction flux ===")
    print("    (true f_rna = 1.000 everywhere — this library has no gDNA at all)\n")
    print(f"{'flanks without a junction':28s} {'n':>5s} {'mass':>12s} {'claimed f_rna':>14s} "
          f"{'f_g':>8s} {'FP mass':>12s} {'share FP':>9s}")
    tot = fp[live].sum()
    for k in (0, 1, 2):
        m = sel & (n_open == k)
        if not m.any():
            continue
        w = mass[m]
        print(f"{('both junctions' if k == 0 else 'one open (TERMINUS)' if k == 1 else 'both open'):28s} "
              f"{int(m.sum()):5d} {mass[m].sum():12,.0f} {np.average(claim[m], weights=w):14.3f} "
              f"{np.average(fg[m], weights=w):8.3f} {fp[m].sum():12,.0f} "
              f"{100 * fp[m].sum() / max(tot, 1):8.1f} %")

    print("\n=== B. Is the shortfall the NON-CROSSING mature mass? ===")
    print("    The graft can only see RNA that crosses a seam. Compare each exon's own unspliced DENSITY")
    print("    with the density the graft claims for it; if the graft is a lower bound, the ratio is < 1")
    print("    and it should worsen as the exon gets longer relative to the fragment length.\n")
    dens = mass / np.maximum(eff, _EPS)
    q = np.percentile(eff[sel], [0, 20, 40, 60, 80, 100])
    print(f"{'exon eff-length bin':26s} {'n':>5s} {'med eff':>9s} {'claimed f_rna':>14s} {'f_g':>8s} "
          f"{'FP %':>7s}")
    for i in range(5):
        b = sel & (eff >= q[i]) & (eff <= q[i + 1] if i == 4 else eff < q[i + 1])
        if not b.any():
            continue
        w = mass[b]
        print(f"[{q[i]:10.0f},{q[i + 1]:10.0f}) {int(b.sum()):5d} {np.median(eff[b]):9.0f} "
              f"{np.average(claim[b], weights=w):14.3f} {np.average(fg[b], weights=w):8.3f} "
              f"{100 * fp[b].sum() / max(mass[b].sum(), 1):7.1f}")
    _ = dens

    print("\n=== C. The same node population on a library that DOES have nascent RNA ===")
    print("    `nrna_none` is the worst stratum (FP 30.7 % of exon mass vs 9.6 % with nascent RNA).")
    print("    With nascent RNA the introns carry mass, get solved by the intron factory, and can emit —")
    print("    so the exon has a SECOND opinion. Without it, the graft is the only voice in the room.")


if __name__ == "__main__":
    main()
