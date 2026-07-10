"""Interrogate the per-strand RNA IMPUTATION around the AMBIG locus (gDNA=0 + nascent, ss0.99).

User's question: the single-strand exon is screaming RNA via strand tilt (u_pos=5/u_neg=622); the AMBIG
node's matching-strand flanks are RNA-rich. The abundant RNA evidence SHOULD dominate and flow OUT to help
the AMBIG node. What are the imputation passes actually doing?

This logs every `_message` call in the LAST sweep pass (gDNA = 6 positional args, RNA = 7), matches each to
its destination node by mass, and prints the per-strand RNA + gDNA messages flowing INTO the POS exon and the
AMBIG node — with the source vs destination densities so the (wrong) direction is visible.

    OMP_NUM_THREADS=1 python scripts/debug/rna_imputation_interrogate.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambig_nascent_trace import build  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402

SC = {-1: "?", 0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)

        log = {"pass": 0, "msgs": []}
        orig_msg = bp._message
        orig_nd = bp.node_densities
        cap = {}
        orig_sweep = bp.node_sweep

        def nd_spy(*a, **k):
            log["pass"] += 1
            log["msgs"].clear()  # keep only the latest pass
            return orig_nd(*a, **k)

        def msg_spy(*args, **k):
            out = orig_msg(*args, **k)
            kind = "gDNA" if len(args) == 6 else "RNA"
            # args = (rho_src, eff_src, eff_dst, mass_dst, rho_dst_cur, varmean[, spliced_dst])
            log["msgs"].append((kind, float(args[0]), float(args[3]), float(args[4]), out[0], out[1]))
            return out

        def sweep_spy(chain, statics, geometry, belief, *a, **k):
            cap["chain"] = chain
            cap["geom"] = geometry
            cap["statics"] = statics
            out = orig_sweep(chain, statics, geometry, belief, *a, **k)
            cap["belief"] = out[0]
            return out

        bp.node_densities = nd_spy
        bp._message = msg_spy
        bp.node_sweep = sweep_spy
        calibrate.__globals__["node_sweep"] = sweep_spy
        res = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
        bp.node_densities = orig_nd
        bp._message = orig_msg
        bp.node_sweep = orig_sweep
        calibrate.__globals__["node_sweep"] = orig_sweep

        chain = cap["chain"]
        geom = cap["geom"]
        st = cap["statics"]
        bel = cap["belief"]
        kind = np.asarray(chain.kind)
        ref_idx = np.asarray(chain.ref_idx)
        is_reg = kind == REGION
        start = np.asarray(ra.start)
        sc_arr = np.asarray(ra.strand_class)
        snap = bp.node_densities(bel, geom)
        kappa = res.rna_sense_frac

        # node table around the AMBIG locus (3000–8000)
        def region_node(s0):
            r = int(np.where(start == s0)[0][0])
            return int(np.where(is_reg & (ref_idx == r))[0][0]), r

        print(f"\n=== RNA imputation around the AMBIG locus (gDNA=0, nascent=30, ss0.99, κ={kappa:.4f}) ===")
        print(f"  {'node':>20} {'cls':>6} {'u_pos':>6} {'u_neg':>6} "
              f"{'f_pos':>6} {'f_neg':>6} {'f_g':>6} {'ρ_pos':>7} {'ρ_neg':>7} {'ρ_g':>6}")
        nodes = []
        for s0 in (1500, 4000, 5000, 6000, 7000):
            hit = np.where(start == s0)[0]
            if hit.size == 0:
                continue
            node, r = region_node(s0)
            nodes.append((s0, node, r))
            print(f"  {f'{s0}-{int(np.asarray(ra.end)[r])}':>20} {SC.get(int(sc_arr[r]), '?'):>6} "
                  f"{st.u_pos[node]:>6.0f} {st.u_neg[node]:>6.0f} "
                  f"{bel.f_pos[node]:>6.3f} {bel.f_neg[node]:>6.3f} {bel.f_g[node]:>6.3f} "
                  f"{snap.rho_pos_left[node]:>7.2f} {snap.rho_neg_left[node]:>7.2f} {snap.rho_g_left[node]:>6.3f}")

        # messages into the POS exon (4000) and the AMBIG (5000), matched by dest mass
        mass_l = np.asarray(geom.mass_left)
        for s0, label in ((4000, "POS exon 4000-5000"), (5000, "AMBIG 5000-6000")):
            node, r = region_node(s0)
            mdst = mass_l[node]
            print(f"\n  --- messages INTO {label} (own f_pos={bel.f_pos[node]:.3f} "
                  f"f_neg={bel.f_neg[node]:.3f} f_g={bel.f_g[node]:.3f}) ---")
            seen = set()
            for knd, rho_src, mass_dst, rho_dst_cur, mu_c, n_src in log["msgs"]:
                if abs(mass_dst - mdst) > 1e-6:
                    continue
                key = (knd, round(rho_src, 4), round(mu_c, 4))
                if key in seen:
                    continue
                seen.add(key)
                print(f"      {knd:>4} msg: ρ_src={rho_src:>7.2f}  ρ_dst_cur={rho_dst_cur:>7.2f}  "
                      f"→ μ_c={mu_c:>6.3f}  N_src={n_src:>6.2f}")
        sc.cleanup()


if __name__ == "__main__":
    main()
