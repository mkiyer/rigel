"""Look at the raw data: the AMBIG region's counts + the two flanking boundaries' per-side counts.

Transcript structure (gDNA=0, nascent=30, ss0.99):
  TA+ exon2 = 4000-6000 (CONTINUOUS), TB- exon1 = 5000-7000 (CONTINUOUS).
  region 4000-5000 = TA only (POS); 5000-6000 = TA∩TB (AMBIG); 6000-7000 = TB only (NEG).
  boundary 5000 = mid-TA-exon (+ crosses) AND TB exon start; boundary 6000 = TA exon end AND mid-TB-exon (- crosses).

Dumps the RAW per-side counts so we can see whether the boundary actually carries the strand-specific RNA
crossing signal, and what the boundary node solves to.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_boundary_counts.py
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
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402

SC = {-1: "?", 0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def vrow(view, i):
    return (float(view.n_unspliced_pos[i]), float(view.n_unspliced_neg[i]),
            float(view.n_spliced_sense[i]), float(view.mass_unspliced[i]), float(view.mass_spliced[i]))


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)
        cs = CalibrationSubstrate.from_payload(pl, ra)
        bs = BoundarySubstrate.from_payload(pl)

        cap = {}
        orig = bp.node_sweep

        def spy(chain, statics, geometry, belief, *a, **k):
            cap["chain"] = chain
            cap["geom"] = geometry
            out = orig(chain, statics, geometry, belief, *a, **k)
            cap["belief"] = out[0]
            return out

        bp.node_sweep = spy
        calibrate.__globals__["node_sweep"] = spy
        res = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
        bp.node_sweep = orig
        calibrate.__globals__["node_sweep"] = orig

        chain = cap["chain"]
        geom = cap["geom"]
        bel = cap["belief"]
        kind = np.asarray(chain.kind)
        ref_idx = np.asarray(chain.ref_idx)
        cleft = np.asarray(chain.left)
        cright = np.asarray(chain.right)
        is_reg = kind == REGION
        start = np.asarray(ra.start)
        end = np.asarray(ra.end)
        sc_arr = np.asarray(ra.strand_class)
        blr = np.asarray(bs.left_region)
        brr = np.asarray(bs.right_region)

        amb_r = int(np.where(start == 5000)[0][0])
        amb_node = int(np.where(is_reg & (ref_idx == amb_r))[0][0])
        lb_node = int(cleft[amb_node])   # boundary at 5000
        rb_node = int(cright[amb_node])  # boundary at 6000

        print(f"\n=== AMBIG region 5000-6000 (κ={res.rna_sense_frac:.4f}) ===")
        up, un, ss_, mu, ms = vrow(cs.contained, amb_r)
        print(f"  strand_class={SC.get(int(sc_arr[amb_r]),'?')}  contained: "
              f"u_pos={up:.0f} u_neg={un:.0f} spliced_sense={ss_:.0f}  mass_unspl={mu:.1f} mass_spl={ms:.1f}")
        print(f"  SOLVED: f_pos={bel.f_pos[amb_node]:.3f} f_neg={bel.f_neg[amb_node]:.3f} f_g={bel.f_g[amb_node]:.3f}")

        for bnode, label in ((lb_node, "boundary 5000 (POS-exon | AMBIG)"),
                             (rb_node, "boundary 6000 (AMBIG | NEG-exon)")):
            bi = int(ref_idx[bnode])
            lr, rr = int(blr[bi]), int(brr[bi])
            print(f"\n=== {label}  [node {bnode}, bidx {bi}] ===")
            print(f"  left_region={lr}({SC.get(int(sc_arr[lr]),'?') if lr>=0 else '-'} "
                  f"{start[lr] if lr>=0 else '-'}-{end[lr] if lr>=0 else '-'})  "
                  f"right_region={rr}({SC.get(int(sc_arr[rr]),'?') if rr>=0 else '-'} "
                  f"{start[rr] if rr>=0 else '-'}-{end[rr] if rr>=0 else '-'})")
            for side, view in (("left side", bs.left), ("right side", bs.right)):
                up, un, ss_, mu, ms = vrow(view, bi)
                print(f"  {side:>10}: u_pos={up:.0f} u_neg={un:.0f} spliced_sense={ss_:.0f}  "
                      f"mass_unspl={mu:.1f} mass_spl={ms:.1f}")
            print(f"  SOLVED belief: f_pos={bel.f_pos[bnode]:.3f} f_neg={bel.f_neg[bnode]:.3f} "
                  f"f_g={bel.f_g[bnode]:.3f}")
            print(f"  geometry: spliced_pos_l={geom.spliced_pos_left[bnode]:.1f} "
                  f"spliced_pos_r={geom.spliced_pos_right[bnode]:.1f} "
                  f"spliced_neg_l={geom.spliced_neg_left[bnode]:.1f} "
                  f"spliced_neg_r={geom.spliced_neg_right[bnode]:.1f}")
            print(f"            eff_rna_l={geom.eff_rna_left[bnode]:.1f} eff_rna_r={geom.eff_rna_right[bnode]:.1f} "
                  f"mass_l={geom.mass_left[bnode]:.1f} mass_r={geom.mass_right[bnode]:.1f}")
        sc.cleanup()


if __name__ == "__main__":
    main()
