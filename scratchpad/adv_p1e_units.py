"""ADVERSARIAL unit / frame checks on P1e's statistic.

 (a) ``1/n_dst`` is the destination's POISSON relative variance of ``M``. It is only that if ``M`` IS the
     count ``n_dst``. `bp_solver.py:584` warns the spliced MASS is ~2x the COUNT; check the ratio the code
     actually uses (``M`` vs ``_n_node``, per node class).

 (b) the INTERGENIC constant offset: delta = -0.57 (gdna300) / -0.70 (gdna100) on 100 % of messages. Is it a
     frame/units artefact (a fixed ratio) or a real density disagreement? Decompose S into its components.

    OMP_NUM_THREADS=1 python scratchpad/adv_p1e_units.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

for name in ("gdna300_ss0.50_present_capOFF", "gdna100_ss0.50_none_capOFF"):
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    print(f"\n{'=' * 110}\n{name}\n{'=' * 110}")
    print("  (a)  M / n_dst   [1 ⇒ mass IS the count ⇒ 1/n_dst is the right Poisson term]")
    is_reg = nf["kind"] == 1  # REGION == 1? printed below for safety
    nd = np.where(is_reg, nf["n_unspl_l"], nf["n_unspl_l"] + nf["n_unspl_r"])
    for cls in ("exon", "boundary", "intron", "intergenic"):
        m = (nf["cls"] == cls) & (nd > 0)
        if not m.any():
            continue
        r = nf["M"][m] / nd[m]
        print(f"    {cls:<12} n={m.sum():>5}  median M/n = {np.median(r):.4f}  "
              f"p10 {np.percentile(r, 10):.4f}  p90 {np.percentile(r, 90):.4f}")

    print("\n  (b)  INTERGENIC destinations: what is S made of?")
    sel = (t["cls"] == "intergenic") & (t["nsup"] > 0)
    if sel.any():
        print(f"    n={sel.sum()}  median delta={np.median(t['delta'][sel]):+.4f}  "
              f"median S/M={np.median(t['S'][sel] / t['M'][sel]):.4f}")
        print(f"    median alpha_g={np.median(t['alpha_g'][sel]):.4f}  "
              f"alpha_p={np.median(t['alpha_p'][sel]):.4f}  alpha_n={np.median(t['alpha_n'][sel]):.4f}")
        print(f"    supplied gDNA on {100 * np.mean(t['sup_g'][sel]):.1f}%  "
              f"supplied RNA on {100 * np.mean(t['sup_r'][sel]):.1f}%  "
              f"lambda-emitting {100 * np.mean(t['lam_emit'][sel]):.1f}%")
        d = t["dst"][sel]
        print(f"    dst E_g/E_r median = {np.median(nf['E_g'][d] / np.maximum(nf['E_r'][d], 1e-9)):.4f}")
        print(f"    dst oracle f_g median = {np.median(nf['fo'][d]):.4f}   "
              f"dst own og*E_g/M median = "
              f"{np.median(nf['og'][d] * nf['E_g'][d] / np.maximum(nf['M'][d], 1e-9)):.4f}")
        src = t["src"][sel]
        print(f"    src class counts: "
              f"{dict(zip(*np.unique(nf['cls'][src], return_counts=True), strict=False))}")
        # the supplied gDNA claim alone, as a share of M
        ag = t["alpha_g"][sel] * t["S"][sel] / np.maximum(t["M"][sel], 1e-9)
        print(f"    a_g = rho_g*E_g/M   median {np.median(ag):.4f}  p10 {np.percentile(ag, 10):.4f}  "
              f"p90 {np.percentile(ag, 90):.4f}")
