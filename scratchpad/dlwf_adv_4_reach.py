"""ADVERSARY #4 — REACH of the hazards found in adv #2, measured on real solver state, plus the
_RHO_ITERS kill-decision FLIP rate (the DL hinge is discrete, so a small iteration-2 reframe shift can
flip a message from kept to killed).

  Q1  what fraction of DL-LIVE destinations receive a PARTIAL message (adv#2 P2/P4: r does not cancel,
      dG/dlog(o) is not -1)?
  Q2  what fraction receive a PEEL (adv#2 P3) or a GRAFT (a junction COUNT in the measurement stream that
      the plan deflates by a COMPOSITION-mismatch estimate)?
  Q3  the entrenchment set: destinations with o_c == 0 on a component whose message precision is > 0.
  Q4  iteration flip: does the kill decision change between _RHO_ITERS 1 and 2?
  Q5  how much of the message channel survives DL, by MASS (the honest 'is pass-0 still a BP solver?').

    OMP_NUM_THREADS=1 python scratchpad/dlwf_adv_4_reach.py
"""

from __future__ import annotations

import dataclasses
import importlib
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
SCRATCH = Path("/Users/mkiyer/proj/rigel/scratchpad/_dlwf_work")
_EPS = 1.0e-9
CONDS = (
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna5_ss_0.50_nrna_present_capture_verystrong",
)


def main():
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from selfsolve_diag import _scan_and_truth

    SCRATCH.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    calmod = importlib.import_module("rigel.calibration.calibrate")
    os.environ["RIGEL_S2T_OFF"] = "1"
    for cond in CONDS:
        inp = _scan_and_truth(SUITE, cond, index, cfg, SCRATCH, SUITE / "_selfsolve_cache")
        dbg = {}
        calmod.calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
        )
        cap = dbg["capture"]
        st = cap["_uni_static"]
        uni = cap["_uni"]
        n = st["M"].shape[0]
        M, E_g = np.asarray(st["M"], float), np.asarray(st["E_g"], float)
        og, op = st["og"], st["op"]
        fg_own = np.asarray(cap["fg_loc"], float)
        tau = np.asarray(cap["_tau0_lam"], float)
        lock = np.asarray(st["struct_lock"], bool)
        solvable = np.asarray(cap["solvable"], bool)
        is_bnd, ex_a = np.asarray(st["is_bnd"], bool), np.asarray(st["is_exon"], bool)
        left, right = np.asarray(st["left"]), np.asarray(st["right"])
        has_own = lock | (tau > _EPS)
        base = solvable & has_own
        print(f"\n{'='*100}\nCONDITION {cond}")
        print(f"[Q0] DL-live destinations (solvable & finite v_own): {int(base.sum())} / {int(solvable.sum())} solvable")

        for tag, src, kg, kp, pg_, pp_ in (
            ("L", left, "ag", "ap", "apg", "app"), ("R", right, "bg", "bp", "bpg", "bpp"),
        ):
            u = uni[-1]
            valid = (src >= 0) & base
            s = np.clip(src, 0, n - 1)
            graft = ex_a & is_bnd[s] & valid
            peel = is_bnd & ex_a[s] & valid
            pg = np.asarray(u[pg_], float) > 0.0
            pp = np.asarray(u[pp_], float) > 0.0
            tg, tp = np.asarray(u[kg], float), np.asarray(u[kp], float)
            partial = valid & (pg ^ pp)
            print(f"[Q1/Q2] msg {tag}: DL-live valid n={int(valid.sum()):5d}"
                  f"  PARTIAL (r survives the pin) {100.0*partial.sum()/max(valid.sum(),1):5.1f}%"
                  f"  PEEL {100.0*peel.sum()/max(valid.sum(),1):5.1f}%"
                  f"  GRAFT (junction COUNT in the measurement stream) {100.0*graft.sum()/max(valid.sum(),1):5.1f}%")
            ent_p = valid & pp & (np.asarray(op, float) <= _EPS)
            ent_g = valid & pg & (np.asarray(og, float) <= _EPS)
            zt_p = valid & pp & (tp <= _EPS)
            print(f"[Q3]    msg {tag}: ENTRENCHED (o_c==0 but a live message on that component): "
                  f"gDNA {int(ent_g.sum()):4d}  RNA+ {int(ent_p.sum()):4d}"
                  f"   |  message density t_c==0 with a live precision: RNA+ {int(zt_p.sum()):4d} "
                  f"({100.0*zt_p.sum()/max((valid&pp).sum(),1):4.1f}% of live RNA+ msgs)")

        # Q4 — iteration flip of the kill decision
        v_own_g = np.where(lock, 0.0, np.where(tau > _EPS, (1 - fg_own) ** 2 / np.maximum(tau, _EPS), np.inf))
        flips = {}
        for tag, kg, pg_ in (("L", "ag", "apg"), ("R", "bg", "bpg")):
            dec = []
            for it in range(len(uni)):
                t = np.asarray(uni[it][kg], float)
                p = np.asarray(uni[it][pg_], float)
                G = np.log(np.maximum(t, _EPS) / np.maximum(og, _EPS))
                vm = np.where(p > 0, 1.0 / np.maximum(p, _EPS), np.inf)
                exc = np.where(np.isfinite(v_own_g), np.maximum(0.0, G * G - vm - v_own_g), 0.0)
                damp = np.where(p > 0, 1.0 / (1.0 + exc * p), 1.0)
                dec.append((p > 0) & base & (damp < 0.5))
            if len(dec) >= 2:
                flips[tag] = float(np.mean(dec[0][base] != dec[-1][base]))
        print(f"[Q4]    kill-decision FLIP rate between _RHO_ITERS 1 and 2: "
              + "  ".join(f"{k}={v*100:.2f}%" for k, v in flips.items()))

        # Q5 — surviving message channel by mass
        for tag, kg, pg_, kp, pp_ in (("L", "ag", "apg", "ap", "app"), ("R", "bg", "bpg", "bp", "bpp")):
            u = uni[-1]
            tot_g = tot_p = sur_g = sur_p = 0.0
            for k_, p_, o_, vo in ((kg, pg_, og, v_own_g),
                                   (kp, pp_, op, np.where(lock, 0.0, np.where(tau > _EPS, fg_own**2 / np.maximum(tau, _EPS), np.inf)))):
                t = np.asarray(u[k_], float)
                p = np.asarray(u[p_], float)
                live = (p > 0) & base
                G = np.log(np.maximum(t, _EPS) / np.maximum(np.asarray(o_, float), _EPS))
                vm = np.where(p > 0, 1.0 / np.maximum(p, _EPS), np.inf)
                exc = np.where(np.isfinite(vo), np.maximum(0.0, G * G - vm - vo), 0.0)
                pn = np.where(p > 0, 1.0 / (1.0 / np.maximum(p, _EPS) + exc), 0.0)
                if k_ == kg:
                    tot_g, sur_g = float(np.sum(p[live] * M[live])), float(np.sum(pn[live] * M[live]))
                else:
                    tot_p, sur_p = float(np.sum(p[live] * M[live])), float(np.sum(pn[live] * M[live]))
            print(f"[Q5]    msg {tag}: mass-weighted message PRECISION surviving DL: "
                  f"gDNA {100.0*sur_g/max(tot_g,_EPS):5.2f}%   RNA+ {100.0*sur_p/max(tot_p,_EPS):5.2f}%")


if __name__ == "__main__":
    main()
