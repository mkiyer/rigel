"""Trace the calibration error to INDIVIDUAL messages — watch a correct local solve get corrupted.

For every region node with oracle truth we compare three states:
  * ``fg_strand`` — the strand likelihood ALONE (no prior, no messages),
  * ``fg_loc``    — the message-free LOCAL solve (strand + prior),
  * ``f_g``       — the FINAL belief (local ⊗ the forward/backward messages),
and attribute the damage: a node is **CORRUPTED** when the local solve was right and the messages made it
wrong. For those nodes we dump the competing precisions — the local precision (what the node uses to defend
itself) vs each incoming message's precision (gDNA / RNA+ / RNA−).

Two hypotheses this tests:
  1. **Strand too imprecise to hold ground** — local precision ≪ message precision at corrupted nodes.
  2. **The degrees-of-freedom double-count** — a SINGLE-STRAND node has ONE free parameter (f_g, with
     f_active = 1−f_g), yet it receives BOTH a gDNA message (on log f_g) AND an RNA message (on log(1−f_g)).
     Both are added to the same ψ and both derive from the same source belief, so their precisions ADD on one
     axis ⇒ messages get ~2× the intended weight and can overrule the strand.

    OMP_NUM_THREADS=1 python scripts/debug/msg_trace.py --cache-dir DIR --condition NAME
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from rigel.calibration import calibrate
from rigel.calibration.node_chain import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-12


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--condition", default="gdna_gdna100_ss_0.99_nrna_none_capture_off")
    ap.add_argument("--good", type=float, default=0.10, help="local |Δf_g| below this = local was CORRECT")
    ap.add_argument("--bad", type=float, default=0.30, help="final |Δf_g| above this = FINAL is WRONG")
    args = ap.parse_args()

    index = TranscriptIndex.load(str(Path(args.suite) / "rigel_index"))
    ra = RegionArrays.from_index(index)
    cfg = PipelineConfig()
    with open(Path(args.cache_dir) / f"{args.condition}.pkl", "rb") as fh:
        inp = pickle.load(fh)

    dbg: dict = {}
    calibrate(payload=inp["payload"], region_arrays=ra, strand_model=inp["strand_model"],
              gdna_fl_pmf=inp["gdna_fl_pmf"], rna_fl_pmf=inp["rna_fl_pmf"],
              config=cfg.calibration, _debug=dbg)
    cap, chain = dbg["capture"], dbg["chain"]

    # oracle per-region truth f_g on the contained unspliced mass → mapped onto region nodes
    p = inp["pools"]
    true_g = p["gdna_pos"] + p["gdna_neg"]
    tot = true_g + p["mat_uns_pos"] + p["nas_uns_pos"] + p["mat_uns_neg"] + p["nas_uns_neg"]
    with np.errstate(invalid="ignore", divide="ignore"):
        tfg_reg = np.where(tot > 0, true_g / np.maximum(tot, _EPS), np.nan)

    is_reg = np.asarray(chain.kind) == REGION
    ref = np.asarray(chain.ref_idx, np.int64)
    R = tfg_reg.shape[0]
    tfg = np.where(is_reg, tfg_reg[np.clip(ref, 0, R - 1)], np.nan)

    M = np.asarray(cap["mass_global"], float)
    fg_str, fg_loc, fg_fin = np.asarray(cap["fg_strand"]), np.asarray(cap["fg_loc"]), np.asarray(cap["f_g"])
    vg_loc = np.asarray(cap["vg_loc"], float)
    prec_g, prec_p, prec_n = (np.asarray(cap[k], float) for k in ("prec_g", "prec_p", "prec_n"))
    free_pos, free_neg = np.asarray(cap["free_pos"], bool), np.asarray(cap["free_neg"], bool)
    single = free_pos ^ free_neg  # single-strand ⇒ ONE degree of freedom
    solvable = np.asarray(cap["solvable"], bool)

    ok = is_reg & solvable & np.isfinite(tfg) & (M > 0)
    e_str, e_loc, e_fin = np.abs(fg_str - tfg), np.abs(fg_loc - tfg), np.abs(fg_fin - tfg)

    def wmae(e, m):
        return float(np.sum(M[m] * e[m]) / max(np.sum(M[m]), 1.0))

    print(f"=== {args.condition} ===")
    print(f"nodes: {int(ok.sum())} scored ({int((ok & single).sum())} single-strand, "
          f"{int((ok & free_pos & free_neg).sum())} AMBIG)")
    print("\n--- mass-weighted |Δf_g| by stage (does the message HELP or HURT?) ---")
    for nm, e in (("strand only", e_str), ("local (strand+prior)", e_loc), ("FINAL (+messages)", e_fin)):
        print(f"  {nm:22s} {wmae(e, ok):.4f}")

    corrupted = ok & (e_loc < args.good) & (e_fin > args.bad)
    rescued = ok & (e_loc > args.bad) & (e_fin < args.good)
    print(f"\n--- message verdict (local<{args.good} → final>{args.bad}) ---")
    print(f"  CORRUPTED by messages: {int(corrupted.sum())} nodes, {M[corrupted].sum():.3g} mass")
    print(f"  RESCUED  by messages: {int(rescued.sum())} nodes, {M[rescued].sum():.3g} mass")

    if corrupted.any():
        c = corrupted
        loc_prec = 1.0 / np.maximum(vg_loc, _EPS)
        msg_tot = prec_g + prec_p + prec_n
        print("\n--- at CORRUPTED nodes: can the node defend itself? ---")
        print(f"  local precision   median {np.median(loc_prec[c]):.3g}")
        print(f"  message precision median {np.median(msg_tot[c]):.3g}  "
              f"(gDNA {np.median(prec_g[c]):.3g} | RNA+ {np.median(prec_p[c]):.3g} | RNA− {np.median(prec_n[c]):.3g})")
        print(f"  msg/local ratio   median {np.median(msg_tot[c] / np.maximum(loc_prec[c], _EPS)):.3g}"
              "   (>1 ⇒ messages OVERRULE the node)")
        cs = c & single
        print("\n--- DOF double-count check (single-strand ⇒ ONE free axis) ---")
        print(f"  corrupted single-strand nodes: {int(cs.sum())}")
        if cs.any():
            both = cs & (prec_g > 0) & ((prec_p > 0) | (prec_n > 0))
            print(f"  of these, receiving BOTH a gDNA and an RNA message on that ONE axis: "
                  f"{int(both.sum())} ({100 * both.sum() / max(cs.sum(), 1):.0f}%)")
            rna = np.maximum(prec_p, prec_n)
            print(f"  their precision split: gDNA {np.median(prec_g[both]):.3g} + "
                  f"RNA {np.median(rna[both]):.3g} = {np.median(prec_g[both] + rna[both]):.3g} "
                  f"vs local {np.median(loc_prec[both]):.3g}")
        w = np.argsort(-(M * corrupted))[:6]
        print("\n--- worst corrupted nodes (by mass) ---")
        print(f"  {'node':>7} {'mass':>9} {'true':>6} {'strand':>7} {'local':>6} {'FINAL':>6} "
              f"{'loc_pr':>8} {'pr_g':>8} {'pr_+':>8} {'pr_−':>8} {'1DOF':>5}")
        for i in w:
            print(f"  {i:7d} {M[i]:9.1f} {tfg[i]:6.3f} {fg_str[i]:7.3f} {fg_loc[i]:6.3f} {fg_fin[i]:6.3f} "
                  f"{loc_prec[i]:8.3g} {prec_g[i]:8.3g} {prec_p[i]:8.3g} {prec_n[i]:8.3g} {str(bool(single[i])):>5}")


if __name__ == "__main__":
    main()
