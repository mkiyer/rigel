"""P-2 RESIDUAL, PART E — WHICH component contaminates the pin's budget?

The `w ≡ 1` control restores 52 % of `mo_g` and 86 % of `mo_p` but only ~11 % of the error, so the
BP-legality gate (`_clean`: no unsupplied component may lend the destination's own density to the budget)
blocks about half the messages the old pin acted on.  Which component does the blocking, and how much
destination mass is behind each?

If the blocked messages are blocked by the **gDNA** arm (`prec_g == 0` while `og > 0`), then the pin was
borrowing the destination's own gDNA density to set a gDNA level — the exact self-confirmation
`pin_derivation.md` §5 puts in closed form — and that part of the residual is the HONEST PRICE of removing
the violation, not a recoverable defect.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_e_gate.py [--cond ...]
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", nargs="*", default=CONDS)
    a = ap.parse_args()
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    print("=== E1. the pin's budget, per MESSAGE: is it BP-legal, and if not, which arm borrows? ===")
    print("   shares are of DESTINATION MASS over messages that carry any precision at all.\n")
    print(f"{'condition':44s} {'msgs':>7s} {'CLEAN':>7s} {'gDNA':>7s} {'RNA+':>7s} {'RNA-':>7s} "
          f"{'|delta| clean':>13s} {'|delta| dirty':>13s}")
    keep = {}
    for cond in a.cond:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        dbg: dict = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
        cap, chain = dbg["capture"], dbg["chain"]
        S = cap["_uni_static"]
        og, op, on = (np.asarray(S[k], float) for k in ("og", "op", "on"))
        mass = np.asarray(cap["mass_global"], float)
        rt, _ = _node_region_type(chain, ra)
        cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        fo = np.where(Gp + Gn + Rp + Rn > _EPS, (Gp + Gn) / np.maximum(Gp + Gn + Rp + Rn, _EPS), np.nan)
        # the two per-message pre-pin states (left msg then right msg)
        pins = cap["_pin"][-2:]
        p3s = cap["_p3"][-2:]
        tw = ncl = 0.0
        blk = np.zeros(3)
        dc, dd, wc, wd = [], [], [], []
        for pn_, p3 in zip(pins, p3s):
            sup = np.stack([pn_["tpg"], pn_["tpp"], pn_["tpn"]], -1) > 0.0
            own = np.stack([og, op, on], -1)
            live = np.asarray(pn_["valid"], bool) & sup.any(-1) & (mass > _EPS)
            borrow = (~sup) & (own > 0.0)
            clean = np.asarray(p3["clean"], bool)
            tw += mass[live].sum()
            ncl += mass[live & clean].sum()
            for j in range(3):
                blk[j] += mass[live & ~clean & borrow[..., j]].sum()
            d = np.asarray(p3["dlt"], float)
            dc.append(np.abs(d)[live & clean]); wc.append(mass[live & clean])
            dd.append(np.abs(d)[live & ~clean]); wd.append(mass[live & ~clean])
        dc, wc, dd, wd = (np.concatenate(x) for x in (dc, wc, dd, wd))
        print(f"{cond[5:]:44s} {tw:7,.0f} {100 * ncl / max(tw, 1):6.1f}% "
              f"{100 * blk[0] / max(tw, 1):6.1f}% {100 * blk[1] / max(tw, 1):6.1f}% "
              f"{100 * blk[2] / max(tw, 1):6.1f}% "
              f"{np.average(dc, weights=wc) if wc.sum() else np.nan:13.4f} "
              f"{np.average(dd, weights=wd) if wd.sum() else np.nan:13.4f}")
        keep[cond] = (cap, chain, cls, fo, mass, og, op, on, pins, p3s)

    print("\n   CLEAN + gDNA + RNA+ + RNA- may exceed 100 %: a message can borrow on more than one arm.")

    print("\n\n=== E2. the BLOCKED messages: what would the pin have done, and to whom? ===")
    print("   'pin factor' = M/S, the multiplier the old pin applied to every supplied component.\n")
    print(f"{'condition':44s} {'class':10s} {'msg mass':>11s} {'med M/S':>8s} {'p10':>7s} {'p90':>7s} "
          f"{'oracle f_g':>11s} {'own f_g':>8s}")
    CL = ("intergenic", "intron", "exon", "boundary")
    for cond, (cap, chain, cls, fo, mass, og, op, on, pins, p3s) in keep.items():
        E_g = np.asarray(cap["_uni_static"]["E_g"], float)
        E_r = np.asarray(cap["_uni_static"]["E_r"], float)
        M = np.asarray(cap["_uni_static"]["M"], float)
        r0 = np.asarray(cap["_uni_static"]["rho_node0"], float)
        for ci in (2, 3):
            ks, ws, fos, fgo = [], [], [], []
            for pn_, p3 in zip(pins, p3s):
                sup = np.stack([pn_["tpg"], pn_["tpp"], pn_["tpn"]], -1) > 0.0
                live = np.asarray(pn_["valid"], bool) & sup.any(-1) & (mass > _EPS)
                sel = live & ~np.asarray(p3["clean"], bool) & (cls == ci)
                if not sel.any():
                    continue
                mc = np.where(sup, np.stack([pn_["tg"], pn_["tp"], pn_["tn"]], -1),
                              np.stack([og, op, on], -1)) * np.stack([E_g, E_r, E_r], -1)
                Ssum = mc.sum(-1)
                k = np.where(Ssum > _EPS, M / np.maximum(Ssum, _EPS), 1.0)
                ks.append(k[sel]); ws.append(mass[sel]); fos.append(fo[sel])
                fgo.append((og / np.maximum(r0, _EPS))[sel])
            if not ks:
                continue
            ks, ws, fos, fgo = (np.concatenate(x) for x in (ks, ws, fos, fgo))
            print(f"{cond[5:]:44s} {CL[ci]:10s} {ws.sum():11,.0f} {np.median(ks):8.3f} "
                  f"{np.quantile(ks, 0.1):7.3f} {np.quantile(ks, 0.9):7.3f} "
                  f"{np.nanmean(fos):11.3f} {np.median(fgo):8.3f}")


if __name__ == "__main__":
    main()
