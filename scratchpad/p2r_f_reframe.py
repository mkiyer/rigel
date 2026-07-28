"""P-2 RESIDUAL, PART F — is the exposed defect the REFRAME applied to the gDNA component?

E1/E2 established that on the BP-legal branch the conservation violation is tiny (|δ| = 0.07) — so no legal
pin recovers the residual — while on the BLOCKED branch the old pin was applying a **×0.648** shrink to
every delivered component of a message whose budget it had completed from the destination's own density.

What that shrink was masking is visible in the numbers: the delivered gDNA LEVEL on the harmed exons is
`e^moG` = 0.42–0.77 against an oracle `f_g` of 0.008–0.043 — **10–50× too large**.

THE HYPOTHESIS.  The reframe `r = ρ_tot(dst)/ρ_tot(src)` exists to cancel the hybrid-capture step.  Between
an INTRON (or boundary) and an EXON in a gDNA-bearing library the total-density ratio is dominated by RNA,
not by capture — so reframing the **gDNA** component by it multiplies the gDNA claim by the destination's
RNA content.  gDNA is uniform, so its correct transfer across a capture-OFF hop is `r = 1`.

THE TEST, and it is sharp:  if the reframe is the culprit, then on the harmed exons

        ρ_g^msg / ρ_g^true(dst)   ≈   r        (and ρ_g^msg/r ≈ ρ_g^true, i.e. the UNREFRAMED claim is right)

FALSIFIER: if the un-reframed claim is no closer to the truth than the reframed one, the reframe is not the
mechanism and the residual has another cause.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_f_reframe.py
"""
from __future__ import annotations

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
CL = ("intergenic", "intron", "exon", "boundary")
CONDS = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
]


def main():
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    print("=== F1. the gDNA claim delivered into EXONS: reframed vs un-reframed, against the truth ===")
    print("   err = mass-weighted mean |log10(claim / true rho_g)|.  Lower is better.\n")
    print(f"{'condition':44s} {'n':>6s} {'msg mass':>11s} {'med r':>7s} {'err REFRAMED':>13s} "
          f"{'err UNREFRAMED':>15s} {'src=intron/bnd':>15s}")
    for cond in CONDS:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        dbg: dict = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
        cap, chain = dbg["capture"], dbg["chain"]
        S = cap["_uni_static"]
        E_g = np.asarray(S["E_g"], float)
        mass = np.asarray(cap["mass_global"], float)
        rt, _ = _node_region_type(chain, ra)
        cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        rg_true = (Gp + Gn) / np.maximum(E_g, _EPS)
        rl = np.asarray(S["rho_l0"], float)
        rr = np.asarray(S["rho_r0"], float)
        acc = {k: [] for k in ("re", "un", "w", "r", "sc")}
        for pn_, (dstf, srcf) in zip(cap["_pin"][-2:], ((rl, rr), (rr, rl))):
            src = np.asarray(pn_["src"], np.int64)
            valid = np.asarray(pn_["valid"], bool)
            tg = np.asarray(pn_["tg"], float)
            tpg = np.asarray(pn_["tpg"], float)
            framed = valid & (srcf[src] > _EPS) & (dstf > _EPS)
            r = np.where(framed, dstf / np.maximum(srcf[src], _EPS), 1.0)
            sel = valid & (tpg > 0.0) & (tg > _EPS) & (cls == 2) & (rg_true > _EPS) & (mass > _EPS)
            if not sel.any():
                continue
            acc["re"].append(np.abs(np.log10(tg[sel] / rg_true[sel])))
            acc["un"].append(np.abs(np.log10(tg[sel] / np.maximum(r[sel], _EPS) / rg_true[sel])))
            acc["w"].append(mass[sel]); acc["r"].append(r[sel])
            acc["sc"].append(np.isin(cls[src][sel], (1, 3)))
        if not acc["re"]:
            continue
        re, un, w, r, sc = (np.concatenate(acc[k]) for k in ("re", "un", "w", "r", "sc"))
        print(f"{cond[5:]:44s} {re.size:6d} {w.sum():11,.0f} {np.median(r):7.2f} "
              f"{np.average(re, weights=w):13.3f} {np.average(un, weights=w):15.3f} "
              f"{np.average(sc, weights=w):14.1%}")

    print("\n\n=== F2. the same, split by whether the hop crosses a capture step ===")
    print("   the reframe is SUPPOSED to carry the capture step; the question is whether it carries")
    print("   anything else.  capture OFF ⇒ the true enrichment ratio is 1, so any r != 1 is spurious.\n")
    print(f"{'condition':44s} {'r p10':>7s} {'r p50':>7s} {'r p90':>7s} "
          f"{'|log r| ':>9s} {'improve if r dropped':>21s}")
    for cond in CONDS:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        dbg: dict = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
        cap, chain = dbg["capture"], dbg["chain"]
        S = cap["_uni_static"]
        E_g = np.asarray(S["E_g"], float)
        mass = np.asarray(cap["mass_global"], float)
        rt, _ = _node_region_type(chain, ra)
        cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        rg_true = (Gp + Gn) / np.maximum(E_g, _EPS)
        rl, rr = np.asarray(S["rho_l0"], float), np.asarray(S["rho_r0"], float)
        rs, ws, im = [], [], []
        for pn_, (dstf, srcf) in zip(cap["_pin"][-2:], ((rl, rr), (rr, rl))):
            src = np.asarray(pn_["src"], np.int64)
            valid = np.asarray(pn_["valid"], bool)
            tg, tpg = np.asarray(pn_["tg"], float), np.asarray(pn_["tpg"], float)
            framed = valid & (srcf[src] > _EPS) & (dstf > _EPS)
            r = np.where(framed, dstf / np.maximum(srcf[src], _EPS), 1.0)
            sel = valid & (tpg > 0.0) & (tg > _EPS) & (cls == 2) & (rg_true > _EPS) & (mass > _EPS)
            if not sel.any():
                continue
            rs.append(r[sel]); ws.append(mass[sel])
            im.append(np.abs(np.log10(tg[sel] / rg_true[sel]))
                      - np.abs(np.log10(tg[sel] / np.maximum(r[sel], _EPS) / rg_true[sel])))
        if not rs:
            continue
        rs, ws, im = np.concatenate(rs), np.concatenate(ws), np.concatenate(im)
        print(f"{cond[5:]:44s} {np.quantile(rs, 0.1):7.2f} {np.median(rs):7.2f} "
              f"{np.quantile(rs, 0.9):7.2f} {np.average(np.abs(np.log10(rs)), weights=ws):9.3f} "
              f"{np.average(im > 0, weights=ws):20.1%}")


if __name__ == "__main__":
    main()
