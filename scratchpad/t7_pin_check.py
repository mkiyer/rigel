"""TARGET step 7 — does the DELIVERED gDNA mode depend on `r` at all?

`_pin_v` rescales each message so Sum_c rho_c*E_c = M(dst). When all three components are supplied that
scaling is proportional to 1/r, so `r` CANCELS EXACTLY and the delivered composition reduces to

    f_g : f_R  =  rho_g^src * E_g^dst  :  rho_R^src * E_r^dst

i.e. the SOURCE's composition re-expressed in the DESTINATION's effective lengths -- with no `r` in it.
(That identity was verified to 1.8e-15 during the M8 work.) If that holds here, then step 4's "the error is
the transport ratio" is an artifact of measuring the PRE-pin claim, no re-estimate of `r` can help those
nodes, and the mode error is the imputation premise itself.

Where it does NOT hold is where `_pin_v` substitutes the node's OWN density for an unsupplied component --
the gDNA-only / RNA-only strata. There `r` survives and is worth attacking.

    OMP_NUM_THREADS=1 python scratchpad/t7_pin_check.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_none_capture_on"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
calmod.calibrate(
    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
    np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
)
chain, cap = dbg["chain"], dbg["capture"]
uni, us = cap["_uni"][-1], cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
T = Gp + Gn + Rp + Rn
fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
mass = np.asarray(cap["mass_global"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
is_exon = (kind == REGION) & (rt == 2)
ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool) & is_exon
err = np.abs(np.asarray(cap["f_g"]) - fo) * mass


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), 1e-12))


# The POST-transport per-side densities the combine actually fused (ag/ap/an from the left, bg/bp/bn right).
ag, ap, an = uni["ag"], uni["ap"], uni["an"]
bg, bp, bn = uni["bg"], uni["bp"], uni["bn"]
# the fused, published gDNA mode
mo_g = uni["mo_g"]
cg = uni["cg"]

print(f"{COND}\n{'=' * 108}")
# (1) is the published mo_g reproduced by the fused density (sanity), and does the pin bind?
pred = lg(cg * E_g / np.maximum(M, _EPS))
print(f"  published mo_g vs log(cg*E_g/M):  max|Δ| = {np.nanmax(np.abs(pred[ok] - mo_g[ok])):.2e}  (sanity)")

# (2) per side: is the delivered composition r-FREE, i.e. equal to the source composition in dst eff-lengths?
for lab, (g_, p_, n_) in (("left", (ag, ap, an)), ("right", (bg, bp, bn))):
    sh_del = g_ * E_g / np.maximum(g_ * E_g + (p_ + n_) * E_r, _EPS)  # delivered gDNA share
    m = ok & (g_ > _EPS) & ((p_ + n_) > _EPS)
    w = mass[m]
    print(f"  {lab:<6} both components supplied: n={int(m.sum()):>4}  "
          f"delivered f_g share vs ORACLE  |Δ|={np.average(np.abs(sh_del[m] - fo[m]), weights=w):.4f}   "
          f"signed={np.average((sh_del - fo)[m], weights=w):+.4f}")

# (3) the pin's binding factor k per stratum — where all components are supplied, r cancels EXACTLY
cpn = uni["cm_p"] + uni["cm_n"]
strat = {
    "both gDNA+RNA supplied": ok & (uni["cm_g"] > _EPS) & (cpn > _EPS),
    "gDNA only": ok & (uni["cm_g"] > _EPS) & (cpn <= _EPS),
    "RNA only": ok & (uni["cm_g"] <= _EPS) & (cpn > _EPS),
}
print(f"\n  {'stratum':<26}{'n':>5}{'ERR':>10}{'share':>8}   -> is `r` cancelled by the pin?")
E = err[ok].sum()
for lab, m in strat.items():
    if not m.any():
        continue
    tag = ("YES — r is irrelevant here" if "both" in lab else
           "NO — the pin fills the missing component from the node's OWN density, so r survives")
    print(f"  {lab:<26}{int(m.sum()):>5}{err[m].sum():>10,.0f}{err[m].sum() / E:>8.1%}   {tag}")

# (4) so how good is the SOURCE composition as a predictor of the destination's? (the imputation premise)
print(f"\n  THE IMPUTATION PREMISE, measured directly (oracle f_g at src vs dst, exon edges):")
li, ri = us["left"], us["right"]
n = len(M)
acc = []
for src_i in (li, ri):
    s = np.clip(src_i, 0, n - 1)
    m = ok & (src_i >= 0) & np.isfinite(fo[s]) & (fo[s] > 1e-6)
    acc.append((np.abs(fo[m] - fo[s][m]), mass[m], np.abs(lg(fo[m]) - lg(fo[s][m]))))
d, w, dl = (np.concatenate([a[i] for a in acc]) for i in range(3))
print(f"    |f_g(dst) − f_g(src)| = {np.average(d, weights=w):.4f}   "
      f"|log f_g(dst) − log f_g(src)| = {np.average(dl, weights=w):.4f}   (n={d.size})")
print("    ^ this is the FLOOR on any message's mode error where the pin cancels r.")
