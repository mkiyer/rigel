"""P1 — root-cause the CONFIDENTLY-WRONG single-strand exons on unstranded x capture-OFF.

378,815 reads (25% of all confidently-wrong mass) with z2|Q1 = 33, against 0-6 in every other regime.

HYPOTHESIS: capture-OFF removes the LAST damping term. On unstranded data tau_own = 0, so DL is inert
(v_own = inf => b_hat^2 = 0, by design). Capture-OFF then makes r ~ 1, so M5's sigma^2_transfer ~ 0 too.
Nothing is left to damp the messages, so they arrive at full precision and pin the node.

The test is PAIRED: capture-off and capture-on are the same genome and the same region partition, so the
same node index is the same piece of DNA. Compare the channel precisions node-by-node.

    OMP_NUM_THREADS=1 python scratchpad/p1_capoff.py
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
PAIR = [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",  # the target: z2|Q1 = 33
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",  # the control: z2|Q1 = 6
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

S = {}
for cond in PAIR:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    kind = np.asarray(chain.kind)
    rt, _ = _node_region_type(chain, ra)
    n = kind.shape[0]
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    S[cond] = dict(
        fo=fo, fg=np.asarray(cap["f_g"]), var=np.asarray(cap["var_g"]),
        mass=np.asarray(cap["mass_global"]), cls=cls, amb=fp & fn,
        solv=np.asarray(cap["solvable"], bool), lv=us["logvar_tot"], tau_own=us["tau_own"],
        c_tau=uni["c_tau"], cm_g=uni["cm_g"], cm_p=uni["cm_p"], cm_n=uni["cm_n"],
        cpg=uni["cpg"], fg_loc=np.asarray(cap["fg_loc"]),
    )

off, on = (S[c] for c in PAIR)
ok = (
    np.isfinite(off["fo"]) & np.isfinite(on["fo"]) & (off["mass"] > _EPS)
    & off["solv"] & (off["cls"] == "exon") & ~off["amb"]
)
err_off = np.abs(off["fg"] - off["fo"])
err_on = np.abs(on["fg"] - on["fo"])
sd_off, sd_on = np.sqrt(np.maximum(off["var"], 0)), np.sqrt(np.maximum(on["var"], 0))
# the SUITE-WIDE confident-quartile threshold (from pass0_error_table.py's saved state), so the population
# here is exactly the one the CWRONG census counts -- not a scenario-local quantile.
_st = np.load("/tmp/pass0_state.npz", allow_pickle=True)
_v = _st["var"]
q1 = float(np.quantile(_v[np.isfinite(_v)], 0.25))
conf = ok & np.isfinite(off["var"]) & (off["var"] <= q1)

print(f"{'=' * 116}\nP1 — single-strand EXONS, unstranded, capture-OFF vs capture-ON (paired: same nodes)"
      f"\n{'=' * 116}")
print(f"  {'stratum':<34}{'n':>6}{'reads':>12}{'|err|':>8}{'sd':>8}{'z2':>8}"
      f"{'c_tau':>9}{'cm_g':>9}{'cm_R':>9}{'sd(log r)':>11}")
for lab, arm, m in (
    ("capture OFF  all exon-single", off, ok), ("capture OFF  confident quartile", off, conf),
    ("capture ON   all exon-single", on, ok), ("capture ON   same nodes as above", on, conf),
):
    e = np.abs(arm["fg"] - arm["fo"])
    w = arm["mass"]
    z2 = float(np.sum(w[m] * e[m] ** 2) / max(np.sum(w[m] * arm["var"][m]), _EPS))
    print(f"  {lab:<34}{int(m.sum()):>6}{w[m].sum():>12,.0f}"
          f"{np.average(e[m], weights=w[m]):>8.4f}"
          f"{np.average(np.sqrt(np.maximum(arm['var'][m], 0)), weights=w[m]):>8.4f}{z2:>8.1f}"
          f"{np.average(arm['c_tau'][m], weights=w[m]):>9.2f}"
          f"{np.average(arm['cm_g'][m], weights=w[m]):>9.2f}"
          f"{np.average(arm['cm_p'][m] + arm['cm_n'][m], weights=w[m]):>9.2f}"
          f"{np.median(np.sqrt(np.maximum(arm['lv'][m], 0))):>11.4f}")

print(f"\n  the CONFIDENT + WRONG exons on capture-OFF ({int(conf.sum())} nodes, "
      f"{off['mass'][conf].sum():,.0f} reads):")
print(f"  {'':<4}{'oracle':>8}{'self':>8}{'solved':>8}{'|err|':>8}{'sd':>8}{'z':>7}|"
      f"{'c_tau':>8}{'cm_g':>8}{'cm_R':>8}{'tau_own':>9}|{'reads':>10}")
idx = np.flatnonzero(conf)
idx = idx[np.argsort(-(err_off[idx] * off["mass"][idx]))][:12]
for i in idx:
    z = err_off[i] / max(sd_off[i], 1e-9)
    print(f"  {i:<4}{off['fo'][i]:>8.4f}{off['fg_loc'][i]:>8.4f}{off['fg'][i]:>8.4f}"
          f"{err_off[i]:>8.4f}{sd_off[i]:>8.4f}{z:>7.1f}|"
          f"{off['c_tau'][i]:>8.2f}{off['cm_g'][i]:>8.2f}{off['cm_p'][i] + off['cm_n'][i]:>8.2f}"
          f"{off['tau_own'][i]:>9.3f}|{off['mass'][i]:>10,.0f}")

print("\n  WHERE THE POSTERIOR PRECISION COMES FROM (mass-weighted, confident+wrong exons):")
for lab, arm, m in (("capture OFF", off, conf), ("capture ON ", on, conf)):
    w = arm["mass"][m]
    tot = arm["c_tau"][m] + arm["cm_g"][m] + arm["cm_p"][m] + arm["cm_n"][m]
    print(f"    {lab}: c_tau {np.average(arm['c_tau'][m], weights=w):>8.2f}"
          f" ({np.average(arm['c_tau'][m] / np.maximum(tot, _EPS), weights=w):>5.1%})"
          f"   cm_g {np.average(arm['cm_g'][m], weights=w):>8.2f}"
          f" ({np.average(arm['cm_g'][m] / np.maximum(tot, _EPS), weights=w):>5.1%})"
          f"   cm_R {np.average(arm['cm_p'][m] + arm['cm_n'][m], weights=w):>8.2f}"
          f" ({np.average((arm['cm_p'][m] + arm['cm_n'][m]) / np.maximum(tot, _EPS), weights=w):>5.1%})"
          f"   tau_own {np.average(arm['tau_own'][m], weights=w):>7.3f}")
