"""TARGET gdna300_ss_0.50_nrna_none_capture_on — step 1: characterize the 1.36 M-read error.

Why this scenario (pass0_error_table.py): largest single error scenario in the suite (1,358,610 reads,
10.8 % of all suite error, mwae 0.1870), single-strand carries 75 % of it, and its stated precision is HONEST
(4.1 % of error on the most-confident quartile) so nothing here is hidden behind false confidence.

Its physics: kappa = 1/2 ⇒ the strand likelihood is FLAT ⇒ tau_own = 0 at EVERY node except introns (the
factory). No nascent RNA ⇒ every intron's true RNA is exactly 0. Capture ON ⇒ large enrichment cliffs. So the
intron factory + message propagation are the ONLY information in the whole solve.

    OMP_NUM_THREADS=1 python scratchpad/t1_char.py
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
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_none_capture_on"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
res = calmod.calibrate(
    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
    np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cc, calib_refit_iters=0), _debug=dbg,
)
chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
uni = cap["_uni"][-1]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
T = G + R
fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
amb = fp & fn
mass = np.asarray(cap["mass_global"])
ship = np.asarray(cap["f_g"])
self_fg = np.asarray(cap["fg_loc"])
solvable = np.asarray(cap["solvable"], bool)
kind = np.asarray(chain.kind)
ridx = np.clip(np.asarray(chain.ref_idx, np.int64), 0, len(ra.signature) - 1)
sig = np.asarray(ra.signature)[ridx]
ep, en = (sig & BIT_EXON_POS) != 0, (sig & BIT_EXON_NEG) != 0
ip, inn = (sig & BIT_INTRON_POS) != 0, (sig & BIT_INTRON_NEG) != 0
ret = (ep & ip) | (en & inn)
ex, it = ep | en, ip | inn
rcls = np.where(kind != REGION, "boundary",
                np.where(ret, "RETAINED", np.where(ex & it, "exon+intron(x)",
                         np.where(ex, "exon", np.where(it, "intron", "intergenic")))))
ok = np.isfinite(fo) & (mass > _EPS) & solvable
err = np.abs(ship - fo) * mass
serr = np.abs(self_fg - fo) * mass


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


# ── the psi replay + channel ablation (bit-exact) ─────────────────────────────────────────────────────────
cp, cn = uni["cp"], uni["cn"]
cR = cp + cn
th_msg = np.arcsin(np.clip(np.where(cR > _EPS, (cp - cn) / np.maximum(cR, _EPS), 0.0), -1, 1))
th_prec = np.where(amb, uni["cm_p"] + uni["cm_n"], 0.0)
Z = np.zeros_like(uni["cm_g"])


def psi(*, lam=True, gdna=True, rna=True, tilt=True, factory=True):
    return _solve_nodes_logodds_all(
        st.u_pos, st.u_neg, fp, fn, st.mass_unspliced, st.mass_spliced,
        kappa=float(res.rna_sense_frac), od_g=float(res.gdna_strand_overdispersion),
        od_r=float(res.rna_strand_overdispersion),
        n_grid=int(cc.sweep_n_grid), L=float(cc.sweep_logodds_window),
        n_tilt=cc.sweep_n_tilt, n_grid_ss=cc.sweep_n_grid_single_strand,
        global_logprior=cap["global_lp"],
        gdna_imp_mode=uni["mo_g"], gdna_imp_prec=uni["cm_g"] if gdna else Z,
        rna_imp_mode=(uni["mo_p"], uni["mo_n"]),
        rna_imp_prec=(uni["cm_p"], uni["cm_n"]) if rna else (Z, Z),
        lam_logprior=cap["intron_prior"] if factory else None,
        lam_imp_mode=uni["lam_msg"], lam_imp_prec=uni["c_tau"] if lam else Z,
        theta_imp_mode=th_msg, theta_imp_prec=th_prec if tilt else Z,
        fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
    ).gdna_frac


base = np.where(solvable, psi(), ship)
print(f"{COND}\n  replay fidelity max|Δ| = {np.max(np.abs(base - ship)):.2e}   "
      f"kappa={float(res.rna_sense_frac):.5f}  mass={mass[ok].sum():,.0f}  "
      f"ERR={err[ok].sum():,.0f}  selfERR={serr[ok].sum():,.0f}")

arms = {"-lam": psi(lam=False), "-gdna": psi(gdna=False), "-rna": psi(rna=False),
        "-tilt": psi(tilt=False), "-factory": psi(factory=False),
        "nomsg": psi(lam=False, gdna=False, rna=False, tilt=False)}
print(f"\n  ABLATION (error reads, whole scenario):  shipped {err[ok].sum():>10,.0f}")
for k, v in arms.items():
    e = (np.abs(np.where(solvable, v, ship) - fo) * mass)[ok].sum()
    print(f"     {k:<10}{e:>12,.0f}   ({e / err[ok].sum() - 1:+.1%})")

# ── where is it, by class × DOF ───────────────────────────────────────────────────────────────────────────
print(f"\n  {'class':<17}{'DOF':<7}{'n':>6}{'reads':>12}{'ERR':>11}{'share':>7}"
      f"{'mwae':>8}{'self':>8}{'signed':>9}|{'tau_own':>9}{'c_tau':>8}{'cm_g':>8}{'cm_p+n':>9}")
for cl in ("exon", "exon+intron(x)", "RETAINED", "intron", "intergenic", "boundary"):
    for lab, dm in (("single", ok & ~amb), ("AMBIG", ok & amb)):
        m = dm & (rcls == cl)
        if m.sum() < 3:
            continue
        w = mass[m]
        print(f"  {cl:<17}{lab:<7}{int(m.sum()):>6}{w.sum():>12,.0f}{err[m].sum():>11,.0f}"
              f"{err[m].sum() / err[ok].sum():>7.1%}{mw(np.abs(ship - fo)[m], w):>8.4f}"
              f"{mw(np.abs(self_fg - fo)[m], w):>8.4f}{mw((ship - fo)[m], w):>+9.4f}|"
              f"{mw(np.asarray(cap['_tau0_lam'])[m], w):>9.3f}{mw(uni['c_tau'][m], w):>8.2f}"
              f"{mw(uni['cm_g'][m], w):>8.2f}{mw((uni['cm_p'] + uni['cm_n'])[m], w):>9.2f}")

# ── the worst individual nodes ────────────────────────────────────────────────────────────────────────────
print(f"\n  WORST 20 NODES by error reads\n  {'node':>6}{'cls':>16}{'dof':>7}{'reads':>10}"
      f"{'orc':>6}{'self':>6}{'ship':>6}|{'ERR':>9}|{'c_tau':>8}{'lamfg':>7}{'cm_g':>8}{'mo_g':>7}"
      f"{'cm_p':>8}{'mo_p':>7}{'cm_n':>8}{'mo_n':>7}|{'nL':>14}{'nR':>14}")
left, right = np.asarray(chain.left), np.asarray(chain.right)
for i in np.argsort(-np.where(ok, err, 0))[:20]:
    lo, ro = int(left[i]), int(right[i])
    print(f"  {i:>6}{rcls[i]:>16}{'AMBIG' if amb[i] else 'single':>7}{mass[i]:>10,.0f}"
          f"{fo[i]:>6.3f}{self_fg[i]:>6.3f}{ship[i]:>6.3f}|{err[i]:>9,.0f}|"
          f"{uni['c_tau'][i]:>8.2f}{1 / (1 + np.exp(-uni['lam_msg'][i])):>7.3f}"
          f"{uni['cm_g'][i]:>8.2f}{min(np.exp(uni['mo_g'][i]), 9.99):>7.3f}"
          f"{uni['cm_p'][i]:>8.2f}{min(np.exp(uni['mo_p'][i]), 9.99):>7.3f}"
          f"{uni['cm_n'][i]:>8.2f}{min(np.exp(uni['mo_n'][i]), 9.99):>7.3f}|"
          f"{(rcls[lo] if lo >= 0 else '-'):>14}{(rcls[ro] if ro >= 0 else '-'):>14}")
