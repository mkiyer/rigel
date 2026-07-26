"""AMBIG study, step 1 — the census by RAW signature class, and the SIMPLEX-BOUND geometry.

HANDOFF_8 §6. AMBIG nodes have `tau_own = 0` by construction (`node_init`: the strand Beta-Binomial is
rank-1, so the Schur-marginal lambda precision of a 2-DOF node is 0 -- the strand constrains only the TILT).

But that is an INTERIOR argument. The strand mixture is p = 1/2 + (kappa-1/2)*d with d = f_+ - f_-, and the
simplex imposes |d| <= 1 - f_g. When one strand's RNA is ~0 the constraint is ACTIVE, d saturates, and f_g IS
identified: f_g = 1 - |d|. `exon+intron(x-strand)` -- one strand's EXON and the other's INTRON -- is exactly
the geometry where that should happen (and in nrna_none the intron strand's RNA is exactly 0).

This measures, per AMBIG class: the oracle f_g, the oracle MINORITY-strand RNA share, the oracle |d|, the
bound 1-|d|, and how self/solved compare -- i.e. is there identifiable information the model is discarding?

    OMP_NUM_THREADS=1 python scratchpad/ambig_1_census.py
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
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def region_class(sig):
    ep, en = (sig & BIT_EXON_POS) != 0, (sig & BIT_EXON_NEG) != 0
    ip, inn = (sig & BIT_INTRON_POS) != 0, (sig & BIT_INTRON_NEG) != 0
    retained = (ep & ip) | (en & inn)
    ex, it = ep | en, ip | inn
    return np.where(retained, "RETAINED", np.where(ex & it, "exon+intron(x)",
                    np.where(ex, "exon", np.where(it, "intron", "intergenic"))))


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    res = calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    # the ORACLE tilt magnitude |d| = |R+ - R-| / T  (convention-free) and the simplex bound it implies
    d_abs = np.where(T > _EPS, np.abs(Rp - Rn) / np.maximum(T, _EPS), np.nan)
    minor = np.where(R > _EPS, np.minimum(Rp, Rn) / np.maximum(R, _EPS), np.nan)
    fg, mass = np.asarray(cap["f_g"]), np.asarray(cap["mass_global"])
    self_fg, tau = np.asarray(cap["fg_loc"]), np.asarray(cap["_tau0_lam"])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    kind = np.asarray(chain.kind)
    ridx = np.clip(np.asarray(chain.ref_idx, np.int64), 0, len(ra.signature) - 1)
    rcls = region_class(np.asarray(ra.signature)[ridx])
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool)
    err = np.abs(fg - fo)
    em = err * mass
    tot = em[ok].sum()
    print(f"\n{'=' * 122}\n{cond[5:]}   kappa={float(res.rna_sense_frac):.3f}   "
          f"suite mwae {mw(err[ok], mass[ok]):.4f}\n{'=' * 122}")
    print(f"  {'class':<17}{'DOF':<7}{'n':>6}{'err-mass':>11}{'share':>7}{'self':>8}{'solved':>8}"
          f"|{'oracle fg':>10}{'|d|':>7}{'1-|d|':>7}{'minorRNA':>10}{'tau_own':>9}")
    for cl in ("exon", "exon+intron(x)", "RETAINED", "intron", "intergenic"):
        for lab, dm in (("single", ok & (kind == REGION) & ~amb), ("AMBIG", ok & (kind == REGION) & amb)):
            m = dm & (rcls == cl)
            if m.sum() < 3:
                continue
            w = mass[m]
            print(f"  {cl:<17}{lab:<7}{int(m.sum()):>6}{em[m].sum():>11,.0f}{em[m].sum() / tot:>7.1%}"
                  f"{mw(np.abs(self_fg - fo)[m], w):>8.4f}{mw(err[m], w):>8.4f}|{mw(fo[m], w):>10.4f}"
                  f"{mw(d_abs[m], w):>7.3f}{mw(1 - d_abs[m], w):>7.3f}{mw(minor[m], w):>10.4f}"
                  f"{mw(tau[m], w):>9.3f}")
    mb = ok & (kind != REGION)
    for lab, dm in (("single", mb & ~amb), ("AMBIG", mb & amb)):
        if dm.sum() < 3:
            continue
        w = mass[dm]
        print(f"  {'BOUNDARY':<17}{lab:<7}{int(dm.sum()):>6}{em[dm].sum():>11,.0f}{em[dm].sum() / tot:>7.1%}"
              f"{mw(np.abs(self_fg - fo)[dm], w):>8.4f}{mw(err[dm], w):>8.4f}|{mw(fo[dm], w):>10.4f}"
              f"{mw(d_abs[dm], w):>7.3f}{mw(1 - d_abs[dm], w):>7.3f}{mw(minor[dm], w):>10.4f}"
              f"{mw(tau[dm], w):>9.3f}")
    # ── is the simplex bound TIGHT? (1-|d|) - f_g == 0 exactly when the minority strand carries no RNA ──
    a = ok & amb
    slack = (1 - d_abs - fo)[a]
    w = mass[a]
    print(f"\n  SIMPLEX BOUND on AMBIG: mean slack (1-|d|)-f_g = {mw(slack, w):.4f}   "
          f"frac slack<0.05: {mw((slack < 0.05).astype(float), w):.1%}   "
          f"frac minorRNA<0.05: {mw((minor[a] < 0.05).astype(float), w):.1%}")
