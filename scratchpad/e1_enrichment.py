"""ENRICHMENT-RATIO INTERROGATION — the owner's boundary_work_notes.md, claim by claim, against the data.

Five measurements:
  1. do the gDNA and RNA fragment-length distributions actually DIFFER in this suite? (the whole
     "density depends on composition" mechanism is proportional to 1/E_g - 1/E_r)
  2. E_g vs E_r per node class — is the FL differential live where it matters (boundaries)?
  3. the node COUNT distribution — how much mass sits where `1/n` is an over-confident stand-in for the
     Poisson log-variance (psi'(n) is the exact one, and psi'(1)/1 = 1.64)
  4. logvar_tot decomposed into its COUNT half and its COMPOSITION half, per class — which one carries
     Var(log r) in practice
  5. the per-fragment FL-MIXTURE Fisher information about f_g, against the strand likelihood's — is the
     fragment length an unused FOURTH information source?

    OMP_NUM_THREADS=1 python scratchpad/e1_enrichment.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    pg = np.asarray(inp["gdna_fl_pmf"], float)
    pr = np.asarray(inp["rna_fl_pmf"], float)
    L = max(pg.size, pr.size)
    pg = np.pad(pg, (0, L - pg.size))
    pr = np.pad(pr, (0, L - pr.size))
    pg = pg / max(pg.sum(), _EPS)
    pr = pr / max(pr.sum(), _EPS)
    ell = np.arange(L, dtype=float)
    mg, mr = float((pg * ell).sum()), float((pr * ell).sum())
    tv = 0.5 * float(np.abs(pg - pr).sum())

    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    isr = kind == REGION
    n = kind.shape[0]
    M, E_g, E_r, lv = us["M"], us["E_g"], us["E_r"], us["logvar_tot"]
    n_node = np.where(isr, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    live = (M > _EPS) & (n_node > 0)

    print(f"\n{'=' * 112}\n{cond[5:]}\n{'=' * 112}")
    print(f"  [1] FRAGMENT-LENGTH DISTRIBUTIONS: gDNA mean {mg:.1f} bp, RNA mean {mr:.1f} bp, "
          f"ratio {mr / max(mg, _EPS):.3f}, total-variation distance {tv:.4f}")

    # [5] per-fragment FL-mixture Fisher information about f_g, at f_g = 1/2, vs the strand's
    mix = 0.5 * pg + 0.5 * pr
    ok = mix > 1e-12
    fisher_fl = float(np.sum((pg[ok] - pr[ok]) ** 2 / mix[ok]))
    # the strand's per-fragment Fisher about f_g, for scale: p = kappa + f_g*(1/2-kappa),
    # I = (dp/df_g)^2 / (p(1-p)) = (1/2-kappa)^2 / (p(1-p)).  Identically 0 at kappa = 1/2.
    def _i_strand(kappa, f_g=0.5):
        p = kappa + f_g * (0.5 - kappa)
        return (0.5 - kappa) ** 2 / (p * (1.0 - p))
    print(f"      per-fragment FL-mixture Fisher information about f_g (at f_g=1/2): I = {fisher_fl:.4f}"
          f"   [for scale, the STRAND: I = {_i_strand(0.99):.4f} at kappa=0.99, "
          f"{_i_strand(0.5):.4f} at kappa=0.50]")

    print(f"\n  [2/3/4] per NODE CLASS ({'mass-weighted' if True else ''}):")
    print(f"  {'class':<12}{'n':>6}{'reads':>11}{'E_g':>9}{'E_r':>9}{'E_g/E_r':>9}"
          f"{'|coef|':>9}{'med n':>8}{'n<=5':>7}{'1/n':>9}{'comp':>9}{'psi/1n':>8}")
    for c in ("exon", "intron", "intergenic", "boundary"):
        m = live & (cls == c)
        if m.sum() < 3:
            continue
        w = M[m]
        eg, er = E_g[m], E_r[m]
        coef = np.abs((1.0 / np.maximum(eg, _EPS) - 1.0 / np.maximum(er, _EPS))
                      / np.maximum(0.5 / np.maximum(eg, _EPS) + 0.5 / np.maximum(er, _EPS), _EPS))
        nn = n_node[m]
        inv_n = 1.0 / np.maximum(nn, _EPS)
        comp = np.maximum(lv[m] - inv_n, 0.0)
        # the exact Poisson log-variance vs the shipped 1/n stand-in
        ratio = polygamma(1, np.maximum(nn, 1.0)) / np.maximum(inv_n, _EPS)
        print(f"  {c:<12}{int(m.sum()):>6}{w.sum():>11,.0f}"
              f"{np.average(eg, weights=w):>9.1f}{np.average(er, weights=w):>9.1f}"
              f"{np.average(eg / np.maximum(er, _EPS), weights=w):>9.3f}"
              f"{np.average(coef, weights=w):>9.3f}"
              f"{np.median(nn):>8.0f}{np.mean(nn <= 5):>7.1%}"
              f"{np.average(inv_n, weights=w):>9.4f}{np.average(comp, weights=w):>9.4f}"
              f"{np.average(ratio, weights=w):>8.3f}")

    # how much of the suite's node MASS sits on low-count nodes (where 1/n is the worst stand-in)
    for thr in (1, 2, 5, 20):
        sel = live & (n_node <= thr)
        print(f"      nodes with n <= {thr:<3d}: {sel.sum():>5} nodes, "
              f"{M[sel].sum() / M[live].sum():>6.2%} of node mass, "
              f"{np.mean(cls[sel] == 'intron') if sel.sum() else 0:>5.0%} of them introns")
