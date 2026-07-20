"""TRAINING-NODE ANALYSIS — after the Phase-1a prior-free first pass, WHICH nodes are accurately solved, and
does the solver's own PRECISION predict that accuracy? This is the design substrate for the gDNA hyperprior
(Goal 2): we may only train the deconvolved-gDNA NPMLE on nodes we can trust, so we must know which those are.

Drives the REAL production `calibrate()` (Phase 1a), then per node compares the solved f_g to the ORACLE and to
the solver's precision (var_gdna = Var(log f_g); LOW var = HIGH precision/confidence). Answers:

  1. by node CLASS (gonly=structural-gDNA / single-strand / AMBIG / boundary) and signature: accuracy + precision.
  2. does precision predict accuracy? — |err| across precision deciles + the rank correlation.
  3. HIGH-PRECISION-WRONG nodes (confident but wrong) — count, mass, where. The danger for the fit.
  4. candidate TRAINING SETS — the mass-weighted accuracy of each (gonly, single, +background, boundaries, all).

    OMP_NUM_THREADS=1 python training_node_analysis.py [--condition C] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node, _solve  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG  # noqa: E402

_EPS = 1e-12
_EXON = BIT_EXON_POS | BIT_EXON_NEG
_INTRON = BIT_INTRON_POS | BIT_INTRON_NEG


def _sig_type(sig):
    """Coarse region signature class from the 4-bit strand/type signature."""
    if sig == 0:
        return "intergenic"
    if sig & _EXON:
        return "exon"
    if sig & _INTRON:
        return "intron"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--condition", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    inp = _scan_and_truth(suite, a.condition, index, cfg, work, cache)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    dbg = _solve(inp, ra, 0)  # Phase 1a (the only solve; refit deferred)
    chain = dbg["chain"]
    cap = dbg["capture"]
    belief = dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fg = np.asarray(belief.f_g, float)
    var_g = np.asarray(belief.var_gdna, float)  # Var(log f_g): LOW = confident
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    kind = np.asarray(chain.kind)
    ref_idx = np.asarray(chain.ref_idx, np.int64)

    tot = G + R
    fo = np.where(tot > _EPS, G / np.maximum(tot, _EPS), np.nan)
    live = (eff > 1e-9 * 1.001) & (mass > _EPS) & np.isfinite(fo) & (tot > _EPS)
    err = np.where(live, np.abs(fg - fo), np.nan)
    cls = np.where(fp & fn, "AMBIG", np.where(fp | fn, "single", "gonly"))
    cls = np.where(~isr, "bndry", cls)
    # region signature type (regions only)
    sig_arr = index.region_df["signature"].to_numpy()
    sigt = np.array([
        _sig_type(int(sig_arr[ref_idx[i]])) if (kind[i] == REGION and ref_idx[i] < sig_arr.shape[0]) else "bndry"
        for i in range(fg.shape[0])
    ])

    def _wmae(mask):
        w = np.where(mask & live, mass, 0.0)
        if w.sum() <= 0:
            return float("nan"), 0.0, 0
        e = float((w * np.where(mask & live, np.abs(fg - fo), 0.0)).sum() / w.sum())
        return e, float(w.sum() / np.where(live, mass, 0.0).sum()), int((mask & live).sum())

    print(f"=== {a.condition}   [Phase-1a first pass]  live={int(live.sum())}  "
          f"true gDNAfrac={G[live].sum()/max(tot[live].sum(),1):.3f} ===\n")

    print("(1) BY CLASS — accuracy + precision (var_gdna: LOW=confident):")
    print(f"{'class':>8} {'n':>6} {'massfr':>7} {'wmae':>6} {'medVar':>7} {'p10Var':>7} {'trueG':>6} {'trueR':>6}")
    for name in ("gonly", "single", "AMBIG", "bndry"):
        m = (cls == name) & live
        e, mf, n = _wmae(cls == name)
        v = var_g[m]
        mv = float(np.median(v)) if v.size else float("nan")
        p10 = float(np.percentile(v, 10)) if v.size else float("nan")
        w = np.where(m, mass, 0.0)
        tg = float(np.where(m & (fo > 0.8), mass, 0).sum() / max(w.sum(), _EPS))
        tr = float(np.where(m & (fo < 0.2), mass, 0).sum() / max(w.sum(), _EPS))
        print(f"{name:>8} {n:>6} {mf:>7.3f} {e:>6.3f} {mv:>7.3f} {p10:>7.3f} {tg:>6.2f} {tr:>6.2f}")

    print("\n(1b) BY REGION SIGNATURE (regions only):")
    for name in ("intergenic", "intron", "exon"):
        e, mf, n = _wmae(sigt == name)
        m = (sigt == name) & live
        mv = float(np.median(var_g[m])) if m.any() else float("nan")
        print(f"  {name:>10}: n={n:>5} massfr={mf:.3f} wmae={e:.3f} medVar={mv:.3f}")

    # (2) precision → accuracy: deciles of var_g over live nodes
    print("\n(2) PRECISION → ACCURACY (live nodes binned by var_gdna decile; low decile = most confident):")
    lv = var_g[live]
    le = np.abs(fg[live] - fo[live])
    lm = mass[live]
    q = np.quantile(lv, np.linspace(0, 1, 11))
    print(f"{'decile':>7} {'varRange':>16} {'wmae':>6} {'massfr':>7} {'n':>6}")
    for d in range(10):
        sel = (lv >= q[d]) & (lv <= q[d + 1]) if d == 9 else (lv >= q[d]) & (lv < q[d + 1])
        w = lm[sel]
        e = float((w * le[sel]).sum() / max(w.sum(), _EPS)) if w.sum() > 0 else float("nan")
        print(f"{d:>7} [{q[d]:>6.3f},{q[d+1]:>6.3f}] {e:>6.3f} {w.sum()/lm.sum():>7.3f} {int(sel.sum()):>6}")
    # rank correlation precision(=−var) vs accuracy(=−err): Spearman via rank
    from scipy.stats import spearmanr
    rho, _ = spearmanr(lv, le)
    print(f"  Spearman(var_gdna, |err|) = {rho:+.3f}   (>0 ⇒ higher var/lower precision ⇒ higher error = GOOD)")

    # (3) high-precision-wrong: most-confident 30% (lowest var) that are still badly wrong
    thr = np.quantile(lv, 0.30)
    conf = live & (var_g <= thr)
    hpw = conf & (np.abs(fg - fo) > 0.3)
    w_conf = np.where(conf, mass, 0.0).sum()
    w_hpw = np.where(hpw, mass, 0.0).sum()
    print(f"\n(3) HIGH-PRECISION-WRONG (var ≤ p30={thr:.3f} AND |err|>0.3): "
          f"n={int(hpw.sum())} massfrac(of confident)={w_hpw/max(w_conf,_EPS):.3f}")
    for name in ("gonly", "single", "AMBIG", "bndry"):
        m = hpw & (cls == name)
        if m.any():
            print(f"    {name:>7}: n={int(m.sum())} mass={np.where(m,mass,0).sum():.0f} "
                  f"signed(fg−fo)={np.where(m,fg-fo,0).sum()/max(m.sum(),1):+.2f}")

    # (4) candidate training sets: mass-weighted accuracy of each (how reliable as gDNA training data)
    print("\n(4) CANDIDATE TRAINING SETS (mass-wt |f_g − oracle| — lower = more reliable to fit the gDNA NPMLE on):")
    sets = {
        "gonly (structural gDNA)": cls == "gonly",
        "intergenic+intron (background)": np.isin(sigt, ["intergenic", "intron"]),
        "single-strand": cls == "single",
        "single + gonly": np.isin(cls, ["single", "gonly"]),
        "boundaries": cls == "bndry",
        "ALL live": np.ones_like(live),
    }
    for name, m in sets.items():
        e, mf, n = _wmae(m)
        print(f"  {name:<32} wmae={e:.3f}  massfr={mf:.3f}  n={n}")

    # ---- plot: var_gdna vs |err|, by class ----
    outdir = Path(a.out) if a.out else Path(os.environ.get("RIGEL_SCRATCH", "/tmp"))
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = {"gonly": "tab:green", "single": "tab:blue", "AMBIG": "tab:red", "bndry": "0.5"}
    for name, c in colors.items():
        m = (cls == name) & live
        ax.scatter(var_g[m], np.abs(fg[m] - fo[m]), s=np.sqrt(np.maximum(mass[m], 1)) * 0.3,
                   c=c, alpha=0.4, label=f"{name} (n={int(m.sum())})")
    ax.set_xlabel("var_gdna = Var(log f_g)   (← more confident        less confident →)")
    ax.set_ylabel("|f_g solved − oracle|")
    ax.set_title(f"Phase-1a: precision vs accuracy per node\n{a.condition}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = outdir / "training_node_analysis.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
