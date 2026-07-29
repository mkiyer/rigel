"""IDENTIFIABILITY CENSUS — which nodes can be solved by messages alone, and which need the hyperprior?

The derived rank test (docs/calibration/archive/message_layer_derivation.md §3):

    free DOF of a node = number of ACTIVE RNA strands
      (unknowns rho_g, rho_+, rho_-  minus the HARD observation constraint  sum_c rho_c E_c = M;
       a structurally-absent strand is pinned to 0 and is not an unknown)

        0 active strands  -> 0 DOF   pure gDNA. rho_g = M/E_g. SOLVED BY MEASUREMENT.
        1 active strand   -> 1 DOF   needs ONE independent constraint.
        2 active strands  -> 2 DOF   (AMBIG) needs TWO.

    constraint supply:
        strand likelihood      rank 1   iff kappa != 1/2 beyond the derived noise floor (tau0_lam > 0).
                                        NB for AMBIG it constrains only the product f_r*tau -> rank 1, not 2.
        an incident edge       rank |C_e| - 1  when the edge's log-enrichment step delta_e is UNKNOWN
                                        (the SHIFT limit; C_e = components that structurally transport)
                                        rank |C_e|      when delta_e is KNOWN (capture off / measured gradient)
        adjacent spliced       rank 1   on that strand's RNA axis (a direct measurement, own delta)
        gDNA hyperprior        rank 1   on the gDNA axis  <- pass-2, the fallback

THE SHARP PREDICTION this script tests: an RNA-ISOLATED node (admits RNA, but no incident edge transports any
RNA strand -- topologically, the exon of a SINGLE-EXON transcript) has two gDNA-only edges. |C_e| = 1, so it is
IDENTIFIED when delta is known (capture OFF: delta == 0) and UNIDENTIFIED when delta is unknown (capture ON).
If the derivation is right, that class's corr(f_g, oracle) must COLLAPSE from capture-off to capture-on while
the ordinary classes do not.

Correlation is computed WITHIN each scenario and then pooled (mass-weighted mean), because a cross-scenario
pooled correlation is confounded by the library gDNA level -- every node tracks it, even an uninformed one.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/identifiability_census.py
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
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(
    d.name
    for d in SUITE.iterdir()
    if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name
)

# per (capture, class) -> list of (per-scenario corr, n, mass, mean|err|, mean precision)
acc: dict[tuple, list] = {}


def wcorr(a, b, w):
    if a.size < 4:
        return float("nan")
    ma = np.average(a, weights=w)
    mb = np.average(b, weights=w)
    va = np.average((a - ma) ** 2, weights=w)
    vb = np.average((b - mb) ** 2, weights=w)
    if va < 1e-12 or vb < 1e-12:
        return float("nan")
    return float(np.average((a - ma) * (b - mb), weights=w) / np.sqrt(va * vb))


for cond in conds:
    cap_state = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)  # PASS-0 ONLY (no hyperprior)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)

    fg = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    tau0 = np.asarray(cap["_tau0_lam"])
    spl = np.asarray(cap["spl_l"]) + np.asarray(cap["spl_r"])
    prec = np.asarray(cap["prec_g"]) + np.asarray(cap["prec_p"]) + np.asarray(cap["prec_n"])
    kind = np.asarray(chain.kind)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    nid = np.arange(kind.shape[0])

    n_active = fp.astype(int) + fn.astype(int)  # == free DOF

    def shared(nb):
        ok = nb >= 0
        a = np.clip(nb, 0, None)
        return np.where(ok, (fp[a] & fp[nid]).astype(int) + (fn[a] & fn[nid]).astype(int), -1)

    sh_l, sh_r = shared(left), shared(right)
    rna_edges = (np.maximum(sh_l, 0) + np.maximum(sh_r, 0))
    spl_adj = np.zeros_like(spl)
    for nb in (left, right):
        ok = nb >= 0
        spl_adj = np.maximum(spl_adj, np.where(ok, spl[np.clip(nb, 0, None)], 0.0))
    spl_adj = np.maximum(spl_adj, spl)

    # a BOUNDARY is split by whether it owns a spliced channel: a SPLICE JUNCTION (spliced > 0, so the
    # RNA routing is MEASURED) vs a SEAM (no spliced -- TSS/TES, exon-exon, or a probe-depleted junction).
    kindname = np.where(kind != REGION, np.where(spl > _EPS, "bnd-junction", "bnd-seam"), "reg")
    cls = np.where(
        n_active == 0, "0DOF pure-gDNA",
        np.where(
            (n_active > 0) & (rna_edges == 0), "RNA-ISOLATED",
            np.where(n_active == 1, "1DOF single-strand", "2DOF AMBIG"),
        ),
    )
    cls = np.char.add(np.char.add(cls, " ("), np.char.add(kindname, ")"))

    ok = np.isfinite(fo) & (mass > _EPS)
    for c in np.unique(cls[ok]):
        m = ok & (cls == c)
        if m.sum() < 4:
            continue
        w = mass[m]
        err = np.abs(fg[m] - fo[m])
        acc.setdefault((cap_state, c), []).append(
            (
                wcorr(fg[m], fo[m], w),
                int(m.sum()),
                float(w.sum()),
                float(np.average(err, weights=w)),
                float(np.average(prec[m], weights=w)),
                # HONESTY diagnostic: corr(message precision, |err|). Must be NEGATIVE -- a node the solver
                # is confident about should be a node it got right. Positive => CONFIDENTLY WRONG.
                wcorr(prec[m], err, w),
            )
        )

print(f"\n{len(conds)} UNSTRANDED (ss=0.50) scenarios, PASS-0 only (no hyperprior).")
print("corr is computed WITHIN scenario (mass-weighted), then averaged -- no library-level confound.\n")
hdr = (
    f"{'capture':>7} {'class':>32} | {'nodes':>6} {'mass%':>6} |"
    f" {'corr':>7} {'sd':>6} | {'|err|':>6} | {'prec':>6} | {'r(pr,err)':>9} | verdict"
)
print(hdr)
print("-" * len(hdr))
order = ["0DOF pure-gDNA", "1DOF single-strand", "2DOF AMBIG", "RNA-ISOLATED"]
for capst in ("OFF", "ON ", "VSTR"):
    tot = sum(sum(r[2] for r in v) for k, v in acc.items() if k[0] == capst)
    if tot <= 0:
        continue
    for base in order:
        for k in sorted(acc, key=str):
            if k[0] != capst or not k[1].startswith(base):
                continue
            v = acc[k]
            cr = np.array([r[0] for r in v])
            cr = cr[np.isfinite(cr)]
            hn = np.array([r[5] for r in v])
            hn = hn[np.isfinite(hn)]
            n = sum(r[1] for r in v)
            ms = sum(r[2] for r in v)
            er = np.average([r[3] for r in v], weights=[r[2] for r in v])
            pr = np.average([r[4] for r in v], weights=[r[2] for r in v])
            c = float(cr.mean()) if cr.size else float("nan")
            sd = float(cr.std()) if cr.size > 1 else float("nan")
            h = float(hn.mean()) if hn.size else float("nan")
            note = "" if np.isnan(c) else ("COIN-FLIP" if abs(c) < 0.15 else ("ok" if c > 0.3 else "weak"))
            if not np.isnan(h) and h > 0.05 and not np.isnan(c) and abs(c) < 0.15:
                note += " +CONF-WRONG"
            print(
                f"{capst:>7} {k[1]:>32} | {n:>6} {100*ms/tot:>5.1f}% |"
                f" {c:>7.3f} {sd:>6.3f} | {er:>6.3f} | {pr:>6.2f} | {h:>9.3f} | {note}"
            )
    print()
