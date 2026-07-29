"""Boundary-rule dissection — WHAT distinguishes the MEANINGFUL boundaries (corr 0.68 with the oracle) from the
COIN-FLIP ones (corr 0.13)? The §6B DOF gate inverts on boundaries (`docs/calibration/archive/solve_gate_design.md` §4): the boundaries it
calls SOLVABLE are coin-flips, those it WITHHOLDS are meaningful. The owner's hypothesis: it is the ENRICHMENT
CLIFF — an intron|exon boundary straddles it and its gDNA estimate is muddled by discordant flank messages (the
depleted-intron message vs the enriched-exon message), NOT a technical DOF/solvability issue.

This pools the UNSTRANDED (ss0.50) ambig_dense_10mb scenarios and, for every boundary node, computes:
  * the flank-type pair (intron|exon / exon|exon / …) — the owner's primary cut;
  * the enrichment gap |mu_proj[L] − mu_proj[R]| (the cliff magnitude across the boundary) + the per-flank
    σ²_transfer that damps each message;
  * the per-flank gDNA message (mode, precision) — the forward scan carries the LEFT-flank message, the backward
    scan the RIGHT-flank; their disagreement is the "muddle" signal;
  * one-sided-spliced presence; the boundary's own strand info (tau0_lam); nfree.
Then it stratifies corr(f_g, oracle) by each. Error is NOT the metric (unsolved f_g is an arbitrary default) —
CORRELATION is (does the forced solve track the oracle?)."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and "0.50" in d.name)

_TNAME = {-1: "term", 0: "intgc", 1: "intron", 2: "exon"}

# per-boundary pooled columns
cols = {k: [] for k in (
    "nid", "fg", "fo", "mass", "ltyp", "rtyp", "gap", "vproj", "s2t_l", "s2t_r",
    "mgl", "pgl", "mgr", "pgr", "t0", "nfree", "spl", "dof")}


def _flank_type(node_type, chain, i, side_arr):
    j = int(side_arr[i])
    return int(node_type[j]) if j >= 0 else -1


for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    node_type, _ = _node_region_type(chain, ra)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"]); mass = np.asarray(cap["mass_global"])
    fp = np.asarray(cap["free_pos"], bool); fn = np.asarray(cap["free_neg"], bool)
    pg, pp, pn = np.asarray(cap["prec_g"]), np.asarray(cap["prec_p"]), np.asarray(cap["prec_n"])
    t0 = np.asarray(cap["_tau0_lam"]); kind = np.asarray(chain.kind)
    # the LIVE enrichment scale + its variance (the retired NPMLE projection's replacement)
    _us = cap["_uni_static"]
    mu_proj = np.log(np.maximum(np.asarray(_us["rho_node0"]), 1e-12))
    var_proj = np.asarray(_us["logvar_tot"])
    a_fwd, b_bwd = cap["a_fwd"], cap["b_bwd"]  # (amg,apg,amp,app,amn,apn) each; fwd=left msg, bwd=right msg
    mgl, pgl = np.asarray(a_fwd[0]), np.asarray(a_fwd[1])  # gDNA message from the LEFT flank
    mgr, pgr = np.asarray(b_bwd[0]), np.asarray(b_bwd[1])  # gDNA message from the RIGHT flank
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    spl = np.asarray(cap["spl_l"]) + np.asarray(cap["spl_r"])
    lam_id = (t0 > _EPS) | (pg > _EPS) | (pp + pn > _EPS)
    th_id = (t0 > _EPS) | (pp > _EPS) | (pn > _EPS)
    nfree = fp.astype(int) + fn.astype(int)
    dof_solv = np.where(nfree >= 2, lam_id & th_id, np.where(nfree == 1, lam_id, False))

    bnd = (kind != REGION) & np.isfinite(fo) & (mass > _EPS)
    for i in np.where(bnd)[0]:
        lt = _flank_type(node_type, chain, i, left)
        rt = _flank_type(node_type, chain, i, right)
        gap_l = mu_proj[i] - (mu_proj[left[i]] if left[i] >= 0 else mu_proj[i])
        gap_r = mu_proj[i] - (mu_proj[right[i]] if right[i] >= 0 else mu_proj[i])
        flank_gap = abs((mu_proj[left[i]] if left[i] >= 0 else mu_proj[i]) -
                        (mu_proj[right[i]] if right[i] >= 0 else mu_proj[i]))
        cols["nid"].append(int(i))
        cols["fg"].append(fg[i]); cols["fo"].append(fo[i]); cols["mass"].append(mass[i])
        cols["ltyp"].append(lt); cols["rtyp"].append(rt)
        cols["gap"].append(flank_gap); cols["vproj"].append(var_proj[i])
        cols["s2t_l"].append(var_proj[i] + gap_l * gap_l)
        cols["s2t_r"].append(var_proj[i] + gap_r * gap_r)
        cols["mgl"].append(mgl[i]); cols["pgl"].append(pgl[i])
        cols["mgr"].append(mgr[i]); cols["pgr"].append(pgr[i])
        cols["t0"].append(t0[i]); cols["nfree"].append(int(nfree[i]))
        cols["spl"].append(spl[i]); cols["dof"].append(bool(dof_solv[i]))

C = {k: np.array(v) for k, v in cols.items()}
tot_mass = C["mass"].sum()
print(f"pooled {len(conds)} unstranded scenarios — {len(C['fg'])} boundary node-observations\n")
print("WARNING: pooled corr is confounded by between-structure signal (Simpson) — see the per-position\n"
      "         across-scenario metric at the bottom, which is the honest one.\n")


def report(title, groups):
    """groups = list of (label, boolean mask). Print corr(fg,oracle) + n + mass% + mean|err| per group."""
    print(f"── {title} " + "─" * max(2, 58 - len(title)))
    print(f"{'group':>22} | {'n':>5} {'mass%':>6} | {'corr':>7} | {'mean|fg-or|':>11} | note")
    for label, m in groups:
        if m.sum() < 3:
            print(f"{label:>22} | {int(m.sum()):>5} {'':>6} | {'—':>7} |")
            continue
        fgg, fog, mg = C["fg"][m], C["fo"][m], C["mass"][m]
        corr = np.corrcoef(fgg, fog)[0, 1] if fgg.std() > 1e-9 and fog.std() > 1e-9 else float("nan")
        err = float(np.average(np.abs(fgg - fog), weights=mg))
        note = ("COIN-FLIP" if (not np.isnan(corr) and abs(corr) < 0.20)
                else "meaningful" if (not np.isnan(corr) and corr > 0.35) else "")
        print(f"{label:>22} | {int(m.sum()):>5} {100*mg.sum()/tot_mass:>5.1f}% | {corr:>7.3f} | "
              f"{err:>11.3f} | {note}")
    print()


# ── CUT 1: the DOF verdict (reproduce the inversion, as the anchor) ──
report("CUT 1 — DOF verdict", [
    ("DOF solvable", C["dof"]),
    ("DOF unsolvable", ~C["dof"]),
])

# ── CUT 2: flank-type pair (the owner's PRIMARY hypothesis) ──
def _pair(a, b):
    x, y = sorted((a, b))
    return f"{_TNAME[x]}|{_TNAME[y]}"
pair = np.array([_pair(int(lf), int(rf)) for lf, rf in zip(C["ltyp"], C["rtyp"])])
upairs = sorted(set(pair), key=lambda p: -np.count_nonzero(pair == p))
report("CUT 2 — flank-type pair", [(p, pair == p) for p in upairs])

# ── CUT 3: enrichment gap between flanks (the cliff magnitude) ──
gap = C["gap"]
qs = np.quantile(gap, [0.0, 0.25, 0.5, 0.75, 1.0])
report("CUT 3 — flank enrichment gap (log-density; quartiles)", [
    (f"gap≤{qs[1]:.2f} (concordant)", gap <= qs[1]),
    (f"{qs[1]:.2f}<gap≤{qs[2]:.2f}", (gap > qs[1]) & (gap <= qs[2])),
    (f"{qs[2]:.2f}<gap≤{qs[3]:.2f}", (gap > qs[2]) & (gap <= qs[3])),
    (f"gap>{qs[3]:.2f} (cliff)", gap > qs[3]),
])

# ── CUT 4: how many flanks actually send a gDNA message (both / one / neither) ──
has_l = C["pgl"] > _EPS
has_r = C["pgr"] > _EPS
report("CUT 4 — flanks sending a gDNA message", [
    ("both flanks", has_l & has_r),
    ("one flank", has_l ^ has_r),
    ("neither flank", ~has_l & ~has_r),
])

# ── CUT 5: among two-message boundaries, do the flank gDNA messages AGREE? ──
both = has_l & has_r
disagree = np.abs(C["mgl"] - C["mgr"])  # |log f_g^left − log f_g^right| (only meaningful where both>0)
if both.sum() > 10:
    dq = np.quantile(disagree[both], [0.5])
    report("CUT 5 — flank gDNA-message disagreement (two-message boundaries)", [
        ("agree  (Δlog f_g ≤ med)", both & (disagree <= dq[0])),
        ("disagree (Δlog f_g > med)", both & (disagree > dq[0])),
    ])

# ── CUT 6: one-sided spliced (mature) present at the boundary? ──
report("CUT 6 — one-sided spliced present", [
    ("has spliced", C["spl"] > _EPS),
    ("no spliced", C["spl"] <= _EPS),
])

# ── CUT 7: the confound — DOF verdict × flank-type (is DOF just a proxy for flank-type?) ──
print("── CUT 7 — DOF verdict × flank-type (the confound) " + "─" * 8)
print(f"{'pair':>16} | {'DOF-solv n/corr':>18} | {'DOF-unsolv n/corr':>18}")
for p in upairs[:6]:
    row = []
    for m in (C["dof"], ~C["dof"]):
        mm = m & (pair == p)
        if mm.sum() < 3:
            row.append(f"{int(mm.sum()):>5} {'—':>8}")
            continue
        c = np.corrcoef(C["fg"][mm], C["fo"][mm])[0, 1] if C["fo"][mm].std() > 1e-9 else float("nan")
        row.append(f"{int(mm.sum()):>5} {c:>8.3f}")
    print(f"{p:>16} | {row[0]:>18} | {row[1]:>18}")

# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE HONEST METRIC — per-position correlation ACROSS scenarios. Hold a boundary FIXED (same chain node id,
# stable across scenarios) and vary the gDNA level: does its solved f_g track its oracle f_g? This removes the
# between-structure (Simpson) confound the pooled corr suffers. A position is "resolving" if corr_i is high.
# Aggregate corr_i by stratum (mean, mass-weighted). Requires ≥ MINOBS scenarios with finite fo AND fo-variance.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════
_MINOBS = 6
nid = C["nid"]
uniq = np.unique(nid)
pos = {}  # nid -> per-position summary
for u in uniq:
    m = nid == u
    fog, fgg = C["fo"][m], C["fg"][m]
    if m.sum() < _MINOBS or fog.std() < 1e-6 or fgg.std() < 1e-6:
        ci = float("nan")
    else:
        ci = np.corrcoef(fgg, fog)[0, 1]
    # a position's static covariates (constant across scenarios; take the modal / first)
    pos[u] = dict(
        corr=ci, mass=C["mass"][m].mean(),
        spl=(C["spl"][m] > _EPS).mean() > 0.5,
        ltyp=int(np.round(np.median(C["ltyp"][m]))), rtyp=int(np.round(np.median(C["rtyp"][m]))),
        dof=(C["dof"][m].mean() > 0.5),
        gap=np.median(C["gap"][m]),
        bothmsg=((C["pgl"][m] > _EPS) & (C["pgr"][m] > _EPS)).mean() > 0.5,
    )
P = {k: np.array([pos[u][k] for u in uniq]) for k in ("corr", "mass", "spl", "ltyp", "rtyp", "dof", "gap", "bothmsg")}
valid = np.isfinite(P["corr"])
print(f"\n{'='*90}\nHONEST per-position metric — {valid.sum()}/{len(uniq)} boundary positions with ≥{_MINOBS} "
      f"informative scenarios\n{'='*90}")


def pos_report(title, groups):
    print(f"── {title} " + "─" * max(2, 66 - len(title)))
    print(f"{'group':>26} | {'n_pos':>6} | {'mean corr_i':>11} | {'mass-wtd corr_i':>15} | note")
    for label, m in groups:
        mm = m & valid
        if mm.sum() < 5:
            print(f"{label:>26} | {int(mm.sum()):>6} | {'—':>11} |")
            continue
        c = P["corr"][mm]
        mc = float(np.mean(c))
        wc = float(np.average(c, weights=P["mass"][mm]))
        note = "RESOLVES" if mc > 0.35 else ("coin-flip" if mc < 0.20 else "")
        print(f"{label:>26} | {int(mm.sum()):>6} | {mc:>11.3f} | {wc:>15.3f} | {note}")
    print()


pos_report("P1 — spliced present", [
    ("has spliced", P["spl"]),
    ("no spliced", ~P["spl"]),
])
ppair = np.array([_pair(int(lf), int(rf)) for lf, rf in zip(P["ltyp"], P["rtyp"])])
pos_report("P2 — flank-type pair", [(p, ppair == p) for p in sorted(set(ppair))])
pos_report("P3 — spliced × intron|exon (control for flank-type)", [
    ("intron|exon, spliced", (ppair == "intron|exon") & P["spl"]),
    ("intron|exon, no spliced", (ppair == "intron|exon") & ~P["spl"]),
])
pos_report("P4 — DOF verdict (per-position)", [
    ("DOF solvable", P["dof"]),
    ("DOF unsolvable", ~P["dof"]),
])
pos_report("P5 — both flanks message", [
    ("both-flank msg", P["bothmsg"]),
    ("not both", ~P["bothmsg"]),
])
