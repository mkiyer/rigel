"""P-2 RESIDUAL, PART A — is the residual an OFF-SIMPLEX RNA claim?

THE HYPOTHESIS (pre-registered before running).  `_pin_v` used to rescale every incoming claim so that
`Σ_c ρ_c·E_c = M`, which made the delivered log-fraction target `mo_c = log(ρ_c·E_c/M)` **≤ 0 by
construction** — a partial RNA-only claim came out at `1/(1+f_g_own) < 1`.  P-2 stopped that rewrite (it was
building the budget from the DESTINATION's own belief, a BP violation on the mode), and with it went the only
thing holding the target inside the simplex.  So for the first time `mo_R > 0` is reachable: the message
asserts MORE RNA fragments than the node observed in total.

ψ consumes it as `−½·p·(log f_R − mo_R)²` with `log f_R ≤ 0`, so an off-simplex target is a monotone pull to
the RNA vertex whose FORCE AT THE VERTEX is `p·mo_R` — unbounded in the size of the over-claim.  Yet beyond
`f_R = 1` the target carries NO composition information: every `mo_R ≥ 0` says exactly "all of it is RNA".

PREDICTIONS
  1. on unstranded × capture-OFF × gDNA-bearing (the regressing stratum) a large share of the ERROR MASS
     sits on nodes with `mo_R > 0`;
  2. those nodes UNDER-call f_g (solved < oracle) — the direction the residual was reported in;
  3. on `gdna_none` the over-claim is just as common but HARMLESS, because there the truth IS the RNA vertex.
     That asymmetry is what a simplex clip would cost, so it must be sized here, not after the A/B;
  4. stranded capture-OFF barely participates (P-2 left it untouched: +0.0000/+0.0001).

FALSIFICATION: if `mo_R > 0` carries little of the error mass on the regressing conditions, the residual is an
IN-simplex over-claim, no bound at the simplex can help, and this whole direction is wrong.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_a_overclaim.py [--conds substr]
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
CLASSES = ("intergenic", "intron", "exon", "boundary")

# the strata the residual is defined on, plus the controls
GROUPS = [
    ("REGRESSING  unstr x capOFF x gDNA", lambda c: "ss_0.50" in c and "capture_off" in c
     and not c.startswith("gdna_none_")),
    ("control     unstr x capON  x gDNA", lambda c: "ss_0.50" in c and "capture_off" not in c
     and not c.startswith("gdna_none_")),
    ("control     stranded x capOFF", lambda c: "ss_0.99" in c and "capture_off" in c),
    ("control     gdna_none (P-2's win)", lambda c: c.startswith("gdna_none_")),
]


def run(cond, index, ra, cfg):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    cap, chain = dbg["capture"], dbg["chain"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"], float)
    mass = np.asarray(cap["mass_global"], float)
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable_mask"], bool)
    # the delivered RNA target, as a log FRACTION of this node's own observed mass
    mp, pp = np.asarray(cap["mode_p"], float), np.asarray(cap["prec_p"], float)
    mn, pn = np.asarray(cap["mode_n"], float), np.asarray(cap["prec_n"], float)
    mg, pg = np.asarray(cap["mode_g"], float), np.asarray(cap["prec_g"], float)
    # the RNA claim ψ actually feels is the per-strand one (both are evaluated against log(1−f_g));
    # take the live strand's mode, precision-weighted where both are live.
    pr = pp + pn
    moR = np.where(pr > _EPS, (pp * mp + pn * mn) / np.maximum(pr, _EPS), -np.inf)
    return dict(cond=cond, fo=fo, fg=fg, mass=mass, cls=cls, ok=ok, moR=moR, pr=pr,
                mg=mg, pg=pg, err=np.abs(fg - fo))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conds", default=None)
    a = ap.parse_args()
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    if a.conds:
        conds = [c for c in conds if a.conds in c]
    res = {c: run(c, index, ra, cfg) for c in conds}

    print("=== A1. WHERE THE ERROR MASS SITS: RNA claims INSIDE vs OUTSIDE the simplex ===")
    print("   over-claim := mo_R > 0, i.e. the message asserts more RNA fragments than the node observed.\n")
    print(f"{'group':38s} {'n_msg':>7s} {'%nodes':>7s} {'%errmass':>9s} {'med mo_R':>9s} "
          f"{'signed err':>11s} {'signed err':>11s}")
    print(f"{'':38s} {'':>7s} {'over':>7s} {'over':>9s} {'over':>9s} {'OVER':>11s} {'in-simplex':>11s}")
    for gname, sel in GROUPS:
        cs = [c for c in conds if sel(c[5:] if c.startswith("gdna_") else c)]
        if not cs:
            continue
        nm = ov = 0
        em_o = em_i = 0.0
        mo_all, se_o, se_i, w_o, w_i = [], [], [], [], []
        for c in cs:
            r = res[c]
            m = r["ok"] & (r["pr"] > _EPS)
            o = m & (r["moR"] > 0.0)
            i_ = m & (r["moR"] <= 0.0)
            nm += int(m.sum())
            ov += int(o.sum())
            em_o += float((r["mass"] * r["err"])[o].sum())
            em_i += float((r["mass"] * r["err"])[i_].sum())
            mo_all.append(r["moR"][o])
            se_o.append((r["fg"] - r["fo"])[o]); w_o.append(r["mass"][o])
            se_i.append((r["fg"] - r["fo"])[i_]); w_i.append(r["mass"][i_])
        mo_all = np.concatenate(mo_all) if mo_all else np.array([])
        se_o, w_o = np.concatenate(se_o), np.concatenate(w_o)
        se_i, w_i = np.concatenate(se_i), np.concatenate(w_i)
        f_o = np.average(se_o, weights=w_o) if w_o.sum() > 0 else np.nan
        f_i = np.average(se_i, weights=w_i) if w_i.sum() > 0 else np.nan
        print(f"{gname:38s} {nm:7d} {100 * ov / max(nm, 1):6.1f}% "
              f"{100 * em_o / max(em_o + em_i, _EPS):8.1f}% "
              f"{np.median(mo_all) if mo_all.size else np.nan:9.3f} {f_o:+11.3f} {f_i:+11.3f}")

    print("\n   signed err = mass-weighted mean (solved f_g - oracle f_g). NEGATIVE = gDNA UNDER-called.")

    print("\n\n=== A2. PER CONDITION — the over-claim's size and its error share ===")
    print(f"{'condition':50s} {'n_msg':>6s} {'%over':>6s} {'%errmass':>9s} {'p50 mo_R':>9s} "
          f"{'p95 mo_R':>9s} {'signed(over)':>13s}")
    for c in conds:
        r = res[c]
        m = r["ok"] & (r["pr"] > _EPS)
        o = m & (r["moR"] > 0.0)
        i_ = m & (r["moR"] <= 0.0)
        eo = float((r["mass"] * r["err"])[o].sum())
        ei = float((r["mass"] * r["err"])[i_].sum())
        mo = r["moR"][o]
        se = np.average((r["fg"] - r["fo"])[o], weights=r["mass"][o]) if o.any() else np.nan
        print(f"{c[5:]:50s} {int(m.sum()):6d} {100 * o.sum() / max(int(m.sum()), 1):5.1f}% "
              f"{100 * eo / max(eo + ei, _EPS):8.1f}% "
              f"{np.median(mo) if mo.size else np.nan:9.3f} "
              f"{np.quantile(mo, 0.95) if mo.size else np.nan:9.3f} {se:+13.3f}")

    print("\n\n=== A3. THE OVER-CLAIMING NODES, by class (regressing stratum only) ===")
    cs = [c for c in conds if "ss_0.50" in c and "capture_off" in c and not c.startswith("gdna_none_")]
    print(f"{'class':12s} {'n over':>7s} {'errmass':>12s} {'mean|err|':>10s} {'med mo_R':>9s} "
          f"{'med oracle f_g':>15s} {'med solved f_g':>15s}")
    for ci, cn in enumerate(CLASSES):
        nn = 0
        em = 0.0
        mo, orc, slv, ers, ws = [], [], [], [], []
        for c in cs:
            r = res[c]
            m = r["ok"] & (r["pr"] > _EPS) & (r["moR"] > 0.0) & (r["cls"] == ci)
            nn += int(m.sum())
            em += float((r["mass"] * r["err"])[m].sum())
            mo.append(r["moR"][m]); orc.append(r["fo"][m]); slv.append(r["fg"][m])
            ers.append(r["err"][m]); ws.append(r["mass"][m])
        if not nn:
            continue
        mo, orc, slv = np.concatenate(mo), np.concatenate(orc), np.concatenate(slv)
        ers, ws = np.concatenate(ers), np.concatenate(ws)
        print(f"{cn:12s} {nn:7d} {em:12,.0f} {np.average(ers, weights=ws):10.3f} "
              f"{np.median(mo):9.3f} {np.median(orc):15.3f} {np.median(slv):15.3f}")


if __name__ == "__main__":
    main()
