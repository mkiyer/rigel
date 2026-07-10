"""Per-node calibration error attribution by INFORMATION SOURCE (priority 1=strand, 2=imputation, 3=global).

For one condition: run production calibrate (capturing the converged belief + the sweep's internal geometry via
a node_sweep wrapper), then for every solvable node decompose the per-node solve into:
  P1  strand-only  : _local_loglik with the strand likelihood + (single-strand) Jeffreys ONLY  -> f_g + var
  P1+P3  +global   : add the count-space global NB log-prior
  P1+P2  +imputation: add the gDNA + per-strand RNA messages (recomputed from the CONVERGED neighbours)
  FULL              : the actual converged f_g (P1+P2+P3)
and compare each to the by-origin ORACLE f_g. Classifies each node + aggregates the mass-weighted error
(Σ (f_g−oracle)·mass) by category, so we can see whether the error is in nodes that SHOULD be strand-solvable
(P1 overruled — a priority bug) vs genuinely-ambiguous AMBIG nodes (that need the P2/P3 tiebreaker).

    OMP_NUM_THREADS=1 python scripts/debug/node_error_attribution.py [condition] [--top N]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.simplex_sweep import _fg_median, _fg_var, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}
_CAP = {}


def _wrap_node_sweep(orig):
    def w(chain, statics, geometry, belief, region_arrays, boundary_substrate, **kw):
        out = orig(chain, statics, geometry, belief, region_arrays, boundary_substrate, **kw)
        _CAP.update(chain=chain, statics=statics, geom=geometry, bsub=boundary_substrate,
                    belief_in=belief, belief_out=out[0], kappa=kw["rna_sense_frac"],
                    od_g=kw.get("gdna_strand_overdispersion", 0.0),
                    od_r=kw.get("rna_strand_overdispersion", 0.0), n_grid=kw["n_grid"])
        return out
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("condition", nargs="?", default="gdna_gdna300_ss_0.99_nrna_none_capture_on")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    index, blob = build_or_load_cache(args.condition, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    calibrate.__globals__["node_sweep"] = _wrap_node_sweep(bp.node_sweep)  # patch calibrate's own namespace
    calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

    ch, st, geom, bsub = _CAP["chain"], _CAP["statics"], _CAP["geom"], _CAP["bsub"]
    bi, bo = _CAP["belief_in"], _CAP["belief_out"]
    kappa, od_g, od_r = _CAP["kappa"], _CAP["od_g"], _CAP["od_r"]
    lat = _simplex_lattice(int(_CAP["n_grid"]))
    fpg, fng, fgg = lat

    # the count-space global, exactly as the sweep built it (seeds use the INIT belief)
    rho_g, gvm, var_mean = bp._gdna_seed_estimate(ch, st, geom, ra, bsub, np.asarray(bi.f_g, float).copy(), kappa)
    kind = np.asarray(ch.kind)
    is_reg = kind == REGION
    EG = (geom.eff_gdna_left, geom.eff_gdna_right)
    MS = (geom.mass_left, geom.mass_right)
    eff_global = np.where(is_reg, EG[0], 0.5 * (EG[0] + EG[1]))
    mass_global = np.where(is_reg, MS[0], MS[0] + MS[1])
    glp = bp._global_logprior(fgg, mass_global, eff_global, rho_g, gvm, var_mean)
    s2 = max(float(gvm.predict(np.array([max(rho_g, _EPS)]))[0]), 0.0)
    mu_glob = np.clip(rho_g * eff_global / np.maximum(mass_global, _EPS), 0.0, 1.0)
    N_glob = np.full_like(mass_global, rho_g * rho_g / max(var_mean + s2, _EPS))  # M-independent scalar

    # oracle f_g per chain node: regions from contained gdna/rna; boundaries from crossing gdna/rna.
    g_reg = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    r_reg = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    bsg, bsr = BoundarySubstrate.from_payload(blob["payload_gdna"]), BoundarySubstrate.from_payload(blob["payload_rna"])
    gb = np.asarray(bsg.left.mass_unspliced, float) + np.asarray(bsg.right.mass_unspliced, float)
    rb = np.asarray(bsr.left.mass_unspliced, float) + np.asarray(bsr.right.mass_unspliced, float)
    refidx = np.asarray(ch.ref_idx)
    o_fg = np.full(ch.n_nodes, np.nan)
    o_fg[is_reg] = g_reg[refidx[is_reg]] / (g_reg[refidx[is_reg]] + r_reg[refidx[is_reg]] + _EPS)
    o_fg[~is_reg] = gb[refidx[~is_reg]] / (gb[refidx[~is_reg]] + rb[refidx[~is_reg]] + _EPS)

    # per-chain-node class/type (regions from signature; boundaries inherit a coarse label from flanks)
    sig = np.asarray(ra.signature)
    scl = np.full(ch.n_nodes, -1)
    tcl = np.full(ch.n_nodes, -1)
    scl[is_reg] = [coarse_strand_from_signature(int(sig[r])) for r in refidx[is_reg]]
    tcl[is_reg] = [coarse_type_from_signature(int(sig[r])) for r in refidx[is_reg]]

    fp, fn = st.free_pos, st.free_neg
    f_full = np.asarray(bo.f_g, float)
    f_pos_c, f_neg_c = np.asarray(bo.f_pos, float), np.asarray(bo.f_neg, float)
    solvable = (fp | fn) & (st.mass_unspliced > 0.0)

    def solve(i, *, glob=False, msg=False):
        gl = glp[i:i + 1] if glob else None
        mg = mn = mp = 0.0
        ng = nn = npz = 0.0
        if msg:  # recompute the dominant gDNA + per-strand messages from the converged neighbours
            for nbr, df, sf in ((np.asarray(ch.left), 0, 1), (np.asarray(ch.right), 1, 0)):
                src = nbr[i]
                if src < 0 or MS[sf][src] <= _EPS or not solvable[src]:
                    continue
                sm = MS[sf][src]
                a, b = bp._message(f_full[src] * sm / max(EG[sf][src], _EPS), EG[sf][src], EG[df][i],
                                   MS[df][i], 0.0, gvm)
                if b > ng:
                    mg, ng = a, b
        psi = _local_loglik(st.u_pos[i:i + 1], st.u_neg[i:i + 1], st.spliced_pos[i:i + 1],
                            st.spliced_neg[i:i + 1], fp[i:i + 1], fn[i:i + 1], kappa, od_g, od_r, lat,
                            strand_obs=st.strand_obs[i:i + 1], global_logprior=gl,
                            gdna_imp_frac=np.array([mg]), gdna_imp_count=np.array([ng]),
                            rna_imp_frac=(np.array([mp]), np.array([mn])),
                            rna_imp_count=(np.array([npz]), np.array([nn])))
        return float(np.clip(_fg_median(psi, fgg)[0], 0, 1)), float(_fg_var(psi, fgg)[0]), mg, ng

    nodes = np.where(solvable & np.isfinite(o_fg))[0]
    rows = []
    for i in nodes:
        fs, vs, _, _ = solve(i, glob=False, msg=False)          # P1 strand-only
        fsg, _, _, _ = solve(i, glob=True, msg=False)           # P1+P3
        fsi, _, mg, ng = solve(i, glob=False, msg=True)         # P1+P2
        m = mass_global[i]
        rows.append(dict(i=int(i), reg=is_reg[i], cls=SC.get(int(scl[i]), "bnd"),
                         typ=TC.get(int(tcl[i]), "-"), mass=m, o=o_fg[i], strand=fs, svar=vs,
                         sg=fsg, si=fsi, full=f_full[i], muN=(mg, ng),
                         err=(f_full[i] - o_fg[i]) * m, strand_err=(fs - o_fg[i]) * m))

    # classification
    STR_TOL, VAR_TOL = 0.12, 0.02

    def classify(r):
        strand_ok = abs(r["strand"] - r["o"]) < STR_TOL and r["svar"] < VAR_TOL
        full_ok = abs(r["full"] - r["o"]) < STR_TOL
        if strand_ok and full_ok:
            return "P1_solved_ok"
        if strand_ok and not full_ok:
            return "P1_OVERRULED"          # strand had it right, full got it wrong
        return "needs_P2P3"                 # strand uninformative/wrong (AMBIG-like)

    for r in rows:
        r["cat"] = classify(r)

    tot = sum(abs(r["err"]) for r in rows)
    print(f"=== {args.condition} : per-node calibration error attribution by info source ===")
    print(f"solvable nodes={len(rows)}  Σ|mass-wtd f_g err|={tot:,.0f}  (rho_global={rho_g:.4g} sigma2_g*={s2:.4g})\n")
    print("ERROR BY CATEGORY (mass-weighted |Δf_g|, signed net):")
    for cat in ["P1_solved_ok", "P1_OVERRULED", "needs_P2P3"]:
        sub = [r for r in rows if r["cat"] == cat]
        ae = sum(abs(r["err"]) for r in sub)
        ne = sum(r["err"] for r in sub)
        print(f"  {cat:>14}: n={len(sub):>4}  |err|={ae:>11,.0f} ({100*ae/max(tot,1):>4.1f}%)  net={ne:>+11,.0f}")
    print("\n  split by node kind × class:")
    for reg in (True, False):
        for cls in ("AMBIG", "NEG", "POS", "NONE", "bnd"):
            sub = [r for r in rows if r["reg"] == reg and r["cls"] == cls]
            if not sub:
                continue
            ae = sum(abs(r["err"]) for r in sub)
            if ae < 500:
                continue
            kind_s = "region" if reg else "bound "
            print(f"    {kind_s} {cls:>6}: n={len(sub):>4} |err|={ae:>11,.0f}  net={sum(r['err'] for r in sub):>+11,.0f}")

    rows.sort(key=lambda r: -abs(r["err"]))
    print(f"\nTOP {args.top} error nodes  (strand=P1-only f_g, svar=strand confidence(var), sg=+global, si=+imput):")
    print(f"{'node':>5} {'kind':>4} {'cls':>5} {'typ':>5} {'mass':>8} {'oracle':>6} {'strand':>6} {'svar':>6} "
          f"{'sg':>6} {'si':>6} {'FULL':>6} {'msgμ/N':>11} {'err':>9} {'cat':>13}")
    for r in rows[:args.top]:
        mg, ng = r["muN"]
        print(f"{r['i']:>5} {'reg' if r['reg'] else 'bnd':>4} {r['cls']:>5} {r['typ']:>5} {r['mass']:>8,.0f} "
              f"{r['o']:>6.3f} {r['strand']:>6.3f} {r['svar']:>6.3f} {r['sg']:>6.3f} {r['si']:>6.3f} "
              f"{r['full']:>6.3f} {mg:>5.2f}/{ng:>5.0f} {r['err']:>+9,.0f} {r['cat']:>13}")


if __name__ == "__main__":
    main()
