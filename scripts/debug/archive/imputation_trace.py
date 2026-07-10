"""Trace WHY imputation fails at the worst flagship calibration-error nodes (strand + imputation focus).

The global is now a gentle fallback, so the remaining AMBIG/boundary error must be resolvable (or not) by the
strand (priority-1) and imputation (priority-2). For each worst node this prints: its strand answer + confidence,
its oracle, and — for each chain neighbour — whether the neighbour is a CONFIDENT gDNA SOURCE (strand-solvable,
correct) or itself ambiguous/under-called, plus the gDNA message it sends (μ, N_src) and the message components
(ρ_src, τ_ρ). Classifies the cause: no-confident-source (cluster of balanced nodes ⇒ genuinely needs the
fallback) vs weak-precision (a good source exists but N_src too small) vs source-undercalled (cascade).

    OMP_NUM_THREADS=1 python scripts/debug/imputation_trace.py [condition] [--top N]
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
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _fg_var, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMB"}
TC = {0: "intg", 1: "intr", 2: "exon"}
_CAP = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("condition", nargs="?", default="gdna_gdna300_ss_0.99_nrna_none_capture_on")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    index, blob = build_or_load_cache(args.condition, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    def wrap(orig):
        def w(chain, st, geom, belief, rga, bsub, **k):
            out = orig(chain, st, geom, belief, rga, bsub, **k)
            _CAP.update(chain=chain, st=st, geom=geom, bi=belief, bo=out[0], bsub=bsub,
                        kappa=k["rna_sense_frac"], odg=k.get("gdna_strand_overdispersion", 0.0),
                        odr=k.get("rna_strand_overdispersion", 0.0), ng=k["n_grid"])
            return out
        return w

    calibrate.__globals__["node_sweep"] = wrap(bp.node_sweep)
    calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    ch, st, geom, bi, bo, bsub = (_CAP[k] for k in ("chain", "st", "geom", "bi", "bo", "bsub"))
    kappa, odg, odr = _CAP["kappa"], _CAP["odg"], _CAP["odr"]
    lat = _simplex_lattice(_CAP["ng"])
    fpg, fng, fgg = lat
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    is_reg = np.asarray(ch.kind) == REGION
    refidx = np.asarray(ch.ref_idx)
    EG = (geom.eff_gdna_left, geom.eff_gdna_right)
    MS = (geom.mass_left, geom.mass_right)
    fp, fn = st.free_pos, st.free_neg
    solvable = (fp | fn) & (st.mass_unspliced > 0.0)
    rho_g, gvm, _ = bp._gdna_seed_estimate(ch, st, geom, ra, bsub, np.asarray(bi.f_g, float).copy(), kappa)

    # oracle f_g per node
    gR = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    rR = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    bsg, bsr = BoundarySubstrate.from_payload(blob["payload_gdna"]), BoundarySubstrate.from_payload(blob["payload_rna"])
    gB = np.asarray(bsg.left.mass_unspliced, float) + np.asarray(bsg.right.mass_unspliced, float)
    rB = np.asarray(bsr.left.mass_unspliced, float) + np.asarray(bsr.right.mass_unspliced, float)
    o_fg = np.where(is_reg, gR[np.clip(refidx, 0, len(gR) - 1)] / (gR[np.clip(refidx, 0, len(gR) - 1)] + rR[np.clip(refidx, 0, len(rR) - 1)] + _EPS),
                    gB[np.clip(refidx, 0, len(gB) - 1)] / (gB[np.clip(refidx, 0, len(gB) - 1)] + rB[np.clip(refidx, 0, len(rB) - 1)] + _EPS))

    sig = np.asarray(ra.signature)
    cls = np.full(ch.n_nodes, -1)
    typ = np.full(ch.n_nodes, -1)
    cls[is_reg] = [coarse_strand_from_signature(int(sig[r])) for r in refidx[is_reg]]
    typ[is_reg] = [coarse_type_from_signature(int(sig[r])) for r in refidx[is_reg]]

    f_full = np.asarray(bo.f_g, float)
    mass_node = np.where(is_reg, MS[0], MS[0] + MS[1])
    err = (f_full - o_fg) * mass_node

    def strand_solve(i):  # P1-only: strand + Jeffreys, no global/imputation
        psi = _local_loglik(st.u_pos[i:i + 1], st.u_neg[i:i + 1], st.spliced_pos[i:i + 1],
                            st.spliced_neg[i:i + 1], fp[i:i + 1], fn[i:i + 1], kappa, odg, odr, lat,
                            strand_obs=st.strand_obs[i:i + 1])
        return float(_fg_median(psi, fgg)[0]), float(_fg_var(psi, fgg)[0])

    def gmsg(src, df, sf, i):  # gDNA message src->i (uses converged f_g[src])
        if src < 0 or MS[sf][src] <= _EPS or not solvable[src]:
            return None
        rho = f_full[src] * MS[sf][src] / max(EG[sf][src], _EPS)
        mu, N = bp._message(rho, EG[sf][src], EG[df][i], MS[df][i], 0.0, gvm)
        return mu, N, rho

    order = np.argsort(-np.abs(err))
    order = [i for i in order if solvable[i] and np.isfinite(o_fg[i])][: args.top]
    print(f"=== {args.condition}: imputation trace for top-{args.top} error nodes (rho_global={rho_g:.3g}) ===")
    print("strand_obs: is the node's OWN strand informative? source = a neighbour that's strand-solvable+correct\n")
    for i in order:
        sfg, svar = strand_solve(i)
        kind_s = "REG" if is_reg[i] else "BND"
        ci = SC.get(int(cls[i]), "-") if is_reg[i] else "bnd"
        print(f"node {i} [{kind_s} {ci} {TC.get(int(typ[i]),'-')}] mass={mass_node[i]:,.0f}  "
              f"oracle={o_fg[i]:.3f} FULL={f_full[i]:.3f} err={err[i]:+,.0f} | "
              f"strand f_g={sfg:.3f}(var{svar:.3f}, obs={bool(st.strand_obs[i])})")
        for nbr, df, sf, lab in ((left, 0, 1, "L"), (right, 1, 0, "R")):
            s = nbr[i]
            if s < 0:
                print(f"    {lab}: terminal")
                continue
            nk = "REG" if is_reg[s] else "BND"
            nc = SC.get(int(cls[s]), "-") if is_reg[s] else "bnd"
            ns, nv = strand_solve(s)
            src_ok = bool(st.strand_obs[s]) and abs(ns - o_fg[s]) < 0.15  # confident strand source?
            m = gmsg(s, df, sf, i)
            msg = f"gDNAmsg μ={m[0]:.3f} N={m[1]:.2f}" if m else "no msg (empty/unsolvable)"
            print(f"    {lab}: node{s}[{nk} {nc}] f_g={f_full[s]:.3f} oracle={o_fg[s]:.3f} "
                  f"strand={ns:.3f}(obs={bool(st.strand_obs[s])}) confident_source={src_ok} | {msg}")
        print()


if __name__ == "__main__":
    main()
