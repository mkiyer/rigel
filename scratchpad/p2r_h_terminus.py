"""P-2 / PASS-0 — IS THE REFRAME BLOW-UP THE MISSING TSS/TES BIT?

Part G traced the gDNA over-call to the reframe (96 % of a +1.56-decade error) and the individual nodes all
look the same:

    node 3087 (exon, true f_g=0.008)  <- RIGHT boundary 3088, pure gDNA (RNA=0), rho_g measured EXACTLY
                                         right (0.2099), r = 24.07/0.2099 = 114.6, delivered 127x too big.
    its LEFT boundary 3086            <- also pure gDNA, but r = only 6.0.

The difference between those two boundaries is structural: `G0316.1` spans 6,934,891-6,946,101 and its last
exon ENDS at 6,946,101 — the right boundary is a transcript **TERMINUS**, so no RNA continues across it and
its total density is pure gDNA, while the exon it feeds is RNA-dominated.  `r = ρ_tot(dst)/ρ_tot(src)` is
then the RNA ratio, not a capture step.

THE HYPOTHESIS: **`r` blows up exactly at boundaries that carry a transcript START or END**, because that is
where the RNA content of the two sides is structurally different — which is the SAME missing bit that
`ω_graft` (P1d) was measured to split ≥30× on, appearing in a second, independent place.

THE TEST: classify every source boundary as TERMINUS (some transcript starts or ends there) vs JUNCTION-ONLY
(only internal splice junctions), and split `log10 r` and the delivered gDNA error by that bit.

FALSIFIER: if `log10 r` is the same at termini and junction-only boundaries, the bit is not what drives it
and the cause is something else (source RNA content for a non-structural reason, or the frame itself).

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_h_terminus.py [--cond ...]
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

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-12
CONDS = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",
]


def boundary_classes(ix):
    """position -> (is_terminus, is_junction) over every transcript in the annotation."""
    term, junc = set(), set()
    for _, t in ix.t_df.iterrows():
        if bool(t["is_nrna"]):
            continue                      # nRNA spans are synthetic rows over the same locus
        ex = ix.get_exon_intervals(int(t["t_index"]))
        if ex is None or len(ex) == 0:
            continue
        term.add((t["ref"], int(ex[0][0])))
        term.add((t["ref"], int(ex[-1][1])))
        for k in range(len(ex) - 1):
            junc.add((t["ref"], int(ex[k][1])))
            junc.add((t["ref"], int(ex[k + 1][0])))
    return term, junc


def run(cond, ix, ra, cfg, term, junc):
    inp = _scan_and_truth(SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    cap, chain = dbg["capture"], dbg["chain"]
    S = cap["_uni_static"]
    E_g = np.asarray(S["E_g"], float)
    rl, rr = np.asarray(S["rho_l0"], float), np.asarray(S["rho_r0"], float)
    mass = np.asarray(cap["mass_global"], float)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    cls = np.where(kind != REGION, 3, rt)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    rho_true = G / np.maximum(E_g, _EPS)
    ridx = np.asarray(chain.ref_idx, np.int64)
    # per-node structural class for BOUNDARY nodes
    is_term = np.zeros(len(kind), bool)
    is_junc = np.zeros(len(kind), bool)
    bdf = ix.boundary_df
    for i in np.nonzero(kind != REGION)[0]:
        row = bdf.iloc[int(ridx[i])]
        key = (row["ref_name"], int(row["position"]))
        is_term[i] = key in term
        is_junc[i] = key in junc
    out = []
    for rec, dstf, srcf in ((cap["_pin"][-2], rl, rr), (cap["_pin"][-1], rr, rl)):
        src = np.asarray(rec["src"], np.int64)
        valid = np.asarray(rec["valid"], bool)
        tg, tpg = np.asarray(rec["tg"], float), np.asarray(rec["tpg"], float)
        framed = valid & (srcf[src] > _EPS) & (dstf > _EPS)
        r = np.where(framed, dstf / np.maximum(srcf[src], _EPS), 1.0)
        sel = valid & (tpg > 0) & (tg > _EPS) & (cls == 2) & (rho_true > _EPS) & (mass > _EPS)
        sel = sel & (cls[src] == 3)      # boundary sources only (100 % of them per part F)
        i = np.nonzero(sel)[0]
        if not i.size:
            continue
        s = src[i]
        out.append(dict(
            logr=np.log10(np.maximum(r[i], _EPS)),
            err=np.log10(tg[i] / np.maximum(rho_true[i], _EPS)),
            w=mass[i],
            term=is_term[s], junc=is_junc[s],
            src_rna_frac=np.where(G[s] + R[s] > _EPS, R[s] / np.maximum(G[s] + R[s], _EPS), np.nan),
            src_face=srcf[s], dst_face=dstf[i],
        ))
    return {k: np.concatenate([o[k] for o in out]) for k in out[0]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", nargs="*", default=CONDS)
    a = ap.parse_args()
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    term, junc = boundary_classes(ix)
    print(f"# annotation: {len(term)} terminus positions, {len(junc)} junction positions, "
          f"{len(term & junc)} both\n")

    print("=== H1. log10(r) AND THE DELIVERED gDNA ERROR, SPLIT BY THE SOURCE BOUNDARY'S STRUCTURAL BIT ===")
    print("    TERMINUS = a transcript starts or ends there (no RNA continues across).\n")
    print(f"{'condition':40s} {'class':16s} {'n':>6s} {'mass%':>6s} {'med log10 r':>12s} "
          f"{'mean log10 r':>13s} {'med err(dec)':>13s} {'med err(x)':>11s}")
    res = {}
    for cond in a.cond:
        d = run(cond, ix, ra, cfg, term, junc)
        res[cond] = d
        tw = d["w"].sum()
        for nm, m in (("TERMINUS", d["term"]),
                      ("junction-only", ~d["term"] & d["junc"]),
                      ("neither", ~d["term"] & ~d["junc"])):
            if not m.any():
                continue
            print(f"{cond[5:]:40s} {nm:16s} {int(m.sum()):6d} {100 * d['w'][m].sum() / tw:5.1f}% "
                  f"{np.median(d['logr'][m]):12.3f} "
                  f"{np.average(d['logr'][m], weights=d['w'][m]):13.3f} "
                  f"{np.median(d['err'][m]):13.3f} {10 ** np.median(d['err'][m]):10.1f}x")
        print()

    print("\n=== H2. THE MECHANISM, DIRECTLY: r vs the SOURCE FACE's RNA content ===")
    print("    r = rho_face(dst)/rho_face(src).  If the source face is pure gDNA it is tiny, so r explodes.")
    print("    Binned by the source boundary's TRUE RNA fraction (oracle).\n")
    d = res[a.cond[0]]
    ok = np.isfinite(d["src_rna_frac"])
    edges = [0.0, 1e-9, 0.05, 0.25, 0.5, 0.75, 0.95, 1.01]
    lab = ["exactly 0", "0-5%", "5-25%", "25-50%", "50-75%", "75-95%", ">95%"]
    print(f"{'source face RNA fraction':26s} {'n':>6s} {'%TERMINUS':>10s} {'med log10 r':>12s} "
          f"{'med src face rho':>17s} {'med delivered err':>18s}")
    for k in range(len(edges) - 1):
        m = ok & (d["src_rna_frac"] >= edges[k]) & (d["src_rna_frac"] < edges[k + 1])
        if not m.any():
            continue
        print(f"{lab[k]:26s} {int(m.sum()):6d} {100 * d['term'][m].mean():9.1f}% "
              f"{np.median(d['logr'][m]):12.3f} {np.median(d['src_face'][m]):17.4f} "
              f"{np.median(d['err'][m]):17.3f}")

    print("\n\n=== H3. HOW MUCH OF THE ERROR MASS IS ON TERMINUS EDGES? ===")
    print(f"{'condition':40s} {'edges TERM':>11s} {'errmass TERM':>13s} {'errmass junc':>13s}")
    for cond in a.cond:
        d = res[cond]
        em = d["w"] * np.abs(d["err"])
        t = d["term"]
        print(f"{cond[5:]:40s} {100 * t.mean():10.1f}% {100 * em[t].sum() / max(em.sum(), _EPS):12.1f}% "
              f"{100 * em[~t & d['junc']].sum() / max(em.sum(), _EPS):12.1f}%")


if __name__ == "__main__":
    main()
