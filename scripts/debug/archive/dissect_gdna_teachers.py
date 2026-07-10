"""Dissect the Phase-2 teacher densities — isolate the geometry bias (D1) + the ρ_resid weighting/level
(D2 follow-ups). Runs the real calibrate() (toy_prod _debug hook) and dumps, per chain node, the density
decomposition vs the oracle:

  * REGION: E_gdna, M_contained, oracle gDNA count, oracle density (gcount/E), solved f_g, solved density
    (f_g·M/E), and solved/oracle ratio. Under UNIFORM gDNA (capture off) every oracle density must equal
    ρ_uniform — any spread is an EFF-LENGTH model bias; any solved≠oracle is a SOLVE bias.
  * BOUNDARY: the crossing density each side presents, vs its flanking regions' contained densities (should
    match under uniform gDNA — a mismatch is the contained-vs-crossing normalization bias).

Usage: python -m scripts.debug.dissect_gdna_teachers [--kappa K] [--gdna G] [--no-capture]
"""
from __future__ import annotations

import argparse
import collections

import numpy as np
import pysam

from scripts.debug.toy_prod import run, SCRATCH  # noqa: F401  (SCRATCH re-exported for convenience)
from rigel.calibration.bp_solver import REGION, BOUNDARY, node_densities, _node_region_type
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG

_RT = {0: "intergenic", 1: "intron", 2: "exon"}


def _default_genes(n_genes: int):
    genes, x, INTRON, GAP, SIZES = [], 5000, 3000, 6000, (300, 1200, 5000)
    for k in range(n_genes):
        s = x
        ex = []
        for sz in SIZES:
            ex.append((s, s + sz)); s = s + sz + INTRON
        genes.append((f"G{k}", "+" if k % 2 == 0 else "-", ex, 100.0 if k % 2 == 0 else 0.0))
        x = s + GAP
    return genes, x + GAP


def _oracle_region_gcount(bam_path, region_arrays):
    """Per-region contained oracle gDNA fragment count (read1 mate+TLEN span fully inside the region).
    Single-ref toy (chr1) → position lookup only, no ref-name mapping."""
    starts = np.asarray(region_arrays.start); ends = np.asarray(region_arrays.end); R = starts.shape[0]
    order = np.argsort(starts)
    s_sorted, e_sorted, idx_sorted = starts[order], ends[order], order
    g = np.zeros(R)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if (rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1
                or rd.mate_is_unmapped or "gdna" not in rd.query_name.lower()):
            continue
        tl = abs(rd.template_length)
        if tl == 0:
            continue
        lo = min(rd.reference_start, rd.next_reference_start); hi = lo + tl
        j = np.searchsorted(s_sorted, lo, side="right") - 1
        if 0 <= j < R and s_sorted[j] <= lo and hi <= e_sorted[j]:
            g[idx_sorted[j]] += 1.0
    bam.close()
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.99)
    ap.add_argument("--n-genes", type=int, default=6)
    ap.add_argument("--gdna", type=float, default=0.3)
    ap.add_argument("--no-capture", action="store_true")
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    genes, glen = _default_genes(args.n_genes)
    dbg: dict = {}
    run("dissect", genes, kappa=args.kappa, n_rna=12000 * args.n_genes, gdna_fraction=args.gdna,
        capture=not args.no_capture, capture_strength=20.0, genome_length=glen, seed=args.seed, _debug=dbg)

    chain = dbg["chain"]; belief = dbg["belief"]; geom = dbg["geometry"]; ra = dbg["region_arrays"]
    reff = np.asarray(dbg["region_eff_len"], dtype=np.float64)
    dens = node_densities(belief, geom)
    node_rtype, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind); ridx = np.asarray(chain.ref_idx)
    fg = np.asarray(belief.f_g)

    gcount = _oracle_region_gcount(dbg["bam_path"], ra)

    print(f"\n=== κ_fit={float(dbg['rna_sense_frac']):.3f}  capture={'off' if args.no_capture else 'on'}  "
          f"gdna={args.gdna}  n_genes={args.n_genes} ===")
    # ---- REGION table ----
    print("\nREGION nodes (oracle density should be UNIFORM under capture-off):")
    print(f"  {'type':>10} {'L':>6} {'E_gdna':>7} {'M':>9} {'gcnt':>7} {'orc_dns':>8} {'solv_fg':>7} "
          f"{'solv_dns':>8} {'s/o':>6}")
    R = node_rtype.shape[0] if False else int(np.max(ridx[kind == REGION])) + 1
    rows = collections.defaultdict(list)
    for nd in np.where(kind == REGION)[0]:
        ri = int(ridx[nd])
        L = float(ra.region_size_bp[ri]); E = float(reff[ri]); M = float(geom.mass_left[nd])
        solv_dns = float(dens.rho_g_left[nd])
        gc = float(gcount[ri]) if gcount is not None else float("nan")
        orc_dns = gc / E if E > 0 else float("nan")
        so = solv_dns / orc_dns if orc_dns > 0 else float("nan")
        rt = _RT.get(int(node_rtype[nd]), "?")
        rows[(rt, int(round(L)))].append((solv_dns, orc_dns))
        print(f"  {rt:>10} {int(L):>6} {E:>7.0f} {M:>9.1f} {gc:>7.0f} {orc_dns:>8.4f} {fg[nd]:>7.3f} "
              f"{solv_dns:>8.4f} {so:>6.2f}")

    print("\n  by (type,size): mean oracle_dns  vs  mean solved_dns  (log; spread should be ~0 off-capture)")
    for key in sorted(rows):
        o = np.array([r[1] for r in rows[key]]); s = np.array([r[0] for r in rows[key]])
        o = o[np.isfinite(o) & (o > 0)]; s = s[np.isfinite(s) & (s > 0)]
        lo = np.log(o).mean() if o.size else float("nan")
        ls = np.log(s).mean() if s.size else float("nan")
        print(f"    {key[0]:>10} L={key[1]:>5}: oracle logρ={lo:+.2f}  solved logρ={ls:+.2f}  Δ={ls - lo:+.2f}")

    # ---- BOUNDARY vs flanking regions ----
    print("\nBOUNDARY crossing density vs flanking region contained density (should match off-capture):")
    left = np.asarray(chain.left); right = np.asarray(chain.right)
    print(f"  {'bnd':>5} {'Lreg_dns':>9} {'bnd_L':>8} {'bnd_R':>8} {'Rreg_dns':>9}")
    shown = 0
    for nd in np.where(kind == BOUNDARY)[0]:
        lb = int(left[nd]); rb = int(right[nd])
        ldns = float(dens.rho_g_right[lb]) if lb >= 0 and kind[lb] == REGION else float("nan")
        rdns = float(dens.rho_g_left[rb]) if rb >= 0 and kind[rb] == REGION else float("nan")
        bL = float(dens.rho_g_left[nd]); bR = float(dens.rho_g_right[nd])
        if np.isfinite(ldns) or np.isfinite(rdns):
            print(f"  {nd:>5} {ldns:>9.4f} {bL:>8.4f} {bR:>8.4f} {rdns:>9.4f}")
            shown += 1
        if shown >= 12:
            break


if __name__ == "__main__":
    main()
