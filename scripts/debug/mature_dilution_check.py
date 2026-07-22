"""EMPIRICAL validation of the exon→IE-boundary MATURE-DILUTION identity, on the grounded full-transcript toy.

THE CLAIM (derivation, 2026-07-22). Capture is nucleic-acid-agnostic (everything at x scales by e(x)):
    intron / boundary-crossing :  (ρ_bg + ν)·e          [mature-free — a spliced fragment cannot cross an
    boundary spliced           :   μ·e                   intron-exon junction contiguously]
    exon unspliced             :  (ρ_bg + ν + μ)·e
Because the boundary's TOTAL (unspliced + spliced) has the SAME component set as the exon, the enrichment
cancels and:
        f_g^boundary / f_g^exon  =  (ρ_bg+ν+μ)/(ρ_bg+ν)  =  (D_B + S_B) / D_B
so the exon→boundary message is the ordinary composition SHIFT plus one additive log term
``log((D_B+S_B)/D_B)`` — enrichment-invariant, zero free constants, and the mature correction uses ONLY the
boundary's own measured spliced/unspliced split (no fragile eff-length extrapolation of within-exon mature).

This script tests BOTH halves against oracle truth, per splice-junction boundary:
  LHS (oracle)   = f_g^B / f_g^E              from the oracle BAM (contiguous crossing vs contained)
  RHS (measured) = (D_B + S_B) / D_B          from the payload geometry
plus the model's premise `mature does not cross contiguously` (oracle mature-crossing fraction ≈ 0).

Usage:  python mature_dilution_check.py [--grid] [--fl-sweep]
"""
from __future__ import annotations
import argparse
import collections
import contextlib
import io
import numpy as np
import pysam

import toy_inject as ti
from rigel.calibration.bp_solver import build_node_geometry
from rigel.calibration.node_chain import REGION
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate

_EPS = 1e-12
EXONS = [(6000, 7000), (11000, 12000), (16000, 17000)]
GENOME = 30000
# the four splice-junction boundary coordinates (exon|intron and intron|exon)
SJ_COORDS = [7000, 11000, 12000, 16000]


def _origin(qname: str) -> str:
    q = qname.lower()
    if q.startswith("gdna"):
        return "gdna"
    if q.startswith("nrna_"):
        return "nascent"
    return "mature"


def oracle_crossing(bam_path, coords):
    """Per coord: contiguously-crossing mass by origin (a block [bs,be) with bs < c < be)."""
    out = {c: collections.Counter() for c in coords}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if rd.is_secondary or rd.is_supplementary or rd.is_unmapped:
            continue
        k = _origin(rd.query_name)
        for bs, be in rd.get_blocks():
            for c in coords:
                if bs < c < be:
                    out[c][k] += 1
    bam.close()
    return out


def oracle_contained(bam_path, ra, ref_names):
    """Per region: contained mass by origin (mate+TLEN span fully inside)."""
    starts, ends, rid = ra.start, ra.end, ra.ref_id
    R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[j] for j in arr]), np.array([ends[j] for j in arr]), np.array(arr))
    name2ref = {n: i for i, n in enumerate(ref_names)}
    out = [collections.Counter() for _ in range(R)]
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if (rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1
                or rd.mate_is_unmapped):
            continue
        rf = name2ref.get(bam.get_reference_name(rd.reference_id))
        if rf is None or rf not in by_ref:
            continue
        tl = abs(rd.template_length)
        if tl == 0:
            continue
        lo = min(rd.reference_start, rd.next_reference_start)
        hi = lo + tl
        s, e, ii = by_ref[rf]
        j = int(np.searchsorted(s, lo, side="right") - 1)
        if 0 <= j < R and s[j] <= lo and hi <= e[j]:
            out[int(ii[j])][_origin(rd.query_name)] += 1
    bam.close()
    return out


def measure(toy, bam_path):
    """Per SJ boundary: oracle LHS, measured RHS, and the mature-crossing premise check."""
    ra, chain, rdf = toy["ra"], toy["chain"], toy["rdf"]
    sub = CalibrationSubstrate.from_payload(toy["payload"], ra)
    bsub = BoundarySubstrate.from_payload(toy["payload"])
    geo = build_node_geometry(chain, sub, bsub, ra, toy["gdna_fl_pmf"], toy["rna_fl_pmf"])
    kind = np.asarray(chain.kind); ref_idx = np.asarray(chain.ref_idx)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    ctype = coarse_type_array(np.asarray(ra.signature)).astype(int)
    starts = rdf["start"].to_numpy(); ends = rdf["end"].to_numpy()

    ml = np.asarray(geo.mass_left, float); mr = np.asarray(geo.mass_right, float)
    egl = np.asarray(geo.eff_gdna_left, float); egr = np.asarray(geo.eff_gdna_right, float)
    spl = (np.asarray(geo.spliced_pos_left, float) + np.asarray(geo.spliced_neg_left, float))
    spr = (np.asarray(geo.spliced_pos_right, float) + np.asarray(geo.spliced_neg_right, float))
    esl = np.asarray(geo.eff_spl_left, float); esr = np.asarray(geo.eff_spl_right, float)

    cross = oracle_crossing(bam_path, SJ_COORDS)
    cont = oracle_contained(bam_path, ra, toy["rdf"].attrs.get("ref_names", ["chr1"]))

    rows = []
    for n in range(chain.n_nodes):
        if kind[n] == REGION:
            continue
        lf, rt = int(left[n]), int(right[n])
        if lf < 0 or rt < 0:
            continue
        # boundary genomic coord = shared edge of its two flank regions
        li, rix = int(ref_idx[lf]), int(ref_idx[rt])
        coord = int(ends[li])
        if coord not in cross or coord != int(starts[rix]):
            continue
        # the adjacent EXON flank (coarse type 2)
        exon_reg = li if ctype[li] == 2 else (rix if ctype[rix] == 2 else -1)
        intron_reg = li if ctype[li] == 1 else (rix if ctype[rix] == 1 else -1)
        if exon_reg < 0 or intron_reg < 0:
            continue  # only exon|intron junctions
        cc = cross[coord]
        g_c, n_c, m_c = cc["gdna"], cc["nascent"], cc["mature"]
        ce = cont[exon_reg]
        g_e, n_e, m_e = ce["gdna"], ce["nascent"], ce["mature"]
        if (g_c + n_c) == 0 or (g_e + n_e + m_e) == 0 or g_e == 0:
            continue
        fgB = g_c / (g_c + n_c)                     # boundary crossing (mature-free) gDNA fraction
        fgE = g_e / (g_e + n_e + m_e)               # exon contained gDNA fraction
        D_B = (ml[n] + mr[n]) / max(egl[n] + egr[n], _EPS)          # unspliced crossing density
        S_B = (spl[n] + spr[n]) / max(esl[n] + esr[n], _EPS)        # spliced (mature) density
        S_B_alt = spl[n] / max(esl[n], _EPS) + spr[n] / max(esr[n], _EPS)  # per-side density sum
        rows.append(dict(coord=coord, fgB=fgB, fgE=fgE, lhs=fgB / max(fgE, _EPS),
                         D_B=D_B, S_B=S_B, S_B_alt=S_B_alt,
                         rhs=(D_B + S_B) / max(D_B, _EPS),
                         rhs_alt=(D_B + S_B_alt) / max(D_B, _EPS),
                         ml=ml[n], mr=mr[n], egl=egl[n], egr=egr[n],
                         spl=spl[n], spr=spr[n], esl=esl[n], esr=esr[n],
                         g_c=g_c, n_c=n_c, g_e=g_e, n_e=n_e, m_e=m_e,
                         mat_cross_frac=m_c / max(g_c + n_c + m_c, _EPS)))
    return rows


def _probe_writer(layout):
    """Probe layouts. 'full' = whole exon; 'centre' = central 50% (DEPLETES junction reads); 'junction' =
    probes centred ON each splice junction (owner: a probe over the junction can enrich the BOUNDARY above
    the exon — the enrichment ratio inverts, which the rule must survive since e cancels)."""
    def w(path, ref, exons, fraction=None):
        with open(path, "w") as f:
            if layout == "junction":
                for i, c in enumerate(SJ_COORDS):
                    s, e = c - 150, c + 150
                    f.write(f"{ref}\t{s}\t{e}\tj{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")
            else:
                frac = 0.5 if layout == "centre" else 1.0
                for i, (s0, e0) in enumerate(exons):
                    wid = int((e0 - s0) * frac)
                    ps = s0 + ((e0 - s0) - wid) // 2
                    pe = ps + wid
                    f.write(f"{ref}\t{ps}\t{pe}\tp{i}\t0\t+\t{ps}\t{pe}\t0\t1\t{pe - ps}\t0\n")
    return w


def run_cell(label, *, capture_on, gdna, nascent, mature, rna_fl=None, gdna_fl=None, n_rna=4000,
             probes="full"):
    _rna_save, _gdna_save, _pw_save = ti.RNA_FL, ti.GDNA_FL, ti._write_probes
    if rna_fl:
        ti.RNA_FL = dict(ti.RNA_FL, **rna_fl)
    if gdna_fl:
        ti.GDNA_FL = dict(ti.GDNA_FL, **gdna_fl)
    ti._write_probes = _probe_writer(probes)
    dbg = {}
    with contextlib.redirect_stdout(io.StringIO()):
        toy = ti.build_toy(f"md_{label}", exons=EXONS, gdna_fraction=gdna, nascent=nascent,
                           mature=mature, capture_on=capture_on, n_rna=n_rna, genome_length=GENOME)
        # rebuild the oracle bam path (build_toy stores it only transiently) — re-derive via scenario dir
    ti.RNA_FL, ti.GDNA_FL, ti._write_probes = _rna_save, _gdna_save, _pw_save
    cand = sorted((ti.SCRATCH / f"inj_md_{label}").glob("*.bam"))
    if not cand:
        return []
    return measure(toy, cand[0])


def report(rows, label):
    if not rows:
        print(f"  {label:>34} | (no usable junctions)")
        return
    lhs = np.array([r["lhs"] for r in rows]); rhs = np.array([r["rhs"] for r in rows])
    mc = np.array([r["mat_cross_frac"] for r in rows])
    rel = (rhs - lhs) / np.maximum(lhs, _EPS)
    print(f"  {label:>34} | n={len(rows):>2} | LHS(oracle)={np.median(lhs):>6.3f} "
          f"RHS(meas)={np.median(rhs):>6.3f} | rel.err={np.median(rel):>+7.1%} "
          f"| mature-cross={np.median(mc):>6.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", action="store_true")
    ap.add_argument("--fl-sweep", action="store_true")
    args = ap.parse_args()

    print("MATURE-DILUTION identity check:  f_g^B/f_g^E  ?=  (D_B+S_B)/D_B")
    print("  (LHS from oracle BAM; RHS from measured payload; mature-cross should be ~0)\n")

    if not args.grid and not args.fl_sweep:
        rows = run_cell("base", capture_on=True, gdna=3.0, nascent=200.0, mature=100.0)
        report(rows, "capture=on g=3.0 n=200 m=100")
        print("\n  RAW per-junction (find the frame error):")
        for r in rows:
            print(f"   coord={r['coord']} | ORACLE cross g={r['g_c']:>6} n={r['n_c']:>6} | "
                  f"exon g={r['g_e']:>6} n={r['n_e']:>6} m={r['m_e']:>6}")
            print(f"     f_gB={r['fgB']:.4f} f_gE={r['fgE']:.4f} LHS={r['lhs']:.4f}   "
                  f"oracle-implied  mu/(rho+nu) = {r['lhs'] - 1:.4f}")
            print(f"     mass  ml={r['ml']:.1f} mr={r['mr']:.1f} | eff_g  egl={r['egl']:.1f} egr={r['egr']:.1f}"
                  f"  -> D_B={r['D_B']:.5f}")
            print(f"     spl   l={r['spl']:.1f} r={r['spr']:.1f} | eff_spl esl={r['esl']:.1f} esr={r['esr']:.1f}"
                  f"  -> S_B={r['S_B']:.5f} (alt per-side {r['S_B_alt']:.5f})")
            print(f"     RHS={r['rhs']:.4f}  RHS_alt={r['rhs_alt']:.4f}")
        return

    if args.grid:
        print("  --- capture x gDNA x nascent x mature ---")
        for cap in (False, True):
            for g in (0.05, 3.0):
                for na in (20.0, 200.0):
                    for ma in (20.0, 200.0):
                        lab = f"cap={'on' if cap else 'off'} g={g} n={na:.0f} m={ma:.0f}"
                        report(run_cell(f"{cap}_{g}_{na}_{ma}", capture_on=cap, gdna=g,
                                        nascent=na, mature=ma), lab)
        print("\n  --- probe layout (capture ON) ---")
        for lay in ("full", "centre", "junction"):
            report(run_cell(f"pr_{lay}", capture_on=True, gdna=3.0, nascent=200.0, mature=100.0,
                            probes=lay), f"probes={lay}")
    if args.fl_sweep:
        print("\n  --- FL sweep (RNA/gDNA fragment length) ---")
        for rl, gl in ((200, 100), (300, 100), (150, 150), (400, 80)):
            report(run_cell(f"fl_{rl}_{gl}", capture_on=True, gdna=3.0, nascent=200.0, mature=100.0,
                            rna_fl=dict(frag_mean=rl, frag_min=max(50, rl - 100), frag_max=rl + 100),
                            gdna_fl=dict(frag_mean=gl)), f"RNA_FL={rl} gDNA_FL={gl}")


if __name__ == "__main__":
    main()
