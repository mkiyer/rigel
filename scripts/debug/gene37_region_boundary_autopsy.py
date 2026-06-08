#!/usr/bin/env python3
"""TEMP DEBUG — region + boundary autopsy of one locus's count prior vs strand vs truth.

For a genomic window (default GENE0037 / locus 36), re-runs scan+calibrate and dumps, with
ORACLE truth, the inputs that feed the joint deconvolution:

  REGIONS:   contained mass, count-observability, the swept count-clue density + count_evidence,
             count_gf (count-prior mean), strand sense/antisense + strand_gf, joint cal_gf,
             and TRUE contained gDNA/RNA + the FL of the contained gDNA.
  BOUNDARIES: exon/intron seam, count-observability, crossing mass + conduit weight, deconvolved
             crossing gDNA, and TRUE crossing gDNA/RNA.
  TRANSPORT: per-region gDNA mass BEFORE vs AFTER priors._transport_boundary_flux (does the
             exon gain / the intron lose boundary smear?).
  FL HYPOTHESIS: per exon, are short gDNA fragments CONTAINED (invisible to the count clue) while
             only long ones CROSS a boundary (observable)?

Usage:
  python scripts/debug/gene37_region_boundary_autopsy.py \
    --index <idx> --bam <oracle.bam> --ref chr_syn --start 1659774 --end 1701056 --gtf <genes.gtf>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pysam

from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.priors import _transport_boundary_flux
from rigel.calibration.signature import TS_NEG, TS_NONE
from rigel.pipeline import scan_and_buffer, _check_region_payload_alignment
from rigel.splice import SpliceType
from rigel.sim.read_name import parse_origin

_STRAND = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def _exon_bp(gtf, ref, starts, ends):
    exons = []
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) > 4 and f[2] == "exon" and f[0] == ref:
            exons.append((int(f[3]) - 1, int(f[4])))
    exons.sort()
    out = np.zeros(len(starts))
    for i, (s, e) in enumerate(zip(starts, ends)):
        for es, ee in exons:
            if ee <= s or es >= e:
                continue
            out[i] += min(e, ee) - max(s, es)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--ref", default="chr_syn")
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    ap.add_argument("--gtf", required=True)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    _stats, strand_models, fl_acc, _buf, payload = scan_and_buffer(
        str(args.bam), index, BamScanConfig(sj_strand_tag="auto")
    )
    strand_models.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    _check_region_payload_alignment(ra, payload)
    fl_models = build_fl_models(
        global_counts=fl_acc.global_model.counts,
        rna_counts=fl_acc.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=fl_acc.max_size,
    )
    gdna_pmf = fl_models.gdna_pmf
    cal = calibrate(payload, ra, strand_models, gdna_pmf, CalibrationConfig())

    sub = CalibrationSubstrate.from_payload(payload, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, gdna_pmf)
    fl_mean = boundary_eff_length(gdna_pmf)
    kappa = cal.rna_sense_frac

    # Strand-cleaned gdna_frac exactly as calibrate() computes it (feeds the count clue).
    ts = np.asarray(ra.strand_class)
    cc = sub.contained
    n_unspl = (cc.n_unspliced_pos + cc.n_unspliced_neg).astype(np.float64)
    sense = np.where(ts == TS_NEG, cc.n_unspliced_neg, cc.n_unspliced_pos).astype(np.float64)
    sfrac = np.where(n_unspl > 0, sense / np.maximum(n_unspl, 1e-9), 0.5)
    denom = 0.5 - kappa
    clean_gf = np.clip((sfrac - kappa) / denom, 0, 1) if abs(denom) > 1e-6 else np.ones_like(sfrac)
    clean_gf = np.where(ts == TS_NONE, 1.0, clean_gf)
    nd = node_gdna_density(sub, ra, region_eff_len, fl_mean, gdna_frac=clean_gf)

    # Boundary-transport: per-region gDNA mass before vs after.
    g_before = cal.mass_gdna_contained + cal.mass_gdna_left + cal.mass_gdna_right
    g_after = _transport_boundary_flux(
        cal.mass_gdna_contained, cal.mass_gdna_left, cal.mass_gdna_right,
        ra.region_size_bp, cal.gdna_boundary_len, np.asarray(ra.ref_id),
    )

    ref_id = index.ref_name_to_id[args.ref]
    win = np.where((ra.ref_id == ref_id) & (ra.end > args.start) & (ra.start < args.end))[0]
    exb = _exon_bp(args.gtf, args.ref, ra.start[win], ra.end[win])

    def rtype(i, k):
        if int(ra.strand_class[i]) == 0:
            return "interg"
        return "EXON" if exb[k] > 0 else "intron"

    # ---- ORACLE truth: classify each window fragment as contained-in-region or crossing-boundary.
    ws, we = ra.start[win], ra.end[win]
    bpos = ws[1:]  # boundary positions = shared region edges (contiguous)
    true_cont_g = np.zeros(len(win)); true_cont_r = np.zeros(len(win))
    fl_cont_g = [[] for _ in win]
    cross_g = np.zeros(max(1, len(win) - 1)); cross_r = np.zeros(max(1, len(win) - 1))
    # FL hypothesis per exon: contained vs crossing among gDNA fragments OVERLAPPING the region
    ov_g_contained = np.zeros(len(win)); ov_g_crossing = np.zeros(len(win))
    with pysam.AlignmentFile(str(args.bam), "rb") as b:
        for r in b:
            if r.is_read2 or r.is_secondary or r.is_supplementary or r.reference_name != args.ref:
                continue
            # Use the full TEMPLATE span (paired-end fragment extent), not just read1, so
            # contained-vs-crossing matches the accumulator's fragment-level attribution.
            fl = abs(r.template_length) or ((r.reference_end or r.reference_start) - r.reference_start)
            s = min(r.reference_start, r.next_reference_start) if r.template_length else r.reference_start
            e = s + fl
            if e <= args.start or s >= args.end:
                continue
            is_g = parse_origin(r.query_name).kind == "gdna"
            # contained region
            j = np.searchsorted(we, s, side="right")
            if 0 <= j < len(win) and ws[j] <= s and e <= we[j]:
                if is_g:
                    true_cont_g[j] += 1; fl_cont_g[j].append(fl); ov_g_contained[j] += 1
                else:
                    true_cont_r[j] += 1
            # crossing boundaries (a fragment may straddle one of the internal edges)
            for bi, p in enumerate(bpos):
                if s < p < e:
                    if is_g:
                        cross_g[bi] += 1
                    else:
                        cross_r[bi] += 1
            # FL hypothesis: gDNA fragments overlapping an exon, contained vs crossing its edges
            if is_g:
                k = np.searchsorted(we, (s + e) // 2, side="right")
                if 0 <= k < len(win) and exb[k] > 0:  # midpoint in an exon region
                    if ws[k] <= s and e <= we[k]:
                        pass  # already counted ov_g_contained via j above if fully inside
                    else:
                        ov_g_crossing[k] += 1

    # ---------------- REGION TABLE ----------------
    cg, crna = cal.mass_gdna_contained, cal.mass_rna_contained
    print(f"\n=== {args.ref}:{args.start:,}-{args.end:,}  REGIONS  (kappa={kappa:.3f}) ===")
    print(f"{'reg':>4} {'type':>6} {'eff':>5} {'obs':>3} | {'cont_mass':>9} {'density':>8} {'cnt_ev':>7} "
          f"{'count_gf':>8} {'strand_gf':>9} {'cal_gf':>6} {'true_gf':>7} | {'t_cont_g':>8} {'t_cont_r':>8} {'FLg':>5}")
    for k, i in enumerate(win):
        cgf = cg[i] / (cg[i] + crna[i]) if (cg[i] + crna[i]) > 0 else float("nan")
        tg, tr = true_cont_g[k], true_cont_r[k]
        tgf = tg / (tg + tr) if (tg + tr) > 0 else float("nan")
        cgf_count = min(nd.density[i] * region_eff_len[i] / max(cc.mass_unspliced[i], 1e-9), 1.0)
        flg = np.mean(fl_cont_g[k]) if fl_cont_g[k] else 0.0
        print(f"{i:>4} {rtype(i,k):>6} {region_eff_len[i]:>5.0f} {'Y' if nd.region_count_observable[i] else 'n':>3} | "
              f"{cc.mass_unspliced[i]:>9.0f} {nd.density[i]:>8.4f} {nd.count_evidence[i]:>7.0f} "
              f"{cgf_count:>8.2f} {clean_gf[i]:>9.2f} {cgf:>6.2f} {tgf:>7.2f} | {tg:>8.0f} {tr:>8.0f} {flg:>5.0f}")

    # ---------------- BOUNDARY TABLE ----------------
    print(f"\n=== BOUNDARIES (between consecutive regions) ===")
    print(f"{'L→R reg':>9} {'L→R type':>14} {'obs':>3} | {'cross_mass':>10} {'cross_gdna':>10} "
          f"{'conduit_w':>9} | {'true_cross_g':>12} {'true_cross_r':>12}")
    for k in range(len(win) - 1):
        i, j = win[k], win[k + 1]
        cross_mass = sub.right.mass_unspliced[i] + sub.left.mass_unspliced[j]
        cross_gd = cal.mass_gdna_right[i] + cal.mass_gdna_left[j]
        flux = float(sub.right.n_unspliced[i] + sub.left.n_unspliced[j])
        w = flux / (flux + 1.0)
        obs = "Y" if nd.boundary_count_observable[i] else "n"
        print(f"{str(i)+'→'+str(j):>9} {rtype(i,k)+'→'+rtype(j,k+1):>14} {obs:>3} | "
              f"{cross_mass:>10.0f} {cross_gd:>10.0f} {w:>9.3f} | {cross_g[k]:>12.0f} {cross_r[k]:>12.0f}")

    # ---------------- TRANSPORT ----------------
    print(f"\n=== BOUNDARY TRANSPORT — per-region gDNA mass before → after _transport_boundary_flux ===")
    print(f"{'reg':>4} {'type':>6} | {'gdna_before':>11} {'gdna_after':>10} {'Δ':>8} | {'true_total_g(cont)':>18}")
    for k, i in enumerate(win):
        print(f"{i:>4} {rtype(i,k):>6} | {g_before[i]:>11.0f} {g_after[i]:>10.0f} "
              f"{g_after[i]-g_before[i]:>+8.0f} | {true_cont_g[k]:>18.0f}")

    # ---------------- FL HYPOTHESIS ----------------
    print(f"\n=== FL HYPOTHESIS (exon regions): are SHORT gDNA fragments contained (invisible)? ===")
    print(f"  gDNA FL pmf mean={float(np.dot(np.arange(gdna_pmf.size), gdna_pmf)):.0f}")
    print(f"{'reg':>4} {'exon_bp':>7} {'reg_span':>8} | {'g_contained':>11} {'g_crossing':>10} {'%contained':>10} {'FL_contained':>12}")
    for k, i in enumerate(win):
        if exb[k] <= 0:
            continue
        nc, nx = ov_g_contained[k], ov_g_crossing[k]
        pct = 100 * nc / (nc + nx) if (nc + nx) else 0.0
        flg = np.mean(fl_cont_g[k]) if fl_cont_g[k] else 0.0
        print(f"{i:>4} {exb[k]:>7.0f} {ra.region_size_bp[i]:>8.0f} | {nc:>11.0f} {nx:>10.0f} {pct:>9.0f}% {flg:>12.0f}")


if __name__ == "__main__":
    main()
