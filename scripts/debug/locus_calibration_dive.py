#!/usr/bin/env python3
"""TEMP DEBUG — region/boundary-level calibration autopsy of a leaking locus.

Re-runs scan + the acyclic calibrator on a condition BAM, then for a target genomic
window dumps, REGION BY REGION:
  - region type (exonic frac, strand class), coords, size, gDNA FL effective length
  - calibration's deconvolved gDNA vs RNA contained mass + gdna_frac
  - the ORACLE truth: true gDNA vs RNA fragments in that region (+ true gdna_frac)
  - the per-region IPR inputs (gdna_geom_len, g, g^2/geom)
and reconstructs the locus gDNA effective length (IPR) to show why it inflates.

Usage:
  python scripts/debug/locus_calibration_dive.py \
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
from rigel.calibration.effective_length import region_eff_length
from rigel.pipeline import scan_and_buffer, _check_region_payload_alignment
from rigel.splice import SpliceType
from rigel.sim.read_name import parse_origin
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import boundary_eff_length
from rigel.calibration.signature import TS_NEG

_STRAND = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}  # TS_* codes


def _exon_bp_per_region(gtf, ref, starts, ends):
    """Bp of annotated exon overlapping each region (the captured/probe-eligible footprint)."""
    exons = []
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) > 4 and f[2] == "exon" and f[0] == ref:
            exons.append((int(f[3]) - 1, int(f[4])))
    exons.sort()
    out = np.zeros(len(starts))
    for i, (s, e) in enumerate(zip(starts, ends)):
        for es, ee in exons:
            if ee <= s:
                continue
            if es >= e:
                break
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
    stats, strand_models, fl_acc, buffer, payload = scan_and_buffer(
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
    cal = calibrate(payload, ra, strand_models, fl_models.gdna_pmf, CalibrationConfig())
    region_eff_len = region_eff_length(ra.region_size_bp, fl_models.gdna_pmf)

    # Recompute the count clue + strand-only estimate to expose the count-vs-strand tug-of-war.
    sub = CalibrationSubstrate.from_payload(payload, ra)
    fl_mean = boundary_eff_length(fl_models.gdna_pmf)
    nd = node_gdna_density(sub, ra, region_eff_len, fl_mean, gdna_frac=None)
    kappa = cal.rna_sense_frac
    cc = sub.contained
    n_unspl = (cc.n_unspliced_pos + cc.n_unspliced_neg).astype(np.float64)
    ts_all = np.asarray(ra.strand_class)
    sense_all = np.where(ts_all == TS_NEG, cc.n_unspliced_neg, cc.n_unspliced_pos).astype(np.float64)
    sfrac = np.where(n_unspl > 0, sense_all / np.maximum(n_unspl, 1e-9), 0.5)
    denom = 0.5 - kappa
    strand_only_gf = np.clip((sfrac - kappa) / denom, 0.0, 1.0) if abs(denom) > 1e-6 else np.ones_like(sfrac)
    m_unspl = np.asarray(cc.mass_unspliced, dtype=np.float64)
    count_gf = np.clip(np.where(m_unspl > 0, nd.density * region_eff_len / np.maximum(m_unspl, 1e-9), 0.0), 0, 1)

    ref_id = index.ref_name_to_id[args.ref]
    in_win = (ra.ref_id == ref_id) & (ra.end > args.start) & (ra.start < args.end)
    idx = np.where(in_win)[0]
    exon_bp = _exon_bp_per_region(args.gtf, args.ref, ra.start[idx], ra.end[idx])

    # ORACLE truth: true gDNA/RNA fragments per region (assign by fragment midpoint).
    win_starts, win_ends = ra.start[idx], ra.end[idx]
    true_g = np.zeros(len(idx))
    true_r = np.zeros(len(idx))
    with pysam.AlignmentFile(str(args.bam), "rb") as b:
        for r in b:  # name-sorted BAM → full pass + position filter
            if r.is_read2 or r.is_secondary or r.is_supplementary:
                continue
            if r.reference_name != args.ref or r.reference_start is None:
                continue
            mid = (r.reference_start + (r.reference_end or r.reference_start)) // 2
            if not (args.start <= mid < args.end):
                continue
            j = np.searchsorted(win_ends, mid, side="right")
            if j >= len(idx) or not (win_starts[j] <= mid < win_ends[j]):
                continue
            if parse_origin(r.query_name).kind == "gdna":
                true_g[j] += 1
            else:
                true_r[j] += 1

    cg = cal.mass_gdna_contained[idx]
    cr = cal.mass_rna_contained[idx]
    geom = cal.gdna_geom_len[idx]
    print(f"\n=== LOCUS {args.ref}:{args.start:,}-{args.end:,}  ({len(idx)} regions) ===")
    print("cal_gf = joint deconv gdna_frac | count_gf = count-clue prior mean | strand_gf = strand-only estimate | true_gf = oracle")
    print(f"{'region':>6} {'exon_bp':>7} {'str':>5} {'obs':>3} | "
          f"{'count_gf':>8} {'strand_gf':>9} {'cal_gf':>6} {'true_gf':>7} | "
          f"{'cal_gdna':>9} {'true_gdna':>9} {'cnt_ev':>7}")
    for k, i in enumerate(idx):
        cgf = cg[k] / (cg[k] + cr[k]) if (cg[k] + cr[k]) > 0 else float("nan")
        tgf = true_g[k] / (true_g[k] + true_r[k]) if (true_g[k] + true_r[k]) > 0 else float("nan")
        obs = "Y" if nd.region_count_observable[i] else "n"
        print(f"{i:>6} {exon_bp[k]:>7.0f} {_STRAND.get(int(ra.strand_class[i]),'?'):>5} {obs:>3} | "
              f"{count_gf[i]:>8.2f} {strand_only_gf[i]:>9.2f} {cgf:>6.2f} {tgf:>7.2f} | "
              f"{cg[k]:>9.0f} {true_g[k]:>9.0f} {nd.count_evidence[i]:>7.0f}")

    # Locus-level IPR reconstruction (contained only; sides add a little).
    G = cg.sum()
    sq = np.where(geom > 0, cg**2 / np.maximum(geom, 1e-9), 0.0).sum()
    span = geom.sum()
    eff = (G + 1.0) ** 2 / (sq + (2 * G + 1.0) / max(span, 1e-9))
    exon_total = exon_bp.sum()
    print(f"\nIPR reconstruction (contained mass): G(cal gDNA)={G:.0f}  Σg²/geom={sq:.1f}  span(Σgeom)={span:.0f}")
    print(f"  → gDNA eff_len ≈ {min(eff, span):.0f}   vs Σexon_bp={exon_total:.0f}  (ratio efflen/exon={min(eff,span)/max(exon_total,1):.2f})")
    print(f"  true gDNA total in window={true_g.sum():.0f}  cal gDNA contained={G:.0f}  (cal/true={G/max(true_g.sum(),1):.2f})")
    print(f"  true RNA total={true_r.sum():.0f}  cal RNA contained={cr.sum():.0f}")


if __name__ == "__main__":
    main()
