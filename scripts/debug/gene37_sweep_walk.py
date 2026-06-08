#!/usr/bin/env python3
"""TEMP DEBUG — node-by-node walk of the Phase-1 density SWEEP across one locus.

Walks region <-> boundary <-> region ... in genomic order and prints, per node, the raw
fractional-accumulator counts (pos/neg unspliced, sense/antisense spliced), the signature,
the DIRECT gDNA observation that feeds the sweep (BEFORE), and the SWEPT density the sweep
produces (AFTER). Boundaries show the conduit weight w = flux/(flux+1).

This exposes density_model.node_gdna_density's "alternating region<->boundary density sweep"
(Phase 1, BEFORE the joint deconvolution — it produces the count-prior density).

Usage:
  python scripts/debug/gene37_sweep_walk.py \
    --index <idx> --bam <oracle.bam> --ref chr_syn --start 1659774 --end 1701056
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.calibration.density_model import node_gdna_density, count_observable_masks
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.signature import TS_NEG, TS_NONE, coarse_type_array
from rigel.pipeline import scan_and_buffer, _check_region_payload_alignment
from rigel.splice import SpliceType

_TYPE = {0: "interg", 1: "intron", 2: "EXON"}


def sig_label(s: int) -> str:
    p = []
    if s & 0x2:
        p.append("ex+")
    if s & 0x1:
        p.append("ex-")
    if s & 0x8:
        p.append("in+")
    if s & 0x4:
        p.append("in-")
    return "|".join(p) if p else "intergenic"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--ref", default="chr_syn")
    ap.add_argument("--start", required=True, type=int)
    ap.add_argument("--end", required=True, type=int)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    _s, strand_models, fl_acc, _b, payload = scan_and_buffer(
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
    gpmf = fl_models.gdna_pmf
    cal = calibrate(payload, ra, strand_models, gpmf, CalibrationConfig())
    sub = CalibrationSubstrate.from_payload(payload, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, gpmf)
    fl_mean = boundary_eff_length(gpmf)
    kappa = cal.rna_sense_frac

    ts = np.asarray(ra.strand_class)
    cc = sub.contained
    n_unspl = (cc.n_unspliced_pos + cc.n_unspliced_neg).astype(np.float64)
    sense = np.where(ts == TS_NEG, cc.n_unspliced_neg, cc.n_unspliced_pos).astype(np.float64)
    sfrac = np.where(n_unspl > 0, sense / np.maximum(n_unspl, 1e-9), 0.5)
    denom = 0.5 - kappa
    clean_gf = np.clip((sfrac - kappa) / denom, 0, 1) if abs(denom) > 1e-6 else np.ones_like(sfrac)
    clean_gf = np.where(ts == TS_NONE, 1.0, clean_gf)
    nd = node_gdna_density(sub, ra, region_eff_len, fl_mean, gdna_frac=clean_gf)
    reg_obs, bnd_obs = count_observable_masks(np.asarray(ra.signature), np.asarray(ra.ref_id))
    ctype = coarse_type_array(np.asarray(ra.signature))

    # DIRECT observations feeding the sweep (BEFORE):
    #   region: strand-cleaned contained gDNA mass (observable regions only)
    #   boundary: raw crossing unspliced mass (observable boundaries only; assumed gDNA)
    reg_mass = np.where(reg_obs, clean_gf * cc.mass_unspliced, 0.0)

    ref_id = index.ref_name_to_id[args.ref]
    win = np.where((ra.ref_id == ref_id) & (ra.end > args.start) & (ra.start < args.end))[0]

    print(f"\n=== SWEEP WALK {args.ref}:{args.start:,}-{args.end:,}  (kappa={kappa:.3f}, fl_mean={fl_mean:.0f}) ===")
    print("counts shown as unspl(+/-) / spl(sense/anti). DIRECT = sweep input; SWEPT density = sweep output.\n")
    for n, i in enumerate(win):
        # REGION node
        up, un = int(cc.n_unspliced_pos[i]), int(cc.n_unspliced_neg[i])
        sp, sa = int(cc.n_spliced_sense[i]), int(cc.n_spliced_antisense[i])
        direct = f"{reg_mass[i]:>8.0f}" if reg_obs[i] else "  (impute)"
        sweptg = nd.density[i] * region_eff_len[i]
        print(f"R{i:<4} {sig_label(int(ra.signature[i])):<9} {_TYPE[ctype[i]]:<6} "
              f"obs={'Y' if reg_obs[i] else 'n'} eff={region_eff_len[i]:>6.0f} "
              f"| unspl({up:>6},{un:>6}) spl({sp:>5},{sa:>5}) mass={cc.mass_unspliced[i]:>8.0f} "
              f"| DIRECT gdna={direct} | SWEPT dens={nd.density[i]:.4f} swept_gdna={sweptg:>7.0f}")
        if n == len(win) - 1:
            break
        j = win[n + 1]
        # BOUNDARY node between i and j (boundary index = i; describes seam i/i+1)
        rup, run_ = int(sub.right.n_unspliced_pos[i]), int(sub.right.n_unspliced_neg[i])
        lup, lun = int(sub.left.n_unspliced_pos[j]), int(sub.left.n_unspliced_neg[j])
        cross_mass = sub.right.mass_unspliced[i] + sub.left.mass_unspliced[j]
        flux = float(sub.right.n_unspliced[i] + sub.left.n_unspliced[j])
        w = flux / (flux + 1.0)
        bdirect = f"{cross_mass:>8.0f}" if bnd_obs[i] else "  (shared-exon: not obs)"
        print(f"  └B{i:<3} {_TYPE[ctype[i]]}->{_TYPE[ctype[j]]:<7} obs={'Y' if bnd_obs[i] else 'n'} "
              f"Rside({rup:>5},{run_:>5}) Lside({lup:>5},{lun:>5}) flux={flux:>7.0f} "
              f"conduit_w={w:.4f} | DIRECT gdna(bnd_mass)={bdirect}")
    print(f"\nGLOBAL seed density (Σ observable gdna / Σ observable len) ≈ {nd.density[win].mean():.4f} "
          f"— note every region's SWEPT density collapses to ~this value.")


if __name__ == "__main__":
    main()
