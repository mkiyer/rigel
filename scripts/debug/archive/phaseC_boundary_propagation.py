"""Phase C prototype: BOUNDARY propagation — verify the flux physics + local enrichment-matched nascent.

The user's physics for a boundary node (a single genomic position with a left and right side):
  - UNSPLICED fragments propagate THROUGH in BOTH directions -> mass on both sides (gDNA + nascent,
    mature-free, since a contiguous unspliced read cannot span a junction).
  - SPLICED fragments ENTER/EXIT on ONE side, ONE direction -> the junction read's body is on the exon
    side only (mature, the junction's strand).

This prototype, on locus 21, (a) confirms that physics from the accumulator side-flux, and (b) shows the
boundary's per-strand NASCENT crossing density is LOCAL and enrichment-matched (low at intron-facing
boundaries, high at exon-facing boundaries) — the fix for the region-chain prototype's distant-carry
failure (region 224, an intron, wrongly inherited an exon seed's nascent).

Usage:  python scripts/debug/phaseC_boundary_propagation.py [condition]
"""
import dataclasses
import sys

import numpy as np

from rigel.calibration.effective_length import boundary_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.index import TranscriptIndex
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
_SIG_BITS = ((BIT_EXON_POS, "E+"), (BIT_EXON_NEG, "E-"), (BIT_INTRON_POS, "I+"), (BIT_INTRON_NEG, "I-"))


def decode(s):
    return "|".join(tag for bit, tag in _SIG_BITS if s & bit) or "."


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_rnd_capture_on"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    rna_fl_mean = boundary_eff_length(fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)

    sig = np.asarray(ra.signature).astype(np.int64)
    df = idx.region_df
    rid = df["region_id"].to_numpy()
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()

    # For boundary (i, i+1): LEFT side = region i's RIGHT view (sub.right[i]); RIGHT side = region i+1's
    # LEFT view (sub.left[i+1]). Unspliced crosses (both sides); spliced enters one side.
    R = sub.right
    L = sub.left

    def view_row(v, i):
        return (float(v.n_unspliced_pos[i]) + float(v.n_unspliced_neg[i]),
                float(v.n_spliced_sense[i]) + float(v.n_spliced_antisense[i]),
                float(v.n_unspliced_neg[i]), float(v.n_spliced_antisense[i]))

    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165))
    print(f"=== {cond}: locus-21 boundary flux (verify physics: unspliced 2-sided, spliced 1-sided) ===")
    print(f"  kappa={kappa:.3f}  rna_fl_mean={rna_fl_mean:.1f}")
    print(f"{'bnd(i|i+1)':>14}{'Lsig':>6}{'Rsig':>6}{'Luns':>7}{'Runs':>7}{'Lspl':>7}{'Rspl':>7}"
          f"{'spl_1side?':>11}{'rhoNasc_neg':>12}")
    for i in sel[:-1]:
        if i + 1 >= len(sig):
            continue
        l_uns, l_spl, l_uns_neg, l_spl_anti = view_row(R, i)       # left side (region i, right view)
        r_uns, r_spl, r_uns_neg, r_spl_anti = view_row(L, i + 1)    # right side (region i+1, left view)
        # spliced one-sidedness: fraction of spliced on the dominant side
        tot_spl = l_spl + r_spl
        one_sided = max(l_spl, r_spl) / tot_spl if tot_spl > 0 else float("nan")
        # boundary per-strand NASCENT density (crude): unspliced crossing is gDNA+nascent (mature-free);
        # the NEG nascent crossing density ~ neg-unspliced / rna_fl_mean (local, enrichment-matched).
        uns_neg = l_uns_neg + r_uns_neg
        rho_nasc_neg = uns_neg / rna_fl_mean
        print(f"{f'{rid[i]}|{rid[i+1]}':>14}{decode(int(sig[i])):>6}{decode(int(sig[i+1])):>6}"
              f"{l_uns:>7.0f}{r_uns:>7.0f}{l_spl:>7.0f}{r_spl:>7.0f}{one_sided:>11.2f}{rho_nasc_neg:>12.2f}")
    print("\n  spl_1side? ~1.0 confirms spliced enters one side only (the user's physics).")
    print("  rhoNasc_neg: boundary-LOCAL neg crossing density — should be LOW at intron-facing boundaries")
    print("  (e.g. around region 224, the GENE0023 intron) and HIGHER at exon-facing ones (region 226).")


if __name__ == "__main__":
    main()
