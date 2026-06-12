"""Phase C prototype v2: the #1 resolution — propagate the seedable strand, node-strand resolves the rest.

User's key insight: an AMBIG region has at most 2 strands; one is always seedable (a nested gene's TSS/TES
boundary exposes the OTHER gene's nascent as a 2-component, strand-solvable node; and a gene with
single-strand stretches seeds directly). So: propagate the seedable strand's nascent (enrichment-matched),
subtract it + the matures, and the residual = gDNA + the nested strand's nascent -> the node's OWN strand
resolves it (gDNA symmetric, nascent tilted).

This validates that claim on the HARD case (nrna_rnd, where v1 failed at region 224, an intron that wrongly
inherited an exon seed's nascent). Enrichment-matching here = per-strand nascent density split by the
region's own exon/intron class for that strand (on-target vs off-target under capture).

For locus 21: GENE0023(-) is seedable (single-strand NEG regions); GENE0024(+) is nested. So propagate
neg-nascent (enrichment-matched), subtract neg-nascent + matures, then the residual's strand resolves the
+ nascent vs gDNA via the moment estimate (the strand posterior is the production refinement):
    nasc+ = (U'_pos - U'_neg)/(2k-1);   gDNA = (U'_pos + U'_neg) - nasc+;   fg = gDNA/U.

Usage:  python scripts/debug/phaseC_propagation_v2.py [condition]
"""
import dataclasses
import sys

import numpy as np
import pysam

from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.mature_density import mature_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_INTRON_NEG, TS_AMBIG, TS_NEG, TS_NONE, TS_POS,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import strand_posterior_gdna_frac
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
SC = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}
MIN_SEED = 200.0


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_rnd_capture_on"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    cfg = PipelineConfig()
    ccfg = cfg.calibration
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    reg_el = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bnd_el = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion
    md = mature_density(sub, ra, e_mu, boundary_eff_length(fl.rna_pmf))

    ts = np.asarray(ra.strand_class)
    sig = np.asarray(ra.signature).astype(np.int64)
    up = sub.contained.n_unspliced_pos.astype(float)
    un = sub.contained.n_unspliced_neg.astype(float)
    U = up + un

    # ---- seed neg-nascent density, split by ENRICHMENT CLASS (the region's own GENE0023 exon/intron) ----
    # single-strand NEG seeds: strand solve -> fg -> neg RNA -> neg nascent density; classify on-target
    # (E-) vs off-target (I-) and take the count-weighted mean per class (a crude global-per-class stand-in
    # for the boundary-LOCAL value; production reads it from the bounding boundaries).
    neg_seed = (ts == TS_NEG) & (U >= MIN_SEED)
    im = np.flatnonzero(neg_seed)
    g_q, _ = strand_posterior_gdna_frac(un[im], up[im], kappa, gdna_strand_overdispersion=od_g,
                                        rna_strand_overdispersion=od_r, n_grid=ccfg.n_grid)
    nasc_neg_seed = np.maximum((1 - g_q) * U[im] - md.mature_neg[im], 0.0) / np.maximum(e_mu[im], 1e-9)
    is_exon_neg = (sig[im] & BIT_EXON_NEG) != 0
    w_seed = U[im]
    on_tgt = (np.sum(nasc_neg_seed[is_exon_neg] * w_seed[is_exon_neg])
              / max(np.sum(w_seed[is_exon_neg]), 1e-9))
    off_tgt = (np.sum(nasc_neg_seed[~is_exon_neg] * w_seed[~is_exon_neg])
               / max(np.sum(w_seed[~is_exon_neg]), 1e-9))
    print(f"=== {cond}: v2 (#1 resolution), locus 21  kappa={kappa:.3f} ===")
    print(f"  neg-nascent density seeds: on-target(E-)={on_tgt:.2f}  off-target(I-)={off_tgt:.2f}")

    # ---- AMBIG solve: subtract enrichment-matched neg-nascent + both matures -> residual -> strand ----
    nasc_neg = np.where((sig & BIT_EXON_NEG) != 0, on_tgt,
                        np.where((sig & BIT_INTRON_NEG) != 0, off_tgt, 0.0)) * e_mu
    seeded_pos = (1 - kappa) * (nasc_neg + md.mature_neg) + kappa * md.mature_pos
    seeded_neg = kappa * (nasc_neg + md.mature_neg) + (1 - kappa) * md.mature_pos
    up2 = np.maximum(up - seeded_pos, 0.0)
    un2 = np.maximum(un - seeded_neg, 0.0)
    resid = up2 + un2
    nasc_pos = (up2 - un2) / (2 * kappa - 1)                       # the residual's + nascent (moment est)
    gdna = np.clip(resid - np.abs(nasc_pos), 0.0, None)
    fg = np.clip(np.where(U > 0, gdna / np.maximum(U, 1e-9), np.nan), 0.0, 1.0)

    # ---- oracle ----
    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    og = np.zeros(len(df))
    om = np.zeros(len(df))
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i])
            (og if parse_origin(r.query_name).kind == "gdna" else om)[rid] += 1
    o_fg = np.where(og + om > 0, og / np.maximum(og + om, 1e-9), np.nan)

    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()
    rid_a = df["region_id"].to_numpy()
    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165) & (ts == TS_AMBIG) & (U > 200))
    print(f"{'rid':>4}{'sigE':>6}{'U':>8}{'nascNeg':>8}{'M+':>7}{'M-':>7}{'oracle_fg':>10}{'v2_fg':>7}")
    for i in sel:
        cls = "E-" if (sig[i] & BIT_EXON_NEG) else ("I-" if (sig[i] & BIT_INTRON_NEG) else ".")
        print(f"{rid_a[i]:>4}{cls:>6}{U[i]:>8.0f}{nasc_neg[i]:>8.0f}{md.mature_pos[i]:>7.0f}"
              f"{md.mature_neg[i]:>7.0f}{o_fg[i]:>10.2f}{fg[i]:>7.2f}")
    err = np.abs(fg[sel] - o_fg[sel])
    print(f"\n  AMBIG mean|v2_fg - oracle_fg| = {np.nanmean(err):.3f}  (v1 region-chain gave ~0.5 on nrna_rnd"
          f"; the depleted-intron 224 was the failure)")


if __name__ == "__main__":
    main()
