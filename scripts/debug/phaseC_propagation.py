"""Phase C prototype: the iterative PROPAGATION deconvolution (region-chain core), validated on locus 21.

Architecture (docs/calibration/phaseC_ambig_resolution_plan.md rev 3): deconvolution is a graph of nodes
solved from SEEDS by propagation to convergence. This prototype validates the CORE — that propagating the
per-strand NASCENT density from single-strand seeds resolves the AMBIG regions as a residual — before
building the full region+boundary `propagation.py`.

Per region, partition unspliced U into {f+, f-, fg}:
  SEED  NONE (intergenic): fg=1.
  SEED  POS/NEG single-strand: the strand posterior gives fg; RNA_s=(1-fg)*U; nascent_s = RNA_s - mature_s;
        emit per-bp nascent density rho_nasc_s = nascent_s / eff_len_rna.
  AMBIG: carry rho_nasc_s from single-strand neighbours WITHIN the strand-s gene body (E_s|I_s run); then
        RNA_s = rho_nasc_s*eff_len_rna + mature_s;  fg = 1 - (RNA+ + RNA-)/U.
        (gDNA = U - RNA = the region's OWN enriched count minus non-depleted RNA -> the depletion fix.)

Iterate sweeps until fg stable (here: seeds -> carry -> AMBIG converges in 1-2 sweeps for nrna_none).

Usage:  python scripts/debug/phaseC_propagation.py [condition]
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
from rigel.calibration.run_fill import runfill_bidirectional
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS, TS_AMBIG, TS_NEG, TS_NONE, TS_POS,
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
POS_BODY = BIT_EXON_POS | BIT_INTRON_POS
NEG_BODY = BIT_EXON_NEG | BIT_INTRON_NEG


def _body_segments(in_body, ref_id):
    """Segment id constant within each contiguous in-body run (per ref) — for strand-gated run-fill."""
    in_body = np.asarray(in_body, dtype=bool)
    ref = np.asarray(ref_id)
    r = in_body.shape[0]
    brk = np.zeros(r, dtype=bool)
    brk[1:] = (in_body[1:] != in_body[:-1]) | (ref[1:] != ref[:-1])
    return np.cumsum(brk)


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
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
    rna_fl_mean = boundary_eff_length(fl.rna_pmf)
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)  # RNA contained eff-len
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd_raw = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd_raw, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion
    md = mature_density(sub, ra, e_mu, rna_fl_mean)

    ts = np.asarray(ra.strand_class)
    sig = np.asarray(ra.signature).astype(np.int64)
    ref_id = np.asarray(ra.ref_id)
    up = sub.contained.n_unspliced_pos.astype(float)
    un = sub.contained.n_unspliced_neg.astype(float)
    U = up + un
    R = len(ts)

    # ---- SEEDS: single-strand regions -> fg via strand posterior -> per-strand nascent density ----
    rho_nasc_pos = np.full(R, np.nan)
    rho_nasc_neg = np.full(R, np.nan)
    fg = np.full(R, np.nan)
    for cls, sense_arr, rho in ((TS_POS, up, rho_nasc_pos), (TS_NEG, un, rho_nasc_neg)):
        m = (ts == cls) & (U > 0)
        idx_m = np.flatnonzero(m)
        if idx_m.size:
            g_q, _ = strand_posterior_gdna_frac(sense_arr[idx_m], (U - sense_arr)[idx_m], kappa,
                gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r, n_grid=ccfg.n_grid)
            fg[idx_m] = g_q
            rna_s = (1.0 - g_q) * U[idx_m]
            mat_s = (md.mature_pos if cls == TS_POS else md.mature_neg)[idx_m]
            nasc = np.maximum(rna_s - mat_s, 0.0)
            rho[idx_m] = np.where(e_mu[idx_m] > 0, nasc / np.maximum(e_mu[idx_m], 1e-9), 0.0)
    fg[ts == TS_NONE] = 1.0
    rho_nasc_pos[ts == TS_NONE] = 0.0
    rho_nasc_neg[ts == TS_NONE] = 0.0

    # ---- RELIABILITY FILTER: a tiny-mass seed cannot estimate nascent (its strand is noise) and would
    # poison the nearest-neighbour carry. Drop low-count seeds before the fill so the SUBSTANTIAL seeds
    # propagate. (Prototype uses a count threshold; production uses a precision-weighted carry — the seed
    # strand info I=N·(2κ−1)² is the weight, no magic number.) ----
    MIN_SEED = 200.0
    rho_nasc_pos = np.where(U >= MIN_SEED, rho_nasc_pos, np.nan)
    rho_nasc_neg = np.where(U >= MIN_SEED, rho_nasc_neg, np.nan)
    rho_nasc_pos[ts == TS_NONE] = 0.0  # intergenic: genuinely zero nascent (a reliable seed)
    rho_nasc_neg[ts == TS_NONE] = 0.0

    # ---- PROPAGATE: carry rho_nasc_s within the strand-s gene body (E_s|I_s contiguous run) ----
    seg_pos = _body_segments((sig & POS_BODY) != 0, ref_id)
    seg_neg = _body_segments((sig & NEG_BODY) != 0, ref_id)
    rho_nasc_pos = runfill_bidirectional(rho_nasc_pos, seg_pos)
    rho_nasc_neg = runfill_bidirectional(rho_nasc_neg, seg_neg)
    rho_nasc_pos = np.where(((sig & POS_BODY) != 0) & np.isfinite(rho_nasc_pos), rho_nasc_pos, 0.0)
    rho_nasc_neg = np.where(((sig & NEG_BODY) != 0) & np.isfinite(rho_nasc_neg), rho_nasc_neg, 0.0)

    # ---- AMBIG solve: RNA_s = nascent_s (carried) + mature_s; fg = 1 - (RNA+ + RNA-)/U ----
    amb = (ts == TS_AMBIG) & (U > 0)
    rna_p = rho_nasc_pos * e_mu + md.mature_pos
    rna_n = rho_nasc_neg * e_mu + md.mature_neg
    fg_amb = np.clip(1.0 - (rna_p + rna_n) / np.maximum(U, 1e-9), 0.0, 1.0)
    fg = np.where(amb, fg_amb, fg)

    # ---- oracle fg (contained-unspliced) ----
    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    og = np.zeros(R)
    om = np.zeros(R)
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
            k = parse_origin(r.query_name).kind
            (og if k == "gdna" else om if k == "mrna" else om)[rid] += 1
    o_fg = np.where(og + om > 0, og / np.maximum(og + om, 1e-9), np.nan)

    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()
    rid_a = df["region_id"].to_numpy()
    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165)
                         & ((sig & (BIT_EXON_POS | BIT_EXON_NEG | BIT_INTRON_POS | BIT_INTRON_NEG)) != 0)
                         & (U > 50))
    print(f"=== {cond}: propagation prototype, locus 21 (kappa={kappa:.3f}) ===")
    print(f"{'rid':>4}{'cls':>6}{'U':>8}{'M+':>7}{'M-':>7}{'rhoN+':>7}{'rhoN-':>7}{'oracle_fg':>10}{'prop_fg':>8}")
    for i in sel:
        print(f"{rid_a[i]:>4}{SC[ts[i]]:>6}{U[i]:>8.0f}{md.mature_pos[i]:>7.0f}{md.mature_neg[i]:>7.0f}"
              f"{rho_nasc_pos[i]:>7.2f}{rho_nasc_neg[i]:>7.2f}{o_fg[i]:>10.2f}{fg[i]:>8.2f}")
    amb_sel = sel[ts[sel] == TS_AMBIG]
    err = np.abs(fg[amb_sel] - o_fg[amb_sel])
    print(f"\n  AMBIG regions: mean|prop_fg - oracle_fg| = {np.nanmean(err):.3f}  (per-node prototype gave"
          f" ~0.27-0.37 vs oracle 0.57-0.98 -> large error)")


if __name__ == "__main__":
    main()
