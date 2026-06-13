"""Phase C prototype: the per-node likelihood deconvolution over the 3-simplex {RNA+, RNA-, gDNA}.

Standalone validation BEFORE touching calibrate (Phase-A discipline). For each AMBIG region of the flagship
locus 21, build the per-node posterior over (r, t) -- r = RNA fraction = 1-f_g, t = tilt = f+ - f- -- from
three likelihoods and report the MAP f_g vs the oracle:

  L_strand(t)   : BetaBinomial(U_pos | U, 1/2 + (kappa-1/2)*t, od)   -- pins the tilt only
  L_mat+(f+)    : soft LOWER bound f+ >= M+/U   (clipped-Gaussian, width = flux Poisson + geometry scatter)
  L_mat-(f-)    : soft LOWER bound f- >= M-/U
  L_gdna(f_g)   : soft LOWER bound f_g >= g_count  (capture-depleted)

f+ = (r+t)/2, f- = (r-t)/2, f_g = 1-r. 2-D grid over (f_g, t) on the valid region; MAP f_g.

Usage:  python scripts/debug/phaseC_node_posterior.py [condition] [geom_cv]
  geom_cv = the geometry-scatter relative width for the mature/gDNA lower bounds (default 0.12, ~Step-0).
"""
import dataclasses
import sys

import numpy as np
import pysam
from scipy.special import gammaln

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.mature_density import mature_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"


def betabinom_logpmf(k, n, mean, od):
    """Beta-Binomial log-pmf with mean in (0,1) and overdispersion od = 1/(a+b+1) in [0,1)."""
    k = np.asarray(k, dtype=np.float64)
    mean = np.clip(np.asarray(mean, dtype=np.float64), 1e-6, 1 - 1e-6)
    if od <= 1e-9:  # binomial limit
        return (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
                + k * np.log(mean) + (n - k) * np.log1p(-mean))
    s = (1.0 - od) / od  # a + b
    a = mean * s
    b = (1.0 - mean) * s
    return (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
            + gammaln(k + a) + gammaln(n - k + b) - gammaln(n + a + b)
            + gammaln(a + b) - gammaln(a) - gammaln(b))


def soft_lower(x, lo, width):
    """Log of a soft one-sided LOWER bound: ~0 for x >= lo, clipped-Gaussian penalty for x < lo."""
    w = max(width, 1e-6)
    d = np.minimum(x - lo, 0.0) / w
    return -0.5 * d * d


def node_map_fg(u_pos, u_neg, kappa, m_pos, m_neg, g_count, *, od, geom_cv, n_grid=120, allow_pos=True,
                allow_neg=True):
    """MAP f_g for one node over the (f_g, t) grid (constraints: f+>=0, f->=0, |t|<=r)."""
    U = u_pos + u_neg
    if U <= 0:
        return float("nan")
    # flux-Poisson + geometry-scatter widths (as FRACTIONS of U), clamped to a small floor
    def width(mass, n_evidence):
        pois = 1.0 / max(n_evidence, 1.0)
        return float(np.sqrt(pois + geom_cv ** 2)) * (mass / U) + 0.02
    w_mp = width(m_pos, m_pos)
    w_mn = width(m_neg, m_neg)
    w_g = width(g_count, g_count) + 0.05  # extra depletion allowance
    fg_grid = np.linspace(0.0, 1.0, n_grid)
    t_grid = np.linspace(-1.0, 1.0, n_grid)
    FG, T = np.meshgrid(fg_grid, t_grid, indexing="ij")
    R = 1.0 - FG
    FP = (R + T) / 2.0
    FN = (R - T) / 2.0
    valid = (FP >= -1e-9) & (FN >= -1e-9)
    if not allow_pos:
        valid &= np.abs(FP) < 1e-9
    if not allow_neg:
        valid &= np.abs(FN) < 1e-9
    p_pos = 0.5 + (kappa - 0.5) * T
    logL = betabinom_logpmf(u_pos, U, p_pos, od)                 # strand: pins t
    logL = logL + soft_lower(FP, m_pos / U, w_mp)                # mature+ lower bound
    logL = logL + soft_lower(FN, m_neg / U, w_mn)                # mature- lower bound
    logL = logL + soft_lower(FG, g_count / U, w_g)               # gDNA lower bound (depleted)
    logL = np.where(valid, logL, -np.inf)
    i = np.unravel_index(np.argmax(logL), logL.shape)
    return float(FG[i])


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    geom_cv = float(sys.argv[2]) if len(sys.argv) > 2 else 0.12
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
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd_raw = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    gstr = fit_gdna_strand_from_substrate(sub, ra, nd_raw, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight)
    od = gstr.gdna_strand_overdispersion
    md = mature_density(sub, ra, e_mu, rna_fl_mean)
    # g_count (depleted gDNA contained fraction → mass): use the raw count density × eff_len
    nd = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    gcount_mass = nd.count_gdna_frac * (sub.contained.n_unspliced_pos + sub.contained.n_unspliced_neg)

    # oracle f_g (contained-unspliced) per region
    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    og = np.zeros(len(df)); om = np.zeros(len(df))
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
            if k == "gdna":
                og[rid] += 1
            elif k == "mrna":
                om[rid] += 1
    o_fg = np.where(og + om > 0, og / np.maximum(og + om, 1e-9), np.nan)

    ts = np.asarray(ra.strand_class)
    s = df["start"].to_numpy(); e = df["end"].to_numpy(); refn = df["ref_name"].to_numpy()
    rid_a = df["region_id"].to_numpy()
    up = sub.contained.n_unspliced_pos.astype(float)
    un = sub.contained.n_unspliced_neg.astype(float)
    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165) & (ts == TS_AMBIG)
                         & (up + un > 50))
    print(f"=== {cond}: per-node likelihood MAP f_g (AMBIG regions, geom_cv={geom_cv}, od={od:.4f}) ===")
    print(f"{'rid':>4}{'U':>8}{'M+':>7}{'M-':>7}{'gcnt':>7}{'oracle_fg':>10}{'MAP_fg':>8}{'prev(splice)':>13}")
    for i in sel:
        fg = node_map_fg(up[i], un[i], kappa, md.mature_pos[i], md.mature_neg[i], gcount_mass[i],
                         od=od, geom_cv=geom_cv)
        print(f"{rid_a[i]:>4}{up[i]+un[i]:>8.0f}{md.mature_pos[i]:>7.0f}{md.mature_neg[i]:>7.0f}"
              f"{gcount_mass[i]:>7.0f}{o_fg[i]:>10.2f}{fg:>8.2f}")
    print("\n  MAP_fg should match oracle_fg (the AMBIG leak regions were stuck at ~0.37 with the splice"
          " fraction; oracle ~0.57-0.60).")


if __name__ == "__main__":
    main()
