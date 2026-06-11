"""Per-region/boundary calibration introspection for the flagship leak (Issue #3).

Re-runs scan + the calibrate() sub-steps and dumps, for each region overlapping the silent
gDNA-only locus 81 (GENE0082, chr_syn:3977582-4017099, all fragments are true gDNA so every node's
true gDNA fraction is 1.0), what the STRAND clue says, what the COUNT clue says, and what the JOINT
deconvolution decided — to find WHY the per-region gDNA is under-called (which then under-contracts
the eff_len and under-sets the count prior).
"""

import numpy as np
from dataclasses import replace as _replace

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.count_dispersion import fit_gdna_count_overdispersion, effective_count
from rigel.calibration.joint_deconv import deconv_regions, deconv_sides, boundary_side_seeds
from rigel.calibration.effective_length import (
    region_eff_length, boundary_side_eff_length, boundary_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import TS_AMBIG, TS_POS, TS_NEG, TS_NONE
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
BAM = f"{SUITE}/gdna_gdna400_ss_0.99_nrna_none_capture_on/sim_oracle.bam"
REF, LO, HI = "chr_syn", 3977582, 4017099

index = TranscriptIndex.load(f"{SUITE}/rigel_index")
cfg = PipelineConfig()
scan = cfg.scan
scan = _replace(scan, sj_strand_tag=_native_detect_sj_tag(BAM))
stats, strand_models, flm, buffer, payload = scan_and_buffer(BAM, index, scan)
strand_models.finalize()

ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
fl_models = build_fl_models(
    global_counts=flm.global_model.counts,
    rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
    gdna_counts=gdna_fl_mass(payload),
    max_size=flm.max_size,
)
gpmf = fl_models.gdna_pmf
cc = cfg.calibration

sub = CalibrationSubstrate.from_payload(payload, ra)
region_eff_len = region_eff_length(ra.region_size_bp, gpmf)
boundary_eff_len = boundary_side_eff_length(gpmf, ra.region_size_bp)
fl_mean = boundary_eff_length(gpmf)
kappa = float(fit_strand_balance(strand_models).rna_sense_frac)
nd = node_gdna_density(sub, ra, region_eff_len, fl_mean, rna_sense_frac=kappa)
gd = fit_gdna_strand_from_substrate(
    sub, ra, nd, boundary_eff_len, rna_sense_frac=kappa,
    prior_overdispersion=overdispersion_for_beta(cc.gdna_strand_prior_alpha_beta),
    prior_weight=cc.gdna_strand_prior_weight,
)
rs = fit_rna_strand_from_substrate(
    sub, rna_sense_frac=kappa,
    prior_overdispersion=overdispersion_for_beta(cc.rna_strand_prior_alpha_beta),
    prior_weight=cc.rna_strand_prior_weight,
)
ts = np.asarray(ra.strand_class)
cseed = nd.region_count_observable & (ts != TS_AMBIG)
ccount = (nd.density * region_eff_len)[cseed]
clen = region_eff_len[cseed]
_bs, b_total, b_weight = boundary_side_seeds(sub, ra, nd, boundary_eff_len, rna_sense_frac=kappa)
xcount = b_weight * b_total
xlen = np.full(xcount.shape, fl_mean)
cd = fit_gdna_count_overdispersion(
    ccount, clen, xcount, xlen,
    prior_alpha=cc.count_overdispersion_prior, prior_weight=cc.count_overdispersion_prior_weight,
)
region_alpha = np.where(nd.region_count_observable, cd.alpha_contained, cd.alpha_crossing)
count_evidence = effective_count(nd.count_support, region_alpha)
regions = deconv_regions(
    sub, ra, nd, count_evidence, rna_sense_frac=kappa,
    gdna_strand_overdispersion=gd.gdna_strand_overdispersion,
    rna_strand_overdispersion=rs.rna_strand_overdispersion,
    deconv_quantile=cc.gdna_deconv_quantile, n_grid=cc.n_grid,
)
left, right = deconv_sides(
    sub, ra, nd, boundary_eff_len, rna_sense_frac=kappa, alpha_crossing=cd.alpha_crossing,
    gdna_strand_overdispersion=gd.gdna_strand_overdispersion,
    rna_strand_overdispersion=rs.rna_strand_overdispersion,
    deconv_quantile=cc.gdna_deconv_quantile, n_grid=cc.n_grid,
)

print(f"kappa(rna_sense_frac)={kappa:.4f}  fl_mean={fl_mean:.1f}  "
      f"od_gdna={gd.gdna_strand_overdispersion:.4g} od_rna={rs.rna_strand_overdispersion:.4g}")
print(f"alpha_contained={cd.alpha_contained:.3f} alpha_crossing={cd.alpha_crossing:.3f}")
TSNAME = {TS_NONE: "NONE", TS_POS: "POS", TS_NEG: "NEG", TS_AMBIG: "AMBIG"}

rid = index.ref_name_to_id[REF]
mask = (np.asarray(ra.ref_id) == rid) & (np.asarray(ra.end) > LO) & (np.asarray(ra.start) < HI)
idxs = np.where(mask)[0]
c = sub.contained
cpos = np.asarray(c.n_unspliced_pos, float); cneg = np.asarray(c.n_unspliced_neg, float)
cmass = np.asarray(c.mass_unspliced, float)

print(f"\n{len(idxs)} regions in locus 81 span. CONTAINED node per region "
      "(true gDNA frac = 1.0 everywhere — silent locus):")
hdr = ("idx", "start", "len", "ts", "obs", "pos", "neg", "sfrac", "dens",
       "cgfrac", "csupp", "cEvid", "gdna_f", "gMass", "rMass")
print("%6s %8s %6s %5s %3s %8s %8s %6s %7s %6s %8s %7s %6s %9s %9s" % hdr)
for i in idxs:
    sden = (cpos[i] + cneg[i])
    sfrac = cpos[i] / sden if sden > 0 else float("nan")
    gf = regions.gdna_frac[i]
    print("%6d %8d %6d %5s %3d %8.0f %8.0f %6.3f %7.3f %6.3f %8.0f %7.2f %6.3f %9.0f %9.0f" % (
        i, int(ra.start[i]), int(ra.end[i] - ra.start[i]), TSNAME.get(int(ts[i]), "?"),
        int(nd.region_count_observable[i]), cpos[i], cneg[i], sfrac, nd.density[i],
        nd.count_gdna_frac[i], nd.count_support[i], count_evidence[i], gf,
        regions.gdna_mass[i], regions.rna_mass[i],
    ))

# Aggregate over the locus's contained nodes
gm = regions.gdna_mass[idxs].sum() + left.gdna_mass[idxs].sum() + right.gdna_mass[idxs].sum()
rm = regions.rna_mass[idxs].sum() + left.rna_mass[idxs].sum() + right.rna_mass[idxs].sum()
print(f"\nlocus-81 deconv totals: gDNA={gm:.0f} RNA={rm:.0f} gdna_frac={gm/(gm+rm):.4f} "
      f"(true=1.0; RNA here is pure leak)")
print(f"  contained: gDNA={regions.gdna_mass[idxs].sum():.0f} RNA={regions.rna_mass[idxs].sum():.0f}")
print(f"  left side: gDNA={left.gdna_mass[idxs].sum():.0f} RNA={left.rna_mass[idxs].sum():.0f}")
print(f"  right side:gDNA={right.gdna_mass[idxs].sum():.0f} RNA={right.rna_mass[idxs].sum():.0f}")

# Boundary SIDES per region (left = right side of boundary (r-1,r); right = left side of (r,r+1)).
for name, view, dec in (("LEFT", sub.left, left), ("RIGHT", sub.right, right)):
    pos = np.asarray(view.n_unspliced_pos, float); neg = np.asarray(view.n_unspliced_neg, float)
    print(f"\n{name} boundary-side nodes (true gDNA frac = 1.0):")
    print("%6s %8s %8s %8s %6s %7s %9s %9s" % (
        "idx", "pos", "neg", "n_side", "sfrac", "gdna_f", "gMass", "rMass"))
    for i in idxs:
        n = pos[i] + neg[i]
        if n <= 0:
            continue
        sf = pos[i] / n
        print("%6d %8.0f %8.0f %8.0f %6.3f %7.3f %9.0f %9.0f" % (
            i, pos[i], neg[i], n, sf, dec.gdna_frac[i], dec.gdna_mass[i], dec.rna_mass[i]))

# Confirm the leaky edges face intergenic + recompute strand_observable per side, then tally
# the DATASET-WIDE leak from count-observable-but-NOT-strand-observable boundary sides whose own
# sfrac says gDNA (|sfrac-0.5|<0.05) — i.e. clearly-gDNA crossing flux defaulting to Jeffreys.
from rigel.calibration.joint_deconv import _left_right_neighbors
print("\nneighbors of the leaky edges:")
for i in (868, 884):
    print(f"  region {i}: ts={TSNAME.get(int(ts[i]))}  "
          f"left_nbr(r-1={i-1})_ts={TSNAME.get(int(ts[i-1]))}  right_nbr(r+1={i+1})_ts={TSNAME.get(int(ts[i+1]))}")

l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
    ts, np.asarray(ra.ref_id), nd.boundary_count_observable)
left_strand_obs = (l_same & (ts == TS_POS) & (ts_prev == TS_POS)) | (l_same & (ts == TS_NEG) & (ts_prev == TS_NEG))
right_strand_obs = (r_same & (ts == TS_POS) & (ts_next == TS_POS)) | (r_same & (ts == TS_NEG) & (ts_next == TS_NEG))

def side_sfrac(view):
    p = np.asarray(view.n_unspliced_pos, float); ng = np.asarray(view.n_unspliced_neg, float)
    n = p + ng
    return np.where(n > 0, p / np.maximum(n, 1), np.nan), n

print("\nDATASET-WIDE boundary-side leak tally (rna_mass on sides whose own sfrac∈[0.45,0.55] ⇒ gDNA):")
tot_leak_strandblind = 0.0; tot_leak_strandobs = 0.0; n_blind = 0
for name, view, dec, sobs, obs in (
    ("left", sub.left, left, left_strand_obs, np.r_[False, nd.boundary_count_observable[:-1]]),
    ("right", sub.right, right, right_strand_obs, nd.boundary_count_observable),
):
    sf, n = side_sfrac(view)
    gdna_like = (np.abs(sf - 0.5) < 0.05) & (n > 50)   # clearly-gDNA by its own strand, non-trivial
    blind = gdna_like & (~sobs)        # strand likelihood skipped on a clearly-gDNA side
    seen = gdna_like & sobs
    tot_leak_strandblind += dec.rna_mass[blind].sum()
    tot_leak_strandobs += dec.rna_mass[seen].sum()
    n_blind += int(blind.sum())
print(f"  strand-BLIND clearly-gDNA sides: n={n_blind}  total rna_mass(leak)={tot_leak_strandblind:.0f}")
print(f"  strand-OBSERVABLE clearly-gDNA sides: total rna_mass(leak)={tot_leak_strandobs:.0f}")
print(f"  [dataset in-locus net leak was ~241,000]")
