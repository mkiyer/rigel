"""Phase-2b: the odds-propagation grid sum-product sweep (`calibration.simplex_sweep`).

(1) Sweep mechanics: order-independence (tree BP); decoupled edges reduce to the per-node local belief.
(2) Toy AMBIG integration: the AMBIG region inherits +odds from its + neighbour and −odds from its −
neighbour → reads RNA when gDNA=0, real gDNA when present.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel.calibration.simplex import _simplex_lattice
from rigel.calibration.simplex_sweep import _log_odds, _sweep_chain


def _toy_psi(m, P, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, 1.0, (m, P))


def test_sweep_order_independent():
    """Forward+backward on a chain is order-independent (a chain is a tree → unique fixed point)."""
    P = _simplex_lattice(12)[0].size
    fp, fn, fg = _simplex_lattice(12)
    lo_pos, lo_neg = _log_odds(fp, fn, fg)
    psi = _toy_psi(5, P, 3)
    qp = np.array([0.2, 0.5, 0.3, 0.4])
    qn = np.array([0.3, 0.2, 0.6, 0.5])
    fwd = _sweep_chain(psi, qp, qn, lo_pos, lo_neg)
    rev = _sweep_chain(psi[::-1], qp[::-1], qn[::-1], lo_pos, lo_neg)
    assert np.allclose(fwd, rev[::-1], atol=1e-9)


def test_sweep_decoupled_reduces_to_local():
    """With all edges decoupled (Q=∞) the edge logφ=0, so the message is a per-node constant and the
    POSTERIOR (softmax of the belief) equals the local softmax(ψ) — no propagation."""
    from scipy.special import logsumexp
    fp, fn, fg = _simplex_lattice(12)
    lo_pos, lo_neg = _log_odds(fp, fn, fg)
    psi = _toy_psi(4, fp.size, 7)
    inf = np.full(3, np.inf)
    belief = _sweep_chain(psi, inf, inf, lo_pos, lo_neg)
    post_b = np.exp(belief - logsumexp(belief, axis=1, keepdims=True))
    post_p = np.exp(psi - logsumexp(psi, axis=1, keepdims=True))
    assert np.allclose(post_b, post_p, atol=1e-9)


# --- toy AMBIG integration ----------------------------------------------------------------------


def _ambig_sweep_fraction(work_dir, *, gdna_abundance, nrna_abundance):
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.calibration.gdna_strand import (
        fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
    )
    from rigel.calibration.density_model import node_gdna_density
    from rigel.calibration.effective_length import (
        boundary_eff_length, boundary_side_eff_length, region_eff_length,
    )
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.signature import TS_AMBIG
    from rigel.calibration.simplex_sweep import deconv_regions_sweep
    from rigel.calibration.strand_balance import fit_strand_balance
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
    from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
    from rigel.splice import SpliceType

    sc = Scenario("ambig_sweep", genome_length=30000, seed=7, work_dir=work_dir)
    sc.add_gene("gA", "+", [{"t_id": "TA", "exons": [(1000, 1500), (4000, 6000)], "abundance": 100}])
    sc.add_gene("gB", "-", [{"t_id": "TB", "exons": [(5000, 7000), (10000, 10500)], "abundance": 100}])
    sc.add_gene("s1", "+", [{"t_id": "S1", "exons": [(12000, 12500), (13500, 14000), (15000, 15500)],
                             "abundance": 120}])
    sc.add_gene("s2", "-", [{"t_id": "S2", "exons": [(17000, 17500), (18500, 19000), (20000, 20500)],
                             "abundance": 120}])
    gd = (GDNAConfig(abundance=gdna_abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
          if gdna_abundance > 0 else None)
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=7),
        gdna_config=gd, nrna_abundance=float(nrna_abundance),
    )
    idx, bam = res.index, str(res.bam_path)
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
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion
    regions = deconv_regions_sweep(sub, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=od_g,
                                   rna_strand_overdispersion=od_r)
    sc.cleanup()
    amb = np.flatnonzero(np.asarray(ra.strand_class) == TS_AMBIG)
    g = np.asarray(regions.gdna_mass)[amb]
    r = np.asarray(regions.rna_mass)[amb]
    return float(g.sum() / max(g.sum() + r.sum(), 1e-9))


def test_sweep_no_false_gdna_from_nascent(tmp_path):
    # gDNA=0 + nascent: the AMBIG inherits high +/− odds from its single-strand exon neighbours → reads RNA.
    frac = _ambig_sweep_fraction(tmp_path / "none", gdna_abundance=0, nrna_abundance=30.0)
    assert frac < 0.25, f"AMBIG gDNA fraction {frac:.3f} too high at gDNA=0+nascent"


def test_sweep_reads_gdna_when_present(tmp_path):
    frac = _ambig_sweep_fraction(tmp_path / "gdna", gdna_abundance=200, nrna_abundance=0.0)
    assert frac > 0.2, f"AMBIG gDNA fraction {frac:.3f} too low with real gDNA present"
