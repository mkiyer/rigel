"""propagation.py — the iterative propagation deconvolution, on the toy overlapping-opposite-strand genome.

Locks the AMBIG behaviour via the production solver directly (cf. test_ambig_scenario, which exercises the
current calibrate path): an AMBIG node must read its nascent/mature RNA as RNA (not gDNA), yet read real
gDNA when present. Same toy scenario; the inputs are assembled the way calibrate will call
propagate_regions.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length, boundary_side_eff_length, region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.gdna_strand import (
    fit_gdna_strand_from_substrate, fit_rna_strand_from_substrate, overdispersion_for_beta,
)
from rigel.calibration.propagation import propagate_regions
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType


def _ambig_gdna_fraction(work_dir, *, gdna_abundance: int, nrna_abundance: float) -> float:
    sc = Scenario("ambig_prop", genome_length=30000, seed=7, work_dir=work_dir)
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
    e_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    nd = node_gdna_density(sub, ra, reg_el, fl_mean, need_count_variance=False)
    od_g = fit_gdna_strand_from_substrate(sub, ra, nd, bnd_el, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.gdna_strand_prior_alpha_beta),
        prior_weight=ccfg.gdna_strand_prior_weight).gdna_strand_overdispersion
    od_r = fit_rna_strand_from_substrate(sub, rna_sense_frac=kappa,
        prior_overdispersion=overdispersion_for_beta(ccfg.rna_strand_prior_alpha_beta),
        prior_weight=ccfg.rna_strand_prior_weight).rna_strand_overdispersion
    res_dec = propagate_regions(
        sub, ra, rna_region_eff_len=e_rna, rna_fl_mean=boundary_eff_length(fl.rna_pmf),
        rna_sense_frac=kappa, gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r,
        count_gdna_frac=nd.count_gdna_frac, n_grid=ccfg.n_grid,
    )
    sc.cleanup()
    ambig = np.flatnonzero(np.asarray(ra.strand_class) == TS_AMBIG)
    assert ambig.size >= 1, "the overlapping pair did not form an AMBIG node"
    g = np.asarray(res_dec.gdna_mass)[ambig]
    r = np.asarray(res_dec.rna_mass)[ambig]
    return float(g.sum() / max(g.sum() + r.sum(), 1e-9))


def test_propagation_no_false_gdna_from_nascent(tmp_path):
    # gDNA=0 + nascent: propagation must read the AMBIG node's RNA as RNA, not gDNA.
    frac = _ambig_gdna_fraction(tmp_path / "none", gdna_abundance=0, nrna_abundance=30.0)
    assert frac < 0.15, f"AMBIG gDNA fraction {frac:.3f} too high at gDNA=0+nascent"


def test_propagation_reads_gdna_when_present(tmp_path):
    # gDNA present, no nascent: the AMBIG node reads substantial gDNA.
    frac = _ambig_gdna_fraction(tmp_path / "gdna", gdna_abundance=200, nrna_abundance=0.0)
    assert frac > 0.2, f"AMBIG gDNA fraction {frac:.3f} too low with real gDNA present"
