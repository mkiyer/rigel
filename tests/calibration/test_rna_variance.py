"""Phase-2a: the splice-junction-pair RNA var~mean (`calibration.rna_variance`).

Validates that internal exons (flanked by two same-strand splice junctions) become var~mean observations
and the LOESS fits a sensible, non-negative σ²_RNA — on a synthetic genome with enough multi-exon genes
to exceed the LOESS minimum.
"""

from __future__ import annotations

import dataclasses

import numpy as np


def test_rna_spliced_variance_fits_from_junction_pairs(tmp_path):
    from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.rna_variance import rna_spliced_variance
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
    from rigel.sim import ReadSimConfig, Scenario
    from rigel.splice import SpliceType

    sc = Scenario("rna_var", genome_length=120000, seed=11, work_dir=tmp_path)
    # 8 three-exon + genes (each has one INTERNAL exon flanked by two + junctions → a var~mean pair),
    # with a spread of abundances so the fit sees a range of mean RNA densities.
    for k in range(8):
        base = 2000 + k * 13000
        ab = 40 + 30 * k  # varying expression → a mean-density range
        sc.add_gene(f"g{k}", "+", [{"t_id": f"T{k}", "exons": [
            (base, base + 500), (base + 1500, base + 2200), (base + 3000, base + 3500)], "abundance": ab}])
    res = sc.build_oracle(
        n_fragments=12000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=None, nrna_abundance=0.0,
    )
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    reg_el_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    fl_mean_rna = boundary_eff_length(fl.rna_pmf)
    rv = rna_spliced_variance(sub, ra, reg_el_rna, fl_mean_rna)

    # internal exons gave same-strand junction-pair observations
    assert len(rv.fit_mu) >= 5, f"too few splice-pair fit points: {len(rv.fit_mu)}"
    assert np.all(rv.fit_strand[np.isfinite(rv.fit_mu)] == 1)  # all + genes
    assert np.all(rv.fit_var[np.isfinite(rv.fit_var)] >= 0.0)
    # the LOESS evaluated σ²_RNA at + exons, non-negative; no − transcripts ⇒ neg curve all NaN
    finite_pos = np.isfinite(rv.sigma2_pos)
    assert finite_pos.sum() > 0
    assert np.all(rv.sigma2_pos[finite_pos] >= 0.0)
    assert not np.isfinite(rv.sigma2_neg).any()
