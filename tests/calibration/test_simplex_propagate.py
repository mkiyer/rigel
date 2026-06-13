"""simplex_propagate.py — forward–backward gDNA-density propagation (increment 3).

Two layers: (1) isolated tests of the RTS smoother (the exact two-sweep BP) — order-independence (the
theory's headline guarantee) and seed-filling; (2) the integration test on the toy overlapping
opposite-strand genome (same scenario as test_propagation), asserting the AMBIG node reads nascent RNA as
RNA yet reads real gDNA when present — now via the simplex pie with the density coupling.
"""

from __future__ import annotations


import numpy as np

from rigel.calibration.simplex_propagate import _rts_smooth


# --- isolated RTS smoother (the two-sweep BP) ----------------------------------------------------


def test_rts_order_independent():
    """The smoothed result is identical whether the chain is read left→right or right→left — the
    tree-BP uniqueness guarantee the user required (non-determinism)."""
    y = np.array([1.0, 0.0, 0.0, 5.0])
    r = np.array([10.0, 0.0, 0.0, 10.0])
    ms, _ = _rts_smooth(y, r, 0.1)
    ms_rev, _ = _rts_smooth(y[::-1], r[::-1], 0.1)
    assert np.allclose(ms, ms_rev[::-1])


def test_rts_fills_between_seeds():
    """An unobserved node between two equal seeds inherits their value (precision-weighted fusion)."""
    y = np.array([2.0, 0.0, 0.0, 2.0])
    r = np.array([100.0, 0.0, 0.0, 100.0])
    ms, ps = _rts_smooth(y, r, 0.001)
    assert np.allclose(ms, 2.0, atol=0.05)
    assert ps[1] > ps[0]  # the unobserved interior is less certain than a seed


def test_rts_observation_dominates_at_low_process_noise():
    """With small process noise a confident observation pins its node (≈ y, low variance)."""
    y = np.array([3.0, 0.0])
    r = np.array([1000.0, 0.0])
    ms, ps = _rts_smooth(y, r, 0.001)
    assert abs(ms[0] - 3.0) < 0.01
    assert ps[0] < ps[1]


def test_rts_per_node_q_matches_scalar_when_constant():
    """A per-node Q array equal to a constant reproduces the scalar-Q result (the array path generalizes
    the scalar path)."""
    y = np.array([1.0, 0.0, 0.0, 4.0])
    r = np.array([20.0, 0.0, 0.0, 20.0])
    ms_s, ps_s = _rts_smooth(y, r, 0.05)
    ms_a, ps_a = _rts_smooth(y, r, np.full(4, 0.05))
    assert np.allclose(ms_s, ms_a) and np.allclose(ps_s, ps_a)


def test_rts_high_q_grows_uncertainty_with_distance():
    """A larger per-hop process variance makes a node far from the only seed MORE uncertain (the seed
    carries less confidently) — the precision-decay-with-distance the coupling encodes."""
    y = np.array([5.0, 0.0, 0.0])
    r = np.array([1000.0, 0.0, 0.0])
    _, low_ps = _rts_smooth(y, r, np.array([0.0, 0.01, 0.01]))
    _, high_ps = _rts_smooth(y, r, np.array([0.0, 100.0, 100.0]))
    assert high_ps[2] > low_ps[2]


def test_rts_two_seeds_q_controls_interior_pull():
    """With two seeds, a larger interior process variance lets each interior node follow its NEAR seed
    more (less averaging toward the far seed)."""
    y = np.array([10.0, 0.0, 0.0, 0.0])
    r = np.array([1000.0, 0.0, 0.0, 1000.0])
    lowq, _ = _rts_smooth(y, r, np.full(4, 0.001))
    highq, _ = _rts_smooth(y, r, np.full(4, 50.0))
    # node 1 is adjacent to the high seed (10); high Q → less pulled toward the far 0-seed → larger.
    assert highq[1] > lowq[1]


# --- integration on the toy AMBIG genome ---------------------------------------------------------


def _ambig_gdna_fraction(work_dir, *, gdna_abundance: int, nrna_abundance: float) -> float:
    """AMBIG contained gDNA fraction via the REAL wired calibrate path (use_propagation=True), so the
    count signal is the splice-upgraded region_count_frac, exactly as production feeds it."""
    import dataclasses as _dc

    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.signature import TS_AMBIG
    from rigel.config import PipelineConfig
    from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
    from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
    from rigel.splice import SpliceType

    sc = Scenario("ambig_simplex", genome_length=30000, seed=7, work_dir=work_dir)
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
    ccfg = _dc.replace(cfg.calibration, use_propagation=True)
    scan = _dc.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    result = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, ccfg)
    sc.cleanup()
    ambig = np.flatnonzero(np.asarray(ra.strand_class) == TS_AMBIG)
    assert ambig.size >= 1, "the overlapping pair did not form an AMBIG node"
    g = np.asarray(result.mass_gdna_contained)[ambig]
    r = np.asarray(result.mass_rna_contained)[ambig]
    return float(g.sum() / max(g.sum() + r.sum(), 1e-9))


def test_simplex_propagate_no_false_gdna_from_nascent(tmp_path):
    # gDNA=0 + nascent: the AMBIG node's RNA must read as RNA, not gDNA. The splice-upgraded count clue
    # (region_splice_gdna_frac) sees the mature/nascent as RNA; the β combine keeps f_g low.
    frac = _ambig_gdna_fraction(tmp_path / "none", gdna_abundance=0, nrna_abundance=30.0)
    assert frac < 0.2, f"AMBIG gDNA fraction {frac:.3f} too high at gDNA=0+nascent"


def test_simplex_propagate_reads_gdna_when_present(tmp_path):
    # gDNA present, no nascent: the AMBIG node reads substantial gDNA.
    frac = _ambig_gdna_fraction(tmp_path / "gdna", gdna_abundance=200, nrna_abundance=0.0)
    assert frac > 0.2, f"AMBIG gDNA fraction {frac:.3f} too low with real gDNA present"
