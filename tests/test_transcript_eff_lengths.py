"""The payload's RNA pmf reaches every transcript's effective length in the EM, along the chain
``payload -> build_fl_models -> FLModels.rna_pmf -> FragmentLengthModel.from_pmf ->
compute_all_transcript_eff_lens -> effective_lengths``.

The shrunk pmf is reused for the EM's transcript effective lengths, so a frame error in the anchor
propagates into every transcript row and not only into calibration. The left-hand side comes out of a
full ``run_pipeline`` while the right-hand side is rebuilt from the payload alone, so a pipeline feeding
its effective lengths from anything else, or from a differently-built model, cannot agree with it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

SEED = 7


@pytest.fixture(scope="module")
def scenario(tmp_path_factory):
    work = tmp_path_factory.mktemp("eff_lens")
    sc = Scenario("eff_lens", genome_length=12000, seed=SEED, work_dir=work / "work")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(500, 900), (1500, 1900), (2500, 2900)], "abundance": 80},
            {"t_id": "t2", "exons": [(500, 900), (2500, 2900)], "abundance": 40},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t3", "exons": [(6000, 6400), (7000, 7400)], "abundance": 50}])
    sim = ReadSimConfig(
        frag_mean=200,
        frag_std=40,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.9,
        seed=SEED,
    )
    result = sc.build_oracle(
        n_rna_fragments=2000,
        sim_config=sim,
        gdna_fraction=0.2,
        gdna_config=GDNAConfig(),
        nrna_abundance=15.0,
    )
    yield result
    sc.cleanup()


def _config():
    # `fractional` — no draw at all, so an equality assertion here measures the lengths and nothing else.
    return PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )


def test_every_transcript_eff_length_comes_from_the_payloads_rna_pmf(scenario):
    """The pipeline's per-transcript effective lengths must be **exactly** what the payload's RNA pmf
    produces. Any other source (a stale model, a default fallback) gives a different array.
    """
    from rigel.pipeline import drain_side_buffer, library_fl_models

    config = _config()
    pr = run_pipeline(scenario.bam_path, scenario.index, config=config)
    shipped = np.asarray(pr.estimator.effective_lengths, dtype=np.float64)

    # Rebuilt from the payload alone, by the same route production takes but from an independent scan.
    #
    # The drain is part of that route: the RNA pmf calibration sees is fitted from the drained tally,
    # so an undrained payload gives a different array on every transcript. Production's own helper is
    # called rather than its three steps repeated, so "the same route" cannot quietly stop being true.
    _, strand_models, _, payload = scan_and_buffer(
        str(scenario.bam_path), scenario.index, config.scan
    )
    payload = drain_side_buffer(
        payload, scenario.index, strand_models, seed=config.second_pass_seed
    )
    fl_models = library_fl_models(payload, scenario.index)
    rna_fl = FragmentLengthModel.from_pmf(fl_models.rna_pmf, fl_models.max_size)
    exonic = scenario.index.t_df["length"].values.astype(np.int64)
    expected = rna_fl.compute_all_transcript_eff_lens(exonic)

    assert shipped.shape == expected.shape
    np.testing.assert_allclose(shipped, expected, rtol=0, atol=0)


def test_the_eff_lengths_move_when_the_rna_pmf_moves(scenario):
    """The equality above is worth nothing unless the quantity is actually sensitive to the pmf.

    Perturbation: a constant array satisfies it, so shifting the RNA pmf must change the effective
    lengths.
    """
    from rigel.pipeline import library_fl_models

    _, _, _, payload = scan_and_buffer(str(scenario.bam_path), scenario.index, _config().scan)
    fl_models = library_fl_models(payload, scenario.index)
    exonic = scenario.index.t_df["length"].values.astype(np.int64)

    base = FragmentLengthModel.from_pmf(
        fl_models.rna_pmf, fl_models.max_size
    ).compute_all_transcript_eff_lens(exonic)

    # Shift the pmf 50 bp longer: longer fragments ⇒ shorter effective lengths, on every transcript.
    shifted_pmf = np.roll(np.asarray(fl_models.rna_pmf, dtype=np.float64), 50)
    shifted_pmf[:50] = 0.0
    shifted = FragmentLengthModel.from_pmf(
        shifted_pmf, fl_models.max_size
    ).compute_all_transcript_eff_lens(exonic)

    assert not np.allclose(base, shifted), "effective lengths are insensitive to the RNA pmf"
    assert np.all(shifted <= base + 1e-9), (
        "a longer fragment-length distribution must not lengthen any effective length"
    )


def test_the_n_observations_guard_on_that_path_cannot_fire(scenario):
    """A finding, pinned so it is not mistaken for a live safety net.

    ``transcript_geometry`` refuses a model with ``rna_fl.n_observations <= 0`` (the gate
    below). On the production path ``rna_fl`` always comes from
    :meth:`FragmentLengthModel.from_pmf`, which sets ``_total_weight = 1.0`` — so ``n_observations`` is
    **exactly 1, always**, whatever went into the pmf, and the refusal cannot fire there.

    That is right: an empty RNA pool is EB-shrunk all the way to the anchor by ``build_fl_models``. It is
    pinned because the refusal's *appearance* of an empty-pool check could hide an empty-pool bug.
    """
    from rigel.pipeline import library_fl_models

    _, _, _, payload = scan_and_buffer(str(scenario.bam_path), scenario.index, _config().scan)
    fl_models = library_fl_models(payload, scenario.index)
    rna_fl = FragmentLengthModel.from_pmf(fl_models.rna_pmf, fl_models.max_size)
    assert rna_fl.n_observations == 1

    # …and it is 1 even for a pmf built from a pool with no observations at all.
    flat = np.full(fl_models.max_size + 1, 1.0 / (fl_models.max_size + 1))
    assert FragmentLengthModel.from_pmf(flat, fl_models.max_size).n_observations == 1


def test_an_rna_length_model_with_no_observations_is_refused_not_defaulted(mini_index):
    """A library with no RNA fragment-length observation has no effective length to compute, and the
    pipeline must say so rather than substitute a constant: production reaches this function with a
    pmf-built model that always carries observations, so a fallback here is a magic number nothing
    can exercise. The refusal names the condition."""
    from rigel.pipeline import transcript_geometry

    with pytest.raises(ValueError, match="no RNA fragment-length observation"):
        transcript_geometry(mini_index, FragmentLengthModel(), None, None, None)
