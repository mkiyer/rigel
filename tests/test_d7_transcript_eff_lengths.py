"""⭐ C2.5 / D7 — the corrected RNA pmf reaches EVERY TRANSCRIPT'S EFFECTIVE LENGTH in the EM.

``docs/FRAGMENT_LENGTH_AUDIT.md`` D7: *"The shrunk pmf is reused for the EM's transcript effective
lengths, so a frame error in the anchor propagates into every transcript's effective length — not only
into calibration."* The audit is explicit that this must be **checked, not assumed** — it is the path
with the largest blast radius and the least visibility, and nothing was asserting it.

The chain under test::

    payload  ──build_fl_models──▶  FLModels.rna_pmf
                                        │  FragmentLengthModel.from_pmf
                                        ▼
                                     rna_fl  ──compute_all_transcript_eff_lens──▶  effective_lengths
                                                                                        │
                                                                            every transcript row in the EM

⚠ Re-deriving the expected value through the *same* helper would be ``CARRY_FORWARD.md`` §3 trap 1 —
a check that cannot fail. What makes this a real gate is that the left-hand side comes out of a full
``run_pipeline`` and the right-hand side is rebuilt from the PAYLOAD ALONE: if the pipeline fed its
effective lengths from anything else, or from a differently-built model, the two cannot agree.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario

SEED = 7


@pytest.fixture(scope="module")
def scenario(tmp_path_factory):
    work = tmp_path_factory.mktemp("d7")
    sc = Scenario("d7", genome_length=12000, seed=SEED, work_dir=work / "work")
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
        n_fragments=2000, sim_config=sim, gdna_fraction=0.2, nrna_abundance=15.0
    )
    yield result
    sc.cleanup()


def _config():
    # ⚠ `map` — the DETERMINISTIC assignment mode. The EM samples the posterior by default, and an
    # equality assertion under `sample` would be measuring the sampler (CLAUDE.md).
    return PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="map", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )


def test_every_transcript_eff_length_comes_from_the_PAYLOADS_rna_pmf(scenario):
    """⭐ D7, verified rather than assumed.

    The pipeline's per-transcript effective lengths must be **exactly** what the payload's RNA pmf
    produces. Any other source — a stale model, the scanner's deleted histogram, a default fallback —
    gives a different array.
    """
    from rigel.calibration.fl import build_fl_models

    pr = run_pipeline(scenario.bam_path, scenario.index, config=_config())
    shipped = np.asarray(pr.estimator.effective_lengths, dtype=np.float64)

    # Rebuilt from the payload alone, by the same route production takes but from an independent scan.
    _, _, _, payload = scan_and_buffer(
        str(scenario.bam_path), scenario.index, _config().scan
    )
    fl_models = build_fl_models(payload)
    rna_fl = FragmentLengthModel.from_pmf(fl_models.rna_pmf, fl_models.max_size)
    exonic = scenario.index.t_df["length"].values.astype(np.int64)
    expected = rna_fl.compute_all_transcript_eff_lens(exonic)

    assert shipped.shape == expected.shape
    np.testing.assert_allclose(shipped, expected, rtol=0, atol=0)


def test_the_eff_lengths_MOVE_when_the_rna_pmf_moves(scenario):
    """⚠ The equality above is worth nothing unless the quantity is actually sensitive to the pmf.

    A constant array would satisfy it. Perturbing the RNA pmf must change the effective lengths, or
    D7's blast radius is imaginary and the check is decoration.
    """
    from rigel.calibration.fl import build_fl_models

    _, _, _, payload = scan_and_buffer(
        str(scenario.bam_path), scenario.index, _config().scan
    )
    fl_models = build_fl_models(payload)
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


def test_the_n_observations_GUARD_ON_THAT_PATH_CANNOT_FIRE(scenario):
    """⚠ A finding, pinned so it is not mistaken for a live safety net.

    ``pipeline`` guards the effective-length computation with ``if rna_fl.n_observations > 0``. On this
    path ``rna_fl`` always comes from :meth:`FragmentLengthModel.from_pmf`, which sets
    ``_total_weight = 1.0`` — so ``n_observations`` is **exactly 1, always**, whatever went into the
    pmf. The guard reads as "fall back to a default mean if there is no RNA data" and it can never
    take that branch.

    ⛔ That is not a bug to fix here — an empty RNA pool is EB-shrunk all the way to the unconditional
    anchor by ``build_fl_models``, which is a better answer than a hard-coded 200 bp mean, so the
    reachable behaviour is the right one. It is recorded because the guard's *appearance* of a
    fallback is what would let a future empty-pool bug hide behind it.
    """
    from rigel.calibration.fl import build_fl_models

    _, _, _, payload = scan_and_buffer(
        str(scenario.bam_path), scenario.index, _config().scan
    )
    fl_models = build_fl_models(payload)
    rna_fl = FragmentLengthModel.from_pmf(fl_models.rna_pmf, fl_models.max_size)
    assert rna_fl.n_observations == 1

    # …and it is 1 even for a pmf built from a pool with no observations at all.
    flat = np.full(fl_models.max_size + 1, 1.0 / (fl_models.max_size + 1))
    assert FragmentLengthModel.from_pmf(flat, fl_models.max_size).n_observations == 1
