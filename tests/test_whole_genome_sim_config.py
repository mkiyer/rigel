"""Tests for whole-genome simulator configuration and abundance helpers."""

import pytest

from rigel.sim.whole_genome import apply_random_nrna_fraction, parse_yaml_config
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _transcript(
    transcript_id: str,
    abundance: float,
    exons: list[tuple[int, int]],
) -> Transcript:
    return Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(start, end) for start, end in exons],
        t_id=transcript_id,
        g_id=f"gene_{transcript_id}",
        abundance=abundance,
    )


def _entity(transcript_id: str, start: int, end: int, t_index: int) -> Transcript:
    """An nRNA entity as the index makes it: single-exon over the clustered span, synthetic."""
    t = Transcript(
        ref="chr1", strand=Strand.POS, exons=[Interval(start, end)], t_id=transcript_id,
        g_id=transcript_id, t_index=t_index, is_nrna=True, is_synthetic=True, abundance=0.0,
    )
    t.length = t.compute_length()
    return t


def test_random_nrna_fraction_pools_expressed_multi_exon_nascent_onto_entities():
    """⭐ Nascent molecules live on the index's nRNA ENTITIES (owner, 2026-08-19): each expressed
    multi-exon transcript contributes ``abundance × ratio`` molecules to the entity its
    ``nrna_t_index`` names, single-exon and silent transcripts contribute nothing, and the
    contributors themselves end with ``nrna_abundance = 0``."""
    multi_a = _transcript("multi_a", 100.0, [(0, 100), (200, 300)])
    multi_b = _transcript("multi_b", 200.0, [(400, 500), (700, 800)])
    single = _transcript("single", 300.0, [(900, 1200)])
    off = _transcript("off", 0.0, [(1300, 1400), (1500, 1600)])
    shared = _entity("N_shared", 0, 800, 4)  # multi_a and multi_b cluster onto one span
    lone = _entity("N_off", 1300, 1600, 5)
    for i, t in enumerate((multi_a, multi_b, single, off)):
        t.t_index = i
    multi_a.nrna_t_index = multi_b.nrna_t_index = 4
    off.nrna_t_index = 5
    transcripts = [multi_a, multi_b, single, off, shared, lone]

    realized_ratio = apply_random_nrna_fraction(
        transcripts,
        (0.5, 0.5),
        eligible_fraction=1.0,
        seed=17,
    )

    assert shared.nrna_abundance == pytest.approx(50.0 + 100.0)
    assert lone.nrna_abundance == 0.0
    assert all(t.nrna_abundance == 0.0 for t in (multi_a, multi_b, single, off))
    assert realized_ratio == pytest.approx(150.0 / 600.0)


def test_a_multi_exon_contributor_with_no_entity_is_a_defect_not_a_skip():
    """⛔ Every multi-exon transcript of a rigel index links to an entity; a transcript list that does
    not is not the index's, and pooling onto a missing entity raises rather than dropping molecules."""
    from rigel.sim.whole_genome import apply_nrna_ratio

    orphan = _transcript("orphan", 100.0, [(0, 100), (200, 300)])
    orphan.t_index = 0
    with pytest.raises(ValueError, match="nascent entity"):
        apply_nrna_ratio([orphan], 0.25)


def test_fragment_share_solves_the_molecular_ratio_from_the_annotation():
    """⭐⭐ **A PANEL STATES THE NASCENT FRAGMENT SHARE; THE MOLECULAR RATIO IS DERIVED** (owner,
    2026-08-19). A nascent ENTITY spans a whole gene while a mature transcript is spliced, so the two
    are NOT interchangeable: with a 10x length ratio a molecular ratio of 0.25 is nowhere near a 25 %
    fragment share, and the panel's meaning would move silently.

    PERTURBATIONS: (a) the solved share must be REACHED — feeding the solved ratio back through the
    weights reproduces the target; (b) the naive reading (ratio = share) must NOT reach it, or the
    solve is doing nothing; (c) share 0 leaves no nascent at all.
    """
    from rigel.sim.whole_genome import (
        apply_nrna_fragment_share,
        apply_nrna_ratio,
        expected_rna_weights,
    )
    from rigel.sim.wgs_config import SimulationParams

    sim = SimulationParams(frag_mean=200, frag_std=40, frag_min=100, frag_max=400)
    # one mature transcript of 2,000 spliced bp; its entity spans 20,000 bp — the 10x that matters
    mature = _transcript("T", 100.0, [(0, 1000), (19000, 20000)])
    mature.t_index = 0
    mature.nrna_t_index = 1
    entity = _entity("N", 0, 20000, 1)
    rows = [mature, entity]

    ratio = apply_nrna_fragment_share(rows, 0.20, sim)
    w_m, w_n = expected_rna_weights(rows, sim)
    assert w_n / (w_m + w_n) == pytest.approx(0.20, abs=1e-9), "(a) the target share must be reached"
    assert ratio == pytest.approx(0.0227, abs=5e-4), "the derived ratio is ~10x below the share"

    apply_nrna_ratio(rows, 0.20)  # (b) the naive reading
    w_m2, w_n2 = expected_rna_weights(rows, sim)
    assert w_n2 / (w_m2 + w_n2) > 0.65, "(b) a molecular ratio of 0.20 gives a MUCH larger share"

    assert apply_nrna_fragment_share(rows, 0.0, sim) == 0.0  # (c)
    assert entity.nrna_abundance == 0.0


def test_fragment_share_refuses_an_unreachable_target():
    """⛔ With no expressed multi-exon transcript there is no nascent opportunity, and a nonzero share
    is unreachable — that must raise, not silently produce a nascent-free library."""
    from rigel.sim.whole_genome import apply_nrna_fragment_share
    from rigel.sim.wgs_config import SimulationParams

    single = _transcript("S", 100.0, [(0, 2000)])
    single.t_index = 0
    with pytest.raises(ValueError, match="unreachable"):
        apply_nrna_fragment_share([single], 0.2, SimulationParams())


def test_the_fl_pmf_is_the_one_the_engine_draws_from():
    """⚠ The share solve integrates over the fragment-length pmf. If that pmf were not the engine's,
    the solved ratio would be right about a distribution nothing samples.

    PERTURBATION: the analytic pmf must match a large SAMPLE from `truncated_normal_frag_lengths`,
    and must NOT match one drawn with a different sd.
    """
    import numpy as np

    from rigel.sim.sampling import truncated_normal_frag_lengths
    from rigel.sim.wgs_config import SimulationParams
    from rigel.sim.whole_genome import fl_pmf

    sim = SimulationParams(frag_mean=206, frag_std=98, frag_min=50, frag_max=500)
    widths, p = fl_pmf(sim)
    rng = np.random.default_rng(0)
    draws = truncated_normal_frag_lengths(rng, 400_000, sim.frag_mean, sim.frag_std,
                                          sim.frag_min, sim.frag_max)
    emp = np.bincount(draws - sim.frag_min, minlength=len(widths))[: len(widths)] / len(draws)
    assert np.abs(emp - p).max() < 2e-3, "the analytic pmf must be the sampler's"
    wrong = fl_pmf(SimulationParams(frag_mean=206, frag_std=40, frag_min=50, frag_max=500))[1]
    assert np.abs(emp - wrong).max() > 5e-3, "a different sd must be distinguishable"


def test_parse_random_fraction_nrna_config(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\n"
        "gtf: /tmp/genes.gtf\n"
        "nrna:\n"
        "  mode: random_fraction\n"
        "  ratio_ranges:\n"
        "    - [0.0, 0.0]\n"
        "    - [0.05, 0.25]\n"
        "  ratio_labels: [zero, low]\n"
        "  eligible_fraction: 0.5\n"
        "  seed: 123\n"
    )

    config = parse_yaml_config(config_path)

    assert config.nrna.mode == "random_fraction"
    assert config.nrna.ratio_ranges == [(0.0, 0.0), (0.05, 0.25)]
    assert config.nrna.ratio_labels == ["zero", "low"]
    assert config.nrna.eligible_fraction == 0.5
    assert config.nrna.seed == 123


def test_random_fraction_requires_ratio_ranges(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\nnrna:\n  mode: random_fraction\n"
    )

    with pytest.raises(ValueError, match="ratio_ranges"):
        parse_yaml_config(config_path)
