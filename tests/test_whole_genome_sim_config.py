"""Tests for whole-genome simulator configuration and abundance helpers."""

import numpy as np
import pytest

from rigel.sim.whole_genome import apply_sparse_nrna, parse_yaml_config
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


def test_sparse_nascent_pools_onto_entities_and_the_SPAN_is_the_unit():
    """⭐⭐ Nascent molecules live on the index's nRNA ENTITIES, and under the SPARSE model the entity
    is also the unit of the on/off draw (owner, 2026-08-22) — one gene span is transcribed or it is
    not, however many isoforms share it. A span whose only contributor is SILENT gets nothing, single-
    exon transcripts contribute nothing, and the contributors themselves end at ``nrna_abundance = 0``.
    ⛔ Drawing per CONTRIBUTOR instead would make a 5-isoform gene 41 % likely to carry nascent at
    ``on_fraction = 0.1``, so the intron slots calibration reads would be four times less sparse than
    configured."""
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

    realized_ratio = apply_sparse_nrna(
        transcripts, (150.0, 150.0), on_fraction=1.0, seed=17
    )

    # the two contributors SHARE one span, so the span is ONE draw — not one per isoform
    assert shared.nrna_abundance == pytest.approx(150.0)
    assert lone.nrna_abundance == 0.0, "a span whose only contributor is silent is not transcribed"
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


def test_parse_sparse_nrna_config(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\n"
        "gtf: /tmp/genes.gtf\n"
        "nrna:\n"
        "  mode: sparse\n"
        "  abundance_ranges:\n"
        "    - [1.0, 10.0]\n"
        "    - [10.0, 1000.0]\n"
        "  ratio_labels: [low, mixed]\n"
        "  on_fraction: 0.5\n"
        "  seed: 123\n"
    )

    config = parse_yaml_config(config_path)

    assert config.nrna.mode == "sparse"
    assert config.nrna.abundance_ranges == [(1.0, 10.0), (10.0, 1000.0)]
    assert config.nrna.ratio_labels == ["low", "mixed"]
    assert config.nrna.on_fraction == 0.5
    assert config.nrna.seed == 123


def test_sparse_requires_abundance_ranges(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\nnrna:\n  mode: sparse\n"
    )

    with pytest.raises(ValueError, match="abundance_ranges"):
        parse_yaml_config(config_path)


def test_a_log_uniform_range_may_not_touch_zero(tmp_path):
    """⛔ A log-uniform draw has no zero end, so ``[0, x]`` is not "no nascent" — it is undefined.
    Absence is expressed by ``on_fraction``, which is the model's own switch."""
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\nnrna:\n  mode: sparse\n"
        "  abundance_ranges:\n    - [0.0, 10.0]\n  ratio_labels: [bad]\n"
    )
    with pytest.raises(ValueError, match="0 < min"):
        parse_yaml_config(config_path)




def _sparse_population(n: int, *, seed: int = 3) -> tuple[list[Transcript], list[Transcript]]:
    """``n`` expressed multi-exon transcripts, one entity each, so an entity's ``nrna_abundance`` reads
    back exactly one draw. Mature abundances span the ladder's own three decades (log-uniform 1 →
    10,000), which is what lets the independence gate below see a correlation if there is one."""
    rng = np.random.default_rng(seed)
    contributors, entities = [], []
    for i in range(n):
        base = i * 10_000
        t = _transcript(f"t{i}", float(10.0 ** rng.uniform(0, 4)), [(base, base + 100), (base + 200, base + 300)])
        t.t_index = i
        t.nrna_t_index = n + i
        contributors.append(t)
        entities.append(_entity(f"N{i}", base, base + 300, n + i))
    return contributors, entities


def test_the_on_fraction_leaves_most_gene_spans_at_EXACTLY_zero():
    """⭐⭐ **SPARSITY IS THE MODEL AND NOTHING GATED IT** — the retired mode's only behavioural test ran
    at full eligibility, where the Bernoulli draw cannot be observed at all.

    Nascent RNA is absent from most gene spans and present in a minority (owner, 2026-08-22), so
    ``on_fraction = 0.1`` must leave ~90 % of eligible spans at nascent EXACTLY zero — not small, zero
    — while the rest carry a real level. The band is the binomial's own five sigma, derived rather than
    tuned.
    """
    n, frac = 400, 0.1
    contributors, entities = _sparse_population(n)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=frac, seed=11)

    on = [e for e in entities if e.nrna_abundance > 0.0]
    off = [e for e in entities if e.nrna_abundance == 0.0]
    assert len(on) + len(off) == n, "a span is either on or exactly zero"
    assert on and off, "zero-inflation is not observable — the Bernoulli draw did nothing"
    sd = (n * frac * (1 - frac)) ** 0.5
    lo, hi = n * frac - 5 * sd, n * frac + 5 * sd
    assert lo <= len(on) <= hi, (
        f"{len(on)}/{n} spans on at on_fraction={frac}: outside the binomial's own "
        f"[{lo:.0f}, {hi:.0f}] five-sigma band, so the draw is not Bernoulli({frac})"
    )


def test_on_fraction_zero_transcribes_NOTHING_and_one_transcribes_EVERY_ELIGIBLE_span():
    """The two closed ends, exactly — a fraction is a fraction only if its endpoints are honoured."""
    contributors, entities = _sparse_population(40)
    assert apply_sparse_nrna(
        contributors + entities, (1.0, 1000.0), on_fraction=0.0, seed=5
    ) == 0.0
    assert all(e.nrna_abundance == 0.0 for e in entities)

    contributors, entities = _sparse_population(40)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=1.0, seed=5)
    assert all(e.nrna_abundance > 0.0 for e in entities)


def test_the_level_is_LOG_uniform_over_the_range_not_linear():
    """⭐ Where nascent is present its level spans decades — very low in some spans, high in others
    (owner, 2026-08-22). ⛔ A LINEAR draw on (1, 1000) would put **90 % of its mass in the top decade**
    and could not express "very low"; a log-uniform draw spreads the ORDERS OF MAGNITUDE evenly, so the
    decade occupancies are near-equal. That is the property gated here, and it is what distinguishes
    the two draws.
    """
    lo, hi = 1.0, 1000.0
    contributors, entities = _sparse_population(600)
    apply_sparse_nrna(contributors + entities, (lo, hi), on_fraction=1.0, seed=7)
    levels = np.array([e.nrna_abundance for e in entities])

    assert levels.min() >= lo and levels.max() <= hi, "levels escaped the range"
    # three decades, so a log-uniform draw puts ~1/3 of the mass in each; a linear draw puts ~0.9 in
    # the last one. The band is generous — it only has to separate 1/3 from 9/10.
    per_decade = [np.mean((levels >= 10.0**k) & (levels < 10.0 ** (k + 1))) for k in range(3)]
    assert max(per_decade) < 0.5, (
        f"decade occupancies {[round(x, 3) for x in per_decade]} — the level is concentrated in one "
        f"decade, which is a linear draw, not a log-uniform one"
    )
    assert min(per_decade) > 0.15, f"a decade is nearly empty: {[round(x, 3) for x in per_decade]}"


def test_the_nascent_LEVEL_IS_INDEPENDENT_OF_THE_MATURE_LEVEL():
    """⛔⛔ **THE DECISION THIS MODE EXISTS FOR** (owner, 2026-08-22): the nascent level must be drawn
    INDEPENDENTLY of the mature level, so ``nascent > mature`` is a real case the tool has to survive.
    The retired ratio modes set ``nascent = mature x ratio``, which makes the two perfectly rank-
    correlated and puts a ceiling under the interesting case.

    Two things are asserted: the rank correlation between the mature and nascent levels is ~0 (a
    ratio model reads 1.0), and spans where nascent EXCEEDS mature actually occur.
    """
    contributors, entities = _sparse_population(500, seed=9)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=1.0, seed=23)

    mature = np.array([c.abundance for c in contributors])
    nascent = np.array([e.nrna_abundance for e in entities])
    rank_m = np.argsort(np.argsort(mature)).astype(float)
    rank_n = np.argsort(np.argsort(nascent)).astype(float)
    rho = float(np.corrcoef(rank_m, rank_n)[0, 1])
    # 500 independent pairs: sd of Spearman's rho under the null is ~1/sqrt(n-1) = 0.045, so five
    # sigma is 0.22 — derived from the sample size, not tuned
    assert abs(rho) < 5.0 / (len(mature) - 1) ** 0.5, (
        f"rank correlation between mature and nascent levels is {rho:.3f} — the nascent level is "
        f"tracking the mature one, which is the ratio model this mode replaced"
    )
    assert (nascent > mature).mean() > 0.1, (
        f"only {(nascent > mature).mean():.1%} of spans carry more nascent than mature — the "
        f"robustness case the owner asked for is not being exercised"
    )


def test_the_on_off_DRAW_IS_PER_SPAN_NOT_PER_ISOFORM():
    """⛔⛔ **THE UNIT OF SPARSITY, GATED WHERE IT IS ACTUALLY OBSERVABLE.** With ``k`` isoforms sharing
    one gene span, a per-CONTRIBUTOR draw gives that span ``1 - (1 - p)^k`` chance of carrying nascent
    — at ``p = 0.1`` and ``k = 5`` that is **41 %**, four times the configured sparsity — while a
    per-SPAN draw gives exactly ``p``. The intron slots calibration reads belong to the span, so the
    per-span reading is the one that makes the configured number mean what it says.

    ⚠ This gate exists because a perturbation to per-contributor drawing passed every other test in
    this file: the shared-span test runs at ``on_fraction = 1.0``, where the two are indistinguishable
    (`TRAPS: perturb-every-gate`).
    """
    n_spans, k, p = 300, 5, 0.1
    rng = np.random.default_rng(4)
    rows: list[Transcript] = []
    entities: list[Transcript] = []
    for sp in range(n_spans):
        base = sp * 100_000
        entity = _entity(f"N{sp}", base, base + 50_000, 10_000 + sp)
        entities.append(entity)
        for j in range(k):
            off = base + j * 1_000
            t = _transcript(f"t{sp}_{j}", float(10.0 ** rng.uniform(0, 4)),
                            [(off, off + 100), (off + 200, off + 300)])
            t.t_index = sp * k + j
            t.nrna_t_index = 10_000 + sp
            rows.append(t)
    apply_sparse_nrna(rows + entities, (1.0, 1000.0), on_fraction=p, seed=31)

    on_frac = float(np.mean([e.nrna_abundance > 0.0 for e in entities]))
    per_isoform = 1.0 - (1.0 - p) ** k  # 0.41 — what a per-contributor draw would give
    sd = (p * (1 - p) / n_spans) ** 0.5
    assert abs(on_frac - p) < 5 * sd, (
        f"{on_frac:.3f} of spans on at on_fraction={p} with {k} isoforms each: the per-span value is "
        f"{p}, a per-isoform draw would give {per_isoform:.3f} — the draw is not per span"
    )


def test_the_SEED_actually_drives_the_draw():
    """⛔ **THE SEED WAS PARSED, ASSERTED AND UNUSED-ABLE.** Every other gate here is statistical, so
    replacing `default_rng(seed)` with `default_rng(42)` passed all of them — and the orchestrator adds
    `seed + nrna_index` precisely so a multi-range sweep gets DIFFERENT on-sets, which would have
    collapsed to one pattern silently (`TRAPS: perturb-every-gate`).

    Two seeds must disagree, and one seed must reproduce itself exactly.
    """
    def draw(seed):
        contributors, entities = _sparse_population(200)
        apply_sparse_nrna(contributors + entities, (1.0, 100.0), on_fraction=0.5, seed=seed)
        return np.array([e.nrna_abundance for e in entities])

    a, b, a_again = draw(1), draw(2), draw(1)
    assert np.array_equal(a, a_again), "the same seed did not reproduce its own draw"
    assert not np.array_equal(a, b), "two different seeds produced identical draws — the seed is inert"
    assert (a > 0).sum() and (b > 0).sum(), "a draw with nothing on cannot distinguish anything"


def test_a_multi_exon_contributor_with_no_ENTITY_raises_in_the_SPARSE_path_too():
    """⛔ The entity-defect rule was gated only through `apply_nrna_ratio`, so replacing the sparse
    path's own raise with `continue` passed every test. A transcript list without entities is not a
    rigel index's, and dropping its molecules silently is the failure the rule exists to stop."""
    orphan = _transcript("orphan", 100.0, [(0, 100), (200, 300)])
    orphan.t_index = 0  # nrna_t_index stays unset/-1: no entity
    with pytest.raises(ValueError, match="nascent entity"):
        apply_sparse_nrna([orphan], (1.0, 100.0), on_fraction=1.0, seed=1)


def test_the_FUNCTION_validates_its_own_range_and_fraction_not_just_the_PARSER():
    """⛔⛔ **THE PARSER IS NOT THE ONLY DOOR.** `suite.py` builds `NRNAConfig` from argv and validates
    neither the range nor the fraction, so `apply_sparse_nrna`'s own guards are the only thing standing
    between `--nrna-abundance-ranges '0,100'` and `log(0) = -inf` inside the draw. Deleting either
    guard passed all 13 gates before this one existed."""
    rows = sum(_sparse_population(4), [])
    for bad, why in [((0.0, 100.0), "zero lower end has no log"),
                     ((-1.0, 10.0), "negative lower end"),
                     ((100.0, 1.0), "inverted endpoints")]:
        with pytest.raises(ValueError, match="0 < lo <= hi|abundance_range"):
            apply_sparse_nrna(list(rows), bad, on_fraction=0.5, seed=1)
    for bad in (1.5, -0.1):
        with pytest.raises(ValueError, match="on_fraction"):
            apply_sparse_nrna(list(rows), (1.0, 10.0), on_fraction=bad, seed=1)
