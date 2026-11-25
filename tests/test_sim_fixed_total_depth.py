"""The simulator's fixed-total depth mode (``simulation.n_total_fragments``).

``resolve_depths`` has two modes. With ``n_total_fragments`` unset (the default) the RNA depth is
fixed and gDNA is added on top at ``n_gdna = rate x n_rna``, which cannot reach the high-gDNA end of
the spectrum: a 98 % gDNA library is ``rate = 49``, so 490 M fragments at a 10 M RNA depth. With it
set, the total is held and ``rate`` decides only the split, as a sequencing run with a fixed budget
does. ``Scenario.build_oracle`` has the same two modes, and refuses a gDNA fraction its mode cannot read. Each gate
carries a perturbation that shows it can fire.
"""

from __future__ import annotations


import pytest

from rigel.sim.wgs_config import SimulationParams


from _sim_grid import resolve_depths

#: A rate ladder from zero to 98 % gDNA. Rates, not fractions: ``f_gdna = rate / (1 + rate)``.
LADDER = [0.0, 0.010101, 0.052632, 0.111111, 0.333333, 1.0, 3.0, 9.0, 49.0]
TOTAL = 10_000_000


def _params(**kw) -> SimulationParams:
    return SimulationParams(n_rna_fragments=10_000_000, **kw)


def test_the_total_is_conserved_exactly_at_every_rung():
    """``n_rna + n_gdna`` equals the requested total on every rung exactly, not to within rounding
    drift: rungs of different depths cannot be compared, which is the purpose of fixing the total.

    Perturbation: the additive path (``n_total_fragments=None``) must not conserve it, or this gate
    passes on a property the additive mode already has.
    """
    sim = _params(n_total_fragments=TOTAL)
    for rate in LADDER:
        d = resolve_depths(sim, gdna_rate=rate)
        assert d.n_rna + d.n_gdna == TOTAL, f"rate {rate}: {d.n_rna} + {d.n_gdna} != {TOTAL}"

    additive = _params(n_total_fragments=None)
    totals = {resolve_depths(additive, gdna_rate=r).total for r in LADDER}
    assert len(totals) > 1, "the additive path already fixes the total; this mode adds nothing"


def test_the_realised_gdna_fraction_tracks_the_rung():
    """The realised gDNA fraction is the requested one, checked against ``rate/(1+rate)``, which is
    what the rate means, not against the label.

    Perturbation: a ladder that produced a constant fraction would still conserve the total, so the
    gate above cannot catch it; this one requires the fractions to span the range they claim.
    """
    sim = _params(n_total_fragments=TOTAL)
    fracs = []
    for rate in LADDER:
        d = resolve_depths(sim, gdna_rate=rate)
        want = rate / (1.0 + rate)
        got = d.n_gdna / d.total
        assert abs(got - want) <= 1.0 / TOTAL, f"rate {rate}: f_gdna {got:.6f} != {want:.6f}"
        fracs.append(got)
    assert fracs[0] == 0.0, "the zero rung must be EXACTLY zero — it is the false-positive control"
    assert fracs[-1] > 0.97, "the top rung must actually reach the high-gDNA regime"
    assert fracs == sorted(fracs), "the ladder must be monotone"


def test_the_nascent_split_is_inside_the_rna_budget_not_on_top(tmp_path):
    """With nascent enabled the total still holds: nascent RNA comes out of the RNA share, not on
    top of it. The nascent entities are rows of the one RNA multinomial, so the gate is on the engine:
    every fragment the RNA draw assigns, mature or nascent, is one of the ``n_rna`` the budget gave it.

    Perturbation: the same draw with the entity's nascent abundance zeroed must put strictly more
    fragments on the mature row, or the gate is vacuous.
    """
    from rigel.sim.annotation import GeneBuilder
    from rigel.sim.genome import MutableGenome
    from rigel.sim.wgs_config import GDNASimConfig
    from rigel.sim.wgs_engine import WholeGenomeSimulator
    from rigel.sim.transcriptome import assign_nascent_to_entities
    from rigel.transcript import Transcript
    from rigel.types import Interval, Strand

    d = resolve_depths(_params(n_total_fragments=TOTAL), gdna_rate=1.0)
    assert d.n_rna + d.n_gdna == TOTAL

    genome = MutableGenome(4000, seed=3, name="chr1")
    builder = GeneBuilder(genome)
    builder.add_gene(
        "g1", "+", [{"t_id": "T1", "exons": [(500, 1000), (2000, 2500)], "abundance": 100.0}]
    )
    fasta = genome.write_fasta(tmp_path)
    mature = builder.get_transcripts()[0]
    mature.t_index = 0
    mature.span_t_index = 1
    entity = Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(500, 2500)],
        t_id="N1",
        g_id="N1",
        t_index=1,
        is_single_exon=True,
        is_synthetic=True,
        abundance=0.0,
    )
    entity.length = entity.compute_length()

    def draw(nascent_per_mature: float) -> tuple[int, int]:
        import numpy as np

        rows = [mature, entity]
        assign_nascent_to_entities(rows, np.array([100.0 * nascent_per_mature, 0.0]))
        sim = WholeGenomeSimulator(
            fasta,
            rows,
            SimulationParams(
                sim_seed=5, frag_mean=100, frag_std=1, frag_min=100, frag_max=100, read_length=50
            ),
            GDNASimConfig(),
            genomic_refs=["chr1"],
        )
        try:
            m, n = sim._accumulate_rna_counts(2000)
        finally:
            sim.close()
        return sum(sum(v.values()) for v in m.values()), sum(sum(v.values()) for v in n.values())

    m_on, n_on = draw(0.25)
    assert m_on + n_on == 2000, "every RNA fragment is one of the budget's n_rna"
    assert n_on > 0 and n_on < 2000
    m_off, n_off = draw(0.0)
    assert n_off == 0 and m_off == 2000 and m_off > m_on


def test_omitting_the_total_gives_the_additive_depths_exactly():
    """The mode is opt-in: a config that omits ``n_total_fragments`` gets the additive depths
    exactly (``n_rna`` is the configured RNA depth and ``n_gdna = round(rate x n_rna)``).

    Perturbation: switching the mode on must move the numbers, or the identity is vacuous.
    """
    additive = _params(n_total_fragments=None)
    for rate in LADDER:
        d = resolve_depths(additive, gdna_rate=rate)
        assert d.n_rna == 10_000_000
        assert d.n_gdna == round(rate * 10_000_000)

    on = resolve_depths(_params(n_total_fragments=TOTAL), gdna_rate=1.0)
    off = resolve_depths(additive, gdna_rate=1.0)
    assert (on.n_rna, on.n_gdna) != (off.n_rna, off.n_gdna)


def test_a_total_too_small_for_the_rung_raises_rather_than_rounding_to_zero():
    """At rate 49 a 100-fragment total leaves 2 RNA fragments; a smaller one rounds the RNA side to
    zero, and the condition would silently stop being an RNA-seq library, so it raises.

    Perturbation: a total that is large enough must not raise.
    """
    with pytest.raises(ValueError, match="RNA"):
        resolve_depths(_params(n_total_fragments=10), gdna_rate=49.0)
    ok = resolve_depths(_params(n_total_fragments=TOTAL), gdna_rate=49.0)
    assert ok.n_rna == 200_000


def test_a_scenario_refuses_a_gdna_fraction_it_would_ignore(tmp_path):
    """``Scenario.build_oracle`` reads ``gdna_fraction`` only on top of a fixed RNA depth (``n_rna_fragments``) and
    with a ``GDNAConfig`` to draw gDNA from. A fixed total (``n_fragments``) is split by the gDNA config's abundance
    instead, so a fraction passed there, or without a gDNA config, would simulate a library with no gDNA from a call
    that asked for some. Both are refused before anything is written.

    Perturbation: the same fraction on a fixed RNA depth with a ``GDNAConfig`` is read, and adds its gDNA fragments.
    """
    from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

    sc = Scenario("frac", genome_length=4000, seed=3, work_dir=tmp_path / "sim")
    sc.add_gene("g1", "+", [{"t_id": "T1", "exons": [(500, 1000), (2000, 2500)], "abundance": 100}])
    sim = ReadSimConfig(
        frag_mean=150, frag_std=20, frag_min=80, frag_max=250, read_length=50, seed=3
    )
    with pytest.raises(ValueError, match="n_rna_fragments"):
        sc.build_oracle(
            n_fragments=100, sim_config=sim, gdna_fraction=0.2, gdna_config=GDNAConfig()
        )
    with pytest.raises(ValueError, match="GDNAConfig"):
        sc.build_oracle(n_rna_fragments=100, sim_config=sim, gdna_fraction=0.2)
    assert not (tmp_path / "sim" / "index").exists(), "a refused call must write nothing"

    result = sc.build_oracle(
        n_rna_fragments=100, sim_config=sim, gdna_fraction=0.2, gdna_config=GDNAConfig()
    )
    assert result.n_simulated == 120
