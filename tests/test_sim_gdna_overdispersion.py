"""The simulator's gDNA strand overdispersion knob.

`gdna_strand_overdispersion` must convert to the internal Beta concentration correctly and produce
per-region sense rates with the overdispersion it asked for, so the suite can generate the
overdispersed-gDNA conditions the Beta-Binomial strand fit is developed against. A knob that
converted wrongly would leave every such condition quietly at the default dispersion.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.sim.reads import GDNAConfig


@pytest.mark.parametrize("bad", [-0.1, 1.0, 1.5])
def test_overdispersion_out_of_range_rejected(bad):
    with pytest.raises(ValueError):
        GDNAConfig(gdna_strand_overdispersion=bad)


def test_scenario_gdna_config_overdispersion_reaches_the_engine():
    """Scenario's GDNAConfig → engine GDNASimConfig: the overdispersion knob flows through, and an
    unset knob or no gDNA config reaches the engine as 0."""
    from rigel.sim.scenario import _to_gdna_sim

    assert _to_gdna_sim(
        GDNAConfig(gdna_strand_overdispersion=0.2)
    ).strand_overdispersion == pytest.approx(0.2)
    assert _to_gdna_sim(GDNAConfig()).strand_overdispersion == 0.0
    assert _to_gdna_sim(None).strand_overdispersion == 0.0


@pytest.mark.parametrize("od", [0.05, 0.2])
def test_whole_genome_gdna_strand_regions_overdispersion(od):
    """WholeGenomeSimulator (the suite simulator) builds per-ref strand regions at target od."""
    from rigel.sim.wgs_engine import WholeGenomeSimulator
    from rigel.transcript import Transcript
    from rigel.types import Interval, Strand

    sim = WholeGenomeSimulator.__new__(WholeGenomeSimulator)
    sim._rng = np.random.default_rng(0)
    sim.transcripts = [
        Transcript(
            ref="chr1",
            strand=Strand.POS,
            exons=[Interval(s, s + 500), Interval(s + 1500, s + 2000)],
        )
        for s in range(1000, 1000 + 60 * 5000, 5000)  # 60 genes → >30 exon-derived regions
    ]
    sim._gdna_refs = ["chr1"]
    sim._gdna_ref_lengths = [400_000]
    sim._gdna_strand_regions = {}
    sim._init_gdna_strand_regions(od)

    boundaries, p_plus = sim._gdna_strand_regions["chr1"]
    assert len(p_plus) > 30
    # Beta(a, a) with a = ½(1−od)/od has mean ½ and variance ¼·od (intra-class correlation).
    assert float(np.mean(p_plus)) == pytest.approx(0.5, abs=0.05)
    assert float(np.var(p_plus)) == pytest.approx(0.25 * od, rel=0.30)


def test_whole_genome_gdna_config_default_no_overdispersion():
    """Default GDNASimConfig has overdispersion 0 (no strand regions ⇒ uniform 50/50)."""
    from rigel.sim.wgs_config import GDNASimConfig

    assert GDNASimConfig().strand_overdispersion == 0.0
