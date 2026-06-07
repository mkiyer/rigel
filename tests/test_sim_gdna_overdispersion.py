"""Simulator gDNA strand overdispersion (docs/em_strand/03 §3).

Confirms the `gdna_strand_overdispersion` knob (a) converts to the internal Beta
concentration correctly and (b) produces per-region sense rates with the intended
overdispersion — so the suite can generate overdispersed-gDNA conditions for the BB fit.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.sim.reads import GDNAConfig


def test_overdispersion_converts_to_beta_concentration():
    # od = 1/(kappa+1)  ⇒  kappa = (1-od)/od
    assert GDNAConfig(gdna_strand_overdispersion=0.2).strand_kappa == pytest.approx(4.0)
    assert GDNAConfig(gdna_strand_overdispersion=0.05).strand_kappa == pytest.approx(19.0)
    # od == 0 ⇒ Binomial (no region partition)
    assert GDNAConfig(gdna_strand_overdispersion=0.0).strand_kappa is None
    # explicit strand_kappa still works when the clear knob is unset
    assert GDNAConfig(strand_kappa=9.0).strand_kappa == 9.0


@pytest.mark.parametrize("bad", [-0.1, 1.0, 1.5])
def test_overdispersion_out_of_range_rejected(bad):
    with pytest.raises(ValueError):
        GDNAConfig(gdna_strand_overdispersion=bad)


@pytest.mark.parametrize("od", [0.05, 0.2])
def test_whole_genome_gdna_strand_regions_overdispersion(od):
    """WholeGenomeSimulator (the suite simulator) builds per-ref strand regions at target od."""
    from rigel.sim.whole_genome import WholeGenomeSimulator
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
    from rigel.sim.whole_genome import GDNASimConfig

    assert GDNASimConfig().strand_overdispersion == 0.0
