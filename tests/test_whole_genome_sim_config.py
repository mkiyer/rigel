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


def test_random_nrna_fraction_only_spikes_expressed_multi_exon_transcripts():
    transcripts = [
        _transcript("multi_a", 100.0, [(0, 100), (200, 300)]),
        _transcript("multi_b", 200.0, [(400, 500), (700, 800)]),
        _transcript("single", 300.0, [(900, 1200)]),
        _transcript("off", 0.0, [(1300, 1400), (1500, 1600)]),
    ]

    realized_ratio = apply_random_nrna_fraction(
        transcripts,
        (0.5, 0.5),
        eligible_fraction=1.0,
        seed=17,
    )

    assert transcripts[0].nrna_abundance == pytest.approx(50.0)
    assert transcripts[1].nrna_abundance == pytest.approx(100.0)
    assert transcripts[2].nrna_abundance == 0.0
    assert transcripts[3].nrna_abundance == 0.0
    assert realized_ratio == pytest.approx(150.0 / 600.0)


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
        "genome: /tmp/genome.fa\n"
        "gtf: /tmp/genes.gtf\n"
        "nrna:\n"
        "  mode: random_fraction\n"
    )

    with pytest.raises(ValueError, match="ratio_ranges"):
        parse_yaml_config(config_path)