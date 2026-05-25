"""Tests for hybrid-capture simulation weighting."""

import pytest

from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.annotation import GeneBuilder
from rigel.sim.genome import MutableGenome
from rigel.sim.whole_genome import (
    GDNASimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    parse_yaml_config,
)
from rigel.sim.suite import design_capture_probe_intervals, write_random_capture_probes
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _transcript(
    transcript_id: str,
    exons: list[tuple[int, int]],
    *,
    strand: Strand = Strand.POS,
    abundance: float = 100.0,
) -> Transcript:
    transcript = Transcript(
        ref="chr1",
        strand=strand,
        exons=[Interval(start, end) for start, end in exons],
        t_id=transcript_id,
        g_id=f"g_{transcript_id}",
        abundance=abundance,
    )
    transcript.compute_length()
    return transcript


def test_transcript_probe_weights_match_overlap_example(tmp_path):
    probes = tmp_path / "probes.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t200\t320\nT\t500\t620\nT\t1500\t1620\n")
    transcript = _transcript("T", [(0, 2000)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(
            probes=str(probes),
            binding_per_base=1.0,
            off_target_weight=1.0,
        ),
        [transcript],
        {"chr1": 2000},
    )

    assert sampler.fragment_weight("mrna", 0, 2000, 1000, 200) == pytest.approx(1.0)
    assert sampler.fragment_weight("mrna", 0, 2000, 1310, 200) == pytest.approx(11.0)
    assert sampler.fragment_weight("mrna", 0, 2000, 1480, 200) == pytest.approx(121.0)


def test_transcript_probe_crossing_junction_penalizes_gdna_and_nrna(tmp_path):
    probes = tmp_path / "junction_probe.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t80\t140\n")
    transcript = _transcript("T", [(100, 200), (400, 500)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(
            probes=str(probes),
            binding_per_base=1.0,
            off_target_weight=1.0,
            gdna_split_penalty=0.2,
        ),
        [transcript],
        {"chr1": 1000},
    )

    assert sampler.fragment_weight("mrna", 0, 200, 80, 60) == pytest.approx(61.0)
    assert sampler.fragment_weight("gdna", "chr1", 1000, 170, 280) == pytest.approx(13.0)
    assert sampler.fragment_weight("nrna", 0, 400, 70, 280) == pytest.approx(13.0)


def test_bed12_probe_projects_to_transcript_and_gdna(tmp_path):
    bed = tmp_path / "probes.bed"
    bed.write_text("chr1\t180\t440\tprobe1\t0\t+\t180\t440\t0\t2\t20,40\t0,220\n")
    transcript = _transcript("T", [(100, 200), (400, 500)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(
            probes=str(bed),
            probe_format="bed12",
            binding_per_base=1.0,
            off_target_weight=1.0,
            gdna_split_penalty=0.2,
        ),
        [transcript],
        {"chr1": 1000},
    )

    assert sampler.fragment_weight("mrna", 0, 200, 80, 60) == pytest.approx(61.0)
    assert sampler.fragment_weight("gdna", "chr1", 1000, 170, 280) == pytest.approx(13.0)


def test_capture_partition_increases_targeted_transcript_weight(tmp_path):
    probes = tmp_path / "probes.tsv"
    probes.write_text("T1\t400\t520\n")
    transcripts = [
        _transcript("T1", [(0, 1000)]),
        _transcript("T2", [(2000, 3000)]),
    ]
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=10.0),
        transcripts,
        {"chr1": 3000},
    )

    targeted = sampler.partition("mrna", 0, 1000, 100)
    untargeted = sampler.partition("mrna", 1, 1000, 100)

    assert targeted > untargeted * 20


def test_overlapping_probes_do_not_stack_capture_strength(tmp_path):
    probes = tmp_path / "overlap.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t100\t220\nT\t140\t260\n")
    transcript = _transcript("T", [(0, 300)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=1.0, off_target_weight=1.0),
        [transcript],
        {"chr1": 300},
    )

    assert sampler.fragment_weight("mrna", 0, 300, 140, 80) == pytest.approx(81.0)

    expected_extra = 0.0
    for start in range(300 - 80 + 1):
        frag_end = start + 80
        first = max(0, min(frag_end, 220) - max(start, 100))
        second = max(0, min(frag_end, 260) - max(start, 140))
        expected_extra += max(first, second)
    assert sampler.partition("mrna", 0, 300, 80) == pytest.approx(221 + expected_extra)


def test_duplicate_isoform_probes_do_not_stack_gdna_strength(tmp_path):
    probes = tmp_path / "duplicate.tsv"
    probes.write_text("transcript_id\tstart\tend\nT1\t20\t140\nT2\t20\t140\n")
    transcripts = [
        _transcript("T1", [(100, 300)]),
        _transcript("T2", [(100, 300)]),
    ]
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=1.0, off_target_weight=1.0),
        transcripts,
        {"chr1": 500},
    )

    assert sampler.fragment_weight("gdna", "chr1", 500, 120, 120) == pytest.approx(121.0)


def test_parse_yaml_capture_config(tmp_path):
    probes = tmp_path / "probes.tsv"
    probes.write_text("T1\t10\t130\n")
    config = tmp_path / "sim.yaml"
    config.write_text(
        f"genome: genome.fa\n"
        f"gtf: annotation.gtf\n"
        f"capture:\n"
        f"  probes: {probes}\n"
        f"  format: transcript\n"
        f"  off_target_weight: 0.5\n"
        f"  binding_per_base: 7\n"
        f"  gdna_split_penalty: 0.15\n"
        f"  min_overlap: 4\n"
    )

    cfg = parse_yaml_config(config)

    assert cfg.capture.probes == str(probes)
    assert cfg.capture.probe_format == "transcript"
    assert cfg.capture.off_target_weight == pytest.approx(0.5)
    assert cfg.capture.binding_per_base == pytest.approx(7.0)
    assert cfg.capture.gdna_split_penalty == pytest.approx(0.15)
    assert cfg.capture.min_overlap == 4


def test_whole_genome_simulator_uses_capture_partition_for_assignment(tmp_path):
    genome = MutableGenome(3000, seed=7, name="chr1")
    builder = GeneBuilder(genome)
    builder.add_gene("g1", "+", [{"t_id": "T1", "exons": [(100, 1100)]}])
    builder.add_gene("g2", "+", [{"t_id": "T2", "exons": [(1500, 2500)]}])
    transcripts = builder.get_transcripts()
    fasta = genome.write_fasta(tmp_path)

    probes = tmp_path / "probes.tsv"
    probes.write_text("T1\t400\t520\n")
    sim = WholeGenomeSimulator(
        fasta,
        transcripts,
        SimulationParams(
            sim_seed=11,
            frag_mean=100,
            frag_std=1,
            frag_min=100,
            frag_max=100,
            read_length=50,
        ),
        GDNASimConfig(),
        capture_config=CaptureConfig(probes=str(probes), binding_per_base=10.0),
    )
    try:
        mrna_counts, _ = sim._accumulate_rna_counts(0, n_mrna=500, n_nrna=0)
    finally:
        sim.close()

    t1_count = sum(mrna_counts.get(0, {}).values())
    t2_count = sum(mrna_counts.get(1, {}).values())
    assert t1_count > t2_count * 20


def test_generated_probe_tiling_is_non_overlapping_and_centered():
    exact = design_capture_probe_intervals(1200, probe_length=120, probe_density=1.0)
    assert exact == [(i * 120, (i + 1) * 120) for i in range(10)]

    with_slack = design_capture_probe_intervals(1319, probe_length=120, probe_density=1.0)
    assert len(with_slack) == 10
    assert all(end - start == 120 for start, end in with_slack)
    assert all(left[1] <= right[0] for left, right in zip(with_slack, with_slack[1:]))
    assert abs(with_slack[0][0] - (1319 - with_slack[-1][1])) <= 1

    half_density = design_capture_probe_intervals(1200, probe_length=120, probe_density=0.5)
    assert len(half_density) == 5
    assert all(left[1] <= right[0] for left, right in zip(half_density, half_density[1:]))


def test_random_mini_genome_probe_writer_targets_only_captured_transcripts(tmp_path):
    transcripts = [
        _transcript("T1", [(0, 1200)]),
        _transcript("T2", [(2000, 2500)]),
        _transcript("T3", [(3000, 3080)]),
    ]
    probes = tmp_path / "capture_probes.tsv"

    result = write_random_capture_probes(
        transcripts,
        probes,
        capture_fraction=1.0,
        probe_length=120,
        probe_density=0.5,
        seed=7,
    )

    rows = [line.split("\t") for line in probes.read_text().splitlines()[1:]]
    assert result.n_transcripts == 3
    assert result.n_eligible == 2
    assert result.n_captured == 2
    assert result.n_probes == 7
    assert {row[0] for row in rows} == {"T1", "T2"}

    by_transcript: dict[str, list[tuple[int, int]]] = {}
    for transcript_id, start, end in rows:
        by_transcript.setdefault(transcript_id, []).append((int(start), int(end)))
    for intervals in by_transcript.values():
        assert all(end - start == 120 for start, end in intervals)
        assert all(left[1] <= right[0] for left, right in zip(intervals, intervals[1:]))
