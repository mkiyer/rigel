"""Hybrid-capture simulation — the weight a probe panel puts on a transcript and on gDNA, the
partition that turns those weights into per-fragment opportunity, and the length selection that
opportunity has to produce. The first block gates the weighting and the panel plumbing: probe-to-
transcript projection for split and BED12 probes, overlapping and duplicate probes that must not
stack, the YAML capture config and its sweep, external panels, the paired condition seed, the probe
group key, and the random probe writer. The second gates `CaptureSampler.partition_array` against a
brute-force enumeration over every start position, together with the array entry point and the bound
on what it caches. The third gates that capture moves the fragment length marginal and not merely the
conditional.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.sim.annotation import GeneBuilder
from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.capture.design import design_capture_probe_intervals, write_random_capture_probes
from rigel.sim.capture.sampler import ProbeInterval
from rigel.sim.genome import MutableGenome, random_dna_array
from rigel.sim.manifest import condition_dir_name
from rigel.sim.orchestrator import capture_paired_condition_seed
from rigel.sim.whole_genome import (
    GDNASimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    parse_yaml_config,
)
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _transcript(
    transcript_id: str,
    exons: list[tuple[int, int]],
    *,
    strand: Strand = Strand.POS,
    abundance: float = 100.0,
    gene_id: str | None = None,
) -> Transcript:
    transcript = Transcript(
        ref="chr1",
        strand=strand,
        exons=[Interval(start, end) for start, end in exons],
        t_id=transcript_id,
        g_id=gene_id or f"g_{transcript_id}",
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


def test_a_split_probe_binds_every_template_through_ONE_contiguous_part(tmp_path):
    """A probe across a splice junction is one contiguous 60-mer only in a transcript that holds the
    junction; everywhere else — gDNA, the NASCENT ENTITY (a single-exon transcript over the span, the
    index's nRNA row), an isoform that lacks the junction — it is two separate parts, and a fragment
    hybridises through one contiguous stretch, so it binds the better part and never the sum. There is
    no penalty beyond that geometry: the part it does hold binds as it would for any molecule
    (owner, 2026-09-19). Here the probe is 20 bases in the first exon and 40 in the second."""
    probes = tmp_path / "sj_probe.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t80\t140\n")
    transcript = _transcript("T", [(100, 200), (400, 500)])
    entity = _transcript("T_nascent", [(100, 500)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=1.0, off_target_weight=1.0),
        [transcript, entity],
        {"chr1": 1000},
    )

    # the spliced transcript holds the whole probe contiguously
    assert sampler.fragment_weight("mrna", 0, 200, 80, 60) == pytest.approx(61.0)
    # gDNA over both parts, across the intron: the better part (40), not their sum
    assert sampler.fragment_weight("gdna", "chr1", 1000, 170, 280) == pytest.approx(41.0)
    # the entity: transcript coordinate 70 of its 400-bp span is genomic 170, the same fragment
    assert sampler.fragment_weight("mrna", 1, 400, 70, 280) == pytest.approx(41.0)
    # there is no per-transcript nascent space: the entity IS a transcript row
    with pytest.raises(ValueError):
        sampler.fragment_weight("nrna", 0, 400, 70, 280)


def test_gdna_and_cdna_bind_the_SAME_half_of_a_split_probe_alike(tmp_path):
    """Hybridisation does not know what a molecule is: a gDNA fragment and a cDNA fragment of an
    isoform that lacks the junction hold the same 20 bases of the probe's first part, and weigh the
    same. Only a molecule holding the junction can bind the probe whole."""
    probes = tmp_path / "sj_probe.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t80\t140\n")
    with_junction = _transcript("T", [(100, 200), (400, 500)])
    without_junction = _transcript("T2", [(100, 200), (600, 700)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=1.0, off_target_weight=1.0),
        [with_junction, without_junction],
        {"chr1": 1000},
    )

    # T2 coordinates 60..120 are genomic 160..200 then 600..620: they hold genomic 180..200
    cdna_half = sampler.fragment_weight("mrna", 1, 200, 60, 60)
    gdna_half = sampler.fragment_weight("gdna", "chr1", 1000, 150, 60)
    assert cdna_half == pytest.approx(21.0)
    assert gdna_half == pytest.approx(cdna_half)


def test_a_probe_reaches_every_transcript_its_genomic_blocks_overlap(tmp_path):
    """A probe maps to compatible transcripts by GENOMIC overlap: any isoform, any gene, either
    strand — and to gDNA. A sibling isoform sharing the probed exon, a single-exon
    transcript of another gene under the probe on the OPPOSITE strand, and a nascent entity spanning it
    all carry the probe; a transcript whose exons miss it does not."""
    probes = tmp_path / "probes.tsv"
    probes.write_text("transcript_id\tstart\tend\nT1\t20\t80\n")  # genomic [120, 180) of T1
    t1 = _transcript("T1", [(100, 200), (400, 500)])
    sibling = _transcript("T2", [(100, 200), (600, 700)])  # shares exon 1, not probed itself
    antisense = _transcript(
        "A", [(150, 300)], strand=Strand.NEG
    )  # another gene, other strand, under the probe
    entity = _transcript("N", [(100, 700)])  # a nascent entity spanning the lot
    far = _transcript("F", [(800, 900)])  # no overlap
    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(probes), binding_per_base=1.0, off_target_weight=1.0),
        [t1, sibling, antisense, entity, far],
        {"chr1": 1000},
    )
    # a 60-bp fragment lying exactly under the probe, in each transcript's own coordinates
    assert sampler.fragment_weight("mrna", 0, 200, 20, 60) == pytest.approx(61.0)  # T1
    assert sampler.fragment_weight("mrna", 1, 200, 20, 60) == pytest.approx(61.0)  # sibling
    # antisense transcript A = genomic [150, 300) on −, so its coordinate x is genomic 300 − x and the
    # probe's overlap with it is genomic [150, 180) = A's [120, 150): a 60-bp fragment at A's 90
    # covers genomic [150, 210) and binds the 30 overlapping bases
    assert sampler.fragment_weight("mrna", 2, 150, 90, 60) == pytest.approx(31.0)
    assert sampler.fragment_weight("mrna", 3, 600, 20, 60) == pytest.approx(61.0)  # entity
    assert sampler.fragment_weight("mrna", 4, 100, 20, 60) == pytest.approx(1.0)  # far: off target
    assert sampler.fragment_weight("gdna", "chr1", 1000, 120, 60) == pytest.approx(61.0)


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
        ),
        [transcript],
        {"chr1": 1000},
    )

    assert sampler.fragment_weight("mrna", 0, 200, 80, 60) == pytest.approx(61.0)
    # gDNA holds the two blocks apart across the intron and binds the better one (40)
    assert sampler.fragment_weight("gdna", "chr1", 1000, 170, 280) == pytest.approx(41.0)


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
        f"  min_overlap: 4\n"
    )

    cfg = parse_yaml_config(config)

    assert cfg.capture.probes == str(probes)
    assert cfg.capture.probe_format == "transcript"
    assert cfg.capture.off_target_weight == pytest.approx(0.5)
    assert cfg.capture.binding_per_base == pytest.approx(7.0)
    assert cfg.capture.min_overlap == 4


def test_a_capture_key_the_loader_does_not_know_is_REFUSED(tmp_path):
    """Ignored, an unknown key would read as applied — a config still carrying the retired
    ``gdna_split_penalty`` would simulate the new physics while its YAML claimed the old."""
    probes = tmp_path / "probes.tsv"
    probes.write_text("T1\t10\t130\n")
    config = tmp_path / "sim.yaml"
    config.write_text(
        f"genome: genome.fa\n"
        f"gtf: annotation.gtf\n"
        f"capture:\n"
        f"  probes: {probes}\n"
        f"  gdna_split_penalty: 0.2\n"
    )
    with pytest.raises(ValueError, match="gdna_split_penalty"):
        parse_yaml_config(config)


def test_parse_yaml_capture_config_sweep(tmp_path):
    probes = tmp_path / "probes.tsv"
    probes.write_text("T1\t10\t130\n")
    config = tmp_path / "sim.yaml"
    config.write_text(
        f"genome: genome.fa\n"
        f"gtf: annotation.gtf\n"
        f"capture:\n"
        f"  off_target_weight: 0.5\n"
        f"  configs:\n"
        f"    - label: off\n"
        f"      enabled: false\n"
        f"    - label: on\n"
        f"      probes: {probes}\n"
        f"      format: transcript\n"
        f"      binding_per_base: 7\n"
    )

    cfg = parse_yaml_config(config)

    assert [scenario.label for scenario in cfg.capture_configs] == ["off", "on"]
    assert cfg.capture_configs[0].config.probes is None
    assert cfg.capture_configs[1].config.probes == str(probes)
    assert cfg.capture_configs[1].config.probe_format == "transcript"
    assert cfg.capture_configs[1].config.off_target_weight == pytest.approx(0.5)
    assert cfg.capture_configs[1].config.binding_per_base == pytest.approx(7.0)


def test_capture_sweep_uses_paired_condition_seed():
    seed = capture_paired_condition_seed(42, "none", 0.99)

    assert seed == capture_paired_condition_seed(42, "none", 0.99)
    assert seed != capture_paired_condition_seed(42, "high", 0.99)
    assert condition_dir_name("none", 0.99, "none", "off") != condition_dir_name(
        "none",
        0.99,
        "none",
        "on",
    )


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
        genomic_refs=[genome.name],
        capture_config=CaptureConfig(probes=str(probes), binding_per_base=10.0),
    )
    try:
        mrna_counts, _ = sim._accumulate_rna_counts(500)
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
    assert result.n_eligible_genes == 2
    assert result.n_captured_genes == 2
    assert result.n_probes == 7
    assert {row[0] for row in rows} == {"T1", "T2"}

    by_transcript: dict[str, list[tuple[int, int]]] = {}
    for transcript_id, start, end in rows:
        by_transcript.setdefault(transcript_id, []).append((int(start), int(end)))
    for intervals in by_transcript.values():
        assert all(end - start == 120 for start, end in intervals)
        assert all(left[1] <= right[0] for left, right in zip(intervals, intervals[1:]))


def test_random_probe_writer_outputs_bed12_for_sj_spanning_probe(tmp_path):
    transcripts = [_transcript("T", [(100, 200), (400, 500)])]
    probes = tmp_path / "capture_probes.tsv"
    bed = tmp_path / "capture_probes.bed"

    result = write_random_capture_probes(
        transcripts,
        probes,
        capture_fraction=1.0,
        probe_length=120,
        probe_density=1.0,
        seed=7,
        bed_path=bed,
    )

    assert result.bed_path == bed
    fields = bed.read_text().strip().split("\t")
    assert fields[:6] == ["chr1", "140", "460", "T:probe_1", "0", "+"]
    assert fields[9] == "2"
    assert fields[10] == "60,60"
    assert fields[11] == "0,260"


def test_random_probe_writer_masks_shared_isoform_sequence_by_abundance(tmp_path):
    transcripts = [
        _transcript("LOW", [(100, 200), (300, 400)], abundance=10.0, gene_id="G"),
        _transcript("HIGH", [(100, 200), (300, 400)], abundance=1000.0, gene_id="G"),
    ]
    probes = tmp_path / "capture_probes.tsv"
    bed = tmp_path / "capture_probes.bed"

    result = write_random_capture_probes(
        transcripts,
        probes,
        capture_fraction=1.0,
        probe_length=100,
        probe_density=1.0,
        seed=7,
        bed_path=bed,
    )

    rows = [line.split("\t") for line in probes.read_text().splitlines()[1:]]
    assert result.n_captured == 2
    assert result.n_eligible_genes == 1
    assert result.n_captured_genes == 1
    assert result.n_probes == 2
    assert {row[0] for row in rows} == {"HIGH"}

    sampler = CaptureSampler.from_config(
        CaptureConfig(probes=str(bed), probe_format="bed12", binding_per_base=1.0),
        transcripts,
        {"chr1": 1000},
    )

    assert sampler.fragment_weight("mrna", 0, 200, 0, 100) > 1.0
    assert sampler.fragment_weight("mrna", 1, 200, 0, 100) > 1.0


def test_random_probe_writer_selects_capture_pool_by_gene(tmp_path):
    transcripts = [
        _transcript("GA.1", [(0, 240)], gene_id="GA"),
        _transcript("GA.2", [(500, 740)], gene_id="GA"),
        _transcript("GB.1", [(1000, 1240)], gene_id="GB"),
        _transcript("GB.2", [(1500, 1740)], gene_id="GB"),
    ]
    probes = tmp_path / "capture_probes.tsv"

    result = write_random_capture_probes(
        transcripts,
        probes,
        capture_fraction=0.5,
        probe_length=120,
        probe_density=1.0,
        seed=7,
    )

    rows = [line.split("\t") for line in probes.read_text().splitlines()[1:]]
    row_genes = {row[0].split(".")[0] for row in rows}
    expected_by_gene = {
        "GA": {"GA.1", "GA.2"},
        "GB": {"GB.1", "GB.2"},
    }

    assert result.n_eligible == 4
    assert result.n_captured == 2
    assert result.n_eligible_genes == 2
    assert result.n_captured_genes == 1
    assert len(row_genes) == 1
    assert {row[0] for row in rows} == expected_by_gene[next(iter(row_genes))]


# ── ``CaptureSampler.partition_array`` against brute-force enumeration over every start ──────────
#
# The capture partition path dominates a capture-on simulation and is therefore optimised, and the
# standing rule is that every divisor is unit-tested against brute-force enumeration. So the
# optimisation is gated on equality with a reference that literally loops over start positions and
# takes the max. That oracle is deliberately naive and slow — an explicit Python loop, no numpy, no
# helper shared with the implementation, because a validator that calls the builder's own helper
# validates nothing. The matrix below covers what the real panel deliberately does not, so the
# optimisation cannot silently encode a property of one panel: probes of unequal length (the panel is
# uniformly 120 bp); probes whose start-ranges do and do not interact, both regimes being live across
# the fragment lengths swept; ``min_overlap > 1``, which zeroes part of each trapezoid; a fragment
# longer than, equal to and flush against the template length; a transcript with no probes, and a
# pool where no transcript has probes.


def brute_force_partition(
    intervals: list[ProbeInterval],
    seq_len: int,
    frag_len: int,
    *,
    off_target_weight: float,
    binding_per_base: float,
    min_overlap: int,
) -> float:
    """The definition, spelled out: for every start, the fragment's best overlap with ONE contiguous
    probe part, summed over starts. A fragment hybridises through one contiguous stretch, so separate
    parts — of different probes, or of one probe split across an intron — never add."""
    total = 0.0
    for start in range(max(0, seq_len - frag_len + 1)):
        total += off_target_weight
        best = 0
        for interval in intervals:
            overlap = min(start + frag_len, interval.end) - max(start, interval.start)
            if overlap < min_overlap or overlap <= 0:
                continue
            best = max(best, overlap)
        total += binding_per_base * best
    return total


def make_sampler(intervals_by_key: dict[int, list[ProbeInterval]], **overrides) -> CaptureSampler:
    """A sampler whose mRNA intervals are injected directly, bypassing probe-file parsing."""
    config = CaptureConfig(
        probes="unused",
        probe_format="transcript",
        off_target_weight=overrides.get("off_target_weight", 1.0),
        binding_per_base=overrides.get("binding_per_base", 10.0),
        min_overlap=overrides.get("min_overlap", 1),
    )
    # Built through the real `__init__`, never `__new__`: only the probe intervals are then injected,
    # which is the one thing this fixture exists to control. Bypassing the constructor and hand-setting
    # a few attributes instead means any state `__init__` later gains is simply ABSENT here, so the
    # whole file fails with `AttributeError` while the sampler works perfectly in production
    # (`TRAPS: a-gate-that-reconstructs` — a fixture that rebuilds its subject tests the rebuild).
    sampler = CaptureSampler(config, transcripts=(), ref_lengths={}, enabled=True)
    sampler._mrna_intervals = intervals_by_key
    return sampler


def iv(start: int, end: int) -> ProbeInterval:
    return ProbeInterval(start=start, end=end)


# ── the layouts. Each is a named shape the optimisation could get wrong differently ────────────────
LAYOUTS: dict[str, list[ProbeInterval]] = {
    "single probe, mid-template": [iv(300, 420)],
    "single probe, flush at start": [iv(0, 120)],
    "single probe, flush at end": [iv(880, 1000)],
    "two probes, far apart": [iv(100, 220), iv(700, 820)],
    "two probes, adjacent (gap 2)": [iv(100, 220), iv(222, 342)],
    "two probes, UNEQUAL length": [iv(100, 160), iv(300, 700)],
    "two parts, OVERLAPPING": [iv(100, 220), iv(160, 280)],
    "three probes, dense": [iv(100, 220), iv(240, 360), iv(380, 500)],
    "a split probe's two separate parts": [iv(100, 160), iv(400, 460)],
    "no probes": [],
}

#: The layout the gDNA space is made of: every probe that spans an intron lands in genomic coordinates as
#: separate parts.
SPLIT_LAYOUT = "a split probe's two separate parts"
BATCHABLE_LAYOUTS = list(LAYOUTS)

FRAGMENT_LENGTHS = [1, 2, 50, 119, 120, 121, 200, 500, 999, 1000, 1001, 2000]
SEQ_LEN = 1000


class TestAgainstBruteForce:
    """Every case goes through `partition_array` — the batched path — not the scalar one."""

    @pytest.mark.parametrize("layout_name", BATCHABLE_LAYOUTS)
    @pytest.mark.parametrize("frag_len", FRAGMENT_LENGTHS)
    def test_partition_array_matches_enumeration(self, layout_name: str, frag_len: int) -> None:
        intervals = LAYOUTS[layout_name]
        sampler = make_sampler({0: intervals} if intervals else {})
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=1,
        )
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        actual = float(sampler.partition_array("mrna", [0], lengths, frag_len)[0])
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9), (
            f"{layout_name} at w={frag_len}: got {actual}, enumeration says {expected}"
        )

    @pytest.mark.parametrize("layout_name", list(LAYOUTS))
    @pytest.mark.parametrize("frag_len", FRAGMENT_LENGTHS)
    def test_the_scalar_path_also_matches_enumeration(
        self, layout_name: str, frag_len: int
    ) -> None:
        """`partition` (one key) is gated against the same enumeration as the batched call."""
        intervals = LAYOUTS[layout_name]
        sampler = make_sampler({0: intervals} if intervals else {})
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=1,
        )
        assert sampler.partition("mrna", 0, SEQ_LEN, frag_len) == pytest.approx(
            expected, rel=1e-12, abs=1e-9
        )

    @pytest.mark.parametrize("min_overlap", [1, 5, 60, 121])
    @pytest.mark.parametrize("frag_len", [50, 200, 500])
    def test_min_overlap_gate_matches_enumeration(self, min_overlap: int, frag_len: int) -> None:
        intervals = LAYOUTS["three probes, dense"]
        sampler = make_sampler({0: intervals}, min_overlap=min_overlap)
        expected = brute_force_partition(
            intervals,
            SEQ_LEN,
            frag_len,
            off_target_weight=1.0,
            binding_per_base=10.0,
            min_overlap=min_overlap,
        )
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        actual = float(sampler.partition_array("mrna", [0], lengths, frag_len)[0])
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9)

    def test_a_multi_key_pool_matches_enumeration_key_by_key(self) -> None:
        """The batched path shares ONE buffer across keys. A key writing outside its own slice, or a
        rank colliding between keys, only shows up with several keys of DIFFERENT lengths at once."""
        names = list(BATCHABLE_LAYOUTS)
        intervals_by_key = {i: LAYOUTS[name] for i, name in enumerate(names)}
        sampler = make_sampler(intervals_by_key)
        lengths = np.array([SEQ_LEN + 137 * i for i in range(len(names))], dtype=np.int64)

        for frag_len in (50, 200, 500, 1200):
            actual = sampler.partition_array("mrna", list(range(len(names))), lengths, frag_len)
            for key, name in enumerate(names):
                expected = brute_force_partition(
                    LAYOUTS[name],
                    int(lengths[key]),
                    frag_len,
                    off_target_weight=1.0,
                    binding_per_base=10.0,
                    min_overlap=1,
                )
                assert actual[key] == pytest.approx(expected, rel=1e-12, abs=1e-9), (
                    f"key {key} ({name}) at w={frag_len}"
                )

    def test_a_split_probes_parts_are_never_SUMMED(self) -> None:
        """A fragment over both parts of a split probe binds the better part, not their sum. This is
        the case the gDNA space is made of, so the gate reads a fragment long enough to hold both."""
        intervals = LAYOUTS[SPLIT_LAYOUT]
        sampler = make_sampler({0: intervals})
        # a 400-base fragment from 80 holds 60 of each part: its weight is one part's, 1 + 10 * 60
        assert sampler.fragment_weight("mrna", 0, SEQ_LEN, 80, 400) == pytest.approx(601.0)
        sampler = make_sampler({0: intervals})
        lengths = np.array([SEQ_LEN], dtype=np.int64)
        for frag_len in (50, 200, 500):
            expected = brute_force_partition(
                intervals,
                SEQ_LEN,
                frag_len,
                off_target_weight=1.0,
                binding_per_base=10.0,
                min_overlap=1,
            )
            assert float(
                sampler.partition_array("mrna", [0], lengths, frag_len)[0]
            ) == pytest.approx(expected, rel=1e-12, abs=1e-9)


class TestPartitionArray:
    """`partition_array` must agree with `partition` element for element — it is the batched form."""

    @pytest.mark.parametrize("frag_len", [50, 200, 500, 1200])
    def test_array_agrees_with_the_scalar_over_a_mixed_pool(self, frag_len: int) -> None:
        names = list(BATCHABLE_LAYOUTS)
        intervals_by_key = {i: LAYOUTS[name] for i, name in enumerate(names)}
        intervals_by_key[len(names)] = []  # a transcript with no probes, in the middle of the pool
        lengths = np.array([SEQ_LEN + 37 * i for i in range(len(intervals_by_key))], dtype=np.int64)
        sampler = make_sampler(intervals_by_key)

        actual = sampler.partition_array("mrna", range(len(lengths)), lengths, frag_len)
        expected = np.array(
            [
                sampler.partition("mrna", key, int(lengths[key]), frag_len)
                for key in range(len(lengths))
            ]
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-9)

    def test_a_key_absent_from_the_call_does_not_leak_into_another(self) -> None:
        """`partition_array` may be called with a SUBSET of keys — the live-abundance filter does
        exactly that. A key left out must not shift the buffer slices of the keys left in."""
        intervals_by_key = {
            0: LAYOUTS["three probes, dense"],
            1: LAYOUTS["two probes, far apart"],
            2: LAYOUTS["single probe, mid-template"],
        }
        sampler = make_sampler(intervals_by_key)
        lengths_all = np.array([SEQ_LEN, SEQ_LEN + 500, SEQ_LEN + 900], dtype=np.int64)
        full = sampler.partition_array("mrna", [0, 1, 2], lengths_all, 200)
        subset = sampler.partition_array("mrna", [0, 2], lengths_all[[0, 2]], 200)
        np.testing.assert_allclose(subset, full[[0, 2]], rtol=1e-12, atol=1e-9)

    def test_a_pool_where_nothing_is_on_panel_is_pure_off_target(self) -> None:
        sampler = make_sampler({})
        lengths = np.array([1000, 500, 250], dtype=np.int64)
        actual = sampler.partition_array("mrna", range(3), lengths, 200)
        np.testing.assert_array_equal(actual, lengths - 200 + 1)

    def test_templates_shorter_than_the_fragment_contribute_nothing(self) -> None:
        sampler = make_sampler({0: LAYOUTS["three probes, dense"]})
        lengths = np.array([150], dtype=np.int64)
        assert sampler.partition_array("mrna", range(1), lengths, 200)[0] == 0.0


class TestNoUnboundedCache:
    """A per-call cache here costs tens of gigabytes and never hits.

    `partition_array` visits each (key, fragment length) pair exactly once, and `sample_starts` is
    driven by a counts dict keyed by exactly that pair, so anything keyed per call is written and
    never read back. This pins that no such state exists.
    """

    def test_no_per_call_state_accumulates_across_fragment_lengths(self) -> None:
        """One dict is allowed to hold one entry per WIDTH, and the rule is stricter than "no
        growth" about the thing that actually costs the memory. A cache keyed per call grows without
        ever being read; `_partition_memo` is read back by every later condition sharing the capture
        panel, and it is bounded because a different template population CLEARS it. So every OTHER
        dict must not grow, and the memo must be bounded by the number of distinct widths rather
        than by the number of calls.
        """
        sampler = make_sampler({i: LAYOUTS["three probes, dense"] for i in range(40)})
        lengths = np.full(40, SEQ_LEN, dtype=np.int64)
        sizes_before = {k: len(v) for k, v in vars(sampler).items() if isinstance(v, dict)}
        widths = list(range(50, 350))
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), lengths, frag_len)
        for name, size in sizes_before.items():
            if name == "_partition_memo":
                continue
            grown = len(getattr(sampler, name)) - size
            assert grown <= 1, (
                f"{name} grew by {grown} entries across 300 fragment lengths. Anything that scales "
                f"with fragment length here is the 38 GB coming back."
            )
        assert len(sampler._partition_memo) == len(widths), (
            "the memo holds exactly one entry per width"
        )

        # REPEATING the same sweep must add NOTHING — that is what makes it a cache rather than a leak
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), lengths, frag_len)
        assert len(sampler._partition_memo) == len(widths), "a repeat sweep must be pure cache hits"

        # and a DIFFERENT template population must not accumulate beside the first
        other = np.full(40, SEQ_LEN + 1, dtype=np.int64)
        for frag_len in widths:
            sampler.partition_array("mrna", range(40), other, frag_len)
        assert len(sampler._partition_memo) == len(widths), (
            "a second population must CLEAR the memo, not accumulate beside it — that is the 38 GB"
        )

    def test_the_memo_returns_what_recomputation_would(self) -> None:
        """A cache that returns a stale or aliased vector is worse than no cache. Every width is
        compared against a freshly-built sampler's answer, and the returned array must be a COPY (a
        caller mutating its result must not poison the memo)."""
        layout = {i: LAYOUTS["three probes, dense"] for i in range(4)}
        cached_sampler = make_sampler(layout)
        lengths = np.full(4, SEQ_LEN, dtype=np.int64)
        for frag_len in (60, 120, 121, 200):
            first = cached_sampler.partition_array("mrna", range(4), lengths, frag_len)
            first[:] = -12345.0  # a caller mutating its result
            again = cached_sampler.partition_array("mrna", range(4), lengths, frag_len)
            fresh = make_sampler(layout).partition_array("mrna", range(4), lengths, frag_len)
            np.testing.assert_allclose(again, fresh, err_msg=f"memo differs at width {frag_len}")

    def test_sample_starts_does_not_accumulate_state_either(self) -> None:
        sampler = make_sampler({0: LAYOUTS["three probes, dense"]})
        rng = np.random.default_rng(0)
        sizes_before = {k: len(v) for k, v in vars(sampler).items() if isinstance(v, dict)}
        for frag_len in range(50, 350):
            sampler.sample_starts("mrna", 0, SEQ_LEN, frag_len, 3, rng)
        for name, size in sizes_before.items():
            grown = len(getattr(sampler, name)) - size
            assert grown <= 1, f"{name} grew by {grown} entries across 300 fragment lengths"


# ── Hybrid capture selects for LENGTH, and it has to do so in the marginal ──────────────────────
#
# Hybrid capture works by hybridisation of probes to sequence. A short fragment presents less
# sequence to hybridise, binds less well, and is captured less efficiently, so capture selects for
# longer fragments — gDNA and RNA alike. Two consequences are not optional: the post-capture
# distributions are the true baseline, because the pre-capture ones describe a library that was never
# sequenced; and capture narrows the gDNA-to-RNA length gap whenever gDNA is the shorter component,
# because the short tail it removes is disproportionately gDNA. The failure these gates guard is
# drawing the length marginal first, capture-blind, then computing a capture-aware per-template
# opportunity ``total_eff(w)`` and discarding it by normalising within each length: capture can then
# move only WHERE a fragment lands, never WHETHER its length survives, and the post-capture truth
# summaries come out identical to the pre-capture ones. What must hold instead is
# ``f_post(w) ∝ f_pre(w) · total_eff(w)``. The selection gates are directional and carry no
# threshold, because choosing a magnitude would be inventing the capture efficiency curve, which the
# no-magic-numbers rule forbids — the sign is physics, the magnitude is whatever ``binding_per_base``
# and the probe length imply. The narrowing gate is conditional on gDNA being the shorter, broader
# component, which is what this fixture is; never assume that ordering in general, because for other
# library types it reverses.


GENOMIC_REFS = ["chrA", "chrB"]

#: Fragment-length parameters taken from a real cfRNA library's pools — the count-weighted mean and
#: sd of each, not chosen round numbers.
RNA_LENGTH = dict(frag_mean=206.0, frag_std=98.0, frag_min=50, frag_max=800)
GDNA_LENGTH = dict(frag_mean=157.0, frag_std=125.0, frag_min=40, frag_max=1000)


def _write_fasta(tmp_path: Path, sequences: dict[str, str]) -> Path:
    path = tmp_path / "genome.fa"
    with open(path, "w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def _panel_transcript(t_id, ref, exons, abundance) -> Transcript:
    transcript = Transcript(
        ref=ref,
        strand=Strand.POS,
        exons=[Interval(start, end) for start, end in exons],
        t_id=t_id,
        g_id=f"g_{t_id}",
        g_name=f"g_{t_id}",
        abundance=abundance,
    )
    transcript.compute_length()
    return transcript


@pytest.fixture(scope="module")
def panel(tmp_path_factory):
    """Two 200 kb genomic references carrying 24 transcripts, half of them on a probe panel.

    Sized so the geometry is the pilot's rather than a toy's: a probe run is a small neighbourhood
    inside a large template, which is the regime in which capture's length preference is a real
    reweighting instead of a boundary effect.
    """
    tmp_path = tmp_path_factory.mktemp("capture_length")
    rng = np.random.default_rng(11)
    lengths = {"chrA": 200_000, "chrB": 200_000}
    sequences = {name: "".join(random_dna_array(n, rng)) for name, n in lengths.items()}
    fasta = _write_fasta(tmp_path, sequences)

    transcripts = []
    for ref_index, ref in enumerate(GENOMIC_REFS):
        for i in range(12):
            start = 5_000 + i * 15_000
            # Two exons, so the transcript is 4,000 nt of spliced sequence over a 9,000 bp span.
            exons = [(start, start + 2_000), (start + 7_000, start + 9_000)]
            transcripts.append(
                _panel_transcript(f"T{ref_index}_{i}", ref, exons, abundance=100.0 * (i + 1))
            )

    # Probes on every other transcript, one 120 bp probe well inside the spliced sequence.
    probes = tmp_path / "probes.tsv"
    rows = ["transcript_id\tstart\tend"]
    for index, transcript in enumerate(transcripts):
        if index % 2 == 0:
            rows.append(f"{transcript.t_id}\t1500\t1620")
    probes.write_text("\n".join(rows) + "\n")
    return fasta, transcripts, str(probes)


def _simulator(fasta, transcripts, probes, *, seed=23):
    return WholeGenomeSimulator(
        fasta,
        transcripts,
        SimulationParams(sim_seed=seed, read_length=100, **RNA_LENGTH),
        GDNASimConfig(**GDNA_LENGTH),
        seed=seed,
        genomic_refs=GENOMIC_REFS,
        capture_config=(
            CaptureConfig(probes=probes, probe_format="transcript", binding_per_base=10.0)
            if probes is not None
            else None
        ),
    )


def _mean(length_counts: dict[int, int]) -> float:
    n = sum(length_counts.values())
    return sum(length * count for length, count in length_counts.items()) / n


def _gdna_lengths(simulator, n: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for (_ref_index, length), count in simulator._accumulate_gdna_counts(n).items():
        counts[length] = counts.get(length, 0) + count
    return counts


def _mrna_lengths(simulator, n: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    mrna, _nrna = simulator._accumulate_rna_counts(n)
    for per_length in mrna.values():
        for length, count in per_length.items():
            counts[length] = counts.get(length, 0) + count
    return counts


@pytest.fixture(scope="module")
def arms(panel):
    """``{"off"|"on": {"gdna": mean, "rna": mean}}`` — capture is the ONLY thing varied.

    Each pool gets its own freshly seeded simulator, and that is what makes this an A/B. A
    simulator's first use of its RNG is the pre-capture length draw, so a fresh instance per pool
    means both arms draw the SAME pre-capture lengths and every difference downstream is
    attributable to capture. Measuring both pools from one instance does not: the gDNA pool consumes
    a different amount of RNG in each arm, so the RNA draw that follows it differs between the arms
    by accident. PERTURBATION: removing the reweighting from the RNA pool entirely, the RNA gate
    still reads a rise and says nothing.
    """
    fasta, transcripts, probes = panel
    n = 20_000
    result = {}
    for label, probe_file in (("off", None), ("on", probes)):
        gdna_sim = _simulator(fasta, transcripts, probe_file)
        try:
            gdna = _gdna_lengths(gdna_sim, n)
        finally:
            gdna_sim.close()
        rna_sim = _simulator(fasta, transcripts, probe_file)
        try:
            rna = _mrna_lengths(rna_sim, n)
        finally:
            rna_sim.close()
        result[label] = {
            "gdna": _mean(gdna),
            "rna": _mean(rna),
            "n_gdna": sum(gdna.values()),
            "n_rna": sum(rna.values()),
        }
    return result


class TestCaptureMovesTheLengthMarginal:
    def test_the_arms_are_the_same_size(self, arms):
        """Every requested fragment is still allocated — the reweighting is not a filter."""
        for label, values in arms.items():
            assert values["n_gdna"] == 20_000, label
            assert values["n_rna"] == 20_000, label

    def test_g_s3_gdna_mean_length_rises_under_capture(self, arms):
        """Strictly greater. Byte-identical arms are the falsification."""
        off, on = arms["off"]["gdna"], arms["on"]["gdna"]
        assert on > off, f"gDNA mean length {off:.2f} (off) -> {on:.2f} (on)"

    def test_capture_also_lengthens_the_rna_pool(self, arms):
        """The principle applies to gDNA and RNA alike — a probe hybridises to either."""
        off, on = arms["off"]["rna"], arms["on"]["rna"]
        assert on > off, f"RNA mean length {off:.2f} (off) -> {on:.2f} (on)"

    def test_g_s5_the_gdna_to_rna_length_gap_narrows_under_capture(self, arms):
        """`mu_g - mu_r` is the only thing that would let fragment length identify composition: at
        equal component means it carries exactly zero information at any depth."""
        gap_off = abs(arms["off"]["gdna"] - arms["off"]["rna"])
        gap_on = abs(arms["on"]["gdna"] - arms["on"]["rna"])
        assert gap_on < gap_off, f"|mu_g - mu_r| {gap_off:.2f} (off) -> {gap_on:.2f} (on)"

    def test_the_short_tail_is_what_capture_removes(self, panel):
        """The mechanism, not just its sign: the deficit is concentrated at the short end."""
        fasta, transcripts, probes = panel
        n = 20_000
        histograms = {}
        for label, probe_file in (("off", None), ("on", probes)):
            simulator = _simulator(fasta, transcripts, probe_file)
            try:
                histograms[label] = _gdna_lengths(simulator, n)
            finally:
                simulator.close()

        def share_below(counts, limit):
            total = sum(counts.values())
            return sum(c for length, c in counts.items() if length < limit) / total

        assert share_below(histograms["on"], 100) < share_below(histograms["off"], 100)


class TestTheConditionalWasAlreadyRight:
    """A regression guard rather than a falsification, and the distinction is the finding.

    "On-target gDNA is longer than off-target" holds even when capture is not length-selective at
    all, because the CONDITIONAL is right either way: at a fixed width, the probability of landing on
    a probe rises with the width. What a non-selective engine is missing is the MARGINAL — it
    normalises ``total_eff(w)`` away, so the length distribution never moves. The marginal gate above
    is therefore the falsification and this one guards the piece that already worked
    (`TRAPS: a-gate-that-already-passed`, `TRAPS: right-conditional-wrong-marginal`).

    An on-target population defined by the START's territory cannot be satisfied by any capture model
    of this form, and reads INVERTED: conditioned on being captured, a fragment whose start is in the
    intron is one that was long enough to REACH the probe, so an intronic start carries weight
    ~ w^2/2 while an exonic start carries ~ p^2/2, flat in w
    (`TRAPS: on-target-by-start-is-geometry`). Ask what a pool's selection rule correlates with, not
    only what the selected fragments are.
    """

    def test_the_on_probe_share_of_capture_mass_rises_with_fragment_length(self, tmp_path):
        """Checked against the closed form, not against a recorded number.

        For a probe ``[b, e)`` sitting far enough inside a template that no start is truncated,
        ``sum_s overlap(s, w) = (e - b) * w`` exactly — every probe base is covered by exactly ``w``
        of the admissible starts. So the capture-weighted opportunity is
        ``off_target * (L - w + 1) + binding_per_base * p * w``, whose probe-attributable share is
        strictly increasing in ``w``. That is the whole of "capture selects for length".
        """
        transcript_length, probe_start, probe_length = 20_000, 8_000, 120
        transcript = _panel_transcript("T", "chrA", [(0, transcript_length)], 100.0)
        probes = tmp_path / "probes.tsv"
        probes.write_text(
            f"transcript_id\tstart\tend\nT\t{probe_start}\t{probe_start + probe_length}\n"
        )
        config = CaptureConfig(
            probes=str(probes),
            probe_format="transcript",
            off_target_weight=1.0,
            binding_per_base=10.0,
        )
        sampler = CaptureSampler.from_config(config, [transcript], {"chrA": transcript_length})

        shares = []
        for width in (60, 120, 240, 480, 800):
            total = sampler.partition("mrna", 0, transcript_length, width)
            baseline = config.off_target_weight * (transcript_length - width + 1)
            expected = baseline + config.binding_per_base * probe_length * width
            assert total == pytest.approx(expected), f"width {width}"
            shares.append((total - baseline) / total)

        assert shares == sorted(shares)
        assert shares[0] < shares[-1]
