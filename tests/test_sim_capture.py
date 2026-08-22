"""Tests for hybrid-capture simulation weighting."""

from types import SimpleNamespace

import pytest

from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.capture.design import design_capture_probe_intervals, write_random_capture_probes
from rigel.sim.annotation import GeneBuilder
from rigel.sim.genome import MutableGenome
from rigel.sim.whole_genome import (
    GDNASimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    parse_yaml_config,
)
from rigel.sim.manifest import condition_dir_name
from rigel.sim.suite import (
    SuiteCaptureSpec,
    _load_suite_config,
    _suite_capture_specs,
    capture_paired_condition_seed,
    capture_probe_group_key,
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


def _suite_capture_args(**overrides) -> SimpleNamespace:
    values = {
        "capture_configs": None,
        "capture_fraction": 0.0,
        "capture_probes": None,
        "capture_probe_format": "auto",
        "probe_length": 120,
        "probe_density": 1.0,
        "capture_off_target_weight": 1.0,
        "capture_binding_per_base": 10.0,
        "capture_gdna_split_penalty": 0.2,
        "capture_min_overlap": 1,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


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


def test_transcript_probe_crossing_sj_penalizes_gdna_and_the_nascent_entity(tmp_path):
    """A probe across a splice junction binds the spliced transcript whole (scale 1), and binds gDNA
    and the NASCENT ENTITY — a single-exon transcript over the span, the index's nRNA row — as two
    separated genomic pieces at ``gdna_split_penalty``. ⭐ The entity is an ordinary transcript in the
    list (owner, 2026-08-19): its capture comes from every probe whose genomic blocks it overlaps, so
    its weight equals the gDNA weight at the same position, which is the physics."""
    probes = tmp_path / "sj_probe.tsv"
    probes.write_text("transcript_id\tstart\tend\nT\t80\t140\n")
    transcript = _transcript("T", [(100, 200), (400, 500)])
    entity = _transcript("T_nascent", [(100, 500)])
    sampler = CaptureSampler.from_config(
        CaptureConfig(
            probes=str(probes),
            binding_per_base=1.0,
            off_target_weight=1.0,
            gdna_split_penalty=0.2,
        ),
        [transcript, entity],
        {"chr1": 1000},
    )

    assert sampler.fragment_weight("mrna", 0, 200, 80, 60) == pytest.approx(61.0)
    assert sampler.fragment_weight("gdna", "chr1", 1000, 170, 280) == pytest.approx(13.0)
    # the entity: transcript coordinate 70 of its 400-bp span is genomic 170, the same fragment
    assert sampler.fragment_weight("mrna", 1, 400, 70, 280) == pytest.approx(13.0)
    # ⛔ there is no per-transcript nascent space any more
    with pytest.raises(ValueError):
        sampler.fragment_weight("nrna", 0, 400, 70, 280)


def test_a_probe_reaches_every_transcript_its_genomic_blocks_overlap(tmp_path):
    """⭐ Owner, 2026-08-19: a probe maps to compatible transcripts by GENOMIC overlap, any isoform, any
    gene, either strand — and to gDNA. A sibling isoform sharing the probed exon, a single-exon
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


def test_suite_capture_config_accepts_top_level_external_bed_panel(tmp_path):
    probes = tmp_path / "panel.bed"
    probes.write_text("chr1\t10\t130\tprobe1\t0\t+\t10\t130\t0\t1\t120\t0\n")
    config = tmp_path / "suite.yaml"
    config.write_text(f"capture:\n  probes: {probes}\n  format: bed12\n  binding_per_base: 7\n")

    values = _load_suite_config(config)
    specs, include_capture_in_names = _suite_capture_specs(_suite_capture_args(**values))

    assert include_capture_in_names is False
    assert len(specs) == 1
    assert specs[0].label == "on"
    assert specs[0].enabled
    assert specs[0].uses_provided_probes
    assert not specs[0].generates_probes
    assert specs[0].probes == str(probes)
    assert specs[0].probe_format == "bed12"
    assert specs[0].binding_per_base == pytest.approx(7.0)


def test_suite_capture_config_can_mix_off_generated_and_external_panels(tmp_path):
    probes = tmp_path / "panel.bed"
    probes.write_text("chr1\t10\t130\tprobe1\t0\t+\t10\t130\t0\t1\t120\t0\n")
    args = _suite_capture_args(
        capture_fraction=0.25,
        capture_configs=[
            {"label": "off", "enabled": False},
            {"label": "generated", "fraction": 0.5},
            {"label": "panel", "probes": str(probes), "format": "bed12"},
        ],
    )

    specs, include_capture_in_names = _suite_capture_specs(args)

    assert include_capture_in_names
    assert [spec.label for spec in specs] == ["off", "generated", "panel"]
    assert not specs[0].enabled
    assert specs[1].generates_probes
    assert not specs[1].uses_provided_probes
    assert specs[1].fraction == pytest.approx(0.5)
    assert specs[2].uses_provided_probes
    assert not specs[2].generates_probes
    assert specs[2].probes == str(probes)
    assert specs[2].probe_format == "bed12"


def test_capture_sweep_uses_paired_condition_seed():
    seed = capture_paired_condition_seed(42, "none", 0.99, "none")

    assert seed == capture_paired_condition_seed(42, "none", 0.99, "none")
    assert seed != capture_paired_condition_seed(42, "high", 0.99, "none")
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


def test_capture_probe_group_key_ignores_binding_energy_not_geometry():
    base = SuiteCaptureSpec(
        label="weak",
        fraction=0.5,
        probe_length=120,
        probe_density=1.0,
        off_target_weight=1.0,
        binding_per_base=5.0,
        gdna_split_penalty=0.2,
        min_overlap=1,
    )
    stronger = SuiteCaptureSpec(
        label="strong",
        fraction=0.5,
        probe_length=120,
        probe_density=1.0,
        off_target_weight=1.0,
        binding_per_base=50.0,
        gdna_split_penalty=0.2,
        min_overlap=1,
    )
    sparser = SuiteCaptureSpec(
        label="sparse",
        fraction=0.5,
        probe_length=120,
        probe_density=0.5,
        off_target_weight=1.0,
        binding_per_base=5.0,
        gdna_split_penalty=0.2,
        min_overlap=1,
    )

    assert capture_probe_group_key(base) == capture_probe_group_key(stronger)
    assert capture_probe_group_key(base) != capture_probe_group_key(sparser)
