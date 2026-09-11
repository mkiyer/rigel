"""gDNA comes from GENOMIC references, from at least two of them, and never runs past the end of one.

Deriving the gDNA reference set from the annotation puts "has an annotation" in place of "is
genomic", and an RNA-only spike-in reference carrying one transcript then qualifies — so the panel
holds gDNA molecules on spike-ins, which cannot exist, on exactly the references whose truth
abundance is the false-positive control (`TRAPS: annotated-is-not-genomic`). The fixture is built so
that the annotation proxy would PASS it: its spike-in reference is ANNOTATED, so the two readings
disagree on it, where a fixture whose RNA-only reference had no annotation would be green with the
defect present. Requiring at least TWO genomic references is not decoration either, because a
single-reference synthetic index hides a reference-id-space mismatch that drops nearly every
fragment inside ``deposit()`` while the golden tests pass
(`TRAPS: one-reference-hides-refid-bugs`), and the gDNA intergenic branch is its own path through
the scanner.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.sim.genome import random_dna_array
from rigel.sim.read_name import parse_origin
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


#: The RNA-only reference's name. Deliberately not "ERCC-*": nothing about the defect is
#: ERCC-specific — any spike-in, decoy or synthetic-transcript reference hits the same boundary.
SPIKE_REF = "SPIKE-1"

GENOMIC_REFS = ["chrA", "chrB"]


def _write_fasta(tmp_path: Path, sequences: dict[str, str]) -> Path:
    path = tmp_path / "genome.fa"
    with open(path, "w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def _transcript(t_id: str, ref: str, exons: list[tuple[int, int]], abundance: float) -> Transcript:
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


@pytest.fixture
def mixed_reference(tmp_path):
    """Two genomic references plus one ANNOTATED RNA-only spike-in reference."""
    rng = np.random.default_rng(7)
    lengths = {"chrA": 40_000, "chrB": 30_000, SPIKE_REF: 900}
    sequences = {name: "".join(random_dna_array(n, rng)) for name, n in lengths.items()}
    fasta = _write_fasta(tmp_path, sequences)
    transcripts = [
        _transcript("A1", "chrA", [(2_000, 4_000), (6_000, 8_000)], 100.0),
        _transcript("A2", "chrA", [(20_000, 24_000)], 50.0),
        _transcript("B1", "chrB", [(3_000, 5_000), (9_000, 11_000)], 100.0),
        # The spike-in: one transcript covering nearly the whole reference, as ERCC does.
        _transcript("SPIKE_T1", SPIKE_REF, [(30, 860)], 200.0),
    ]
    return fasta, transcripts, lengths


def _simulate_gdna(tmp_path, fasta, transcripts, *, genomic_refs, n_gdna, seed=19):
    simulator = WholeGenomeSimulator(
        fasta,
        transcripts,
        SimulationParams(
            sim_seed=seed,
            frag_mean=157,
            frag_std=125,
            frag_min=50,
            frag_max=800,
            read_length=100,
        ),
        GDNASimConfig(frag_mean=157, frag_std=125, frag_min=40, frag_max=1000),
        seed=seed,
        genomic_refs=genomic_refs,
    )
    try:
        _, _, bam_path = simulator.simulate_and_write(
            tmp_path / "out",
            n_rna=0,
            n_gdna=n_gdna,
            oracle_bam=True,
            n_workers=1,
        )
    finally:
        simulator.close()

    origins = []
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for read in bam:
            if read.is_read1:
                origins.append(parse_origin(read.query_name))
    return origins


class TestGenomicReferenceSet:
    def test_g_s1_no_gdna_on_an_rna_only_reference(self, tmp_path, mixed_reference):
        """An absolute count, from the read names, and the target is exactly 0."""
        fasta, transcripts, _lengths = mixed_reference
        origins = _simulate_gdna(
            tmp_path, fasta, transcripts, genomic_refs=GENOMIC_REFS, n_gdna=4_000
        )

        gdna = [origin for origin in origins if origin.kind == "gdna"]
        assert len(gdna) == 4_000, (
            "every requested gDNA fragment must be written, none silently lost"
        )
        on_spike = [origin for origin in gdna if origin.ref == SPIKE_REF]
        assert on_spike == [], f"{len(on_spike)} gDNA fragments on an RNA-only reference"

    def test_g_s2_every_genomic_reference_carries_gdna(self, tmp_path, mixed_reference):
        """At least two genomic references, each with a non-zero deposit census."""
        fasta, transcripts, _lengths = mixed_reference
        origins = _simulate_gdna(
            tmp_path, fasta, transcripts, genomic_refs=GENOMIC_REFS, n_gdna=4_000
        )

        per_ref: dict[str, int] = {}
        for origin in origins:
            if origin.kind == "gdna":
                per_ref[str(origin.ref)] = per_ref.get(str(origin.ref), 0) + 1

        assert len(GENOMIC_REFS) >= 2
        assert set(per_ref) == set(GENOMIC_REFS), f"gDNA reference census: {per_ref}"
        assert all(count > 0 for count in per_ref.values()), per_ref

    def test_g_s6_no_gdna_fragment_runs_past_its_own_reference(self, tmp_path, mixed_reference):
        """G-S6 — a regression guard. A 400 bp molecule cannot come from a 273 bp contig."""
        fasta, transcripts, lengths = mixed_reference
        origins = _simulate_gdna(
            tmp_path, fasta, transcripts, genomic_refs=GENOMIC_REFS, n_gdna=4_000
        )

        violations = [
            origin
            for origin in origins
            if origin.kind == "gdna" and int(origin.end) > lengths[str(origin.ref)]
        ]
        assert violations == [], f"{len(violations)} gDNA fragments run past their reference end"


class TestTheClassificationIsAnInput:
    """The set is stated, never inferred — and a mis-statement must be loud, not silent."""

    def test_an_annotated_reference_left_out_of_the_set_gets_no_gdna(
        self, tmp_path, mixed_reference
    ):
        """The spike-in is annotated, so an annotation-derived reference set would include it."""
        fasta, transcripts, _lengths = mixed_reference
        simulator = WholeGenomeSimulator(
            fasta,
            transcripts,
            SimulationParams(sim_seed=3),
            GDNASimConfig(),
            genomic_refs=GENOMIC_REFS,
        )
        try:
            assert simulator.genomic_refs == tuple(GENOMIC_REFS)
            assert SPIKE_REF not in simulator.genomic_refs
            assert SPIKE_REF in {t.ref for t in transcripts}
        finally:
            simulator.close()

    def test_an_unknown_reference_name_raises(self, tmp_path, mixed_reference):
        """A typo must not silently produce zero gDNA, which is a request answered with nothing."""
        fasta, transcripts, _lengths = mixed_reference
        with pytest.raises(ValueError, match="not in the reference"):
            WholeGenomeSimulator(
                fasta,
                transcripts,
                SimulationParams(sim_seed=3),
                GDNASimConfig(),
                genomic_refs=["chrA", "chr22"],
            )

    def test_requesting_gdna_with_no_genomic_reference_raises(self, tmp_path, mixed_reference):
        """Silently writing zero of a requested five million fragments is unrepresentable."""
        fasta, transcripts, _lengths = mixed_reference
        simulator = WholeGenomeSimulator(
            fasta,
            transcripts,
            SimulationParams(sim_seed=3),
            GDNASimConfig(),
            genomic_refs=[],
        )
        try:
            assert simulator._accumulate_gdna_counts(0) == {}
            with pytest.raises(ValueError, match="no genomic reference"):
                simulator._accumulate_gdna_counts(1_000)
        finally:
            simulator.close()
