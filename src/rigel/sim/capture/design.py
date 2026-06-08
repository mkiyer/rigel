"""Hybrid-capture probe *design* — generate synthetic capture panels for the simulator.

Tile probes over captured transcripts (transcript-coordinate TSV + genomic BED12), avoiding
already-placed genomic probes. This is the design side of hybrid capture; the runtime sampler that
*consumes* a panel is :mod:`capture.sampler`, and the panel config is :mod:`capture.config`.
(Moved out of ``suite.py`` — design is capture logic, not suite orchestration; the suite-config
layer that wires its parameters in stays in ``suite.py``.)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ...transcript import Transcript
from ...types import Strand
from ..bam import transcript_to_genomic_blocks
from ..intervals import merge_intervals, project_genomic_block_to_transcript

__all__ = [
    "CaptureProbeDesignResult",
    "DesignedGenomicProbe",
    "design_capture_probe_intervals",
    "design_capture_probe_intervals_in_open_regions",
    "write_random_capture_probes",
]


@dataclass(frozen=True)
class CaptureProbeDesignResult:
    """Summary of generated random-transcriptome capture probes."""

    path: Path
    n_transcripts: int
    n_eligible: int
    n_captured: int
    n_probes: int
    bed_path: Path | None = None
    n_eligible_genes: int = 0
    n_captured_genes: int = 0


@dataclass(frozen=True)
class DesignedGenomicProbe:
    """Generated probe projected into genomic BED12 block coordinates."""

    ref: str
    strand: Strand
    blocks: tuple[tuple[int, int], ...]


def design_capture_probe_intervals(
    transcript_length: int,
    probe_length: int = 120,
    probe_density: float = 1.0,
) -> list[tuple[int, int]]:
    """Return non-overlapping, centered probe intervals on transcript coordinates."""
    if probe_length <= 0:
        raise ValueError("probe_length must be > 0")
    if not 0.0 <= probe_density <= 1.0:
        raise ValueError("probe_density must be between 0 and 1")
    if transcript_length < probe_length or probe_density <= 0:
        return []

    max_probes = int(transcript_length) // int(probe_length)
    n_probes = int(np.floor(max_probes * probe_density))
    n_probes = max(1, n_probes)
    total_probe_bases = n_probes * probe_length
    total_gap = int(transcript_length) - total_probe_bases
    gap = total_gap / (n_probes + 1)

    intervals: list[tuple[int, int]] = []
    previous_end = 0
    for i in range(n_probes):
        start = int(round(gap * (i + 1) + probe_length * i))
        start = max(start, previous_end)
        start = min(start, int(transcript_length) - probe_length)
        end = start + probe_length
        if start < previous_end:
            raise RuntimeError("probe tiling produced overlapping probes")
        intervals.append((start, end))
        previous_end = end
    return intervals


def _subtract_masked_intervals(
    transcript_length: int,
    masked_intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    open_intervals: list[tuple[int, int]] = []
    cursor = 0
    for start, end in merge_intervals(masked_intervals):
        start = max(0, min(int(transcript_length), int(start)))
        end = max(0, min(int(transcript_length), int(end)))
        if end <= start:
            continue
        if cursor < start:
            open_intervals.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < transcript_length:
        open_intervals.append((cursor, int(transcript_length)))
    return open_intervals


def _project_genomic_probe_to_transcript(
    transcript: Transcript,
    probe: DesignedGenomicProbe,
) -> list[tuple[int, int]]:
    if transcript.ref != probe.ref or transcript.strand != probe.strand:
        return []
    projected: list[tuple[int, int]] = []
    for start, end in probe.blocks:
        intervals = project_genomic_block_to_transcript(transcript, start, end)
        if intervals is None:
            return []
        projected.extend(intervals)
    return projected


def _masked_probe_intervals(
    transcript: Transcript,
    existing_probes: list[DesignedGenomicProbe],
) -> list[tuple[int, int]]:
    masked: list[tuple[int, int]] = []
    for probe in existing_probes:
        masked.extend(_project_genomic_probe_to_transcript(transcript, probe))
    return merge_intervals(masked)


def design_capture_probe_intervals_in_open_regions(
    transcript_length: int,
    probe_length: int,
    probe_density: float,
    masked_intervals: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Tile probes only in unmasked transcript-space intervals.

    Existing genomic probes are treated as already covering their projected
    transcript coordinates.  The requested density is therefore a soft target:
    fragmented open space can make the realized design sparser than requested.
    """
    if not masked_intervals:
        return design_capture_probe_intervals(
            transcript_length,
            probe_length=probe_length,
            probe_density=probe_density,
        )

    open_intervals = _subtract_masked_intervals(transcript_length, masked_intervals)
    candidates: list[tuple[int, int]] = []
    max_probes = 0
    for open_start, open_end in open_intervals:
        open_len = open_end - open_start
        max_probes += open_len // probe_length
        for start, end in design_capture_probe_intervals(
            open_len,
            probe_length=probe_length,
            probe_density=1.0,
        ):
            candidates.append((open_start + start, open_start + end))

    if max_probes <= 0 or not candidates:
        return []
    n_probes = int(np.floor(max_probes * probe_density))
    n_probes = min(len(candidates), max(1, n_probes))
    if n_probes >= len(candidates):
        return candidates

    step = len(candidates) / n_probes
    selected_indices = [int((i + 0.5) * step) for i in range(n_probes)]
    return [candidates[index] for index in selected_indices]


def _bed12_row_from_blocks(
    ref: str,
    strand: Strand,
    name: str,
    blocks: list[tuple[int, int]],
) -> str | None:
    if not blocks:
        return None
    chrom_start = min(block_start for block_start, _ in blocks)
    chrom_end = max(block_end for _, block_end in blocks)
    block_sizes = [block_end - block_start for block_start, block_end in blocks]
    block_starts = [block_start - chrom_start for block_start, _ in blocks]
    return "\t".join([
        str(ref),
        str(chrom_start),
        str(chrom_end),
        name,
        "0",
        strand.to_str(),
        str(chrom_start),
        str(chrom_end),
        "0",
        str(len(blocks)),
        ",".join(str(size) for size in block_sizes),
        ",".join(str(offset) for offset in block_starts),
    ])


def _bed12_row_for_probe(
    transcript: Transcript,
    start: int,
    end: int,
    probe_index: int,
) -> str | None:
    if transcript.ref is None or transcript.t_id is None:
        return None
    blocks = transcript_to_genomic_blocks(start, end, transcript)
    return _bed12_row_from_blocks(
        str(transcript.ref),
        transcript.strand,
        f"{transcript.t_id}:probe_{probe_index}",
        blocks,
    )


def write_random_capture_probes(
    transcripts: list[Transcript],
    path: Path,
    *,
    capture_fraction: float,
    probe_length: int = 120,
    probe_density: float = 1.0,
    seed: int = 42,
    bed_path: Path | None = None,
) -> CaptureProbeDesignResult:
    """Select captured gene groups and write transcript-coordinate TSV plus BED12."""
    if not 0.0 <= capture_fraction <= 1.0:
        raise ValueError("capture_fraction must be between 0 and 1")
    if probe_length <= 0:
        raise ValueError("probe_length must be > 0")
    if not 0.0 <= probe_density <= 1.0:
        raise ValueError("probe_density must be between 0 and 1")
    bed_path = bed_path or path.with_suffix(".bed")
    path.parent.mkdir(parents=True, exist_ok=True)
    bed_path.parent.mkdir(parents=True, exist_ok=True)

    eligible: list[tuple[int, Transcript, int]] = []
    for idx, transcript in enumerate(transcripts):
        length = int(transcript.length or transcript.compute_length())
        if length >= probe_length and transcript.t_id is not None and transcript.ref is not None:
            eligible.append((idx, transcript, length))

    eligible_by_gene: dict[str, list[tuple[int, Transcript, int]]] = {}
    for item in eligible:
        idx, transcript, _ = item
        gene_key = str(transcript.g_id or transcript.t_id or idx)
        eligible_by_gene.setdefault(gene_key, []).append(item)
    eligible_genes = list(eligible_by_gene)

    if capture_fraction <= 0.0 or probe_density <= 0.0 or not eligible:
        path.write_text("transcript_id\tstart\tend\n")
        bed_path.write_text("")
        return CaptureProbeDesignResult(
            path,
            len(transcripts),
            len(eligible),
            0,
            0,
            bed_path,
            len(eligible_genes),
            0,
        )

    n_captured_genes = int(np.floor(len(eligible_genes) * capture_fraction))
    n_captured_genes = min(len(eligible_genes), max(1, n_captured_genes))
    rng = np.random.default_rng(seed)
    selected_positions = set(rng.choice(len(eligible_genes), size=n_captured_genes, replace=False))
    selected_gene_keys = {
        gene_key for pos, gene_key in enumerate(eligible_genes) if pos in selected_positions
    }
    selected = [
        item
        for gene_key in eligible_genes
        if gene_key in selected_gene_keys
        for item in eligible_by_gene[gene_key]
    ]
    selected.sort(
        key=lambda item: (
            -float(item[1].abundance or 0.0),
            str(item[1].t_id),
            item[0],
        ),
    )

    n_probes = 0
    existing_probes: list[DesignedGenomicProbe] = []
    with open(path, "w") as handle, open(bed_path, "w") as bed_handle:
        handle.write("transcript_id\tstart\tend\n")
        for _, transcript, length in selected:
            masked = _masked_probe_intervals(transcript, existing_probes)
            for start, end in design_capture_probe_intervals_in_open_regions(
                length,
                probe_length=probe_length,
                probe_density=probe_density,
                masked_intervals=masked,
            ):
                blocks = transcript_to_genomic_blocks(start, end, transcript)
                if not blocks or transcript.ref is None:
                    continue
                n_probes += 1
                probe_name = f"{transcript.t_id}:probe_{n_probes}"
                bed_row = _bed12_row_from_blocks(
                    str(transcript.ref),
                    transcript.strand,
                    probe_name,
                    blocks,
                )
                if bed_row is None:
                    continue
                handle.write(f"{transcript.t_id}\t{start}\t{end}\n")
                bed_handle.write(f"{bed_row}\n")
                existing_probes.append(DesignedGenomicProbe(
                    str(transcript.ref),
                    transcript.strand,
                    tuple(blocks),
                ))

    return CaptureProbeDesignResult(
        path,
        len(transcripts),
        len(eligible),
        len(selected),
        n_probes,
        bed_path,
        len(eligible_genes),
        n_captured_genes,
    )
