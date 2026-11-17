"""Design a hybrid-capture probe panel over a suite reference: the panel's density step.

Over a run of flat regions a uniform scenario cannot tell "the messages work" from "the global prior reached
it". A panel that covers some gene groups and not others gives a sharp density step: the captured / uncaptured
boundary is a cliff in gDNA density with real transcripts on both sides. The whole-genome simulator requires a
panel when capture is enabled and does not design one (`require_probes_when_enabled`), so this is a
prerequisite of the capture axis. `capture_fraction` is a scenario parameter like `gdna_rate`: the share of the
annotation that is on-panel.

A random `capture_fraction` of the gene groups is captured, and their transcripts are tiled most abundant first,
each probe skipping the spans an earlier probe already covers. The panel is written twice: as a transcript-coordinate
TSV and as a genomic BED12 beside it (`<out>.bed`); `rigel.sim.capture.CaptureSampler` reads either.

    python scripts/sim/design_suite_probes.py --gtf genes.gtf -o panel.tsv --capture-fraction 0.5
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from rigel.sim.bam import transcript_to_genomic_blocks
from rigel.sim.intervals import clip_blocks_to_transcript, merge_intervals
from rigel.transcript import Transcript
from rigel.types import Strand

#: The default bait length in bp: a standard hybrid-capture bait.
PROBE_LENGTH_BP = 120


@dataclass(frozen=True)
class CaptureProbeDesignResult:
    """Summary of a designed panel."""

    path: Path
    n_transcripts: int
    n_eligible: int
    n_captured: int
    n_probes: int
    bed_path: Path
    n_eligible_genes: int
    n_captured_genes: int


@dataclass(frozen=True)
class DesignedGenomicProbe:
    """A placed probe in genomic BED12 block coordinates."""

    ref: str
    strand: Strand
    blocks: tuple[tuple[int, int], ...]


# ── tiling one transcript ────────────────────────────────────────────────────────────────────────


def design_capture_probe_intervals(
    transcript_length: int,
    probe_length: int = PROBE_LENGTH_BP,
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
    n_probes = max(1, int(np.floor(max_probes * probe_density)))
    gap = (int(transcript_length) - n_probes * probe_length) / (n_probes + 1)

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


def _open_intervals(transcript_length: int, masked: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The stretches of ``[0, transcript_length)`` that ``masked`` leaves uncovered."""
    open_intervals: list[tuple[int, int]] = []
    cursor = 0
    for start, end in merge_intervals(masked):
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


def _tile_open_regions(
    transcript_length: int,
    probe_length: int,
    probe_density: float,
    masked: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Tile probes only in the transcript-space stretches ``masked`` leaves open. The density is a soft target:
    fragmented open space can make the design sparser than requested."""
    if not masked:
        return design_capture_probe_intervals(
            transcript_length, probe_length=probe_length, probe_density=probe_density
        )

    candidates: list[tuple[int, int]] = []
    max_probes = 0
    for open_start, open_end in _open_intervals(transcript_length, masked):
        max_probes += (open_end - open_start) // probe_length
        for start, end in design_capture_probe_intervals(
            open_end - open_start, probe_length=probe_length, probe_density=1.0
        ):
            candidates.append((open_start + start, open_start + end))

    if max_probes <= 0 or not candidates:
        return []
    n_probes = min(len(candidates), max(1, int(np.floor(max_probes * probe_density))))
    if n_probes >= len(candidates):
        return candidates
    step = len(candidates) / n_probes
    return [candidates[int((i + 0.5) * step)] for i in range(n_probes)]


def _covered_by(transcript: Transcript, probes: list[DesignedGenomicProbe]) -> list[tuple[int, int]]:
    """The transcript-space spans the placed ``probes`` already cover: each probe that lies wholly in the
    transcript's exons, on its strand."""
    covered: list[tuple[int, int]] = []
    for probe in probes:
        if transcript.ref != probe.ref or transcript.strand != probe.strand:
            continue
        landed = clip_blocks_to_transcript(transcript, probe.blocks)
        if sum(end - start for start, end in landed) == sum(end - start for start, end in probe.blocks):
            covered.extend(landed)
    return merge_intervals(covered)


def _bed12_row(ref: str, strand: Strand, name: str, blocks: list[tuple[int, int]]) -> str:
    chrom_start = min(block_start for block_start, _ in blocks)
    chrom_end = max(block_end for _, block_end in blocks)
    sizes = ",".join(str(block_end - block_start) for block_start, block_end in blocks)
    offsets = ",".join(str(block_start - chrom_start) for block_start, _ in blocks)
    fields = (ref, chrom_start, chrom_end, name, 0, strand.to_str(), chrom_start, chrom_end, 0,
              len(blocks), sizes, offsets)
    return "\t".join(str(field) for field in fields)


# ── the panel ────────────────────────────────────────────────────────────────────────────────────


def _eligible_by_gene(transcripts: list[Transcript], probe_length: int) -> dict[str, list]:
    """``(index, transcript, length)`` of every transcript a probe fits, grouped by gene in first-seen order."""
    by_gene: dict[str, list] = {}
    for idx, transcript in enumerate(transcripts):
        length = int(transcript.length or transcript.compute_length())
        if length >= probe_length and transcript.t_id is not None and transcript.ref is not None:
            gene = str(transcript.g_id or transcript.t_id or idx)
            by_gene.setdefault(gene, []).append((idx, transcript, length))
    return by_gene


def _captured(by_gene: dict[str, list], capture_fraction: float, seed: int) -> tuple[list, int]:
    """The transcripts of a random ``capture_fraction`` of the gene groups, most abundant first, and the number of
    groups captured."""
    genes = list(by_gene)
    n_genes = min(len(genes), max(1, int(np.floor(len(genes) * capture_fraction))))
    positions = set(np.random.default_rng(seed).choice(len(genes), size=n_genes, replace=False))
    captured = [item for pos, gene in enumerate(genes) if pos in positions for item in by_gene[gene]]
    captured.sort(key=lambda item: (-float(item[1].abundance or 0.0), str(item[1].t_id), item[0]))
    return captured, n_genes


def _write_probes(selected: list, path: Path, bed_path: Path, probe_length: int, probe_density: float) -> int:
    """Tile each selected transcript in turn around the probes already placed; return the number written."""
    n_probes = 0
    placed: list[DesignedGenomicProbe] = []
    with open(path, "w") as handle, open(bed_path, "w") as bed_handle:
        handle.write("transcript_id\tstart\tend\n")
        for _, transcript, length in selected:
            masked = _covered_by(transcript, placed)
            for start, end in _tile_open_regions(length, probe_length, probe_density, masked):
                blocks = transcript_to_genomic_blocks(start, end, transcript)
                if not blocks:
                    continue
                n_probes += 1
                name = f"{transcript.t_id}:probe_{n_probes}"
                handle.write(f"{transcript.t_id}\t{start}\t{end}\n")
                bed_handle.write(_bed12_row(str(transcript.ref), transcript.strand, name, blocks) + "\n")
                placed.append(DesignedGenomicProbe(str(transcript.ref), transcript.strand, tuple(blocks)))
    return n_probes


def write_random_capture_probes(
    transcripts: list[Transcript],
    path: Path,
    *,
    capture_fraction: float,
    probe_length: int = PROBE_LENGTH_BP,
    probe_density: float = 1.0,
    seed: int = 42,
) -> CaptureProbeDesignResult:
    """Capture a random ``capture_fraction`` of genes and write their probes as a transcript-coordinate TSV at
    ``path`` and a BED12 beside it, tiling the most abundant transcripts first."""
    if not 0.0 <= capture_fraction <= 1.0:
        raise ValueError("capture_fraction must be between 0 and 1")
    if probe_length <= 0:
        raise ValueError("probe_length must be > 0")
    if not 0.0 <= probe_density <= 1.0:
        raise ValueError("probe_density must be between 0 and 1")
    bed_path = path.with_suffix(".bed")
    if bed_path == path:
        # the BED12 is written beside the TSV, so a `.bed` target is its own sibling: two handles on one file
        raise ValueError(f"the panel's target is its TSV, and {path} is a .bed: the BED12 is written beside it")
    path.parent.mkdir(parents=True, exist_ok=True)

    by_gene = _eligible_by_gene(transcripts, probe_length)
    n_eligible = sum(len(items) for items in by_gene.values())
    if capture_fraction <= 0.0 or probe_density <= 0.0 or not n_eligible:
        path.write_text("transcript_id\tstart\tend\n")
        bed_path.write_text("")
        return CaptureProbeDesignResult(path, len(transcripts), n_eligible, 0, 0, bed_path, len(by_gene), 0)

    selected, n_genes = _captured(by_gene, capture_fraction, seed)
    n_probes = _write_probes(selected, path, bed_path, probe_length, probe_density)
    return CaptureProbeDesignResult(
        path, len(transcripts), n_eligible, len(selected), n_probes, bed_path, len(by_gene), n_genes
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--capture-fraction", type=float, required=True, help="share of GENE GROUPS on panel")
    ap.add_argument("--probe-length", type=int, default=PROBE_LENGTH_BP)
    ap.add_argument("--probe-density", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    transcripts = Transcript.read_gtf(str(args.gtf), parse_mode="warn-skip")
    result = write_random_capture_probes(
        transcripts,
        args.out,
        capture_fraction=args.capture_fraction,
        probe_length=args.probe_length,
        probe_density=args.probe_density,
        seed=args.seed,
    )
    print(f"transcripts       {len(transcripts):,}")
    print(f"eligible          {result.n_eligible:,} over {result.n_eligible_genes:,} gene groups")
    print(f"gene groups ON    {result.n_captured_genes:,}  "
          f"({100 * result.n_captured_genes / max(result.n_eligible_genes, 1):.1f} %)")
    print(f"probes            {result.n_probes:,}")
    print(f"written           {result.path}")


if __name__ == "__main__":
    main()
