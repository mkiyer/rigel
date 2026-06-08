"""Hybrid-capture runtime sampler + probe loading.

``CaptureSampler`` builds a sparse probe representation (per transcript / pre-mRNA / gDNA) and
answers capture-aware effective lengths, fragment-start sampling, and per-fragment weights. The
model is described in :mod:`capture.config`. Probe *design* (generating synthetic panels) lives in
:mod:`capture.design`.
"""

from __future__ import annotations

import gzip
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from ...transcript import Transcript
from ...types import Strand
from ..bam import transcript_to_genomic_blocks
from ..intervals import merge_intervals, project_genomic_blocks_to_transcript
from .config import CaptureConfig

logger = logging.getLogger(__name__)

__all__ = ["CaptureSampler", "WeightedInterval"]


@dataclass(frozen=True, slots=True)
class WeightedInterval:
    """A probe interval on one coordinate axis with a binding multiplier."""

    start: int
    end: int
    scale: float = 1.0
    probe_group: int = 0

    @property
    def length(self) -> int:
        return self.end - self.start


class CaptureSampler:
    """Sparse capture-aware partition functions and start sampling."""

    def __init__(
        self,
        config: CaptureConfig,
        transcripts: Sequence[Transcript],
        ref_lengths: Mapping[str, int],
        *,
        enabled: bool,
    ):
        self.config = config
        self.transcripts = list(transcripts)
        self.ref_lengths = dict(ref_lengths)
        self._enabled = enabled

        self._tx_lengths = np.array(
            [int(t.length or t.compute_length()) for t in transcripts],
            dtype=np.int64,
        )
        self._premrna_lengths = np.array(
            [int(t.end - t.start) for t in transcripts],
            dtype=np.int64,
        )
        self._tx_id_to_index = {
            str(t.t_id): idx for idx, t in enumerate(transcripts) if t.t_id is not None
        }
        self._tx_by_ref: dict[str, list[int]] = defaultdict(list)
        for idx, t in enumerate(transcripts):
            if t.ref is not None:
                self._tx_by_ref[str(t.ref)].append(idx)

        self._mrna_intervals: dict[int, list[WeightedInterval]] = defaultdict(list)
        self._nrna_intervals: dict[int, list[WeightedInterval]] = defaultdict(list)
        self._gdna_intervals: dict[str, list[WeightedInterval]] = defaultdict(list)
        self._next_probe_group = 1

        self._mass_cache: dict[tuple[str, int | str, int, int], tuple[np.ndarray, np.ndarray]] = {}

    @classmethod
    def disabled(
        cls,
        transcripts: Sequence[Transcript] = (),
        ref_lengths: Mapping[str, int] | None = None,
    ) -> "CaptureSampler":
        """Return a disabled sampler with uniform sampling semantics."""
        return cls(
            CaptureConfig(),
            transcripts,
            ref_lengths or {},
            enabled=False,
        )

    @classmethod
    def from_config(
        cls,
        config: CaptureConfig | None,
        transcripts: Sequence[Transcript],
        ref_lengths: Mapping[str, int],
    ) -> "CaptureSampler":
        """Build a sampler from a probe file, or return disabled if absent."""
        if config is None or not config.probes:
            return cls.disabled(transcripts, ref_lengths)
        _validate_config(config)
        sampler = cls(config, transcripts, ref_lengths, enabled=True)
        sampler._load_probe_file(Path(config.probes))
        logger.info(
            "Hybrid capture enabled: %d mRNA targets, %d nRNA targets, %d gDNA refs",
            len(sampler._mrna_intervals),
            len(sampler._nrna_intervals),
            len(sampler._gdna_intervals),
        )
        return sampler

    @property
    def enabled(self) -> bool:
        return self._enabled

    def partition(
        self,
        space: str,
        key: int | str,
        seq_len: int,
        frag_len: int,
    ) -> float:
        """Return the capture-aware effective length for one template."""
        eff_len = int(seq_len) - int(frag_len) + 1
        if eff_len <= 0:
            return 0.0

        baseline = self.config.off_target_weight * eff_len
        if not self.enabled or self.config.binding_per_base <= 0:
            return float(baseline)

        _, weights = self._extra_landscape(space, key, seq_len, frag_len)
        return float(baseline + self.config.binding_per_base * weights.sum())

    def partition_array(
        self,
        space: str,
        keys: Iterable[int | str],
        lengths: np.ndarray,
        frag_len: int,
    ) -> np.ndarray:
        """Vector of capture-aware effective lengths for a pool."""
        eff = np.maximum(0, lengths.astype(np.int64) - int(frag_len) + 1)
        result = eff.astype(np.float64) * self.config.off_target_weight
        if not self.enabled or self.config.binding_per_base <= 0:
            return result

        for pos, key in enumerate(keys):
            if eff[pos] <= 0 or not self._get_intervals(space, key):
                continue
            result[pos] = self.partition(space, key, int(lengths[pos]), frag_len)
        return result

    def sample_starts(
        self,
        space: str,
        key: int | str,
        seq_len: int,
        frag_len: int,
        count: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Sample fragment starts from the sparse capture landscape."""
        eff_len = int(seq_len) - int(frag_len) + 1
        if count <= 0:
            return np.empty(0, dtype=np.int64)
        if eff_len <= 0:
            return np.empty(0, dtype=np.int64)

        if not self.enabled or self.config.binding_per_base <= 0:
            return rng.integers(0, eff_len, size=count, dtype=np.int64)

        extra_starts, raw_weights = self._extra_landscape(space, key, seq_len, frag_len)
        raw_total = float(raw_weights.sum())
        baseline_mass = self.config.off_target_weight * eff_len
        extra_mass = self.config.binding_per_base * raw_total
        total_mass = baseline_mass + extra_mass
        if raw_total <= 0 or total_mass <= 0:
            return rng.integers(0, eff_len, size=count, dtype=np.int64)

        starts = np.empty(count, dtype=np.int64)
        baseline_mask = rng.random(count) < (baseline_mass / total_mass)
        n_baseline = int(baseline_mask.sum())
        if n_baseline > 0:
            starts[baseline_mask] = rng.integers(
                0,
                eff_len,
                size=n_baseline,
                dtype=np.int64,
            )

        extra_positions = np.where(~baseline_mask)[0]
        if len(extra_positions) == 0:
            return starts

        starts[extra_positions] = rng.choice(
            extra_starts,
            size=len(extra_positions),
            p=raw_weights / raw_total,
        )

        return starts

    def fragment_weight(
        self,
        space: str,
        key: int | str,
        seq_len: int,
        frag_start: int,
        frag_len: int,
    ) -> float:
        """Return the unnormalized sampling weight for a concrete fragment."""
        eff_len = int(seq_len) - int(frag_len) + 1
        if eff_len <= 0 or frag_start < 0 or frag_start >= eff_len:
            return 0.0
        group_scores: dict[int, float] = {}
        frag_end = frag_start + frag_len
        for interval in self._get_intervals(space, key):
            overlap = min(frag_end, interval.end) - max(frag_start, interval.start)
            if overlap >= self.config.min_overlap:
                group_scores[interval.probe_group] = (
                    group_scores.get(interval.probe_group, 0.0) + interval.scale * overlap
                )
        score = max(group_scores.values(), default=0.0)
        return self.config.off_target_weight + self.config.binding_per_base * score

    # -- Probe loading -----------------------------------------------------

    def _load_probe_file(self, path: Path) -> None:
        lines = _read_probe_lines(path)
        if not lines:
            raise ValueError(f"capture probe file is empty: {path}")
        probe_format = _detect_probe_format(lines, self.config.probe_format)
        if probe_format == "transcript":
            self._load_transcript_probes(lines, path)
        elif probe_format == "bed12":
            self._load_bed12_probes(lines, path)
        else:
            raise ValueError(f"unsupported capture probe format: {probe_format!r}")

    def _load_transcript_probes(self, lines: list[str], path: Path) -> None:
        header = _transcript_header(_split_fields(lines[0]))
        start_idx = 1 if header is not None else 0
        n_loaded = 0
        for row_number, line in enumerate(lines[start_idx:], start=start_idx + 1):
            fields = _split_fields(line)
            if header is None:
                if len(fields) < 3:
                    raise ValueError(f"{path}:{row_number}: expected transcript_id start end")
                transcript_id, start_text, end_text = fields[:3]
            else:
                transcript_id = fields[header["transcript_id"]]
                start_text = fields[header["start"]]
                end_text = fields[header["end"]]
            self._add_transcript_probe(
                transcript_id,
                int(start_text),
                int(end_text),
                path,
                row_number,
            )
            n_loaded += 1
        logger.info("Loaded %d transcript-coordinate capture probes from %s", n_loaded, path)

    def _load_bed12_probes(self, lines: list[str], path: Path) -> None:
        n_loaded = 0
        for row_number, line in enumerate(lines, start=1):
            fields = _split_fields(line)
            if len(fields) < 12:
                raise ValueError(f"{path}:{row_number}: expected 12 BED fields")
            if not _is_int(fields[1]) or not _is_int(fields[2]):
                continue
            ref = fields[0]
            chrom_start = int(fields[1])
            strand = Strand.from_str(fields[5]) if fields[5] in {"+", "-"} else Strand.NONE
            block_count = int(fields[9])
            block_sizes = _parse_bed_list(fields[10])
            block_starts = _parse_bed_list(fields[11])
            if block_count != len(block_sizes) or block_count != len(block_starts):
                raise ValueError(f"{path}:{row_number}: BED12 blockCount mismatch")
            blocks = [
                (chrom_start + rel_start, chrom_start + rel_start + size)
                for size, rel_start in zip(block_sizes, block_starts)
                if size > 0
            ]
            if not blocks:
                continue
            self._add_bed12_probe(ref, strand, blocks)
            n_loaded += 1
        logger.info("Loaded %d BED12 capture probes from %s", n_loaded, path)

    def _add_transcript_probe(
        self,
        transcript_id: str,
        start: int,
        end: int,
        path: Path,
        row_number: int,
    ) -> None:
        t_idx = self._tx_id_to_index.get(transcript_id)
        if t_idx is None:
            logger.debug("Skipping probe for unknown transcript %s", transcript_id)
            return
        tx_len = int(self._tx_lengths[t_idx])
        if start < 0 or end <= start or end > tx_len:
            raise ValueError(
                f"{path}:{row_number}: invalid probe interval for {transcript_id}: "
                f"[{start}, {end}) outside transcript length {tx_len}"
            )

        probe_group = self._new_probe_group()
        self._mrna_intervals[t_idx].append(WeightedInterval(start, end, 1.0, probe_group))
        transcript = self.transcripts[t_idx]
        blocks = transcript_to_genomic_blocks(start, end, transcript)
        split_scale = self._split_scale(blocks)
        self._add_gdna_blocks(str(transcript.ref), blocks, split_scale, probe_group)
        self._add_nrna_blocks(t_idx, blocks, split_scale, probe_group)

    def _add_bed12_probe(
        self,
        ref: str,
        strand: Strand,
        blocks: list[tuple[int, int]],
    ) -> None:
        split_scale = self._split_scale(blocks)
        probe_group = self._new_probe_group()
        self._add_gdna_blocks(ref, blocks, split_scale, probe_group)

        for t_idx in self._tx_by_ref.get(ref, []):
            transcript = self.transcripts[t_idx]
            if strand != Strand.NONE and transcript.strand != strand:
                continue
            projected = project_genomic_blocks_to_transcript(transcript, blocks)
            if projected is None:
                continue
            for start, end in merge_intervals(projected):
                self._mrna_intervals[t_idx].append(WeightedInterval(start, end, 1.0, probe_group))
            self._add_nrna_blocks(t_idx, blocks, split_scale, probe_group)

    def _add_gdna_blocks(
        self,
        ref: str,
        blocks: Sequence[tuple[int, int]],
        scale: float,
        probe_group: int,
    ) -> None:
        ref_len = self.ref_lengths.get(ref)
        if ref_len is None:
            return
        for start, end in blocks:
            interval = _clip_interval(start, end, ref_len, scale, probe_group)
            if interval is not None:
                self._gdna_intervals[ref].append(interval)

    def _add_nrna_blocks(
        self,
        t_idx: int,
        blocks: Sequence[tuple[int, int]],
        scale: float,
        probe_group: int,
    ) -> None:
        transcript = self.transcripts[t_idx]
        seq_len = int(self._premrna_lengths[t_idx])
        for start, end in blocks:
            mapped = _genomic_block_to_premrna_interval(transcript, start, end)
            if mapped is None:
                continue
            interval = _clip_interval(mapped[0], mapped[1], seq_len, scale, probe_group)
            if interval is not None:
                self._nrna_intervals[t_idx].append(interval)

    def _new_probe_group(self) -> int:
        group = self._next_probe_group
        self._next_probe_group += 1
        return group

    def _split_scale(self, blocks: Sequence[tuple[int, int]]) -> float:
        if len(blocks) <= 1:
            return 1.0
        return self.config.gdna_split_penalty

    # -- Sparse math -------------------------------------------------------

    def _get_intervals(self, space: str, key: int | str) -> list[WeightedInterval]:
        if space == "mrna":
            return self._mrna_intervals.get(int(key), [])
        if space == "nrna":
            return self._nrna_intervals.get(int(key), [])
        if space == "gdna":
            return self._gdna_intervals.get(str(key), [])
        raise ValueError(f"unknown capture coordinate space: {space!r}")

    def _extra_landscape(
        self,
        space: str,
        key: int | str,
        seq_len: int,
        frag_len: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return start positions and best single-probe overlap weights."""
        cache_key = (space, key, int(seq_len), int(frag_len))
        cached = self._mass_cache.get(cache_key)
        if cached is not None:
            return cached

        start_chunks: list[np.ndarray] = []
        weight_chunks: list[np.ndarray] = []
        group_chunks: list[np.ndarray] = []
        for interval in self._get_intervals(space, key):
            starts, weights = self._local_overlap_weights(seq_len, frag_len, interval)
            if len(starts) == 0:
                continue
            scaled_weights = weights * interval.scale
            keep = scaled_weights > 0
            if np.any(keep):
                start_chunks.append(starts[keep])
                weight_chunks.append(scaled_weights[keep])
                group_chunks.append(
                    np.full(int(np.sum(keep)), interval.probe_group, dtype=np.int64),
                )

        if not start_chunks:
            empty_starts = np.empty(0, dtype=np.int64)
            empty_weights = np.empty(0, dtype=np.float64)
            self._mass_cache[cache_key] = (empty_starts, empty_weights)
            return empty_starts, empty_weights

        all_starts = np.concatenate(start_chunks)
        all_weights = np.concatenate(weight_chunks)
        all_groups = np.concatenate(group_chunks)

        order = np.lexsort((all_groups, all_starts))
        sorted_starts = all_starts[order]
        sorted_groups = all_groups[order]
        sorted_weights = all_weights[order]
        start_probe_group_starts = np.flatnonzero(
            np.r_[
                True,
                (sorted_starts[1:] != sorted_starts[:-1])
                | (sorted_groups[1:] != sorted_groups[:-1]),
            ],
        )
        probe_starts = sorted_starts[start_probe_group_starts]
        probe_weights = np.add.reduceat(sorted_weights, start_probe_group_starts)

        order = np.argsort(probe_starts)
        sorted_probe_starts = probe_starts[order]
        sorted_probe_weights = probe_weights[order]
        start_group_starts = np.flatnonzero(
            np.r_[True, sorted_probe_starts[1:] != sorted_probe_starts[:-1]],
        )
        starts = sorted_probe_starts[start_group_starts]
        weights = np.maximum.reduceat(sorted_probe_weights, start_group_starts)

        self._mass_cache[cache_key] = (starts, weights)
        return starts, weights

    def _local_overlap_weights(
        self,
        seq_len: int,
        frag_len: int,
        interval: WeightedInterval,
    ) -> tuple[np.ndarray, np.ndarray]:
        eff_len = int(seq_len) - int(frag_len) + 1
        if eff_len <= 0:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64)
        lo = max(0, interval.start - int(frag_len) + 1)
        hi = min(eff_len, interval.end)
        if hi <= lo:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64)
        starts = np.arange(lo, hi, dtype=np.int64)
        overlaps = np.minimum(starts + int(frag_len), interval.end) - np.maximum(
            starts,
            interval.start,
        )
        if self.config.min_overlap > 1:
            overlaps = np.where(overlaps >= self.config.min_overlap, overlaps, 0)
        weights = overlaps.astype(np.float64)
        return starts, weights


def _validate_config(config: CaptureConfig) -> None:
    if config.off_target_weight < 0:
        raise ValueError("capture.off_target_weight must be >= 0")
    if config.binding_per_base < 0:
        raise ValueError("capture.binding_per_base must be >= 0")
    if config.gdna_split_penalty < 0:
        raise ValueError("capture.gdna_split_penalty must be >= 0")
    if config.min_overlap < 1:
        raise ValueError("capture.min_overlap must be >= 1")


def _read_probe_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [
            line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")
        ]


def _split_fields(line: str) -> list[str]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) == 1:
        fields = line.split()
    return fields


def _detect_probe_format(lines: list[str], configured: str) -> str:
    fmt = configured.lower().replace("-", "_")
    aliases = {
        "auto": "auto",
        "transcript": "transcript",
        "transcript_tsv": "transcript",
        "tsv": "transcript",
        "bed12": "bed12",
        "bed": "bed12",
    }
    if fmt not in aliases:
        raise ValueError(
            "capture.probe_format must be one of auto, transcript, transcript_tsv, bed12"
        )
    fmt = aliases[fmt]
    if fmt != "auto":
        return fmt

    fields = _split_fields(lines[0])
    if _transcript_header(fields) is not None:
        return "transcript"
    if len(fields) >= 12 and _is_int(fields[1]) and _is_int(fields[2]) and _is_int(fields[9]):
        return "bed12"
    return "transcript"


def _transcript_header(fields: Sequence[str]) -> dict[str, int] | None:
    lowered = [field.strip().lower() for field in fields]
    tx_col = _find_col(lowered, ("transcript_id", "tx_id", "t_id", "target_id"))
    start_col = _find_col(lowered, ("start", "tx_start", "transcript_start"))
    end_col = _find_col(lowered, ("end", "tx_end", "transcript_end"))
    if tx_col is None or start_col is None or end_col is None:
        return None
    return {"transcript_id": tx_col, "start": start_col, "end": end_col}


def _find_col(fields: Sequence[str], names: Sequence[str]) -> int | None:
    for name in names:
        try:
            return fields.index(name)
        except ValueError:
            continue
    return None


def _is_int(text: str) -> bool:
    try:
        int(text)
    except ValueError:
        return False
    return True


def _parse_bed_list(text: str) -> list[int]:
    return [int(part) for part in text.rstrip(",").split(",") if part]


def _clip_interval(
    start: int,
    end: int,
    length: int,
    scale: float,
    probe_group: int,
) -> WeightedInterval | None:
    start = max(0, int(start))
    end = min(int(length), int(end))
    if end <= start or scale <= 0:
        return None
    return WeightedInterval(start, end, float(scale), int(probe_group))


# Interval helpers (merge_intervals / project_genomic_block(s)_to_transcript) live in
# rigel.sim.intervals — shared with suite.py's capture-probe design.


def _genomic_block_to_premrna_interval(
    transcript: Transcript,
    block_start: int,
    block_end: int,
) -> tuple[int, int] | None:
    if block_start < transcript.start or block_end > transcript.end or block_end <= block_start:
        return None
    if transcript.strand == Strand.NEG:
        premrna_len = transcript.end - transcript.start
        return premrna_len - (block_end - transcript.start), premrna_len - (
            block_start - transcript.start
        )
    return block_start - transcript.start, block_end - transcript.start
