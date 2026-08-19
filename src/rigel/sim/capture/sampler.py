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
from ..intervals import merge_intervals
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
        #: ⭐ Per reference, transcripts sorted by START with the running maximum END beside them, so a
        #: probe finds the transcripts it overlaps by two bisects instead of scanning the reference.
        #: ⛔ Not an optimisation to taste: a probe now maps to EVERY overlapping transcript, so the
        #: naive scan is O(probes x transcripts) — measured **131 s** for 13,824 probes against 15,669
        #: transcripts, per condition, before this (2026-08-19).
        self._tx_index: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for ref, idxs in self._tx_by_ref.items():
            starts = np.array([transcripts[i].start for i in idxs], dtype=np.int64)
            ends = np.array([transcripts[i].end for i in idxs], dtype=np.int64)
            order = np.argsort(starts, kind="mergesort")
            starts, ends = starts[order], ends[order]
            arr = np.array(idxs, dtype=np.int64)[order]
            # suffix maximum of END: transcripts at or after position j whose end exceeds a query
            suffix_max_end = np.maximum.accumulate(ends[::-1])[::-1]
            self._tx_index[ref] = (starts, ends, arr, suffix_max_end)

        #: ⭐ ONE transcript space (owner, 2026-08-19): a probe is mapped to every transcript whose
        #: exons its GENOMIC blocks overlap — any gene, any isoform, any strand — and to gDNA by the same
        #: overlap. The nascent entities are single-exon transcripts in this list, so a probe on another
        #: gene's exon inside an entity's span reaches the entity exactly as it reaches the gDNA there.
        #: ⛔ There is no per-transcript "nrna" space any more: scoping a pre-mRNA's capture to its own
        #: transcript's probes under-enriched nascent against gDNA by up to 6x (`hop_currency.py`).
        self._mrna_intervals: dict[int, list[WeightedInterval]] = defaultdict(list)
        self._gdna_intervals: dict[str, list[WeightedInterval]] = defaultdict(list)
        self._next_probe_group = 1

        #: Probe layout per space, flattened once. Independent of fragment length, so it is built on
        #: first use and never rebuilt. A one-tuple holding ``None`` means "this space needs the
        #: per-key path" (a probe group with more than one interval); see `_flat_probes`.
        self._flat_probe_cache: dict[str, tuple] = {}
        #: ⭐⭐ **PARTITION MEMO — the capture-aware effective length depends ONLY on the probe panel,
        #: the templates and the fragment width.** It knows nothing about abundances, gDNA rate, strand
        #: specificity or nascent level, so every condition sharing a capture panel computes the SAME
        #: vectors. ⛔ Measured 2026-08-19: 0.45 s per width over 169,399 probe intervals, ~450 widths
        #: ⇒ **260 s per condition**, recomputed identically for each one; shared across conditions the
        #: second and later ones cost **0.2 s**.
        #: ⛔⛔ **IT IS BOUNDED BY CONSTRUCTION, AND THAT IS NOT OPTIONAL** — its ancestor `_mass_cache`
        #: reached **38 GB** by growing per call and never being read back
        #: (`tests/test_sim_capture_partition.py::TestNoUnboundedCache`). This one holds ONE
        #: ``(space, keys, lengths)`` population at a time: a call with a different population CLEARS
        #: it, so the size is at most the number of distinct fragment widths (~450 vectors, ~56 MB at
        #: 15,669 templates) and every entry is one that a later condition reads back.
        self._partition_memo: dict[tuple, np.ndarray] = {}
        self._partition_memo_population: tuple | None = None

    def _transcripts_overlapping(self, ref: str, lo: int, hi: int) -> np.ndarray:
        """Indices of the transcripts on ``ref`` whose span intersects ``[lo, hi)``."""
        idx = self._tx_index.get(ref)
        if idx is None:
            return np.empty(0, dtype=np.int64)
        starts, ends, arr, _suffix = idx
        # candidates start before `hi`; among those keep the ones ending after `lo`
        j = int(np.searchsorted(starts, hi, side="left"))
        if j == 0:
            return np.empty(0, dtype=np.int64)
        keep = ends[:j] > lo
        return arr[:j][keep]

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
            "Hybrid capture enabled: %d transcript targets (incl. nascent entities), %d gDNA refs",
            len(sampler._mrna_intervals),
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

        # ⭐ Prefer the batched path even for one key. `_extra_landscape` materialises and lexsorts a
        # per-start array over EVERY probe on the key, which in the gDNA space means every probe on a
        # whole chromosome: profiled at 6,476,550 `_local_overlap_weights` calls and 110.6 s for one
        # small run, with `_mass_cache` holding those arrays to 37.2 GB. The batched path covers only
        # the merged probe neighbourhoods and caches nothing.
        lengths = np.array([int(seq_len)], dtype=np.int64)
        return float(self.partition_array(space, [key], lengths, frag_len)[0])

    def _flat_probes(self, space: str) -> tuple[np.ndarray, ...] | None:
        """Flatten a space's probe intervals ONCE — the layout does not depend on fragment length.

        Returns ``(keys, start, end, scale, probe_group)`` in ``(key, start)`` order, which is what
        lets `partition_array` find merged probe runs with a single prefix maximum.

        ⚠ ``group`` is carried through rather than collapsed: a group is one physical probe, and a
        probe split across exons must have its pieces **summed** before the max across probes.
        `partition_array` slots on ``(group, piece)`` for exactly that reason.
        """
        cached = self._flat_probe_cache.get(space)
        if cached is not None:
            return cached[0]

        intervals_by_key = {
            "mrna": self._mrna_intervals,
            "gdna": self._gdna_intervals,
        }[space]

        keys: list[int | str] = []
        starts: list[int] = []
        ends: list[int] = []
        scales: list[float] = []
        groups: list[int] = []
        for key, intervals in intervals_by_key.items():
            for interval in sorted(intervals, key=lambda i: (i.start, i.end)):
                keys.append(key)
                starts.append(interval.start)
                ends.append(interval.end)
                scales.append(interval.scale)
                groups.append(interval.probe_group)

        flat = (
            keys,
            np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64),
            np.asarray(scales, dtype=np.float64),
            np.asarray(groups, dtype=np.int64),
        )
        self._flat_probe_cache[space] = (flat,)
        return flat

    def partition_array(
        self,
        space: str,
        keys: Iterable[int | str],
        lengths: np.ndarray,
        frag_len: int,
    ) -> np.ndarray:
        """Vector of capture-aware effective lengths for a pool, computed for every key at once.

        ⭐ **This is the whole cost of a capture-on simulation.** Profiled before this was batched:
        123.6 s of a 127.9 s run, 1,052,172 scalar `partition` calls, and `_mass_cache` grown to 18 GB
        — because the caller loops over every DISTINCT FRAGMENT LENGTH and the old body looped over
        every transcript inside that, doing ten numpy operations on ~300-element arrays each time. The
        arithmetic was never the cost; the per-call overhead was.

        ⚠ **The deleted `ambig_dense_10mb` suite had `frag_std: 0` — ONE distinct length.** A realistic
        distribution draws ~540, so this path is ~540x hotter the moment the suite stops being degenerate
        in exactly the way that made it unable to test anything length-dependent.

        ⚠ **Nothing is cached here, deliberately.** Each ``(key, fragment length)`` pair is visited
        exactly once per call, so the old ``_mass_cache`` writes were pure growth with zero reuse —
        `tests/test_sim_capture_partition.py::TestNoUnboundedCache` pins that.
        """
        keys = list(keys)
        lengths = np.asarray(lengths, dtype=np.int64)
        width = int(frag_len)
        eff = np.maximum(0, lengths - width + 1)
        result = eff.astype(np.float64) * self.config.off_target_weight
        if not self.enabled or self.config.binding_per_base <= 0:
            return result
        population = (space, tuple(keys), lengths.tobytes())
        if population != self._partition_memo_population:
            # a different template population: the stored vectors can never be read again
            self._partition_memo.clear()
            self._partition_memo_population = population
        memo_key = width
        cached = self._partition_memo.get(memo_key)
        if cached is not None:
            return cached.copy()

        flat = self._flat_probes(space)
        flat_keys, starts, ends, scales, probe_groups = flat
        if not flat_keys:
            self._partition_memo[memo_key] = result
            return result.copy()

        key_to_pos = {key: pos for pos, key in enumerate(keys)}
        positions = np.fromiter(
            (key_to_pos.get(key, -1) for key in flat_keys), dtype=np.int64, count=len(flat_keys)
        )
        landscape = self._run_landscape(
            positions, starts, ends, scales, probe_groups, eff, int(frag_len)
        )
        if landscape is None:
            self._partition_memo[memo_key] = result
            return result.copy()
        buffer, run_offset, _run_start, run_first, kept_positions = landscape

        run_sums = np.add.reduceat(buffer, run_offset[:-1])
        sums = np.bincount(kept_positions[run_first], weights=run_sums, minlength=len(keys))
        result = result + self.config.binding_per_base * sums
        self._partition_memo[memo_key] = result
        return result.copy()

    def _run_landscape(self, positions, starts, ends, scales, probe_groups, eff, width):
        """The capture landscape over MERGED PROBE RUNS — the one hot computation, shared by both users.

        ⭐ `partition_array` reduces it to a per-key sum; `_extra_landscape` reads out its nonzero
        positions. Both used to have their own implementation, and the second one — a Python loop over
        every probe on a key plus a lexsort of the concatenation — was **231 s and 38 GB** of a gDNA run
        on its own, because in that space a key is a whole chromosome.

        Returns ``(buffer, run_offset, run_start, run_first, positions)`` or ``None`` if nothing is live.
        The buffer holds ``w(s)`` — the best single probe GROUP's total overlap — over the concatenated
        runs, so run ``r`` covers template positions ``[run_start[r], run_start[r] + run_len[r])``.
        """
        key_eff = np.where(positions >= 0, eff[positions], 0)
        lo = np.maximum(0, starts - width + 1)
        hi = np.minimum(key_eff, ends)
        alive = (positions >= 0) & (hi > lo)
        if not alive.any():
            return None
        positions, starts, ends, scales, probe_groups, lo, hi = (
            positions[alive],
            starts[alive],
            ends[alive],
            scales[alive],
            probe_groups[alive],
            lo[alive],
            hi[alive],
        )
        counts = hi - lo

        # ── merged runs: the buffer covers PROBE NEIGHBOURHOODS, never whole templates ─────────────
        # ⛔ Sizing it by template length is unusable in the gDNA space, where a template is a whole
        # chromosome. Probes arrive grouped by key and sorted by start, so `lo` is non-decreasing within
        # a key and a run ends where the next `lo` clears every `hi` seen so far in that key.
        key_rank = np.concatenate(([0], np.cumsum(positions[1:] != positions[:-1]))).astype(
            np.int64
        )
        stride = int(hi.max()) + 1
        # ⭐ Offsetting by the key's dense rank makes a GLOBAL prefix max come out per-key correct: every
        # entry of key k sits above every entry of key k-1, so it cannot carry across a key boundary.
        running_hi = np.maximum.accumulate(hi + key_rank * stride) - key_rank * stride
        opens_run = np.empty(len(lo), dtype=bool)
        opens_run[0] = True
        opens_run[1:] = (key_rank[1:] != key_rank[:-1]) | (lo[1:] > running_hi[:-1])

        run_id = np.cumsum(opens_run) - 1
        run_first = np.flatnonzero(opens_run)
        run_start = lo[run_first]
        run_end = np.maximum.reduceat(hi, run_first)
        run_offset = np.concatenate(([0], np.cumsum(run_end - run_start)))
        buffer = np.zeros(int(run_offset[-1]), dtype=np.float64)
        base = run_offset[run_id] - run_start[run_id]

        # ── slots: a probe GROUP's rank within its run, and a piece's rank within its group ────────
        # ⛔ Ranking within the KEY does not work in the gDNA space, where one chromosome carries every
        # probe. Within a RUN it is a handful, and two probes sharing a slot are then always in
        # different runs — disjoint buffer ranges, which keeps the scatter duplicate-free.
        order = np.lexsort((probe_groups, run_id))
        ordered_run, ordered_group = run_id[order], probe_groups[order]
        opens_group = np.empty(len(order), dtype=bool)
        opens_group[0] = True
        opens_group[1:] = (ordered_run[1:] != ordered_run[:-1]) | (
            ordered_group[1:] != ordered_group[:-1]
        )
        group_first = np.flatnonzero(opens_group)
        group_sizes = np.diff(np.r_[group_first, len(order)])
        group_run = ordered_run[group_first]
        opens_group_run = np.empty(len(group_first), dtype=bool)
        opens_group_run[0] = True
        opens_group_run[1:] = group_run[1:] != group_run[:-1]
        run_first_group = np.maximum.accumulate(
            np.where(opens_group_run, np.arange(len(group_first)), 0)
        )
        group_slot = np.empty(len(order), dtype=np.int64)
        piece_slot = np.empty(len(order), dtype=np.int64)
        group_slot[order] = np.repeat(np.arange(len(group_first)) - run_first_group, group_sizes)
        piece_slot[order] = np.arange(len(order)) - np.repeat(group_first, group_sizes)

        def scatter(selection):
            """Buffer indices and capture weights for one slot's probes."""
            s_lo, s_counts = lo[selection], counts[selection]
            s_start, s_end, s_scale = starts[selection], ends[selection], scales[selection]
            total = int(s_counts.sum())
            segment = np.repeat(np.arange(len(s_counts)), s_counts)
            within = np.arange(total) - np.repeat(np.cumsum(s_counts) - s_counts, s_counts)
            position = s_lo[segment] + within
            overlap = np.minimum(position + width, s_end[segment]) - np.maximum(
                position, s_start[segment]
            )
            if self.config.min_overlap > 1:
                overlap = np.where(overlap >= self.config.min_overlap, overlap, 0)
            weights = np.maximum(overlap, 0).astype(np.float64) * s_scale[segment]
            return base[selection][segment] + position, weights

        group_buffer = np.zeros_like(buffer)
        n_pieces = int(piece_slot.max()) + 1
        for slot in range(int(group_slot.max()) + 1):
            in_slot = group_slot == slot
            touched = []
            for piece in range(n_pieces):
                selection = in_slot & (piece_slot == piece)
                if not selection.any():
                    continue
                index, weights = scatter(selection)
                # ⚠ A probe GROUP is one physical probe; pieces split across exons must be SUMMED at a
                # shared start before the max across probes. Indices are unique within one
                # (slot, piece), so `+=` accumulates across pieces without `np.add.at`.
                group_buffer[index] += weights
                touched.append(index)
            if not touched:
                continue
            index = np.concatenate(touched)
            buffer[index] = np.maximum(buffer[index], group_buffer[index])
            group_buffer[index] = 0.0

        return buffer, run_offset, run_start, run_first, positions

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
        transcript = self.transcripts[t_idx]
        blocks = transcript_to_genomic_blocks(start, end, transcript)
        self._add_genomic_probe(str(transcript.ref), blocks, probe_group)

    def _add_bed12_probe(
        self,
        ref: str,
        strand: Strand,
        blocks: list[tuple[int, int]],
    ) -> None:
        # ⚠ ``strand`` is parsed and ignored: the library is DNA at capture time, so a probe binds the
        # molecules of either strand that carry its sequence (owner, 2026-08-19).
        del strand
        self._add_genomic_probe(ref, blocks, self._new_probe_group())

    def _add_genomic_probe(
        self,
        ref: str,
        blocks: Sequence[tuple[int, int]],
        probe_group: int,
    ) -> None:
        """One probe, as GENOMIC blocks, mapped by overlap to gDNA and to EVERY transcript on the
        reference whose exons it touches — the owner's ruling, 2026-08-19.

        Per transcript each block is CLIPPED to the exons it overlaps and projected into transcript
        coordinates — overlap, not containment: a fragment carries whatever probe bases it holds, and
        the weight is the fragment's overlap with those, exactly as it is for gDNA at an exon edge.
        The landed pieces keep the full scale when they sit contiguously in transcript coordinates;
        pieces SEPARATED by an intron — a nascent entity's span, like gDNA — take ``gdna_split_penalty``,
        the model's one price for a probe that cannot hybridise as a contiguous whole.
        """
        self._add_gdna_blocks(ref, blocks, self._split_scale(blocks), probe_group)
        lo = min(b[0] for b in blocks)
        hi = max(b[1] for b in blocks)
        for t_idx in self._transcripts_overlapping(ref, lo, hi):
            transcript = self.transcripts[int(t_idx)]
            landed = _clip_blocks_to_transcript(transcript, blocks)
            if not landed:
                continue
            merged = merge_intervals(landed)
            scale = 1.0 if len(merged) == 1 else self.config.gdna_split_penalty
            for start, end in merged:
                self._mrna_intervals[int(t_idx)].append(
                    WeightedInterval(start, end, scale, probe_group)
                )

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
        """Return start positions and best single-probe-GROUP overlap weights, ascending by position.

        ⛔ **This was 231 s and 38 GB of a gDNA run.** It looped in Python over every probe on the key
        and lexsorted the concatenation — and in the gDNA space a key is a whole chromosome carrying
        every probe, so one call did 281 `_local_overlap_weights` calls and a 76 k-element sort **to
        place a mean of about four fragments**. Profiled: 22,838 calls, 6,409,160 inner calls, 21 s of
        `argsort` alone.

        It now reads out `_run_landscape`, the same vectorised computation `partition_array` uses, so
        the two can no longer disagree and the Python loop is gone.

        ⚠ **Nothing is cached.** `sample_starts` asks for each `(key, fragment length)` pair exactly
        once — the caller iterates a counts dict keyed by exactly that — so every `_mass_cache` entry
        was written and never read. That was the 38 GB.
        """
        width = int(frag_len)
        eff_len = int(seq_len) - width + 1
        empty = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64))
        if eff_len <= 0:
            return empty

        intervals = self._get_intervals(space, key)
        if not intervals:
            return empty
        ordered = sorted(intervals, key=lambda i: (i.start, i.end))
        n = len(ordered)
        landscape = self._run_landscape(
            np.zeros(n, dtype=np.int64),
            np.fromiter((i.start for i in ordered), dtype=np.int64, count=n),
            np.fromiter((i.end for i in ordered), dtype=np.int64, count=n),
            np.fromiter((i.scale for i in ordered), dtype=np.float64, count=n),
            np.fromiter((i.probe_group for i in ordered), dtype=np.int64, count=n),
            np.array([eff_len], dtype=np.int64),
            width,
        )
        if landscape is None:
            return empty
        buffer, run_offset, run_start, _run_first, _positions = landscape

        # Runs are disjoint and ascending within a key, so concatenating their positions is already
        # sorted — which is the contract `sample_starts` relies on.
        run_len = np.diff(run_offset)
        position = np.repeat(run_start, run_len) + (
            np.arange(len(buffer)) - np.repeat(run_offset[:-1], run_len)
        )
        nonzero = buffer > 0
        return position[nonzero], buffer[nonzero]

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


def _clip_blocks_to_transcript(
    transcript: Transcript, blocks: Sequence[tuple[int, int]]
) -> list[tuple[int, int]]:
    """The parts of genomic ``blocks`` that fall inside ``transcript``'s exons, in transcript
    coordinates (strand-oriented). ⚠ Unlike `intervals.project_genomic_block_to_transcript` this does
    NOT require a block to be fully exonic — capture binds by overlap, so a probe hanging off an exon
    edge still binds the bases that are there."""
    tx_len = int(transcript.length or transcript.compute_length())
    out: list[tuple[int, int]] = []
    for block_start, block_end in blocks:
        consumed = 0
        for exon in transcript.exons:
            lo = max(block_start, exon.start)
            hi = min(block_end, exon.end)
            if lo < hi:
                tx_start = consumed + (lo - exon.start)
                tx_end = consumed + (hi - exon.start)
                if transcript.strand == Strand.NEG:
                    tx_start, tx_end = tx_len - tx_end, tx_len - tx_start
                out.append((tx_start, tx_end))
            consumed += exon.end - exon.start
    return out
