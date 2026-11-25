"""Hybrid-capture runtime: probe loading, and the questions the read engine asks of a panel.

``CaptureSampler`` holds one probe panel as weighted intervals in two coordinate spaces, ``"transcript"`` keyed by
transcript row (nascent entities included) and ``"gdna"`` keyed by reference. For a template and a fragment width it
gives the capture-aware effective length (:meth:`~CaptureSampler.partition_array`), draws fragment starts
(:meth:`~CaptureSampler.sample_starts`) and weighs one fragment (:meth:`~CaptureSampler.fragment_weight`). Without a
probe file the sampler is disabled: uniform starts and plain effective lengths.

A fragment weighs ``off_target_weight + binding_per_base * max(strength * overlap)``, the maximum over the contiguous
probe parts in its own template, ``overlap`` the fragment's overlap with a part and ``strength`` the part's probe's.
A probe binds at a strength its own sequence gives, :func:`gc_strength` of its GC content, 1 for the best-binding
probe; it is read once from the reference when the panel loads. A molecule hybridises through one contiguous stretch,
so a probe split across a splice junction is whole in a transcript holding the junction and two parts elsewhere
(gDNA, a nascent span, an isoform without the junction), where a fragment binds the better part, never their sum;
the split probe is one probe, its sequence its blocks joined, and its parts share its strength. Each probe maps by
genomic overlap, not containment, onto gDNA and onto every transcript whose exons it touches, on either strand,
since the library is DNA at capture.

Bases the reference does not know (``N``) are masked per template: a fragment overlapping one is never drawn and adds
no opportunity. The effective length and the start distribution read one vectorised computation over merged probe
runs (:meth:`~CaptureSampler._run_landscape`), so they cannot disagree, and its buffer is sized by probe
neighbourhoods, not template length, which keeps whole-chromosome gDNA templates tractable.
"""

from __future__ import annotations

import gzip
import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from ...observable import true_runs, unknown_bases
from ...transcript import Transcript
from ...types import Strand
from ..bam import transcript_to_genomic_blocks
from ..intervals import clip_blocks_to_transcript, merge_intervals
from .config import CaptureConfig

logger = logging.getLogger(__name__)

__all__ = ["CaptureSampler", "ProbeInterval", "gc_strength"]


def gc_strength(gc, config: CaptureConfig):
    """A probe's binding strength relative to the best-binding probe, from the GC fraction of its own sequence:
    ``2 ** -((gc - gc_optimum) / gc_half_distance) ** 2``.

    Capture falls away at both GC extremes, so the simplest form right in shape is one peak at ``gc_optimum`` that
    falls symmetrically on a log scale, halving at ``gc_half_distance`` on either side, a quarter-strength probe
    lying ``sqrt(2)`` half distances away. Scalar in, float out; an array in, an array out. ``gc_half_distance =
    inf`` gives exactly 1 everywhere: every probe binds alike.
    """
    out = np.exp2(
        -np.square((np.asarray(gc, dtype=np.float64) - config.gc_optimum) / config.gc_half_distance)
    )
    return float(out) if out.ndim == 0 else out


def _merged(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Half-open intervals ascending by start, overlapping ones merged: ``int64[n, 2]``."""
    if lo.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    run_hi = np.maximum.accumulate(hi)
    opens = np.empty(lo.size, dtype=bool)
    opens[0] = True
    opens[1:] = lo[1:] > run_hi[:-1]
    first = np.flatnonzero(opens)
    ends = np.maximum.reduceat(hi, first)
    return np.stack([lo[first], ends], axis=1)


@dataclass(frozen=True, slots=True)
class ProbeInterval:
    """One contiguous part of a probe on one coordinate axis — the unit a fragment binds through — with its probe's
    binding strength (:func:`gc_strength`), which every part of the probe shares."""

    start: int
    end: int
    strength: float


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
        self._tx_id_to_index = {
            str(t.t_id): idx for idx, t in enumerate(transcripts) if t.t_id is not None
        }
        # a duplicate the index folded into a kept transcript (`--collapse-duplicate-transcripts`) is that
        # transcript: a probe naming it binds the kept one, the same molecule
        for idx, t in enumerate(transcripts):
            for alias in t.aliases:
                self._tx_id_to_index.setdefault(str(alias), idx)
        self._tx_by_ref: dict[str, list[int]] = defaultdict(list)
        for idx, t in enumerate(transcripts):
            if t.ref is not None:
                self._tx_by_ref[str(t.ref)].append(idx)
        #: Per reference, the transcripts' starts, ends and indices, sorted by start.
        self._tx_index: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for ref, idxs in self._tx_by_ref.items():
            starts = np.array([transcripts[i].start for i in idxs], dtype=np.int64)
            ends = np.array([transcripts[i].end for i in idxs], dtype=np.int64)
            order = np.argsort(starts, kind="mergesort")
            starts, ends = starts[order], ends[order]
            arr = np.array(idxs, dtype=np.int64)[order]
            self._tx_index[ref] = (starts, ends, arr)

        #: One transcript space: a probe maps to every transcript whose exons it overlaps, and to gDNA by the same
        #: overlap. A nascent entity is a single-exon row here, so a probe on another gene's exon inside its span
        #: reaches it as it reaches the gDNA there. Nascent capture must not be scoped to its own transcript's
        #: probes, or nascent RNA would be under-enriched against the gDNA beside it.
        self._transcript_intervals: dict[int, list[ProbeInterval]] = defaultdict(list)
        self._gdna_intervals: dict[str, list[ProbeInterval]] = defaultdict(list)
        #: Each loaded probe's binding strength, in the probe file's order.
        self._probe_strengths: list[float] = []

        #: Probe layout per space, flattened on first use (`_flat_probes`); it does not depend on fragment length.
        self._flat_probe_cache: dict[str, tuple] = {}
        #: Partition memo, keyed by fragment width. The effective length depends only on the panel, the templates and
        #: the width, so every condition sharing a panel reads the same vectors. It must stay bounded: it holds one
        #: ``(space, keys, lengths)`` population and clears when the population changes, so it holds at most one
        #: entry per distinct width.
        self._partition_memo: dict[tuple, np.ndarray] = {}
        self._partition_memo_population: tuple | None = None
        #: Per template, the runs of bases the reference does not know: ``{space: {key: int64[n, 2]}}``, sorted,
        #: disjoint ``[start, end)`` in the template's own coordinates. A fragment overlapping one is never drawn.
        self._unknown: dict[str, dict[int | str, np.ndarray]] = {}

    # -- Unknown sequence --------------------------------------------------

    @staticmethod
    def _norm(space: str, key: int | str) -> int | str:
        return int(key) if space == "transcript" else str(key)

    def mask_unknown(self, space: str, key: int | str, runs) -> None:
        """Install one template's unknown runs, ``[start, end)`` pairs in any order, replacing any before; they are
        kept sorted and merged. A changed mask clears the partition memo; an unchanged one keeps it."""
        runs = np.asarray(runs, dtype=np.int64).reshape(-1, 2)
        if len(runs):
            if (runs[:, 1] <= runs[:, 0]).any() or (runs[:, 0] < 0).any():
                raise ValueError(
                    f"an unknown run is a non-empty [start, end) at or after 0: {runs.tolist()}"
                )
            order = np.argsort(runs[:, 0], kind="stable")
            runs = _merged(runs[order, 0], runs[order, 1])
        masks = self._unknown.setdefault(space, {})
        k = self._norm(space, key)
        old = masks.get(k)
        if (old is None and len(runs) == 0) or (old is not None and np.array_equal(old, runs)):
            return
        if len(runs):
            masks[k] = runs
        else:
            del masks[k]
        self._partition_memo.clear()
        self._partition_memo_population = None

    def unknown(self, space: str, key: int | str) -> np.ndarray:
        """One template's unknown runs, ``int64[n, 2]`` (empty when it has none)."""
        runs = self._unknown.get(space, {}).get(self._norm(space, key))
        return np.empty((0, 2), dtype=np.int64) if runs is None else runs

    def mask_unknown_sequence(self, fasta, gdna_refs) -> None:
        """Mask every template from ``fasta`` (a ``pysam.FastaFile``): each reference in ``gdna_refs`` at its own
        ``N`` runs, and each transcript row at the runs of its spliced template, read 5′ to 3′."""
        gdna_refs = {str(r) for r in gdna_refs}
        for ref in fasta.references:
            rows = self._tx_by_ref.get(ref, ())
            if ref not in gdna_refs and not rows:
                continue
            isn = unknown_bases(fasta.fetch(ref))
            runs = true_runs(isn)
            if ref in gdna_refs:
                self.mask_unknown("gdna", ref, runs)
            n_before = np.concatenate(([0], np.cumsum(isn)))
            for t_idx in rows:
                t = self.transcripts[t_idx]
                if n_before[t.end] == n_before[t.start]:
                    self.mask_unknown("transcript", t_idx, ())
                    continue
                template = np.concatenate([isn[e.start : e.end] for e in t.exons])
                if t.strand == Strand.NEG:
                    template = template[::-1]
                self.mask_unknown("transcript", t_idx, true_runs(template))

    def _masked_among(self, space: str, keys) -> list[tuple[int, int | str]]:
        """``(position, key)`` of each template in ``keys`` that has unknown runs."""
        masks = self._unknown.get(space)
        if not masks:
            return []
        if space == "transcript":
            karr = np.asarray(keys, dtype=np.int64)
            hit = np.flatnonzero(
                np.isin(karr, np.fromiter(masks, dtype=np.int64, count=len(masks)))
            )
            return [(int(i), int(karr[i])) for i in hit]
        return [(i, str(k)) for i, k in enumerate(keys) if str(k) in masks]

    def blocked_starts(
        self, space: str, key: int | str, eff_len: int, width: int
    ) -> np.ndarray | None:
        """The starts in ``[0, eff_len)`` whose fragment of ``width`` overlaps an unknown base, as merged
        ``[start, end)`` runs; ``None`` when there are none."""
        runs = self._unknown.get(space, {}).get(self._norm(space, key))
        if runs is None or eff_len <= 0:
            return None
        lo = np.maximum(runs[:, 0] - int(width) + 1, 0)
        hi = np.minimum(runs[:, 1], eff_len)
        keep = hi > lo
        if not keep.any():
            return None
        return _merged(lo[keep], hi[keep])

    def admissible_counts(self, space: str, keys, lengths: np.ndarray, frag_len: int) -> np.ndarray:
        """Per template, the number of starts whose fragment of ``frag_len`` covers no unknown base."""
        lengths = np.asarray(lengths, dtype=np.int64)
        width = int(frag_len)
        out = np.maximum(0, lengths - width + 1)
        for i, key in self._masked_among(space, keys):
            blocked = self.blocked_starts(space, key, int(out[i]), width)
            if blocked is not None:
                out[i] -= int((blocked[:, 1] - blocked[:, 0]).sum())
        return out

    def _zero_blocked(
        self, buffer, run_offset, run_start, run_keys, eff, keys, space, width
    ) -> None:
        """Zero the landscape at every start whose fragment overlaps an unknown base. Runs are disjoint and ascending
        within a key, so each blocked stretch meets a contiguous slice of the key's runs."""
        run_end = run_start + np.diff(run_offset)
        for kp, key in self._masked_among(space, keys):
            sel = np.flatnonzero(run_keys == kp)
            if not sel.size:
                continue
            blocked = self.blocked_starts(space, key, int(eff[kp]), width)
            if blocked is None:
                continue
            rs, re, ro = run_start[sel], run_end[sel], run_offset[sel]
            for b_lo, b_hi in blocked.tolist():
                for i in range(
                    int(np.searchsorted(re, b_lo, side="right")), int(np.searchsorted(rs, b_hi))
                ):
                    lo, hi = max(int(rs[i]), b_lo), min(int(re[i]), b_hi)
                    if hi > lo:
                        buffer[int(ro[i]) + lo - int(rs[i]) : int(ro[i]) + hi - int(rs[i])] = 0.0

    def _uniform_starts(
        self, space: str, key: int | str, eff_len: int, frag_len: int, size: int, rng
    ) -> np.ndarray:
        """Uniform draws over the admissible starts: the k-th admissible start, k uniform."""
        blocked = self.blocked_starts(space, key, eff_len, frag_len)
        if blocked is None:
            return rng.integers(0, eff_len, size=size, dtype=np.int64)
        free = np.concatenate(([0], blocked.ravel(), [eff_len])).reshape(-1, 2)
        free = free[free[:, 1] > free[:, 0]]
        before = np.concatenate(([0], np.cumsum(free[:, 1] - free[:, 0])))
        k = rng.integers(0, int(before[-1]), size=size, dtype=np.int64)
        seg = np.searchsorted(before, k, side="right") - 1
        return free[seg, 0] + (k - before[seg])

    def _transcripts_overlapping(self, ref: str, lo: int, hi: int) -> np.ndarray:
        """Indices of the transcripts on ``ref`` whose span intersects ``[lo, hi)``."""
        idx = self._tx_index.get(ref)
        if idx is None:
            return np.empty(0, dtype=np.int64)
        starts, ends, arr = idx
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
        fasta=None,
    ) -> "CaptureSampler":
        """Build a sampler from a probe file, or return disabled if absent. ``fasta`` (a ``pysam.FastaFile``) is the
        reference each probe's sequence, and so its binding strength, is read from; it is required unless the
        config has no GC effect (``gc_half_distance = inf``)."""
        if config is None or not config.probes:
            return cls.disabled(transcripts, ref_lengths)
        _validate_config(config)
        if fasta is None and not math.isinf(config.gc_half_distance):
            raise ValueError(
                "a probe binds at the strength its own GC content gives, so the capture sampler reads each probe's "
                "sequence: pass the genome the panel lies on (`fasta`)"
            )
        sampler = cls(config, transcripts, ref_lengths, enabled=True)
        sampler._load_probe_file(Path(config.probes), fasta)
        strengths = sampler.probe_strengths()
        logger.info(
            "Hybrid capture enabled: %d transcript targets (incl. nascent entities), %d gDNA refs; %d probes, "
            "binding strength 5/50/95 %%: %s",
            len(sampler._transcript_intervals),
            len(sampler._gdna_intervals),
            len(strengths),
            " / ".join(f"{q:.3f}" for q in np.percentile(strengths, [5, 50, 95]))
            if len(strengths)
            else "-",
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
        """The capture-aware effective length for one template: :meth:`partition_array` for one key."""
        if int(seq_len) - int(frag_len) + 1 <= 0:
            return 0.0
        lengths = np.array([int(seq_len)], dtype=np.int64)
        return float(self.partition_array(space, [key], lengths, frag_len)[0])

    def _flat_probes(self, space: str) -> tuple[np.ndarray, ...] | None:
        """A space's probe parts as ``(keys, start, end, strength)`` in ``(key, start)`` order, which lets
        `_run_landscape` find merged runs with one prefix maximum. A fragment binds the part with the largest
        strength times overlap, so a part carries its probe's strength and no other identity. Built once per space.
        """
        cached = self._flat_probe_cache.get(space)
        if cached is not None:
            return cached[0]

        intervals_by_key = {
            "transcript": self._transcript_intervals,
            "gdna": self._gdna_intervals,
        }[space]

        keys: list[int | str] = []
        starts: list[int] = []
        ends: list[int] = []
        strengths: list[float] = []
        for key, intervals in intervals_by_key.items():
            for interval in sorted(intervals, key=lambda i: (i.start, i.end)):
                keys.append(key)
                starts.append(interval.start)
                ends.append(interval.end)
                strengths.append(interval.strength)

        flat = (
            keys,
            np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64),
            np.asarray(strengths, dtype=np.float64),
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
        """Capture-aware effective lengths of a pool of templates at one fragment width, every key at once.

        The caller loops over every distinct width, so the work is batched per width, not per template. Results are
        memoised per width within one template population (``_partition_memo``).
        """
        keys = list(keys)
        lengths = np.asarray(lengths, dtype=np.int64)
        width = int(frag_len)
        eff = np.maximum(0, lengths - width + 1)
        result = (
            self.admissible_counts(space, keys, lengths, width).astype(np.float64)
            * self.config.off_target_weight
        )
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
        flat_keys, starts, ends, strengths = flat
        if not flat_keys:
            self._partition_memo[memo_key] = result
            return result.copy()

        key_to_pos = {key: pos for pos, key in enumerate(keys)}
        positions = np.fromiter(
            (key_to_pos.get(key, -1) for key in flat_keys), dtype=np.int64, count=len(flat_keys)
        )
        landscape = self._run_landscape(positions, starts, ends, strengths, eff, int(frag_len))
        if landscape is None:
            self._partition_memo[memo_key] = result
            return result.copy()
        buffer, run_offset, run_start, run_first, kept_positions = landscape
        self._zero_blocked(
            buffer, run_offset, run_start, kept_positions[run_first], eff, keys, space, width
        )

        run_sums = np.add.reduceat(buffer, run_offset[:-1])
        sums = np.bincount(kept_positions[run_first], weights=run_sums, minlength=len(keys))
        result = result + self.config.binding_per_base * sums
        self._partition_memo[memo_key] = result
        return result.copy()

    def _run_landscape(self, positions, starts, ends, strengths, eff, width):
        """The capture landscape over merged probe runs: `partition_array` sums it per key and `binding_landscape` reads
        its nonzero positions, so the effective length and the start distribution cannot disagree.

        Returns ``(buffer, run_offset, run_start, run_first, positions)``, or ``None`` if nothing is live. The buffer
        holds ``w(s)``, the largest strength times overlap with one contiguous probe part, over the concatenated runs:
        run ``r`` covers template positions ``[run_start[r], run_start[r] + run_len[r])``.
        """
        key_eff = np.where(positions >= 0, eff[positions], 0)
        lo = np.maximum(0, starts - width + 1)
        hi = np.minimum(key_eff, ends)
        alive = (positions >= 0) & (hi > lo)
        if not alive.any():
            return None
        positions, starts, ends, strengths, lo, hi = (
            positions[alive],
            starts[alive],
            ends[alive],
            strengths[alive],
            lo[alive],
            hi[alive],
        )
        counts = hi - lo

        # Merged runs: the buffer covers probe neighbourhoods, never whole templates. Probes arrive grouped by key and
        # sorted by start, so a run ends where the next `lo` clears every `hi` seen so far in that key.
        key_rank = np.concatenate(([0], np.cumsum(positions[1:] != positions[:-1]))).astype(
            np.int64
        )
        stride = int(hi.max()) + 1
        # Offsetting by the key's dense rank keeps a global prefix max from carrying across a key boundary.
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

        # Slots: a part's rank within its run. Two parts sharing a slot are always in different runs, whose buffer
        # ranges are disjoint, so each slot's scatter is duplicate-free.
        slot = np.arange(len(lo)) - run_first[run_id]

        def scatter(selection):
            """Buffer indices and capture weights for one slot's parts."""
            s_lo, s_counts = lo[selection], counts[selection]
            s_start, s_end = starts[selection], ends[selection]
            total = int(s_counts.sum())
            segment = np.repeat(np.arange(len(s_counts)), s_counts)
            within = np.arange(total) - np.repeat(np.cumsum(s_counts) - s_counts, s_counts)
            position = s_lo[segment] + within
            overlap = np.minimum(position + width, s_end[segment]) - np.maximum(
                position, s_start[segment]
            )
            if self.config.min_overlap > 1:
                overlap = np.where(overlap >= self.config.min_overlap, overlap, 0)
            weights = np.maximum(overlap, 0) * strengths[selection][segment]
            return base[selection][segment] + position, weights

        # A fragment binds through one part, the one with the largest strength times overlap: the max across parts,
        # never their sum.
        for rank in range(int(slot.max()) + 1):
            index, weights = scatter(slot == rank)
            buffer[index] = np.maximum(buffer[index], weights)

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

        n_free = int(self.admissible_counts(space, [key], np.array([int(seq_len)]), frag_len)[0])
        if n_free <= 0:
            raise ValueError(
                f"{count} fragments of {frag_len} bp asked of {space} template {key!r}, where every start covers an "
                "unknown base: its opportunity is 0, so it can be allocated none"
            )
        if not self.enabled or self.config.binding_per_base <= 0:
            return self._uniform_starts(space, key, eff_len, frag_len, count, rng)

        extra_starts, raw_weights = self.binding_landscape(space, key, seq_len, frag_len)
        raw_total = float(raw_weights.sum())
        baseline_mass = self.config.off_target_weight * n_free
        extra_mass = self.config.binding_per_base * raw_total
        total_mass = baseline_mass + extra_mass
        if raw_total <= 0 or total_mass <= 0:
            return self._uniform_starts(space, key, eff_len, frag_len, count, rng)

        starts = np.empty(count, dtype=np.int64)
        baseline_mask = rng.random(count) < (baseline_mass / total_mass)
        n_baseline = int(baseline_mask.sum())
        if n_baseline > 0:
            starts[baseline_mask] = self._uniform_starts(
                space, key, eff_len, frag_len, n_baseline, rng
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
        """The unnormalised sampling weight of one concrete fragment: the weight model written fragment by fragment,
        the reference the vectorised landscape (:meth:`partition_array`, :meth:`sample_starts`) is checked against.
        The simulator itself never calls it."""
        eff_len = int(seq_len) - int(frag_len) + 1
        if eff_len <= 0 or frag_start < 0 or frag_start >= eff_len:
            return 0.0
        blocked = self.blocked_starts(space, key, eff_len, frag_len)
        if (
            blocked is not None
            and ((blocked[:, 0] <= frag_start) & (frag_start < blocked[:, 1])).any()
        ):
            return 0.0
        if not self.enabled or self.config.binding_per_base <= 0:
            return self.config.off_target_weight
        frag_end = frag_start + frag_len
        score = 0.0
        for interval in self.probe_intervals(space, key):
            overlap = min(frag_end, interval.end) - max(frag_start, interval.start)
            if overlap >= self.config.min_overlap:
                score = max(score, interval.strength * overlap)
        return self.config.off_target_weight + self.config.binding_per_base * score

    # -- Probe loading -----------------------------------------------------

    def _load_probe_file(self, path: Path, fasta) -> None:
        lines = _read_probe_lines(path)
        if not lines:
            raise ValueError(f"capture probe file is empty: {path}")
        probe_format = _detect_probe_format(lines, self.config.probe_format)
        if probe_format == "transcript":
            self._load_transcript_probes(lines, path, fasta)
        elif probe_format == "bed12":
            self._load_bed12_probes(lines, path, fasta)
        else:
            raise ValueError(f"unsupported capture probe format: {probe_format!r}")

    def _load_transcript_probes(self, lines: list[str], path: Path, fasta) -> None:
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
                fasta,
            )
            n_loaded += 1
        logger.info("Loaded %d transcript-coordinate capture probes from %s", n_loaded, path)

    def _load_bed12_probes(self, lines: list[str], path: Path, fasta) -> None:
        n_loaded = 0
        for row_number, line in enumerate(lines, start=1):
            fields = _split_fields(line)
            if len(fields) < 12:
                raise ValueError(f"{path}:{row_number}: expected 12 BED fields")
            if not _is_int(fields[1]) or not _is_int(fields[2]):
                if row_number == 1:  # a header, as the transcript format allows on its first row
                    continue
                raise ValueError(f"{path}:{row_number}: chromStart and chromEnd must be integers")
            ref = fields[0]
            chrom_start = int(fields[1])
            block_count = int(fields[9])
            block_sizes = _parse_bed_list(fields[10])
            block_starts = _parse_bed_list(fields[11])
            if block_count != len(block_sizes) or block_count != len(block_starts):
                raise ValueError(f"{path}:{row_number}: BED12 blockCount mismatch")
            if min(block_sizes) < 1:
                raise ValueError(f"{path}:{row_number}: a BED12 block holds no base")
            blocks = [
                (chrom_start + rel_start, chrom_start + rel_start + size)
                for size, rel_start in zip(block_sizes, block_starts)
            ]
            self._add_genomic_probe(
                ref, blocks, self._probe_strength(ref, blocks, fasta, f"{path}:{row_number}")
            )
            n_loaded += 1
        logger.info("Loaded %d BED12 capture probes from %s", n_loaded, path)

    def _add_transcript_probe(
        self,
        transcript_id: str,
        start: int,
        end: int,
        path: Path,
        row_number: int,
        fasta,
    ) -> None:
        t_idx = self._tx_id_to_index.get(transcript_id)
        if t_idx is None:
            raise ValueError(
                f"{path}:{row_number}: a probe on {transcript_id}, which is neither a transcript of the index nor a "
                "duplicate it folded into one: the panel was designed on another annotation"
            )
        tx_len = int(self._tx_lengths[t_idx])
        if start < 0 or end <= start or end > tx_len:
            raise ValueError(
                f"{path}:{row_number}: invalid probe interval for {transcript_id}: "
                f"[{start}, {end}) outside transcript length {tx_len}"
            )

        transcript = self.transcripts[t_idx]
        ref = str(transcript.ref)
        blocks = transcript_to_genomic_blocks(start, end, transcript)
        self._add_genomic_probe(
            ref, blocks, self._probe_strength(ref, blocks, fasta, f"{path}:{row_number}")
        )

    def _probe_strength(
        self, ref: str, blocks: Sequence[tuple[int, int]], fasta, where: str
    ) -> float:
        """One probe's binding strength, :func:`gc_strength` of the GC fraction of its sequence, its genomic blocks
        joined, over the bases the reference knows (``N`` is unknown, as for the mask). Without a reference there is
        no GC effect (``from_config`` refuses otherwise), and every probe binds at 1. A probe that lies wholly on
        unknown bases is refused: no fragment over it can ever be drawn, so the panel was designed on another
        reference."""
        if fasta is None:
            strength = 1.0
        else:
            if ref not in fasta.references:
                raise ValueError(f"{where}: a probe on {ref}, which the reference does not hold")
            sequence = "".join(fasta.fetch(ref, int(lo), int(hi)) for lo, hi in blocks).upper()
            known = len(sequence) - sequence.count("N")
            if known == 0:
                raise ValueError(
                    f"{where}: a probe on {ref} {list(blocks)} lies wholly on unknown bases (N), where no fragment "
                    "can be drawn: the panel was designed on another reference"
                )
            strength = gc_strength((sequence.count("G") + sequence.count("C")) / known, self.config)
        self._probe_strengths.append(strength)
        return strength

    def _add_genomic_probe(
        self, ref: str, blocks: Sequence[tuple[int, int]], strength: float
    ) -> None:
        """Map one probe, as genomic blocks, by overlap to gDNA and to every transcript whose exons it touches, every
        part at the probe's one ``strength``.

        Each template receives the probe's contiguous parts in its own coordinates: a transcript holding the
        junction a probe spans receives it whole, since its landed pieces merge, while gDNA, a nascent span and an
        isoform without the junction receive the separate parts.
        """
        self._add_gdna_blocks(ref, merge_intervals(blocks), strength)
        lo = min(b[0] for b in blocks)
        hi = max(b[1] for b in blocks)
        for t_idx in self._transcripts_overlapping(ref, lo, hi):
            transcript = self.transcripts[int(t_idx)]
            landed = clip_blocks_to_transcript(transcript, blocks)
            if not landed:
                continue
            for start, end in merge_intervals(landed):
                self._transcript_intervals[int(t_idx)].append(ProbeInterval(start, end, strength))

    def _add_gdna_blocks(
        self, ref: str, blocks: Sequence[tuple[int, int]], strength: float
    ) -> None:
        ref_len = self.ref_lengths.get(ref)
        if ref_len is None:
            return
        for start, end in blocks:
            interval = _clip_interval(start, end, ref_len, strength)
            if interval is not None:
                self._gdna_intervals[ref].append(interval)

    def probe_strengths(self) -> np.ndarray:
        """Each loaded probe's binding strength, in the probe file's order."""
        return np.asarray(self._probe_strengths, dtype=np.float64)

    def probed_bases(self, space: str, keys: Iterable[int | str]) -> np.ndarray:
        """Per template, how many of its bases lie under a probe part and are known to the reference: the panel's
        geometry, blind to how strongly each probe binds."""
        keys = list(keys)
        out = np.zeros(len(keys), dtype=np.int64)
        if not self.enabled:
            return out
        for pos, key in enumerate(keys):
            parts = self.probe_intervals(space, key)
            if not parts:
                continue
            lo = np.fromiter((p.start for p in parts), dtype=np.int64, count=len(parts))
            hi = np.fromiter((p.end for p in parts), dtype=np.int64, count=len(parts))
            order = np.argsort(lo, kind="stable")
            covered = _merged(lo[order], hi[order])
            n = int((covered[:, 1] - covered[:, 0]).sum())
            for u_lo, u_hi in self.unknown(space, key).tolist():
                n -= int(
                    np.maximum(
                        np.minimum(covered[:, 1], u_hi) - np.maximum(covered[:, 0], u_lo), 0
                    ).sum()
                )
            out[pos] = n
        return out

    def probe_intervals(self, space: str, key: int | str) -> list[ProbeInterval]:
        """The probe parts mapped onto one template: a transcript row in ``"transcript"``, a reference in ``"gdna"``."""
        if space == "transcript":
            return self._transcript_intervals.get(int(key), [])
        if space == "gdna":
            return self._gdna_intervals.get(str(key), [])
        raise ValueError(f"unknown capture coordinate space: {space!r}")

    def binding_landscape(
        self,
        space: str,
        key: int | str,
        seq_len: int,
        frag_len: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Start positions and their best single-part weights (strength times overlap), ascending by position, read
        out of `_run_landscape`. Not cached: `sample_starts` asks for each ``(key, fragment length)`` pair once.
        """
        width = int(frag_len)
        eff_len = int(seq_len) - width + 1
        empty = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64))
        if eff_len <= 0:
            return empty

        intervals = self.probe_intervals(space, key)
        if not intervals:
            return empty
        ordered = sorted(intervals, key=lambda i: (i.start, i.end))
        n = len(ordered)
        landscape = self._run_landscape(
            np.zeros(n, dtype=np.int64),
            np.fromiter((i.start for i in ordered), dtype=np.int64, count=n),
            np.fromiter((i.end for i in ordered), dtype=np.int64, count=n),
            np.fromiter((i.strength for i in ordered), dtype=np.float64, count=n),
            np.array([eff_len], dtype=np.int64),
            width,
        )
        if landscape is None:
            return empty
        buffer, run_offset, run_start, _run_first, _positions = landscape
        self._zero_blocked(
            buffer,
            run_offset,
            run_start,
            np.zeros(len(run_start), dtype=np.int64),
            np.array([eff_len]),
            [key],
            space,
            width,
        )

        # Runs are disjoint and ascending, so the concatenated positions are sorted.
        run_len = np.diff(run_offset)
        position = np.repeat(run_start, run_len) + (
            np.arange(len(buffer)) - np.repeat(run_offset[:-1], run_len)
        )
        nonzero = buffer > 0
        return position[nonzero], buffer[nonzero]


def _validate_config(config: CaptureConfig) -> None:
    if config.off_target_weight < 0:
        raise ValueError("capture.off_target_weight must be >= 0")
    if config.binding_per_base < 0:
        raise ValueError("capture.binding_per_base must be >= 0")
    if config.min_overlap < 1:
        raise ValueError("capture.min_overlap must be >= 1")
    if not 0.0 <= config.gc_optimum <= 1.0:
        raise ValueError("capture.gc_optimum is a GC fraction, in [0, 1]")
    if not config.gc_half_distance > 0.0:
        raise ValueError("capture.gc_half_distance must be > 0 (inf: no GC effect)")


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
            "capture.probe_format must be one of auto, transcript, transcript_tsv (or tsv), bed12 (or bed)"
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


def _clip_interval(start: int, end: int, length: int, strength: float) -> ProbeInterval | None:
    start = max(0, int(start))
    end = min(int(length), int(end))
    if end <= start:
        return None
    return ProbeInterval(start, end, strength)
