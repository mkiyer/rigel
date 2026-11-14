"""rigel.locus — the EM's loci: connected components of transcripts linked by shared fragments.

* :class:`MultiLocus` — one connected component; the unit the EM is run on.
* :class:`Locus` — one contiguous interval of a multi-locus's merged transcript footprint.
  ``rigel.calibration.priors`` overlaps these intervals with the calibration regions to project
  calibration's per-region numbers onto the EM loci.
* :func:`build_multi_loci` — the connected-component partition, one :class:`MultiLocus` per
  component.

The EM's two pseudocounts are formed per locus by ``rigel.pipeline.em_pseudocounts``, from
calibration's two per-locus counts, its gDNA fragments and every fragment (``rigel.calibration.priors``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from .index import TranscriptIndex
from .native import connected_components as _cc_native
from .scored_fragments import ScoredFragments


# ---------------------------------------------------------------------------
# Locus / MultiLocus dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Locus:
    """One contiguous interval of a :class:`MultiLocus`'s merged transcript footprint.

    Most multi-loci have exactly one. A component whose transcripts do not overlap (linked
    through multimapping fragments, say, or spanning several references) carries several.
    """

    ref: str
    ref_id: int
    start: int
    end: int

    @property
    def span(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, slots=True)
class MultiLocus:
    """A connected component of transcripts linked by shared fragments.

    The unit the EM is run on. Its merged transcript footprint is one or more
    :class:`Locus` intervals.

    Attributes
    ----------
    multi_locus_id : int
        Sequential label (0-based).
    transcript_indices : np.ndarray
        int32 — global transcript indices in this multi-locus.
    unit_indices : np.ndarray
        int32 — EM unit indices (rows in the global CSR) in this multi-locus, in BAM order.
    span_bp : int
        Total merged genomic footprint (bp), ``sum(l.span for l in loci)``: the loci table's ``locus_span_bp``.
    loci : tuple[Locus, ...]
        The contiguous intervals, sorted by the ``ref`` column's categorical code (not ``ref_id``), then ``start``.
    """

    multi_locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray
    span_bp: int
    loci: tuple[Locus, ...]


# ---------------------------------------------------------------------------
# MultiLocus builder: C++ union-find connected components
# ---------------------------------------------------------------------------


def _in_bam_order(units: np.ndarray, frag_ids: np.ndarray) -> np.ndarray:
    """``units`` reordered by their fragments' ``frag_id``, which is BAM order; always a fresh array."""
    return units[np.argsort(frag_ids[units], kind="stable")]


class _TranscriptSpans(NamedTuple):
    """Every transcript's span, for merging a component's footprint. ``ref_codes`` are the ``ref`` column's
    categorical codes, which drive the merge's sort and compare; they are not the canonical ref-id space
    (``index.ref_name_to_id``), which ``Locus.ref_id`` must carry because ``calibration.priors`` bins locus blocks by
    ref id against ``RegionArrays``, so ``ref_ids`` maps a code to it."""

    ref_codes: np.ndarray
    starts: np.ndarray
    ends: np.ndarray
    ref_names: np.ndarray
    ref_ids: np.ndarray

    @classmethod
    def from_index(cls, index: TranscriptIndex) -> "_TranscriptSpans":
        ref_cat = index.t_df["ref"].cat
        ref_names = ref_cat.categories.values
        return cls(
            ref_codes=ref_cat.codes.values,
            starts=index.t_df["start"].values,
            ends=index.t_df["end"].values,
            ref_names=ref_names,
            ref_ids=np.array(
                [index.ref_name_to_id[str(name)] for name in ref_names], dtype=np.int32
            ),
        )

    def footprint(self, t_idx: np.ndarray) -> tuple[tuple[Locus, ...], int]:
        """The transcripts' merged footprint: their spans sorted by ``(ref code, start)`` and merged where they
        overlap, as :class:`Locus` intervals, and its total span in bp."""
        rc, ss, ee = self.ref_codes[t_idx], self.starts[t_idx], self.ends[t_idx]
        order = np.lexsort((ss, rc))
        merged: list[tuple[int, int, int]] = []  # (ref_code, start, end)
        span = 0
        prev_rc, prev_s, prev_e = int(rc[order[0]]), int(ss[order[0]]), int(ee[order[0]])
        for j in order[1:]:
            rj, sj, ej = int(rc[j]), int(ss[j]), int(ee[j])
            if rj != prev_rc or sj > prev_e:
                merged.append((prev_rc, prev_s, prev_e))
                span += prev_e - prev_s
                prev_rc, prev_s, prev_e = rj, sj, ej
            elif ej > prev_e:
                prev_e = ej
        merged.append((prev_rc, prev_s, prev_e))
        span += prev_e - prev_s
        loci = tuple(
            Locus(
                ref=str(self.ref_names[rcode]),
                ref_id=int(self.ref_ids[rcode]),
                start=s,
                end=e,
            )
            for rcode, s, e in merged
        )
        return loci, span


def build_multi_loci(
    em_data: ScoredFragments,
    index: TranscriptIndex,
) -> list[MultiLocus]:
    """Build the multi-loci, the connected components of transcripts linked by fragments, by C++ union-find.

    Parameters
    ----------
    em_data : ScoredFragments
        Global EM data (annotated-transcript and synthetic-span candidates + per-unit metadata).
    index : TranscriptIndex
        Reference index.

    Returns
    -------
    list[MultiLocus]
    """
    if em_data.n_units == 0 or len(em_data.t_indices) == 0:
        return []

    # Per-component transcript and unit lists in CSR form, each sorted ascending.
    n_comp, comp_t_offsets, comp_t_indices, comp_u_offsets, comp_u_indices = _cc_native(
        em_data.offsets,
        em_data.t_indices,
        np.int32(index.num_transcripts),
    )
    spans = _TranscriptSpans.from_index(index)
    multi_loci: list[MultiLocus] = []
    for lid in range(n_comp):
        t_idx = comp_t_indices[comp_t_offsets[lid] : comp_t_offsets[lid + 1]].copy()
        loci, span = spans.footprint(t_idx)
        multi_loci.append(
            MultiLocus(
                multi_locus_id=lid,
                transcript_indices=t_idx,
                # Units are put in BAM order, not unit-index order: a unit's row in the scan buffer depends on the
                # order worker threads finished their chunks, while `frag_id` is the fragment's own identity. The
                # EM's sampled draws and float sums consume units in this order, so it must not vary from run to
                # run. The order must come from identity, never content: a content key ties on duplicate
                # fragments, which would then all draw alike.
                unit_indices=_in_bam_order(
                    comp_u_indices[comp_u_offsets[lid] : comp_u_offsets[lid + 1]],
                    em_data.frag_ids,
                ),
                span_bp=span,
                loci=loci,
            )
        )
    return multi_loci
