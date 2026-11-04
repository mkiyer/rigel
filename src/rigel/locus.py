"""rigel.locus — the EM's loci: connected components of transcripts linked by shared fragments.

* :class:`MultiLocus` — one connected component; the unit the EM is run on.
* :class:`Locus` — one contiguous interval of a multi-locus's merged transcript footprint.
  ``rigel.calibration.priors`` overlaps these intervals with the calibration regions to project
  calibration's per-region numbers onto the EM loci.
* :func:`build_multi_loci` — the connected-component partition, one :class:`MultiLocus` per
  component.

The EM's two pseudocounts are formed per locus by ``rigel.pipeline.em_pseudocounts``, from
calibration's per-locus gDNA count (``rigel.calibration.priors``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from .native import connected_components as _cc_native

from .scored_fragments import ScoredFragments

from .index import TranscriptIndex


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
    n_transcripts = index.num_transcripts
    offsets = em_data.offsets
    t_indices = em_data.t_indices
    n_units = em_data.n_units

    if n_units == 0 or len(t_indices) == 0:
        return []

    # Per-component transcript and unit lists in CSR form, each sorted ascending.
    n_comp, comp_t_offsets, comp_t_indices, comp_u_offsets, comp_u_indices = _cc_native(
        offsets,
        t_indices,
        np.int32(n_transcripts),
    )

    # Units are put in BAM order (`_in_bam_order`), not unit-index order: a unit's row in the scan buffer depends on
    # the order worker threads finished their chunks, while `frag_id` is the fragment's own identity. The EM's sampled
    # draws and float sums consume units in this order, so it must not vary from run to run. The order must come from
    # identity, never content: a content key ties on duplicate fragments, which would then all draw alike.

    t_starts_all = index.t_df["start"].values
    t_ends_all = index.t_df["end"].values

    # The "ref" column's categorical codes drive the merge loop's sort and compare. They are not the canonical ref-id
    # space (`index.ref_name_to_id`), which `Locus.ref_id` must carry: `calibration.priors` bins locus blocks by
    # ref_id against RegionArrays.
    ref_cat = index.t_df["ref"].cat
    _ref_names = ref_cat.categories.values
    _ref_codes = ref_cat.codes.values
    _cat_to_canonical_ref_id = np.array(
        [index.ref_name_to_id[str(name)] for name in _ref_names],
        dtype=np.int32,
    )

    multi_loci: list[MultiLocus] = []
    for lid in range(n_comp):
        t_lo = comp_t_offsets[lid]
        t_hi = comp_t_offsets[lid + 1]
        t_idx = comp_t_indices[t_lo:t_hi].copy()

        # Sort the transcript intervals by (ref_code, start) and merge overlaps into the genomic footprint.
        rc = _ref_codes[t_idx]
        ss = t_starts_all[t_idx]
        ee = t_ends_all[t_idx]
        order = np.lexsort((ss, rc))

        merged: list[tuple[int, int, int]] = []  # (ref_code, start, end)
        span = 0
        prev_rc = int(rc[order[0]])
        prev_s = int(ss[order[0]])
        prev_e = int(ee[order[0]])
        for k in range(1, len(order)):
            j = order[k]
            rj, sj, ej = int(rc[j]), int(ss[j]), int(ee[j])
            if rj != prev_rc or sj > prev_e:
                merged.append((prev_rc, prev_s, prev_e))
                span += prev_e - prev_s
                prev_rc, prev_s, prev_e = rj, sj, ej
            else:
                if ej > prev_e:
                    prev_e = ej
        merged.append((prev_rc, prev_s, prev_e))
        span += prev_e - prev_s

        loci_tuple = tuple(
            Locus(
                ref=str(_ref_names[rcode]),
                ref_id=int(_cat_to_canonical_ref_id[rcode]),
                start=s,
                end=e,
            )
            for rcode, s, e in merged
        )

        multi_loci.append(
            MultiLocus(
                multi_locus_id=lid,
                transcript_indices=t_idx,
                unit_indices=_in_bam_order(
                    comp_u_indices[comp_u_offsets[lid] : comp_u_offsets[lid + 1]],
                    em_data.frag_ids,
                ),
                span_bp=span,
                loci=loci_tuple,
            )
        )

    return multi_loci
