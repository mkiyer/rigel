"""rigel.locus — Locus graph construction.

* :class:`Locus` — one contiguous genomic interval (the calibration
  estimation unit; introduced in M5).
* :class:`MultiLocus` — one connected component of transcripts linked
  by shared fragments (the unit the EM is run on; was named ``Locus``
  prior to M5).
* :func:`build_multi_loci` — connected-component partitioning of
  transcripts producing per-component :class:`MultiLocus` records.

Per-locus EM Dirichlet priors are assembled by ``rigel.calibration.priors``.
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
    """One contiguous genomic interval — the calibration-estimation unit.

    A :class:`MultiLocus` is composed of one or more :class:`Locus`
    intervals; most ``MultiLocus``es have exactly one, but paralog
    clusters spanning multiple references carry several.
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

    The unit the EM is run on.  Composed of one or more :class:`Locus`
    intervals (the calibration-estimation unit).

    Attributes
    ----------
    multi_locus_id : int
        Sequential label (0-based).
    transcript_indices : np.ndarray
        int32 — global transcript indices in this multi-locus.
    unit_indices : np.ndarray
        int32 — EM unit indices (rows in global CSR) belonging to
        this multi-locus.
    gdna_span : int
        Total merged genomic footprint (bp).  Equal to
        ``sum(l.span for l in loci)`` (precomputed cache for the EM
        hot path).
    loci : tuple[Locus, ...]
        The contiguous intervals composing this multi-locus, sorted
        ascending by ``(ref_id, start)``.
    """

    multi_locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray
    gdna_span: int
    loci: tuple[Locus, ...]


# ---------------------------------------------------------------------------
# MultiLocus builder: C++ union-find connected components
# ---------------------------------------------------------------------------


def build_multi_loci(
    em_data: ScoredFragments,
    index: TranscriptIndex,
) -> list[MultiLocus]:
    """Build multi-loci as connected components of transcripts linked
    by fragments.

    Uses C++ union-find (disjoint-set with path compression and union
    by rank) for fast component detection.

    Parameters
    ----------
    em_data : ScoredFragments
        Global EM data (mRNA + nRNA candidates + per-unit metadata).
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

    # C++ union-find returns per-component transcript and unit lists
    # in CSR form (offsets + flat index arrays), already sorted ascending.
    n_comp, comp_t_offsets, comp_t_indices, comp_u_offsets, comp_u_indices = _cc_native(
        offsets,
        t_indices,
        np.int32(n_transcripts),
    )

    # Pre-extract transcript coordinates as numpy arrays.
    t_starts_all = index.t_df["start"].values
    t_ends_all = index.t_df["end"].values

    # Integer ref codes for fast sort/compare within merge loops.
    # The "ref" column is already categorical; use its codes directly
    # instead of an O(N log N) np.unique sort on 457K string objects.
    #
    # NOTE: pandas categorical codes are NOT the canonical resolver/BAM
    # ref-id space (which is defined by ``index.ref_lengths`` insertion
    # order via ``index.ref_name_to_id``).  We store both: ``_ref_codes``
    # is used only for fast intra-loop sorting/comparison; ``Locus.ref_id``
    # must carry the canonical id so downstream consumers
    # (``calibration.priors._project_regions_to_loci``, which bins locus
    # blocks by ref_id against RegionArrays) match the correct contig.
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

        # Sort transcript intervals by (ref_code, start), then merge
        # overlapping spans to compute the genomic footprint.
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
                unit_indices=comp_u_indices[comp_u_offsets[lid] : comp_u_offsets[lid + 1]].copy(),
                gdna_span=max(span, 1),
                loci=loci_tuple,
            )
        )

    return multi_loci


