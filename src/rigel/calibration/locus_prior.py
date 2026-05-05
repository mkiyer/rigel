"""Per-locus prior helpers for the EM solver.

Currently exposes :func:`build_prior_weight_rna`, the canonical
constructor for the ``prior_weight_rna`` vector consumed by
``run_locus_em_native`` / ``batch_locus_em_partitioned``.

The component layout for the EM is ``[t0, t1, ..., t_{n_t-1}, gDNA]``
(length ``n_t + 1``).  The gDNA entry is ignored by the solver — it is
routed through ``alpha_gdna`` — but the array is sized to
``n_components`` for layout symmetry.

When per-transcript nRNA components (M6/M8) land, this helper will gain
a ``nrna_weight`` parameter that maps ``[mRNA_t0, nRNA_t0, ..., gDNA]``
component indices onto the right multipliers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..scored_fragments import Locus  # noqa: F401  (renamed to MultiLocus in Part B)


def build_prior_weight_rna(
    multi_locus,
    em_data=None,  # reserved for M6: per-transcript nRNA component info
    *,
    nrna_weight: float = 0.0,  # noqa: ARG001 (reserved for M6/M8)
) -> np.ndarray:
    """Construct the per-component nRNA-suppression weight vector.

    Parameters
    ----------
    multi_locus
        The locus the weight vector is being built for; only its
        ``transcript_indices`` length is read here.
    em_data
        Reserved.  Will be used by M6 to identify which components are
        synthetic nRNA shadows.
    nrna_weight
        Reserved.  Will be applied to nRNA components by M6/M8.

    Returns
    -------
    np.ndarray
        ``float32`` array of length ``n_components = n_transcripts + 1``,
        all entries set to ``1.0``.  Passing this to the EM is
        bit-identical to passing ``None``.
    """
    n_t = int(multi_locus.transcript_indices.shape[0])
    return np.ones(n_t + 1, dtype=np.float32)
