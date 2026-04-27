"""rigel.locus — Locus graph construction and per-locus EM priors.

* :func:`build_loci` — connected-component partitioning of transcripts
  linked by shared fragments, producing per-locus :class:`Locus` records.
* :func:`compute_locus_priors_from_partitions` — SRD v1 per-locus
  Dirichlet prior (``α_gdna``, ``α_rna``) computed from per-fragment
  posteriors on already-scored partitions; consumed by the C++ EM solver
  as a persistent M-step pseudocount and warm-start ratio.
"""

from __future__ import annotations

import numpy as np
from .native import connected_components as _cc_native

from .scored_fragments import Locus, LocusPartition, ScoredFragments

from .index import TranscriptIndex


# ---------------------------------------------------------------------------
# Locus builder: C++ union-find connected components
# ---------------------------------------------------------------------------


def build_loci(
    em_data: ScoredFragments,
    index: TranscriptIndex,
) -> list[Locus]:
    """Build loci as connected components of transcripts linked by fragments.

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
    list[Locus]
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
    ref_cat = index.t_df["ref"].cat
    _ref_names = ref_cat.categories.values
    _ref_codes = ref_cat.codes.values

    loci = []
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

        merged = []
        span = 0
        prev_rc = int(rc[order[0]])
        prev_s = int(ss[order[0]])
        prev_e = int(ee[order[0]])
        for k in range(1, len(order)):
            j = order[k]
            rj, sj, ej = int(rc[j]), int(ss[j]), int(ee[j])
            if rj != prev_rc or sj > prev_e:
                merged.append((_ref_names[prev_rc], prev_s, prev_e))
                span += prev_e - prev_s
                prev_rc, prev_s, prev_e = rj, sj, ej
            else:
                if ej > prev_e:
                    prev_e = ej
        merged.append((_ref_names[prev_rc], prev_s, prev_e))
        span += prev_e - prev_s

        loci.append(
            Locus(
                locus_id=lid,
                transcript_indices=t_idx,
                unit_indices=comp_u_indices[comp_u_offsets[lid] : comp_u_offsets[lid + 1]].copy(),
                gdna_span=max(span, 1),
                merged_intervals=merged,
            )
        )

    return loci


# ---------------------------------------------------------------------------
# SRD v1 per-locus Dirichlet priors from per-fragment posteriors
# ---------------------------------------------------------------------------


# Default Dirichlet evidence strength.  Matches the v5 ``gdna_prior_c_base``
# default; kept as a module constant rather than a config knob per the SRD
# simplicity constraint.  Promote to ``CalibrationConfig`` only if Phase 5
# benchmarks show sensitivity.
_C_BASE_DEFAULT = 5.0

# Floor on ``pi_pool`` so the log-prior never blows up when the calibrator
# returns 0 (QC-fail libraries) or 1 (degenerate gDNA-only inputs).
_PI_FLOOR = 1e-6


def compute_locus_priors_from_partitions(
    partitions: dict[int, LocusPartition],
    loci: list[Locus],
    pi_pool: float,
    *,
    c_base: float = _C_BASE_DEFAULT,
) -> tuple[np.ndarray, np.ndarray]:
    """SRD v1 per-locus Dirichlet priors from per-fragment posteriors.

    For each fragment ``f`` in each per-locus partition, compute the
    2-class gDNA posterior

    .. math::

        \\gamma_f = \\sigma\\!\\Big(
            \\log\\tfrac{\\pi_{\\text{pool}}}{1-\\pi_{\\text{pool}}}
            + \\log p_{\\text{gDNA}}(f)
            - \\log\\max_t p_{\\text{RNA}}(f \\mid t)
        \\Big)

    using the calibrated ``gdna_log_liks`` and per-candidate RNA
    ``log_liks`` already stored in the partition.  Aggregate to a
    per-locus mean ``\\gamma_\\ell`` and emit symmetric Dirichlet
    pseudocounts ``\\alpha_{gDNA}[\\ell] = c_{\\text{base}}\\gamma_\\ell``,
    ``\\alpha_{RNA}[\\ell] = c_{\\text{base}}(1 - \\gamma_\\ell)``.

    Parameters
    ----------
    partitions
        Per-locus partitions keyed by ``locus_id``.  Each partition's
        ``log_liks`` (per-candidate RNA log-likelihood) and
        ``gdna_log_liks`` (per-unit gDNA log-likelihood) must already
        reflect the calibrated ``gdna_fl_model`` / ``rna_fl_model``.
    loci
        Locus list defining iteration order; the returned arrays are
        indexed positionally (``[i]`` corresponds to ``loci[i]``).
    pi_pool
        Library-wide gDNA prior from
        :class:`~rigel.calibration.CalibrationResult.pi_pool`.
    c_base
        Dirichlet evidence strength.  Default 5.0.

    Returns
    -------
    alpha_gdna : np.ndarray, shape (n_loci,), float64
    alpha_rna : np.ndarray, shape (n_loci,), float64
    """
    n_loci = len(loci)
    alpha_gdna = np.zeros(n_loci, dtype=np.float64)
    alpha_rna = np.zeros(n_loci, dtype=np.float64)

    pi = float(np.clip(pi_pool, _PI_FLOOR, 1.0 - _PI_FLOOR))
    log_prior_g = np.log(pi) - np.log1p(-pi)

    for li, locus in enumerate(loci):
        part = partitions.get(locus.locus_id)
        if part is None or part.n_units == 0:
            continue

        offsets = np.asarray(part.offsets, dtype=np.intp)
        log_liks = np.asarray(part.log_liks, dtype=np.float64)
        gdna_log_liks = np.asarray(part.gdna_log_liks, dtype=np.float64)

        # Per-unit max RNA log-lik via reduceat over CSR offsets.
        max_rna = np.maximum.reduceat(log_liks, offsets[:-1])

        # γ_f = σ(z); use clipped z to avoid overflow in exp().
        z = log_prior_g + gdna_log_liks - max_rna
        z = np.clip(z, -50.0, 50.0)
        gamma_f = 1.0 / (1.0 + np.exp(-z))
        gamma_l = float(gamma_f.mean())

        alpha_gdna[li] = c_base * gamma_l
        alpha_rna[li] = c_base * (1.0 - gamma_l)

    return alpha_gdna, alpha_rna
