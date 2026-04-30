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


# Default Dirichlet evidence strength for the per-locus (gDNA, RNA)
# prior derived from per-fragment SRD posteriors.  Empirically tuned;
# kept as a module constant rather than a config knob per the SRD
# simplicity constraint.
_C_BASE_DEFAULT = 10.0

# Floor on ``pi_pool`` so the log-prior never blows up when the calibrator
# returns 0 (QC-fail libraries) or 1 (degenerate gDNA-only inputs).
_PI_FLOOR = 1e-6


def compute_locus_priors_from_partitions(
    partitions: dict[int, LocusPartition],
    loci: list[Locus],
    pi_pool: float,
    *,
    c_base: float | None = None,  # retained for back-compat (vestigial; see Notes)
) -> tuple[np.ndarray, np.ndarray]:
    """SRD v3 per-locus Dirichlet priors from per-fragment posteriors,
    restricted to gDNA-eligible units.

    For each unit ``f`` in each per-locus partition that carries a finite
    ``gdna_log_lik`` (i.e. unspliced units that are scored against the
    gDNA hypothesis), compute the 2-class gDNA posterior

    .. math::

        \\gamma_f = \\sigma\\!\\Big(
            \\log\\tfrac{\\pi_{\\text{pool}}}{1-\\pi_{\\text{pool}}}
            + \\log p_{\\text{gDNA}}(f)
            - \\log\\max_t p_{\\text{RNA}}(f \\mid t)
        \\Big)

    and emit per-locus pseudocounts equal to the *sum* of those
    posteriors over the eligible subset:

    .. math::

        \\alpha_{gDNA}[\\ell] = \\sum_{f \\in E_\\ell} \\gamma_f , \\qquad
        \\alpha_{RNA}[\\ell]  = \\sum_{f \\in E_\\ell} (1 - \\gamma_f) ,

    where ``E_ℓ`` is the gDNA-eligible (unspliced) subset of locus
    ``ℓ``.  This expresses the prior as the *expected count* of gDNA
    fragments inside the eligible subset.

    Notes
    -----
    Why restrict to gDNA-eligible units?  The previous SRD-v1 form
    ``α_gDNA = c_base · mean_f γ_f`` averaged over **all** units in the
    locus.  Spliced / exonic units have ``gdna_log_lik = -inf`` ⇒
    ``γ_f = 0``, which dragged the mean toward zero whenever the locus
    was dominated by exonic mRNA (the typical case).  The resulting
    ``α_gDNA`` was tiny relative to ``α_RNA``, and because intronic
    gDNA-vs-nRNA per-fragment likelihoods are nearly degenerate
    (especially for unstranded data), the locus EM defaulted to the
    prior and siphoned ~95-99% of intronic gDNA into the synthetic
    nRNA component.  See ``docs/calibration/srd_v3_phase2_oracle_findings.md``
    and ``scripts/debug/nrna_siphon_minimal.py`` (experiments E5–E7)
    for the diagnostic chain.

    Why is the formulation **asymmetric** in scope (mean over eligible
    only) but symmetric in scale (a single ``c_base``)?  Because the
    calibrator only estimates the gDNA share of the eligible pool —
    averaging into spliced/exonic units (whose ``γ_f = 0`` by
    construction) is what diluted ``α_gDNA`` to ineffectual values.
    Using the eligible-only mean preserves the *ratio* the calibrator
    intends.  Keeping the *strength* at a small constant (``c_base``)
    avoids double-counting the per-fragment evidence that already
    shapes ``γ_f`` — sum-formulations effectively fold in the data
    twice and over-pull the EM toward gDNA at heavy contamination,
    breaking the paralog tie-breaker and identical-sequence loci.

    The default ``c_base = 10.0`` (raised from SRD v1's 5.0) reflects
    a compromise calibrated against both the synthetic siphon harness
    (``scripts/debug/nrna_siphon_minimal.py`` E7) and the scenario
    test suite.  Higher values resolve the synthetic siphon more
    aggressively but inflate the nRNA prior in pure-mRNA loci because
    of the C++ side's coverage-proportional ``α_RNA`` allocator
    (``em_solver.cpp:838``).  Resolving the residual siphon below
    ~50% requires changes to that allocator (alpha_rna should *not*
    flow into the synthetic nRNA component when there is no
    independent gDNA-eligible evidence for it); see TODO in
    ``docs/calibration/srd_v3_phase2_nrna_siphon.md``.

    Parameters
    ----------
    partitions
        Per-locus partitions keyed by ``locus_id``.
    loci
        Locus list defining iteration order.
    pi_pool
        Library-wide gDNA prior from
        :class:`~rigel.calibration.CalibrationResult.pi_pool`.
    c_base
        Per-locus prior strength (default 10.0).

    Returns
    -------
    alpha_gdna : np.ndarray, shape (n_loci,), float64
    alpha_rna  : np.ndarray, shape (n_loci,), float64
    """
    n_loci = len(loci)
    alpha_gdna = np.zeros(n_loci, dtype=np.float64)
    alpha_rna = np.zeros(n_loci, dtype=np.float64)

    pi = float(np.clip(pi_pool, _PI_FLOOR, 1.0 - _PI_FLOOR))
    log_prior_g = np.log(pi) - np.log1p(-pi)

    cb = _C_BASE_DEFAULT if c_base is None else float(c_base)

    for li, locus in enumerate(loci):
        part = partitions.get(locus.locus_id)
        if part is None or part.n_units == 0:
            continue

        offsets = np.asarray(part.offsets, dtype=np.intp)
        log_liks = np.asarray(part.log_liks, dtype=np.float64)
        gdna_log_liks = np.asarray(part.gdna_log_liks, dtype=np.float64)

        eligible = np.isfinite(gdna_log_liks)
        if not eligible.any():
            # No eligible (unspliced) units: emit the calibrated symmetric
            # prior so EM can still proceed if needed.
            alpha_gdna[li] = cb * pi
            alpha_rna[li] = cb * (1.0 - pi)
            continue

        # Per-unit max RNA log-lik via reduceat over CSR offsets.
        max_rna = np.maximum.reduceat(log_liks, offsets[:-1])

        gdna_e = gdna_log_liks[eligible]
        max_rna_e = max_rna[eligible]
        # Clip the diff before adding the prior to dodge ±inf
        # (units whose only candidate is gDNA → max_rna = -inf).
        diff = np.where(np.isfinite(max_rna_e), gdna_e - max_rna_e, 50.0)
        z = log_prior_g + diff
        z = np.clip(z, -50.0, 50.0)
        gamma_e = 1.0 / (1.0 + np.exp(-z))

        # Eligible-only mean preserves the calibrator's intended ratio
        # without diluting against spliced/exonic units (where γ_f = 0).
        gamma_l = float(gamma_e.mean())
        alpha_gdna[li] = cb * gamma_l
        alpha_rna[li] = cb * (1.0 - gamma_l)

    return alpha_gdna, alpha_rna
