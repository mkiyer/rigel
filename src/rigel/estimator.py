"""The abundance estimator: per-transcript and per-gene counts from the deterministic path and the locus EM.

The EM starts from a coverage-weighted seed and adds a gDNA-vs-RNA prior per locus: ``a_g`` on
gDNA and ``a_r`` on the RNA pool, shared among the RNA components in proportion to their current EM count unless
``rna_prior_weight`` names the shares (``apply_grouped_prior_update`` in em_solver.cpp). Production passes no
weights, so the distribution within the RNA pool is the data's.

``ScoredFragments`` and ``LocusPartition`` live in ``rigel.scored_fragments``, ``MultiLocus`` in ``rigel.locus``,
and the per-locus prior is assembled in ``rigel.calibration.priors``.
"""

import logging

import numpy as np
import pandas as pd

from .config import EMConfig, TranscriptGeometry
from .native import (
    ASSIGN_FRACTIONAL,
    ASSIGN_SAMPLE,
    EM_CLASSES,
)
from .native import batch_locus_em_partitioned as _batch_locus_em_partitioned
from .splice import NUM_SPLICE_STRAND_COLS, SPLICED_COLS

logger = logging.getLogger(__name__)

#: ``EMConfig``'s mode names → the codes the native EM defines (it refuses any other code).
ASSIGNMENT_MODE_CODES = {"fractional": ASSIGN_FRACTIONAL, "sample": ASSIGN_SAMPLE}


# ======================================================================
# Helpers
# ======================================================================


def _aggregate_to_gene(t_to_g: np.ndarray, n_genes: int, values: np.ndarray) -> np.ndarray:
    """Sum per-transcript *values* to gene level using np.add.at."""
    shape = (n_genes,) if values.ndim == 1 else (n_genes, values.shape[1])
    out = np.zeros(shape, dtype=np.float64)
    np.add.at(out, t_to_g, values)
    return out


def _gene_locus_ids(
    t_to_g: np.ndarray,
    n_genes: int,
    locus_id_per_transcript: np.ndarray,
    t_mrna: np.ndarray,
) -> np.ndarray:
    """For each gene, the locus_id of its highest-count transcript.

    The caller passes each transcript's total count as ``t_mrna``.
    """
    g_locus_id = np.full(n_genes, -1, dtype=np.int32)
    g_best_mrna = np.full(n_genes, -1.0, dtype=np.float64)
    valid = locus_id_per_transcript >= 0
    if not valid.any():
        return g_locus_id
    v_g = t_to_g[valid]
    v_mrna = t_mrna[valid]
    v_lid = locus_id_per_transcript[valid]
    # Descending count order; the first write per gene wins
    order = np.argsort(-v_mrna, kind="stable")
    for i in order:
        g = int(v_g[i])
        if v_mrna[i] > g_best_mrna[g]:
            g_best_mrna[g] = v_mrna[i]
            g_locus_id[g] = int(v_lid[i])
    return g_locus_id


#: Strand code (rigel.types.Strand: NONE=0, POS=1, NEG=2, AMBIGUOUS=3) → display glyph.
_STRAND_STR = np.array([".", "+", "-", "."])


def _strand_glyph(strand_codes: np.ndarray) -> np.ndarray:
    """Map an int strand-code array to ``+`` / ``-`` / ``.`` glyphs."""
    return _STRAND_STR[np.clip(np.asarray(strand_codes), 0, 3)]


def _tpm(
    counts: np.ndarray, eff_len: np.ndarray, denom_counts: np.ndarray | None = None
) -> np.ndarray:
    """Transcripts per million from count + effective-length arrays.

    If ``denom_counts`` is provided, TPM is normalized by the rpk sum of
    that alternate pool (used for ``tpm_total_rna``, where annotated
    counts are normalized over the annotated + synthetic pool).

    All inputs must share the same length; the denominator uses
    ``eff_len`` in both the numerator and denominator rpk.
    """
    rpk = counts / np.maximum(eff_len, 1.0)
    denom_rpk = rpk if denom_counts is None else denom_counts / np.maximum(eff_len, 1.0)
    denom = denom_rpk.sum()
    return (rpk / denom * 1e6) if denom > 0 else np.zeros_like(rpk)


def _span_children_counts(
    span_t_idx: np.ndarray,
    t_totals: np.ndarray,
    n_transcripts: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate multi-exon children to their span entity.

    For every transcript ``i`` with ``span_t_idx[i] >= 0`` and
    ``span_t_idx[i] != i`` (i.e., a *proper* child — not the self-reference
    of a single-exon transcript that covers a span), accumulate
    ``t_totals[i]`` into the entity's row.

    Returns
    -------
    n_children : int32[n_transcripts]
        Number of proper children per transcript.
    children_count : float64[n_transcripts]
        Sum of children's totals per transcript.
    """
    child_mask = (span_t_idx >= 0) & (span_t_idx != np.arange(n_transcripts))
    parents = span_t_idx[child_mask]
    n_children = np.bincount(parents, minlength=n_transcripts).astype(np.int32)
    children_count = np.bincount(
        parents,
        weights=t_totals[child_mask],
        minlength=n_transcripts,
    ).astype(np.float64)
    return n_children, children_count


def _per_locus(values, n_loci: int, name: str, *, fill: float) -> np.ndarray:
    """A per-locus float64 array for the native EM: *values*, checked to hold one entry per locus, or *fill* at every
    locus when ``None``."""
    if values is None:
        return np.full(n_loci, fill, dtype=np.float64)
    values = np.ascontiguousarray(values, dtype=np.float64)
    if values.shape[0] != n_loci:
        raise ValueError(f"{name} length {values.shape[0]} != n_loci {n_loci}.")
    return values


def _rna_prior_shares(weight, n_transcripts: int) -> np.ndarray:
    """``rna_prior_weight`` as the native EM reads it: per transcript, flat and whole, never indexed by locus; an empty
    array shares the RNA prior by current evidence."""
    if weight is None:
        return np.zeros(0, dtype=np.float64)
    weight = np.ascontiguousarray(weight, dtype=np.float64)
    if weight.shape[0] != n_transcripts:
        raise ValueError(
            f"rna_prior_weight length {weight.shape[0]} != n_transcripts "
            f"{n_transcripts}: the array is indexed per transcript, not per locus."
        )
    if not np.all(np.isfinite(weight)) or np.any(weight < 0.0):
        raise ValueError(
            "rna_prior_weight must be finite and non-negative — it is a share of a "
            "pseudo-fragment count."
        )
    return weight


def _gene_effective_lengths(
    t_to_g: np.ndarray, n_genes: int, t_counts: np.ndarray, t_eff: np.ndarray
) -> np.ndarray:
    """Each gene's effective length: its transcripts' effective lengths weighted by their counts, the unweighted mean
    for a gene with no count, floored at 1."""
    g_eff_num = _aggregate_to_gene(t_to_g, n_genes, t_counts * t_eff)
    g_eff_den = _aggregate_to_gene(t_to_g, n_genes, t_counts)
    g_eff_sum = _aggregate_to_gene(t_to_g, n_genes, t_eff)
    g_eff_cnt = _aggregate_to_gene(t_to_g, n_genes, np.ones_like(t_eff))
    with np.errstate(divide="ignore", invalid="ignore"):
        g_eff_mean = np.where(g_eff_cnt > 0, g_eff_sum / g_eff_cnt, 1.0)
        g_eff_len = np.where(g_eff_den > 0, g_eff_num / g_eff_den, g_eff_mean)
    return np.maximum(g_eff_len, 1.0)


#: The nascent table's columns, in order (:meth:`AbundanceEstimator.get_nascent_quant_df`).
_NASCENT_COLUMNS = (
    "nascent_id",
    "effective_length",
    "em_effective_length",
    "locus_id",
    "count",
    "n_annotated",
    "count_annotated",
    "tpm_total_rna",
)

#: The loci table's columns, in order (:meth:`AbundanceEstimator.get_loci_df`).
_LOCI_COLUMNS = (
    "locus_id",
    "locus_span_bp",
    "n_transcripts",
    "n_nascent",
    "n_genes",
    "n_em_fragments",
    "count_unambig",
    "count_annotated",
    "count_nascent",
    "count_gdna",
    "count",
    "fraction_gdna",
    "count_gdna_prior",
    "count_rna_prior",
    "has_gdna",
    "gdna_em_effective_length",
)


def _locus_row(
    r: dict,
    members: np.ndarray,
    t_total: np.ndarray,
    t_unambig: np.ndarray,
    is_synthetic: np.ndarray,
) -> dict:
    """One locus's row of the loci table, from its EM result *r* and its member transcripts.

    ``r["rna_total"]`` (em_solver.cpp, assign_posteriors) is the posterior count of every non-gDNA component in the
    locus, the annotated transcripts and the synthetic spans together. The table's counts do not overlap::

        count_nascent   = the synthetic spans' counts
        count_annotated = rna_total - count_nascent + unambig (annotated only)
        count           = count_annotated + count_nascent + count_gdna
    """
    syn_in_locus = members[is_synthetic[members]]
    annot_in_locus = members[~is_synthetic[members]]
    count_nascent = float(t_total[syn_in_locus].sum())
    n_nascent = int(syn_in_locus.size)
    # Unambiguous counts for annotated transcripts in this locus (splice-confirmed fragments that bypass EM)
    count_unambig = float(t_unambig[annot_in_locus].sum())
    rna_total = float(r["rna_total"])  # EM: annotated transcripts + synthetic spans
    annotated_em = max(rna_total - count_nascent, 0.0)  # the annotated transcripts' EM part
    count_annotated = annotated_em + count_unambig
    count_gdna = float(r["count_gdna"])
    count = count_annotated + count_nascent + count_gdna
    return {
        "locus_id": r["locus_id"],
        "locus_span_bp": r["locus_span_bp"],
        "n_transcripts": int(members.size) - n_nascent,
        "n_nascent": n_nascent,
        "n_genes": r["n_genes"],
        "n_em_fragments": r["n_em_fragments"],
        "count_unambig": count_unambig,
        "count_annotated": count_annotated,
        "count_nascent": count_nascent,
        "count_gdna": count_gdna,
        "count": count,
        "fraction_gdna": count_gdna / count if count > 0 else 0.0,
        "count_gdna_prior": float(r["count_gdna_prior"]),
        "count_rna_prior": float(r["count_rna_prior"]),
        "has_gdna": bool(r["has_gdna"]),
        "gdna_em_effective_length": float(r["gdna_em_effective_length"]),
    }


# ======================================================================
# AbundanceEstimator
# ======================================================================


class AbundanceEstimator:
    """Bayesian fragment abundance estimator with locus-level EM.

    Accumulates per-transcript and per-gene fragment counts from both the
    deterministic path (an unambiguous SPLICED_ANNOT fragment, assigned with no EM)
    and the probabilistic (EM) path.
    Estimates transcript abundances and each locus's gDNA count.

    Locus-level EM architecture::

        Per locus:
        [0, n_t)       - transcript components (annotated transcripts + synthetic nRNA spans)
        [n_t]          - gDNA: one component per locus

    Parameters
    ----------
    num_transcripts : int
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    geometry : TranscriptGeometry or None
        The effective lengths.  None gives every transcript length 1 (minimal construction in tests).
    """

    def __init__(
        self,
        num_transcripts: int,
        *,
        em_config: EMConfig | None = None,
        geometry: TranscriptGeometry | None = None,
        is_synthetic: np.ndarray | None = None,
    ):
        if em_config is None:
            em_config = EMConfig()
        self.em_config = em_config
        self.num_transcripts = num_transcripts
        self._rng = np.random.default_rng(em_config.seed)
        self._synthetic_mask = (
            np.asarray(is_synthetic, dtype=bool)
            if is_synthetic is not None
            else np.zeros(num_transcripts, dtype=bool)
        )

        # --- Effective lengths: from TranscriptGeometry, or 1 per transcript ---
        if geometry is not None:
            self._t_eff_len_output = np.maximum(
                np.asarray(geometry.effective_lengths, dtype=np.float64),
                1.0,
            )
            # The EM's yields, unfloored: a yield of 0 is a component that cannot emit (the C++ reads it
            # so) and every positive yield enters as it is — the E-step decides on ratios within a locus.
            em = np.asarray(geometry.effective_lengths_em, dtype=np.float64)
            if not np.all(np.isfinite(em)) or np.any(em < 0.0):
                raise ValueError("effective_lengths_em must be finite and non-negative")
            self._t_eff_len_em = em
        else:
            self._t_eff_len_output = np.ones(num_transcripts, dtype=np.float64)
            self._t_eff_len_em = self._t_eff_len_output

        self.unambig_counts = np.zeros((num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        self.em_counts = np.zeros((num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64)

        # --- gDNA: per locus, not per gene ---
        # Total gDNA count as assigned by locus EM (sum across all loci)
        self._gdna_em_total = 0.0
        # The intergenic term's mass: the multimappers' gDNA the EM put at a place outside every locus
        self._gdna_intergenic_em_total = 0.0

        # Per-locus results, one dict per locus (``pipeline.run_locus_em_partitioned``).
        self.locus_results: list[dict] = []

        # Per-transcript locus_id assignment (-1 = no locus).
        self.locus_id_per_transcript = np.full(num_transcripts, -1, dtype=np.int32)

        # Per-transcript confidence tracking:
        self._em_posterior_sum: np.ndarray | None = None
        self._em_n_assigned: np.ndarray | None = None

        # Per-locus profiling stats (populated if emit_locus_stats=True)
        self.locus_stats: list[dict] | None = None

    @property
    def effective_lengths(self) -> np.ndarray:
        """Per-transcript effective lengths (read-only)."""
        return self._t_eff_len_output

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA count across all loci."""
        return self._gdna_em_total

    @property
    def gdna_intergenic_em_count(self) -> float:
        """The EM's gDNA at intergenic places, outside every locus: intergenic gDNA beside the templates no hit of
        which resolves."""
        return self._gdna_intergenic_em_total

    @property
    def synthetic_count(self) -> float:
        """Total count on the synthetic spans: their unambiguous and EM counts, summed."""
        mask = self._synthetic_mask
        if not mask.any():
            return 0.0
        return float(self.t_counts[mask].sum())

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unambig + EM).

        Shape (N_t, NUM_SPLICE_STRAND_COLS) = (N_t, 8).
        """
        return self.unambig_counts + self.em_counts

    # ------------------------------------------------------------------
    # Batch locus EM (partition-native)
    # ------------------------------------------------------------------

    def run_batch_locus_em_partitioned(
        self,
        partition_tuples: list,
        locus_transcript_indices: list,
        gdna_prior_count: np.ndarray,
        *,
        rna_prior_count: np.ndarray | None = None,
        rna_prior_weight: np.ndarray | None = None,
        gdna_eff_len: np.ndarray | None = None,
        intergenic_rate: float = 0.0,
        em_iterations: int = EMConfig.iterations,
        em_convergence_delta: float = EMConfig.convergence_delta,
        emit_locus_stats: bool = False,
        emit_assignments: bool = False,
        emit_class_mass: bool = False,
    ) -> tuple[float, np.ndarray, np.ndarray, ...]:
        """Run locus-level EM from pre-partitioned data.

        Parameters
        ----------
        partition_tuples : list[tuple]
            List of 8-tuples, one per locus, containing partition arrays.
        locus_transcript_indices : list[np.ndarray]
            List of int32 transcript index arrays, one per locus.
        gdna_prior_count : np.ndarray
            float64 — additive aggregate gDNA alpha per locus.
        rna_prior_count : np.ndarray, optional
            float64 — additive aggregate RNA alpha per locus; the pipeline always
            passes it, and ``None`` means zeros.
        rna_prior_weight : np.ndarray, optional
            float64[n_transcripts] — how the per-locus RNA prior is shared among
            a locus's RNA components. ``None`` (the default) shares it in proportion
            to the evidence each component already carries. It is per transcript,
            like ``t_eff_lens``: passed to the C++ flat and whole, never indexed by
            ``multi_locus_id``, and the C++ refuses any length but ``n_transcripts`` or 0.
        gdna_eff_len : np.ndarray, optional
            float64 array of length ``n_loci``: the capture-contracted effective
            length of each locus's gDNA component (``LocusPriors.gdna_eff_len``,
            uncontracted off capture); the pipeline always passes it, and ``None``
            means ones.
        em_iterations, em_convergence_delta
            EM algorithm parameters.
        emit_locus_stats : bool
            If True, collect per-locus profiling stats.
        emit_assignments : bool
            If True, return per-unit winner_tid / posterior / n_candidates.
        emit_class_mass : bool
            If True, return every unit's posterior mass by class, ``(n_units, 4)`` float64 with the
            columns ``native.EM_CLASSES`` names (annotated transcripts, synthetic spans, the locus's
            gDNA, the intergenic term), from the posterior the assignment keeps or draws from: each
            row sums to 1 (0 for a unit no component can emit), and under fractional assignment
            the columns add up to the pools. Its rows are the units in the per-unit annotation's
            order, each locus's units in its partition's order. It changes no count.

        Returns
        -------
        tuple
            (total_gdna_em, locus_rna, locus_gdna[, winner_tid, winner_post, n_candidates][,
            class_mass]): the EM's total gDNA count, then per locus its RNA count (annotated
            transcripts plus synthetic nRNA spans) and its gDNA count.
        """
        n_transcripts = self.num_transcripts
        n_loci = len(partition_tuples)
        gdna_eff_len = _per_locus(gdna_eff_len, n_loci, "gdna_eff_len", fill=1.0)
        rna_prior_count = _per_locus(rna_prior_count, n_loci, "rna_prior_count", fill=0.0)
        t_rna_prior_weight = _rna_prior_shares(rna_prior_weight, n_transcripts)
        if self._em_posterior_sum is None:
            self._em_posterior_sum = np.zeros(n_transcripts, dtype=np.float64)
            self._em_n_assigned = np.zeros(n_transcripts, dtype=np.float64)

        rng_seed = int(self._rng.integers(0, 2**63))

        if not (np.isfinite(intergenic_rate) and intergenic_rate >= 0.0):
            raise ValueError(
                f"intergenic_rate must be finite and non-negative; got {intergenic_rate!r}."
            )
        # The EM's per-transcript effective lengths. The C++ EM subtracts log L_t per component in the E-step and in
        # assign_posteriors, in place of any per-fragment length correction.
        t_eff_lens = np.ascontiguousarray(self._t_eff_len_em, dtype=np.float64)
        (
            total_gdna_em,
            locus_rna,
            locus_gdna,
            locus_stats_raw,
            *assignments,  # each unit's winner, its posterior and its candidate count
            total_intergenic_em,
            class_mass,
        ) = _batch_locus_em_partitioned(
            partition_tuples,
            locus_transcript_indices,
            np.ascontiguousarray(gdna_prior_count, dtype=np.float64),
            rna_prior_count,
            gdna_eff_len,
            self.unambig_counts,
            t_eff_lens,
            t_rna_prior_weight,
            self.em_counts,
            self._em_posterior_sum,
            self._em_n_assigned,
            em_iterations,
            em_convergence_delta,
            self.em_config.mode == "vbem",
            ASSIGNMENT_MODE_CODES[self.em_config.assignment_mode],
            rng_seed,
            n_transcripts,
            NUM_SPLICE_STRAND_COLS,
            self.em_config.n_threads,
            emit_locus_stats,
            emit_assignments,
            float(intergenic_rate),
            np.ascontiguousarray(self._synthetic_mask, dtype=np.uint8)
            if emit_class_mass
            else np.zeros(0, dtype=np.uint8),
        )
        self._gdna_intergenic_em_total += float(total_intergenic_em)

        if emit_locus_stats:
            if self.locus_stats is None:
                self.locus_stats = []
            self.locus_stats.extend(locus_stats_raw)

        result = (total_gdna_em, np.asarray(locus_rna), np.asarray(locus_gdna))
        if emit_assignments:
            result = (*result, *(np.asarray(a) for a in assignments))
        if emit_class_mass:
            no_units = np.zeros((0, len(EM_CLASSES)), dtype=np.float64)
            result = (*result, no_units if class_mass is None else np.asarray(class_mass))
        return result

    # ------------------------------------------------------------------
    # Confidence / posterior metrics
    # ------------------------------------------------------------------

    def posterior_mean(self) -> np.ndarray:
        """Per-transcript mean posterior at assignment, ``Σp² / Σp`` over the EM units with a candidate at the
        transcript (in sample mode the units drawn to it), ``p`` each unit's posterior there. ``NaN`` where the EM
        assigned it nothing: a transcript whose every fragment is deterministic, or that no unit gave any posterior
        (a component the EM left at zero takes none)."""
        n_assigned = self._em_n_assigned
        if n_assigned is None:
            return np.full(self.num_transcripts, np.nan)
        post_sum = self._em_posterior_sum
        with np.errstate(divide="ignore", invalid="ignore"):
            result = np.where(
                n_assigned > 0,
                post_sum / n_assigned,
                np.nan,
            )
        return result

    # ------------------------------------------------------------------
    # Output - primary counts
    # ------------------------------------------------------------------

    def get_counts_df(self, index) -> pd.DataFrame:
        """Primary transcript-level abundance estimates: one row per annotated transcript.

        Every reference transcript is a row, single-exon included. The synthetic spans are not:
        they are the rows of the nascent table (:meth:`get_nascent_quant_df`).

        Columns
        -------
        transcript_id, gene_id, gene_name, gene_type : identifiers + biotype
        ref, strand, start, end, length : genomic locus + spliced length
        effective_length : FL-marginal effective length floored at 1 (TPM's length)
        em_effective_length : the EM's effective length
        locus_id : int32, EM locus (-1 if no locus)
        nascent_id : str, the transcript's synthetic span ("." when it has none: a
            single-exon transcript, or a multi-exon one whose span an annotated
            single-exon transcript covers)
        is_basic, is_mane : GENCODE flags
        is_single_exon : bool, the transcript has one exon
        count : fragment count assigned to this transcript
        count_unambig : uniquely-assigned count
        count_em : EM-assigned count
        count_spliced : spliced fragment count
        count_nascent : the synthetic span's own count (0 when ``nascent_id`` is ".")
        tpm : transcripts per million (annotated transcripts only)
        tpm_total_rna : TPM normalized over all RNA (annotated transcripts + synthetic spans)
        posterior_mean : mean posterior at EM assignment (:meth:`posterior_mean`), NaN where the EM assigned the
            transcript nothing
        """
        total = self.t_counts
        unique = self.unambig_counts

        t_total = total.sum(axis=1)
        count_unambig = unique.sum(axis=1)
        count_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)
        count_em = self.em_counts.sum(axis=1)

        # A synthetic span's count is the nascent table's: 0 here, and its row is dropped below.
        is_synthetic = index.t_df["is_synthetic"].values
        is_single_exon = index.t_df["is_single_exon"].values
        count = np.where(is_synthetic, 0.0, t_total)
        pmean = self.posterior_mean()

        # TPM:
        #   - ``tpm``           — annotated-only (synthetic counts forced to 0).
        #   - ``tpm_total_rna`` — annotated-only numerator, but normalized over
        #     annotated + synthetic so values are directly comparable to the
        #     nascent table's ``tpm_total_rna``.
        eff = self._t_eff_len_output
        em_eff = self._t_eff_len_em
        tpm = _tpm(count, eff)
        tpm_total_rna = _tpm(count, eff, denom_counts=t_total)

        # The nascent span: the span entity ``span_t_index`` names, when that entity is a
        # synthetic span. An annotated single-exon transcript that covers a multi-exon span is the
        # entity there (it names itself), and it is annotated RNA with a row of its own here, so a
        # transcript whose span it covers has no nascent span. The column is guaranteed by the
        # index schema (see TranscriptIndex.load).
        t_ids = index.t_df["t_id"].values
        span = index.t_df["span_t_index"].values
        has_nascent = (span >= 0) & is_synthetic[span.clip(0)]
        nascent_id = np.where(has_nascent, t_ids[span.clip(0)], ".")
        count_nascent = np.where(has_nascent, t_total[span.clip(0)], 0.0)

        df = pd.DataFrame(
            {
                "transcript_id": t_ids,
                "gene_id": index.t_df["g_id"].values,
                "gene_name": index.t_df["g_name"].values,
                "gene_type": index.t_df["g_type"].values,
                "ref": index.t_df["ref"].values,
                "strand": _strand_glyph(index.t_df["strand"].to_numpy()),
                "start": index.t_df["start"].values,
                "end": index.t_df["end"].values,
                "length": index.t_df["length"].values,
                "effective_length": eff,
                "em_effective_length": em_eff,
                "locus_id": self.locus_id_per_transcript,
                "nascent_id": nascent_id,
                "is_basic": index.t_df["is_basic"].values,
                "is_mane": index.t_df["is_mane"].values,
                "is_single_exon": is_single_exon,
                "count": count,
                "count_unambig": count_unambig,
                "count_em": count_em,
                "count_spliced": count_spliced,
                "count_nascent": count_nascent,
                "tpm": tpm,
                "tpm_total_rna": tpm_total_rna,
                "posterior_mean": pmean,
            }
        )
        # Drop the synthetic spans (they are the nascent table's rows)
        df = df[~is_synthetic].reset_index(drop=True)
        # Sort by transcript_id for deterministic output
        df = df.sort_values("transcript_id").reset_index(drop=True)
        return df

    def get_gene_counts_df(self, index) -> pd.DataFrame:
        """Primary gene-level abundance estimates.

        A gene's counts sum its annotated transcripts, single-exon included. The synthetic spans
        are gene-neutral (one span may model transcripts of more than one gene): they are the rows
        of the nascent table (:meth:`get_nascent_quant_df`), never attributed to a gene here.

        Columns
        -------
        gene_id, gene_name, gene_type : identifiers + biotype
        ref, strand, start, end : genomic locus
        n_transcripts : int, number of annotated transcripts
        locus_id : int32, primary EM locus for this gene (-1 if none)
        effective_length : abundance-weighted mean effective length
        count : sum of annotated transcript counts
        count_unambig : uniquely-assigned count
        count_em : EM-assigned count
        count_spliced : spliced fragment count
        tpm : transcripts per million
        """
        t_to_g = index.t_to_g_arr
        n_genes = index.num_genes
        # The annotated transcripts only: every synthetic span is its own gene row, so no aggregate puts its count on
        # an annotated gene, and the synthetic rows are dropped below
        annotated_mask = (~index.t_df["is_synthetic"].values).astype(np.float64)
        t_counts_all = self.t_counts.sum(axis=1)
        columns = {
            "n_transcripts": _aggregate_to_gene(t_to_g, n_genes, annotated_mask).astype(int),
            # the locus of the gene's highest-count transcript
            "locus_id": _gene_locus_ids(
                t_to_g, n_genes, self.locus_id_per_transcript, t_counts_all
            ),
            "effective_length": _gene_effective_lengths(
                t_to_g, n_genes, t_counts_all, self._t_eff_len_output
            ),
            "count": _aggregate_to_gene(t_to_g, n_genes, t_counts_all * annotated_mask),
            "count_unambig": _aggregate_to_gene(t_to_g, n_genes, self.unambig_counts).sum(axis=1),
            "count_em": _aggregate_to_gene(
                t_to_g, n_genes, self.em_counts.sum(axis=1) * annotated_mask
            ),
            "count_spliced": _aggregate_to_gene(t_to_g, n_genes, self.t_counts)[
                :, list(SPLICED_COLS)
            ].sum(axis=1),
        }
        annotated_gene = ~index.g_df["is_synthetic"].to_numpy()
        columns = {name: values[annotated_gene] for name, values in columns.items()}
        # Per-gene biotype (g_df carries no g_type): the first g_type per gene, keyed by the transcript→gene map
        g_type = (
            pd.Series(index.t_df["g_type"].to_numpy())
            .groupby(np.asarray(t_to_g))
            .first()
            .reindex(range(n_genes))
            .to_numpy()
        )[annotated_gene]

        df = pd.DataFrame(
            {
                "gene_id": index.g_df["g_id"].values[annotated_gene],
                "gene_name": index.g_df["g_name"].values[annotated_gene],
                "gene_type": g_type,
                "ref": index.g_df["ref"].values[annotated_gene],
                "strand": _strand_glyph(index.g_df["strand"].to_numpy()[annotated_gene]),
                "start": index.g_df["start"].values[annotated_gene],
                "end": index.g_df["end"].values[annotated_gene],
                **columns,
                "tpm": _tpm(columns["count"], columns["effective_length"]),
            }
        )
        # Sort by gene_id for deterministic output
        return df.sort_values("gene_id").reset_index(drop=True)

    # ------------------------------------------------------------------
    # Output - nascent RNA table
    # ------------------------------------------------------------------

    def get_nascent_quant_df(self, index) -> pd.DataFrame:
        """Nascent RNA table: one row per synthetic span.

        A synthetic span is the transcript the index adds to model the unspliced RNA of the
        multi-exon transcripts whose merged span it is, and those transcripts name it in
        ``span_t_index``. A merged span that an annotated single-exon transcript covers gets no
        synthetic span: that transcript is annotated RNA, a row of the transcript table, and is no
        row here. A synthetic span that no multi-exon transcript names is not a row either.

        Columns
        -------
        nascent_id : str, the synthetic span's transcript ID
        effective_length : FL-marginal effective length floored at 1 (TPM's length)
        em_effective_length : the EM's effective length
        locus_id : int32, EM locus
        count : float, the span's own fragment count
        n_annotated : int, the annotated transcripts the span models (the ones that name it)
        count_annotated : float, the sum of their counts
        tpm_total_rna : float, TPM over all RNA (comparable to the transcript table's)
        """
        t_df = index.t_df
        is_synthetic = t_df["is_synthetic"].values

        if not is_synthetic.any():
            return pd.DataFrame(columns=list(_NASCENT_COLUMNS))

        t_total = self.t_counts.sum(axis=1)
        # The span → named-by aggregation once (span_t_index is in the index schema), then the synthetic rows that at
        # least one multi-exon transcript names
        n_annotated_per_span, count_annotated_per_span = _span_children_counts(
            t_df["span_t_index"].values,
            t_total,
            n_transcripts=len(t_df),
        )
        idx = np.where(is_synthetic)[0]
        idx = idx[n_annotated_per_span[idx] > 0]
        counts = t_total[idx]
        eff = self._t_eff_len_output[idx]

        # Total RNA TPM, normalized over the annotated + synthetic rpk pool so it is comparable to the transcript
        # table's tpm_total_rna; inline, because the denominator is the full pool, not this array
        rpk = counts / np.maximum(eff, 1.0)
        total_rpk = float((t_total / np.maximum(self._t_eff_len_output, 1.0)).sum())
        tpm_total_rna = (rpk / total_rpk * 1e6) if total_rpk > 0 else np.zeros_like(rpk)

        df = pd.DataFrame(
            {
                "nascent_id": t_df["t_id"].values[idx],
                "effective_length": eff,
                "em_effective_length": self._t_eff_len_em[idx],
                "locus_id": self.locus_id_per_transcript[idx],
                "count": counts,
                "n_annotated": n_annotated_per_span[idx],
                "count_annotated": count_annotated_per_span[idx],
                "tpm_total_rna": tpm_total_rna,
            },
            columns=list(_NASCENT_COLUMNS),
        )
        return df.sort_values("nascent_id").reset_index(drop=True)

    # ------------------------------------------------------------------
    # Output - locus-level breakdown
    # ------------------------------------------------------------------

    def get_loci_df(self, index) -> pd.DataFrame:
        """Locus-level output: one row per EM locus, its count by population.

        Returns an empty DataFrame (with correct columns) if no loci
        were processed.

        Columns
        -------
        locus_id : int, sequential locus identifier
        locus_span_bp : int, merged genomic footprint in base pairs
        n_transcripts : int, annotated transcripts in the locus
        n_nascent : int, synthetic spans in the locus
        n_genes : int, annotated genes
        n_em_fragments : int, ambiguous fragments entering EM
        count_unambig : float, unambiguous (splice-confirmed) fragments on the annotated transcripts
        count_annotated : float, the annotated transcripts' count (EM + unambig)
        count_nascent : float, the synthetic spans' count
        count_gdna : float, gDNA count from the locus EM
        count : float, count_annotated + count_nascent + count_gdna
        fraction_gdna : float, count_gdna / count
        count_gdna_prior : float, per-locus gDNA Dirichlet prior scalar
        count_rna_prior : float, per-locus RNA Dirichlet prior scalar
        has_gdna : bool, whether the locus admits a gDNA component
        gdna_em_effective_length : float
            Exposure-adjusted EM effective length for the locus gDNA component.
        """
        if not self.locus_results:
            return pd.DataFrame(columns=list(_LOCI_COLUMNS))

        t_total = self.t_counts.sum(axis=1)
        t_unambig = self.unambig_counts.sum(axis=1)
        is_synthetic = index.t_df["is_synthetic"].values
        locus_ids = self.locus_id_per_transcript
        # Transcripts grouped by locus once, by a stable sort that keeps each locus's members in index order, so
        # every sum adds the values a per-locus mask would select, in the same order
        order = np.argsort(locus_ids, kind="stable")
        sorted_ids = locus_ids[order]
        # Every locus's member range in two vectorised searches rather than one call per locus
        lids = np.array([r["locus_id"] for r in self.locus_results], dtype=sorted_ids.dtype)
        first = np.searchsorted(sorted_ids, lids, side="left")
        last = np.searchsorted(sorted_ids, lids, side="right")
        rows = [
            _locus_row(r, order[lo:hi], t_total, t_unambig, is_synthetic)
            for r, lo, hi in zip(self.locus_results, first, last)
        ]
        return pd.DataFrame(rows, columns=list(_LOCI_COLUMNS))
