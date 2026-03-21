"""rigel.estimator — Bayesian fragment abundance estimator.

``AbundanceEstimator`` manages per-transcript and per-gene fragment
count accumulation combined with Bayesian abundance estimation via
locus-level EM.

For unambiguously-mapping fragments, assignment is deterministic.
For ambiguous fragments, transcript abundances are estimated via
per-locus EM with a coverage-weighted One Virtual Read (OVR) prior,
then counts are accumulated from the converged posterior.

Data containers (``ScoredFragments``, ``Locus``, ``LocusEMInput``)
live in ``rigel.scored_fragments``.  Calibration and prior computation
live in ``rigel.calibration`` and ``rigel.locus``.
"""

import logging

import numpy as np
import pandas as pd

from .config import EMConfig, TranscriptGeometry
from .native import batch_locus_em as _batch_locus_em
from .scored_fragments import ScoredFragments
from .splice import (
    SpliceType,
    SpliceStrandCol,
    NUM_SPLICE_STRAND_COLS,
    SPLICED_COLS,
)

logger = logging.getLogger(__name__)


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
    """For each gene, assign the locus_id of its highest-mRNA transcript."""
    g_locus_id = np.full(n_genes, -1, dtype=np.int32)
    g_best_mrna = np.full(n_genes, -1.0, dtype=np.float64)
    valid = locus_id_per_transcript >= 0
    if not valid.any():
        return g_locus_id
    v_g = t_to_g[valid]
    v_mrna = t_mrna[valid]
    v_lid = locus_id_per_transcript[valid]
    # Process in descending mRNA order; first write per gene wins.
    order = np.argsort(-v_mrna, kind="stable")
    for i in order:
        g = int(v_g[i])
        if v_mrna[i] > g_best_mrna[g]:
            g_best_mrna[g] = v_mrna[i]
            g_locus_id[g] = int(v_lid[i])
    return g_locus_id


# ======================================================================
# AbundanceEstimator
# ======================================================================


class AbundanceEstimator:
    """Bayesian fragment abundance estimator with locus-level EM.

    Accumulates per-transcript and per-gene fragment counts from both
    deterministic (unambig) and probabilistic (EM) assignment paths.
    Estimates transcript-level mRNA, nRNA, and gDNA abundances.

    Locus-level EM architecture::

        Per locus:
        [0, n_t)       - mRNA transcript components
        [n_t, 2*n_t)   - nRNA shadow per transcript
        [2*n_t]        - gDNA: ONE shadow per locus

    Fragment routing:

    - SPLICED_ANNOT + unambig -> deterministic mRNA (no EM)
    - All other -> enter locus EM
      - Spliced fragments: mRNA + nRNA candidates only
      - Unspliced fragments: mRNA + nRNA + gDNA candidates

    Parameters
    ----------
    num_transcripts : int
    em_config : EMConfig or None
        EM algorithm configuration.  Defaults to ``EMConfig()``.
    geometry : TranscriptGeometry or None
        Pre-computed transcript/gene geometry.  None for minimal
        construction (e.g. in tests).
    """

    def __init__(
        self,
        num_transcripts: int,
        *,
        num_nrna: int = 0,
        t_to_nrna: np.ndarray | None = None,
        em_config: EMConfig | None = None,
        geometry: TranscriptGeometry | None = None,
    ):
        if em_config is None:
            em_config = EMConfig()
        self.em_config = em_config
        self.num_transcripts = num_transcripts
        self.num_nrna = num_nrna if num_nrna > 0 else num_transcripts
        self._t_to_nrna = t_to_nrna  # global transcript → global nRNA mapping
        self._rng = np.random.default_rng(em_config.seed)

        # nRNA shadow base index (for global CSR component numbering).
        self.nrna_base_index = num_transcripts

        # --- Geometry: from TranscriptGeometry or defaults ---
        if geometry is not None:
            self._t_to_g = np.asarray(geometry.t_to_g, dtype=np.int32)
            self._transcript_spans = np.asarray(geometry.transcript_spans, dtype=np.float64)
            self._exonic_lengths = np.asarray(geometry.exonic_lengths, dtype=np.float64)
            self._t_eff_len = np.maximum(
                np.asarray(geometry.effective_lengths, dtype=np.float64),
                1.0,
            )
        else:
            self._t_to_g = None
            self._transcript_spans = None
            self._exonic_lengths = None
            self._t_eff_len = np.ones(num_transcripts, dtype=np.float64)

        self.unambig_counts = np.zeros((num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64)
        self.em_counts = np.zeros((num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64)

        # Per-transcript high-confidence EM counts (posterior >= threshold).
        self.em_high_conf_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # Per-nRNA shadow initialization and EM counts.
        self.nrna_em_counts = np.zeros(self.num_nrna, dtype=np.float64)

        # --- gDNA: locus-level, NOT per-gene ---
        # Total gDNA count as assigned by locus EM (sum across all loci)
        self._gdna_em_total = 0.0

        # Per-locus results (locus_id, ref, n_transcripts, n_genes,
        # n_em_fragments, mrna, nrna, gdna, gdna_init).
        self.locus_results: list[dict] = []

        # Per-transcript locus_id assignment (-1 = no locus).
        self.locus_id_per_transcript = np.full(num_transcripts, -1, dtype=np.int32)

        # Per-transcript gDNA locus attribution (for reporting).
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # --- Pre-EM strand accumulators (for gDNA init via EB) ---
        # Per-nRNA: all same-strand UNSPLICED fragments.
        # Accumulated during scan pass for empirical Bayes gDNA prior.
        # Locus-level totals are derived by summing over locus.transcript_indices.
        self.transcript_unspliced_sense = np.zeros(self.num_nrna, dtype=np.float64)
        self.transcript_unspliced_antisense = np.zeros(self.num_nrna, dtype=np.float64)

        # --- Pre-EM intronic accumulators (for nRNA init) ---
        self.transcript_intronic_sense = np.zeros(self.num_nrna, dtype=np.float64)
        self.transcript_intronic_antisense = np.zeros(self.num_nrna, dtype=np.float64)

        # Per-transcript confidence tracking:
        self._em_posterior_sum: np.ndarray | None = None
        self._em_n_assigned: np.ndarray | None = None

    @property
    def effective_lengths(self) -> np.ndarray:
        """Per-transcript effective lengths (read-only)."""
        return self._t_eff_len

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA count across all loci."""
        return self._gdna_em_total

    @property
    def nrna_em_count(self) -> float:
        """Total EM-assigned nRNA count."""
        return float(self.nrna_em_counts.sum())

    @property
    def gdna_contamination_rate(self) -> float:
        """Fraction of all fragments attributed to gDNA."""
        total_rna = float(self.unambig_counts.sum() + self.em_counts.sum())
        total = total_rna + self._gdna_em_total
        if total == 0:
            return 0.0
        return self._gdna_em_total / total

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unambig + EM), shape (N_t, 8)."""
        return self.unambig_counts + self.em_counts

    @property
    def t_high_conf_counts(self) -> np.ndarray:
        """High-confidence transcript counts (unambig + EM above threshold)."""
        return self.unambig_counts + self.em_high_conf_counts

    # ------------------------------------------------------------------
    # Batch locus EM — Phase 2A: single C++ call for all loci
    # ------------------------------------------------------------------

    def run_batch_locus_em(
        self,
        loci: list,
        em_data: ScoredFragments,
        index,
        locus_gammas: np.ndarray,
        *,
        em_iterations: int = 1000,
        em_convergence_delta: float = 1e-6,
        confidence_threshold: float = 0.95,
    ) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        """Run locus-level EM for ALL loci in a single C++ call.

        Replaces the Python per-locus for-loop
        (build_locus_em_data → run_locus_em → assign_locus_ambiguous)
        with one batched C++ call to ``batch_locus_em()``.

        Parameters
        ----------
        loci : list[Locus]
            All loci to process.
        em_data : ScoredFragments
            Global CSR data.
        index : TranscriptIndex
            Reference index.
        locus_gammas : np.ndarray
            float64[n_loci] — calibration gDNA fraction ∈ [0, 1] per locus.
        em_iterations, em_convergence_delta, confidence_threshold
            EM algorithm parameters.

        Returns
        -------
        tuple[float, np.ndarray, np.ndarray, np.ndarray]
            (total_gdna_em, locus_mrna, locus_nrna, locus_gdna)
        """
        n_loci = len(loci)
        n_transcripts = self.num_transcripts

        # Flatten locus definitions into contiguous CSR arrays.
        locus_t_offsets = np.empty(n_loci + 1, dtype=np.int64)
        locus_u_offsets = np.empty(n_loci + 1, dtype=np.int64)
        locus_t_offsets[0] = 0
        locus_u_offsets[0] = 0
        for i, locus in enumerate(loci):
            locus_t_offsets[i + 1] = locus_t_offsets[i] + len(locus.transcript_indices)
            locus_u_offsets[i + 1] = locus_u_offsets[i] + len(locus.unit_indices)

        total_t = int(locus_t_offsets[-1])
        total_u = int(locus_u_offsets[-1])
        locus_t_flat = np.empty(total_t, dtype=np.int32)
        locus_u_flat = np.empty(total_u, dtype=np.int32)
        for i, locus in enumerate(loci):
            ts = int(locus_t_offsets[i])
            te = int(locus_t_offsets[i + 1])
            us = int(locus_u_offsets[i])
            ue = int(locus_u_offsets[i + 1])
            locus_t_flat[ts:te] = locus.transcript_indices
            locus_u_flat[us:ue] = locus.unit_indices

        # Per-transcript geometry arrays
        t_starts = index.t_df["start"].values.astype(np.int64)
        t_ends = index.t_df["end"].values.astype(np.int64)
        t_lengths = index.t_df["length"].values.astype(np.int64)

        # Prepare mutable output accumulators
        # (em_counts, em_high_conf_counts already on self)
        if self._em_posterior_sum is None:
            self._em_posterior_sum = np.zeros(n_transcripts, dtype=np.float64)
            self._em_n_assigned = np.zeros(n_transcripts, dtype=np.float64)

        # Map assignment_mode string to int for C++: 0=fractional, 1=map, 2=sample
        _ASSIGNMENT_MODE_MAP = {"fractional": 0, "map": 1, "sample": 2}
        assignment_mode_int = _ASSIGNMENT_MODE_MAP[self.em_config.assignment_mode]

        # Derive a per-batch RNG seed from the Python-level RNG
        rng_seed = int(self._rng.integers(0, 2**63))

        total_gdna_em, locus_mrna, locus_nrna, locus_gdna = _batch_locus_em(
            # Global CSR
            em_data.offsets,
            em_data.t_indices,
            em_data.log_liks,
            em_data.coverage_weights,
            em_data.tx_starts,
            em_data.tx_ends,
            em_data.count_cols,
            em_data.is_spliced.view(np.uint8),
            em_data.gdna_log_liks,
            em_data.genomic_footprints,
            em_data.locus_t_indices,
            em_data.locus_count_cols,
            np.int32(em_data.nrna_base_index),
            # Locus definitions
            locus_t_offsets,
            locus_t_flat,
            locus_u_offsets,
            locus_u_flat,
            np.ascontiguousarray(locus_gammas, dtype=np.float64),
            # Pre-computed union genomic footprints
            np.array([loc.gdna_span for loc in loci], dtype=np.int64),
            # Per-transcript data
            self.unambig_counts,
            t_starts,
            t_ends,
            t_lengths,
            self._transcript_spans,
            self._exonic_lengths,
            # nRNA mapping
            index.t_to_nrna_arr,
            index.num_nrna,
            # Mutable output accumulators
            self.em_counts,
            self.em_high_conf_counts,
            self.nrna_em_counts,
            self.gdna_locus_counts,
            self._em_posterior_sum,
            self._em_n_assigned,
            # EM config
            em_iterations,
            em_convergence_delta,
            self.em_config.prior_pseudocount,
            self.em_config.mode == "vbem",
            confidence_threshold,
            assignment_mode_int,
            self.em_config.assignment_min_posterior,
            rng_seed,
            n_transcripts,
            NUM_SPLICE_STRAND_COLS,
            self.em_config.n_threads,
        )

        return (
            total_gdna_em,
            np.asarray(locus_mrna),
            np.asarray(locus_nrna),
            np.asarray(locus_gdna),
        )

    # ------------------------------------------------------------------
    # Confidence / posterior metrics
    # ------------------------------------------------------------------

    def posterior_mean(self) -> np.ndarray:
        """Per-transcript mean posterior at assignment."""
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
    # Helpers
    # ------------------------------------------------------------------

    def _nrna_per_transcript(self) -> np.ndarray:
        """Fan out per-nRNA EM counts to per-transcript (equal share)."""
        if self._t_to_nrna is not None:
            fan_counts = np.bincount(
                self._t_to_nrna,
                minlength=self.num_nrna,
            ).astype(np.float64)
            fan_counts = np.maximum(fan_counts, 1.0)
            return self.nrna_em_counts[self._t_to_nrna] / fan_counts[self._t_to_nrna]
        return self.nrna_em_counts

    # ------------------------------------------------------------------
    # Output - primary counts
    # ------------------------------------------------------------------

    def get_counts_df(self, index) -> pd.DataFrame:
        """Primary transcript-level abundance estimates.

        Columns
        -------
        transcript_id, gene_id, gene_name : identifiers
        locus_id : int32, EM locus (-1 if no locus)
        effective_length : bias-corrected effective length
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_high_conf : high-confidence mRNA (unambig + EM >= threshold)
        mrna_spliced : spliced mRNA fragments
        nrna : nascent RNA count (EM-assigned only)
        rna_total : mrna + nrna
        tpm : transcripts per million (mRNA-based)
        posterior_mean : mean posterior at EM assignment
        """
        total = self.t_counts
        unique = self.unambig_counts

        mrna = total.sum(axis=1)
        mrna_unambig = unique.sum(axis=1)
        mrna_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)
        mrna_em = self.em_counts.sum(axis=1)
        mrna_high_conf = self.t_high_conf_counts.sum(axis=1)
        nrna = self._nrna_per_transcript()
        rna_total = mrna + nrna
        pmean = self.posterior_mean()

        # TPM: mRNA-based, using effective lengths
        eff = self._t_eff_len
        rpk = mrna / eff  # reads per kilobase (un-scaled)
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        df = pd.DataFrame(
            {
                "transcript_id": index.t_df["t_id"].values,
                "gene_id": index.t_df["g_id"].values,
                "gene_name": index.t_df["g_name"].values,
                "locus_id": self.locus_id_per_transcript,
                "effective_length": eff,
                "mrna": mrna,
                "mrna_unambig": mrna_unambig,
                "mrna_em": mrna_em,
                "mrna_high_conf": mrna_high_conf,
                "mrna_spliced": mrna_spliced,
                "nrna": nrna,
                "rna_total": rna_total,
                "tpm": tpm,
                "posterior_mean": pmean,
            }
        )
        return df

    def get_gene_counts_df(self, index) -> pd.DataFrame:
        """Primary gene-level abundance estimates.

        Columns
        -------
        gene_id, gene_name : identifiers
        locus_id : int32, primary EM locus for this gene (-1 if none)
        effective_length : abundance-weighted mean effective length
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_high_conf : high-confidence mRNA
        mrna_spliced : spliced mRNA fragments
        nrna : nascent RNA count
        rna_total : mrna + nrna
        tpm : transcripts per million (mRNA-based)
        """
        t_to_g = index.t_to_g_arr
        total = self.t_counts
        unique = self.unambig_counts
        n_genes = index.num_genes

        g_total = _aggregate_to_gene(t_to_g, n_genes, total)
        g_unambig = _aggregate_to_gene(t_to_g, n_genes, unique)

        mrna = g_total.sum(axis=1)
        mrna_unambig = g_unambig.sum(axis=1)
        mrna_spliced = g_total[:, list(SPLICED_COLS)].sum(axis=1)

        mrna_em_arr = _aggregate_to_gene(t_to_g, n_genes, self.em_counts.sum(axis=1))
        mrna_hc_arr = _aggregate_to_gene(t_to_g, n_genes, self.t_high_conf_counts.sum(axis=1))

        # nRNA per gene (EM-assigned only)
        nrna = _aggregate_to_gene(t_to_g, n_genes, self._nrna_per_transcript())
        rna_total = mrna + nrna

        # Gene effective length: abundance-weighted mean of transcript
        # effective lengths.  For genes with zero counts, use the
        # unweighted mean of transcript effective lengths.
        t_eff = self._t_eff_len
        t_counts_flat = (self.unambig_counts + self.em_counts).sum(axis=1)
        g_eff_num = _aggregate_to_gene(t_to_g, n_genes, t_counts_flat * t_eff)
        g_eff_den = _aggregate_to_gene(t_to_g, n_genes, t_counts_flat)
        g_eff_sum = _aggregate_to_gene(t_to_g, n_genes, t_eff)
        g_eff_cnt = _aggregate_to_gene(t_to_g, n_genes, np.ones_like(t_eff))
        with np.errstate(divide="ignore", invalid="ignore"):
            g_eff_mean = np.where(g_eff_cnt > 0, g_eff_sum / g_eff_cnt, 1.0)
            g_eff_len = np.where(g_eff_den > 0, g_eff_num / g_eff_den, g_eff_mean)
        g_eff_len = np.maximum(g_eff_len, 1.0)

        # Gene-level locus_id: assign the locus of the transcript with
        # the highest mRNA count.
        g_locus_id = _gene_locus_ids(
            t_to_g, n_genes, self.locus_id_per_transcript, total.sum(axis=1)
        )

        # TPM: mRNA-based at gene level
        rpk = mrna / g_eff_len
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        df = pd.DataFrame(
            {
                "gene_id": index.g_df["g_id"].values,
                "gene_name": index.g_df["g_name"].values,
                "locus_id": g_locus_id,
                "effective_length": g_eff_len,
                "mrna": mrna,
                "mrna_unambig": mrna_unambig,
                "mrna_em": mrna_em_arr,
                "mrna_high_conf": mrna_hc_arr,
                "mrna_spliced": mrna_spliced,
                "nrna": nrna,
                "rna_total": rna_total,
                "tpm": tpm,
            }
        )
        return df

    # ------------------------------------------------------------------
    # Output - locus-level breakdown
    # ------------------------------------------------------------------

    def get_loci_df(self) -> pd.DataFrame:
        """Locus-level output with gDNA/nRNA/mRNA breakdown.

        Returns an empty DataFrame (with correct columns) if no loci
        were processed.

        Columns
        -------
        locus_id : int, sequential locus identifier
        ref : str, primary reference name
        n_transcripts : int
        n_genes : int
        n_em_fragments : int, ambiguous fragments entering EM
        mrna : float, mRNA count from locus EM
        nrna : float, nRNA count from locus EM
        gdna : float, gDNA count from locus EM
        gdna_rate : float, gdna / (mrna + nrna + gdna)
        gdna_init : float, calibration locus γ (gDNA fraction)
        """
        cols = [
            "locus_id",
            "ref",
            "n_transcripts",
            "n_genes",
            "n_em_fragments",
            "mrna",
            "nrna",
            "gdna",
            "gdna_rate",
            "gdna_init",
        ]
        if not self.locus_results:
            return pd.DataFrame(columns=cols)

        rows = []
        for r in self.locus_results:
            total = r["mrna"] + r["nrna"] + r["gdna"]
            rate = r["gdna"] / total if total > 0 else 0.0
            rows.append(
                {
                    "locus_id": r["locus_id"],
                    "ref": r.get("ref", ""),
                    "n_transcripts": r["n_transcripts"],
                    "n_genes": r["n_genes"],
                    "n_em_fragments": r["n_em_fragments"],
                    "mrna": r["mrna"],
                    "nrna": r["nrna"],
                    "gdna": r["gdna"],
                    "gdna_rate": rate,
                    "gdna_init": r.get("gdna_init", 0.0),
                }
            )
        return pd.DataFrame(rows, columns=cols)

    # ------------------------------------------------------------------
    # Output - detail (long format QC breakdown)
    # ------------------------------------------------------------------

    def get_detail_df(self, index) -> pd.DataFrame:
        """Detailed counts in long format for QC."""
        t_ids = index.t_df["t_id"].values
        g_ids = index.t_df["g_id"].values

        frames = []
        for source_name, counts in (
            ("unambig", self.unambig_counts),
            ("em", self.em_counts),
        ):
            cat_counts = np.zeros(
                (self.num_transcripts, len(SpliceType)),
                dtype=np.float64,
            )
            for cat in SpliceType:
                sense_col = SpliceStrandCol.from_category(cat, False)
                anti_col = SpliceStrandCol.from_category(cat, True)
                cat_counts[:, int(cat)] = counts[:, sense_col] + counts[:, anti_col]

            nz_t, nz_cat = np.nonzero(cat_counts)
            if len(nz_t) == 0:
                continue
            cat_order = [c.name.lower() for c in SpliceType]
            frames.append(
                pd.DataFrame(
                    {
                        "transcript_id": t_ids[nz_t],
                        "gene_id": g_ids[nz_t],
                        "category": pd.Categorical(
                            [SpliceType(c).name.lower() for c in nz_cat],
                            categories=cat_order,
                        ),
                        "source": source_name,
                        "count": cat_counts[nz_t, nz_cat],
                    }
                )
            )

        if not frames:
            return pd.DataFrame(
                columns=[
                    "transcript_id",
                    "gene_id",
                    "category",
                    "source",
                    "count",
                ]
            )
        return pd.concat(frames, ignore_index=True)
