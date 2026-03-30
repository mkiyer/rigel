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
from .native import batch_locus_em_partitioned as _batch_locus_em_partitioned
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
    Estimates transcript-level mRNA and gDNA abundances.

    Locus-level EM architecture::

        Per locus:
        [0, n_t)       - transcript components (mRNA + synthetic nRNA)
        [n_t]          - gDNA: ONE shadow per locus

    Fragment routing:

    - SPLICED_ANNOT + unambig -> deterministic mRNA (no EM)
    - All other -> enter locus EM
      - Spliced fragments: transcript candidates only
      - Unspliced fragments: transcript + gDNA candidates

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
        em_config: EMConfig | None = None,
        geometry: TranscriptGeometry | None = None,
        is_synthetic_nrna: np.ndarray | None = None,
    ):
        if em_config is None:
            em_config = EMConfig()
        self.em_config = em_config
        self.num_transcripts = num_transcripts
        self._rng = np.random.default_rng(em_config.seed)
        self._synthetic_mask = (
            np.asarray(is_synthetic_nrna, dtype=bool)
            if is_synthetic_nrna is not None
            else np.zeros(num_transcripts, dtype=bool)
        )

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

        # --- gDNA: locus-level, NOT per-gene ---
        # Total gDNA count as assigned by locus EM (sum across all loci)
        self._gdna_em_total = 0.0

        # Per-locus results (locus_id, ref, n_transcripts, n_genes,
        # n_em_fragments, mrna, gdna, gdna_init).
        self.locus_results: list[dict] = []

        # Per-transcript locus_id assignment (-1 = no locus).
        self.locus_id_per_transcript = np.full(num_transcripts, -1, dtype=np.int32)

        # Per-transcript gDNA locus attribution (for reporting).
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_SPLICE_STRAND_COLS), dtype=np.float64
        )

        # Per-transcript confidence tracking:
        self._em_posterior_sum: np.ndarray | None = None
        self._em_n_assigned: np.ndarray | None = None

        # Per-locus profiling stats (populated if emit_locus_stats=True)
        self.locus_stats: list[dict] | None = None

    @property
    def effective_lengths(self) -> np.ndarray:
        """Per-transcript effective lengths (read-only)."""
        return self._t_eff_len

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA count across all loci."""
        return self._gdna_em_total

    @property
    def gdna_contamination_rate(self) -> float:
        """Fraction of all fragments attributed to gDNA."""
        total_rna = float(self.unambig_counts.sum() + self.em_counts.sum())
        total = total_rna + self._gdna_em_total
        if total == 0:
            return 0.0
        return self._gdna_em_total / total

    @property
    def nrna_em_count(self) -> float:
        """Total nRNA count (sum of EM + unambig counts for synthetic nRNA transcripts)."""
        mask = self._synthetic_mask
        if not mask.any():
            return 0.0
        return float(self.t_counts[mask].sum())

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unambig + EM), shape (N_t, 8)."""
        return self.unambig_counts + self.em_counts

    # ------------------------------------------------------------------
    # Batch locus EM (partition-native)
    # ------------------------------------------------------------------

    def run_batch_locus_em_partitioned(
        self,
        partition_tuples: list,
        locus_transcript_indices: list,
        locus_gammas: np.ndarray,
        gdna_spans: np.ndarray,
        index,
        *,
        em_iterations: int = 1000,
        em_convergence_delta: float = 1e-6,
        emit_locus_stats: bool = False,
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """Run locus-level EM from pre-partitioned data.

        Parameters
        ----------
        partition_tuples : list[tuple]
            List of 12-tuples, one per locus, containing partition arrays.
        locus_transcript_indices : list[np.ndarray]
            List of int32 transcript index arrays, one per locus.
        locus_gammas : np.ndarray
            float64 — calibration gDNA fraction per locus being processed.
        gdna_spans : np.ndarray
            int64 — pre-computed union genomic footprints per locus.
        index : TranscriptIndex
            Reference index.
        em_iterations, em_convergence_delta
            EM algorithm parameters.
        emit_locus_stats : bool
            If True, collect per-locus profiling stats.

        Returns
        -------
        tuple[float, np.ndarray, np.ndarray]
            (total_gdna_em, locus_mrna, locus_gdna)
        """
        n_transcripts = self.num_transcripts
        t_lengths = index.t_df["length"].values.astype(np.int64)

        if self._em_posterior_sum is None:
            self._em_posterior_sum = np.zeros(n_transcripts, dtype=np.float64)
            self._em_n_assigned = np.zeros(n_transcripts, dtype=np.float64)

        _ASSIGNMENT_MODE_MAP = {"fractional": 0, "map": 1, "sample": 2}
        assignment_mode_int = _ASSIGNMENT_MODE_MAP[self.em_config.assignment_mode]
        rng_seed = int(self._rng.integers(0, 2**63))

        total_gdna_em, locus_mrna, locus_gdna, locus_stats_raw = (
            _batch_locus_em_partitioned(
                partition_tuples,
                locus_transcript_indices,
                np.ascontiguousarray(locus_gammas, dtype=np.float64),
                np.ascontiguousarray(gdna_spans, dtype=np.int64),
                self.unambig_counts,
                t_lengths,
                self.em_counts,
                self.gdna_locus_counts,
                self._em_posterior_sum,
                self._em_n_assigned,
                em_iterations,
                em_convergence_delta,
                self.em_config.prior_pseudocount,
                self.em_config.mode == "vbem",
                assignment_mode_int,
                self.em_config.assignment_min_posterior,
                rng_seed,
                n_transcripts,
                NUM_SPLICE_STRAND_COLS,
                self.em_config.n_threads,
                emit_locus_stats,
            )
        )

        if emit_locus_stats:
            if self.locus_stats is None:
                self.locus_stats = []
            self.locus_stats.extend(locus_stats_raw)

        return (
            total_gdna_em,
            np.asarray(locus_mrna),
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

    # ------------------------------------------------------------------
    # Output - primary counts
    # ------------------------------------------------------------------

    def get_counts_df(self, index) -> pd.DataFrame:
        """Primary transcript-level abundance estimates (annotated only).

        Columns
        -------
        transcript_id, gene_id, gene_name : identifiers
        effective_length : bias-corrected effective length
        locus_id : int32, EM locus (-1 if no locus)
        nrna_id : str, associated nRNA entity transcript ID ("." if none)
        is_basic, is_mane, is_nascent_equiv : flags
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_spliced : spliced mRNA fragments
        tpm : transcripts per million (mRNA-based)
        posterior_mean : mean posterior at EM assignment
        """
        total = self.t_counts
        unique = self.unambig_counts

        t_total = total.sum(axis=1)
        mrna_unambig = unique.sum(axis=1)
        mrna_spliced = total[:, list(SPLICED_COLS)].sum(axis=1)
        mrna_em = self.em_counts.sum(axis=1)

        # nRNA: for synthetic nRNA transcripts, their transcript count IS the
        # nRNA count; their mrna count is 0. For annotated transcripts, nrna = 0.
        is_synthetic = index.t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in index.t_df.columns else np.zeros(self.num_transcripts, dtype=bool)
        mrna = np.where(is_synthetic, 0.0, t_total)
        pmean = self.posterior_mean()

        # TPM: mRNA-based, using effective lengths
        eff = self._t_eff_len
        rpk = mrna / eff
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        # Build nrna_id: lookup nrna_t_index → transcript ID
        t_ids = index.t_df["t_id"].values
        nrna_t_idx = index.t_df["nrna_t_index"].values if "nrna_t_index" in index.t_df.columns else np.full(len(t_ids), -1, dtype=int)
        nrna_id = np.where(
            nrna_t_idx >= 0,
            t_ids[nrna_t_idx.clip(0)],
            ".",
        )
        # Fix entries where nrna_t_index == -1 (clip(0) forced index 0)
        nrna_id = np.where(nrna_t_idx >= 0, nrna_id, ".")

        df = pd.DataFrame(
            {
                "transcript_id": t_ids,
                "gene_id": index.t_df["g_id"].values,
                "gene_name": index.t_df["g_name"].values,
                "effective_length": eff,
                "locus_id": self.locus_id_per_transcript,
                "nrna_id": nrna_id,
                "is_basic": index.t_df["is_basic"].values,
                "is_mane": index.t_df["is_mane"].values,
                "is_nascent_equiv": index.t_df["is_nascent_equiv"].values,
                "mrna": mrna,
                "mrna_unambig": mrna_unambig,
                "mrna_em": mrna_em,
                "mrna_spliced": mrna_spliced,
                "tpm": tpm,
                "posterior_mean": pmean,
            }
        )
        # Filter out synthetic nRNA transcripts (they go to nrna_quant)
        df = df[~is_synthetic].reset_index(drop=True)
        # Sort by transcript_id for deterministic output
        df = df.sort_values("transcript_id").reset_index(drop=True)
        return df

    def get_gene_counts_df(self, index) -> pd.DataFrame:
        """Primary gene-level abundance estimates.

        Columns
        -------
        gene_id, gene_name : identifiers
        n_transcripts : int, number of annotated transcripts
        locus_id : int32, primary EM locus for this gene (-1 if none)
        effective_length : abundance-weighted mean effective length
        mrna : total mRNA count (unambig + EM)
        mrna_unambig : uniquely-assigned mRNA
        mrna_em : EM-assigned mRNA
        mrna_spliced : spliced mRNA fragments
        tpm : transcripts per million (mRNA-based)
        """
        t_to_g = index.t_to_g_arr
        total = self.t_counts
        unique = self.unambig_counts
        n_genes = index.num_genes

        g_total = _aggregate_to_gene(t_to_g, n_genes, total)
        g_unambig = _aggregate_to_gene(t_to_g, n_genes, unique)

        g_all = g_total.sum(axis=1)
        mrna_unambig = g_unambig.sum(axis=1)
        mrna_spliced = g_total[:, list(SPLICED_COLS)].sum(axis=1)

        mrna_em_arr = _aggregate_to_gene(t_to_g, n_genes, self.em_counts.sum(axis=1))

        # nRNA per gene: sum counts of synthetic nRNA transcripts in each gene
        is_synthetic = index.t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in index.t_df.columns else np.zeros(len(t_to_g), dtype=bool)
        t_counts_all = (self.unambig_counts + self.em_counts).sum(axis=1)
        t_nrna = np.where(is_synthetic, t_counts_all, 0.0)
        nrna = _aggregate_to_gene(t_to_g, n_genes, t_nrna)
        mrna = g_all - nrna

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

        # n_transcripts (annotated only, excluding synthetics)
        annotated_mask = (~is_synthetic).astype(np.float64)
        n_annotated = _aggregate_to_gene(t_to_g, n_genes, annotated_mask).astype(int)

        # TPM: mRNA-based at gene level
        rpk = mrna / g_eff_len
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        df = pd.DataFrame(
            {
                "gene_id": index.g_df["g_id"].values,
                "gene_name": index.g_df["g_name"].values,
                "n_transcripts": n_annotated,
                "locus_id": g_locus_id,
                "effective_length": g_eff_len,
                "mrna": mrna,
                "mrna_unambig": mrna_unambig,
                "mrna_em": mrna_em_arr,
                "mrna_spliced": mrna_spliced,
                "tpm": tpm,
            }
        )
        # Sort by gene_id for deterministic output
        df = df.sort_values("gene_id").reset_index(drop=True)
        return df

    # ------------------------------------------------------------------
    # Output - nRNA entity table
    # ------------------------------------------------------------------

    def get_nrna_counts_df(self, index) -> pd.DataFrame:
        """Nascent RNA entity table (synthetics + nascent equivalents).

        Columns
        -------
        nrna_id : str, transcript ID for the nRNA entity
        effective_length : bias-corrected effective length
        locus_id : int32, EM locus
        is_synthetic : bool, True for RIGEL-generated synthetics
        n_contributing_transcripts : int, annotated transcripts merged
        count : float, nascent RNA fragment count
        tpm : float, TPM (nRNA-pool normalized)
        """
        t_df = index.t_df
        is_syn = t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in t_df.columns else np.zeros(self.num_transcripts, dtype=bool)
        is_equiv = t_df["is_nascent_equiv"].values if "is_nascent_equiv" in t_df.columns else np.zeros(self.num_transcripts, dtype=bool)
        mask = is_syn | is_equiv

        if not mask.any():
            return pd.DataFrame(columns=[
                "nrna_id", "effective_length",
                "locus_id", "is_synthetic", "n_contributing_transcripts",
                "count", "tpm",
            ])

        idx = np.where(mask)[0]
        t_total = self.t_counts.sum(axis=1)

        # For synthetics, count = their total. For nascent-equiv, count = 0
        # (their fragments are counted as mRNA in quant.tsv).
        counts = np.where(is_syn[mask], t_total[idx], 0.0)

        eff = self._t_eff_len[idx]

        # nRNA-pool TPM
        rpk = counts / np.maximum(eff, 1.0)
        rpk_sum = rpk.sum()
        tpm = (rpk / rpk_sum * 1e6) if rpk_sum > 0 else np.zeros_like(rpk)

        n_contrib = t_df["nrna_n_contributors"].values[idx] if "nrna_n_contributors" in t_df.columns else np.zeros(len(idx), dtype=int)

        df = pd.DataFrame(
            {
                "nrna_id": t_df["t_id"].values[idx],
                "effective_length": eff,
                "locus_id": self.locus_id_per_transcript[idx],
                "is_synthetic": is_syn[mask],
                "n_contributing_transcripts": n_contrib,
                "count": counts,
                "tpm": tpm,
            }
        )
        df = df.sort_values("nrna_id").reset_index(drop=True)
        return df

    # ------------------------------------------------------------------
    # Output - locus-level breakdown
    # ------------------------------------------------------------------

    def get_loci_df(self, index=None) -> pd.DataFrame:
        """Locus-level output with three-pool breakdown.

        Returns an empty DataFrame (with correct columns) if no loci
        were processed.

        Columns
        -------
        locus_id : int, sequential locus identifier
        locus_span_bp : int, merged genomic footprint in base pairs
        n_transcripts : int, total transcripts in locus
        n_annotated_transcripts : int, annotated transcripts only
        n_nrna_entities : int, synthetic nRNA transcripts in locus
        n_genes : int
        n_em_fragments : int, ambiguous fragments entering EM
        mrna : float, mRNA count from locus EM
        nrna : float, nRNA count from locus EM
        gdna : float, gDNA count from locus EM
        total : float, mrna + nrna + gdna
        gdna_rate : float, gdna / total
        gdna_prior : float, calibration-derived gDNA prior (γ)
        """
        cols = [
            "locus_id",
            "locus_span_bp",
            "n_transcripts",
            "n_annotated_transcripts",
            "n_nrna_entities",
            "n_genes",
            "n_em_fragments",
            "mrna",
            "nrna",
            "gdna",
            "total",
            "gdna_rate",
            "gdna_prior",
        ]
        if not self.locus_results:
            return pd.DataFrame(columns=cols)

        # Pre-compute per-transcript nRNA counts and synthetic mask
        t_total = self.t_counts.sum(axis=1)
        if index is not None:
            is_syn = index.t_df["is_synthetic_nrna"].values if "is_synthetic_nrna" in index.t_df.columns else np.zeros(self.num_transcripts, dtype=bool)
        else:
            is_syn = np.zeros(self.num_transcripts, dtype=bool)
        locus_ids = self.locus_id_per_transcript

        rows = []
        for r in self.locus_results:
            lid = r["locus_id"]
            in_locus = locus_ids == lid
            syn_in_locus = in_locus & is_syn
            nrna = float(t_total[syn_in_locus].sum())
            n_syn = int(syn_in_locus.sum())
            n_annot = int(in_locus.sum()) - n_syn

            total = r["mrna"] + nrna + r["gdna"]
            rate = r["gdna"] / total if total > 0 else 0.0
            rows.append(
                {
                    "locus_id": lid,
                    "locus_span_bp": r.get("locus_span_bp", 0),
                    "n_transcripts": r["n_transcripts"],
                    "n_annotated_transcripts": n_annot,
                    "n_nrna_entities": n_syn,
                    "n_genes": r["n_genes"],
                    "n_em_fragments": r["n_em_fragments"],
                    "mrna": r["mrna"],
                    "nrna": nrna,
                    "gdna": r["gdna"],
                    "total": total,
                    "gdna_rate": rate,
                    "gdna_prior": r.get("gdna_init", 0.0),
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
