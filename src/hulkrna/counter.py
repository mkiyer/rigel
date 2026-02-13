"""
hulkrna.counter — Read count accumulation with EM-based Bayesian assignment.

``ReadCounter`` accumulates discrete read counts into transcript
arrays.  Each fragment contributes exactly **one count** (1.0) to a
single transcript, preserving the biological interpretation of one
sequenced molecule = one count.

For uniquely-mapping fragments (single gene, single transcript,
NH = 1), assignment is deterministic.  For ambiguous fragments —
whether isoform-ambiguous (single gene, multiple transcripts),
gene-ambiguous (multiple genes), or multimapped (NH > 1) — transcript
abundances are first estimated via Expectation-Maximization (EM), then
each fragment is assigned stochastically using the converged
posteriors.

This unified EM approach eliminates order-dependent bias between
ambiguity classes.  All ambiguous fragments compete for transcript
probability mass simultaneously, producing consistent results
regardless of the relative proportions of isoform-ambiguous,
gene-ambiguous, and multimapped fragments.  This is the standard
approach used by RSEM, Salmon, and Kallisto.

Hierarchical per-gene gDNA shadow model
-----------------------------------------
Each gene has its own gDNA "shadow" component that competes with
transcripts during EM.  Shadow priors are initialized via empirical
Bayes shrinkage using the "Iceberg" model:

    θ_global = (n_intergenic + 2 × Σ antisense) / n_total
    shadow_init_g = 2 × antisense_g + κ × depth_g × θ_global

Genomic DNA distributes equally on both strands, so observed
antisense reads are the tip of the iceberg — an equal mass of
sense gDNA hides among RNA reads.  The global rate θ_global uses
*all* gDNA evidence (intergenic + 2× antisense) rather than
intergenic alone, which would underestimate the true rate in modern
annotations (GENCODE v44+) where intergenic space is minimal.

The per-gene shrinkage term ``κ × depth_g × θ_global`` scales by
gene depth, so deeply-sequenced genes get proportionally larger
gDNA floors than lowly-sequenced genes.

EM algorithm
------------
1. **Initialize** transcript abundance estimates (θ) from unique
   counts plus a Dirichlet pseudocount (α).  Shadow components
   are initialized from ``shadow_init``.
2. **E-step**: For each ambiguous fragment, compute the posterior
   probability of each candidate transcript or shadow:
   ``P(t | frag) ∝ θ_t × P_strand(t) × P_insert(frag)``
3. **M-step**: Update θ by accumulating fractional counts.
4. **Repeat** until convergence or iteration limit.
5. **Assign**: Use converged posteriors for final stochastic
   assignment (one discrete count per fragment).

All vectorized operations use segment-wise NumPy ``reduceat`` for
efficient batched computation across variable-length candidate sets
stored in CSR (Compressed Sparse Row) format.
"""

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .types import Strand
from .categories import (
    AssignmentSource,
    CountCategory,
    CountStrand,
    CountType,
    NUM_COUNT_TYPES,
)

logger = logging.getLogger(__name__)


# Pre-computed lookup tables for sparse output
_CT_CATEGORY_NAMES = [
    CountCategory(i // 3).name.lower() for i in range(NUM_COUNT_TYPES)
]
_CT_STRAND_NAMES = [
    CountStrand(i % 3).name.lower() for i in range(NUM_COUNT_TYPES)
]
_CATEGORY_ORDER = [c.name.lower() for c in CountCategory]
_STRAND_ORDER = [s.name.lower() for s in CountStrand]
_SOURCE_ORDER = [s.name.lower() for s in AssignmentSource]


# ======================================================================
# EMData — pre-computed CSR arrays for vectorized EM
# ======================================================================


# Sentinel transcript index for the gDNA pseudo-component in EM.
# Set dynamically by ReadCounter based on num_transcripts.
GDNA_INDEX: int = -1  # Placeholder; overridden per-counter instance


@dataclass(slots=True)
class EMData:
    """Pre-computed CSR arrays for vectorized EM.

    Stores per-candidate log-likelihoods and count type indices in
    Compressed Sparse Row format.  Each "unit" is either a single
    ambiguous fragment or a multimapper molecule (all alignments of
    one molecule grouped together).

    Per-gene gDNA shadow components are included as extra candidates
    (with ``t_indices == gdna_base_index + gene_index``) for
    non-SPLICED_ANNOT fragments.  Each gene has its own shadow index
    in the range ``[gdna_base_index, gdna_base_index + num_genes)``.

    ``locus_t_indices`` and ``locus_count_types`` track the best
    transcript candidate per unit for gDNA locus attribution.

    Attributes
    ----------
    offsets : np.ndarray
        int64[n_units + 1] — CSR offsets into flat arrays.
    t_indices : np.ndarray
        int32[n_candidates] — candidate transcript or shadow indices.
        Shadow indices are ``gdna_base_index + gene_index``.
    log_liks : np.ndarray
        float64[n_candidates] — log(P_strand × P_insert) per candidate.
        These are the **data likelihoods** (fixed across EM iterations).
    count_types : np.ndarray
        uint8[n_candidates] — CountType column index per candidate.
    locus_t_indices : np.ndarray
        int32[n_units] — best transcript index per unit (highest
        log_lik among non-shadow candidates).  Used to attribute
        gDNA assignments to their genomic locus.  -1 if no
        transcript candidate exists.
    locus_count_types : np.ndarray
        uint8[n_units] — count type for the locus transcript.
    n_units : int
        Number of ambiguous units (fragments + multimapper molecules).
    n_candidates : int
        Total number of (unit, candidate) entries across all units.
    gdna_base_index : int
        First shadow index.  Shadow indices occupy the range
        ``[gdna_base_index, gdna_base_index + num_genes)``.
    """

    offsets: np.ndarray
    t_indices: np.ndarray
    log_liks: np.ndarray
    count_types: np.ndarray
    locus_t_indices: np.ndarray
    locus_count_types: np.ndarray
    n_units: int
    n_candidates: int
    gdna_base_index: int


# ======================================================================
# ReadCounter
# ======================================================================


class ReadCounter:
    """Accumulates discrete read counts into transcript arrays.

    Maintains two separate count arrays for provenance tracking:

    - ``unique_counts``: deterministic assignments from unambiguous
      fragments (1 gene, 1 transcript, NH = 1).
    - ``em_counts``: stochastic assignments from EM-resolved
      ambiguous fragments.

    Count arrays have shape ``(N, 12)`` where the 12 columns
    correspond to ``CountType`` values (4 categories × 3 strands).

    Hierarchical per-gene gDNA shadows
    -----------------------------------
    Each gene has its own shadow component in the θ vector at
    positions ``[gdna_base_index, gdna_base_index + num_genes)``.
    Shadow priors (``shadow_init``) are computed externally via
    empirical Bayes shrinkage and passed in at construction time.

    Parameters
    ----------
    num_transcripts : int
        Number of transcripts in the index.
    num_genes : int
        Number of genes in the index.
    seed : int or None
        Random seed for reproducibility (default None = non-deterministic).
    alpha : float
        Dirichlet pseudocount added to transcript counts before
        computing EM priors (default 1.0).
    shadow_init : np.ndarray or None
        Per-gene shadow initialization values, shape ``(num_genes,)``.
        If None, all shadows start with 0 (equivalent to no gDNA prior
        beyond the Dirichlet alpha).
    """

    def __init__(
        self,
        num_transcripts: int,
        num_genes: int,
        *,
        seed=None,
        alpha: float = 1.0,
        shadow_init: np.ndarray | None = None,
    ):
        self.num_transcripts = num_transcripts
        self.num_genes = num_genes
        self.alpha = alpha
        self.rng = np.random.default_rng(seed)

        # Base index for per-gene shadow components in EM theta vector.
        # Shadows occupy indices [gdna_base_index, gdna_base_index + num_genes).
        self.gdna_base_index = num_transcripts

        self.unique_counts = np.zeros(
            (num_transcripts, NUM_COUNT_TYPES), dtype=np.float64
        )
        self.em_counts = np.zeros(
            (num_transcripts, NUM_COUNT_TYPES), dtype=np.float64
        )

        # --- Per-gene gDNA shadow tracking ---
        # shadow_init[g] encodes the EB-shrunk prior for gene g's shadow.
        if shadow_init is not None:
            self.shadow_init = np.asarray(shadow_init, dtype=np.float64)
        else:
            self.shadow_init = np.zeros(num_genes, dtype=np.float64)

        # Per-gene EM-assigned gDNA counts.
        self.gdna_em_counts = np.zeros(num_genes, dtype=np.float64)

        # Per-transcript locus attribution of gDNA assignments.
        # When a fragment is assigned to a gene's shadow by EM, the
        # transcript with the highest data likelihood is recorded here.
        self.gdna_locus_counts = np.zeros(
            (num_transcripts, NUM_COUNT_TYPES), dtype=np.float64
        )

        # Converged transcript abundance from EM (set by run_em)
        self._converged_theta: np.ndarray | None = None

    # --- Backward-compatible scalar properties ---

    @property
    def gdna_index(self) -> int:
        """Backward-compatible alias for gdna_base_index."""
        return self.gdna_base_index

    @property
    def gdna_unique_count(self) -> float:
        """Total shadow_init (backward-compatible scalar)."""
        return float(self.shadow_init.sum())

    @property
    def gdna_em_count(self) -> float:
        """Total EM-assigned gDNA (backward-compatible scalar)."""
        return float(self.gdna_em_counts.sum())

    @property
    def gdna_total(self) -> float:
        """Total gDNA count (shadow_init + EM-assigned)."""
        return self.gdna_unique_count + self.gdna_em_count

    @property
    def gdna_contamination_rate(self) -> float:
        """Fraction of all fragments attributed to gDNA.

        Returns 0.0 if no fragments were counted.
        """
        total_rna = float(self.unique_counts.sum() + self.em_counts.sum())
        total = total_rna + self.gdna_total
        if total == 0:
            return 0.0
        return self.gdna_total / total

    @property
    def t_counts(self) -> np.ndarray:
        """Total transcript counts (unique + EM), shape ``(N_t, 12)``."""
        return self.unique_counts + self.em_counts

    # ------------------------------------------------------------------
    # Strand classification
    # ------------------------------------------------------------------

    @staticmethod
    def classify_strand(
        exon_strand,
        gene_strand: int,
        strand_model,
    ) -> CountStrand:
        """Classify the strand relationship using the trained model.

        Uses the strand model's learned probability to determine
        sense/antisense/ambiguous.  For highly stranded libraries
        (specificity > 0.8), fragments matching the gene strand
        are SENSE and mismatching are ANTISENSE.  For weakly stranded
        or unstranded libraries, most fragments become AMBIGUOUS.

        Parameters
        ----------
        exon_strand : Strand
            Combined alignment strand of the fragment.
        gene_strand : int
            Annotated strand of the candidate gene.
        strand_model : StrandModel
            Trained strand model from the BAM scan.

        Returns
        -------
        CountStrand
        """
        p = strand_model.strand_likelihood(exon_strand, Strand(gene_strand))
        if p > 0.8:
            return CountStrand.SENSE
        elif p < 0.2:
            return CountStrand.ANTISENSE
        else:
            return CountStrand.AMBIGUOUS

    # ------------------------------------------------------------------
    # Unique assignment
    # ------------------------------------------------------------------

    def assign_unique(self, resolved, index, strand_model) -> None:
        """Deterministic assignment for a truly unique fragment.

        One gene, one transcript, NH = 1 → count += 1 to
        ``unique_counts``.

        Parameters
        ----------
        resolved : ResolvedFragment or BufferedFragment
            Resolution output — must be single-gene, single-transcript.
        index : HulkIndex
            The loaded reference index.
        strand_model : StrandModel
            Trained strand model.
        """
        t_inds = resolved.t_inds

        if len(t_inds) == 0:
            return

        t_idx = int(next(iter(t_inds)))
        g_idx = int(index.t_to_g_arr[t_idx])
        gene_strand = int(index.g_to_strand_arr[g_idx])
        count_strand = self.classify_strand(
            resolved.exon_strand, gene_strand, strand_model
        )
        count_type = int(resolved.count_cat) * 3 + int(count_strand)
        self.unique_counts[t_idx, count_type] += 1.0

    # ------------------------------------------------------------------
    # EM iteration (vectorized)
    # ------------------------------------------------------------------

    def run_em(self, em_data: EMData, *, em_iterations: int = 10) -> None:
        """Run vectorized EM to estimate transcript + per-gene shadow abundances.

        The abundance vector θ has ``num_transcripts + num_genes``
        elements.  Positions ``[0, num_transcripts)`` are transcripts;
        positions ``[gdna_base_index, gdna_base_index + num_genes)``
        are per-gene gDNA shadows.

        Starting from unique-count-derived priors (transcripts) and
        shadow_init-derived priors (per-gene shadows) with Dirichlet
        smoothing, iteratively refines abundance estimates using all
        ambiguous fragments simultaneously.

        Parameters
        ----------
        em_data : EMData
            Pre-computed CSR arrays from the pipeline's buffer scan.
        em_iterations : int
            Maximum number of EM iterations (default 10).
            Iteration stops early if the total absolute change in
            θ drops below 1e-6.
        """
        n_total = self.num_transcripts + self.num_genes
        gdna_base = self.gdna_base_index

        # Initialize from unique counts + pseudocount
        unique_totals = np.zeros(n_total, dtype=np.float64)
        unique_totals[:self.num_transcripts] = self.unique_counts.sum(axis=1)
        unique_totals[gdna_base:gdna_base + self.num_genes] = self.shadow_init

        theta = unique_totals + self.alpha
        theta /= theta.sum()

        if em_data.n_units == 0 or em_data.n_candidates == 0:
            self._converged_theta = theta
            return

        offsets = em_data.offsets
        t_indices = em_data.t_indices
        log_liks = em_data.log_liks
        seg_lengths = np.diff(offsets)

        for it in range(em_iterations):
            # E-step: posterior ∝ likelihood × prior
            log_posteriors = log_liks + np.log(theta[t_indices] + 1e-300)

            # Segment-wise log-sum-exp normalization
            seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
            log_posteriors -= np.repeat(seg_max, seg_lengths)
            posteriors = np.exp(log_posteriors)
            seg_sum = np.add.reduceat(posteriors, offsets[:-1])
            posteriors /= np.repeat(seg_sum, seg_lengths)

            # M-step: accumulate fractional counts per component
            em_totals = np.zeros(n_total, dtype=np.float64)
            np.add.at(em_totals, t_indices, posteriors)

            # Update theta
            theta_new = unique_totals + em_totals + self.alpha
            theta_new /= theta_new.sum()

            delta = np.abs(theta_new - theta).sum()
            theta = theta_new

            logger.debug(f"  EM iteration {it + 1}: delta={delta:.6f}")

            if delta < 1e-6:
                logger.debug(
                    f"  EM converged after {it + 1} iterations"
                )
                break

        self._converged_theta = theta
        gdna_theta_sum = theta[gdna_base:gdna_base + self.num_genes].sum()
        shadow_init_sum = self.shadow_init.sum()
        logger.info(
            f"  gDNA theta (sum) = {gdna_theta_sum:.6f}, "
            f"shadow_init (sum) = {shadow_init_sum:.0f} / "
            f"{unique_totals.sum():.0f}"
        )

    def assign_ambiguous(self, em_data: EMData) -> None:
        """Stochastic assignment of ambiguous fragments using converged posteriors.

        Each ambiguous unit contributes exactly one count: either to a
        transcript (written to ``em_counts``) or to a per-gene gDNA
        shadow (written to ``gdna_em_counts`` and ``gdna_locus_counts``).

        Uses vectorized cumsum + searchsorted for efficient batched
        sampling.

        Must be called after :meth:`run_em`.

        Parameters
        ----------
        em_data : EMData
            Same pre-computed data used for :meth:`run_em`.
        """
        if em_data.n_units == 0 or em_data.n_candidates == 0:
            return

        theta = self._converged_theta
        if theta is None:
            raise RuntimeError(
                "run_em() must be called before assign_ambiguous()"
            )

        offsets = em_data.offsets
        t_indices = em_data.t_indices
        log_liks = em_data.log_liks
        count_types = em_data.count_types
        seg_lengths = np.diff(offsets)
        n_units = em_data.n_units

        # Compute final posteriors using converged theta
        log_posteriors = log_liks + np.log(theta[t_indices] + 1e-300)

        # Segment-wise softmax normalization
        seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
        log_posteriors -= np.repeat(seg_max, seg_lengths)
        posteriors = np.exp(log_posteriors)
        seg_sum = np.add.reduceat(posteriors, offsets[:-1])
        posteriors /= np.repeat(seg_sum, seg_lengths)

        # --- Vectorized stochastic choice via cumsum + searchsorted ---
        global_cumsum = np.cumsum(posteriors)

        # Per-segment offset into global cumsum
        start_vals = np.zeros(n_units, dtype=np.float64)
        if n_units > 1:
            start_vals[1:] = global_cumsum[offsets[1:-1] - 1]

        # Segment totals for scaling random values
        seg_totals = global_cumsum[offsets[1:] - 1] - start_vals

        # Generate random values ∈ [0, seg_total) per unit
        u_rand = self.rng.random(n_units) * seg_totals

        # Shift to global cumsum space and searchsorted
        u_rand_global = u_rand + start_vals
        chosen_flat = np.searchsorted(
            global_cumsum, u_rand_global, side="left"
        )

        # Clip to valid range (floating-point safety)
        chosen_flat = np.clip(chosen_flat, offsets[:-1], offsets[1:] - 1)

        # Separate gDNA assignments from transcript assignments
        chosen_t = t_indices[chosen_flat]
        chosen_ct = count_types[chosen_flat]

        gdna_base = self.gdna_base_index
        gdna_mask = chosen_t >= gdna_base
        rna_mask = ~gdna_mask

        # --- Transcript assignments → em_counts ---
        if rna_mask.any():
            np.add.at(
                self.em_counts,
                (chosen_t[rna_mask], chosen_ct[rna_mask]),
                1.0,
            )

        # --- gDNA assignments → per-gene gdna_em_counts + gdna_locus_counts ---
        if gdna_mask.any():
            gene_indices = chosen_t[gdna_mask] - gdna_base
            np.add.at(self.gdna_em_counts, gene_indices, 1.0)
            # Attribute gDNA assignments to locus transcripts
            locus_t = em_data.locus_t_indices[gdna_mask]
            locus_ct = em_data.locus_count_types[gdna_mask]
            # Filter out units with no locus transcript (locus_t == -1)
            valid = locus_t >= 0
            if valid.any():
                np.add.at(
                    self.gdna_locus_counts,
                    (locus_t[valid], locus_ct[valid]),
                    1.0,
                )

    # ------------------------------------------------------------------
    # Output — wide format
    # ------------------------------------------------------------------

    def get_t_counts_df(self) -> pd.DataFrame:
        """Return total transcript counts as a wide (N_t × 12) DataFrame."""
        return pd.DataFrame(
            self.t_counts, columns=list(CountType.columns())
        )

    def get_g_counts_df(self, t_to_g_arr: np.ndarray) -> pd.DataFrame:
        """Return total gene counts as a wide (N_g × 12) DataFrame.

        Gene counts are derived by summing transcript counts per gene.

        Parameters
        ----------
        t_to_g_arr : np.ndarray
            Mapping array ``t_index → g_index``.
        """
        g_counts = np.zeros(
            (self.num_genes, NUM_COUNT_TYPES), dtype=np.float64
        )
        np.add.at(g_counts, t_to_g_arr, self.t_counts)
        return pd.DataFrame(
            g_counts, columns=list(CountType.columns())
        )

    # ------------------------------------------------------------------
    # Output — sparse long format
    # ------------------------------------------------------------------

    def get_sparse_t_counts_df(self) -> pd.DataFrame:
        """Return transcript counts in sparse long format.

        Each row represents a nonzero ``(transcript, category, strand,
        source)`` combination with an integer count.

        Returns
        -------
        pd.DataFrame
            Columns: ``transcript_idx``, ``category``, ``strand``,
            ``source``, ``count``.
        """
        frames = []
        for source, counts in (
            (AssignmentSource.UNIQUE, self.unique_counts),
            (AssignmentSource.EM, self.em_counts),
        ):
            nz_t, nz_ct = np.nonzero(counts)
            if len(nz_t) == 0:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "transcript_idx": nz_t.astype(np.int32),
                        "category": pd.Categorical(
                            [_CT_CATEGORY_NAMES[c] for c in nz_ct],
                            categories=_CATEGORY_ORDER,
                        ),
                        "strand": pd.Categorical(
                            [_CT_STRAND_NAMES[c] for c in nz_ct],
                            categories=_STRAND_ORDER,
                        ),
                        "source": source.name.lower(),
                        "count": counts[nz_t, nz_ct].astype(np.int64),
                    }
                )
            )

        if not frames:
            return pd.DataFrame(
                columns=[
                    "transcript_idx",
                    "category",
                    "strand",
                    "source",
                    "count",
                ]
            )
        return pd.concat(frames, ignore_index=True)

    # ------------------------------------------------------------------
    # Output — gDNA locus breakdown
    # ------------------------------------------------------------------

    def get_t_gdna_df(self) -> pd.DataFrame:
        """Transcript-level gDNA locus counts, shape ``(N_t, 12)``.

        When a fragment is assigned to gDNA by EM, the transcript with
        the highest data likelihood is recorded here.  This shows how
        many gDNA fragments originated from each transcript's genomic
        locus, broken down by category × strand.
        """
        return pd.DataFrame(
            self.gdna_locus_counts,
            columns=list(CountType.columns()),
        )

    def get_g_gdna_df(self, t_to_g_arr: np.ndarray) -> pd.DataFrame:
        """Gene-level gDNA locus counts, shape ``(N_g, 12)``.

        Aggregates ``gdna_locus_counts`` by gene.
        """
        g_gdna = np.zeros(
            (self.num_genes, NUM_COUNT_TYPES), dtype=np.float64
        )
        np.add.at(g_gdna, t_to_g_arr, self.gdna_locus_counts)
        return pd.DataFrame(
            g_gdna, columns=list(CountType.columns()),
        )

    def gdna_summary(self) -> dict:
        """Summary dict of gDNA contamination estimates."""
        return {
            "gdna_unique_count": float(self.gdna_unique_count),
            "gdna_em_count": float(self.gdna_em_count),
            "gdna_total": float(self.gdna_total),
            "gdna_contamination_rate": float(
                round(self.gdna_contamination_rate, 6)
            ),
            "rna_unique_total": float(self.unique_counts.sum()),
            "rna_em_total": float(self.em_counts.sum()),
        }

    def get_sparse_g_counts_df(
        self, t_to_g_arr: np.ndarray
    ) -> pd.DataFrame:
        """Return gene counts in sparse long format.

        Each row represents a nonzero ``(gene, category, strand,
        source)`` combination with an integer count.

        Parameters
        ----------
        t_to_g_arr : np.ndarray
            Mapping array ``t_index → g_index``.

        Returns
        -------
        pd.DataFrame
            Columns: ``gene_idx``, ``category``, ``strand``,
            ``source``, ``count``.
        """
        frames = []
        for source, t_counts in (
            (AssignmentSource.UNIQUE, self.unique_counts),
            (AssignmentSource.EM, self.em_counts),
        ):
            g_counts = np.zeros(
                (self.num_genes, NUM_COUNT_TYPES), dtype=np.float64
            )
            np.add.at(g_counts, t_to_g_arr, t_counts)

            nz_g, nz_ct = np.nonzero(g_counts)
            if len(nz_g) == 0:
                continue
            frames.append(
                pd.DataFrame(
                    {
                        "gene_idx": nz_g.astype(np.int32),
                        "category": pd.Categorical(
                            [_CT_CATEGORY_NAMES[c] for c in nz_ct],
                            categories=_CATEGORY_ORDER,
                        ),
                        "strand": pd.Categorical(
                            [_CT_STRAND_NAMES[c] for c in nz_ct],
                            categories=_STRAND_ORDER,
                        ),
                        "source": source.name.lower(),
                        "count": g_counts[nz_g, nz_ct].astype(np.int64),
                    }
                )
            )

        if not frames:
            return pd.DataFrame(
                columns=[
                    "gene_idx",
                    "category",
                    "strand",
                    "source",
                    "count",
                ]
            )
        return pd.concat(frames, ignore_index=True)
