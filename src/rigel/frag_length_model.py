"""
rigel.frag_length_model — Fragment length distribution model.

Learns the fragment length distribution from fragment-to-transcript
mappings produced by Pass 1.  Fragments where all candidate
transcripts yield the *same* fragment length contribute to training;
ambiguous fragments are deferred to Bayesian quantification where the
learned distribution provides a likelihood term.

The distribution is stored as a histogram (float vector) indexed
by fragment length in bases.  Sizes >= ``max_size`` are clamped into
a single overflow bin.

Typical RNA-seq fragment lengths range from 100–500 bp.  The default
``max_size=2000`` accommodates the vast majority of libraries.

Fragments whose length exceeds ``max_size`` receive exponential tail
decay: each additional base-pair beyond the boundary incurs a fixed
log-probability penalty (``_TAIL_DECAY_LP ≈ log(0.99) ≈ −0.01``
per bp), so very long fragments are penalised heavily rather than
assigned the same probability as the overflow bin.
"""

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Per-bp log-probability penalty applied to fragment lengths that exceed
# max_size.  log(0.99) ≈ −0.01005 per bp, giving:
#   100 bp over → −1.0 extra log penalty  (≈ 2.7× less likely)
#   500 bp over → −5.0                    (≈ 150× less likely)
#  1000 bp over → −10.0                   (≈ 22 000× less likely)
_TAIL_DECAY_LP: float = math.log(0.99)

#: Default maximum fragment length tracked individually.
#: Sizes >= this go into a single overflow bin.  Must match
#: ``BamScanConfig.max_frag_length`` (config.py) for consistency
#: between standalone and production usage.
DEFAULT_MAX_FRAG_SIZE: int = 1000

logger = logging.getLogger(__name__)


@dataclass
class FragmentLengthModel:
    """Fragment length distribution model.

    The histogram has ``max_size + 1`` bins:
      - ``counts[k]`` for ``0 <= k < max_size``: fragments with
        fragment length exactly *k*.
      - ``counts[max_size]``: overflow bin for fragment lengths >= max_size.

    Training
    --------
    Call :meth:`observe` with each unambiguous fragment length.  After
    training, :meth:`log_likelihood` returns the log-probability
    of a given fragment length under the learned distribution.

    Attributes
    ----------
    max_size : int
        Maximum fragment length tracked individually (sizes >= this
        go into the overflow bin).  Queries above this receive
        exponential tail decay (see ``_TAIL_DECAY_LP``).
    counts : np.ndarray
        Float64 histogram of shape ``(max_size + 1,)``.
    """

    max_size: int = DEFAULT_MAX_FRAG_SIZE
    counts: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        if self.counts is None:
            self.counts = np.zeros(self.max_size + 1, dtype=np.float64)
            self._total_weight: float = 0.0
        else:
            self._total_weight: float = float(self.counts.sum())
        self._log_prob: np.ndarray | None = None
        self._finalized: bool = False

    @property
    def n_observations(self) -> int:
        """Number of observations (derived from total histogram weight)."""
        return int(self._total_weight)

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def from_counts(
        cls,
        counts,
        max_size: int | None = None,
    ) -> "FragmentLengthModel":
        """Create a finalized model from a pre-set histogram.

        Useful for tests and scenarios that need a model with a known
        distribution without going through the observe/finalize cycle.

        Parameters
        ----------
        counts : array-like
            Histogram of fragment length frequencies.  Index *k*
            corresponds to fragment length *k*.  The last element is
            the overflow bin.  If shorter than ``max_size + 1``, the
            remaining bins are zero-filled.
        max_size : int or None
            Maximum fragment length tracked individually.  If *None*,
            inferred as ``len(counts) - 1``.

        Returns
        -------
        FragmentLengthModel
            A finalized model ready for ``log_likelihood()`` calls.
        """
        counts = np.asarray(counts, dtype=np.float64)
        if max_size is None:
            max_size = len(counts) - 1
        model = cls(max_size=max_size)
        n = min(len(counts), max_size + 1)
        model.counts[:n] = counts[:n]
        model._total_weight = float(model.counts.sum())
        model.finalize()
        return model

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def observe(self, frag_length: int, weight: float = 1.0) -> None:
        """Record one fragment length observation.

        Parameters
        ----------
        frag_length : int
            Fragment fragment length (must be >= 0).
        weight : float
            Observation weight (default 1.0).
        """
        idx = min(max(frag_length, 0), self.max_size)
        self.counts[idx] += weight
        self._total_weight += weight

    def observe_batch(self, frag_lengths: "np.ndarray") -> None:
        """Record a batch of fragment length observations (vectorized).

        Parameters
        ----------
        frag_lengths : np.ndarray
            Integer array of fragment lengths.  All weights are 1.0.
        """
        lengths = np.asarray(frag_lengths, dtype=np.intp)
        clamped = np.clip(lengths, 0, self.max_size)
        counts = np.bincount(clamped, minlength=self.max_size + 1)
        self.counts += counts[:self.max_size + 1].astype(np.float64)
        self._total_weight += float(len(lengths))

    # ------------------------------------------------------------------
    # Distribution properties
    # ------------------------------------------------------------------

    @property
    def total_weight(self) -> float:
        """Sum of all histogram weights (may differ from n_observations
        if non-unit weights are used)."""
        return self._total_weight

    @property
    def mean(self) -> float:
        """Weighted mean fragment length.

        The overflow bin is treated as having fragment length = max_size.
        Returns 0.0 if no observations.
        """
        total = self.total_weight
        if total == 0:
            return 0.0
        indices = np.arange(self.max_size + 1, dtype=np.float64)
        return float(np.dot(indices, self.counts) / total)

    @property
    def std(self) -> float:
        """Weighted standard deviation.

        Returns 0.0 if no observations.
        """
        total = self.total_weight
        if total == 0:
            return 0.0
        indices = np.arange(self.max_size + 1, dtype=np.float64)
        mu = np.dot(indices, self.counts) / total
        var = np.dot((indices - mu) ** 2, self.counts) / total
        return float(np.sqrt(var))

    @property
    def median(self) -> float:
        """Weighted median fragment length.

        Returns 0.0 if no observations.
        """
        total = self.total_weight
        if total == 0:
            return 0.0
        cumsum = np.cumsum(self.counts)
        idx = np.searchsorted(cumsum, total / 2.0)
        return float(min(idx, self.max_size))

    @property
    def mode(self) -> int:
        """Fragment length with the highest count.

        Returns 0 if no observations.
        """
        if self.total_weight == 0:
            return 0
        return int(np.argmax(self.counts))

    # ------------------------------------------------------------------
    # Finalization (call after training, before scoring)
    # ------------------------------------------------------------------

    def finalize(self) -> None:
        """Pre-compute log-likelihood lookup table for fast scoring.

        Builds ``_log_prob`` array so ``log_likelihood()`` becomes a
        single array index instead of 2× ``np.log`` per call.

        Also caches the tail-decay base value so that queries beyond
        ``max_size`` can be answered with a single multiply-add.
        """
        total = self._total_weight
        if total == 0:
            n = self.max_size + 1
            self._log_prob = np.full(n, -np.log(n), dtype=np.float64)
        else:
            self._log_prob = (
                np.log(self.counts + 1.0)
                - np.log(total + self.max_size + 1)
            )
        self._tail_base: float = float(self._log_prob[self.max_size])
        self._finalized = True

    # ------------------------------------------------------------------
    # Likelihood for Bayesian quantification
    # ------------------------------------------------------------------

    def log_likelihood(self, frag_length: int) -> float:
        """Log-probability of a fragment length under the learned distribution.

        Uses add-one (Laplace) smoothing to avoid -inf for unseen sizes.
        Fragment lengths beyond ``max_size`` receive exponential tail
        decay: ``log_prob[max_size] + (frag_length − max_size) × _TAIL_DECAY_LP``.

        Parameters
        ----------
        frag_length : int
            Query fragment length.

        Returns
        -------
        float
            Log-probability (natural log).
        """
        if frag_length > self.max_size:
            # Exponential tail decay beyond the histogram
            if self._finalized:
                base = self._tail_base
            else:
                total = self._total_weight
                if total == 0:
                    base = -math.log(self.max_size + 1)
                else:
                    base = float(
                        np.log(self.counts[self.max_size] + 1.0)
                        - np.log(total + self.max_size + 1)
                    )
            return base + (frag_length - self.max_size) * _TAIL_DECAY_LP

        idx = max(frag_length, 0)
        if self._finalized:
            return float(self._log_prob[idx])
        total = self._total_weight
        if total == 0:
            return -math.log(self.max_size + 1)
        # Laplace smoothing
        return float(
            np.log(self.counts[idx] + 1.0)
            - np.log(total + self.max_size + 1)
        )

    # ------------------------------------------------------------------
    # eCDF-based effective length computation
    # ------------------------------------------------------------------

    def _normalized_probs(self) -> np.ndarray:
        """Return Laplace-smoothed probability vector, shape (max_size+1,)."""
        total = self.total_weight
        if total == 0:
            n = self.max_size + 1
            return np.full(n, 1.0 / n, dtype=np.float64)
        return (self.counts + 1.0) / (total + self.max_size + 1)

    # ------------------------------------------------------------------
    # Analytical transcript effective length (salmon-style eCDF)
    # ------------------------------------------------------------------

    def _build_eff_len_cache(self) -> tuple[np.ndarray, np.ndarray]:
        """Precompute cumulative CDF and moment arrays for vectorized
        effective length computation.

        Returns (cdf, cmom) each of shape (max_size+1,) where:
            cdf[k] = sum_{l=0}^{k} P(l)
            cmom[k] = sum_{l=0}^{k} l * P(l)
        """
        probs = self._normalized_probs()  # shape (max_size+1,)
        cdf = np.cumsum(probs)
        l_vals = np.arange(len(probs), dtype=np.float64)
        cmom = np.cumsum(probs * l_vals)
        return cdf, cmom

    def compute_all_transcript_eff_lens(
        self, lengths: np.ndarray,
    ) -> np.ndarray:
        """Vectorized analytical effective length for an array of
        transcript exonic lengths.

        Uses precomputed cumulative CDF and moment arrays so
        each transcript requires only two table lookups and a
        multiply-subtract::

            eff_len = (L + 1) × (CDF[k] − P(0)) − CMOM[k]

        where ``k = min(L, max_size)`` and ``CDF / CMOM`` are
        precomputed cumulative sums of the fragment length
        distribution.

        This is the recommended entry point for computing effective
        lengths for all transcripts in a single call.

        Parameters
        ----------
        lengths : np.ndarray
            1-D array of transcript exonic (spliced) lengths.

        Returns
        -------
        np.ndarray
            float64 array of effective lengths, each ≥ 1.0.
        """
        lengths = np.asarray(lengths, dtype=np.int64)
        probs = self._normalized_probs()
        p0 = float(probs[0])
        cdf, cmom = self._build_eff_len_cache()

        # Clamp to valid range for the histogram
        k = np.clip(lengths, 0, self.max_size)

        # eff = (L + 1) * sum_{l=1}^{k} P(l) - sum_{l=1}^{k} l*P(l)
        #     = (L + 1) * (CDF[k] - P(0)) - CMOM[k]
        eff = (lengths + 1).astype(np.float64) * (cdf[k] - p0) - cmom[k]

        # Overflow bin: for transcripts longer than max_size, fragments
        # in the overflow bin can still map.  Same treatment as the
        # existing compute_exonic_eff_len method.
        overflow_mask = lengths > self.max_size
        if overflow_mask.any():
            p_overflow = float(probs[self.max_size])
            eff[overflow_mask] += (
                p_overflow * (lengths[overflow_mask] - self.max_size)
            )

        # Floor at 1.0 (same as salmon: prevents log(0) for very short
        # transcripts that are shorter than virtually all fragments)
        np.maximum(eff, 1.0, out=eff)
        return eff

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of the fragment length model.

        Includes summary statistics and a trimmed histogram (bins with
        zero counts at the tails are omitted).
        """
        # Find the range of non-zero bins for compact output
        nonzero = np.nonzero(self.counts)[0]
        if len(nonzero) > 0:
            lo, hi = int(nonzero[0]), int(nonzero[-1])
            hist_range = [lo, hi]
            hist_values = [float(v) for v in self.counts[lo:hi + 1]]
        else:
            hist_range = [0, 0]
            hist_values = []

        return {
            "summary": {
                "n_observations": int(self.n_observations),
                "total_weight": float(round(self.total_weight, 2)),
                "max_size": int(self.max_size),
                "mean": float(round(self.mean, 2)),
                "std": float(round(self.std, 2)),
                "median": float(round(self.median, 2)),
                "mode": int(self.mode),
            },
            "histogram": {
                "range": hist_range,
                "values": hist_values,
            },
        }



# ---------------------------------------------------------------------------
# FragmentLengthModels — per-category container
# ---------------------------------------------------------------------------

class FragmentLengthModels:
    """Container for per-category, global, and intergenic fragment length models.

    Wraps one ``FragmentLengthModel`` per ``SpliceType`` plus a global
    model (all categories) and an intergenic model (fragments with no
    gene overlap).

    The :meth:`observe` method routes each observation to the global
    model plus the appropriate category or intergenic model.
    """

    def __init__(self, max_size: int = DEFAULT_MAX_FRAG_SIZE):
        from .splice import SpliceType

        self.max_size = max_size
        self.global_model = FragmentLengthModel(max_size=max_size)
        self.intergenic = FragmentLengthModel(max_size=max_size)
        self.category_models: dict = {
            cat: FragmentLengthModel(max_size=max_size)
            for cat in SpliceType
        }
        # Combined gDNA fragment length model:
        # After mix_models(), holds intergenic + weighted unspliced gDNA.
        self.gdna_model = FragmentLengthModel(max_size=max_size)
        # RNA fragment length model:
        # After mix_models(), holds spliced + weighted unspliced RNA.
        self.rna_model = FragmentLengthModel(max_size=max_size)
        # Per-strand unspliced histograms (raw, pre-mixing):
        self.unspliced_same_strand = FragmentLengthModel(max_size=max_size)
        self.unspliced_opp_strand = FragmentLengthModel(max_size=max_size)

    def finalize(self) -> None:
        """Cache derived values on all sub-models for fast scoring.

        Call after all ``observe()`` calls are complete (and after
        ``mix_models()`` if applicable) and before any
        ``log_likelihood()`` calls during EM scoring.
        """
        self.global_model.finalize()
        self.intergenic.finalize()
        self.gdna_model.finalize()
        self.rna_model.finalize()
        self.unspliced_same_strand.finalize()
        self.unspliced_opp_strand.finalize()
        for m in self.category_models.values():
            m.finalize()

    def mix_models(self, s_rna: float, p_r1_sense: float) -> None:
        """Build gDNA and RNA distributions via strand-weighted mixing.

        Uses the strand-weighted histogram mixing approach: unspliced
        fragments at genic loci are split into sense/antisense pools,
        and the known strand specificity of RNA is used to estimate the
        gDNA fraction in each pool.

        Must be called AFTER ``StrandModels.finalize()`` (so strand
        specificity is known) and BEFORE ``self.finalize()`` (so the
        mixed models get their log-probability lookup tables cached).

        Parameters
        ----------
        s_rna : float
            Strand specificity (0.5 = unstranded, 1.0 = perfectly
            stranded).  From ``strand_models.strand_specificity``.
        p_r1_sense : float
            P(read 1 aligns in gene-sense direction).
            From ``strand_models.exonic_spliced.p_r1_sense``.
        """
        from .splice import SpliceType

        # 1. Map same/opp strand to sense/antisense
        if p_r1_sense >= 0.5:
            H_sense = self.unspliced_same_strand.counts.copy()
            H_anti = self.unspliced_opp_strand.counts.copy()
        else:
            H_sense = self.unspliced_opp_strand.counts.copy()
            H_anti = self.unspliced_same_strand.counts.copy()

        # 2. Total counts in each strand bin
        S = float(H_sense.sum())
        A = float(H_anti.sum())

        # 3. Estimate total RNA and gDNA in the unspliced pool
        denom = 2.0 * s_rna - 1.0
        if denom > 1e-6:
            R = max(0.0, (S - A) / denom)
            G = (S + A) - R
        else:
            # No strand information — cannot separate
            R = 0.0
            G = 0.0

        # 4. Per-strand Bayesian posterior P(gDNA | strand)
        W_sense = min(max((G * 0.5 / S) if S > 0 else 0.0, 0.0), 1.0)
        W_anti = min(max((G * 0.5 / A) if A > 0 else 0.0, 0.0), 1.0)

        # 5. Directional certainty penalty
        #
        # A pool with posterior gDNA fraction W is a mixture of gDNA and
        # RNA.  Applying a uniform scalar W to the pool's histogram
        # preserves the pool's shape, not just the gDNA component's
        # shape — a pool that is 80% RNA contributes RNA-shaped noise
        # to the gDNA model.
        #
        # The directional certainty penalty ensures each pool only
        # contributes to the model where it is the majority component:
        #   gDNA certainty: c(W) = max(2W − 1, 0)   (active when W > 0.5)
        #   RNA  certainty: c(W) = max(1 − 2W, 0)   (active when W < 0.5)
        #
        # Effective weights per pool:
        #   w_eff_gDNA = W · max(2W − 1, 0)
        #   w_eff_RNA  = (1 − W) · max(1 − 2W, 0)
        #
        # At W = 1.0: gDNA certainty = 1.0 (pure gDNA, full weight).
        # At W = 0.5: both zero (no information → excluded from both).
        # At W = 0.0: RNA certainty = 1.0 (pure RNA, full weight).
        C_gdna_sense = max(2.0 * W_sense - 1.0, 0.0)
        C_gdna_anti = max(2.0 * W_anti - 1.0, 0.0)
        C_rna_sense = max(1.0 - 2.0 * W_sense, 0.0)
        C_rna_anti = max(1.0 - 2.0 * W_anti, 0.0)

        # 6. Build mixed histograms with certainty-penalized weights
        H_intergenic = self.intergenic.counts
        H_spliced = self.category_models[SpliceType.SPLICED_ANNOT].counts

        gdna_counts = (
            H_intergenic
            + (W_sense * C_gdna_sense) * H_sense
            + (W_anti * C_gdna_anti) * H_anti
        )
        rna_counts = (
            H_spliced
            + ((1.0 - W_sense) * C_rna_sense) * H_sense
            + ((1.0 - W_anti) * C_rna_anti) * H_anti
        )

        # 7. Write into models (replace counts and total_weight)
        self.gdna_model.counts = gdna_counts
        self.gdna_model._total_weight = float(gdna_counts.sum())

        self.rna_model.counts = rna_counts
        self.rna_model._total_weight = float(rna_counts.sum())

        # 8. Log diagnostics
        logger.info(
            f"Fragment length mixing: S={S:.0f} sense, A={A:.0f} anti, "
            f"s_rna={s_rna:.4f}, W_sense={W_sense:.4f}, "
            f"W_anti={W_anti:.4f}"
        )
        logger.info(
            f"  Certainty: C_gdna=[{C_gdna_sense:.3f}, {C_gdna_anti:.3f}], "
            f"C_rna=[{C_rna_sense:.3f}, {C_rna_anti:.3f}]"
        )
        logger.info(
            f"  gDNA model: {self.gdna_model.total_weight:.0f} "
            f"weighted obs, mean={self.gdna_model.mean:.1f}"
        )
        logger.info(
            f"  RNA model: {self.rna_model.total_weight:.0f} "
            f"weighted obs, mean={self.rna_model.mean:.1f}"
        )

    @property
    def n_observations(self) -> int:
        """Total observations across all categories (delegated to global)."""
        return self.global_model.n_observations

    def observe(
        self,
        frag_length: int,
        splice_type=None,
        weight: float = 1.0,
    ) -> None:
        """Record one fragment length observation.

        Parameters
        ----------
        frag_length : int
            Fragment fragment length (must be >= 0).
        splice_type : SpliceType or None
            The fragment's count category, or ``None`` for intergenic.
        weight : float
            Observation weight (default 1.0).
        """
        self.global_model.observe(frag_length, weight)
        if splice_type is None:
            self.intergenic.observe(frag_length, weight)
            # Intergenic fragments are ~100% gDNA → train gDNA model
            self.gdna_model.observe(frag_length, weight)
        else:
            self.category_models[splice_type].observe(frag_length, weight)

    def observe_batch(
        self,
        frag_lengths: "np.ndarray",
        splice_types: "np.ndarray",
    ) -> None:
        """Record a batch of fragment length observations (vectorized).

        Parameters
        ----------
        frag_lengths : np.ndarray
            Integer array of fragment lengths.
        splice_types : np.ndarray
            Integer array of splice type values.  Dispatches to
            category sub-models based on value.
        """
        from .splice import SpliceType
        lengths = np.asarray(frag_lengths, dtype=np.intp)
        stypes = np.asarray(splice_types, dtype=np.intp)
        # Global model gets all observations
        self.global_model.observe_batch(lengths)
        # Per-category sub-models
        for cat in SpliceType:
            mask = stypes == int(cat)
            if mask.any():
                self.category_models[cat].observe_batch(lengths[mask])

    def observe_intergenic_batch(
        self,
        frag_lengths: "np.ndarray",
    ) -> None:
        """Record a batch of intergenic fragment lengths (vectorized).

        Parameters
        ----------
        frag_lengths : np.ndarray
            Integer array of intergenic fragment lengths.
        """
        lengths = np.asarray(frag_lengths, dtype=np.intp)
        self.global_model.observe_batch(lengths)
        self.intergenic.observe_batch(lengths)
        self.gdna_model.observe_batch(lengths)

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of all fragment length models."""
        from .splice import SpliceType

        d: dict = {"global": self.global_model.to_dict()}
        d["intergenic"] = self.intergenic.to_dict()
        d["gdna"] = self.gdna_model.to_dict()
        d["rna"] = self.rna_model.to_dict()
        d["unspliced_same_strand"] = self.unspliced_same_strand.to_dict()
        d["unspliced_opp_strand"] = self.unspliced_opp_strand.to_dict()
        for cat in SpliceType:
            d[cat.name.lower()] = self.category_models[cat].to_dict()
        return d

    def write_json(self, path) -> None:
        """Write all fragment length models to a single JSON file."""
        import json

        path = Path(path)
        with open(path, "w") as fh:
            json.dump(
                {"frag_length_models": self.to_dict()},
                fh,
                indent=2,
            )
        logger.info(f"Wrote fragment length models to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of all fragment length models."""
        from .splice import SpliceType

        logger.info(
            f"Fragment length models: {self.global_model.n_observations:,} "
            f"total observations"
        )
        if self.global_model.n_observations > 0:
            logger.info(
                f"  Global: mean={self.global_model.mean:.1f}, "
                f"std={self.global_model.std:.1f}, "
                f"median={self.global_model.median:.1f}, "
                f"mode={self.global_model.mode}"
            )
        for cat in SpliceType:
            m = self.category_models[cat]
            if m.n_observations > 0:
                logger.info(
                    f"  {cat.name}: {m.n_observations:,} obs, "
                    f"mean={m.mean:.1f}, mode={m.mode}"
                )
        if self.intergenic.n_observations > 0:
            logger.info(
                f"  INTERGENIC: {self.intergenic.n_observations:,} obs, "
                f"mean={self.intergenic.mean:.1f}, mode={self.intergenic.mode}"
            )
        if self.gdna_model.total_weight > 0:
            logger.info(
                f"  gDNA (mixed): {self.gdna_model.total_weight:.0f} weighted obs, "
                f"mean={self.gdna_model.mean:.1f}, mode={self.gdna_model.mode}"
            )
        if self.rna_model.total_weight > 0:
            logger.info(
                f"  RNA (mixed): {self.rna_model.total_weight:.0f} weighted obs, "
                f"mean={self.rna_model.mean:.1f}, mode={self.rna_model.mode}"
            )
        if self.unspliced_same_strand.n_observations > 0:
            logger.info(
                f"  Unspliced same-strand: "
                f"{self.unspliced_same_strand.n_observations:,} obs, "
                f"mean={self.unspliced_same_strand.mean:.1f}, "
                f"mode={self.unspliced_same_strand.mode}"
            )
        if self.unspliced_opp_strand.n_observations > 0:
            logger.info(
                f"  Unspliced opp-strand: "
                f"{self.unspliced_opp_strand.n_observations:,} obs, "
                f"mean={self.unspliced_opp_strand.mean:.1f}, "
                f"mode={self.unspliced_opp_strand.mode}"
            )
