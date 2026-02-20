"""
hulkrna.insert_model — Insert size distribution model.

Learns the insert size distribution from fragment-to-transcript
mappings produced by Pass 1.  Fragments where all candidate
transcripts yield the *same* insert size contribute to training;
ambiguous fragments are deferred to Bayesian counting where the
learned distribution provides a likelihood term.

The distribution is stored as a histogram (float vector) indexed
by insert size in bases.  Sizes >= ``max_size`` are clamped into
a single overflow bin.

Typical RNA-seq insert sizes range from 100–500 bp.  The default
``max_size=1000`` accommodates the vast majority of libraries.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class InsertSizeModel:
    """Insert size distribution model.

    The histogram has ``max_size + 1`` bins:
      - ``counts[k]`` for ``0 <= k < max_size``: fragments with
        insert size exactly *k*.
      - ``counts[max_size]``: overflow bin for insert sizes >= max_size.

    Training
    --------
    Call :meth:`observe` with each unambiguous insert size.  After
    training, :meth:`log_likelihood` returns the log-probability
    of a given insert size under the learned distribution.

    Attributes
    ----------
    max_size : int
        Maximum insert size tracked individually (sizes >= this
        go into the overflow bin).
    counts : np.ndarray
        Float64 histogram of shape ``(max_size + 1,)``.
    n_observations : int
        Total number of observations added.
    """

    max_size: int = 1000
    counts: np.ndarray = field(default=None, repr=False)
    n_observations: int = 0

    def __post_init__(self):
        if self.counts is None:
            self.counts = np.zeros(self.max_size + 1, dtype=np.float64)
            self._total_weight: float = 0.0
        else:
            self._total_weight: float = float(self.counts.sum())
        self._log_prob: np.ndarray | None = None
        self._finalized: bool = False

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def observe(self, insert_size: int, weight: float = 1.0) -> None:
        """Record one insert size observation.

        Parameters
        ----------
        insert_size : int
            Fragment insert size (must be >= 0).
        weight : float
            Observation weight (default 1.0).
        """
        idx = min(max(insert_size, 0), self.max_size)
        self.counts[idx] += weight
        self.n_observations += 1
        self._total_weight += weight

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
        """Weighted mean insert size.

        The overflow bin is treated as having insert size = max_size.
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
        """Weighted median insert size.

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
        """Insert size with the highest count.

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
        self._finalized = True

    # ------------------------------------------------------------------
    # Likelihood for Bayesian counting
    # ------------------------------------------------------------------

    def log_likelihood(self, insert_size: int) -> float:
        """Log-probability of an insert size under the learned distribution.

        Uses add-one (Laplace) smoothing to avoid -inf for unseen sizes.

        Parameters
        ----------
        insert_size : int
            Query insert size.

        Returns
        -------
        float
            Log-probability (natural log).
        """
        total = self._total_weight
        idx = min(max(insert_size, 0), self.max_size)
        if self._finalized:
            return float(self._log_prob[idx])
        if total == 0:
            # Uniform prior when no data
            return -np.log(self.max_size + 1)
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

    def compute_exonic_eff_len(self, exon_lengths: list[int] | np.ndarray) -> float:
        """eCDF-weighted effective length for gDNA fragments on exons.

        For each exon of length *E*, the expected number of gDNA
        fragment start positions that produce a fragment fully
        contained within the exon is::

            Σ_{l=1}^{E}  P(l) × (E − l + 1)

        where *P(l)* is the Laplace-smoothed gDNA insert-size
        probability.  Summing across all exons gives the total
        exonic effective length for the gene.

        This is exact under the assumption that gDNA fragments
        are uniformly distributed across the genome.

        Parameters
        ----------
        exon_lengths : list[int] or np.ndarray
            Lengths of the (merged) exons for one gene.

        Returns
        -------
        float
            Total exonic effective length (≥ 1.0).
        """
        probs = self._normalized_probs()  # P(l) for l=0..max_size
        total = 0.0
        for E in exon_lengths:
            E = int(E)
            if E <= 0:
                continue
            # l ranges from 1..min(E, max_size)
            L_max = min(E, self.max_size)
            # positions[l] = E - l + 1  for l = 1..L_max
            l_vals = np.arange(1, L_max + 1)
            positions = E - l_vals + 1  # (E, E-1, ..., E-L_max+1)
            total += float(np.dot(probs[l_vals], positions))
            # Overflow bin: fragments of length > max_size can only
            # fit if E > max_size.  They get 1 position each
            # (the exon must fully contain them).
            if E > self.max_size:
                total += float(probs[self.max_size]) * (E - self.max_size)
        return max(total, 1.0)

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of the insert size model.

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

    def write_json(self, path: Path | str) -> None:
        """Write the insert size model to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output JSON file path.
        """
        import json

        path = Path(path)
        d = self.to_dict()

        with open(path, "w") as fh:
            json.dump(
                {"insert_size_model": d},
                fh,
                indent=2,
            )

        logger.info(f"Wrote insert size model to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of the learned insert size model."""
        logger.info(
            f"Insert size model: {self.n_observations:,} observations, "
            f"mean={self.mean:.1f}, std={self.std:.1f}, "
            f"median={self.median:.1f}, mode={self.mode}"
        )


# ---------------------------------------------------------------------------
# InsertSizeModels — per-category container
# ---------------------------------------------------------------------------

class InsertSizeModels:
    """Container for per-category, global, and intergenic insert size models.

    Wraps one ``InsertSizeModel`` per ``SpliceType`` plus a global
    model (all categories) and an intergenic model (fragments with no
    gene overlap).

    The :meth:`observe` method routes each observation to the global
    model plus the appropriate category or intergenic model.
    """

    def __init__(self, max_size: int = 1000):
        from .categories import SpliceType

        self.max_size = max_size
        self.global_model = InsertSizeModel(max_size=max_size)
        self.intergenic = InsertSizeModel(max_size=max_size)
        self.category_models: dict = {
            cat: InsertSizeModel(max_size=max_size)
            for cat in SpliceType
        }
        # Combined gDNA insert size model:
        # Trained from intergenic + antisense-genic fragments.
        # Provides a richer gDNA insert profile than intergenic alone
        # (which can be sparse for small genomes / targeted panels).
        self.gdna_model = InsertSizeModel(max_size=max_size)

    def finalize(self) -> None:
        """Cache derived values on all sub-models for fast scoring.

        Call after all ``observe()`` calls are complete and before
        any ``log_likelihood()`` calls during EM scoring.
        """
        self.global_model.finalize()
        self.intergenic.finalize()
        self.gdna_model.finalize()
        for m in self.category_models.values():
            m.finalize()

    @property
    def n_observations(self) -> int:
        """Total observations across all categories (delegated to global)."""
        return self.global_model.n_observations

    def observe(
        self,
        insert_size: int,
        splice_type=None,
        weight: float = 1.0,
    ) -> None:
        """Record one insert size observation.

        Parameters
        ----------
        insert_size : int
            Fragment insert size (must be >= 0).
        splice_type : SpliceType or None
            The fragment's count category, or ``None`` for intergenic.
        weight : float
            Observation weight (default 1.0).
        """
        self.global_model.observe(insert_size, weight)
        if splice_type is None:
            self.intergenic.observe(insert_size, weight)
            # Intergenic fragments are ~100% gDNA → train gDNA model
            self.gdna_model.observe(insert_size, weight)
        else:
            self.category_models[splice_type].observe(insert_size, weight)

    def observe_gdna(self, insert_size: int, weight: float = 1.0) -> None:
        """Record an insert size observation for the gDNA model only.

        Used for antisense genic fragments (exonic + intronic) which
        are enriched for gDNA but already routed to their category
        model via :meth:`observe`.
        """
        self.gdna_model.observe(insert_size, weight)

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of all insert size models."""
        from .categories import SpliceType

        d: dict = {"global": self.global_model.to_dict()}
        d["intergenic"] = self.intergenic.to_dict()
        d["gdna"] = self.gdna_model.to_dict()
        for cat in SpliceType:
            d[cat.name.lower()] = self.category_models[cat].to_dict()
        return d

    def write_json(self, path) -> None:
        """Write all insert size models to a single JSON file."""
        import json
        from pathlib import Path

        path = Path(path)
        with open(path, "w") as fh:
            json.dump(
                {"insert_size_models": self.to_dict()},
                fh,
                indent=2,
            )
        logger.info(f"Wrote insert size models to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of all insert size models."""
        from .categories import SpliceType

        logger.info(
            f"Insert size models: {self.global_model.n_observations:,} "
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
        if self.gdna_model.n_observations > 0:
            logger.info(
                f"  gDNA (combined): {self.gdna_model.n_observations:,} obs, "
                f"mean={self.gdna_model.mean:.1f}, mode={self.gdna_model.mode}"
            )
