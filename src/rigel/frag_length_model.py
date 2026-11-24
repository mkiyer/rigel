"""rigel.frag_length_model — the fragment-length distribution.

:meth:`FragmentLengthModel.from_pmf` builds the scoring and effective-length model from a pmf that
``calibration.fl.build_fl_models`` derives from the accumulator payload; ``FragmentLengthModel(counts=...)`` wraps a raw
histogram for the QC report's statistics. The histogram is indexed by length in bases, with lengths ``>= max_size`` in
one overflow bin. The scorer penalises each base beyond ``max_size`` by a fixed log-probability (``TAIL_DECAY_LP`` in
``native/constants.h``), so a very long fragment scores below the overflow bin.
"""

from dataclasses import dataclass, field

import numpy as np

from .config import BamScanConfig


@dataclass
class FragmentLengthModel:
    """A fragment-length distribution over ``max_size + 1`` bins: ``counts[k]`` holds length *k* for
    ``k < max_size``, and ``counts[max_size]`` every length ``>= max_size``.

    :meth:`from_pmf` gives the finalized scoring model (the scorer reads ``_log_prob`` and ``_tail_base``); the
    ``counts=`` constructor gives the raw QC view, whose statistics read the histogram itself.

    Attributes
    ----------
    max_size : int
        The overflow bin's index; the scorer decays lengths above it.
    counts : np.ndarray
        Float64 histogram of shape ``(max_size + 1,)``.
    """

    max_size: int = BamScanConfig.max_frag_length
    counts: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        if self.counts is None:
            self.counts = np.zeros(self.max_size + 1, dtype=np.float64)
            self._total_weight: float = 0.0
        else:
            self._total_weight: float = float(self.counts.sum())
        self._log_prob: np.ndarray | None = None
        self._prob: np.ndarray | None = None
        # the log-probability the scorer decays a length past ``max_size`` from; 0 until finalized
        self._tail_base: float = 0.0
        # Lazy cache for _build_eff_len_cache(), populated on the first call; only a finalized model has one.
        self._cdf_cache: np.ndarray | None = None
        self._cmom_cache: np.ndarray | None = None

    @property
    def n_observations(self) -> int:
        """The total histogram weight, rounded. Rounded rather than truncated because a length pool divided by its
        own opportunity holds fractional weights, and truncating would report one fragment fewer than it holds.
        ``total_weight`` is the unrounded quantity."""
        return int(round(self._total_weight))

    # ------------------------------------------------------------------
    # Factory methods
    # ------------------------------------------------------------------

    @classmethod
    def from_pmf(
        cls,
        pmf,
        max_size: int | None = None,
    ) -> "FragmentLengthModel":
        """A finalized model from a pmf built upstream. Mass beyond ``max_size`` joins the overflow bin, and every
        bin is floored at the smallest positive float so its log is finite."""
        prob_in = np.asarray(pmf, dtype=np.float64)
        if prob_in.ndim != 1:
            raise ValueError("FragmentLengthModel.from_pmf: pmf must be one-dimensional.")
        if max_size is None:
            max_size = len(prob_in) - 1
        if max_size < 0:
            raise ValueError(
                f"FragmentLengthModel.from_pmf: max_size must be >= 0; got {max_size}."
            )
        if np.any(~np.isfinite(prob_in)) or np.any(prob_in < 0.0):
            raise ValueError("FragmentLengthModel.from_pmf: pmf must be finite and non-negative.")

        prob = np.zeros(max_size + 1, dtype=np.float64)
        n = min(len(prob_in), max_size + 1)
        prob[:n] = prob_in[:n]
        if len(prob_in) > max_size + 1:
            prob[max_size] += float(np.sum(prob_in[max_size + 1 :], dtype=np.float64))
        total = float(np.sum(prob, dtype=np.float64))
        if total <= 0.0:
            raise ValueError("FragmentLengthModel.from_pmf: pmf total mass must be > 0.")

        tiny = np.finfo(np.float64).tiny
        prob = np.maximum(prob / total, tiny)
        prob /= float(np.sum(prob, dtype=np.float64))

        model = cls(max_size=max_size)
        model.counts[:] = prob
        model._total_weight = 1.0
        model._prob = prob
        model._log_prob = np.log(prob)
        model._tail_base = float(model._log_prob[model.max_size])
        return model

    # ------------------------------------------------------------------
    # Distribution properties
    # ------------------------------------------------------------------

    @property
    def total_weight(self) -> float:
        """Sum of the histogram weights (unrounded :attr:`n_observations`)."""
        return self._total_weight

    @property
    def pmf(self) -> np.ndarray:
        """The finalized pmf of a :meth:`from_pmf` model, shape ``(max_size + 1,)``: the same vector the effective
        lengths use, so a consumer stays consistent with them. The model's own array; do not mutate."""
        return self._normalized_probs()

    @property
    def mean(self) -> float:
        """Weighted mean of the in-range fragment lengths; 0.0 when there are none.

        The overflow bin is right-censored, so counting it as a point mass at ``max_size`` would bias the mean toward
        the ceiling; it is excluded, and reported on its own (:attr:`overflow_count`).
        """
        w = self.counts[: self.max_size]
        total = float(w.sum())
        if total <= 0.0:
            return 0.0
        idx = np.arange(self.max_size, dtype=np.float64)
        return float(np.dot(idx, w) / total)

    @property
    def std(self) -> float:
        """Weighted standard deviation of the in-range fragment lengths, about :attr:`mean`; 0.0 when there are
        none."""
        w = self.counts[: self.max_size]
        total = float(w.sum())
        if total <= 0.0:
            return 0.0
        idx = np.arange(self.max_size, dtype=np.float64)
        mu = np.dot(idx, w) / total
        var = np.dot((idx - mu) ** 2, w) / total
        return float(np.sqrt(max(var, 0.0)))

    @property
    def median(self) -> float:
        """Weighted median fragment length over the full distribution; 0.0 when empty. A rank statistic survives
        right-censoring while the overflow share is below one half, so the overflow bin is kept here."""
        total = self.total_weight
        if total == 0:
            return 0.0
        cumsum = np.cumsum(self.counts)
        idx = np.searchsorted(cumsum, total / 2.0)
        return float(min(idx, self.max_size))

    @property
    def mode(self) -> int:
        """The in-range fragment length with the highest weight; 0 when there are none. The overflow bin holds the
        whole censored tail, so it is excluded: it can be the largest bin without being a peak."""
        w = self.counts[: self.max_size]
        if float(w.sum()) <= 0.0:
            return 0
        return int(np.argmax(w))

    @property
    def overflow_count(self) -> float:
        """Weight in the ``>= max_size`` overflow (right-censored) bin."""
        return float(self.counts[self.max_size])

    # ------------------------------------------------------------------
    # Effective length (salmon-style eCDF)
    # ------------------------------------------------------------------

    def _normalized_probs(self) -> np.ndarray:
        """The finalized probability vector, shape ``(max_size + 1,)``; a raw-count model refuses."""
        if self._prob is None:
            raise ValueError(
                "this FragmentLengthModel holds raw counts only; pmf() and the effective lengths need "
                "a model built by FragmentLengthModel.from_pmf"
            )
        return self._prob

    def _build_eff_len_cache(self) -> tuple[np.ndarray, np.ndarray]:
        """The cumulative pmf and first moment, memoized on a finalized model.

        Returns ``(cdf, cmom)``, each of shape ``(max_size + 1,)``::

            cdf[k]  = sum_{l=0}^{k} P(l)
            cmom[k] = sum_{l=0}^{k} l * P(l)
        """
        if self._cdf_cache is not None:
            return self._cdf_cache, self._cmom_cache
        probs = self._normalized_probs()  # shape (max_size+1,)
        cdf = np.cumsum(probs)
        l_vals = np.arange(len(probs), dtype=np.float64)
        cmom = np.cumsum(probs * l_vals)
        self._cdf_cache, self._cmom_cache = cdf, cmom
        return cdf, cmom

    def compute_all_transcript_eff_lens(
        self,
        lengths: np.ndarray,
        *,
        min_value: float = 1.0,
    ) -> np.ndarray:
        """The effective length of each transcript exonic length, from the cumulative sums of the pmf::

            eff_len = (L + 1) × (CDF[k] − P(0)) − CMOM[k],   k = min(L, max_size)

        Parameters
        ----------
        lengths : np.ndarray
            1-D array of transcript exonic (spliced) lengths.
        min_value : float, default 1.0
            Floor on the result. ``1.0`` keeps a short transcript's log finite (salmon's convention); pass ``0.0``
            where the value is a Poisson-rate denominator, so a zero-span region contributes zero.

        Returns
        -------
        np.ndarray
            float64 effective lengths, each ≥ ``min_value``.
        """
        lengths = np.asarray(lengths, dtype=np.int64)
        probs = self._normalized_probs()
        p0 = float(probs[0])
        cdf, cmom = self._build_eff_len_cache()

        k = np.clip(lengths, 0, self.max_size)

        # eff = (L + 1) * sum_{l=1}^{k} P(l) - sum_{l=1}^{k} l*P(l)
        #     = (L + 1) * (CDF[k] - P(0)) - CMOM[k]
        eff = (lengths + 1).astype(np.float64) * (cdf[k] - p0) - cmom[k]

        np.maximum(eff, min_value, out=eff)
        return eff
