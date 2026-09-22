"""
rigel.frag_length_model — Fragment length distribution model.

Two shapes, and nothing trains it. :meth:`FragmentLengthModel.from_pmf` builds the scoring and
effective-length model from a pmf that ``calibration.fl.build_fl_models`` derived from the
accumulator payload; ``FragmentLengthModel(counts=...)`` wraps a raw histogram for the QC report's
summary statistics.

The distribution is stored as a histogram (float vector) indexed by fragment length in bases. Sizes
>= ``max_size`` are clamped into a single overflow bin.

Fragments whose length exceeds ``max_size`` receive exponential tail decay in the scorer: each
additional base-pair beyond the boundary incurs a fixed log-probability penalty (the scorer's
``TAIL_DECAY_LP = log(0.99) ≈ −0.01`` per bp, ``native/constants.h``), so very long fragments are
penalised heavily rather than assigned the same probability as the overflow bin.
"""

from dataclasses import dataclass, field

import numpy as np

#: Default maximum fragment length tracked individually.
#: Sizes >= this go into a single overflow bin.  Must match
#: ``BamScanConfig.max_frag_length`` (config.py) for consistency
#: between standalone and production usage.
DEFAULT_MAX_FRAG_SIZE: int = 1000

# One pseudo-observation spread across the full FL support is enough to keep
# unseen lengths finite without pulling short libraries toward the midpoint of
# the histogram.
_UNSEEN_FL_SMOOTHING_ESS: float = 1.0


@dataclass
class FragmentLengthModel:
    """Fragment length distribution model.

    The histogram has ``max_size + 1`` bins:
      - ``counts[k]`` for ``0 <= k < max_size``: fragments with
        fragment length exactly *k*.
      - ``counts[max_size]``: overflow bin for fragment lengths >= max_size.

    :meth:`from_pmf` gives the finalized scoring model (``_log_prob`` and ``_tail_base`` are what
    the scorer reads); the ``counts=`` constructor gives the raw QC view, whose statistics read the
    histogram itself.

    Attributes
    ----------
    max_size : int
        Maximum fragment length tracked individually (sizes >= this
        go into the overflow bin).  Queries above this receive
        exponential tail decay in the scorer.
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
        self._prob: np.ndarray | None = None
        self._stats_use_prob: bool = False
        self._finalized: bool = False
        # Lazy cache for _build_eff_len_cache(), populated on the first call on a finalized model.
        self._cdf_cache: np.ndarray | None = None
        self._cmom_cache: np.ndarray | None = None

    @property
    def n_observations(self) -> int:
        """Number of observations, derived from the total histogram weight.

        ROUNDED, not truncated, and that matters because some histograms are fractional. A length
        pool divided by its own opportunity holds non-integer weights by construction (the de-tilt
        preserves the pool's total but redistributes it), so truncating would report one fewer
        fragment than the histogram beside it holds, and a reader has no way to tell a rounding
        artefact from a lost fragment. ``total_weight`` is the unrounded quantity.
        """
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
        """Create a finalized model from an already-normalized PMF.

        This is intended for calibrated scoring surfaces whose probabilities
        have already been regularized and normalized upstream, so it adds no
        unseen-bin reserve.
        """
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
        model._stats_use_prob = True
        model._tail_base = float(model._log_prob[model.max_size])
        model._finalized = True
        model._cdf_cache = None
        model._cmom_cache = None
        return model

    # ------------------------------------------------------------------
    # Distribution properties
    # ------------------------------------------------------------------

    @property
    def total_weight(self) -> float:
        """Sum of all histogram weights (may differ from n_observations
        if non-unit weights are used)."""
        return self._total_weight

    @property
    def pmf(self) -> np.ndarray:
        """Public probability vector, shape ``(max_size + 1,)``.

        Returns the finalized pmf of a :meth:`from_pmf` model, or a smoothed
        estimate from the raw counts of a QC view. This is the
        *same* vector used internally by
        :meth:`compute_all_transcript_eff_lens` and the analytical
        moment helpers, so consumers that need the underlying PMF (e.g.
        the calibration boundary-crossing exposure) stay consistent
        with effective-length math without reaching into private state.

        Always returns a freshly normalized ``float64`` array; do not
        mutate.
        """
        return self._normalized_probs()

    @property
    def mean(self) -> float:
        """Weighted mean of the *in-range* fragment lengths.

        The ``>= max_size`` overflow bin is right-censored — its members' exact
        lengths are unknown, only that they exceed ``max_size`` — so treating it
        as a point mass at ``max_size`` would bias the mean toward the ceiling.
        It is therefore EXCLUDED; this reports the mean of the resolved lengths.
        The censored count is surfaced separately (``to_dict()`` → ``overflow``).
        Returns 0.0 if there are no in-range observations.
        """
        w = self._stats_weights()[: self.max_size]
        total = float(w.sum())
        if total <= 0.0:
            return 0.0
        idx = np.arange(self.max_size, dtype=np.float64)
        return float(np.dot(idx, w) / total)

    @property
    def std(self) -> float:
        """Weighted standard deviation of the in-range fragment lengths.

        Excludes the ``>= max_size`` overflow bin, consistently with :meth:`mean`
        (the spread is reported about the in-range mean). Returns 0.0 if there
        are no in-range observations.
        """
        w = self._stats_weights()[: self.max_size]
        total = float(w.sum())
        if total <= 0.0:
            return 0.0
        idx = np.arange(self.max_size, dtype=np.float64)
        mu = np.dot(idx, w) / total
        var = np.dot((idx - mu) ** 2, w) / total
        return float(np.sqrt(max(var, 0.0)))

    @property
    def median(self) -> float:
        """Weighted median fragment length (over the full distribution).

        Unlike :meth:`mean` / :meth:`mode`, the median is a rank statistic and
        stays identifiable under right-censoring as long as the overflow fraction
        is below 50%, so the ``>= max_size`` bin is retained here (it correctly
        places the censored tail on the high side). Returns 0.0 if no observations.
        """
        probs = self._stats_probs()
        if probs is not None:
            cumsum = np.cumsum(probs)
            idx = np.searchsorted(cumsum, 0.5)
            return float(min(idx, self.max_size))
        total = self.total_weight
        if total == 0:
            return 0.0
        cumsum = np.cumsum(self.counts)
        idx = np.searchsorted(cumsum, total / 2.0)
        return float(min(idx, self.max_size))

    @property
    def mode(self) -> int:
        """In-range fragment length with the highest weight.

        Excludes the ``>= max_size`` overflow bin: that bin aggregates the whole
        censored tail, so it can be the single largest bin without representing a
        genuine peak (see :meth:`mean`). Returns 0 if no in-range observations.
        """
        w = self._stats_weights()[: self.max_size]
        if float(w.sum()) <= 0.0:
            return 0
        return int(np.argmax(w))

    def _stats_probs(self) -> np.ndarray | None:
        """Finalized probability vector for summary statistics, when informative."""
        if self._finalized and self._prob is not None and self._stats_use_prob:
            return self._prob
        return None

    def _stats_weights(self) -> np.ndarray:
        """Weight vector for summary statistics: the finalized posterior-predictive
        PMF when informative, else the raw counts."""
        probs = self._stats_probs()
        return probs if probs is not None else self.counts

    @property
    def overflow_count(self) -> float:
        """Weight in the ``>= max_size`` overflow (right-censored) bin."""
        return float(self.counts[self.max_size])

    # ------------------------------------------------------------------
    # eCDF-based effective length computation
    # ------------------------------------------------------------------

    def _normalized_probs(self) -> np.ndarray:
        """Return the smoothed probability vector, shape (max_size+1,)."""
        if self._finalized and self._prob is not None:
            return self._prob
        total = self.total_weight
        if total == 0:
            n = self.max_size + 1
            return np.full(n, 1.0 / n, dtype=np.float64)
        n = self.max_size + 1
        smoothing_per_bin = _UNSEEN_FL_SMOOTHING_ESS / n
        return (self.counts + smoothing_per_bin) / (total + _UNSEEN_FL_SMOOTHING_ESS)

    # ------------------------------------------------------------------
    # Analytical transcript effective length (salmon-style eCDF)
    # ------------------------------------------------------------------

    def _build_eff_len_cache(self) -> tuple[np.ndarray, np.ndarray]:
        """Precompute cumulative CDF and moment arrays for vectorized
        effective length computation.

        Returns (cdf, cmom) each of shape (max_size+1,) where:
            cdf[k] = sum_{l=0}^{k} P(l)
            cmom[k] = sum_{l=0}^{k} l * P(l)

        A finalized model is immutable, so its cumulative arrays are
        computed once and memoized; a raw QC view recomputes on every call.
        """
        if self._finalized and self._cdf_cache is not None:
            return self._cdf_cache, self._cmom_cache
        probs = self._normalized_probs()  # shape (max_size+1,)
        cdf = np.cumsum(probs)
        l_vals = np.arange(len(probs), dtype=np.float64)
        cmom = np.cumsum(probs * l_vals)
        if self._finalized:
            self._cdf_cache = cdf
            self._cmom_cache = cmom
        return cdf, cmom

    def compute_all_transcript_eff_lens(
        self,
        lengths: np.ndarray,
        *,
        min_value: float = 1.0,
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
        min_value : float, default 1.0
            Lower floor applied to the result. Default ``1.0`` matches
            salmon's transcript-scoring convention (prevents log(0)
            for very short transcripts). Pass ``0.0`` when the value
            is used as a Poisson-rate denominator (gDNA density
            estimation), where a zero-span region must contribute
            zero effective length.

        Returns
        -------
        np.ndarray
            float64 array of effective lengths, each ≥ ``min_value``.
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

        np.maximum(eff, min_value, out=eff)
        return eff

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of the fragment length model.

        Includes summary statistics and a trimmed histogram (bins with
        zero counts at the tails are omitted). ``mean`` / ``std`` / ``mode``
        describe the in-range lengths; the right-censored ``>= max_size`` tail is
        reported separately under ``summary.overflow`` (count + fraction).
        """
        # Find the range of non-zero bins for compact output
        nonzero = np.nonzero(self.counts)[0]
        if len(nonzero) > 0:
            lo, hi = int(nonzero[0]), int(nonzero[-1])
            hist_range = [lo, hi]
            hist_values = [float(v) for v in self.counts[lo : hi + 1]]
        else:
            hist_range = [0, 0]
            hist_values = []

        total = self.total_weight
        overflow = self.overflow_count

        return {
            "summary": {
                "n_observations": int(self.n_observations),
                "total_weight": float(round(total, 2)),
                "max_size": int(self.max_size),
                "mean": float(round(self.mean, 2)),
                "std": float(round(self.std, 2)),
                "median": float(round(self.median, 2)),
                "mode": int(self.mode),
                # Right-censored fragments (length >= max_size), excluded from
                # mean/std/mode above; reported here for QC.
                "overflow": {
                    "count": float(round(overflow, 2)),
                    "fraction": float(round(overflow / total, 6)) if total > 0 else 0.0,
                },
            },
            "histogram": {
                "range": hist_range,
                "values": hist_values,
            },
        }


# There is deliberately no plural container of raw per-SpliceType histograms here: mixing a genomic
# footprint for one subset of fragments with a transcript-space length for a disjoint one gives an
# anchor measured differently from the pools it anchors.
#
# `FragmentLengthModel` above is the scoring and effective-length model, built by `from_pmf` from a
# pmf that `calibration.fl.build_fl_models` derives from the accumulator payload. The per-splice-type
# QC counts come from the scanner's own census (`rigel.splice.census_field`).
