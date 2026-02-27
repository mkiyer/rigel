"""hulkrna.bias — Positional bias model with prefix-sum acceleration.

Provides ``BiasProfile`` for per-transcript positional bias with O(1)
fragment weight queries and O(L) effective-length computation via
prefix sums.

The uniform model (current default) collapses to exactly the standard
effective-length correction::

    log P(f | t) = −log(max(t_len − l_f + 1, 1))

Future biological bias profiles (GC content, 3' degradation, etc.)
can be dropped into the same ``BiasProfile`` container and will
automatically flow through the EM via ``_apply_bias_correction``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class BiasProfile:
    """Per-transcript positional bias profile with prefix-sum acceleration.

    The prefix-sum array ``B`` has length ``L + 1`` where ``L`` is the
    transcript length.  ``B[0] = 0`` and ``B[i] = b[0] + ... + b[i−1]``
    so that ``sum(b[s..e−1]) = B[e] − B[s]`` in O(1).

    Attributes
    ----------
    prefix_sum : np.ndarray
        float64[L+1] — cumulative sum with leading zero.
    """

    prefix_sum: np.ndarray

    # ---- factories ---------------------------------------------------

    @staticmethod
    def uniform(length: int) -> BiasProfile:
        """Uniform bias: ``b[i] = 1.0`` for all positions.

        ``prefix_sum = [0, 1, 2, ..., L]``.
        """
        return BiasProfile(
            prefix_sum=np.arange(length + 1, dtype=np.float64),
        )

    # ---- queries -----------------------------------------------------

    @property
    def length(self) -> int:
        """Transcript length L (number of bases)."""
        return len(self.prefix_sum) - 1

    def fragment_weight(self, tx_start: int, tx_end: int) -> float:
        """Mean bias over the fragment footprint ``[tx_start, tx_end)``.

        Returns ``W(f|t) = (B[e] − B[s]) / (e − s)`` after clipping
        ``s = max(tx_start, 0)`` and ``e = min(tx_end, L)`` so that
        out-of-bounds fragments (e.g. near transcript ends) return
        graceful values instead of raising index errors.

        Under uniform bias the result is always 1.0 regardless of
        clipping.
        """
        L = self.length
        s = max(tx_start, 0)
        e = min(tx_end, L)
        span = e - s
        if span <= 0:
            return 0.0
        return float(
            (self.prefix_sum[e] - self.prefix_sum[s]) / span
        )

    def effective_length(self, frag_len: int) -> float:
        """Bias-aware effective length for a given fragment length.

        .. math::

            \\tilde{L}_t(l_f) = \\frac{1}{l_f}
                \\sum_{v=0}^{L - l_f}
                    \\bigl(B[v + l_f] - B[v]\\bigr)

        Under uniform bias this equals ``max(L − l_f + 1, 1)``.

        Returns at least 1.0 to prevent ``log(0)``.
        """
        L = self.length
        n_valid = L - frag_len + 1
        if n_valid <= 0:
            return 1.0
        B = self.prefix_sum
        # Vectorised difference of shifted prefix sums
        eff = float(
            (B[frag_len : frag_len + n_valid] - B[:n_valid]).sum()
            / frag_len
        )
        return max(eff, 1.0)


# ------------------------------------------------------------------
# Convenience builders
# ------------------------------------------------------------------


def build_uniform_profiles(t_lengths: np.ndarray) -> list[BiasProfile]:
    """Build a uniform (null) ``BiasProfile`` for each transcript.

    Parameters
    ----------
    t_lengths : np.ndarray
        int or float array of transcript lengths.

    Returns
    -------
    list[BiasProfile]
        One uniform profile per transcript.
    """
    return [BiasProfile.uniform(int(tl)) for tl in t_lengths]
