"""Public result container produced by :func:`rigel.calibration.calibrate_gdna`."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..frag_length_model import FragmentLengthModel


@dataclass(frozen=True)
class CalibrationResult:
    """Continuous gDNA deconvolution results."""

    region_e_gdna: np.ndarray
    region_n_total: np.ndarray
    gdna_fl_model: FragmentLengthModel | None
    lambda_gdna: float
    strand_specificity: float

    region_gamma: np.ndarray | None = None
    region_gamma_strand: np.ndarray | None = None
    mu_R: float | None = None
    sigma_R: float | None = None
    mixing_pi: float | None = None
    mixing_pi_soft: float | None = None
    kappa_G: float = 0.0
    kappa_R: float = 0.0
    em_n_iter: int = 0
    em_converged: bool = False
    n_eligible: int = 0
    n_soft: int = 0
    n_spliced_hard: int = 0

    def to_summary_dict(self) -> dict:
        total_e_gdna = float(self.region_e_gdna.sum())
        total_n = float(self.region_n_total.sum())
        d: dict = {
            "lambda_gdna": _round_or_none(self.lambda_gdna, 8),
            "mu_R": _round_or_none(self.mu_R, 4),
            "sigma_R": _round_or_none(self.sigma_R, 4),
            "mixing_pi": _round_or_none(self.mixing_pi, 4),
            "mixing_pi_soft": _round_or_none(self.mixing_pi_soft, 4),
            "kappa_G": round(float(self.kappa_G), 4),
            "kappa_R": round(float(self.kappa_R), 4),
            "em_n_iter": int(self.em_n_iter),
            "em_converged": bool(self.em_converged),
            "n_eligible": int(self.n_eligible),
            "n_soft": int(self.n_soft),
            "n_spliced_hard": int(self.n_spliced_hard),
            "total_expected_gdna": round(total_e_gdna, 1),
            "gdna_fraction": round(total_e_gdna / max(total_n, 1.0), 4),
            "strand_specificity": round(self.strand_specificity, 4),
            "gdna_fl_mean": None,
            "gdna_fl_observations": 0,
        }
        if self.gdna_fl_model is not None:
            d["gdna_fl_mean"] = round(self.gdna_fl_model.mean, 2)
            d["gdna_fl_observations"] = self.gdna_fl_model.n_observations
        return d


def _round_or_none(x, digits):
    if x is None:
        return None
    return round(float(x), digits)
