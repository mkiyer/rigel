"""rigel.calibration._diagnostics \u2014 fractional payload diagnostics.

Replaces the legacy 8-mask integer breakdown. Carries:

* Fragment-level integer counters (sum to ``n_total`` when supplied).
* Fractional mass summaries from ``region_counts``, ``channel_mass``,
  ``signature_mass``.
* FL pool totals keyed by name.

``to_summary_dict()`` returns a JSON-safe nested dict suitable for
``summary.json``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .fractional_evidence import is_exon_any, is_intergenic, is_intron_only
from .scan_payload import CalibrationScanPayload
from .signature import FL_POOL_NAMES, N_CHANNELS, channel_tuple


__all__ = ["Diagnostics"]


_COMPARTMENT_NAMES = ("CONTAINED", "BOUNDARY_LEFT", "BOUNDARY_RIGHT")
_SPLICE_NAMES = ("UNSPLICED", "SPLICED")
_STRAND_NAMES = ("POS", "NEG")
_COARSE_NAMES = ("INTERGENIC", "INTRON", "EXON")


@dataclass(frozen=True, slots=True)
class Diagnostics:
    """Fractional payload summary."""

    # ---- exact integer counters --------------------------------------------
    n_observed: int
    n_excluded_multimap: int
    n_excluded_chimera: int
    n_excluded_artifact: int
    n_excluded_strand_ambig: int
    n_unobserved: int
    n_unannotated_ref: int  # subset of n_observed
    n_fl_unavailable: int  # subset of n_observed

    # ---- fractional mass summaries -----------------------------------------
    total_region_mass: float
    mass_by_coarse_class: Mapping[str, float]  # INTERGENIC, INTRON, EXON
    mass_by_compartment: Mapping[str, float]  # CONTAINED, BOUNDARY_LEFT, BOUNDARY_RIGHT
    mass_by_splice: Mapping[str, float]  # UNSPLICED, SPLICED
    mass_by_strand: Mapping[str, float]  # POS, NEG
    mass_by_signature: np.ndarray  # float64[16]

    # ---- FL pool totals ----------------------------------------------------
    fl_pool_total: Mapping[str, float]  # 6 named pools

    def total(self) -> int:
        """Sum of all fragment-level exclusion + observation counters."""
        return (
            self.n_observed
            + self.n_excluded_multimap
            + self.n_excluded_chimera
            + self.n_excluded_artifact
            + self.n_excluded_strand_ambig
            + self.n_unobserved
        )

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "n_observed": int(self.n_observed),
            "n_excluded_multimap": int(self.n_excluded_multimap),
            "n_excluded_chimera": int(self.n_excluded_chimera),
            "n_excluded_artifact": int(self.n_excluded_artifact),
            "n_excluded_strand_ambig": int(self.n_excluded_strand_ambig),
            "n_unobserved": int(self.n_unobserved),
            "n_unannotated_ref": int(self.n_unannotated_ref),
            "n_fl_unavailable": int(self.n_fl_unavailable),
            "total_region_mass": float(self.total_region_mass),
            "mass_by_coarse_class": {k: float(v) for k, v in self.mass_by_coarse_class.items()},
            "mass_by_compartment": {k: float(v) for k, v in self.mass_by_compartment.items()},
            "mass_by_splice": {k: float(v) for k, v in self.mass_by_splice.items()},
            "mass_by_strand": {k: float(v) for k, v in self.mass_by_strand.items()},
            "mass_by_signature": [float(v) for v in self.mass_by_signature],
            "fl_pool_total": {k: float(v) for k, v in self.fl_pool_total.items()},
        }

    @classmethod
    def from_payload(
        cls,
        payload: CalibrationScanPayload,
        *,
        signature: np.ndarray | None = None,
    ) -> "Diagnostics":
        """Build a Diagnostics from a payload (optionally with region signatures).

        ``signature`` is the per-region uint8 signature array. When omitted,
        ``mass_by_coarse_class`` is populated from ``signature_mass``
        (global signature marginals, equivalent to the per-region split).
        """
        ch = payload.channel_mass  # float64[12]

        # Marginals over the 12-channel index by unpacking (comp, splice, strand).
        comp_totals = {name: 0.0 for name in _COMPARTMENT_NAMES}
        splice_totals = {name: 0.0 for name in _SPLICE_NAMES}
        strand_totals = {name: 0.0 for name in _STRAND_NAMES}
        for c in range(N_CHANNELS):
            comp_idx, splice_idx, strand_idx = channel_tuple(c)
            mass = float(ch[c])
            comp_totals[_COMPARTMENT_NAMES[comp_idx]] += mass
            splice_totals[_SPLICE_NAMES[splice_idx]] += mass
            strand_totals[_STRAND_NAMES[strand_idx]] += mass

        sig_mass = np.asarray(payload.signature_mass, dtype=np.float64)
        if signature is not None:
            sig = np.asarray(signature, dtype=np.uint8)
            if sig.shape != (payload.n_regions,):
                raise ValueError(
                    f"Diagnostics.from_payload: signature shape {sig.shape} does not match "
                    f"payload n_regions ({payload.n_regions})."
                )
            if np.any(sig >= sig_mass.size):
                raise ValueError(
                    "Diagnostics.from_payload: signature contains values outside [0, 15]."
                )
            row_mass = payload.region_counts.sum(axis=1, dtype=np.float64)
            observed_sig_mass = np.bincount(sig, weights=row_mass, minlength=sig_mass.size)[
                : sig_mass.size
            ]
            tol = max(1e-3, 1e-5 * max(float(row_mass.sum()), 1.0))
            if not np.allclose(observed_sig_mass, sig_mass, rtol=0.0, atol=tol):
                max_diff = float(np.max(np.abs(observed_sig_mass - sig_mass)))
                raise ValueError(
                    "Diagnostics.from_payload: payload signature_mass does not match "
                    f"region_counts grouped by region signature (max_diff={max_diff}, tol={tol})."
                )
        # mass_by_coarse_class from signature_mass marginals.
        sig_idx = np.arange(sig_mass.size, dtype=np.uint8)
        coarse = {
            "INTERGENIC": float(sig_mass[is_intergenic(sig_idx)].sum()),
            "INTRON": float(sig_mass[is_intron_only(sig_idx)].sum()),
            "EXON": float(sig_mass[is_exon_any(sig_idx)].sum()),
        }
        fl_pool_total = {
            FL_POOL_NAMES[i]: float(payload.fl_pool_total[i]) for i in range(len(FL_POOL_NAMES))
        }

        return cls(
            n_observed=int(payload.n_observed),
            n_excluded_multimap=int(payload.n_excluded_multimap),
            n_excluded_chimera=int(payload.n_excluded_chimera),
            n_excluded_artifact=int(payload.n_excluded_artifact),
            n_excluded_strand_ambig=int(payload.n_excluded_strand_ambig),
            n_unobserved=int(payload.n_unobserved),
            n_unannotated_ref=int(payload.n_unannotated_ref),
            n_fl_unavailable=int(payload.n_fl_unavailable),
            total_region_mass=float(payload.region_counts.sum(dtype=np.float64)),
            mass_by_coarse_class=coarse,
            mass_by_compartment=comp_totals,
            mass_by_splice=splice_totals,
            mass_by_strand=strand_totals,
            mass_by_signature=sig_mass,
            fl_pool_total=fl_pool_total,
        )
