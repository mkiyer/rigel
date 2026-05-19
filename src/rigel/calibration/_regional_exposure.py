"""rigel.calibration._regional_exposure — per-region gDNA exposure weights.

Hybrid-capture gDNA is not uniformly exposed across the genome. This
module learns a per-region relative gDNA exposure weight ``A_r in (0, 1]``
from the conservative calibration evidence already collected by v6 and
provides vectorized lookups for downstream scorers.

See ``docs/calibration/gdna_exposure_plan_v3.md`` for the full design.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from ._arrays import PayloadArrays, RegionArrays
from ._exposure import boundary_crossing_exposure
from ._orient import ORIENT_OPP, ORIENT_SAME, StrandSummary
from .density_global import (
    GlobalDensityTable,
    GlobalGdnaDensity,
    _gdna_count_moment,
    _strand_identifiable_rows,
    l_eff_contained,
    strand_correction_usable,
)
from ..frag_length_model import FragmentLengthModel
from .regions import RegionType

__all__ = [
    "ExposureMode",
    "RegionalGdnaExposure",
    "RegionalWeightApplicationStats",
    "LOG_A_FLOOR",
    "LOG_RHO_CLIP_NATS",
    "REFERENCE_QUANTILE",
]

ExposureMode = Literal["uniform", "regional"]

REFERENCE_QUANTILE = 0.95
LOG_A_FLOOR = float(np.log(1.0e-4))
LOG_RHO_FLOOR = float(np.log(np.finfo(np.float64).tiny))
LOG_RHO_CLIP_NATS = 8.0
SPREAD_EPS = 1.0e-12
Z_Q95 = 1.6448536269514722

_CLASS_INTERGENIC = "INTERGENIC"
_CLASS_INTRON = "INTRON"
_CLASS_EXON = "EXON-INTRON"
_CLASS_ORDER = (_CLASS_INTERGENIC, _CLASS_INTRON, _CLASS_EXON)
_INT64_MIN = np.iinfo(np.int64).min


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class RegionalWeightApplicationStats:
    """Per-call counters for `_apply_unit_gdna_weights`.

    Emitted in ``summary.json["regional_exposure"]["application"]`` so the
    user can audit how many units actually received a weight versus how
    many were skipped for each documented reason.
    """

    n_units_seen: int = 0
    n_units_weighted: int = 0
    n_units_no_gdna: int = 0
    n_units_missing_midpoint: int = 0
    n_units_cross_ref_skipped: int = 0

    def to_dict(self) -> dict[str, int]:
        return {
            "n_units_seen": int(self.n_units_seen),
            "n_units_weighted": int(self.n_units_weighted),
            "n_units_no_gdna": int(self.n_units_no_gdna),
            "n_units_missing_midpoint": int(self.n_units_missing_midpoint),
            "n_units_cross_ref_skipped": int(self.n_units_cross_ref_skipped),
        }


@dataclass(frozen=True, slots=True)
class RegionalGdnaExposure:
    """Per-region gDNA exposure weight ``A_r`` and lookup table.

    Arrays are aligned to ``RegionArrays`` sorted (ref_id, start) order.
    The region geometry is copied here so the lookup does not need to
    materialize a separate ``RegionIndexPy``.
    """

    rho_hat: np.ndarray
    log_weight: np.ndarray
    weight: np.ndarray
    mode: ExposureMode
    rho_ref: float
    n_at_floor: int
    per_class: dict[str, dict[str, float]] = field(default_factory=dict)

    ref_offsets: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    ref_id: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int32))
    start: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    end: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))

    # -- constructors ------------------------------------------------------

    @classmethod
    def uniform(cls, region_arrays: RegionArrays) -> "RegionalGdnaExposure":
        """Identity exposure: ``A_r == 1`` and ``log A_r == 0`` for all rows."""
        R = int(region_arrays.start.size)
        zeros = np.zeros(R, dtype=np.float64)
        ones = np.ones(R, dtype=np.float64)
        return cls(
            rho_hat=zeros,
            log_weight=zeros.copy(),
            weight=ones,
            mode="uniform",
            rho_ref=0.0,
            n_at_floor=0,
            per_class={},
            ref_offsets=np.asarray(region_arrays.ref_offsets, dtype=np.int32).copy(),
            ref_id=np.asarray(region_arrays.ref_id, dtype=np.int32).copy(),
            start=np.asarray(region_arrays.start, dtype=np.int64).copy(),
            end=np.asarray(region_arrays.end, dtype=np.int64).copy(),
        )

    @classmethod
    def build(
        cls,
        region_arrays: RegionArrays,
        payload_arrays: PayloadArrays,
        global_densities: GlobalDensityTable,
        gdna_fl: FragmentLengthModel,
        *,
        strand_summary: StrandSummary | None = None,
        splicing_anchor_tolerance: int = 0,
        enabled: bool = True,
    ) -> "RegionalGdnaExposure":
        """Build regional exposure or return uniform when disabled."""
        if not enabled:
            return cls.uniform(region_arrays)

        R = int(region_arrays.start.size)
        if R == 0:
            return cls.uniform(region_arrays)

        # Per-region exposure E_r and conservative count Y_r per class.
        spans = (region_arrays.end - region_arrays.start).astype(np.int64, copy=False)
        b_cross = boundary_crossing_exposure(
            gdna_fl, splicing_anchor_tolerance=int(splicing_anchor_tolerance)
        )

        type_arr = region_arrays.type
        is_intergenic = type_arr == int(RegionType.INTERGENIC)
        is_intron = type_arr == int(RegionType.INTRON)
        is_exon = type_arr == int(RegionType.EXON)

        E = np.zeros(R, dtype=np.float64)
        # Intergenic + intron containment.
        ig_in_mask = is_intergenic | is_intron
        if ig_in_mask.any():
            E[ig_in_mask] = l_eff_contained(spans[ig_in_mask], gdna_fl)
        # Exon boundary flux (eligibility filter: any side that is allowed
        # to contribute to boundary crossing).
        bf_left = region_arrays.bf_left.astype(bool, copy=False)
        bf_right = region_arrays.bf_right.astype(bool, copy=False)
        eligible_left = (is_exon & bf_left).astype(np.int64)
        eligible_right = (is_exon & bf_right).astype(np.int64)
        sides = eligible_left + eligible_right
        if b_cross > 0.0 and is_exon.any():
            E[is_exon] = sides[is_exon].astype(np.float64) * b_cross

        # Strand-corrected count moments per channel.
        Y = np.zeros(R, dtype=np.float64)
        Y[is_intergenic] = payload_arrays.intergenic_per_region[is_intergenic].astype(
            np.float64, copy=False
        )

        # Set up strand-correction once.
        strand_active = (
            strand_summary is not None and strand_correction_usable(strand_summary)
        )
        ssc = (
            float(strand_summary.signed_strand_contrast) if strand_summary is not None else 0.0
        )
        identifiable = _strand_identifiable_rows(region_arrays.strand)

        # Intron.
        in_raw = payload_arrays.intron_by_orient.sum(axis=1).astype(np.float64, copy=False)
        if strand_active:
            same = payload_arrays.intron_by_orient[:, ORIENT_SAME]
            opp = payload_arrays.intron_by_orient[:, ORIENT_OPP]
            corrected = _gdna_count_moment(same, opp, signed_strand_contrast=ssc)
            in_corrected = np.where(identifiable, np.maximum(corrected, 0.0), in_raw)
        else:
            in_corrected = in_raw
        Y[is_intron] = in_corrected[is_intron]

        # Exon boundary (raw is per-region: u_left * 1_L + u_right * 1_R).
        ul = payload_arrays.u_left_by_orient
        ur = payload_arrays.u_right_by_orient
        ex_raw = (
            eligible_left * ul.sum(axis=1) + eligible_right * ur.sum(axis=1)
        ).astype(np.float64, copy=False)
        if strand_active:
            same_ex = (
                eligible_left * ul[:, ORIENT_SAME] + eligible_right * ur[:, ORIENT_SAME]
            )
            opp_ex = (
                eligible_left * ul[:, ORIENT_OPP] + eligible_right * ur[:, ORIENT_OPP]
            )
            corrected_ex = _gdna_count_moment(same_ex, opp_ex, signed_strand_contrast=ssc)
            ex_corrected = np.where(identifiable, np.maximum(corrected_ex, 0.0), ex_raw)
        else:
            ex_corrected = ex_raw
        Y[is_exon] = ex_corrected[is_exon]

        # Per-class rho_global + kappa.
        class_density: dict[str, GlobalGdnaDensity] = {
            _CLASS_INTERGENIC: global_densities.intergenic,
            _CLASS_INTRON: global_densities.intron,
            _CLASS_EXON: global_densities.exon_intron,
        }
        class_mask = {
            _CLASS_INTERGENIC: is_intergenic,
            _CLASS_INTRON: is_intron,
            _CLASS_EXON: is_exon,
        }

        # EB-shrunk per-region rho_hat.
        rho_hat = np.zeros(R, dtype=np.float64)
        per_class_summary: dict[str, dict[str, float]] = {}
        signal_per_class: dict[str, float] = {}
        rho_ref_per_class: dict[str, float] = {}
        for cname in _CLASS_ORDER:
            d = class_density[cname]
            mask = class_mask[cname]
            if not mask.any():
                per_class_summary[cname] = _empty_class_summary(d)
                signal_per_class[cname] = 0.0
                rho_ref_per_class[cname] = float(d.rho)
                continue
            rho_global = float(d.rho)
            kappa = float(d.kappa.value)
            Yc = Y[mask]
            Ec = E[mask]
            denom = Ec + kappa
            with np.errstate(invalid="ignore", divide="ignore"):
                rho_c = np.where(denom > 0.0, (Yc + kappa * rho_global) / denom, rho_global)
            rho_hat[mask] = rho_c

            # Class reference (weighted Q95).
            rho_ref_c = _weighted_quantile(rho_c, Ec, REFERENCE_QUANTILE, fallback=rho_global)
            rho_ref_per_class[cname] = float(rho_ref_c)

            # Auto-uniform signal.
            if rho_global <= 0.0 or Ec.sum() <= 0.0:
                obs_spread = 0.0
                null_spread = 0.0
                signal = 0.0
                q05 = q50 = q95 = float(rho_global)
            else:
                log_floor = float(np.log(rho_global)) - LOG_RHO_CLIP_NATS
                log_rho = np.log(np.maximum(rho_c, max(np.exp(log_floor), np.exp(LOG_RHO_FLOOR))))
                q50_log = _weighted_quantile(log_rho, Ec, 0.5, fallback=float(np.log(rho_global)))
                q95_log = _weighted_quantile(log_rho, Ec, 0.95, fallback=float(np.log(rho_global)))
                obs_spread = float(max(q95_log - q50_log, 0.0))
                # Delta-method null variance.
                with np.errstate(invalid="ignore", divide="ignore"):
                    var_log = Ec / (rho_global * (Ec + kappa) ** 2)
                w_total = float(Ec.sum())
                sigma2 = float((var_log * Ec).sum() / w_total) if w_total > 0.0 else 0.0
                null_spread = Z_Q95 * float(np.sqrt(max(sigma2, 0.0)))
                signal = float(
                    np.clip(
                        (obs_spread - null_spread) / max(obs_spread, SPREAD_EPS), 0.0, 1.0
                    )
                )
                q05 = float(_weighted_quantile(rho_c, Ec, 0.05, fallback=rho_global))
                q50 = float(_weighted_quantile(rho_c, Ec, 0.50, fallback=rho_global))
                q95 = float(_weighted_quantile(rho_c, Ec, 0.95, fallback=rho_global))
            signal_per_class[cname] = signal
            per_class_summary[cname] = {
                "rho_global": float(rho_global),
                "rho_q05": q05,
                "rho_q50": q50,
                "rho_q95": q95,
                "rho_ref_class": float(rho_ref_c),
                "observed_log_spread": float(obs_spread),
                "null_log_spread": float(null_spread),
                "signal": float(signal),
                "kappa": float(kappa),
                "n_regions": int(mask.sum()),
            }

        rho_ref = max(rho_ref_per_class.values()) if rho_ref_per_class else 0.0
        if rho_ref <= 0.0:
            out = cls.uniform(region_arrays)
            return cls(
                rho_hat=rho_hat,
                log_weight=out.log_weight,
                weight=out.weight,
                mode="uniform",
                rho_ref=0.0,
                n_at_floor=0,
                per_class=per_class_summary,
                ref_offsets=out.ref_offsets,
                ref_id=out.ref_id,
                start=out.start,
                end=out.end,
            )

        # If every class signal is zero, no attenuation: return uniform but
        # preserve the diagnostic per_class summary.
        if all(s <= 0.0 for s in signal_per_class.values()):
            out = cls.uniform(region_arrays)
            return cls(
                rho_hat=rho_hat,
                log_weight=out.log_weight,
                weight=out.weight,
                mode="uniform",
                rho_ref=float(rho_ref),
                n_at_floor=0,
                per_class=per_class_summary,
                ref_offsets=out.ref_offsets,
                ref_id=out.ref_id,
                start=out.start,
                end=out.end,
            )

        # Per-class signal-attenuated log A_r.
        log_rho_ref = float(np.log(rho_ref))
        log_weight = np.zeros(R, dtype=np.float64)
        for cname in _CLASS_ORDER:
            mask = class_mask[cname]
            if not mask.any():
                continue
            signal = signal_per_class[cname]
            if signal <= 0.0:
                continue
            rho_c = rho_hat[mask]
            with np.errstate(invalid="ignore", divide="ignore"):
                raw = np.where(rho_c > 0.0, np.log(rho_c) - log_rho_ref, LOG_A_FLOOR)
            raw = np.minimum(raw, 0.0)
            log_weight[mask] = np.maximum(signal * raw, LOG_A_FLOOR)
        log_weight = np.maximum(log_weight, LOG_A_FLOOR)
        n_at_floor = int(np.sum(log_weight <= LOG_A_FLOOR + 1e-15))
        weight = np.exp(log_weight)

        return cls(
            rho_hat=rho_hat,
            log_weight=log_weight,
            weight=weight,
            mode="regional",
            rho_ref=float(rho_ref),
            n_at_floor=n_at_floor,
            per_class=per_class_summary,
            ref_offsets=np.asarray(region_arrays.ref_offsets, dtype=np.int32).copy(),
            ref_id=np.asarray(region_arrays.ref_id, dtype=np.int32).copy(),
            start=np.asarray(region_arrays.start, dtype=np.int64).copy(),
            end=np.asarray(region_arrays.end, dtype=np.int64).copy(),
        )

    # -- lookup -------------------------------------------------------------

    def log_weights_for_positions(
        self,
        ref_ids: np.ndarray,
        positions: np.ndarray,
    ) -> np.ndarray:
        """Vectorized ``log A(ref, pos)``; uniform mode returns zeros."""
        N = int(ref_ids.size)
        out = np.zeros(N, dtype=np.float64)
        if self.mode == "uniform" or N == 0:
            return out
        ref_ids = np.asarray(ref_ids, dtype=np.int32)
        positions = np.asarray(positions, dtype=np.int64)
        n_refs = int(self.ref_offsets.size - 1) if self.ref_offsets.size else 0
        order = np.argsort(ref_ids, kind="stable")
        ref_sorted = ref_ids[order]
        pos_sorted = positions[order]
        # Group boundaries by ref id.
        unique_refs, group_starts = np.unique(ref_sorted, return_index=True)
        group_ends = np.r_[group_starts[1:], N]
        for ref, g0, g1 in zip(unique_refs.tolist(), group_starts.tolist(), group_ends.tolist()):
            if ref < 0 or ref >= n_refs:
                continue  # leave at identity
            r_lo = int(self.ref_offsets[ref])
            r_hi = int(self.ref_offsets[ref + 1])
            if r_hi <= r_lo:
                continue
            ref_starts = self.start[r_lo:r_hi]
            ref_ends = self.end[r_lo:r_hi]
            qpos = pos_sorted[g0:g1]
            # searchsorted: rightmost start <= pos.
            idx = np.searchsorted(ref_starts, qpos, side="right") - 1
            valid = idx >= 0
            if not valid.any():
                continue
            hit = valid.copy()
            hit_idx = np.where(valid, idx, 0)
            hit &= qpos < ref_ends[hit_idx]
            if not hit.any():
                continue
            local = hit_idx[hit] + r_lo
            out_slice = np.zeros(qpos.size, dtype=np.float64)
            out_slice[hit] = self.log_weight[local]
            out[order[g0:g1]] = out_slice
        return out

    def weighted_length_on_ref(self, ref_id: int, start: int, end: int) -> float:
        """Integral of ``A(x)`` over ``[start, end)`` on ``ref_id``.

        Region partitions are gap-free by construction for any contig that
        has at least one region. A gap inside ``[start, end)`` therefore
        indicates a corrupted ``RegionArrays`` and raises ``RuntimeError``.
        Contigs with zero regions (e.g., FASTA reference absent from the
        GTF) return identity exposure ``end - start``.
        """
        if end <= start:
            return 0.0
        if self.mode == "uniform":
            return float(end - start)
        n_refs = int(self.ref_offsets.size - 1) if self.ref_offsets.size else 0
        if ref_id < 0 or ref_id >= n_refs:
            return float(end - start)
        r_lo = int(self.ref_offsets[ref_id])
        r_hi = int(self.ref_offsets[ref_id + 1])
        if r_hi <= r_lo:
            return float(end - start)
        ref_starts = self.start[r_lo:r_hi]
        ref_ends = self.end[r_lo:r_hi]
        # First region containing/after `start`.
        i = int(np.searchsorted(ref_starts, start, side="right") - 1)
        if i < 0:
            i = 0
        total = 0.0
        cursor = int(start)
        end_int = int(end)
        while cursor < end_int and i < ref_starts.size:
            rs = int(ref_starts[i])
            re = int(ref_ends[i])
            if rs > cursor:
                # Gap detected inside a covered contig: data invariant
                # violation, surface it loudly rather than silently
                # weighting at 1.0.
                raise RuntimeError(
                    f"RegionalGdnaExposure: gap-free invariant violated on "
                    f"ref_id={ref_id}: cursor={cursor} < region_start={rs} "
                    f"(query=[{start}, {end}))."
                )
            seg_end = min(re, end_int)
            seg_start = max(rs, cursor)
            if seg_end > seg_start:
                total += float(seg_end - seg_start) * float(self.weight[r_lo + i])
            cursor = re
            i += 1
        if cursor < end_int:
            # Query extends past the last region; treat the tail as identity
            # (this is the "contig has regions but query spills beyond end"
            # case — most likely a coordinate-clipping bug upstream, but
            # behaviorally equivalent to no information).
            total += float(end_int - cursor)
        return total

    # -- summary -----------------------------------------------------------

    def to_summary_dict(self) -> dict:
        return {
            "mode": self.mode,
            "rho_ref": float(self.rho_ref),
            "n_regions": int(self.rho_hat.size),
            "n_at_floor": int(self.n_at_floor),
            "log_a_floor": float(LOG_A_FLOOR),
            "per_class": {
                cname: {k: float(v) for k, v in summary.items()}
                for cname, summary in self.per_class.items()
            },
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    q: float,
    *,
    fallback: float,
) -> float:
    """Lower-step weighted quantile.

    Returns the smallest value whose cumulative weight reaches ``q * total``.
    Falls back to ``fallback`` when total weight is zero.
    """
    v = np.asarray(values, dtype=np.float64).ravel()
    w = np.asarray(weights, dtype=np.float64).ravel()
    keep = np.isfinite(v) & (w > 0.0)
    if not keep.any():
        return float(fallback)
    v = v[keep]
    w = w[keep]
    order = np.argsort(v, kind="stable")
    v = v[order]
    w = w[order]
    cw = np.cumsum(w)
    target = q * cw[-1]
    idx = int(np.searchsorted(cw, target, side="left"))
    if idx >= v.size:
        idx = v.size - 1
    return float(v[idx])


def _empty_class_summary(d: GlobalGdnaDensity) -> dict[str, float]:
    rho = float(d.rho)
    return {
        "rho_global": rho,
        "rho_q05": rho,
        "rho_q50": rho,
        "rho_q95": rho,
        "rho_ref_class": rho,
        "observed_log_spread": 0.0,
        "null_log_spread": 0.0,
        "signal": 0.0,
        "kappa": float(d.kappa.value),
        "n_regions": 0,
    }
