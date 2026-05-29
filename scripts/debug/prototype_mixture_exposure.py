"""Prototype: mixture-responsibility per-region exposure vs. gate vs. no-gate.

Dumps per-region inputs from the paralog scenario at several (gdna, seed) settings,
then runs three estimators offline and compares their `omega` outputs.

Estimators:
  A) MAD + no gate (current Phase 5 + MAD, the regression we're trying to fix)
  B) MAD + p_unx >= 0.80 hard gate (Phase 4a baseline, the working version)
  C) Soft mixture: gamma_r = p(background | log_raw, fit) via single E-step,
     used as the per-region weight on log(omega). No threshold anywhere.

Each estimator returns omega_r per region. We then evaluate:
  - max |omega - 1| over "borderline" regions (the failure mode)
  - agreement on actual paralog regions (p_unx ~= 0 -> omega == 1)
  - agreement on clean background regions (p_unx ~= 1 -> matches gate behavior)
"""
from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "tests" / "scenarios_aligned"))
from conftest import (  # noqa: E402
    GDNAConfig,
    PIPELINE_SEED,
    ReadSimConfig,
    SIM_SEED,
    run_pipeline,
)

from rigel.config import BamScanConfig, EMConfig, PipelineConfig  # noqa: E402
from rigel.sim import Scenario  # noqa: E402

logging.basicConfig(level=logging.WARNING)
log = logging.getLogger("mixture-proto")
log.setLevel(logging.INFO)

_MAD_TO_SIGMA = 1.4826


@dataclass
class RegionInputs:
    """Inputs to the exposure estimator, one row per region."""

    region_bp: np.ndarray
    gdna_mass: np.ndarray
    rna_mass: np.ndarray
    p_unx: np.ndarray
    fused_info: np.ndarray
    support: np.ndarray  # unspliced_counts
    log_raw_ratio: np.ndarray
    v_obs: np.ndarray
    rho0: float


def dump_inputs(gdna: int, seed: int, work_root: Path) -> RegionInputs:
    work_dir = work_root / f"g{gdna}_s{seed}"
    if work_dir.exists():
        import shutil
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    sc = Scenario("paralog_proto", genome_length=12000, seed=SIM_SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1000)], "abundance": 100}])
    sc.add_gene("g2", "+", [{"t_id": "t2", "exons": [(5000, 5500)], "abundance": 100}])
    sc.genome.edit(5000, sc.genome[500:1000])
    sc.add_gene("g_helper", "+", [
        {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50},
    ])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0}])

    sim = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                       read_length=100, strand_specificity=1.0, seed=SIM_SEED)
    gdna_cfg = (GDNAConfig(abundance=gdna, frag_mean=350, frag_std=100,
                          frag_min=100, frag_max=1000) if gdna > 0 else None)
    res = sc.build(n_fragments=3000, sim_config=sim, gdna_config=gdna_cfg, nrna_abundance=0)
    cfg = PipelineConfig(em=EMConfig(seed=seed),
                        scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True))
    pr = run_pipeline(res.bam_path, res.index, config=cfg)

    rc = pr.calibration.region_calibration
    exp = rc.region_exposure
    mass = rc.region_unspliced_mass
    return RegionInputs(
        region_bp=np.asarray(mass.region_size_bp, dtype=np.float64),
        gdna_mass=np.asarray(mass.gdna_mass, dtype=np.float64),
        rna_mass=np.asarray(mass.rna_mass, dtype=np.float64),
        p_unx=np.asarray(rc.p_unexpressed, dtype=np.float64),
        fused_info=np.asarray(mass.density_information + mass.strand_information,
                              dtype=np.float64),
        support=np.asarray(mass.unspliced_counts, dtype=np.float64),
        log_raw_ratio=np.asarray(exp.log_raw_ratio, dtype=np.float64),
        v_obs=np.asarray(exp.v_obs, dtype=np.float64),
        rho0=float(exp.rho0),
    )


def _weighted_mad(values: np.ndarray, weights: np.ndarray) -> float:
    w = np.asarray(weights, dtype=np.float64)
    mask = w > 0
    if not np.any(mask):
        return 0.0
    abs_dev = np.abs(values[mask])
    order = np.argsort(abs_dev, kind="stable")
    cum = np.cumsum(w[mask][order], dtype=np.float64)
    if cum[-1] <= 0:
        return 0.0
    idx = int(np.searchsorted(cum, 0.5 * cum[-1], side="left"))
    return float(abs_dev[order][min(idx, abs_dev.size - 1)])


def estimator_A_mad_no_gate(rin: RegionInputs) -> np.ndarray:
    """Current Phase 5 + MAD: weight = fused_info * p_unx, no gate."""
    weights = rin.fused_info * rin.p_unx
    weights = np.where(rin.support > 0, weights, 0.0)
    mad = _weighted_mad(rin.log_raw_ratio, weights)
    tau2 = max((_MAD_TO_SIGMA * mad) ** 2, 0.0)
    if tau2 <= 0:
        return np.ones_like(rin.log_raw_ratio)
    shrink = tau2 / (tau2 + rin.v_obs)
    log_omega = shrink * rin.log_raw_ratio
    return np.exp(np.clip(log_omega, np.log(1e-6), np.log(1e6)))


def estimator_B_mad_gate(rin: RegionInputs, threshold: float = 0.80) -> np.ndarray:
    """MAD + hard gate at p_unx >= threshold (the working Phase 4a behavior)."""
    pool = (rin.support > 0) & (rin.p_unx >= threshold)
    weights = rin.fused_info * rin.p_unx * pool
    mad = _weighted_mad(rin.log_raw_ratio, weights)
    tau2 = max((_MAD_TO_SIGMA * mad) ** 2, 0.0)
    if tau2 <= 0:
        return np.ones_like(rin.log_raw_ratio)
    shrink = tau2 / (tau2 + rin.v_obs)
    shrink = np.where(pool, shrink, 0.0)
    log_omega = shrink * rin.log_raw_ratio
    return np.exp(np.clip(log_omega, np.log(1e-6), np.log(1e6)))


def estimator_C_mixture(rin: RegionInputs, n_iter: int = 6) -> np.ndarray:
    """Soft mixture: log_raw ~ p_bg * N(0, tau2 + v_obs) + (1-p_bg) * Wide(0, sigma_c^2).

    The 'wide' contamination component is a broad zero-mean Gaussian with a
    pre-specified large variance — it acts as an outlier-absorber. Per-region
    responsibility gamma_r = p(background | log_raw, fit) is then used as both
    (a) the weight when re-estimating tau2 and (b) the multiplier on shrink_weight
    in the final omega. No threshold is used anywhere.

    Prior on the background-membership: pi_r = p_unx (from the EM model).
    """
    # Initialize tau2 from MAD on p_unx-weighted data, like estimator A.
    init_w = rin.fused_info * rin.p_unx
    init_w = np.where(rin.support > 0, init_w, 0.0)
    mad = _weighted_mad(rin.log_raw_ratio, init_w)
    tau2 = max((_MAD_TO_SIGMA * mad) ** 2, 1e-6)

    # Wide-tailed contamination component: fixed large variance.
    # Choose sigma_c^2 = 100 * tau2_init (well outside the bg distribution).
    sigma_c2 = max(100.0 * tau2, 25.0)  # log-ratio of ~5 -> ratio of ~150x

    # Only consider regions with positive support; others can't contribute.
    valid = (rin.support > 0) & (rin.fused_info > 0) & (rin.region_bp > 0)
    pi_prior = np.clip(rin.p_unx, 1e-6, 1.0 - 1e-6)

    x = rin.log_raw_ratio
    # We treat v_obs as the per-region sampling noise variance under the bg model.
    # log_raw | background ~ N(0, tau2 + v_obs).  log_raw | contamination ~ N(0, sigma_c2).
    gamma = np.where(valid, pi_prior, 0.0)
    for _ in range(n_iter):
        var_bg = tau2 + rin.v_obs
        # Gaussian densities (unnormalized const cancels in responsibility).
        log_p_bg = -0.5 * (np.log(var_bg) + (x * x) / var_bg)
        log_p_c = -0.5 * (np.log(sigma_c2) + (x * x) / sigma_c2)
        # Posterior responsibility for background.
        log_num = np.log(pi_prior) + log_p_bg
        log_den_c = np.log(1.0 - pi_prior) + log_p_c
        m = np.maximum(log_num, log_den_c)
        gamma = np.exp(log_num - m) / (np.exp(log_num - m) + np.exp(log_den_c - m))
        gamma = np.where(valid, gamma, 0.0)
        # M-step: re-estimate tau2 from gamma-weighted MAD.
        w = rin.fused_info * gamma
        mad = _weighted_mad(x, w)
        tau2_new = max((_MAD_TO_SIGMA * mad) ** 2, 1e-6)
        if abs(tau2_new - tau2) / max(tau2, 1e-6) < 1e-4:
            tau2 = tau2_new
            break
        tau2 = tau2_new

    # Per-region omega: weight the linear shrinkage by gamma_r (responsibility).
    shrink = gamma * tau2 / (tau2 + rin.v_obs)
    log_omega = shrink * x
    return np.exp(np.clip(log_omega, np.log(1e-6), np.log(1e6)))


def evaluate(gdna: int, seed: int, rin: RegionInputs) -> dict:
    omega_A = estimator_A_mad_no_gate(rin)
    omega_B = estimator_B_mad_gate(rin)
    omega_C = estimator_C_mixture(rin)

    # Categorize regions by p_unx tier.
    paralog = (rin.p_unx < 0.05) & (rin.rna_mass > 0)
    background = (rin.p_unx >= 0.80) & (rin.rna_mass == 0)
    borderline = (rin.p_unx >= 0.05) & (rin.p_unx < 0.80)

    def _dev(o: np.ndarray, mask: np.ndarray) -> float:
        if not np.any(mask):
            return 0.0
        return float(np.max(np.abs(np.log(o[mask]))))

    return {
        "gdna": gdna,
        "seed": seed,
        "n_regions": int(rin.region_bp.size),
        "n_paralog": int(paralog.sum()),
        "n_border": int(borderline.sum()),
        "n_bg": int(background.sum()),
        # Max |log omega| in each category (smaller = closer to omega=1 = less risk)
        "A_para_devlog": _dev(omega_A, paralog),
        "B_para_devlog": _dev(omega_B, paralog),
        "C_para_devlog": _dev(omega_C, paralog),
        "A_border_devlog": _dev(omega_A, borderline),
        "B_border_devlog": _dev(omega_B, borderline),
        "C_border_devlog": _dev(omega_C, borderline),
        "A_bg_devlog": _dev(omega_A, background),
        "B_bg_devlog": _dev(omega_B, background),
        "C_bg_devlog": _dev(omega_C, background),
    }


def print_per_region(gdna: int, seed: int, rin: RegionInputs) -> None:
    omega_A = estimator_A_mad_no_gate(rin)
    omega_B = estimator_B_mad_gate(rin)
    omega_C = estimator_C_mixture(rin)
    print(f"\n=== gdna={gdna} seed={seed}: per-region omega comparison ===")
    print(f"{'idx':>3} {'bp':>5} {'gdna':>6} {'rna':>5} {'p_unx':>6} {'log_raw':>8} "
          f"{'v_obs':>9} {'omega_A':>9} {'omega_B':>9} {'omega_C':>9}")
    for i in range(rin.region_bp.size):
        print(f"{i:>3d} {rin.region_bp[i]:>5.0f} {rin.gdna_mass[i]:>6.1f} "
              f"{rin.rna_mass[i]:>5.1f} {rin.p_unx[i]:>6.3f} "
              f"{rin.log_raw_ratio[i]:>8.3f} {rin.v_obs[i]:>9.3g} "
              f"{omega_A[i]:>9.3g} {omega_B[i]:>9.3g} {omega_C[i]:>9.3g}")


def main() -> None:
    work_root = Path("/tmp/mixture_proto")
    work_root.mkdir(exist_ok=True)

    cases = [(0, PIPELINE_SEED), (20, PIPELINE_SEED), (20, 1), (20, 7),
             (100, PIPELINE_SEED), (100, 1)]
    rows = []
    for gdna, seed in cases:
        try:
            rin = dump_inputs(gdna, seed, work_root)
            row = evaluate(gdna, seed, rin)
            rows.append(row)
            log.info(
                "gdna=%3d seed=%3d  para_devlog A=%.2f B=%.2f C=%.2f  "
                "border_devlog A=%.2f B=%.2f C=%.2f  bg_devlog A=%.2f B=%.2f C=%.2f",
                gdna, seed,
                row["A_para_devlog"], row["B_para_devlog"], row["C_para_devlog"],
                row["A_border_devlog"], row["B_border_devlog"], row["C_border_devlog"],
                row["A_bg_devlog"], row["B_bg_devlog"], row["C_bg_devlog"],
            )
        except Exception as exc:
            log.error("gdna=%d seed=%d failed: %r", gdna, seed, exc)
            raise

    # Dump per-region view for the most informative case.
    rin = dump_inputs(20, PIPELINE_SEED, work_root)
    print_per_region(20, PIPELINE_SEED, rin)
    rin = dump_inputs(100, PIPELINE_SEED, work_root)
    print_per_region(100, PIPELINE_SEED, rin)


if __name__ == "__main__":
    main()
