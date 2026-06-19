"""Phase 0 — bimodal-node characterization (does the trough-mean projection P2 actually bite?).

Constructs a SINGLE balanced-count AMBIG node off-chain, runs the CURRENT `_local_loglik`, and characterizes
the per-component lattice marginals: are they bimodal/ridge-like? where do the MEAN (the message/relay value
for f± today, via `_axis_mean`) and the MEDIAN (the f_g readout, via `_fg_median`) land relative to the
posterior mass? This establishes the baseline the Phase-1 message/readout operator split must beat, and tests
the adversarial P2 concern empirically on our code.

    OMP_NUM_THREADS=1 python scripts/debug/bimodal_node.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

from rigel.calibration.simplex import _simplex_lattice
from rigel.calibration.simplex_sweep import _local_loglik, _fg_median, _axis_mean, _fg_var
from scipy.special import logsumexp

N_GRID = 60


def _marginal(psi_row, axis_g):
    """Posterior marginal over the unique lattice values of one axis: returns (values, prob)."""
    post = np.exp(psi_row - logsumexp(psi_row))
    vals = np.round(axis_g, 6)
    uniq = np.unique(vals)
    prob = np.array([post[vals == v].sum() for v in uniq])
    return uniq, prob


def _modality(vals, prob):
    """Crude bimodality flag: count interior local maxima of the (smoothed) marginal above 5% of the peak."""
    p = prob / max(prob.max(), 1e-12)
    peaks = 0
    for i in range(1, len(p) - 1):
        if p[i] >= p[i - 1] and p[i] >= p[i + 1] and p[i] > 0.05:
            peaks += 1
    # also count edge peaks (mass piled at a vertex)
    edge = (p[0] > 0.5) + (p[-1] > 0.5)
    return peaks, int(edge)


def characterize(label, *, u_pos, u_neg, kappa, od_g, od_r, strand_obs, global_mu, global_tau):
    fp_g, fn_g, fg_g = _simplex_lattice(N_GRID)
    lattice = (fp_g, fn_g, fg_g)

    def arr(x):
        return np.array([float(x)])

    psi = _local_loglik(
        arr(u_pos), arr(u_neg), arr(0.0), arr(0.0), arr(True).astype(bool), arr(True).astype(bool),
        float(kappa), float(od_g), float(od_r), lattice,
        strand_obs=np.array([strand_obs]),
        global_mu=(None if global_mu is None else arr(global_mu)),
        global_tau=(0.0 if global_mu is None else float(global_tau)),
    )
    row = psi[0]
    med_g = _fg_median(psi, fg_g)[0]
    mean_g = _axis_mean(psi, fg_g)[0]
    mean_p = _axis_mean(psi, fp_g)[0]
    mean_n = _axis_mean(psi, fn_g)[0]
    var_g = _fg_var(psi, fg_g)[0]
    print(f"\n  [{label}]  u=({u_pos:.0f},{u_neg:.0f}) kappa={kappa} od=({od_g},{od_r}) "
          f"strand_obs={strand_obs} global_mu={global_mu}")
    for nm, axis in (("f_g", fg_g), ("f_pos", fp_g), ("f_neg", fn_g)):
        vals, prob = _marginal(row, axis)
        peaks, edge = _modality(vals, prob)
        # report mass at the two ends + the interior trough (min in the middle third)
        mid = prob[(vals > 0.33) & (vals < 0.67)]
        shape = "BIMODAL/edge-piled" if (edge >= 1 or peaks >= 2) else "unimodal"
        print(f"      {nm:>5}: {shape:>18}  mass[0]={prob[0]:.3f} mass[1]={prob[-1]:.3f} "
              f"mid-mass={mid.sum():.3f}")
    print(f"      => f_g: median={med_g:.3f} mean={mean_g:.3f} var={var_g:.4f}  | "
          f"MESSAGE f_pos(mean)={mean_p:.3f} f_neg(mean)={mean_n:.3f}")
    # the trough test: is the MEAN message landing where posterior mass is LOW?
    for nm, axis, m in (("f_pos", fp_g, mean_p), ("f_neg", fn_g, mean_n)):
        vals, prob = _marginal(row, axis)
        i = int(np.argmin(np.abs(vals - m)))
        peak = prob.max()
        if prob[i] < 0.5 * peak and (prob[0] > peak * 0.8 or prob[-1] > peak * 0.8):
            print(f"      !! TROUGH-MEAN on {nm}: mean={m:.3f} sits at posterior prob {prob[i]:.3f} "
                  f"(peak {peak:.3f}) — message propagates a state the data rejects")


if __name__ == "__main__":
    print("=" * 96)
    print("BIMODAL-NODE CHARACTERIZATION — balanced AMBIG node, current _local_loglik (P2 test)")
    print("=" * 96)
    # 1) UNSTRANDED balanced AMBIG, zero-gDNA global: the flat-strand ridge case
    characterize("unstranded, global_mu=0", u_pos=100, u_neg=100, kappa=0.5,
                 od_g=0.0, od_r=0.0, strand_obs=False, global_mu=0.0, global_tau=10.0)
    # 2) STRANDED balanced AMBIG (counts balanced despite kappa=0.99 => genuinely 50/50 sense/antisense)
    characterize("stranded kappa=.99, global_mu=0", u_pos=100, u_neg=100, kappa=0.99,
                 od_g=0.0, od_r=0.0, strand_obs=False, global_mu=0.0, global_tau=10.0)
    # 3) STRANDED with overdispersion (realistic)
    characterize("stranded kappa=.99 od=.1, global_mu=0", u_pos=100, u_neg=100, kappa=0.99,
                 od_g=0.1, od_r=0.1, strand_obs=False, global_mu=0.0, global_tau=10.0)
    # 4) high-gDNA global (B-step destination flavor): does the mean track the global or trough?
    characterize("stranded kappa=.99, global_mu=0.8", u_pos=100, u_neg=100, kappa=0.99,
                 od_g=0.0, od_r=0.0, strand_obs=False, global_mu=0.8, global_tau=10.0)
    # 5) strongly + node (NOT balanced) as a control — should be unimodal, mean≈median
    characterize("control: pos-dominant 190:10", u_pos=190, u_neg=10, kappa=0.99,
                 od_g=0.0, od_r=0.0, strand_obs=False, global_mu=0.0, global_tau=10.0)
