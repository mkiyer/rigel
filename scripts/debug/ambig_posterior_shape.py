"""Test the multi-root hypothesis: is an AMBIG node's strand posterior over f_g multi-modal?

Reconstructs the τ-marginal λ-posterior P(f_g) (exactly as _solve_ambig_logodds) from a node's raw
unspliced strand counts, with NO prior and NO messages — the pure strand likelihood. Prints the curve
so we can see whether a balanced node has two peaks (f_g≈0 gDNA-free two-strand-RNA, f_g≈1 gDNA) or is
flat/unimodal. Also sweeps synthetic balanced counts at several depths.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_posterior_shape.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.calibration.simplex import _mixture_strand_loglik  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
index, blob = build_or_load_cache(COND, False)
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
kappa = float(cal.rna_sense_frac)
od_g = float(cal.gdna_strand_overdispersion)
od_r = float(cal.rna_strand_overdispersion)
print(f"kappa(rna_sense_frac)={kappa:.4f}  od_g={od_g:.4f}  od_r={od_r:.4f}")

sub_f = CalibrationSubstrate.from_payload(blob["payload_full"], ra)
upos_all = np.asarray(sub_f.contained.n_unspliced_pos, float)
uneg_all = np.asarray(sub_f.contained.n_unspliced_neg, float)


def posterior_fg(upos, uneg, n_grid=121, n_tilt=81):
    """τ-marginal λ-posterior over f_g (strand-only, + λ Jacobian), replicating _solve_ambig_logodds."""
    N = upos + uneg
    lam = np.linspace(-8, 8, n_grid)
    fg = 1.0 / (1.0 + np.exp(-lam))              # (K,)
    tau = np.linspace(-1.0, 1.0, n_tilt)          # (Kt,)
    fgc = fg[:, None]                              # (K,1)
    fpos = (1.0 - fgc) * (1.0 + tau[None, :]) / 2.0   # (K,Kt)
    fneg = (1.0 - fgc) * (1.0 - tau[None, :]) / 2.0
    psi = _mixture_strand_loglik(np.array([[upos]]), N, fg[None, :, None],
                                 fpos[None], fneg[None], kappa, od_g, od_r)[0]  # (K,Kt)
    log_jac = np.log(np.maximum(fg * (1.0 - fg), 1e-12))    # (K,)
    full = psi + log_jac[:, None]
    # marginalize tau (logsumexp), then normalize over lam
    m = full.max(axis=1)
    marg = m + np.log(np.exp(full - m[:, None]).sum(axis=1))   # (K,)
    marg = marg - marg.max()
    P = np.exp(marg); P = P / P.sum()
    return fg, P


def describe(fg, P, tag):
    # modes: local maxima
    modes = [i for i in range(1, len(P) - 1) if P[i] > P[i - 1] and P[i] >= P[i + 1]]
    cdf = np.cumsum(P); med = fg[np.searchsorted(cdf, 0.5)]
    mean = float((fg * P).sum())
    mode_str = ", ".join(f"f_g={fg[i]:.3f}(p={P[i]:.3f})" for i in sorted(modes, key=lambda i: -P[i])[:3])
    # coarse histogram of mass in f_g bins
    bins = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.01]
    hist = [P[(fg >= bins[b]) & (fg < bins[b + 1])].sum() for b in range(len(bins) - 1)]
    hstr = " ".join(f"[{bins[b]:.1f}-{bins[b+1]:.1f}]={hist[b]:.2f}" for b in range(len(hist)))
    print(f"{tag}: median={med:.3f} mean={mean:.3f}  n_modes={len(modes)}  peaks[{mode_str}]")
    print(f"     mass: {hstr}")


print("\n=== real AMBIG-exon nodes (strand-only posterior over f_g) ===")
for reg in [242, 224, 672, 1080, 1083, 231]:
    up, un = upos_all[reg], uneg_all[reg]
    fg, P = posterior_fg(up, un)
    describe(fg, P, f"reg {reg} (upos={up:.0f} uneg={un:.0f}, bal={up/(up+un):.3f})")

print("\n=== synthetic perfectly-balanced nodes at increasing depth ===")
for N in [50, 500, 5000, 20000]:
    fg, P = posterior_fg(N / 2, N / 2)
    describe(fg, P, f"balanced N={N}")

print("\n=== synthetic strand-skewed nodes (stranded RNA signature) ===")
for frac in [0.6, 0.8, 0.95]:
    N = 20000
    fg, P = posterior_fg(N * frac, N * (1 - frac))
    describe(fg, P, f"skew pos_frac={frac} N={N}")
