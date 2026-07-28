"""MC validation of the P1e CONSERVATION SURPRISE laws — the same shape as `message_variance_mc.py`.

Everything is validated against the ESTIMATOR'S ERROR VS TRUTH, never against a posterior's own spread.

A message delivers per-component densities rho_c with the destination's effective lengths E_c, so it ASSERTS
S = sum_c rho_c E_c against the observed mass M.  The delivered MODE (what psi consumes) is
mo_c = log(rho_c E_c / M), whose error vs truth is exactly the component's log error eps_c, so the whole
statistic can be simulated in mode space.

    P1e-1  THE NULL.  delta = log(M/S).  On a matched reframe the message is built as
           rho_c^msg = rho_c^true e^{eps_c} * rho_tot(dst)/rho_tot(src), and rho_tot(dst) = M*B(f_g) carries
           M EXACTLY -- so S is proportional to M and the destination's Poisson count CANCELS from delta:
                Var(delta) = alpha^T Sigma alpha       (NOT + 1/n_dst)
           The count only survives to the extent that the destination's face density is NOT proportional to
           M, i.e. through its SPLICED part: with q = rho_spl/(M*B),
                d delta / d log M = q/(1+q)   =>   Var(delta) = alpha^T Sigma alpha + (q/(1+q))^2 / n_dst.
           P1e-1a is the region-destination case (q = 0), P1e-1b the boundary case.

    P1e-2  DL RECOVERY against the corrected null: b2 = max(0, delta^2 - alpha^T Sigma alpha - nu) is an
           unbiased method-of-moments estimate of the unmodelled excess E[(alpha^T eps_extra)^2], with the
           same +0.4839*(sd2+nu) positive bias at b = 0 that M7 has (the OVER-damping, safe direction).

    P1e-3  THE INFLATION IS HONEST.  The observation is ONE scalar, so the coherent attribution of the
           excess is along the model's own conditional direction s = Sigma alpha (E[eps | alpha^T eps = t]
           = t*s/sd2), i.e. Sigma' = Sigma + b2 * s s^T / sd2^2, giving alpha^T Sigma' alpha = delta^2
           exactly and, per component, Delta v_c = b2 * s_c^2 / sd2^2 and
           Delta Var(lambda) = b2 * (s_g - s_R)^2 / sd2^2.
           (a) a COMMON-SCALE unmodelled error must be priced onto the LEVELS and must NOT reach lambda;
           (b) a COMPONENT-SPECIFIC unmodelled error must reach lambda -- but only if Sigma's own diag(w)
               part is non-degenerate.  With w = 0 the model cannot attribute it and stays OVER-CONFIDENT on
               lambda: that failure is asserted, not hidden.

    P1e-4  THE PIN IS INFORMATION.  The weighted rescale a = delta*s/sd2 is exactly the Gaussian conditional
           mean given alpha^T eps = -delta, and the post-pin error covariance is exactly the conditional
           (Schur) covariance Sigma - s s^T / sd2.  The COMMON-factor pin (today's `_pin_v`) is the
           projection (I - 1 alpha^T) and its residual covariance is larger in the Loewner order.  This is
           why the DEDUCTION must NOT be taken when the null is in doubt: taking it under a false null
           leaves the message OVER-CONFIDENT.

    OMP_NUM_THREADS=1 python scratchpad/p1e_mc.py [--draws N]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260726)
FAILS = []


def _report(name, pred, emp, tol=0.05):
    pred, emp = float(pred), float(emp)
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    ok = rel <= tol
    if not ok:
        FAILS.append(name)
    print(f"  {'OK ' if ok else '***'} {name:<62} pred {pred:12.6g}   emp {emp:12.6g}   rel {rel:7.2%}")
    return ok


def _assert(name, cond, note=""):
    if not cond:
        FAILS.append(name)
    print(f"  {'OK ' if cond else '***'} {name:<62} {note}")
    return cond


def _sigma(s_cm2, w):
    """Sigma = s_cm^2 * 11^T + diag(w) — the shipped `conservation_rescale` error model."""
    C = len(w)
    return s_cm2 * np.ones((C, C)) + np.diag(np.asarray(w, float))


def _draw(Sigma, N):
    L = np.linalg.cholesky(Sigma + 1e-15 * np.eye(Sigma.shape[0]))
    return rng.standard_normal((N, Sigma.shape[0])) @ L.T


# ─────────────────────────────────────────────────────────────────────────────────────────────────────
def null_region(N, *, f_true, E, s_cm2, w, n_dst):
    """P1e-1a — REGION destination: the reframe carries M exactly, so delta is M-FREE.

    Var(delta) = alpha^T Sigma alpha, and the sketch's ``+ 1/n_dst`` is a subtraction of a noise that is
    NOT THERE — which UNDER-states b2 and therefore UNDER-damps (the expensive direction)."""
    Sigma = _sigma(s_cm2, w)
    eps = _draw(Sigma, N)
    alpha0 = np.asarray(f_true, float)                      # the true mass shares = alpha at eps = 0
    M = rng.poisson(n_dst, N).astype(float) / n_dst         # the destination's own Poisson mass, relative
    # the claim: rho_c^msg = rho_c^true e^{eps_c} * (M/M_true)   (the reframe's rho_tot(dst) = M*B)
    m = alpha0[None, :] * np.exp(eps) * M[:, None]          # per-component asserted mass / M_true
    S = m.sum(axis=1)
    delta = np.log(M) - np.log(S)
    alpha = m / S[:, None]
    sd2 = float(np.mean(np.einsum("nc,cd,nd->n", alpha, Sigma, alpha)))
    ok = _report("P1e-1a  Var(delta) = alpha^T Sigma alpha   (region dst, M cancels)",
                 sd2, float(np.var(delta)), tol=0.05)
    bad = sd2 + 1.0 / n_dst
    print(f"      the sketch's +1/n_dst would predict {bad:.6g} -> z2 = {np.var(delta) / bad:5.3f} "
          f"(<1 = under-charges b2 = UNDER-damps)")
    del E
    return ok


def null_boundary(N, *, f_true, s_cm2, w, n_dst, q):
    """P1e-1b — BOUNDARY destination: the face density is M*B + rho_spl, so only the fraction 1/(1+q) of
    log M reaches S and the surviving count leg is (q/(1+q))^2 / n_dst — NOT 1/n_dst."""
    Sigma = _sigma(s_cm2, w)
    eps = _draw(Sigma, N)
    alpha0 = np.asarray(f_true, float)
    M = rng.poisson(n_dst, N).astype(float) / n_dst
    scale = (M + q) / (1.0 + q)                            # (M*B + rho_spl)/(M_true*B + rho_spl)
    m = alpha0[None, :] * np.exp(eps) * scale[:, None]
    S = m.sum(axis=1)
    delta = np.log(M) - np.log(S)
    alpha = m / S[:, None]
    sd2 = float(np.mean(np.einsum("nc,cd,nd->n", alpha, Sigma, alpha)))
    c = q / (1.0 + q)
    pred = sd2 + c * c / n_dst
    ok = _report(f"P1e-1b  Var(delta) = a^T Sigma a + (q/(1+q))^2/n_dst   [q={q:g}]",
                 pred, float(np.var(delta)), tol=0.06)
    print(f"      c = d delta/d log M = {c:.4f};  the flat +1/n_dst would predict "
          f"{sd2 + 1.0 / n_dst:.6g} (z2 = {np.var(delta) / (sd2 + 1.0 / n_dst):5.3f})")
    return ok


# ─────────────────────────────────────────────────────────────────────────────────────────────────────
def dl_recovery(N, *, f_true, s_cm2, w, n_dst, q, b_extra, direction):
    """P1e-2 — the DerSimonian-Laird recovery of the unmodelled excess against the CORRECTED null."""
    Sigma = _sigma(s_cm2, w)
    eps = _draw(Sigma, N)
    u = np.asarray(direction, float)
    xi = rng.normal(0.0, b_extra, N) if b_extra > 0 else np.zeros(N)
    eps = eps + xi[:, None] * u[None, :]
    alpha0 = np.asarray(f_true, float)
    M = rng.poisson(n_dst, N).astype(float) / n_dst
    scale = (M + q) / (1.0 + q)
    m = alpha0[None, :] * np.exp(eps) * scale[:, None]
    S = m.sum(axis=1)
    delta = np.log(M) - np.log(S)
    alpha = m / S[:, None]
    sd2 = np.einsum("nc,cd,nd->n", alpha, Sigma, alpha)
    nu = (q / (1.0 + q)) ** 2 / n_dst
    b2 = np.maximum(0.0, delta**2 - sd2 - nu)
    truth = b_extra**2 * float(np.mean((alpha @ u) ** 2))
    if b_extra > 0:
        # the linearisation delta ~ -alpha^T eps is second-order accurate, so a large extra OVER-recovers
        # by a few percent -- the OVER-damping (safe) direction.
        return _report(f"P1e-2   E[b2] recovers (alpha^T u)^2 * b^2   [b={b_extra:g}, u={direction}]",
                       truth, float(np.mean(b2)), tol=0.08)
    pred = 0.4839 * float(np.mean(sd2 + nu))
    return _report("P1e-2   b = 0 floor = 0.4839*(sd2+nu)  (the SAFE positive bias)",
                   pred, float(np.mean(b2)), tol=0.06)


# ─────────────────────────────────────────────────────────────────────────────────────────────────────
def inflation_honest(N, *, f_true, s_cm2, w, n_dst, b_extra, direction, tag):
    """P1e-3 — is the inflation HONEST about the delivered mode error, and along WHICH direction?

    One scalar observation identifies ONE parameter, so the direction is a modelling assertion.  Three
    candidates, all satisfying alpha^T Sigma' alpha = delta^2 exactly:
        RANK1_s   Sigma + b2 * s s^T / sd2^2         s = Sigma alpha (the Gaussian conditional direction)
        COMMON    Sigma + b2 * 1 1^T / (alpha^T 1)^2 = Sigma + b2 * 11^T   (a pure SCALE error)
        BIRGE     (delta^2/sd2) * Sigma              (the quasi-likelihood / over-dispersion rescale)
    z2_c = E[eps_c^2]/E[v_c'] per LEVEL and z2_lam = E[(eps_g-eps_R)^2]/E[Var(lambda)'].
    z2 = 1 honest, < 1 conservative, > 1 OVER-CONFIDENT (the failure the project cares about)."""
    Sigma = _sigma(s_cm2, w)
    eps0 = _draw(Sigma, N)
    u = np.asarray(direction, float)
    xi = rng.normal(0.0, b_extra, N)
    eps = eps0 + xi[:, None] * u[None, :]
    alpha0 = np.asarray(f_true, float)
    M = rng.poisson(n_dst, N).astype(float) / n_dst
    m = alpha0[None, :] * np.exp(eps) * M[:, None]
    S = m.sum(axis=1)
    delta = np.log(M) - np.log(S)
    alpha = m / S[:, None]
    s_vec = alpha @ Sigma                                     # s = Sigma alpha, per draw
    sd2 = np.einsum("nc,nc->n", alpha, s_vec)
    b2 = np.maximum(0.0, delta**2 - sd2)                      # q = 0 here => nu = 0
    one = np.ones_like(s_vec)
    a1 = alpha.sum(axis=1)
    err2 = [float(np.mean(eps[:, c] ** 2)) for c in (0, 1)]
    err_lam = float(np.mean((eps[:, 0] - eps[:, 1]) ** 2))
    v_lam0 = float(w[0] + w[1])                               # the common part cancels from lambda
    print(f"    [{tag}]  extra b={b_extra:g} along {direction};  E[b2] = {np.mean(b2):.4g};  "
          f"true E[err^2] g/R/lam = {err2[0]:.4g}/{err2[1]:.4g}/{err_lam:.4g}")
    ok = True
    for law, dv, dlam in (
        ("RANK1_s", b2[:, None] * s_vec**2 / sd2[:, None] ** 2,
         b2 * (s_vec[:, 0] - s_vec[:, 1]) ** 2 / sd2**2),
        ("COMMON ", b2[:, None] * one / a1[:, None] ** 2, np.zeros_like(b2)),
        ("BIRGE  ", (b2 / sd2)[:, None] * np.array([Sigma[0, 0], Sigma[1, 1]])[None, :],
         (b2 / sd2) * v_lam0),
    ):
        z = []
        for c in (0, 1):
            z.append(err2[c] / (float(Sigma[c, c]) + float(np.mean(dv[:, c]))))
        zl = err_lam / (v_lam0 + float(np.mean(dlam))) if (v_lam0 + np.mean(dlam)) > 0 else float("inf")
        good = max(z[0], z[1], zl) <= 1.15
        ok &= _assert(f"P1e-3 {tag} :: {law}", good,
                      f"z2  gDNA {z[0]:7.2f}   RNA {z[1]:7.2f}   SPLIT lambda {zl:10.2f}")
    return ok


# ─────────────────────────────────────────────────────────────────────────────────────────────────────
def pin_is_information(N, *, f_true, s_cm2, w, n_dst):
    """P1e-4 — the WEIGHTED pin is the Gaussian conditional mean and its residual covariance is exactly
    the Schur complement Sigma - s s^T/sd2.  The COMMON-factor pin is (I - 1 alpha^T), strictly worse."""
    Sigma = _sigma(s_cm2, w)
    eps = _draw(Sigma, N)
    alpha0 = np.asarray(f_true, float)
    M = rng.poisson(n_dst, N).astype(float) / n_dst
    m = alpha0[None, :] * np.exp(eps) * M[:, None]
    S = m.sum(axis=1)
    delta = np.log(M) - np.log(S)
    alpha = m / S[:, None]
    s_vec = alpha @ Sigma
    sd2 = np.einsum("nc,nc->n", alpha, s_vec)
    # mode errors: mo_c - truth = eps_c (M cancels from the mode), plus the pin's correction
    e_w = eps + delta[:, None] * s_vec / sd2[:, None]      # weighted pin
    e_c = eps + delta[:, None]                            # common-factor pin
    s0 = Sigma @ alpha0
    sd20 = float(alpha0 @ s0)
    post = Sigma - np.outer(s0, s0) / sd20
    P = np.eye(len(w)) - np.outer(np.ones(len(w)), alpha0)
    post_c = P @ Sigma @ P.T
    ok = True
    for c, nm in ((0, "gDNA"), (1, "RNA")):
        ok &= _report(f"P1e-4  weighted pin residual Var({nm}) = (Sigma - ss^T/sd2)_cc",
                      post[c, c], float(np.var(e_w[:, c])), tol=0.06)
    for c, nm in ((0, "gDNA"), (1, "RNA")):
        ok &= _report(f"P1e-4  common pin  residual Var({nm}) = (P Sigma P^T)_cc",
                      post_c[c, c], float(np.var(e_c[:, c])), tol=0.06)
    ok &= _assert("P1e-4  the weighted pin is <= the common pin on every component",
                  bool(np.all(np.diag(post) <= np.diag(post_c) + 1e-12)),
                  f"diag {np.round(np.diag(post), 5)} vs {np.round(np.diag(post_c), 5)}")
    # taking the DEDUCTION when the null is FALSE leaves the message over-confident
    xi = rng.normal(0.0, 0.5, N)
    eps2 = eps + xi[:, None] * np.array([1.0, 1.0])[None, :]
    m2 = alpha0[None, :] * np.exp(eps2) * M[:, None]
    S2 = m2.sum(axis=1)
    d2 = np.log(M) - np.log(S2)
    a2 = m2 / S2[:, None]
    s2 = a2 @ Sigma
    sd2b = np.einsum("nc,nc->n", a2, s2)
    e2 = eps2 + d2[:, None] * s2 / sd2b[:, None]
    z = float(np.var(e2[:, 0])) / post[0, 0]
    _assert("P1e-4  under a FALSE null the deduction leaves z2 > 1 (=> do not deduct)", z > 1.2,
            f"z2(gDNA, deducted) = {z:.2f}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# MC validation of the P1e CONSERVATION SURPRISE   draws={N:,}\n")
    f2 = [0.62, 0.38]  # a gDNA-rich exon claim: alpha_g, alpha_R

    print("═══ P1e-1  THE NULL: what noise is actually in delta ═══\n")
    null_region(N, f_true=f2, E=[2100.0, 1850.0], s_cm2=0.02, w=[0.05, 0.004], n_dst=400)
    print()
    null_region(N, f_true=f2, E=[2100.0, 1850.0], s_cm2=0.0, w=[0.30, 0.004], n_dst=60)
    print()
    for q in (0.5, 3.0):
        null_boundary(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=120, q=q)
    print()

    print("\n═══ P1e-2  DL recovery of the unmodelled excess, against the CORRECTED null ═══\n")
    dl_recovery(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400, q=0.0,
                b_extra=0.0, direction=[1.0, 1.0])
    dl_recovery(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400, q=0.0,
                b_extra=0.60, direction=[1.0, 1.0])
    dl_recovery(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400, q=0.0,
                b_extra=0.60, direction=[1.0, 0.0])
    dl_recovery(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=120, q=1.5,
                b_extra=0.60, direction=[1.0, 0.0])

    print("\n\n═══ P1e-3  is the rank-1 inflation HONEST about the delivered mode error? ═══\n")
    inflation_honest(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400,
                     b_extra=0.60, direction=[1.0, 1.0], tag="COMMON-scale extra")
    print()
    inflation_honest(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400,
                     b_extra=0.60, direction=[1.0, 0.0], tag="gDNA-ONLY extra")
    print()
    inflation_honest(N, f_true=f2, s_cm2=0.0, w=[0.30, 0.004], n_dst=400,
                     b_extra=0.60, direction=[1.0, 0.0], tag="gDNA-ONLY, w-dominant Sigma")
    print()
    inflation_honest(N, f_true=f2, s_cm2=0.30, w=[1e-9, 1e-9], n_dst=400,
                     b_extra=0.60, direction=[1.0, 0.0], tag="gDNA-ONLY, DEGENERATE Sigma (w=0)")

    print("\n\n═══ P1e-4  the pin is INFORMATION: the conditional covariance ═══\n")
    pin_is_information(N, f_true=f2, s_cm2=0.02, w=[0.05, 0.004], n_dst=400)

    print(f"\n\n{'ALL PASS' if not FAILS else str(len(FAILS)) + ' FAILURES: ' + ', '.join(FAILS)}")


if __name__ == "__main__":
    main()
