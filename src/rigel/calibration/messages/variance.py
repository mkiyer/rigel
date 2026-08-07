"""The message POLICY's variance toolbox — pure arithmetic, no chain and no solver state.

       Gate: ``tests/calibration/test_enrichment_frame.py``

⭐ This was ``calibration/enrichment_frame.py``. It moved under ``messages/`` because every function in it
prices a step a MESSAGE takes, so it belongs to the policy rather than to the backbone — the backbone owns
the loop, the combine, the write-back and the five assertions, and none of them needs a variance.

⚠ **ONE function crosses the layer and it is deliberate**: :func:`count_logvar` is also imported by
``node_init``, which builds the per-slot self-solve the backbone hands the policy as its seed. That is the
Poisson log-count variance and it has exactly ONE home (TRAPS: a-test-that-redefines) — the same ``1/n`` had been written
out twice, and unifying it is what the panel-wide **39 %** belongs to. It stays here, imported.


Hybrid capture multiplies *every* component of a node's density by the same efficiency ``e(x)``, so it is
visible in the **total density and only in the total** (a composition *ratio* ``k = ρ_g/ρ_R`` cancels ``e``).
A message carries a per-base density, so it is transportable across an enrichment cliff only after being
scaled into the destination node's frame by the *measured* ratio of total densities. This module owns that
arithmetic and nothing else: it is a pure function of masses, compositions, and effective lengths — no chain,
no substrate, no solver state — so its facts are pinned by closed-form unit tests
(``tests/calibration/test_enrichment_frame.py``) that cannot drift with solver behaviour.

 (the framework, §1 the total-density
pivot, §2 the bounding lemma, §4 the step-wise junction solve) and (the
junction instance + §5b the r₂/r₁ asymmetry). The solver wiring that consumes these is phase TRAPS: annotated-is-not-genomic (behind a
flag); this module is frame-agnostic by construction (arrays in, arrays out), so the per-face-vs-node-level
frame question TRAPS: annotated-is-not-genomic must settle cannot invalidate anything here.

The design principle behind the API shape: **transport ``k``, never ``f_g``.** ``k`` is enrichment-free and
component-set-shared; ``f_g`` depends on the node's own crossing-vs-contained effective lengths, so copying it
across a frame change is a known past defect (§2 of the junction doc). And the composition assumption is
carried as a **variance**, not a bool gate: the bounding lemma (§2) says a totally wrong composition still
pins ``ρ_tot`` to within the effective-length ratio (1.04–1.5× for normal nodes, 4×+ for short regions), which
is a continuous quantity — :func:`composition_logvar` derives it with no tuned threshold, so a short region
down-weights itself automatically, with no ``L ≲ fl_mean`` cut anywhere in the module.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import log_ndtr, polygamma, zeta

__all__ = [
    "composition_logvar",
    # ── the pass-0 message-VARIANCE laws ──
    "graft_frame_logvar",
    "graft_frame_logvar_scalar",
    "peel_rna_logvar",
    "peel_continue_share",
    "peel_continue_share_scalar",
    "peel_share_logvar",
    "residual_level",
    "residual_level_scalar",
    "transfer_logvar",
    "mismatch_gap",
    "mismatch_deflate",
]

_EPS = 1.0e-12
_INV_EPS = 1.0 / _EPS
_HALF_LOG_2PI = 0.5 * np.log(2.0 * np.pi)  # folded once; the same bits as the inline expression


def _f(x):
    return np.asarray(x, dtype=np.float64)


# ── scalar twins of the three numpy selectors, for the SCALAR path below ───────────────────────────────
# They exist because `np.maximum`/`np.minimum` PROPAGATE nan (`(a > b || isnan(a)) ? a : b`) where a bare
# `x if x > c else c` would silently swallow it, and because `np.clip` KEEPS x inside the interval (so it
# preserves −0.0) rather than composing max-then-min. Both differences are reachable, both are bugs if
# mismatched, and both are one token wide. Bitwise-verified against numpy over every (x, constant) pair
# this module forms, including ±inf / nan / ±0.0 (`test_enrichment_frame.py`).
def _fmax(x, c):
    return x if x > c or x != x else c


def _fmin(x, c):
    return x if x < c or x != x else c


def _fclip01(x):
    return x if x != x else (0.0 if x < 0.0 else (1.0 if x > 1.0 else x))


def _exp(x):
    """``np.exp`` on a scalar — **NOT** ``math.exp``, and the difference is load-bearing.

    ⚠ ``np.exp`` and ``math.exp`` are NOT bit-identical, because ``exp`` is not correctly rounded and the
    two use different implementations: ``math.exp`` is the platform libm, while numpy carries its own
    AVX2/AVX-512 kernel on x86-64. They agree on arm64 (where numpy's float64 ``exp`` *is* libm) and differ
    by 1 ULP on x86-64 Linux — so a scalar twin written with ``math.exp`` passes its bit-identity test on a
    Mac and fails it in CI. That is exactly what happened (`test_residual_level_scalar_...`, corner
    ``(1.0, 1e-13, 0.4, 210.0, 260.0, 1e6)``).

    1 ULP is not a rounding detail here: :func:`residual_level` subtracts two nearly-equal normal pdfs
    (``pdf_a − pdf_b``, magnification ~1.2e4 at that corner) and then forms ``1 − φ + σ·d`` (a second
    cancellation, ``−83.0 + 83.4999…``), so one ULP in ``exp`` moves the result by ~1e6 ULP.

    Using ``np.exp`` on both sides makes the twins identical **by construction on every platform**: numpy's
    ``exp`` is size-independent (bulk == per-element == 0-d == Python float, verified over 1e5 draws), so
    the array form's vectorised call and this scalar call return the same bits.

    ⚠ **Any transcendental added to a scalar twin must come from numpy, not ``math``**, for this reason.
    ``sqrt`` is exempt — IEEE-754 requires it correctly rounded, so the two agree by standard."""
    return float(np.exp(x))


def _log(x):
    """``np.log`` on a scalar — see :func:`_exp`. Same hazard, same reason, same rule."""
    return float(np.log(x))


def count_logvar(count) -> np.ndarray:
    """``Var(log ρ)`` for a Poisson rate seen as ``count`` events over an opportunity — **exactly**, at
    every count including zero. ⭐⭐ **THE ONE HOME**, imported by every consumer.

    A rate seen as ``a`` events over exposure ``E``, under the SAME Jeffreys prior ψ is built on
    (`simplex_logodds._JEFFREYS_REF` derives it as Jeffreys for a Poisson rate), has posterior
    ``Gamma(a + ½, E)`` whose log has variance ``trigamma(a + ½)`` — independent of ``E``, since the
    opportunity moves the location and cannot sharpen the claim.

    This IS the ``1/n`` it replaces: ``trigamma(a + ½) → 1/a``, agreeing to <0.1 % for ``a ≥ 10``. The
    whole difference is at small counts, and at ``a = 0`` it is ``π²/2 = 4.93`` (sd 2.22 nats) instead of
    ``∞``.

    ⛔⛔ **THE SAME ASYMPTOTE WAS WRITTEN OUT TWICE AND BOTH COPIES BROKE THE SAME WAY** (TRAPS: a-zero-count-is-a-measurement,
    TRAPS: a-ratio-cannot-carry-zero). In ``own_precision`` it made a zero-count object emit nothing; in :func:`composition_logvar` it
    made ``Var(log ρ_tot) = ∞``, hence ``σ²_transfer = ∞``, hence ``1/(1/p + ∞) = 0`` — annihilating every
    message that object sent, on all three streams. Fixing one and not the other is why the repair moved
    the refit solve and not pass-0. ⭐ It lives HERE, in the leaf module, so there is one definition and
    nothing to keep in step."""
    return polygamma(1, np.asarray(count, np.float64) + 0.5)


def composition_logvar(f_g, E_g, E_r, var_fg, n):
    """``Var(log ρ_tot)`` — the composition assumption carried as a **variance**, not a bool gate.

    Two independent sources of uncertainty in a total density (the
    bounding lemma, made continuous):

    * **counting** — ``M`` is a Poisson count, so ``Var(log M) = trigamma(M + ½)``
      (:func:`count_logvar`; ``→ 1/n``, but FINITE at ``n = 0``, which ``1/n`` was not);
    * **composition** — ``ρ_tot = M·B(f_g)`` with ``B = f_g/E_g + (1−f_g)/E_r``, and ``dB/df_g = 1/E_g − 1/E_r``,
      so by the delta method ``Var(log ρ_tot | composition) = [(1/E_g − 1/E_r)/B]²·Var(f_g)``.

    ::

        Var(log ρ_tot) = trigamma(M + ½)  +  [ (1/E_g − 1/E_r) / B ]² · Var(f_g)

    **No tuned constant.** The bracket is the bounding lemma's factor expressed as a coefficient: it is ≈0.036
    for a long contained region (composition barely matters), ≈0.4 for a boundary crossing (E_g≠E_r), and
    diverges as ``E_g`` collapses on a short region — so a short region down-weights *itself* as an enrichment
    reference exactly where §2 says it must be excluded, with no ``L ≲ fl_mean`` threshold. At the structural
    ``f_g = 1`` corner (intergenic / seam, ``Var(f_g) = 0``) it reduces to the pure counting precision ``1/n``,
    which is why ``M/E_g`` and the blended estimator are bit-identical there (E0 finding 2)."""
    fg, eg, er = _f(f_g), _f(E_g), _f(E_r)
    vfg, nn = _f(var_fg), _f(n)
    inv_g, inv_r = 1.0 / np.maximum(eg, _EPS), 1.0 / np.maximum(er, _EPS)
    B = fg * inv_g + (1.0 - fg) * inv_r
    coeff = (inv_g - inv_r) / np.maximum(B, _EPS)
    comp = coeff * coeff * np.maximum(vfg, 0.0)
    return count_logvar(nn) + comp


def graft_frame_logvar(r):
    """the-graft-frame-variance — the GRAFT's **frame-mislift** log-variance: ``(log r)²`` on the MEASURED spliced component.

    the-reframe-scale-variance sets ``σ²_transfer = 0`` on a graft edge because ``r`` is common-mode across the matched set ``{g, R}``
    and cancels in the composition. That is true of the IMPUTED continue ``ρ_ν``, which travels with ``ρ_g``
    from the same source — but it is **false of the grafted ``ρ_μ``**. The delivered share is

        f_g : f_R  =  ρ_g^src·E_g^dst : (ρ_ν^src + ρ_μ)·E_r^dst

    (``r`` and the ``_pin_v`` factor multiply every component alike and cancel exactly — verified to 1.8e-15
    against the shipped arrays). So ``ρ_μ`` is being ratioed against the **source's** gDNA density. But a
    spliced fragment's two blocks lie in the flanking EXONS, so ``ρ_μ = S/E_spl`` is measured in the
    DESTINATION exon's frame, not the boundary's — structurally, the same kind of geometric fact as the
    ``free_pos``/``free_neg`` strand gates. Measured on the suite (oracle densities, no solver in the loop):
    ``ρ_R(exon)/ρ_spl(bnd)`` = 1.02–1.86 and is capture-INVARIANT, while the exon↔boundary gDNA step goes
    1.03 → 6.1–6.8 under capture.

    So the grafted component has **no matched gDNA partner to cancel ``r`` against** — exactly the-reframe-scale-variance's peel /
    partial-anchor case, where ``σ²_transfer`` is load-bearing. The un-modelled frame step is what the message
    gets wrong, its magnitude is the step itself, and the model's own estimate of that step is ``r``; charging
    ``(log r)²`` is the method-of-moments second moment of a single observation of it — the same logic as the-mismatch-deflation's
    ``b̂²``, and with no tuned constant. It is identically 0 when there is no frame change (``r = 1``), which is
    why the shipped model is EXACT off-capture (measured required correction ``log c`` = +0.009/−0.008/+0.054)
    and only bites where capture opens a step.

    ⚠ This is NOT the retired ``σ²_cliff = (log r)²`` proxy, which charged the whole cliff on **matched**
    reframes — where the composition genuinely IS preserved across the enrichment step, so it over-damped.
    Here the component set is *unmatched by construction*, so the cliff really is un-cancelled error.

    Applied to the spliced measurement's OWN precision (``1/n_spl ⊕ (log r)²``), so the-composition-variance's share weighting
    ``w_μ²`` arises implicitly from the inverse-variance fusion with the correctly-framed ``ρ_ν`` arm — a
    graft that is a minority of the RNA is damped proportionately less.

    Calibration (delivered ``λ`` vs oracle, graft edges into exons): ``z2 = E[Δ²]/E[v]`` 58–310 → **2.1–3.8**,
    consistently across capture off and on."""
    rr = _f(r)
    lr = np.log(np.where(rr > _EPS, rr, 1.0))
    return lr * lr


def graft_frame_logvar_scalar(r):
    """Scalar twin of :func:`graft_frame_logvar` — see that docstring for the model.

    ⚠ TWIN: mirror any change into both. Pinned bit-for-bit by ``test_enrichment_frame.py``."""
    lr = _log(r) if r > _EPS else 0.0  # r ≤ _EPS (or nan) ⇒ no frame step ⇒ log 1 = 0
    return lr * lr


def peel_rna_logvar(v_log_rho_R, s2_transfer, v_mu, u):
    """the-peel-share-variance — the PEEL (exon→boundary) RNA-continue message log-variance. The boundary receives only what
    CONTINUES: ``ρ_ν = ρ_R(x)/r − ρ_μ`` — a **difference** (an absolute measured density is subtracted, so
    enrichment does NOT cancel). The delta method gives u-weighted (≥1) terms::

        Var(log ρ_ν) = u²·Var(log T) + (u−1)²·v_μ ,   Var(log T) = Var(log ρ_R(x)) + σ²_transfer ,
        u = T/ρ_ν = 1/(fraction that continues) ≥ 1

    A difference DESTROYS precision (subtracting near-equal numbers) — the mirror of the graft's convex weights;
    ``σ²_transfer`` is LOAD-BEARING here (~85–92% of the variance).

    ⚠ The linearization is valid only for ``ε = √(fraction continuing) ≲ 0.15`` (``u ≲ 3``); beyond it (>p75 of
    real junctions) it UNDER-states the variance (over-confident 27–40%). ``u`` MUST therefore gate the peel's
    precision as a per-junction weight, and ``ρ_ν < 0`` is a PRIOR truncation, not an emission gate (handoff §6).
    MC 1–3% in-regime."""
    uu = np.asarray(u, np.float64)
    vT = np.asarray(v_log_rho_R, np.float64) + np.asarray(s2_transfer, np.float64)
    return uu * uu * vT + (uu - 1.0) ** 2 * np.asarray(v_mu, np.float64)


def peel_continue_share(rho_nu, rho_mu):
    """the-continuing-share — the fraction of a seam's RNA that CONTINUES unspliced: ``w = ρ_ν/(ρ_ν + ρ_μ)``.

    This is the object that retires the peel's SUBTRACTION. What continues past a junction is a *share* of the
    RNA at the seam, and a share is **enrichment-free**: capture multiplies the continuing and the splicing
    channels alike (``ρ_ν = e·c_ν``, ``ρ_μ = e·c_μ``), so ``e`` cancels identically inside ``w``. Both inputs
    are taken in the BOUNDARY's own frame — its solved unspliced-RNA density and its measured spliced density.

    **Why a share and not a difference.** The old peel formed ``ρ_ν = ρ_R(x)·r − ρ_μ``, a difference of two
    absolute densities, and a difference does not commute with a scale error. With ``u = A/ρ_ν ≥ 1`` (the
    reciprocal of the continuing fraction), a systematic log-scale error ``δ`` in the reframed source arrives
    as ``log(u·e^δ − (u−1))`` — which is ``u·δ`` in the small-``δ`` limit and is measured at 1.77× / 2.39× /
    5.01× for ``u = 2 / 3 / 10`` at ``δ = 0.30``. Through this share it arrives as ``δ``, EXACTLY, for every
    ``u`` and every ``δ``: a scaling commutes with a scaling. That matters because the exon-facing reframe
    error is irreducible — the boundary samples a ``fl_mean`` window around a point while the exon samples its
    interior, so with mid-exon probes the two sit at genuinely different capture (measured 0.4–1.3 nats), and
    ``u``'s p75 on real junctions is ≈ 3. MC: `message_variance_mc.py` the-continuing-share, exact to 1e-12.

    Degenerate limits, both structural: no spliced flux (``ρ_μ = 0``) ⇒ ``w = 1``, nothing splices away and
    nothing is peeled; no RNA at the seam at all (``ρ_ν + ρ_μ = 0``) ⇒ ``w = 1``, there is nothing to
    apportion and the caller's own gates decide. ``w`` is always in ``[0, 1]`` — the peel can no longer go
    negative, which retires the zero-truncation defect (a fully-consumed peel used to emit ``ρ_ν = 0``, "no
    RNA continues past here", at a live precision)."""
    nu, mu = _f(rho_nu), _f(rho_mu)
    tot = nu + mu
    return np.where(tot > _EPS, nu / np.maximum(tot, _EPS), 1.0)


def peel_continue_share_scalar(rho_nu, rho_mu):
    """Scalar twin of :func:`peel_continue_share` — see that docstring for the model.

    ⚠ TWIN: mirror any change into both. Pinned bit-for-bit by ``test_enrichment_frame.py``."""
    tot = rho_nu + rho_mu
    return rho_nu / tot if tot > _EPS else 1.0


def peel_share_logvar(w_mu, v_nu, v_mu):
    """the-continuing-share — the log-variance the continuing SHARE contributes: ``w_μ²·(v_ν + v_μ)``.

    Delta method on ``log w = log ρ_ν − log(ρ_ν + ρ_μ)``: both partials are ``±w_μ``, the SPLICED share
    ``w_μ = 1 − w``, so::

        Var(log w) = w_μ²·( Var(log ρ_ν^B) + Var(log ρ_μ^B) ) ,   v_μ = 1/n_spl (measured, count-only)

    and the delivered message variance is ``Var(log ρ_R(x)) + σ²_transfer + Var(log w)``. The weights are
    **CONVEX** (``w_μ ≤ 1``) — the exact mirror of the-composition-variance's graft SUM — where the-peel-share-variance's difference carried ``u ≥ 1`` and
    amplified. So this move takes the peel out of the ill-conditioned regime entirely; there is no ``ε``/``u``
    validity limit to respect and no over-confidence tail. MC: `message_variance_mc.py` the-continuing-share, rel 0.2–1.0 %."""
    wm = _f(w_mu)
    return wm * wm * (_f(v_nu) + _f(v_mu))


def residual_level(mass, n_mass, rho_g, E_g, E_r, v_g):
    """the-residual-level — a node's own RNA density read off its **observed mass** against an imputed gDNA **density**.

    This is the generic DENSITY DECONVOLUTION (`density_deconv`, of which the intron factory is the intron
    special case) with the gDNA density prior supplied by a NEIGHBOUR instead of by the intergenic pool. It is
    the only gDNA/RNA split available at a seam with no factory within reach and no strand — i.e. at
    ``exon|exon`` boundaries (97 % have no factory neighbour) and at every seam of a low-gDNA library — and it
    is count-zero-information-legal for the same reason the factory is: the *information* is the imputed gDNA
    DENSITY, and the count only converts that density into a composition::

        φ   = ρ_g·E_g / M              the share of the crossing the imputed gDNA claim accounts for
        f_R = 1 − φ                    the RNA share of the crossing
        σ_f = φ·√v_g                   the gDNA claim's own log-variance, entering through the RATIO

    **The BOUNDS are applied on the FRACTION, not on the mass**, because that is where they are exact. Doing
    the truncation on the mass (``0 ≤ X ≤ M`` with ``X ~ N(M − ρ_g E_g, M²/n + (ρ_g E_g)²v_g)``) uses the
    OBSERVED ``M`` as a bound on a quantity whose true bound is ``M_true``, and that breaks the pure-RNA
    limit: at ``ρ_g = 0`` the answer is ``f_R = 1`` exactly — the count cancels out of the ratio — where the
    mass form leaves an ``M/√n`` residual pressed against the upper bound and under-claims the RNA by
    ``0.8/√n``. Here ``σ_f → 0`` with ``φ``, as it must.

    **``f_R`` is TRUNCATED TO ``[0, 1]``, not clipped**, and both bounds are structural: an RNA share cannot be
    negative, and it cannot exceed the whole crossing (the gDNA share it is subtracted from is itself
    non-negative). The estimator is the moment-matched two-sided truncated normal ``f_R | 0 ≤ f_R ≤ 1`` with
    ``f_R ~ N(1−φ, σ_f²)``; ``ρ_ν = E[f_R]·M/E_r``. Each bound earns its keep and neither is optional:

    * **the LOWER bound** is what lets the estimator say something at an RNA-free seam. A naive ``max(·, 0)``
      would report density 0 at infinite variance — "no opinion" — exactly where the node in fact holds a
      strong one: that its RNA is below its own noise floor.
    * **the UPPER bound** is what stops it inventing one. The imputed gDNA claim routinely arrives at
      ``√v_g ≈ 1.0–1.2 nats`` (measured at exon→boundary edges under capture), which makes ``σ_f`` of order 1;
      a one-sided positive part then returns ``E[f_R] ≈ 0.8·σ_f``, i.e. *"most of my mass is RNA"*, asserted
      **out of pure ignorance** and at a confident-looking ``k ≈ 2``. With the upper bound the same ignorance
      degrades to its correct limit — ``σ_f ≫ 1`` ⇒ ``f_R ~ Uniform(0,1)`` ⇒ ``E = ½`` at ``k = 3`` — a wide
      claim the fuse out-weighs with any real evidence, rather than a narrow wrong one that swamps it.

    **The log-variance is the TRIGAMMA of the level's effective fragment count**, not the delta method::

        k = E[f_R]² / Var[f_R]                    the RNA share's effective COUNT — how many fragments of RNA
        Var(log ρ_ν) = ψ'(k) + 1/(n·E[f_R]²)      the node has evidence for, in exactly the ``Var(log) = 1/n``
                                                  sense used everywhere else in the model, ⊕ the LEVEL's own
                                                  Poisson count

    The second term is the count, propagated with its exact Jacobian: ``ρ_ν·E_r = M − ρ_g E_g``, so
    ``∂log ρ_ν/∂log M = 1/f_R`` — a count error is amplified by the reciprocal of the RNA share, because the
    level is a small difference of two large numbers. Dropping that amplification (writing a bare ``+1/n``)
    under-states the variance by exactly the M/f_R covariance, MC-caught at 19 % on the ``f_R = 0.9`` arm. The
    linearization point is the TRUNCATED mean, not ``1−φ``, which is what keeps it finite as ``φ → 1``.

    i.e. moment-match a Gamma (the count-posterior family) and take its exact log-variance. This is EXACT in
    both limits where the answer is known and the delta method is not: ``k → ∞`` gives ``ψ'(k) → 1/k``, the
    delta/Poisson limit; ``k → 1`` (the deep sub-noise tail, where the truncated normal *is* an exponential)
    gives ``ψ'(1) = π²/6``, where the delta method returns 1 and is over-confident by 1.6×. **``k ≥ 1``
    always**, so the level can never carry more than one fragment's worth of unearned confidence — an exact
    limit of the truncation, not a floor.

    The three operating limits, and they are the reason this exists:
    * **``ρ_g E_g ≪ M``** (a low-gDNA library — where RNA is the entire signal and the peel's old
      no-evidence default silenced the channel outright): ``ρ_ν → M/E_r`` at ``Var → 1/n``. A MEASUREMENT.
    * **``ρ_g E_g → M`` with a precise gDNA claim** (an RNA-free seam): ``ρ_ν → 0``, and the log-variance
      grows without bound — but the LINEAR statement stays tight (``sd(ρ_ν) ∝ M/(E_r√n)``, i.e. *"below a few
      percent of my mass"*), which is exactly the information a near-zero level carries and exactly why the
      consumer fuses in linear space. This is what reproduces "intronic unspliced fragments are gDNA until
      proven otherwise" as a measurement rather than a rule.
    * **an imprecise gDNA claim** (``σ_f ≫ 1``): ``ρ_ν → M/(2E_r)`` at ``k = 3``. Uninformative, and declared
      uninformative — a wide claim the fuse out-weighs with any real evidence.

    Returns ``(ρ_ν, Var(log ρ_ν), Var(ρ_ν))`` — the LINEAR variance as well, because near zero that is the
    coordinate the information lives in: "0.02 ± 0.04" and "log ρ = −3.9 ± 2.0" are the same linearization,
    but only the first can be fused with another estimate in a way that lets a confident near-zero claim pull
    the level DOWN. A log-space (geometric) fuse of positive modes cannot reach zero at all.

    ``(0, +inf, +inf)`` where no level exists — no mass, no count, or a gDNA claim carrying no precision
    (``v_g = +inf``: "supplied" is a statement about PRECISION, never about the density's value — the same
    test the λ-emission gate makes)."""
    m_, n_, rg_ = _f(mass), _f(n_mass), _f(rho_g)
    eg_, er_, vg_ = _f(E_g), _f(E_r), _f(v_g)
    ok = (m_ > _EPS) & (n_ > 0.0) & (er_ > _EPS) & np.isfinite(vg_)
    inv_n = np.where(ok, 1.0 / np.maximum(n_, _EPS), 0.0)
    # substitute a finite v_g BEFORE any product — np.where evaluates both branches, so a masked-out
    # ``v_g = +inf`` against ``φ = 0`` would still form ``0·inf = nan`` (the standing trap).
    vg_fin = np.where(np.isfinite(vg_), np.maximum(vg_, 0.0), 0.0)
    phi = np.where(ok, np.maximum(rg_, 0.0) * eg_ / np.maximum(m_, _EPS), 0.0)
    # The claim's log-variance is linearized at the ADMISSIBLE point ``min(φ, 1)``, not at ``φ`` itself. ``φ``
    # is a share and cannot exceed 1; an estimate that does (routine — a relayed claim asserts more reads of a
    # component than the node sequenced on 52–71 % of nodes) is "at least all of it", and its relative
    # uncertainty applies to the largest admissible value. Linearizing at ``φ`` instead makes ``σ_f`` grow
    # WITH the over-claim, so the [0,1] window widens faster than the mean leaves it and the estimator INVERTS
    # — more imputed gDNA returning more RNA, reaching ``f_R = 1`` at ``φ ≈ 100``. This is the projection onto
    # the admissible set, not a clamp on the answer: the mean still uses the full ``φ``.
    sig = np.minimum(phi, 1.0) * np.sqrt(vg_fin)
    # the standardized bounds of ``0 ≤ f_R ≤ 1``: β = φ/σ ≥ 0 always, α = β − 1/σ. The two branches below
    # each subtract two SAME-SIDE normal tails, so neither loses precision to cancellation.
    tiny = sig <= _EPS
    ss = np.where(tiny, 1.0, sig)
    beta = np.where(tiny, 1.0, phi / ss)
    alpha = beta - 1.0 / ss
    pdf_a = np.exp(-0.5 * alpha * alpha - _HALF_LOG_2PI)
    pdf_b = np.exp(-0.5 * beta * beta - _HALF_LOG_2PI)
    Z = np.where(
        alpha >= 0.0,
        np.exp(log_ndtr(-alpha)) - np.exp(log_ndtr(-beta)),
        np.exp(log_ndtr(beta)) - np.exp(log_ndtr(alpha)),
    )
    # Two limits where ``Z`` underflows and the ratio would be 0/0. Both are supplied exactly:
    #   * ``α > 0`` — the whole untruncated mass sits BELOW the lower bound (the gDNA claim over-explains the
    #     crossing). On [0,1] the density is then an exponential with scale ``σ/α = σ²/(φ−1)``: ``E = V^½``,
    #     so ``k = 1`` — one fragment's worth of evidence, the least this estimator can ever claim.
    #   * otherwise — ``σ_f`` swamps the whole unit interval: ``f_R ~ Uniform(0,1)``, ``E = ½``, ``k = 3``.
    bad = ~(Z > _EPS) & ~tiny
    lo_tail = bad & (alpha > 0.0)
    e_tail = np.where(lo_tail, ss / np.where(lo_tail, alpha, 1.0), 0.0)
    Zs = np.where(Z > _EPS, Z, 1.0)
    d = (pdf_a - pdf_b) / Zs
    f_hi = np.clip(1.0 - phi, 0.0, 1.0)
    f_R = np.clip(
        np.where(
            tiny,
            f_hi,
            np.where(lo_tail, e_tail, np.where(bad, 0.5, 1.0 - phi + ss * d)),
        ),
        0.0,
        1.0,
    )
    v_f = np.where(
        tiny,
        0.0,
        np.where(
            lo_tail,
            e_tail * e_tail,
            np.where(
                bad,
                1.0 / 12.0,
                ss * ss * np.maximum(1.0 + (alpha * pdf_a - beta * pdf_b) / Zs - d * d, 0.0),
            ),
        ),
    )
    k = np.maximum(
        f_R * f_R / np.maximum(v_f, _EPS), 1.0
    )  # the RNA share's effective fragment count
    # ψ'(k) — as ``zeta(2, k)``, which is what ``scipy.special.polygamma(1, k)`` reduces to exactly (its
    # ``(−1)^(n+1)`` and ``Γ(n+1)`` prefactors are both 1.0 at n=1) without also evaluating ψ(k) and
    # discarding it. Bitwise-identical; the equality is pinned in `test_enrichment_frame.py`.
    v_log = np.where(
        ok,
        zeta(2.0, np.minimum(k, _INV_EPS))
        + inv_n / np.maximum(f_R * f_R, _EPS),  # the count, with its exact 1/f_R Jacobian
        np.inf,
    )
    rho = np.where(ok, f_R * m_ / np.maximum(er_, _EPS), 0.0)
    # mask the ∞ BEFORE the product (np.where evaluates both branches: ``0·inf = nan``, the standing trap).
    return rho, v_log, np.where(ok, rho * rho * np.where(ok, v_log, 0.0), np.inf)


def residual_level_scalar(mass, n_mass, rho_g, E_g, E_r, v_g):
    """Scalar twin of :func:`residual_level` — see that docstring for the model and the derivation.

    Same arithmetic in the same association order, with every ``np.where`` turned back into the branch it
    encodes. That is where the speed is: the array form must evaluate BOTH arms of all 21 selectors, so a
    node with no level still pays for ``log_ndtr`` ×4 and ``ζ(2,·)``, and every one of those ops costs
    0.5–0.7 µs on a 0-d array against ~0.02 µs for the float expression. 25× on the main branch, 223× on
    the no-level early-out — and the relay calls this once per node per direction, sequentially, which is
    the bulk of a genome-scale calibration.

    ⚠ TWIN of :func:`residual_level`: mirror any change into both. The equality is pinned bit-for-bit,
    over every branch and every ±inf/nan/±0.0 corner, by ``test_enrichment_frame.py``."""
    m_, n_, rg_ = float(mass), float(n_mass), float(rho_g)
    eg_, er_, vg_ = float(E_g), float(E_r), float(v_g)
    if not (m_ > _EPS and n_ > 0.0 and er_ > _EPS and math.isfinite(vg_)):
        return 0.0, math.inf, math.inf  # no level exists
    inv_n = 1.0 / _fmax(n_, _EPS)
    phi = _fmax(rg_, 0.0) * eg_ / _fmax(m_, _EPS)
    sig = _fmin(phi, 1.0) * math.sqrt(_fmax(vg_, 0.0))
    if sig <= _EPS:  # a certain gDNA claim ⇒ the share is exactly what it leaves over
        f_R, v_f = _fclip01(1.0 - phi), 0.0
    else:
        beta = phi / sig
        alpha = beta - 1.0 / sig
        if (
            alpha >= 0.0
        ):  # subtract SAME-SIDE tails — neither branch loses precision to cancellation
            Z = _exp(log_ndtr(-alpha)) - _exp(log_ndtr(-beta))
        else:
            Z = _exp(log_ndtr(beta)) - _exp(log_ndtr(alpha))
        if not Z > _EPS:  # Z underflowed — the two exact limits the docstring supplies
            if alpha > 0.0:  # the claim over-explains the crossing: an exponential on [0,1], k = 1
                e_tail = sig / alpha
                f_R, v_f = _fclip01(e_tail), e_tail * e_tail
            else:  # σ_f swamps the unit interval: f_R ~ Uniform(0,1), k = 3
                f_R, v_f = 0.5, 1.0 / 12.0
        else:
            pdf_a = _exp(-0.5 * alpha * alpha - _HALF_LOG_2PI)
            pdf_b = _exp(-0.5 * beta * beta - _HALF_LOG_2PI)
            d = (pdf_a - pdf_b) / Z
            f_R = _fclip01(1.0 - phi + sig * d)
            v_f = sig * sig * _fmax(1.0 + (alpha * pdf_a - beta * pdf_b) / Z - d * d, 0.0)
    fr2 = f_R * f_R
    k = _fmin(_fmax(fr2 / _fmax(v_f, _EPS), 1.0), _INV_EPS)
    v_log = float(zeta(2.0, k)) + inv_n / _fmax(fr2, _EPS)
    rho = f_R * m_ / _fmax(er_, _EPS)
    return rho, v_log, rho * rho * v_log


def transfer_logvar(logvar_tot_dst, logvar_tot_src, graft):
    """the-reframe-scale-variance — ``σ²_transfer = Var(log r)``, the enrichment-ratio uncertainty that damps a message across a capture
    cliff: ``Var(log r) = Var(log ρ_tot^dst) + Var(log ρ_tot^src)`` (each from :func:`composition_logvar`).

    DIRECTION-dependent: on the **graft** the reframe ``r`` is common-mode across the matched component set and
    CANCELS in the composition (return 0 — applying it there is the double-count the density-uniformity proxy
    committed); on the **peel** / partial-anchor message it is LOAD-BEARING. This one law replaces the retired
    ``var_proj[dst] + (μ_proj[dst]−μ_proj[src])²`` proxy and covers BOTH the relay and the combine. MC 0.02–0.27%.

    ⭐⭐ **``s`` CAN NO LONGER BE ``+∞``, and that is why there is no second case here.** It used to be,
    at any zero-count slot, because :func:`composition_logvar`'s counting term was ``1/n``; ``1/(1/p + ∞)``
    then annihilated all three streams of every message such a slot sent — including the MEASUREMENT
    stream, which never multiplies by ``r`` at all. That was one bug in two places, and it is fixed at the
    source: :func:`count_logvar`. ⛔ An ``~isfinite`` guard was briefly added HERE instead and is deleted —
    it treated the symptom, and it silently made every genuinely-unscaled hop free. TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero."""
    g = np.asarray(graft, bool)
    s = np.asarray(logvar_tot_dst, np.float64) + np.asarray(logvar_tot_src, np.float64)
    return np.where(g, 0.0, s)


def mismatch_gap(rho_msg, rho_own):
    """the-mismatch-deflation — the observed composition gap ``G`` between a message and the destination's own belief, and the
    CONTRADICTED mask. Both densities are per-component and in the destination's frame, already normalized to
    the destination's mass, so ``E_c`` and ``M_dst`` cancel and ``G = log(ρ^msg/ρ^own) = mo^msg − mo^own`` is a
    pure statement about the SPLIT.

    ``contradicted`` marks exactly one side being zero. A message with ``ρ^msg = 0`` asserts the component is
    ABSENT (``log f_c → −∞``) — an infinitely large gap against a destination whose own evidence says it is
    present, so its honest mismatch variance is ``+∞`` rather than whatever a numerical floor would make it.
    Returning it as a MASK is what keeps the zero-test ``_EPS`` a zero-test instead of promoting it to a model
    constant that sets the surviving precision. Both sides zero ⇒ ``G = 0``: nobody has an opinion.

    Returns ``(gap, contradicted)``."""
    t, o = _f(rho_msg), _f(rho_own)
    live_t, live_o = t > _EPS, o > _EPS
    both = live_t & live_o
    gap = np.where(both, np.log(np.maximum(t, _EPS) / np.maximum(o, _EPS)), 0.0)
    return gap, live_t ^ live_o


def mismatch_deflate(precision, gap, contradicted, var_own):
    """the-mismatch-deflation — deflate a message precision by the DerSimonian–Laird between-source mismatch variance::

        b̂² = max(0, G² − v_msg − v_own)          v_msg = 1/precision
        p_effective = 1 / (v_msg + b̂²) = 1 / max(v_msg, G² − v_own)

    The three regimes fall out of ``var_own`` alone, with no gate and no threshold:

    * ``var_own = +inf`` — the destination has NO composition evidence (τ_λ = 0: every AMBIG node, and every
      node on unstranded data) ⇒ ``b̂² = 0`` ⇒ the message passes at full precision. This is deliberate and it
      is what makes the term safe: where messages are the only information, there is no own opinion to
      contradict them, and inventing one would be the count-votes-composition regression.
    * ``var_own`` finite — a message that AGREES is barely touched; one that CONFLICTS is killed.
    * ``var_own = 0`` — composition CERTAIN (a structural pure-gDNA anchor) ⇒ the full ``G²`` is charged.

    The closed form ``p_eff = 1/max(v_msg, G²−v_own)`` states the safety property exactly: **a message can only
    out-weigh the destination's own belief when it agrees with it to within ``√2·σ_own``**, at every node and
    every composition — "pass-0 must be weak and correctable" as an inequality rather than a tuning knob.

    ``precision = 0`` (no message) stays 0; a CONTRADICTED claim goes to 0 wherever the destination has any
    evidence at all. Never returns a nan or an ``∞`` (the ``+inf`` var_own is substituted BEFORE the
    subtraction, so there is no ``inf − inf``). At ``b = 0`` the estimator is positively biased by
    ``E[max(0, χ²₁−1)]·(v_msg+v_own) ≈ 0.484·(v_msg+v_own)`` — the OVER-damping direction, and harmless because
    a message that agrees with the own belief moves the fused mode nowhere."""
    p = _f(precision)
    g = _f(gap)
    vo = _f(var_own)
    known = np.isfinite(vo)
    v_msg = np.where(p > 0.0, 1.0 / np.maximum(p, _EPS), 0.0)
    b2 = np.where(known, np.maximum(0.0, g * g - v_msg - np.where(known, vo, 0.0)), 0.0)
    out = np.where((p > 0.0) & (b2 > 0.0), 1.0 / np.maximum(v_msg + b2, _EPS), p)
    return np.where(np.asarray(contradicted, bool) & known, 0.0, out)


def graft_premise_logvar(flux_a, flux_b, var_a, var_b):
    """⭐ P1d — the GRAFT's PREMISE log-variance, measured from an exon's two flanking seams.

    **The event it prices.** The graft hands an exon the RNA at one seam — ``ρ_ν + ρ_μ`` — as *the* exon's
    RNA density. Every molecule counted there is in the exon, but the exon may also hold molecules that
    never touch that seam: ones that reach it by the other flank, or that start or end inside it. So what
    the graft knows is an **inequality**, ``ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)``, and it uses it as an equality.
    Nothing else in the ledger prices that. **the-graft-frame-variance comes closest and does not cover it**: the-graft-frame-variance charges
    ``(log r)²``, i.e. it assumes the only reason a seam and an exon differ is the capture step between
    them — so off-capture ``r = 1`` and the-graft-frame-variance charges exactly zero, while a seam and a region still differ.
    ``1/n_spl`` does not cover it either, and never will however large it grows: a count says *how many I
    counted*, this says *whether what I counted speaks for the exon*.

    **The estimator.** The exon's two flanking seams are two INDEPENDENT statements of the same claim —
    measured, not assumed: ``corr(log φ_left, log φ_right)`` is **0.017 raw, −0.036 after the Poisson part**
    at capture-OFF over 14 conditions (0.30 under capture, where the two do share a frame error). So the
    part of their squared gap that their own noise cannot explain IS the premise variance, halved because
    each seam carries half of the gap::

        v_premise = max(0, d² − v_a − v_b) / 2 ,      d = log(ρ_μ^a / ρ_μ^b)

    the same **method of moments** the calibrator already uses for ``κ`` and for both strand
    overdispersions — a fitted quantity, not a tuned constant, and there is no coefficient to choose: the
    ``/2`` is the measured independence and the truncation at 0 is the method's own ("no detectable
    premise error"). Returns ``(per_edge, pooled)``.

    ⚠⚠ **`ω_graft` IS A DEBT, NOT A MODEL.** It partially compensates for a FAILURE IN THE STRUCTURAL
    REPRESENTATION: the region/boundary map has no TSS/TES, so the solver cannot tell a splice junction from
    a transcript terminus — and that distinction is the whole of the effect this term prices (`ω̂` 1.7–1.9 at
    terminus boundaries vs 0.04–0.06 at junction-only ones, a ≥30× split, with 20.8 % of edges carrying
    71.7 % of the error). One library-wide average is standing in for a bimodal quantity. It works because
    over-charging a variance is cheap and under-charging is expensive, **not because it is right**, and it is
    expected to be FRAGILE ON REAL DATA (fitted on ~200 exons at 30–50 % SE, on a Poisson-by-construction
    suite whose isoforms are only nested truncations plus exon skips). **It MUST be re-derived as a per-class
    quantity the moment TSS/TES enter the region map (P1g)** — same equation, one scalar per structural class;
    the partial-pooling block below is the plug-in point.

    **⭐ The POOLED value is what is used — the per-edge one is returned for diagnostics only.** Both halves
    of that were wrong in the first landing (2026-07-26) and are corrected here:

    * **Statistically**, ``d²`` from ONE pair is a single draw of a χ²₁ scaled by the true variance, so its
      own coefficient of variation is √2. A per-edge "measurement" of a variance from one pair is mostly
      noise, and it both over- and UNDER-charges around the right mean — the under-charging half doing the
      damage, because it REPLACES the population value on the ~48 % of edges where a pair exists. Measured:
      per-edge+pooled 870,245 confidently-wrong reads at ``z2`` 11.12; pooled alone **762,000 at 8.98**, with
      exon-single ``z2`` 5.68 → **2.46**. A derived PARTIAL-POOLING variant (James–Stein, with the shrinkage
      weight ``B = τ²_b/(τ²_b + Var(ω̂_i))`` and ``Var(ω̂_i) = (2ω+v_i)²/2`` exact from the χ²₁) was also
      implemented and **deleted**: its own weight comes out **B = 0.82–0.89** — the between-node heterogeneity
      is real and dominant — and it STILL loses (765,281 / 9.23), because ``ω_i`` is right-skewed (the top
      decile of pairs carries 78–92 % of ``Σd²``) so any per-node point estimate charges the median node ≈0,
      and the loss is asymmetric: over-charging only widens a message, under-charging leaves a node
      confidently wrong. **The population mean is right precisely because the loss is asymmetric.**
      Do not rebuild the per-node form.
    * **Structurally** (owner, 2026-07-26), a per-edge form makes the message from the LEFT seam carry a
      variance computed from the RIGHT seam's counts, so a non-adjacent node's data reaches the destination
      twice — a real BP violation. With the pooled form no message's precision depends on anything but its
      own edge and one library constant, which is exactly the standing ``κ`` and both strand overdispersions
      already have.

    **Two things about the fit that were asked and are settled by measurement (2026-07-26):**

    * **A plain average, not an inverse-variance (DerSimonian–Laird) one.** DL is the efficient estimator of
      a *common* effect; there is no common effect here (`τ²_between` = 1.35–3.0 against a sampling variance
      of 0.26–0.78). What must be charged to a message is the EXPECTED squared premise error of a
      randomly-drawn edge, i.e. the arithmetic mean `E[ω_i]` — which is exactly what the plain moment
      difference estimates and DL does not. Measured, the two land within ~1 % of each other anyway
      (confidently-wrong 762,000 vs 755,727; DL slightly worse on error mass), so this is a correctness
      argument, not a scoring one.
    * **The value is fitted from the CURRENT belief, weakly.** The two fluxes are lifted by frames read off
      `ρ_tot(f_cur)`, so the fit is recomputed 3× per sweep — measured 0.2968 / 0.2968 / 0.2937, i.e. it
      moves <2 % as the belief updates. The count and `logvar_tot` legs of the noise subtraction are fixed
      data. The feedback loop this creates can only ever WIDEN a message, so it cannot manufacture
      confidence.

    The *population* the scalar is fitted on is still exons that have two live seams — and that population
    is not representative in the owner's sense: an exon whose flank is a transcript terminus or an
    exon↔exon boundary has no second seam and never enters the fit, and it is measurably the WORSE-behaved
    half (premise 0.60 with one seam, 2.69 with none, against 0.48 with two). So the fit UNDER-states:
    the count-zero-information-safe direction, and the "weak and correctable" one.

    **What the term does NOT claim.** It does not claim the second seam is a good proxy for the exon's RNA on
    any individual exon. It uses the seam pairs only to estimate ONE library-level number — the typical size
    of the graft's premise failure — because that failure is an UNDER-claim, and an under-claim has **no
    local signature**: a seam accounting for less mass than the exon holds is indistinguishable from an exon
    with more gDNA, which is count-zero-information exactly. The over-claim direction *is* locally visible
    (``claim/obs``) and belongs to P1e.

    **What the caller must pass.** Both fluxes ALREADY LIFTED into the destination's frame, so a capture
    step common to the two seams cancels out of ``d`` and only a genuine abundance difference is charged;
    and each ``var`` must carry every noise source the model knows — the seam's spliced COUNT (never its
    mass) **⊕ its lift's own scale sampling** (the-reframe-scale-variance's source leg; the destination's leg is common to both
    lifts and cancels in ``d``). Method of moments books as premise error every noise it fails to subtract,
    and omitting the lift term inflates the fit in proportion to gDNA depth (the frame is read off ``ρ_tot``,
    which gDNA makes noisier) — measurably the worst place to over-charge, because that is exactly where the
    RNA claim is a near-exact measurement and the gDNA claim is the wrong one.

    **It cannot double-count with the-mismatch-deflation** for the reason proves in general: DL's
    ``max()`` means a newly-modelled variance *replaces* the part of the gap the-mismatch-deflation was absorbing as unexplained
    drift, one for one, until ``b̂²`` hits its floor at 0.

    ⚠ The suite this was measured on is **Poisson by construction**, so the premise variance here is purely
    structural/annotation-driven. Real libraries add overdispersion on top of every term, which makes this
    an under-estimate rather than an over-estimate — but the SHAPE (terminus vs junction-only) must be
    re-measured on a real annotation before it is trusted quantitatively."""
    fa, fb = _f(flux_a), _f(flux_b)
    ok = (fa > _EPS) & (fb > _EPS)
    d = np.log(np.maximum(fa, _EPS)) - np.log(np.maximum(fb, _EPS))
    va, vb = _f(var_a), _f(var_b)
    v = np.where(np.isfinite(va), va, 0.0) + np.where(np.isfinite(vb), vb, 0.0)
    per = np.where(ok, np.maximum(0.0, d * d - v) * 0.5, 0.0)
    if not np.any(ok):
        return per, 0.0
    d2, vv = d[ok] ** 2, v[ok]
    return per, float(max(0.0, float(d2.mean()) - float(vv.mean())) * 0.5)
