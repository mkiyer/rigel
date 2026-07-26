"""Enrichment-frame primitives — the generalized enrichment-ratio solver's pure arithmetic (phase E1).

Hybrid capture multiplies *every* component of a node's density by the same efficiency ``e(x)``, so it is
visible in the **total density and only in the total** (a composition *ratio* ``k = ρ_g/ρ_R`` cancels ``e``).
A message carries a per-base density, so it is transportable across an enrichment cliff only after being
scaled into the destination node's frame by the *measured* ratio of total densities. This module owns that
arithmetic and nothing else: it is a pure function of masses, compositions, and effective lengths — no chain,
no substrate, no solver state — so its facts are pinned by closed-form unit tests
(``tests/calibration/test_enrichment_frame.py``) that cannot drift with solver behaviour.

Derivation: ``docs/calibration/enrichment_ratio_generalization.md`` (the framework, §1 the total-density
pivot, §2 the bounding lemma, §4 the step-wise junction solve) and ``junction_enrichment_scaling.md`` (the
junction instance + §5b the r₂/r₁ asymmetry). The solver wiring that consumes these is phase E2 (behind a
flag); this module is frame-agnostic by construction (arrays in, arrays out), so the per-face-vs-node-level
frame question E2 must settle cannot invalidate anything here.

The design principle behind the API shape: **transport ``k``, never ``f_g``.** ``k`` is enrichment-free and
component-set-shared; ``f_g`` depends on the node's own crossing-vs-contained effective lengths, so copying it
across a frame change is a known past defect (§2 of the junction doc). And the composition assumption is
carried as a **variance**, not a bool gate: the bounding lemma (§2) says a totally wrong composition still
pins ``ρ_tot`` to within the effective-length ratio (1.04–1.5× for normal nodes, 4×+ for short regions), which
is a continuous quantity — :func:`composition_logvar` derives it with no tuned threshold, so a short region
down-weights itself automatically. The only genuine bool is :func:`gdna_fallback_admissible`, the §5
*structural* refusals where the 100 %-gDNA default is qualitatively wrong (AMBIG, exon↔exon, no junction,
retained intron) and a wrong frame is worse than no message.
"""

from __future__ import annotations

import numpy as np
from scipy.special import log_ndtr, polygamma

__all__ = [
    "total_density",
    "k_from_belief",
    "f_g_from_k",
    "boundary_unspliced_from_k",
    "enrichment_ratio",
    "composition_logvar",
    "gdna_fallback_admissible",
    "reframe_density",
    "density_mode_logfrac",
    # ── the pass-0 message-VARIANCE laws (message_variance_derivation.md) ──
    "transport_seed_logvar",
    "graft_rna_logvar",
    "graft_frame_logvar",
    "peel_rna_logvar",
    "peel_continue_share",
    "peel_share_logvar",
    "residual_level",
    "conservation_rescale",
    "transfer_logvar",
    "message_precision",
    "mismatch_gap",
    "mismatch_deflate",
]

_EPS = 1.0e-12

# Numerical-resolution knobs for :func:`conservation_rescale`'s root solve, NOT model constants: the bracket
# doubles from ±1 so 48 expansions reach |mu| ~ 1e14, and 96 bisections then resolve it to ~1e-14 of the
# bracket. Both are "enough to be exact in float64", not tunables — the answer is the unique root either way.
_RESCALE_BRACKET = 48
_RESCALE_BISECT = 96
_EXP_LIM = 700.0  # the float64 overflow guard on exp()


def _f(x):
    return np.asarray(x, dtype=np.float64)


def total_density(mass, f_g, E_g, E_r):
    """``ρ_tot = M·[ f_g/E_g + (1−f_g)/E_r ]`` — nucleic acid per base, the quantity an enrichment ratio is
    formed from (``enrichment_ratio_generalization.md`` §1).

    The bracket is a **convex combination** of ``1/E_g`` and ``1/E_r`` (the composition-weighted harmonic
    blend of the two effective lengths), so ``ρ_tot`` is *not* ``M/E_anything`` — that is exactly why it
    depends on composition, and why :func:`composition_logvar` exists. It degenerates to ``M/E_g`` at
    ``f_g = 1`` (an intergenic region / seam, structurally pure gDNA) and to ``M/E_r`` at ``f_g = 0``.
    """
    m, fg, eg, er = _f(mass), _f(f_g), _f(E_g), _f(E_r)
    return m * (fg / np.maximum(eg, _EPS) + (1.0 - fg) / np.maximum(er, _EPS))


def k_from_belief(f_g, E_g, E_r):
    """``k = ρ_g/ρ_R`` from a node's composition and its own effective lengths — the **enrichment-free**
    content ratio (``junction_enrichment_scaling.md`` §2).

    ``ρ_c = f_c·M/E_c`` ⇒ ``k = (f_g/E_g)/((1−f_g)/E_r) = f_g·E_r / ((1−f_g)·E_g)``. The mass ``M`` cancels,
    and so does the capture efficiency ``e`` (component-agnostic at pass-0), which is what makes ``k``
    transportable between component-set-matched nodes where ``f_g`` is not. Returns ``+inf`` at ``f_g = 1``
    (pure gDNA)."""
    fg, eg, er = _f(f_g), _f(E_g), _f(E_r)
    num = fg * er
    den = (1.0 - fg) * eg
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > _EPS, num / np.where(den > _EPS, den, 1.0), np.inf)


def f_g_from_k(k, E_g, E_r):
    """Re-form ``f_g`` in the destination's frame from a transported ``k`` (the inverse of
    :func:`k_from_belief`): ``f_g = k·E_g / (k·E_g + E_r)`` (``junction_enrichment_scaling.md`` §2/§3.1).

    This is the correct way to move a composition across a frame change — transport ``k``, then re-form
    ``f_g`` from the destination's own ``E_g, E_r``. Copying ``f_g`` directly is the known defect."""
    kk, eg, er = _f(k), _f(E_g), _f(E_r)
    finite = np.isfinite(kk)
    # substitute a finite k BEFORE the divide — k = +inf would form inf/inf = nan even though it is masked
    # out (trap #6). The pure-gDNA limit f_g → 1 is supplied explicitly.
    kk_fin = np.where(finite, kk, 0.0)
    keg = kk_fin * eg
    den = keg + er
    out = np.where(den > _EPS, keg / np.maximum(den, _EPS), 1.0)
    return np.where(finite, out, 1.0)


def boundary_unspliced_from_k(k, mass, E_g, E_r):
    """The step-wise resolve of a junction boundary's unspliced composition from a transferred ``k``
    (``enrichment_ratio_generalization.md`` §4 step 1; ``junction_enrichment_scaling.md`` §3.1).

    A junction boundary's unspliced crossing shares its component set ``{g, R_cont}`` with the flanking
    intron, so ``k(B) = k(I)`` transfers, and the boundary's own mass identity ``M_B = ρ_g·E_g + ρ_R·E_r``
    closes the system in the boundary's **own frame** — no enrichment enters (``M_B`` is observed, ``k`` is
    enrichment-free)::

        ρ_R = M / ( k·E_g + E_r ) ,   ρ_g = k·ρ_R

    Returns ``(ρ_g, ρ_R)``. At ``k = +inf`` (the 100 %-gDNA fallback, §4 step 2) this yields
    ``ρ_R = 0, ρ_g = M/E_g`` — the pure-gDNA density in the boundary's frame, resolving the *frame* without
    asserting a belief (the caller supplies zero precision)."""
    kk, m, eg, er = _f(k), _f(mass), _f(E_g), _f(E_r)
    # finite-k branch: rho_R = M/(k*E_g+E_r); infinite-k (pure gDNA) branch: rho_R = 0, rho_g = M/E_g
    finite = np.isfinite(kk)
    # substitute a finite k BEFORE any product — np.where evaluates both branches, so leaving inf in would
    # form inf·0 = nan in the pure-gDNA branch (trap #6). The infinite-k value is supplied explicitly below.
    kk_fin = np.where(finite, kk, 0.0)
    den = kk_fin * eg + er
    rho_R = np.where(finite, m / np.maximum(den, _EPS), 0.0)
    rho_g = np.where(finite, kk_fin * rho_R, m / np.maximum(eg, _EPS))
    return rho_g, rho_R


def enrichment_ratio(rho_tot_dst, rho_tot_src, var_log_dst=0.0, var_log_src=0.0):
    """The enrichment ratio that scales a source density **into the destination's frame** — a single
    orientation, ``r = ρ_tot(dst)/ρ_tot(src)``, so a transported density is ``ρ_dst = ρ_src·r``
    (``enrichment_ratio_generalization.md`` §3.4; the docs write ``r₁ = ρ(E)/ρ(B)`` and ``r₂ = ρ(B)/ρ(I)`` in
    opposite orientations — this canonicalizes to one direction to make orientation bugs impossible).

    Because ``r`` compares two total densities that share their component set, ``e(dst)/e(src)`` is isolated —
    that is the whole content of ``r`` (§3.3). The two totals are ratios of *solved* or partially-solved
    quantities, so ``r`` carries their uncertainty; the log-variance is additive for a ratio of independent
    quantities::

        Var(log r) = Var(log ρ_tot(dst)) + Var(log ρ_tot(src))

    Each input variance comes from :func:`composition_logvar`. Carrying it is what stops E2 from repeating the
    §13.6h failure — an uncertain scale factor applied as exact. Returns ``(r, var_log_r)``."""
    rd, rs = _f(rho_tot_dst), _f(rho_tot_src)
    r = np.where((rd > _EPS) & (rs > _EPS), rd / np.maximum(rs, _EPS), np.nan)
    var_log_r = _f(var_log_dst) + _f(var_log_src)
    return r, var_log_r


def composition_logvar(f_g, E_g, E_r, var_fg, n):
    """``Var(log ρ_tot)`` — the composition assumption carried as a **variance**, not a bool gate.

    Two independent sources of uncertainty in a total density (``enrichment_ratio_generalization.md`` §2, the
    bounding lemma, made continuous):

    * **counting** — ``M`` is a Poisson count, so ``Var(log M) = 1/n``;
    * **composition** — ``ρ_tot = M·B(f_g)`` with ``B = f_g/E_g + (1−f_g)/E_r``, and ``dB/df_g = 1/E_g − 1/E_r``,
      so by the delta method ``Var(log ρ_tot | composition) = [(1/E_g − 1/E_r)/B]²·Var(f_g)``.

    ::

        Var(log ρ_tot) = 1/n  +  [ (1/E_g − 1/E_r) / B ]² · Var(f_g)

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
    count = np.where(nn > 0.0, 1.0 / np.maximum(nn, _EPS), np.inf)
    return count + comp


def reframe_density(rho_c_src, rho_tot_dst, rho_tot_src):
    """Scale a source component density into the destination's enrichment frame — the ONE operation of the
    unified solver (`unified_solver_design.md` §2). ``ρ_c(dst) = ρ_c(src) · ρ_tot(dst)/ρ_tot(src)``.

    The ratio ``r = ρ_tot(dst)/ρ_tot(src) = e(dst)/e(src)`` (an enrichment ratio, robust to the unknown
    composition by the bounding lemma), so the reframed density is at the destination's enrichment. Applied to
    **every** component of a message by the SAME ``r``. It cancels when the message is normalized by the
    destination's own mass (:func:`density_mode_logfrac`), so within a single matched message the reframe is
    invisible; across two messages from different-enrichment sources it is load-bearing (the commensuration)."""
    rc, rd, rs = _f(rho_c_src), _f(rho_tot_dst), _f(rho_tot_src)
    r = np.where(rs > _EPS, rd / np.maximum(rs, _EPS), 1.0)
    return rc * r


def density_mode_logfrac(rho_c_dst_frame, E_c_dst, mass_dst):
    """The message factor on ``log f_c`` — the density mode, ÷ the destination's OWN observed mass
    (`unified_solver_design.md` §2): ``log( ρ_c(dst-frame) · E_c^dst / M_dst )``.

    Enrichment-correct **once ``ρ_c`` is reframed** (:func:`reframe_density`): the reframe carries
    ``e(src)→e(dst)`` and dividing by ``M_dst`` (which carries ``e(dst)``) cancels it, leaving the
    enrichment-free fraction. ÷``M_dst`` (not the imputed total) is what makes it universal — a partial source
    (gDNA only) gives ``f_g = ρ_g·E_g/M_dst < 1`` correctly, where the shift's imputed-total normalizer would
    assert ``f_g = 1``. This is the single factor that retires the shift/density mode split. Returns the log
    (a per-component ``imp_mode`` for the ψ solve); the domain guard floors the argument at ``_EPS``."""
    rc, ec, md = _f(rho_c_dst_frame), _f(E_c_dst), _f(mass_dst)
    return np.log(np.maximum(rc * ec / np.maximum(md, _EPS), _EPS))


def gdna_fallback_admissible(is_ambig, mature_crosses, has_junction, retained_intron):
    """The one genuine bool: is the 100 %-gDNA-at-zero-precision fallback structurally *admissible* on this
    boundary's unspliced crossing (``enrichment_ratio_generalization.md`` §5; ``junction_enrichment_scaling.md``
    §5 K7)?

    The fallback assumes the unspliced crossing is ``{g, R_cont}`` with negligible contiguous RNA. §2's
    bounding lemma does **not** cover a component set that is *qualitatively* different, so where the crossing
    carries something else the fallback must be **refused** (the caller emits nothing rather than a wrong
    frame). The four structural refusals — every one a constant-free boolean already available upstream:

    * ``is_ambig`` — opposite-strand overlap: the crossing carries another transcript's RNA (``free_pos & free_neg``);
    * ``mature_crosses`` — exon↔exon: contiguous mature crosses, the unspliced pool is RNA-rich (``mrna_active_s``);
    * ``~has_junction`` — no splice junction: no routing split, RNA simply continues (``S_B = 0``);
    * ``retained_intron`` — the contiguous RNA is not scarce (the intron's own solve).

    Admissible ⟺ a genuine junction with none of the RNA-rich exceptions. Returns bool, same shape."""
    ambig = np.asarray(is_ambig, dtype=bool)
    mature = np.asarray(mature_crosses, dtype=bool)
    junction = np.asarray(has_junction, dtype=bool)
    retained = np.asarray(retained_intron, dtype=bool)
    return junction & ~ambig & ~mature & ~retained


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# The pass-0 message-VARIANCE laws (``docs/calibration/message_variance_derivation.md``).
#
# Built on the ``τ_λ`` FOUNDATION (``variance_foundation_proposal.md``): a source component's message log-variance
# is its composition variance ``Var(log f_c)`` (the Schur-marginal Jacobian) ⊕ its Poisson sampling ``1/n``,
# transported across an edge by a SUM (the GRAFT, share-weighted, convex ≤1) or a DIFFERENCE (the PEEL,
# u-weighted, ≥1), damped by the enrichment-transfer variance ``σ²_transfer = Var(log r)`` — DIRECTION-dependent
# (~0 on the graft where the reframe cancels, load-bearing on the peel/anchor). The ÷M_dst normalizer is
# common-mode and cancels, so the ψ message precision is simply ``1/Var(log ρ_c)`` — NO destination Jacobian,
# NO ``1/n_dst``. Every law MC-validated (``scratchpad/message_variance_mc.py``; independently re-derived +
# adversarially verified, workflow ``wf_c952640d``) to <1% in-regime. Variances live in log space (additive);
# convert to a ψ precision only at the handoff (:func:`message_precision`).
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────


def transport_seed_logvar(var_log_f, n):
    """M1 — the per-component SOURCE log-variance a belief carries when it becomes a message: composition ⊕
    sampling, ``Var(log ρ_c) = Var(log f_c) + 1/n``.

    Count-zero-information, made precise: the count sets ONLY the total precision ``1/n``; the gDNA-vs-RNA split
    precision is ``Var(log f_c)`` (from ``τ_λ``), NEVER the count. The MEASURED spliced is the
    ``var_log_f = 0`` (composition-CERTAIN) case ⇒ pure ``1/n_spl``. ``n = 0`` (no count) ⇒ ``+inf`` (no
    message). This is the *transport* seed — distinct from the *fusion* weight ``1/Var(log f_c)`` (composition
    only, no ``1/n``); the ``1/n`` enters HERE, not in the fusion (handoff §7, the composition⟂sampling split).

    ⚠ ``Var(log f_c)`` is the FOUNDATION Jacobian; near the ``f_c→0`` wall a small, high-CV minority component's
    log-variance is UNDER-stated (over-confident 36–92% — ``message_variance_derivation.md`` §4 guard). Not
    corrected here (a foundation caveat); watch the minority-arm message in the A/B."""
    v = np.asarray(var_log_f, np.float64)
    nn = np.asarray(n, np.float64)
    count = np.where(nn > 0.0, 1.0 / np.maximum(nn, _EPS), np.inf)
    return v + count


def graft_rna_logvar(v_nu, v_mu, w_nu, w_mu):
    """M2 — the GRAFT (boundary→exon) RNA message log-variance. The exon receives the boundary's WHOLE RNA
    flux, so the RNA density is a **sum** ``ρ_R = ρ_ν + ρ_μ`` (imputed-continue + measured-splice). The delta
    method on a sum gives the **share-weighted** (convex, ≤1) combination — the item-E rule, a minority
    component contributes quadratically LITTLE::

        Var(log ρ_R) = w_ν²·v_ν + w_μ²·v_μ ,     w_ν = ρ_ν/ρ_R,  w_μ = ρ_μ/ρ_R,  w_ν+w_μ = 1

    ``v_ν`` is the continue seed (``Var(log f_R)+1/n``), ``v_μ = 1/n_spl`` (measured). Enrichment CANCELS on the
    matched-set graft, so ``σ²_transfer`` is added as 0 there (:func:`transfer_logvar`). The gDNA grafts alone
    (its message log-variance is just its transport seed). MC <1%."""
    wn, wm = np.asarray(w_nu, np.float64), np.asarray(w_mu, np.float64)
    return wn * wn * np.asarray(v_nu, np.float64) + wm * wm * np.asarray(v_mu, np.float64)


def graft_frame_logvar(r):
    """M8 — the GRAFT's **frame-mislift** log-variance: ``(log r)²`` on the MEASURED spliced component.

    M5 sets ``σ²_transfer = 0`` on a graft edge because ``r`` is common-mode across the matched set ``{g, R}``
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

    So the grafted component has **no matched gDNA partner to cancel ``r`` against** — exactly M5's peel /
    partial-anchor case, where ``σ²_transfer`` is load-bearing. The un-modelled frame step is what the message
    gets wrong, its magnitude is the step itself, and the model's own estimate of that step is ``r``; charging
    ``(log r)²`` is the method-of-moments second moment of a single observation of it — the same logic as M7's
    ``b̂²``, and with no tuned constant. It is identically 0 when there is no frame change (``r = 1``), which is
    why the shipped model is EXACT off-capture (measured required correction ``log c`` = +0.009/−0.008/+0.054)
    and only bites where capture opens a step.

    ⚠ This is NOT the retired ``σ²_cliff = (log r)²`` proxy, which charged the whole cliff on **matched**
    reframes — where the composition genuinely IS preserved across the enrichment step, so it over-damped.
    Here the component set is *unmatched by construction*, so the cliff really is un-cancelled error.

    Applied to the spliced measurement's OWN precision (``1/n_spl ⊕ (log r)²``), so M2's share weighting
    ``w_μ²`` arises implicitly from the inverse-variance fusion with the correctly-framed ``ρ_ν`` arm — a
    graft that is a minority of the RNA is damped proportionately less.

    Calibration (delivered ``λ`` vs oracle, graft edges into exons): ``z2 = E[Δ²]/E[v]`` 58–310 → **2.1–3.8**,
    consistently across capture off and on."""
    rr = _f(r)
    lr = np.log(np.where(rr > _EPS, rr, 1.0))
    return lr * lr


def peel_rna_logvar(v_log_rho_R, s2_transfer, v_mu, u):
    """M3 — the PEEL (exon→boundary) RNA-continue message log-variance. The boundary receives only what
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
    """M10 — the fraction of a seam's RNA that CONTINUES unspliced: ``w = ρ_ν/(ρ_ν + ρ_μ)``.

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
    ``u``'s p75 on real junctions is ≈ 3. MC: `message_variance_mc.py` M10b, exact to 1e-12.

    Degenerate limits, both structural: no spliced flux (``ρ_μ = 0``) ⇒ ``w = 1``, nothing splices away and
    nothing is peeled; no RNA at the seam at all (``ρ_ν + ρ_μ = 0``) ⇒ ``w = 1``, there is nothing to
    apportion and the caller's own gates decide. ``w`` is always in ``[0, 1]`` — the peel can no longer go
    negative, which retires the zero-truncation defect (a fully-consumed peel used to emit ``ρ_ν = 0``, "no
    RNA continues past here", at a live precision)."""
    nu, mu = _f(rho_nu), _f(rho_mu)
    tot = nu + mu
    return np.where(tot > _EPS, nu / np.maximum(tot, _EPS), 1.0)


def peel_share_logvar(w_mu, v_nu, v_mu):
    """M10 — the log-variance the continuing SHARE contributes: ``w_μ²·(v_ν + v_μ)``.

    Delta method on ``log w = log ρ_ν − log(ρ_ν + ρ_μ)``: both partials are ``±w_μ``, the SPLICED share
    ``w_μ = 1 − w``, so::

        Var(log w) = w_μ²·( Var(log ρ_ν^B) + Var(log ρ_μ^B) ) ,   v_μ = 1/n_spl (measured, count-only)

    and the delivered message variance is ``Var(log ρ_R(x)) + σ²_transfer + Var(log w)``. The weights are
    **CONVEX** (``w_μ ≤ 1``) — the exact mirror of M2's graft SUM — where M3's difference carried ``u ≥ 1`` and
    amplified. So this move takes the peel out of the ill-conditioned regime entirely; there is no ``ε``/``u``
    validity limit to respect and no over-confidence tail. MC: `message_variance_mc.py` M10, rel 0.2–1.0 %."""
    wm = _f(w_mu)
    return wm * wm * (_f(v_nu) + _f(v_mu))


def residual_level(mass, n_mass, rho_g, E_g, E_r, v_g):
    """M11 — a node's own RNA density read off its **observed mass** against an imputed gDNA **density**.

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
    pdf_a = np.exp(-0.5 * alpha * alpha - 0.5 * np.log(2.0 * np.pi))
    pdf_b = np.exp(-0.5 * beta * beta - 0.5 * np.log(2.0 * np.pi))
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
    k = np.maximum(f_R * f_R / np.maximum(v_f, _EPS), 1.0)  # the RNA share's effective fragment count
    v_log = np.where(
        ok,
        polygamma(1, np.minimum(k, 1.0 / _EPS))
        + inv_n / np.maximum(f_R * f_R, _EPS),  # the count, with its exact 1/f_R Jacobian
        np.inf,
    )
    rho = np.where(ok, f_R * m_ / np.maximum(er_, _EPS), 0.0)
    # mask the ∞ BEFORE the product (np.where evaluates both branches: ``0·inf = nan``, the standing trap).
    return rho, v_log, np.where(ok, rho * rho * np.where(ok, v_log, 0.0), np.inf)


def conservation_rescale(mass, rho, eff, var, var_common, own):
    """M12 — restore ``Σ_c ρ_c·E_c = M`` by moving each component in proportion to how badly it is known.

    A claim about a node asserts a fragment count, ``S = Σ_c ρ_c·E_c``, and the node observed ``M``. That the
    two agree is an IDENTITY, not an approximation. When they do not, the claim must be projected onto the
    constraint — and the only question is **which component gives ground**. Today's operator (`_pin_v`) scales
    every component by the same ``k = M/S``, which punishes a component that is right in order to accommodate
    one that is wrong; measured on the failing arms, that is a direct spliced-RNA measurement (claimed/true
    1.0×) being shrunk 43 % to make room for an imputed gDNA level that is 47× too large.

    **The derivation.** Write the claim's error in log space, ``ρ_c = ρ_c^true·exp(ε_c)``, and take the
    correction ``a_c`` (``ρ_c → ρ_c·exp(a_c)``) that is most likely under the error covariance ``Σ``::

        minimise  aᵀΣ⁻¹a   subject to   Σ_c m_c·exp(a_c) = M ,    m_c = ρ_c·E_c

    Linearising the constraint (``α_c = m_c/S``, ``δ = log(M/S)`` ⇒ ``αᵀa = δ``) and solving the Lagrangian
    gives ``a = δ·Σα/(αᵀΣα)``. **Only the DIRECTION ``s = Σα`` comes from the error model.** The magnitude is
    then re-solved EXACTLY along that direction — ``find μ with Σ_c m_c·exp(μ·s_c) = M`` — because the
    constraint is an identity and deserves to hold exactly rather than to first order. ``g(μ)`` is
    non-decreasing for ``s ≥ 0``, so the root is unique and a bracketed bisection reaches it with no step
    control, no cap and no tuned constant.

    **The error model.** Every component of a message shares the reframe and the source's own count, so those
    are common-mode; whatever is left in a component's variance is its own::

        Σ = σ_cm²·11ᵀ + diag(w) ,   w_c = max(0, v_c − σ_cm²)   ⇒   s_c = σ_cm² + α_c·w_c

    **Two exact limits, and they are the justification:**

    * ``w → 0`` — the only error is the shared frame ⇒ ``s_c`` is constant ⇒ every component moves by the same
      factor. **That is exactly ``_pin_v``**: the shipped operator is the zero-independent-variance limit of
      this one, not a different operator.
    * ``σ_cm² → 0`` and one component's ``w_c → 0`` — that component is a MEASUREMENT, not an imputation ⇒
      ``s_c = 0`` ⇒ it does not move at all and the others absorb the whole residual.

    Capture-OFF drives the second limit (``r ≈ 1`` ⇒ ``σ_cm² ≈ 0``) exactly where the grafted spliced count
    makes the RNA arm a measurement; capture-ON drives the first (``Var(log r)`` large). The regime switch is
    not a rule — it is what the variances already say.

    **The rank-1 composition term is deliberately NOT in Σ, and that is a result.** A source's own ``f_g``
    error moves its components in opposite directions (``u = (f_R, −f_g)``, the shipped
    `node_init.own_composition_logvar`), and ``uᵀα = ∂log S/∂λ`` vanishes whenever the two nodes' effective-
    length RATIOS agree — so under a strict rank-1 model a conservation violation says nothing about
    composition and this operator degenerates to the common rescale. The violations actually observed
    (1.24–1.31×) far exceed the eff-length bound on composition drift (×1.04 contained, ×1.50 at a boundary
    crossing), so the error is **component-specific**, which rank-1 forbids and ``diag(w)`` admits. It can be
    component-specific because the RNA arm receives an independent additive measurement — the graft — that the
    gDNA arm does not.

    **"Supplied" is a statement about PRECISION** (``var`` finite), never about the density's value — the same
    test the λ-emission gate makes. A component the claim does not supply contributes the node's OWN density
    to the mass budget and **does not move**, which is what keeps a partial claim partial (a message carrying
    gDNA only still delivers ``f_g < 1``); rescaling all components blindly instead regresses capture-OFF 3.6×
    (`scratchpad/derive_2_relay_pin.py`).

    Arguments are ``mass`` ``(n,)``, ``var_common`` ``(n,)``, and ``rho`` / ``eff`` / ``var`` / ``own``
    ``(n, C)`` over the components. Returns the corrected ``rho`` ``(n, C)``. Unchanged where the claim
    already balances, where there is no mass, or where nothing may move."""
    M = _f(mass)
    r, E, v, o = _f(rho), _f(eff), _f(var), _f(own)
    vc = np.maximum(_f(var_common), 0.0)[..., None]
    supplied = np.isfinite(v)
    m = np.where(supplied, r, o) * E  # `_pin_v`'s partial-claim budget substitution
    S = m.sum(axis=-1)
    ok = (S > _EPS) & (M > _EPS)
    alpha = m / np.maximum(S, _EPS)[..., None]
    # the component's OWN variance: what is left once the shared part is removed. Mask the +inf BEFORE the
    # subtraction — an unsupplied component is zeroed by `supplied` anyway, but `inf - finite` would leak.
    w = np.maximum(np.where(supplied, v, 0.0) - vc, 0.0)
    s = np.where(supplied, vc + alpha * w, 0.0)
    # nothing may move (no supplied component, or a model with no uncertainty anywhere): fall back to the
    # common factor over whatever IS supplied — the σ_cm² → 0, w → 0 corner, which is `_pin_v`.
    dead = ~(s > _EPS).any(axis=-1)
    s = np.where(dead[..., None] & supplied, 1.0, s)
    tgt = np.where(ok, M, S)

    def _tot(mu):
        return np.sum(m * np.exp(np.clip(mu[..., None] * s, -_EXP_LIM, _EXP_LIM)), axis=-1)

    lo = np.full(S.shape, -1.0)
    hi = np.full(S.shape, 1.0)
    for _ in range(_RESCALE_BRACKET):  # geometric expansion until the root is straddled
        need_lo, need_hi = _tot(lo) > tgt, _tot(hi) < tgt
        if not (need_lo.any() or need_hi.any()):
            break
        lo = np.where(need_lo, lo * 2.0, lo)
        hi = np.where(need_hi, hi * 2.0, hi)
    for _ in range(_RESCALE_BISECT):  # bisection: g is non-decreasing in mu, so this is exact
        mid = 0.5 * (lo + hi)
        up = _tot(mid) < tgt
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    mu = np.where(ok, 0.5 * (lo + hi), 0.0)
    return r * np.exp(np.clip(mu[..., None] * s, -_EXP_LIM, _EXP_LIM))


def transfer_logvar(logvar_tot_dst, logvar_tot_src, graft):
    """M5 — ``σ²_transfer = Var(log r)``, the enrichment-ratio uncertainty that damps a message across a capture
    cliff: ``Var(log r) = Var(log ρ_tot^dst) + Var(log ρ_tot^src)`` (each from :func:`composition_logvar`).

    DIRECTION-dependent: on the **graft** the reframe ``r`` is common-mode across the matched component set and
    CANCELS in the composition (return 0 — applying it there is the double-count the density-uniformity proxy
    committed); on the **peel** / partial-anchor message it is LOAD-BEARING. This one law replaces the retired
    ``var_proj[dst] + (μ_proj[dst]−μ_proj[src])²`` proxy and covers BOTH the relay and the combine. MC 0.02–0.27%."""
    g = np.asarray(graft, bool)
    s = np.asarray(logvar_tot_dst, np.float64) + np.asarray(logvar_tot_src, np.float64)
    return np.where(g, 0.0, s)


def message_precision(var_log_rho_msg):
    """M4 — the ÷M_dst conversion: the ψ message precision from a component's message log-variance. Because
    ``ρ_tot^dst = M_dst·B_dst`` cancels ``M_dst`` in the reframe ``r``, the destination mass is common-mode and
    drops from the delivered composition mode AND its variance ⇒ the precision is simply ``1/Var(log ρ_c^msg)`` —
    **NO destination Jacobian, NO ``1/n_dst``** (verified to machine precision, invariant under a 100× M_dst
    spread). Returns 0 for a ``+inf`` / non-positive variance (no message), never a nan or ``∞``."""
    v = np.asarray(var_log_rho_msg, np.float64)
    ok = np.isfinite(v) & (v > _EPS)
    return np.where(ok, 1.0 / np.maximum(v, _EPS), 0.0)


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M7 — the CROSS-CLIFF COMPOSITION-MISMATCH variance (`SESSION_2026_07_24_HANDOFF_5.md` §4).
#
# M1–M5 price a message's SAMPLING error. They do not price the imputation PREMISE — "my neighbour and I share a
# composition" — being false, and that premise is the message's whole content. The delivered mode error splits
# EXACTLY (MC-validated to machine precision, `scripts/debug/message_variance_mc.py` M7a/M7b) into
#
#     mo_c − log f_c^dst,true  =  log(s_c^src/s_c^dst,true)  +  log(r̂/r_true)   ⇒   σ² = v_msg + b_c²
#
# — a composition-SHARE mismatch plus the reframe's own scale noise (M5's σ²_transfer, already inside v_msg).
# The bias b_c is a POPULATION quantity (the third information source), unavailable prior-free — but the
# destination holds an INDEPENDENT estimate of the same composition: its own message-free self-solve. Two
# estimators of one quantity is a two-study random-effects meta-analysis, and its between-study variance has a
# closed-form method-of-moments estimator (DerSimonian & Laird 1986) that introduces NO tuned constant.
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────


def mismatch_gap(rho_msg, rho_own):
    """M7 — the observed composition gap ``G`` between a message and the destination's own belief, and the
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
    """M7 — deflate a message precision by the DerSimonian–Laird between-source mismatch variance::

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
    Nothing else in the ledger prices that. **M8 comes closest and does not cover it**: M8 charges
    ``(log r)²``, i.e. it assumes the only reason a seam and an exon differ is the capture step between
    them — so off-capture ``r = 1`` and M8 charges exactly zero, while a seam and a region still differ.
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

    **⭐ The POOLED value is what is used — the per-edge one is returned for diagnostics only.** Both halves
    of that were wrong in the first landing (2026-07-26) and are corrected here:

    * **Statistically**, ``d²`` from ONE pair is a single draw of a χ²₁ scaled by the true variance, so its
      own coefficient of variation is √2. A per-edge "measurement" of a variance from one pair is mostly
      noise, and it both over- and UNDER-charges around the right mean — the under-charging half doing the
      damage, because it REPLACES the population value on the ~48 % of edges where a pair exists. Measured:
      per-edge+pooled 870,245 confidently-wrong reads at ``z2`` 11.12; pooled alone **762,000 at 8.98**, with
      exon-single ``z2`` 5.68 → **2.46**. Shrinking all the way to the pooled mean is simply the better
      estimator, and this is the ordinary empirical-Bayes result, not a tuning choice.
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
    mass) **⊕ its lift's own scale sampling** (M5's source leg; the destination's leg is common to both
    lifts and cancels in ``d``). Method of moments books as premise error every noise it fails to subtract,
    and omitting the lift term inflates the fit in proportion to gDNA depth (the frame is read off ``ρ_tot``,
    which gDNA makes noisier) — measurably the worst place to over-charge, because that is exactly where the
    RNA claim is a near-exact measurement and the gDNA claim is the wrong one.

    **It cannot double-count with M7** for the reason `variance_ledger.md` §2.2 proves in general: DL's
    ``max()`` means a newly-modelled variance *replaces* the part of the gap M7 was absorbing as unexplained
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
        return per, 0.0, per
    d2, vv = d[ok] ** 2, v[ok]
    pooled = float(max(0.0, float(d2.mean()) - float(vv.mean())) * 0.5)
    # ── PARTIAL POOLING (the shrunk per-node estimate; the caller chooses whether to use it) ─────────────
    # ``ω_i`` from one pair is a χ²₁ draw, so its sampling variance is known EXACTLY and needs no constant:
    # ``d_i² ~ (2ω + v_i)·χ²₁`` ⇒ ``Var(d_i²) = 2(2ω+v_i)²`` ⇒ ``Var(ω̂_i) = (2ω+v_i)²/2``. The BETWEEN-node
    # variance is what the population's spread has left over after that, and the James–Stein / hierarchical
    # weight follows: ``B_i = τ²_b / (τ²_b + Var(ω̂_i))``. Every piece is estimated from the same data.
    # ``B → 0`` (full pooling) exactly when the population is explained by sampling noise alone; ``B → 1``
    # (per-node) when the between-node spread dominates. No knob decides which — the data does.
    v_samp = (2.0 * pooled + vv) ** 2 * 0.5
    tau_b = max(0.0, float(per[ok].var()) - float(v_samp.mean()))
    b = tau_b / np.maximum(tau_b + v_samp, _EPS)
    shrunk = np.full(per.shape, pooled)
    shrunk[ok] = np.maximum(0.0, pooled + b * (per[ok] - pooled))
    return per, pooled, shrunk
