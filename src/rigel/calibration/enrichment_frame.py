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
    "peel_rna_logvar",
    "transfer_logvar",
    "message_precision",
]

_EPS = 1.0e-12


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
