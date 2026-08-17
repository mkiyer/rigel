"""ψ's REFERENCE MEAN, SET FROM THE ANNOTATION — the claim, its strength, and the asymmetry it exposed.

`~rigel.calibration.simplex_logodds.structural_reference_location` sets ψ's per-slot reference mean where
**no annotated MATURE transcript is continuous across the position**, and leaves the shipped neutral ½
everywhere else. ``CalibrationConfig.structural_reference`` is its switch.

**THE TWO HALVES ARE CHOSEN BY DIFFERENT ARGUMENTS AND THE FILE IS ORGANISED AROUND THAT.**

⭐ **The LOCATION is annotation-derived**: *nascent RNA is zero unless proven otherwise* (owner) — so
presume gDNA at the four structural classes. The 16-condition ladder ranks this correctly, and it is the
design target.

⭐⭐ **The STRENGTH is ONE PSEUDO-OBSERVATION, and the ladder CANNOT rank it.** ``m = (a+1)/(a+b+1) =
0.75`` exactly — the mean one gDNA observation takes ``Beta(a,b)`` to. Since ``strength = logit(m)``, that
is ``log 3`` nats, i.e. a **3:1** claim, refutable at ~1.5 fragments. ⛔ It replaces ``m = σ(L)``, the
lattice's own width, which asserted **9.31 nats (10,000:1)** and measured WORSE THAN NO PRIOR AT ALL at
being refuted. ⚠ The ladder holds ``nrna = 0``, so it scores only the DELIVER obligation, where more nats
is monotonically better — which is exactly why it preferred the worst option
(`TRAPS: a-single-level-panel-cannot-see-a-constant`, met on a strength rather than a level).

⭐⭐⭐ **AND A REFUTABILITY TEST IS ONLY VALID IF THE REFUTATION CHANNEL IS IN THE FIXTURE.** The two
mechanisms that can prove nascent RNA present are ① strand asymmetry inside an intron (`i_strand`, dead at
κ=½) and ② intron-vs-intergenic density (`tau_fac`, alive unstranded, ``intron_factory = True`` in
production). ⛔ An earlier measurement used a chain with **no intergenic regions**, so
`fit_intron_background` was uninformative, ② was silently absent, and the prior looked catastrophic. On
`_gene_in_intergenic`, where ② reads **τ_fac = 161.4** at every intron, the same prior YIELDS. Every
refutability gate here ASSERTS the channel is present before measuring.

⛔ **Two further things a handoff got wrong, corrected by measurement rather than argument:** ``τ_λ`` never
collapsed to zero (4.21e-05 printed at four decimals, and the predicate still fires), and its fall at a
pinned slot is ~98 % the ``[f(1−f)]²`` Jacobian — so the "asymmetry" repair that would feed a prior's
curvature into it was built and REFUSED. See
:func:`test_the_reference_does_not_contribute_to_the_data_s_own_lambda_precision`.

⭐ Each gate carries its own perturbation (`TRAPS: perturb-every-gate`).
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest
from scipy.special import expit

from rigel.types import Strand

from rigel.config import CalibrationConfig
from rigel.calibration.calibrate import _build_intron_prior
from rigel.calibration.density_deconv import density_factor_precision, fit_intron_background
from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.messages.head import HeadPolicy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION
from rigel.calibration.region_geometry import init_beliefs
from rigel.calibration.region_init import build_region_init, strand_evidence
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    mrna_active_strands,
)
from rigel.calibration.simplex_logodds import (
    CompositionPriors,
    _JEFFREYS_REF,
    _NEUTRAL_LOCATION,
    _location_term,
    _logodds_grid,
    _solve_regions_logodds_all,
    structural_reference_location,
)
from rigel.calibration.sweep import solve_chain

from _synthetic import make_chain_parts

_L = 10.0

#: slots of :func:`_mature_exon_chain`. ``N E N E N E N E N`` — 9 slots, regions at the even ones — so the
#: expressed exon is slot 4 and its two flanking PURE-gDNA introns (truth ``f_g = 1``) are slots 2 and 6.
#: ⚠ Deliberately the same fixture `test_sweep.py` uses, because its
#: ``test_mature_no_nascent_hallucination_in_introns`` is the test that found this.
MX_EXON, MX_INTRONS = 4, [2, 6]


def _delta_pmf(length, size=None):
    p = np.zeros((size or length) + 1, dtype=np.float64)
    p[length] = 1.0
    return p


#: the fixtures' shared geometry, so a chain states only what it varies
_SIZE = 2000.0
_GDNA_FL, _RNA_FL = _delta_pmf(300), _delta_pmf(200)
_EFF_G = float(contained_eff_length(np.full(1, _SIZE), _GDNA_FL)[0])
_EFF_R = float(contained_eff_length(np.full(1, _SIZE), _RNA_FL)[0])
_CROSS_G = float(
    crossing_eff_length(_GDNA_FL, np.full(1, UNBOUNDED_REACH), np.full(1, UNBOUNDED_REACH))[0]
)
_CROSS_R = float(
    crossing_eff_length(_RNA_FL, np.full(1, UNBOUNDED_REACH), np.full(1, UNBOUNDED_REACH))[0]
)
#: `intergenic | exon | INTRON | exon | INTRON | exon | intergenic` — per-REGION masks
_IS_EXON = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0])
_IS_INTRON = np.array([0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0])


def _mature_exon_chain(*, rho_g=0.5, rho_m=1.0, kappa=0.95):
    """``exon+ | intron+ | EXON+ | intron+ | exon+`` — a pure-MATURE expressed gene with NO nascent, so the
    two introns hold gDNA and nothing else and their truth is exactly ``f_g = 1``."""
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    length = 2000.0
    unb = np.full(1, UNBOUNDED_REACH)
    eff_g = float(contained_eff_length(np.full(1, length), gdna_fl)[0])
    eff_r = float(contained_eff_length(np.full(1, length), rna_fl)[0])
    cross_g = float(crossing_eff_length(gdna_fl, unb, unb)[0])
    cross_r = float(crossing_eff_length(rna_fl, unb, unb)[0])
    g_half = rho_g * eff_g / 2.0
    mat = rho_m * eff_r
    is_exon = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
    parts = make_chain_parts(
        [BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS],
        region_size_bp=length,
        region_pos=g_half + mat * is_exon,
        region_neg=g_half,
        boundary_pos=rho_g * cross_g / 2.0,
        boundary_neg=rho_g * cross_g / 2.0,
        boundary_spliced=0.0,
        sj=[
            (0, 2, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, rho_m * cross_r),
            (2, 4, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, rho_m * cross_r),
        ],
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    belief = init_beliefs(
        parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _nascent_chain(*, rho_g, rho_m, rho_n, kappa=0.5):
    """The same five regions, but with NASCENT RNA at density ``rho_n`` inside the introns — the one thing
    the 16-condition ladder cannot produce, because it holds ``nrna = 0`` on every row.

    Returns ``(parts, belief, truth_per_slot)``; the truth is derived from the numbers deposited, not
    assumed."""
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    length = 2000.0
    unb = np.full(1, UNBOUNDED_REACH)
    eff_g = float(contained_eff_length(np.full(1, length), gdna_fl)[0])
    eff_r = float(contained_eff_length(np.full(1, length), rna_fl)[0])
    cross_g = float(crossing_eff_length(gdna_fl, unb, unb)[0])
    cross_r = float(crossing_eff_length(rna_fl, unb, unb)[0])
    is_exon = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
    mature, nascent = rho_m * eff_r * is_exon, rho_n * eff_r * (1.0 - is_exon)
    parts = make_chain_parts(
        [BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS],
        region_size_bp=length,
        region_pos=rho_g * eff_g / 2.0 + mature + nascent,
        region_neg=rho_g * eff_g / 2.0,
        boundary_pos=rho_g * cross_g / 2.0 + rho_n * cross_r,  # nascent CROSSES exon|intron
        boundary_neg=rho_g * cross_g / 2.0,
        boundary_spliced=0.0,
        sj=[
            (0, 2, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, rho_m * cross_r),
            (2, 4, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, rho_m * cross_r),
        ],
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    truth = np.empty(9)
    reg_g, reg_r = np.full(5, rho_g * eff_g), mature + nascent
    bnd_g, bnd_r = np.full(4, rho_g * cross_g), np.full(4, rho_n * cross_r)
    truth[0::2] = reg_g / np.maximum(reg_g + reg_r, 1e-30)
    truth[1::2] = bnd_g / np.maximum(bnd_g + bnd_r, 1e-30)
    belief = init_beliefs(
        parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60
    )
    return parts, belief, truth


def _own(args, *, structural: bool):
    """The per-slot message-free SELF-SOLVE, with and without the structural location — no sweep, no
    relay, so ``τ_λ`` is read where it is built rather than after anything can mask it."""
    chain, statics, geometry, belief, _ra = args
    loc = structural_reference_location(statics, _L) if structural else None
    return build_region_init(
        chain,
        statics,
        geometry,
        kappa=0.95,
        od_g=0.0,
        od_r=0.0,
        n_gdna_obs=10000.0,
        n_rna_obs=10000.0,
        n_grid=60,
        logodds_window=_L,
        n_tilt=None,
        n_grid_ss=256,
        belief=belief,
        priors=CompositionPriors(location=loc),
    )


def _sweep(args, *, structural: bool, policy):
    chain, statics, geometry, belief, region_arrays = args
    return solve_chain(
        chain,
        statics,
        geometry,
        belief,
        region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10000.0,
        n_gdna_obs=10000.0,
        n_grid=60,
        logodds_window=_L,
        structural_reference=structural,
        policy=policy,
    )


# ── the CLAIM: what the builder asserts, and where it says nothing ──────────────────────────────────


def test_the_location_is_neutral_exactly_where_mature_rna_can_be():
    """⭐⭐ The builder touches ONLY the annotation-determined slots. Where mature RNA can be, it returns
    the neutral ½ — at which :func:`_location_term` is identically constant and drops out of ψ — so the
    population the gDNA landscape is fitted to serve is left entirely alone.

    ⛔ The perturbation is the other half: at the slots it DOES speak about, the value must not be ½, or
    the builder is inert and every number below is about dead code."""
    _chain, statics, *_ = _mature_exon_chain()
    m = structural_reference_location(statics, _L)
    mature = np.asarray(statics.mrna_active_pos, bool) | np.asarray(statics.mrna_active_neg, bool)
    assert np.all(m[mature] == 0.5)
    assert np.all(m[~mature] > 0.5)
    # ⛔ and on this fixture BOTH populations are non-empty, so neither assertion is vacuous
    assert mature.any() and (~mature).any()


def _strength_of(m: float, *, L: float = 60.0) -> float:
    """The location term's full RANGE over the λ interval, in nats — how much the prior can move ψ. Read
    off the shipped term rather than restated, so it cannot drift from it. ``L`` defaults wide so the
    identity is measured on the term itself; pass the real window to see the lattice cap bite."""
    lam, _ = _logodds_grid(200001, L)
    term = _location_term(lam, np.array([m]))[0]
    return float(term.max() - term.min())


def test_the_location_on_the_log_odds_scale_is_the_strength_in_nats():
    """⭐⭐⭐ **THE IDENTITY THE WHOLE DESIGN TURNS ON: ``strength = logit(m)``, exactly.** The term is
    ``−log[(1−m)f_g + m(1−f_g)]``; at ``f_g → 1`` it is ``−log(1−m)`` and at ``f_g → 0`` it is ``−log m``,
    so its full range is ``log(m/(1−m))``. **The location written on the λ scale IS its strength in nats.**

    ⭐ That is what makes the strength choosable without a constant: write down the weight the prior is
    entitled to and read off the location.

    ⚠ The one exception is the lattice's own top point ``m = σ(L)``, where ``f_g`` cannot exceed ``m`` and
    both halves of the bracket coincide at ``2m(1−m)`` — there the range saturates at ``L − log 2``. That
    is the SHIPPED-BEFORE value and it is why the cap is a cap and not the choice."""
    for k in (0.5, 1.0, 2.0, 3.0, 5.0):
        assert _strength_of(float(expit(k))) == pytest.approx(k, rel=1e-6), k
    # ⛔ the perturbation: at the lattice's top point the identity SATURATES, which is the whole reason
    #   `σ(L)` was the wrong thing to choose. The saturated value is the closed form
    #   ``log[(m²+(1−m)²)/(2m(1−m))]``, i.e. ``L − log 2`` to within ``e^(−L)``.
    for L in (5.0, 10.0, 20.0):
        m = float(expit(L))
        closed = float(np.log((m * m + (1 - m) ** 2) / (2 * m * (1 - m))))
        assert _strength_of(m, L=L) == pytest.approx(closed, rel=1e-9), L
        assert _strength_of(m, L=L) == pytest.approx(L - np.log(2.0), abs=2.0 * np.exp(-L) + 1e-7)
        assert _strength_of(m, L=L) < L


def test_the_strength_is_one_pseudo_observation_derived_from_the_reference_s_own_exponents():
    """⛔⛔⛔ **THE PRIOR IS WORTH EXACTLY ONE PSEUDO-OBSERVATION OF gDNA, AND THAT IS WRITTEN AS THE MEAN
    ONE WOULD PRODUCE** — ``Beta(a,b) → Beta(a+1,b)``, mean ``(a+1)/(a+b+1)`` = **0.75** at
    ``a = b = _JEFFREYS_REF``. No new number: it is the reference's own exponents, and it moves if they do.

    ⚠ **Units matter here and an earlier draft got them wrong.** ``strength = logit(m)``, so a strength is
    a LOG-ODDS RANGE in nats while ``a + b`` is a PSEUDO-COUNT. Setting ``m = σ(a+b)`` equated the two —
    numerically close (0.731 vs 0.75) but an analogy, not a derivation. The route above stays in
    composition-mean units throughout.

    ⛔ **Both replace ``m = σ(L)``, the lattice's own width, which was measured WORSE THAN NO PRIOR AT
    ALL.** Swept against both obligations — DELIVER where the claim is true, YIELD where evidence refutes
    it — ``Σ|f_g − truth|`` over a depth ladder and both κ::

        no prior              deliver 1.2286   refute 0.3946   total 1.6231
           0.69 nats          deliver 0.9477   refute 0.2606   total 1.2083
        ⭐ log 3 (m = 0.75)    deliver 0.7849   refute 0.3293   total 1.1142   ← SHIPPED
           1.50 nats          deliver 0.6476   refute 0.5346   total 1.1822
        ⛔ 9.31 (σ(L))         deliver 0.0037   refute 2.0247   total 2.0285   ← worst of seven

    ⭐ **The optimum is BROAD — 0.69 → 1.50 nats all within 6 % — which is what makes a DERIVED value safe
    rather than lucky.** ⚠ The 16-condition ladder cannot rank this and that is why it picked the worst
    row: it holds ``nrna = 0``, so it scores the DELIVER column alone, where more nats is monotonically
    better."""

    class _S:
        mrna_active_pos = np.array([True, False])
        mrna_active_neg = np.array([False, False])

    a = b = _JEFFREYS_REF
    expected = (a + 1.0) / (a + b + 1.0)
    assert expected == 0.75  # exactly, at the shipped exponents
    for L in (5.0, 10.0, 20.0):
        m = float(structural_reference_location(_S, L)[1])
        assert m == expected
        # ⭐ the odds are e^strength, so 0.75 IS a 3:1 claim
        assert _strength_of(m) == pytest.approx(float(np.log(3.0)), rel=1e-6), L
        # ⛔ and it does NOT track the window — that is the property being removed
    # ⭐ the cap survives as a CAP: a window narrower than the claim binds
    narrow = float(structural_reference_location(_S, 0.25)[1])
    assert narrow == pytest.approx(float(expit(0.25)), rel=1e-12)
    assert narrow < expected
    # ⛔ the perturbation: the neutral half is untouched by any of this
    assert structural_reference_location(_S, 10.0)[0] == 0.5


def _overturn(loc, u_pos, u_neg, *, kappa=0.99):
    """``|f_g(with the location) − f_g(without)|`` on one single-strand slot — how much of the prior
    survives the data."""

    def s(location):
        return float(
            np.asarray(
                _solve_regions_logodds_all(
                    np.array([u_pos]),
                    np.array([u_neg]),
                    np.array([True]),
                    np.array([False]),
                    np.array([u_pos + u_neg]),
                    np.array([0.0]),
                    kappa=kappa,
                    od_g=0.0,
                    od_r=0.0,
                    n_grid=60,
                    n_grid_ss=256,
                    priors=None
                    if location is None
                    else CompositionPriors(location=np.array([location])),
                ).gdna_frac
            )[0]
        )

    return abs(s(loc) - s(None))


def test_the_overturn_budget_holds_only_where_the_likelihood_is_informative():
    """⭐ The strength converts into FRAGMENTS — ``strength/log(2κ)`` = **13.6** at ``L = 10, κ = 0.99`` —
    which is what makes it a domain statement rather than an opaque epsilon.

    ⚠ **The budget assumes every fragment delivers ``log(2κ)`` nats ON THE λ AXIS, and near the vertex it
    does not** — the strand Fisher information is ``I ∝ [f_g(1−f_g)]²``, which collapses exactly where this
    prior points. That is why the strength cannot be chosen from the lattice: at the old ``σ(L)`` the
    residual reached **0.0058 at 5,000 fragments**, 370× the budget's depth, in that flat regime.

    ⭐ **At the declared ``a + b`` weight the budget is 1.46 fragments and the pathology is gone.**
    Measured on a pure-sense slot: **0.2112** at 1 fragment, 0.0212 by 10, 0.0019 by 100 — and in the flat
    5 %-antisense regime the residual never exceeds **0.006**, i.e. below one ``K = 60`` grid step (0.085).
    A prior that cannot move the answer by a representable amount cannot outvote anything."""
    budget = lambda kappa: (2.0 * _JEFFREYS_REF) / np.log(2.0 * kappa)  # noqa: E731
    assert budget(0.99) == pytest.approx(1.46, abs=0.02)
    # a WEAKER library needs MORE fragments to say the same thing — monotone in κ, and the direction is
    # the half a swapped ratio would get wrong
    assert budget(0.95) > budget(0.99)
    with np.errstate(divide="ignore"):
        assert not np.isfinite(budget(0.5))

    m = float(expit(2.0 * _JEFFREYS_REF))
    # ⭐ audible at one fragment, spent within ten — a one-pseudo-observation prior behaving like one
    assert _overturn(m, 1.0, 0.0) > 0.15
    assert _overturn(m, 10.0, 0.0) < 0.03
    assert _overturn(m, 100.0, 0.0) < 0.005
    # ⛔ the perturbation, and the whole reason the strength moved: in the FLAT regime the old lattice-wide
    #    location still shifted a 5,000-fragment slot, and the declared-weight one stays under a grid step
    old = float(expit(_L))
    assert _overturn(old, 4750.0, 250.0) > 0.003, _overturn(old, 4750.0, 250.0)
    assert _overturn(m, 4750.0, 250.0) < 0.085, _overturn(m, 4750.0, 250.0)


def test_the_strength_is_independent_of_the_window_and_of_the_clamp():
    """⭐ Because the strength is now the reference's own weight rather than the lattice's, the builder is
    **invariant to ``L``** over every window a config could hold — which also puts it far away from
    :func:`_location_term`'s ``_EPS`` clamp.

    ⚠ That clamp is why ``m = σ(L)`` could not have stayed: above ``L ≈ 20.72``, ``1 − σ(L)`` falls under
    it and the strength saturates at ``log(1/2ε) = 20.723`` forever — measured **20.709 at L = 25** against
    a claimed 24.307, and **20.723 at L = 37** against 36.307. The old rule's own justification broke on a
    wide window; this one has nothing to break."""

    class _S:
        mrna_active_pos = np.array([True, False])
        mrna_active_neg = np.array([False, False])

    ref = structural_reference_location(_S, 10.0)
    # ⭐ invariant on every window WIDER than the claim itself (log 3 = 1.0986 nats)
    for L in (1.5, 5.0, 20.0, 25.0, 40.0, 100.0):
        assert np.array_equal(structural_reference_location(_S, L), ref), L
    # ⛔ and BELOW it the lattice cap correctly binds — the claim is trimmed to what the grid can hold
    assert float(structural_reference_location(_S, 1.0)[1]) < float(ref[1])
    # ⛔ the perturbation: the OLD rule was not invariant, and above the clamp it silently saturated at
    #   ``logit(1 − _EPS)`` FOREVER — 20.7233 at L = 25 and at L = 37 alike, against a claimed 24.31/36.31
    saturated = _strength_of(float(expit(25.0)))
    assert saturated == pytest.approx(20.7233, abs=0.01)
    assert _strength_of(float(expit(37.0))) == pytest.approx(saturated, abs=1e-9)
    assert saturated < 25.0 - np.log(2.0)


# ── the MASK: which slots the claim is made at, over the WHOLE signature space ──────────────────────


def _boundary_shape(sig_l: int, sig_r: int) -> str:
    """The owner's four classes, plus the two shapes that would WIDEN them — named so a gate can say which
    one an annotation produced rather than merely that the count moved."""
    ex = lambda s: bool(s & (BIT_EXON_POS | BIT_EXON_NEG))  # noqa: E731
    if sig_l == 0 or sig_r == 0:
        return "gene edge"
    if ex(sig_l) != ex(sig_r):
        return "one flank exonic"
    if not ex(sig_l):
        return "intron|intron"
    return "exon|exon, opposite strands"


def test_the_mask_is_exactly_the_four_classes_plus_two_named_wideners():
    """⛔⛔ **EXHAUSTIVE OVER THE WHOLE SIGNATURE SPACE, NOT OVER ONE ANNOTATION.** ``¬mrna_active`` is
    what the builder asserts on, and the owner's justification — *"no annotated transcription there, so
    assume 100 % gDNA until proven otherwise"* — is about FOUR classes: intergenic REGION, intron REGION,
    gene-edge BOUNDARY, one-flank-exonic BOUNDARY.

    ⭐ Two further shapes satisfy the predicate and are **absent from the human index measured**
    (``1,312 + 2,620 + 9,805 + 19,610 = 33,347`` of 70,176 slots, no remainder) — an ``intron|intron``
    BOUNDARY, and an ``exon|exon`` BOUNDARY whose flanks are exons of OPPOSITE strands. ⛔ **Neither is
    impossible on another annotation, and the opposite-strand one has annotated transcription on BOTH
    sides, so the justification would not hold for it.** This enumerates all ``N_SIGNATURES²`` flank pairs
    so the population is stated rather than assumed, and NAMES the two so a future annotation producing one
    is identified rather than merely counted.

    ⚠ The perturbation is that both wideners are REACHED here: if the enumeration could not construct them,
    this test would be asserting a partition of the empty set."""
    sigs = range(N_SIGNATURES)
    # a REGION's mask is its OWN exon bits, so ¬mature ⟺ no exon bit at all — intergenic or intron-only
    for s in sigs:
        mp, mn = mrna_active_strands(np.array([s], dtype=np.uint8))
        mature = bool(mp[0] or mn[0])
        assert mature == bool(s & (BIT_EXON_POS | BIT_EXON_NEG)), s
        if not mature:
            assert (s == 0) or bool(s & (BIT_INTRON_POS | BIT_INTRON_NEG)), s

    # a BOUNDARY's mask is the per-strand AND over its two flanks — mature must CROSS, not merely be near
    seen: dict[str, int] = {}
    for sl, sr in itertools.product(sigs, sigs):
        lp, ln = mrna_active_strands(np.array([sl], dtype=np.uint8))
        rp, rn = mrna_active_strands(np.array([sr], dtype=np.uint8))
        if bool((lp[0] and rp[0]) or (ln[0] and rn[0])):
            continue  # mature crosses ⇒ the builder says nothing here
        shape = _boundary_shape(sl, sr)
        seen[shape] = seen.get(shape, 0) + 1
    assert set(seen) == {
        "gene edge",
        "one flank exonic",
        "intron|intron",
        "exon|exon, opposite strands",
    }, seen
    assert all(v > 0 for v in seen.values()), seen


# ── ⛔⛔⛔ τ_λ IS THE DATA'S INFORMATION AND THE REFERENCE MAY NOT ADD TO IT ──────────────────────────


def test_the_reference_does_not_contribute_to_the_data_s_own_lambda_precision():
    """⛔⛔⛔ **THE RULING, AND IT IS THE OPPOSITE OF THE ONE THIS FILE WAS FIRST WRITTEN WITH.** ``τ_λ`` is
    the **DATA's** Fisher information on λ. The reference's location is a prior, it carries no count, and
    it must contribute NOTHING here.

    ⭐ The design this replaces argued the other way: a slot the reference pins sits at the grid's edge
    where the strand term vanishes, so ``τ_λ`` falls **3,227×** (0.1358 → 4.21e-05) exactly when the slot's
    belief becomes most certain; the location term is a λ-factor like ``intron_prior``, so let it pay its
    curvature too (``density_factor_precision``, which does read **0.743** off it). ⛔ **It was built,
    measured and REFUSED, on three independent findings — this gate is what stops it coming back:**

    1. **the 3,227× is not a loss.** ``I_strand ∝ [f_g(1−f_g)]²/(4p(1−p))``, and that Jacobian alone
       predicts **3,154×** between ``f_g = 0.98576`` and ``0.99975`` — ~98 % of it. The strand likelihood
       genuinely IS that flat on λ at the point the prior chose. Nothing moved anywhere.
    2. **it is a BOOLEAN GATE FLIP, not a contribution.** At the vertex ``Var(log f_g) = (1−f_g)²/τ`` is
       ~8e-08, so :func:`own_precision` saturates at the COUNT ceiling: τ = 0.029 and τ = 1e6 both return
       850.44 against a ceiling of 850.50. Only ``τ > 0`` does any work, and what it releases is 850
       fragments of count precision — from a prior worth one.
    3. **it credits DATA-FREE slots.** The location carries no count, so an ``n = 0`` slot goes
       ``prec_g`` 0 → 0.2026. ⚠ And the structural reference's own safety argument is that all four of its
       classes are EMPTY in a zero-gDNA library — precisely the population that would gain manufactured
       evidence. ``intron_prior`` is not the precedent it looks like: its NegBinom curvature is
       count-derived and self-limits on thin data.

    ⚠ ``has_own_composition_evidence`` is also 0.8.0's own denominator (`solvability_audit`), so flipping
    it on 33,347 slots silently redefines what "confidently wrong" counts.
    """
    args = _mature_exon_chain()
    chain, statics, geometry, belief, _ra = args
    off, on = _own(args, structural=False), _own(args, structural=True)

    # ⛔ THE GATE: τ_λ is EXACTLY the data's own sources at the location's own solve point — nothing added.
    i_strand, _lock = strand_evidence(
        np.asarray(geometry.unspliced_count, np.float64)[:, 0],
        np.asarray(geometry.unspliced_count, np.float64)[:, 1],
        on.f_g,
        kappa=0.95,
        od_g=0.0,
        od_r=0.0,
        n_gdna_obs=10000.0,
        n_rna_obs=10000.0,
        is_region=np.asarray(chain.kind) == REGION,
        locked=~(
            (np.asarray(statics.free_pos, bool) | np.asarray(statics.free_neg, bool))
            & (np.asarray(geometry.unspliced_count, np.float64).sum(axis=1) > 0.0)
        ),
    )
    single = np.asarray(statics.free_pos, bool) ^ np.asarray(statics.free_neg, bool)
    assert np.array_equal(on.tau_lam, np.where(single, i_strand, 0.0)), (
        "τ_λ carries something other than the strand + factory data terms — if the reference's "
        "location has been given a contribution here, read this test's docstring before restoring it."
    )
    # ⭐ and the FALL is the Jacobian, so it is not evidence of anything missing
    intr = np.asarray(MX_INTRONS)
    a, b = float(off.f_g[intr][0]), float(on.f_g[intr][0])
    jac = ((a * (1 - a)) ** 2) / ((b * (1 - b)) ** 2)
    ratio = float(off.tau_lam[intr][0] / on.tau_lam[intr][0])
    assert ratio == pytest.approx(jac, rel=0.05), (ratio, jac)
    # ⛔ the perturbation: the refused contribution is REAL and would be visible here, so this gate is
    #   not passing because the quantity is unreachable
    lam, _ = _logodds_grid(60, _L)
    refused = density_factor_precision(
        _location_term(lam, structural_reference_location(statics, _L)), lam
    )
    # ⚠ 0.0329 at the declared ``a + b`` weight (it was 0.743 at the old lattice-wide location) — small,
    #   but the defect was never the magnitude: it is that ANY nonzero value flips
    #   `has_own_composition_evidence` and releases the full COUNT precision.
    assert float(refused[intr][0]) > 0.01, refused
    assert not np.array_equal(on.tau_lam, on.tau_lam + refused)


def test_a_neutral_location_is_a_bit_identical_no_op_through_psi():
    """⛔⛔ **PROPERTY 1 IS ENFORCED IN FLOAT64, NOT MERELY TRUE ON PAPER — and the gap was 0.0845.**

    At ``m = ½`` the bracket is the constant ½ in exact arithmetic, but ``logaddexp`` leaves the row with
    ``ptp = 2.22e-16`` over three distinct values at every grid size. `_posterior_median_fg` reads the grid
    point where the CDF first reaches ½, so at a slot whose posterior is exactly symmetric that rounding
    tips the knife-edge and ``f_g`` moves **0.5423 → 0.4577** — one full grid step, at a slot the
    structural reference claims to say nothing about, and the ~37,000 mature slots are exactly the
    population the gDNA landscape is fitted to serve.

    ⭐ :func:`_location_term` therefore returns an EXACT zero at the neutral location. ⛔ The perturbation
    is the near-constant the fix replaces: it must NOT be bit-identical, or there was nothing to fix."""
    lam, _ = _logodds_grid(60, _L)
    assert np.all(_location_term(lam, np.array([_NEUTRAL_LOCATION]))[0] == 0.0)

    def solve(loc, u_pos, u_neg, *, allow_neg=True, kappa=0.5):
        return float(
            np.asarray(
                _solve_regions_logodds_all(
                    np.array([u_pos]),
                    np.array([u_neg]),
                    np.array([True]),
                    np.array([allow_neg]),
                    np.array([u_pos + u_neg]),
                    np.array([0.0]),
                    kappa=kappa,
                    od_g=0.0,
                    od_r=0.0,
                    n_grid=60,
                    n_grid_ss=256,
                    priors=None if loc is None else CompositionPriors(location=np.array([loc])),
                ).gdna_frac
            )[0]
        )

    # ⭐ the knife-edge case: a BALANCED AMBIG slot at κ=½ — an exactly symmetric posterior
    for u, kw in ((100.0, {}), (3.0, {}), (5.0, {"allow_neg": False, "kappa": 0.95})):
        base = solve(None, u, u if kw.get("allow_neg", True) else 0.0, **kw)
        neutral = solve(_NEUTRAL_LOCATION, u, u if kw.get("allow_neg", True) else 0.0, **kw)
        assert neutral == base, (u, kw, base, neutral)
    # ⛔ the perturbation: a location a hair off neutral — the same magnitude of claim the float64
    #   rounding used to make by accident — DOES move that slot, so the equalities above are the exact
    #   zero doing its job and not the knife-edge being unreachable
    nudged = solve(_NEUTRAL_LOCATION + 1e-12, 100.0, 100.0)
    assert nudged != solve(None, 100.0, 100.0), nudged


def test_a_pure_gdna_intron_rises_when_the_prior_says_it_is_pure_gdna():
    """⭐⭐ **THE CONSEQUENCE, ON THE SHIPPED PATH.** The two introns hold gDNA and nothing else, so their
    truth is exactly ``f_g = 1``; a prior asserting ``m = σ(L)`` there can only move them UP. Measured
    under the shipped `SilentPolicy`: **0.9858 → 0.9998**.

    ⛔ The perturbation is the expressed EXON in the same chain: it is mostly RNA, the builder says
    nothing about it, and it must not move at all — a location leaking onto mature-crossing slots, or a
    neutral row that is not bit-neutral, would both show up here."""
    args = _mature_exon_chain()
    off = _sweep(args, structural=False, policy=SilentPolicy())
    on = _sweep(args, structural=True, policy=SilentPolicy())
    intr = np.asarray(MX_INTRONS)
    assert np.all(on.f_g[intr] > off.f_g[intr]), (off.f_g[intr], on.f_g[intr])
    assert float(on.f_g[MX_EXON]) == float(off.f_g[MX_EXON])


@pytest.mark.xfail(
    strict=True,
    reason="⛔ A RELAY DEFECT THE REFERENCE EXPOSES, NOT ONE IT CAUSES — and `message_propagation` is "
    "OFF, so this is a STUDY configuration. Under `HeadPolicy` the same intron goes 0.9006 → 0.7661, "
    "the wrong way, and it is the λ-message that does it (nulling `lam_channel` restores 0.9006 "
    "exactly; `cm_g`/`cm_p` stay 0). The reference's claim is CORRECT at every slot it makes one — true "
    "f_g is 1.0000 there — and the relay then transports that correct claim ACROSS an exon↔intron "
    "boundary, where the mature-RNA population changes. The composition licence knows about transcript "
    "TERMINI (`terminus_flank_gain`) and not about `mrna_active` flipping, which is exactly the "
    "predicate that says the RNA population differs across this hop. ⛔ Do not close this by softening "
    "the prior (measured worse) or by topping up τ_λ (refused — see "
    "`test_the_reference_does_not_contribute_to_the_data_s_own_lambda_precision`).",
)
def test_the_relay_does_not_carry_a_structural_claim_across_a_population_change():
    """⛔⛔ The relay-ON twin of the gate above, written as an executable record of a PROVEN defect rather
    than as a bound to be widened.

    ⚠ Before the reference existed this pathway had never been exercised at ``κ = ½``: measured,
    `HeadPolicy`/OFF is **bit-identical to `SilentPolicy`/OFF on all 9 slots**, because every slot's
    ``τ_λ`` is 0 there and no message carries precision. The reference supplies the first confident
    neighbour the unstranded relay has ever had."""
    args = _mature_exon_chain()
    off = _sweep(args, structural=False, policy=HeadPolicy())
    on = _sweep(args, structural=True, policy=HeadPolicy())
    intr = np.asarray(MX_INTRONS)
    assert np.all(on.f_g[intr] > off.f_g[intr]), (off.f_g[intr], on.f_g[intr])


# ── ⛔⛔⛔ WHY THE DEFAULT IS OFF — the claim is FALSE where nascent RNA lives ───────────────────────


def _gene_in_intergenic(*, rho_g, rho_m, rho_n, kappa, depth=1.0):
    """``intergenic | exon | INTRON | exon | INTRON | exon | intergenic`` — 7 regions, 13 slots, introns
    at 4 and 8.

    ⛔⛔ **THE INTERGENIC FLANKS ARE LOAD-BEARING AND THEIR ABSENCE INVALIDATED AN EARLIER MEASUREMENT.**
    `fit_intron_background` pools INTERGENIC regions only, so on a chain without them it returns
    uninformative, `_build_intron_prior` returns ``None``, and the intron-vs-intergenic density mechanism
    is **silently absent** — leaving the prior with nothing to be refuted by and making it look
    catastrophic. Production always has intergenic regions."""
    g_half = depth * rho_g * _EFF_G / 2.0
    mature = depth * rho_m * _EFF_R * _IS_EXON
    nascent = depth * rho_n * _EFF_R * _IS_INTRON  # sense-strand RNA inside the introns
    parts = make_chain_parts(
        [0, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, 0],
        region_size_bp=_SIZE,
        region_pos=g_half + mature + nascent,
        region_neg=g_half,
        boundary_pos=depth * rho_g * _CROSS_G / 2.0,
        boundary_neg=depth * rho_g * _CROSS_G / 2.0,
        boundary_spliced=0.0,
        sj=[
            (1, 3, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, depth * rho_m * _CROSS_R),
            (3, 5, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, depth * rho_m * _CROSS_R),
        ],
        gdna_fl=_GDNA_FL,
        rna_fl=_RNA_FL,
    )
    reg_g = np.full(7, depth * rho_g * _EFF_G)
    truth_intron = float(reg_g[2] / max(reg_g[2] + (mature + nascent)[2], 1e-30))
    belief = init_beliefs(
        parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60
    )
    return parts, belief, truth_intron


def _solve_with_factory(parts, belief, kappa, *, structural):
    """The production shape: `intron_factory` is ON by default, so mechanism ② is live."""
    eff_g = contained_eff_length(np.full(7, _SIZE), _GDNA_FL)
    bg = fit_intron_background(parts.substrate, parts.region_arrays, eff_g, include_introns=False)
    assert bg.informative, "the intergenic pool is uninformative — mechanism ② would be absent"
    ip = _build_intron_prior(
        parts.chain, parts.substrate, parts.region_arrays, eff_g, CalibrationConfig(), bg=bg
    )
    return solve_chain(
        parts.chain,
        parts.statics,
        parts.geometry,
        belief,
        parts.region_arrays,
        rna_sense_frac=kappa,
        n_rna_obs=1e5,
        n_gdna_obs=1e5,
        n_grid=60,
        logodds_window=_L,
        n_grid_ss=256,
        intron_prior=ip,
        structural_reference=structural,
        policy=SilentPolicy(),
    )


def test_the_prior_yields_to_the_density_mechanism_on_unstranded_data():
    """⭐⭐⭐ **THE GATE THAT RANKS THE STRENGTH, AND THE PANEL CANNOT WRITE IT.** *Nascent RNA is zero
    unless proven otherwise* is the correct prior; the obligation it carries is that the two mechanisms
    which CAN prove otherwise actually overturn it.

    ② is the intron-vs-intergenic density factor, it is ``intron_factory = True`` in production, and — the
    part that matters — **it works UNSTRANDED**, where the strand channel is identically dead. Measured on
    this fixture it carries ``τ_fac = 161.4`` at every intron slot.

    With the prior at its declared ``a + b`` weight the answer is the evidence's, not the prior's:
    at ``ρ_n = 0.25`` (truth 0.6539) the solve reads within 0.02 of the no-prior answer, and at
    ``ρ_n = 1.0`` (truth 0.3208) it is the same to 4 dp.

    ⛔ The perturbation is the CONTROL in the same fixture: at ``nascent = 0`` (truth exactly 1.0) the
    prior must still DELIVER, and it must beat the no-prior solve — otherwise "yields" would just mean
    "does nothing"."""
    for kappa in (0.5, 0.99):
        for rho_n in (0.25, 1.0):
            p, b, truth = _gene_in_intergenic(rho_g=0.5, rho_m=1.0, rho_n=rho_n, kappa=kappa)
            off = float(np.mean(_solve_with_factory(p, b, kappa, structural=False).f_g[[4, 8]]))
            on = float(np.mean(_solve_with_factory(p, b, kappa, structural=True).f_g[[4, 8]]))
            assert abs(on - truth) < abs(off - truth) + 0.03, (kappa, rho_n, truth, off, on)
            assert abs(on - truth) < 0.06, (kappa, rho_n, truth, on)
    # ⛔ the perturbation: where the claim is TRUE the prior must never make things WORSE, and on a
    #   thin slot — where it is supposed to be worth having — it must strictly help. ⚠ On a slot the
    #   factory has already pinned, a one-pseudo-observation prior correctly cannot move a grid point.
    for kappa in (0.5, 0.99):
        for depth, strict in ((0.01, True), (1.0, False)):
            p, b, truth = _gene_in_intergenic(
                rho_g=0.5, rho_m=1.0, rho_n=0.0, kappa=kappa, depth=depth
            )
            assert truth == 1.0
            off = float(np.mean(_solve_with_factory(p, b, kappa, structural=False).f_g[[4, 8]]))
            on = float(np.mean(_solve_with_factory(p, b, kappa, structural=True).f_g[[4, 8]]))
            assert on > off if strict else on >= off, (kappa, depth, off, on)


@pytest.mark.parametrize("rho_n", (0.25, 1.0), ids=("nascent-quarter", "nascent-full"))
def test_with_no_refutation_channel_at_all_the_prior_is_the_only_voice(rho_n):
    """⚠ **THE DEGENERATE CASE, KEPT AS A BOUNDARY RATHER THAN AS AN ALARM.** This chain has NO
    intergenic regions, so `fit_intron_background` has no pool and mechanism ② cannot exist; at ``κ = ½``
    mechanism ① is dead by derivation too. **Neither refutation channel exists**, and there the prior is
    the only voice and wins — which is correct Bayes, not a defect.

    ⭐ It is here because an earlier session measured exactly this fixture, found the prior "2.10–5.74×
    worse with nascent", and concluded the claim was unsafe. The conclusion was an artefact of the missing
    intergenic flanks: on `_gene_in_intergenic`, where ② is live at ``τ_fac = 161.4``, the same prior
    YIELDS to the same nascent. ⛔ **A refutability test is only valid if the refutation channel is in the
    fixture** — check it, do not assume it.

    ⚠ What survives from that measurement is the honest residual: where a slot has no channel at all, the
    prior's accuracy is load-bearing. That is why its STRENGTH is ``a + b`` and not the lattice's width."""
    pinned = np.array([1, 2, 3, 5, 6, 7])

    def err(rho_n_, structural):
        parts, belief, truth = _nascent_chain(rho_g=0.5, rho_m=1.0, rho_n=rho_n_)
        out = solve_chain(
            parts.chain,
            parts.statics,
            parts.geometry,
            belief,
            parts.region_arrays,
            rna_sense_frac=0.5,
            n_rna_obs=10000.0,
            n_gdna_obs=10000.0,
            n_grid=60,
            logodds_window=_L,
            structural_reference=structural,
            policy=SilentPolicy(),
        )
        mass = np.asarray(parts.geometry.unspliced_count, np.float64).sum(axis=1)
        return float(np.sum(np.abs(np.asarray(out.f_g) - truth)[pinned] * mass[pinned]))

    # ⭐ the precondition this test is ABOUT: with no intergenic regions the channel does not exist
    eff_g = contained_eff_length(np.full(5, _SIZE), _GDNA_FL)
    p0, _b0, _t0 = _nascent_chain(rho_g=0.5, rho_m=1.0, rho_n=rho_n)
    assert not fit_intron_background(
        p0.substrate, p0.region_arrays, eff_g, include_introns=False
    ).informative, "this fixture is supposed to have NO refutation channel; it now has one"
    # ⛔ with no channel the prior STANDS, so the error is whatever the claim's own accuracy is. Measured
    #   at the declared weight: still a net win at ρ_n = 0.25 (160 vs 523) and a net loss by ρ_n = 1.0
    #   (2,389 vs 766) — the crossover is where the claim stops being approximately true.
    win = err(rho_n, True) < err(rho_n, False)
    assert win == (rho_n <= 0.25), (rho_n, err(rho_n, False), err(rho_n, True))
    # ⛔ the perturbation: at nascent = 0 the claim is exactly true and the prior helps at every level
    assert err(0.0, True) < 0.7 * err(0.0, False), (err(0.0, False), err(0.0, True))


# ── the SWITCH: off is not "a neutral term", it is NO term ──────────────────────────────────────────


@pytest.mark.parametrize("policy", (SilentPolicy, HeadPolicy), ids=("silent", "head"))
def test_the_flag_off_is_bit_identical_through_the_whole_solve(policy):
    """⛔⛔ **THE STAGE-1 CONTRACT, AT THE SWITCH RATHER THAN AT THE TERM.**
    ``tests/calibration/test_reference_location.py`` pins ``location=None`` ⇒ no term inside ψ; this pins
    that the shipped ``structural_reference = False`` reaches that state through `solve_chain`, on every
    output array.

    ⭐ The perturbation is the ON run in the same call: it must differ, or the flag is unreachable."""
    args = _mature_exon_chain()
    a = _sweep(args, structural=False, policy=policy())
    b = _sweep(args, structural=False, policy=policy())
    fields = ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg")
    for f in fields:
        assert np.array_equal(getattr(a, f), getattr(b, f)), f
    on = _sweep(args, structural=True, policy=policy())
    assert not np.array_equal(on.f_g, a.f_g)


def test_the_builder_needs_only_the_annotation_and_the_window():
    """⭐ It reads no count, no density and no deconvolved array — which is what makes it admissible at
    pass-0, where the gDNA landscape does not exist yet. Two chains that differ ONLY in their fragment
    counts must produce the identical location.

    ⛔ The perturbation: a chain that differs in its ANNOTATION must produce a different one."""
    a = structural_reference_location(_mature_exon_chain(rho_g=0.5, rho_m=1.0)[1], _L)
    b = structural_reference_location(_mature_exon_chain(rho_g=9.0, rho_m=0.01)[1], _L)
    assert np.array_equal(a, b)
    # ⛔ the perturbation — an all-exon chain has nothing to assert on
    parts = make_chain_parts([BIT_EXON_POS, BIT_EXON_POS, BIT_EXON_POS], region_size_bp=2000.0)
    assert np.all(structural_reference_location(parts.statics, _L)[[0, 2, 4]] == 0.5)
    # …and an all-intron one asserts on every slot
    parts = make_chain_parts(
        [BIT_INTRON_NEG, BIT_INTRON_NEG, BIT_INTRON_NEG], region_size_bp=2000.0
    )
    assert np.all(structural_reference_location(parts.statics, _L) > 0.5)
