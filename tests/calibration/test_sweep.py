"""The belief-propagation sweep (`sweep.solve_chain`) and the beliefs it starts from.

⭐ **Every fixture here is on the S5.e axes** — ``_synthetic.make_chain_parts``, i.e. a region axis, a
contiguous-boundary axis with ``k − 1`` entries per reference and **no terminal slots**, and a sj axis
whose boundaries state their own ``(src, dst, strand)``. The per-FACE fixtures this file used to carry are
gone; ``RegionGeometry``'s own gate is ``test_region_geometry.py``, written from scratch against enumerated
start positions.

⚠ **Two things the old shape hid, and their tests say so in place**: a reference terminal was a
data-free boundary SLOT that could be G1-locked and emit structural all-gDNA into its neighbour
(`test_gdna_sweep_zero_gdna_pin_and_monotone`), and the mature flux at an intron↔exon boundary had to be
placed by hand rather than derived from a sj's endpoints (`_mature_exon_chain`).
"""

from __future__ import annotations


import functools

import numpy as np
import pytest

from rigel.types import Strand

from rigel.calibration.messages.relay import RelayPolicy
from rigel.calibration.sweep import solve_chain

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.region_geometry import init_beliefs

from _synthetic import make_chain_parts
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    N_SIGNATURES,
    mrna_active_strands,
    nrna_active_strands,
)


#: ⚠ These gates exercise HEADPOLICY's operators, so the policy is named EXPLICITLY. ``solve_chain``
#: defaults to ``SilentPolicy``, which sends nothing — every assertion below would then be vacuous, which
#: is TRAPS: could-the-arm-have-fired exactly ("check the arm COULD have changed something").
region_sweep = functools.partial(solve_chain, policy=RelayPolicy())


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def test_init_zero_gdna_introns_via_strand():
    # The P1 gate: a zero-gDNA library. 3 regions: intergenic | intron+ | AMBIG (one ref).
    # intergenic gDNA = strand-symmetric; intron+ RNA = strongly sense-tilted (κ=0.95); AMBIG = symmetric.
    parts = make_chain_parts(
        [0, BIT_INTRON_POS, BIT_EXON_POS | BIT_EXON_NEG],
        region_size_bp=[1000.0, 2000.0, 800.0],
        region_pos=[50.0, 95.0, 50.0],
        region_neg=[50.0, 5.0, 50.0],
    )
    b = init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)

    # ⭐ the chain is N E N E N, so the regions are at 0, 2, 4 — there are no terminal slots.
    rid = [0, 2, 4]
    fg = b.f_g[rid]
    # intergenic: locked gDNA sink {0,0,1}, all precision locked (var 0).
    assert fg[0] == 1.0
    assert b.var_gdna[0] == 0.0 and b.var_pos[0] == 0.0 and b.var_neg[0] == 0.0
    # intron+ (zero gDNA): the strand tilt alone drives f_g → 0; the − axis is locked (var 0), + & g finite.
    assert fg[1] < 0.15
    assert b.var_neg[2] == 0.0 and np.isfinite(b.var_gdna[2]) and np.isfinite(b.var_pos[2])
    # AMBIG: unresolved by strand → {0,0,1} default at MAX (inf) variance for the sweep to resolve.
    assert fg[2] == 1.0
    assert np.isinf(b.var_gdna[4]) and np.isinf(b.var_pos[4]) and np.isinf(b.var_neg[4])


def test_init_boundary_continuity_gate():
    # 1 ref, 2 regions (exon+ | intron+) → ONE boundary between them. ⭐ The two terminal boundary slots the
    # predecessor asserted on do not exist: a reference with k regions owns k-1 boundaries, so there is nothing
    # before the first region or after the last to be a sink.
    parts = make_chain_parts(
        [BIT_EXON_POS, BIT_INTRON_POS],
        region_size_bp=[1000.0, 2000.0],
        region_pos=[80.0, 40.0],
        region_neg=[4.0, 30.0],
        # the crossing: sense-tilted unspliced (κ=0.95 ⇒ +) + a certified-RNA (spliced) floor
        boundary_pos=[90.0],
        boundary_neg=[5.0],
        boundary_spliced=[50.0],
    )
    b = init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)
    # slots: N0=0, E0=1, N1=2.
    # E0 (ex+→in+): +strand continuous (G2+) ⇒ the strand tilt resolves f_g → 0; − axis locked (var 0).
    assert b.f_g[1] < 0.15
    assert b.var_neg[1] == 0.0 and np.isfinite(b.var_gdna[1])


def test_init_tss_boundary_is_black_hole():
    # intergenic | exon+ : the internal boundary is a TSS (intergenic↔exon) ⇒ continuity blocks RNA ⇒ sink.
    # the TSS-crossing fragments are sense-tilted, but continuity must STILL block RNA (the black hole).
    parts = make_chain_parts(
        [0, BIT_EXON_POS],
        region_size_bp=[1000.0, 2000.0],
        region_pos=[50.0, 80.0],
        region_neg=[50.0, 4.0],
        boundary_pos=[90.0],
        boundary_neg=[5.0],
    )
    b = init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)
    # slot 1 is the TSS boundary: a locked gDNA sink despite the sense tilt (all precision locked at 0).
    assert b.f_g[1] == 1.0 and b.var_gdna[1] == 0.0 and b.var_pos[1] == 0.0 and b.var_neg[1] == 0.0


def test_precision_state_count_resolution():
    """The log-density solver's precision state is ``Var(log f_g)`` — the message currency (TRAPS: two-gaussians-one-latent). It reflects
    EVIDENCE: a region with more fragments (same composition) resolves its log-density sharper, so a lower
    ``Var(log f_g)`` ⇒ a more confident message. (In LOG space a confident ``f_g→0`` region has WIDE variance —
    a near-zero gDNA density carries little reliable gDNA-density information to impute, the
    "zero-density-is-not-a-measurement" principle — so the lattice's linear ``Var(f_g)`` ordering does not
    carry over.) A region with no fragments reports zero variance."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    kappa = 0.99
    z = np.zeros(2)
    # Same single-strand + composition at two evidence levels: region 1 has 20× the counts of region 0.
    u_pos = np.array([20.0, 400.0])
    u_neg = np.array([20.0, 400.0])
    allow_pos = np.array([True, True])
    allow_neg = np.array([False, False])  # both single-strand
    mass = u_pos + u_neg
    d = _solve_regions_logodds_all(
        u_pos,
        u_neg,
        allow_pos,
        allow_neg,
        mass,
        z,
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_grid=60,
    )
    assert d.gdna_frac_var is not None
    # p̂=0.5 at κ=0.99 ⇒ the fragments look unstranded ⇒ the mean channel points at the gDNA mode f_g=1.
    # Under the count-zero-info variance freeze the count enters as PRECISION: more evidence sharpens that
    # signal, so the higher-count region resolves NEARER the mode with a lower Var(log f_g) — it is not pinned
    # count-independent (that was the pre-freeze f_g-dependent normalizer artifact).
    assert d.gdna_frac[0] > 0.85 and d.gdna_frac[1] > 0.85  # both gDNA-dominant
    assert d.gdna_frac[1] >= d.gdna_frac[0]  # more count ⇒ nearer the mode
    assert d.gdna_frac_var[1] < d.gdna_frac_var[0]  # more fragments ⇒ sharper
    # all per-component variances are present, finite, non-negative for active regions.
    for v in (d.gdna_frac_var, d.rna_pos_frac_var, d.rna_neg_frac_var):
        assert np.all(np.isfinite(v)) and np.all(v >= 0.0)
    # a no-fragment region is inactive ⇒ zero variance on every component.
    d0 = _solve_regions_logodds_all(
        np.array([0.0]),
        np.array([0.0]),
        np.array([True]),
        np.array([True]),
        np.array([0.0]),
        np.array([0.0]),
        kappa=kappa,
        od_g=0.2,
        od_r=0.1,
        n_grid=60,
    )
    assert d0.gdna_frac_var[0] == 0.0 and d0.rna_pos_frac_var[0] == 0.0


def _factor1_uniform_rho():
    """The factor-1 bedrock fixture: a UNIFORM-gDNA chain intergenic | AMBIG | intergenic, every object's
    count laid down as ``rho x its own placements`` (rho = 0.5). Returns the per-SLOT gDNA density after
    the sweep.

    ⭐ **One density per slot, not a (left, right) pair.** The predecessor returned two arrays because a
    boundary's two sides had different divisors; a 0-bp boundary has one. The invariant being tested is
    unchanged: lay down a uniform field and the solver must read it back.
    """
    rho = 0.5
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    region_eff = contained_eff_length(np.full(3, 1000.0), gdna_fl)  # [701, 701, 701]
    boundary_eff = float(
        crossing_eff_length(gdna_fl, np.full(1, UNBOUNDED_REACH), np.full(1, UNBOUNDED_REACH))[0]
    )
    region_count, boundary_count = rho * region_eff, rho * boundary_eff
    parts = make_chain_parts(
        [0, BIT_EXON_POS | BIT_EXON_NEG, 0],
        region_size_bp=1000.0,
        region_pos=region_count / 2,
        region_neg=region_count / 2,
        boundary_pos=boundary_count / 2,
        boundary_neg=boundary_count / 2,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    belief = init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.7, n_grid=40)
    final = region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        belief,
        parts.region_arrays,
        rna_sense_frac=0.7,
        n_grid=40,
    )
    # gDNA density = f_g x count / E_gdna (the formula the sweep inlines).
    count = np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)
    return final.f_g * count / np.asarray(parts.geometry.eff_gdna, float)


def test_gdna_sweep_factor1_intergenic_anchors():
    """The factor-1 bedrock, anchors: on a UNIFORM-gDNA chain the strand/signature-locked intergenic regions
    read back ρ EXACTLY — the measured-gDNA anchor invariant, which every phase must hold."""
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    interg = [0, 4]  # the chain is N E N E N, so the two intergenic regions are slots 0 and 4
    assert np.allclose(rho_g[interg], rho, atol=0.02)


def test_gdna_sweep_factor1_ambig_recovery():
    """The factor-1 bedrock, AMBIG region: a balanced AMBIG region between two ρ=0.5 anchors must read back ρ=0.5.

    ⚠ **This is the minimal reproduction of the PRIOR-FREE AMBIG WEAKNESS, and it passes with little room.**
    An AMBIG region has ``τ_own = 0`` (the strand likelihood constrains only the tilt), so it has no composition
    evidence of its own; all it gets is its neighbours' messages at their honest count precision, and ψ's
    uninformative Jeffreys reference deliberately holds it off the ``f_g = 1`` vertex until the data earn it.
    The shortfall is WEIGHT, not a wrong mode — it shrinks monotonically with depth, and the trained gDNA
    hyperprior is what replaces it (HANDOFF_6 §3).

    It was an ``xfail`` at ``ρ_g = 0.3914`` (21.7 % low). **Measured 2026-07-27: 0.45476, i.e. 9.0 % low —
    |err| 0.0452 against this test's 0.05 bound**, so the marker is gone and this is now a live guard. Expect
    it to be the first thing that trips on an AMBIG-facing change, and read a failure as "how much weight does
    a message deliver to a region with no evidence of its own", not as a tolerance nuisance. Do NOT attack it
    with more damping, and do not widen the bound without deriving what the residual SHOULD be."""
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    ambig = 2  # the AMBIG region slot
    assert np.allclose(rho_g[ambig], rho, atol=0.05)


def test_interior_anchor_is_immovable_and_produces_no_nan():
    """The `struct_lock` interior-anchor regression (HANDOFF_5 §6). A composition-CERTAIN region has
    ``Var(log f_c) = 0``, so any code path that forms a fusion weight as ``1/Var`` produces ``∞`` and cascades
    a nan through the whole chain. Pin both halves of the contract on the factor-1 chain, whose two intergenic
    REGIONs are exactly such anchors sitting INTERIOR to the chain (each has a live neighbour):

    1. **no nan anywhere** — beliefs and variances stay finite (``∞`` is the honest 'unsolved' state and is
       allowed on a variance; nan never is);
    2. **the anchor is IMMOVABLE** — it reads back the true ρ exactly despite receiving messages from an AMBIG
       neighbour that is itself wrong by 22%. Note what does the work: the anchor's own ``pg_own = n`` in the
       relay fuse, NOT the DL ``v_own = 0`` branch (which is inert at the combine because a struct_lock region is
       never `solvable`, so its ψ output is discarded)."""
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    assert not np.any(np.isnan(rho_g)), rho_g
    assert np.all(np.isfinite(rho_g)), rho_g
    assert np.allclose(rho_g[[0, 4]], rho, atol=1e-9)  # exact, not merely close


def test_gdna_emits_across_tss_tes_boundary():
    """Structural-gate regression: gDNA is genomically continuous, so the
    gene-boundary boundaries (TSS/TES) flanking a SINGLE-EXON gene must RELAY a gDNA message into it from the
    intergenic regions beyond — even though neither RNA strand is continuous across those boundaries. Before the
    fix the gDNA message was gated by RNA strand-continuity (`solvable`), so such a gene (both flanks
    intergenic) was a no-relay region, solving on its own local belief alone.

    Conversely, the intergenic flank is structurally RNA-free and emits ZERO RNA authority: the exon receives
    no +/− RNA message from it (a region's confidence about its OWN all-gDNA state grants no authority over a
    neighbour's RNA). The assertions lock the two halves of the three-term emission gate.
    """
    rho = 0.5
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    L, unb = 1000.0, np.full(1, UNBOUNDED_REACH)
    region_count = rho * float(contained_eff_length(np.full(1, L), gdna_fl)[0])
    boundary_count = rho * float(crossing_eff_length(gdna_fl, unb, unb)[0])
    parts = make_chain_parts(  # intergenic | exon+ | intergenic
        [0, BIT_EXON_POS, 0],
        region_size_bp=L,
        region_pos=region_count / 2,
        region_neg=region_count / 2,
        boundary_pos=boundary_count / 2,
        boundary_neg=boundary_count / 2,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    cap = {}
    final = region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.7, n_grid=40),
        parts.region_arrays,
        rna_sense_frac=0.7,
        n_grid=40,
        _capture=cap,
    )
    exon = 2  # the single-exon gene, flanked on both sides by TSS/TES boundaries (chain N E N E N)
    # THE FIX — the exon receives a gDNA relay across the boundary (incoming precision > 0). Pre-fix: 0 (no relay).
    assert cap["prec_g"][exon] > 0.0, (
        "single-exon gene got NO gDNA relay across the TSS/TES boundary"
    )
    # The intergenic flanks emit ZERO RNA authority: the exon receives no +/− RNA message from them.
    assert cap["prec_p"][exon] == 0.0 and cap["prec_n"][exon] == 0.0
    # State ⊥ messages: the intergenic regions stay locked all-gDNA (confident own-state, ignore all inputs).
    assert final.f_g[0] == 1.0 and final.f_g[4] == 1.0


def test_gdna_sweep_zero_gdna_pin_and_monotone():
    # ⚠ Was `xfail` as "pre-existing known-red" while σ²_transfer was identically 0 on this seedless chain:
    # the AMBIG region leant gDNA at ≈0.564 and the strand-solved introns were dragged up to ≈0.564 by the
    # directly-adjacent terminal G1 locks, whose messages were then UNDAMPED. The derived the-reframe-scale-variance σ²_transfer
    # (`Var(log r)` from `composition_logvar`, computed per boundary from counts and eff-lengths) damps them, and
    # measured 2026-07-27 all three regions are back under the 0.50 bound. Marker removed; live guard again.
    # A pure-RNA chain intron+ | AMBIG(in+|in−) | intron−. The AMBIG starts at the all-gDNA init f_g=1; the
    # global (driven to ~0 by the RNA introns) + the RNA-neighbour messages must pull the phantom gDNA down,
    # monotonically.
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    # sense-tilted RNA (κ=0.95): the + intron aligns genome+, the − intron genome−. The two boundaries carry
    # the same tilt as the regions they separate. ⭐ Two boundaries, not four: there are no terminal slots, so
    # the "directly-adjacent terminal G1 lock" the retired xfail blamed no longer exists at all.
    parts = make_chain_parts(
        [BIT_INTRON_POS, BIT_INTRON_POS | BIT_INTRON_NEG, BIT_INTRON_NEG],
        region_size_bp=2000.0,
        region_pos=[95.0, 50.0, 5.0],
        region_neg=[5.0, 50.0, 95.0],
        boundary_pos=[40.0, 2.0],
        boundary_neg=[2.0, 40.0],
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    chain, st, geom, region_arrays = (
        parts.chain,
        parts.statics,
        parts.geometry,
        parts.region_arrays,
    )
    belief = init_beliefs(chain, geom, st, rna_sense_frac=0.95, n_grid=40)
    assert belief.f_g[2] == 1.0  # AMBIG starts all-gDNA
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10000.0,  # library sample sizes so the stranded (κ=0.95) intron seeds fire (τ noise floor)
        n_gdna_obs=10000.0,
        n_grid=40,
    )
    # The AMBIG phantom is pulled DOWN from its all-gDNA init (1.0) toward RNA. This chain is the WORST
    # case for a balanced AMBIG region: it is an ARTIFICIAL all-RNA chain (intron+|AMBIG|intron−) with NO
    # intergenic structural seeds, so the gDNA prior has almost nothing to anchor a zero-gDNA baseline. The
    # AMBIG region's strand is balanced (both strands live), so the strand likelihood is DEGENERATE — a
    # balanced count is equally consistent with gDNA and with balanced ±RNA — and the neutral (f_g, τ)
    # reference measure parsimoniously leans a balanced count toward gDNA, deferring to the prior for the
    # level. With no seeds the (weak) floor + the intron RNA-imputation only pull it to ~0.44 here; on real
    # libraries the intergenic zero-count seeds make the prior decisive (the gdna_none capture-on benchmark
    # shows ~0 false gDNA). Still pulled well below the all-gDNA init and RNA-leaning.
    assert final.f_g[3] < 0.50
    # single-strand introns: the decisive strand wins and the floor DEFERS (a hyperprior cannot overrule a
    # region's own strand evidence) → they stay RNA-leaning, well below their all-gDNA init.
    #
    # ⭐ **S5.e removed the confound this comment used to be about.** The predecessor read ~0.44 rather
    # than ~0.22 because the chain's two TERMINAL boundary slots were G1-locked and emitted their
    # structural all-gDNA into the flanking introns — an artefact of an artificial chain, since an intron
    # running to the chain boundary with no intergenic flank cannot occur in a real annotation. **Those slots
    # no longer exist**: a reference with k regions owns k−1 boundaries, so there is nothing beyond the outer
    # regions to emit anything. The invariant the test protects is unchanged; the artefact is gone.
    assert final.f_g[0] < 0.50 and final.f_g[4] < 0.50


# density-Gaussian message form: two-sided pull + emergent deference --------


def test_density_message_two_sided_mode_not_vertex():
    """A LOG-fraction gDNA message (mode=log 0.2, strong prec) on a balanced AMBIG region (flat strand) pulls
    f_g TOWARD 0.2 — two-sided by construction (a Gaussian on log f_g, no boundary wall), not to the f_g=1
    vertex. The log-density log-odds solver's message form."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    z = np.zeros(1)
    # AMBIG region, balanced counts ⇒ the strand is flat (κ=0.5); only the message shapes f_g.
    d = _solve_regions_logodds_all(
        np.array([50.0]),
        np.array([50.0]),
        np.array([True]),
        np.array([True]),
        np.array([100.0]),
        z,
        kappa=0.5,
        od_g=0.0,
        od_r=0.0,
        n_grid=80,
        gdna_imp_mode=np.array([np.log(0.2)]),
        gdna_imp_prec=np.array([200.0]),
    )
    fg = float(d.gdna_frac[0])
    assert abs(fg - 0.2) < 0.05, fg


def test_density_message_defers_to_decisive_strand():
    """Emergent deference: a WEAK gDNA message (prec=3) trying to pull f_g→0.9 must lose to a decisive
    single-strand region's ~1000-fragment strand likelihood — f_g stays ≈0 (the honest precision blend means a
    weak message cannot override the data; no log-wall to force it off zero)."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    z = np.zeros(1)
    d = _solve_regions_logodds_all(
        np.array([1000.0]),
        np.array([5.0]),
        np.array([True]),
        np.array([False]),
        np.array([1005.0]),
        z,
        kappa=0.99,
        od_g=0.0,
        od_r=0.0,
        n_grid=80,
        gdna_imp_mode=np.array([np.log(0.9)]),
        gdna_imp_prec=np.array([3.0]),
    )
    fg = float(d.gdna_frac[0])
    assert fg < 0.1, fg


# --- mature absorption: the spliced mass "absorbs" the imputed RNA, leaving only NASCENT ---------------
# The RNA message src→dst is
#   ρ = src_nascent/E_r + SP[sf][src]/E_spl_src − SP[df][dst]/E_spl_dst.
# The dst-face term subtracts the mature a sj boundary measures, so a pure-mature exon imputes
# ≈0 nascent into the intron beyond it — no wholesale nascent hallucination.


#: slots of the ``_mature_exon_chain`` fixture. The chain is ``N E N E N E N E N`` — 9 slots, regions at
#: the even ones — so the exon under test is slot 4 and its two flanking introns are slots 2 and 6.
MX_EXON, MX_INTRONS = 4, [2, 6]


def _mature_exon_chain(*, spliced: bool, rho_g=0.5, rho_m=1.0, kappa=0.95, spl_scale=1.0):
    """``exon+ | intron+ | EXON+ | intron+ | exon+`` — a pure-MATURE expressed gene with NO nascent.

    ⭐ **Five regions, not three, and the extra two are load-bearing.** The predecessor put the mature
    flux on the two intron↔exon *boundaries* by hand, because the old accumulator attributed a splice to
    the region's boundary. A sj now states its own ``(src, dst)``, so it has to HAVE endpoints: the
    sj over intron ``n1`` runs ``n0 → n2`` and the one over ``n3`` runs ``n2 → n4``, and
    `build_region_geometry` places their flux on the boundaries they leave and enter. The exon under test
    (``n2``) ends up with mature flux on both its flanking boundaries — which is what the old fixture asserted
    by construction, now derived from the graph instead.

    Physically consistent: every exon's contained unspliced is balanced gDNA + sense (+) mature; the
    introns and every boundary carry balanced gDNA only. ⭐ **`boundary_spliced` is 0 everywhere, and that is a
    measured fact rather than a convenience** — mature RNA never crosses an exon↔intron boundary (0 of 1,146
    boundaries over 7 conditions). It skips the intron as a sj, never as
    a contiguous crossing.
    """
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    L = 2000.0
    unb = np.full(1, UNBOUNDED_REACH)
    Eg = float(contained_eff_length(np.full(1, L), gdna_fl)[0])  # contained gDNA placements
    Er = float(contained_eff_length(np.full(1, L), rna_fl)[0])  # contained RNA placements
    cross_g = float(crossing_eff_length(gdna_fl, unb, unb)[0])  # 299
    cross_r = float(crossing_eff_length(rna_fl, unb, unb)[0])  # 199
    g_half = rho_g * Eg / 2.0  # per-strand contained gDNA count (balanced)
    mat = rho_m * Er  # contained mature count (+ strand only)
    is_exon = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
    # ``spl_scale`` < 1 models a CAPTURE-DEPLETED sj: sj-spanning reads are only partially
    # captured, so the sj UNDER-reports the exon's true mature density ⇒ the boundary→exon mature
    # MEASUREMENT disagrees with the exon's own (confident) unspliced belief. Used by the silencing test.
    j_count = rho_m * cross_r * spl_scale
    sj = (
        [
            (0, 2, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, j_count),
            (2, 4, Strand.POS, UNBOUNDED_REACH, UNBOUNDED_REACH, j_count),
        ]
        if spliced
        else []
    )
    parts = make_chain_parts(
        [BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS],
        region_size_bp=L,
        region_pos=g_half + mat * is_exon,
        region_neg=g_half,
        boundary_pos=rho_g * cross_g / 2.0,
        boundary_neg=rho_g * cross_g / 2.0,
        boundary_spliced=0.0,
        sj=sj,
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    belief = init_beliefs(
        parts.chain, parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _sweep(args, kappa=0.95, n_rna_obs=10000.0, n_gdna_obs=10000.0):
    chain, st, geom, belief, ra = args
    cap = {}
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        ra,
        rna_sense_frac=kappa,
        # The τ strand seed needs the library sample sizes to size its overdispersion noise floor
        # ¼·(1/N + ω); a strongly-stranded fixture (κ=0.95) fires only when N is supplied (the default 0 ⇒
        # ∞ floor ⇒ gated). Large N here ⇒ floor≈0 ⇒ the (2κ−1)² strand seed fires at full strength.
        n_rna_obs=n_rna_obs,
        n_gdna_obs=n_gdna_obs,
        n_grid=60,
        _capture=cap,
    )
    return final, cap


def test_mature_no_nascent_hallucination_in_introns():
    """The owner's red boundary: a pure-mature expressed exon (nascent = 0) must NOT manufacture wholesale nascent
    in its flanking pure-gDNA introns; the introns stay gDNA (truth ``f_g = 1``).

    Was ``xfail`` at ``f_g ≈ 0.82`` — the exon's ~95 %-mature unspliced payload leaking in as nascent, because
    the RNA-total factor does not subtract mature. **Measured 2026-07-27: 0.9271**, comfortably past this
    test's 0.85 bound, so the marker is gone. The residual 0.073 is still the nascent-factory gap
    (``ρ_nascent = ρ_RNA − ρ_mature``, intron-baselined); tighten this bound when that lands rather than
    treating 0.85 as the target."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_introns = fin_m.f_g[MX_INTRONS]
    assert np.all(fg_introns > 0.85), fg_introns


# NOTE: `test_mature_absorption_lowers_nascent_message_into_sj` was RETIRED when the mature-crossing gate
# landed (Phase 4). It asserted the exon→TRAPS: measure-the-ceiling-first +RNA message FIRES (`app[b1] > 0`) so its absorption term could
# lower the imputed nascent; the gate now blocks that boundary entirely (the exon may not manufacture nascent into
# its intron-side sj), so the message no longer exists to absorb. Its replacement is
# `test_exon_does_not_manufacture_nascent_into_intron` (same fixture, same boundary, inverted assertion). The
# B→exon MEASUREMENT + absorption path it half-covered is still guarded by the two `test_mature_measurement_*`
# tests below, which the gate leaves untouched.


def test_mature_measurement_recovers_exon_rna():
    """The companion direction (unchanged B→exon MEASUREMENT): the same chain's expressed exon is
    correctly recovered as mostly RNA (its true f_g ≈ ρ_g·E_g/(ρ_g·E_g+ρ_m·E_r) ≈ 0.32), driven by the
    + strand tilt + the mature measurement — so the absorption does not starve the exon of its own RNA."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_exon = float(fin_m.f_g[MX_EXON])
    assert fg_exon < 0.45, fg_exon  # truth ≈0.32; comfortably RNA-dominated, not pinned to gDNA


def test_mature_measurement_disagreement_silenced():
    """BUG #2 regression: the mature MEASUREMENT message must be DISAGREEMENT-SILENCED like every other RNA
    message (the old exemption applied it at full COUNT precision). Under capture, sj-spanning reads are
    only partially captured, so the B→exon mature density UNDER-reports the exon's true RNA → the measurement
    DISAGREES with the exon's own confident belief. Un-silenced it dragged f_pos down → phantom gDNA by simplex
    complement (−gDNA flagship +0.04→+0.018). Here: a depleted sj genuinely lowers the message target,
    yet the exon's gDNA fraction stays unchanged vs a consistent sj — the disagreeing measurement was
    down-weighted, not applied whole (measured: the delivered precision collapses 188.9 → 1.24, 152×).

    ⚠ **This was a `strict=True` xfail on an OPEN ITEM, and the open item is RESOLVED.** The defect was a
    real cross-component coupling: `_boundary_spliced_mass_increment` folded the SPLICED density into the
    mature-inclusive boundary projection used for σ²_transfer on exon↔boundary boundaries, so depleting the
    spliced channel moved σ²_transfer and hence the attenuation of the **gDNA** relay — components that
    should not touch. Assertion (2) failed at |Δf_g| = 0.0612 against its 0.05 bound. σ²_transfer is now the
    derived the-reframe-scale-variance `Var(log r)` (`composition_logvar`, per boundary from counts and eff-lengths) and that projection
    is out of the path: **measured 2026-07-27, |Δf_g| = 0.000000 — exactly zero, not merely inside the
    bound.**

    Neither empirical bound was widened (the retired marker warned against exactly that). What changed is the
    STIMULUS: the depletion is now 10× (`spl_scale=0.1`) rather than 4×, because at 4× the message target now
    moves 0.2635 nats against assertion (1)'s 0.3 — a precondition on "is the disagreement big enough to be
    worth silencing", not a guarantee this test protects. A 10× depletion moves it 0.6277 and is the harder
    test."""
    ex = MX_EXON
    fin_ok, cap_ok = _sweep(_mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=1.0))
    fin_lo, cap_lo = _sweep(_mature_exon_chain(spliced=True, rho_m=4.0, spl_scale=0.1))
    # (1) the depleted sj really did lower the +RNA message target into the exon (a genuine disagreement)…
    assert cap_lo["mode_p"][ex] < cap_ok["mode_p"][ex] - 0.3, (
        cap_lo["mode_p"][ex],
        cap_ok["mode_p"][ex],
    )
    # (2) …yet the exon's gDNA fraction barely moves (silenced — pre-fix the low measurement inflated f_g).
    assert abs(float(fin_lo.f_g[ex]) - float(fin_ok.f_g[ex])) < 0.05, (
        fin_lo.f_g[ex],
        fin_ok.f_g[ex],
    )
    # (3) and the exon stays RNA-dominated, not pulled toward phantom gDNA.
    assert float(fin_lo.f_g[ex]) < 0.45, fin_lo.f_g[ex]


def test_tau_gag_fix_spliced_sj_emits_when_unstranded():
    """τ-GAG REGRESSION ( §Phase B, 2026-07-21). On UNSTRANDED data
    (κ=½ ⇒ the strand Fisher info ``I_strand`` is identically 0), a splice-junction boundary still carries
    motif-stranded spliced (mature-RNA) fragments — a DIRECT measurement, independent of strand. That
    measurement MUST reach the exon. The bug: the τ-evidence emission gate (keyed on ``I_strand``+``I_struct``
    only, NOT the spliced count) silenced it, and the spliced-precision credit — which lives *inside* the gated
    block — never fired (52% of sj). The fix opens RNA emission on spliced presence while keeping the
    deconvolution PREDICTION τ-gated.

    Pins both halves: (1) a spliced sj DELIVERS a +RNA message to its exon even unstranded; (2) the same
    chain with the spliced REMOVED delivers zero +RNA authority (a vacuous unstranded region manufactures no
    phantom RNA — the deconvolution stays gated). This exact pair fails on the pre-fix gated code."""
    ex = MX_EXON
    fin_spl, cap_spl = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # (1) THE FIX: the spliced (mature) MEASUREMENT reaches the exon with the strand silent (κ=½).
    assert cap_spl["prec_p"][ex] > 0.0, cap_spl["prec_p"][ex]
    # (2) the vacuous control (no spliced, no strand): ZERO +RNA authority — no phantom manufactured.
    assert cap_no["prec_p"][ex] == 0.0, cap_no["prec_p"][ex]
    # (3) the real mature measurement moves the exon TOWARD RNA — never toward phantom gDNA.
    assert float(fin_spl.f_g[ex]) < float(fin_no.f_g[ex]), (fin_spl.f_g[ex], fin_no.f_g[ex])


def test_tau_gag_fix_deconvolution_prediction_stays_gated():
    """The safety half of the τ-gag fix: unblocking the spliced MEASUREMENT must NOT unblock the deconvolution
    PREDICTION on a vacuous source (that is the phantom the τ-precision exists to kill). On the unstranded
    no-spliced chain, the boundary→exon +RNA precision is exactly 0 (asserted above); here we pin that even the
    gDNA message a vacuous boundary sends carries no manufactured composition confidence beyond the honest
    structural/count evidence — the exon's solved f_g stays near the uninformative reference, not driven to a
    confident vertex by a phantom."""
    ex = MX_EXON
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # No spliced + no strand ⇒ the exon has no composition evidence of its own; it must not be pinned to a
    # confident vertex by the messages (the phantom would drive it to ~1). It stays mid-range (reference-led).
    assert 0.2 < float(fin_no.f_g[ex]) < 0.8, fin_no.f_g[ex]


def test_strand_overdispersion_prior_default_is_near_binomial():
    """BUG #1 regression: the shipped default strand-overdispersion prior must be the NEAR-BINOMIAL null
    (α=β=14 ⇒ od₀≈0.034), NOT the old over-conservative 0.143 (α=β=3) that widened the gDNA Beta-Binomial
    and erased its specificity at its own mean ½.

    ⭐ 2026-07-28: the four ``CalibrationConfig`` prior fields were collapsed into the two asserted module
    constants next to the estimator, so this now asserts on those. The assertion itself is UNCHANGED and
    still binds — it is the reason ``od₀ = od_max/2 = 0.1`` (derived independently as the max-entropy mean
    of the ceiling-bounded prior) was REJECTED: at a = 4.5 it fails this test, and it was measured to
    collapse a balanced pure-gDNA region's strand log-evidence 305.4 → 113.9 nats."""
    from rigel.calibration.gdna_strand import (
        _CEIL_ALPHA_BETA,
        _MAX_OVERDISPERSION,
        _PRIOR_ALPHA_BETA,
        _PRIOR_OVERDISPERSION,
        overdispersion_for_beta,
    )

    assert _PRIOR_ALPHA_BETA == 14.0
    assert _PRIOR_OVERDISPERSION < 0.05  # the near-binomial null — BUG #1's fix
    assert overdispersion_for_beta(3.0) > 0.14  # the old default was ~4× more overdispersed
    # the ceiling is the admissibility clamp, and the prior must sit strictly inside it
    assert _CEIL_ALPHA_BETA == 2.0
    assert _MAX_OVERDISPERSION == pytest.approx(0.2)
    assert 0.0 < _PRIOR_OVERDISPERSION < _MAX_OVERDISPERSION


def test_strand_overdispersion_prior_weight_is_derived_not_asserted():
    """The prior's weight ``W`` must be DERIVED from the two asserted constants, in the data's own
    information units — it used to be an asserted ``30`` in *seed-region* units, which is the wrong currency
    for a second moment (a 1-fragment seed carries no information about a correlation between fragments).

    ``W = 1/τ²`` with ``τ²`` the variance of the maximum-entropy distribution on ``[0, od_max]`` with mean
    ``od₀`` — the least-committal prior given exactly what we assert."""
    from rigel.calibration.gdna_strand import (
        _MAX_OVERDISPERSION,
        _PRIOR_INFORMATION,
        _PRIOR_OVERDISPERSION,
    )

    # bracketed by the two distribution-shape extremes on the same support and mean
    two_point = _PRIOR_OVERDISPERSION * (_MAX_OVERDISPERSION - _PRIOR_OVERDISPERSION)
    assert 1.0 / two_point < _PRIOR_INFORMATION < 1.0 / _PRIOR_OVERDISPERSION**2 * 1.5
    assert _PRIOR_INFORMATION == pytest.approx(909.1, rel=1e-3)


def test_null_information_reduces_to_pair_count_at_symmetric_mean():
    """``I = 1/Var(od_mom)|₀`` must equal the PAIR COUNT ``Σ n(n−1)/2`` exactly at μ = ½ (the gDNA case),
    and must NOT be substituted by the pair count away from it (measured ``I/pairs`` = 0.05–0.14 at the RNA
    fit's κ, i.e. pairs overstate RNA information 7–20×)."""
    import numpy as np

    from rigel.calibration.gdna_strand import _null_information

    n = np.array([1.0, 2.0, 2.0, 10.0, 100.0])
    pairs = float((n * (n - 1.0) / 2.0).sum())
    assert _null_information(n, 0.25) == pytest.approx(pairs, rel=1e-12)
    # a singleton contributes nothing, so dropping it changes nothing
    assert _null_information(n[1:], 0.25) == pytest.approx(pairs, rel=1e-12)
    # away from ½ the information is strictly LESS than the pair count
    assert _null_information(n, 0.01 * 0.99) < pairs


def test_pure_gdna_region_confident_at_near_binomial_od():
    """BUG #1 mechanism (unit): a pure-gDNA single-strand region has EXACT 50/50 per-strand counts, which the
    strand mixture (gDNA mean ½, RNA mean κ≠½) must read as gDNA — f_g≈1. At the near-binomial od (the fixed
    default) it does; at the old inflated od=0.143 the widened gDNA BB loses specificity at ½ and the region is
    dragged toward the RNA/gDNA boundary (f_g well below 1). A pure-RNA control (+frac=κ) stays f_g≈0 at both."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    # κ=0.7 (intermediate strand): gDNA mean ½ is near enough to the RNA mean that the gDNA BB width matters —
    # exactly where the inflated prior does its damage (and where the toy battery regressed pre-fix).
    def solve(u_pos, u_neg, od):
        z = np.zeros(1)
        n = float(u_pos + u_neg)
        return float(
            _solve_regions_logodds_all(
                np.array([float(u_pos)]),
                np.array([float(u_neg)]),
                np.array([True]),
                np.array([False]),
                np.array([n]),
                z,
                kappa=0.7,
                od_g=od,
                od_r=od,
                n_grid=80,
            ).gdna_frac[0]
        )

    # pure gDNA (exact 50/50, truth f_g=1): near-binomial od → confidently gDNA; inflating od monotonically
    # under-calls it. The harm is REAL but far smaller than it used to look: this assertion previously needed
    # only od=0.143 to force a >0.15 collapse, because ψ then also carried the improper `+0.5·λ` ramp, which
    # AMPLIFIED od harm by fighting the strand near the vertex. With ψ bare the strand speaks cleanly and the
    # solver is materially more od-robust (0.9945 → 0.9683 at od=0.143, vs a collapse below 0.844 before) —
    # i.e. the old numbers were further from the truth, not closer. Assert the physics, not the artifact.
    fg = [solve(500, 500, od) for od in (0.034, 0.143, 0.4)]
    assert fg[0] > 0.8, fg
    assert fg[0] > fg[1] > fg[2], fg  # monotone: inflating od always degrades the gDNA call
    assert fg[2] < fg[0] - 0.15, fg  # and at a materially inflated od the damage is large
    # pure RNA (+frac = κ = 0.7): near-binomial stays RNA-dominated; the inflated prior's symmetric harm is
    # MORE false gDNA on RNA too (it pulls every region toward ½). (At this intermediate κ the gDNA/RNA means
    # are close, so a small residual f_g is inherent — the point is near-binomial is cleaner.)
    rna_near = solve(700, 300, 0.034)
    rna_infl = solve(700, 300, 0.143)
    assert rna_near < 0.25, rna_near
    assert rna_infl > rna_near, (rna_infl, rna_near)


# ---------------------------------------------------------------------------
# RNA-message routing after the mature-crossing gate was DISMANTLED
#
# Only the STRUCTURAL per-strand `free_s` continuity gate remains: each RNA strand's density flows wherever that
# strand is continuous on BOTH endpoints (intron↔exon in either direction, intron→boundary, boundary→exon), and
# gDNA flows genomically. The asymmetric `send_s = mrna_active_s[dst] or not mrna_active_s[src]` gate that used
# to silence exon→intron RNA is GONE. On the `_mature_exon_chain` fixture (intron+ | exon+ | intron+, chain
# B0 R0 B1 R1 B2 R2 B3) EVERY continuous-strand boundary now fires — including the formerly-silenced exon R1→B1
# (backward) and exon R1→B2 (forward). That re-opens the mature leak into the introns; the honest σ²_transfer
# precision + the nascent factory (ρ_nascent = ρ_RNA − ρ_mature) are what will counter it (see §6 of the doc).
# The `mrna_active_*` mask itself stays computed in the statics (the nascent factory will consume it).
# ---------------------------------------------------------------------------

# chain region ids for the mature-exon fixture (intergenic|intron R0|B1|exon R1|B2|intron R2|...):
_R1_EXON = 3  # the expressed exon R1
_B1 = 2  # intron→exon sj; its right neighbour (backward src) is R1
_B2 = 4  # exon→intron sj; its left neighbour (forward src) is R1


def test_intron_relays_nascent_into_exon_both_directions():
    """The structural-continuity guard: +RNA (nascent) must keep flowing along the +strand-continuous chain in
    BOTH directions (intron R0→B1→exon R1 forward; intron R2→B2→exon R1 backward). The unified relay fuses each
    region's own belief with the transported neighbour, so a live +RNA precision (`fwd_pp`/`bwd_pp` > 0) at these
    regions is the relay firing. Guards against a regression that would delete the intron→exon nascent relay."""
    _, cap = _sweep(_mature_exon_chain(spliced=True))
    uni = cap["_uni_static"]
    fpp, bpp = uni["fwd_pp"], uni["bwd_pp"]  # forward / backward +RNA precision after the relay
    # +strand-continuous chain ⇒ the fused +RNA precision is live at the sj and the exon, both directions
    assert fpp[_B1] > 0.0, fpp[
        _B1
    ]  # forward relay reaches the intron→exon sj TRAPS: measure-the-ceiling-first
    assert bpp[_B2] > 0.0, bpp[
        _B2
    ]  # backward relay reaches the exon→intron sj TRAPS: score-against-truth
    assert fpp[_R1_EXON] > 0.0 and bpp[_R1_EXON] > 0.0  # the exon receives +RNA from both flanks


def test_mrna_active_matches_same_strand_exon_rule():
    """The `mrna_active_strands` mature-presence mask (no longer wired into the emission gate, but kept in the
    statics for the coming nascent factory `ρ_nascent = ρ_RNA − ρ_mature`) is exactly the user's rule: mature is
    present on strand s across a boundary iff the SAME-STRANDED exon bit is set on BOTH flanks. Intron bits never
    qualify; `EX+EX- | EX+EX-` passes on BOTH strands. Enumerate all 16×16 signature pairs (a boundary's two
    flanks) and check `mrna_active_strands` against that predicate, plus the subsumption `mrna_active_s ⇒
    nrna_active_s` (mature ⇒ nascent). Pure, no sweep."""
    sigs = np.arange(N_SIGNATURES, dtype=np.int64)
    for sl in sigs:
        for sr in sigs:
            mrp_l, mrn_l = mrna_active_strands(np.array([sl]))
            mrp_r, mrn_r = mrna_active_strands(np.array([sr]))
            # a boundary's per-strand mature-crossing = AND of the two flanks' own exon bits
            mrp = bool(mrp_l[0] and mrp_r[0])
            mrn = bool(mrn_l[0] and mrn_r[0])
            exp_pos = bool((sl & BIT_EXON_POS) and (sr & BIT_EXON_POS))
            exp_neg = bool((sl & BIT_EXON_NEG) and (sr & BIT_EXON_NEG))
            assert mrp == exp_pos, (sl, sr, mrp, exp_pos)
            assert mrn == exp_neg, (sl, sr, mrn, exp_neg)
            # subsumption: mature-active ⇒ nascent-active (an exon carries nascent too), per strand
            nrp_l, nrn_l = nrna_active_strands(np.array([sl]))
            nrp_r, nrn_r = nrna_active_strands(np.array([sr]))
            nrp = bool(nrp_l[0] and nrp_r[0])
            nrn = bool(nrn_l[0] and nrn_r[0])
            assert not mrp or nrp, (sl, sr)  # mrp ⇒ nrp
            assert not mrn or nrn, (sl, sr)

    # the user's headline case: overlapping opposite-strand exons on both flanks ⇒ mature passes on BOTH strands
    both = BIT_EXON_POS | BIT_EXON_NEG
    mrp, mrn = mrna_active_strands(np.array([both]))
    mrp2, mrn2 = mrna_active_strands(np.array([both]))
    assert bool(mrp[0] and mrp2[0]) and bool(mrn[0] and mrn2[0])
    # and an exon+intron mixed flank does NOT block the exon's own strand (+ passes; − is intron→intron ⇒ no)
    mixed = BIT_EXON_POS | BIT_INTRON_NEG  # exon on +, intron on −
    mp_l, mn_l = mrna_active_strands(np.array([mixed]))
    mp_r, mn_r = mrna_active_strands(np.array([mixed]))
    assert bool(mp_l[0] and mp_r[0])  # + strand: exon|exon ⇒ mature passes
    assert not bool(mn_l[0] and mn_r[0])  # − strand: intron|intron ⇒ no mature


def test_sweep_finite_over_extreme_configs():
    """F: no nan/inf reaches the fold. The real region_sweep over spliced/±, stranded/unstranded, and extreme
    gDNA/mature densities (pure-gDNA, pure-RNA, empty, tiny, huge) — every final fraction is finite & in range,
    every variance is ≥0 (∞ = the honest 'unsolved' state is allowed; nan is not), every emitted message
    mode/precision is finite."""
    for spliced in (True, False):
        for kappa in (0.5, 0.95):
            for rho_g, rho_m in [(0.5, 1.0), (0.0, 1.0), (2.0, 0.0), (1e-6, 1e-6), (1e4, 1e4)]:
                cfg = dict(spliced=spliced, kappa=kappa, rho_g=rho_g, rho_m=rho_m)
                final, cap = _sweep(_mature_exon_chain(**cfg), kappa=kappa)
                for nm in ("f_g", "f_pos", "f_neg"):
                    v = np.asarray(getattr(final, nm))
                    assert np.all(np.isfinite(v)), (cfg, nm, v)
                    assert np.all(v >= -1e-9) and np.all(v <= 1.0 + 1e-9), (cfg, nm, v)
                for nm in ("var_gdna", "var_pos", "var_neg"):
                    v = np.asarray(getattr(final, nm))
                    assert not np.any(np.isnan(v)) and np.all(v >= -1e-12), (
                        cfg,
                        nm,
                        v,
                    )  # ∞ ok, nan not
                # every unified message mode/precision (relay + combine) is finite — nothing nan reaches the ψ solve
                for nm in ("mode_g", "prec_g", "mode_p", "prec_p", "mode_n", "prec_n"):
                    assert np.all(np.isfinite(np.asarray(cap[nm]))), (cfg, nm)
                uni = cap["_uni_static"]
                for nm in ("fwd_g", "fwd_p", "fwd_n", "fwd_pg", "fwd_pp", "fwd_pn"):
                    assert np.all(np.isfinite(np.asarray(uni[nm]))), (cfg, nm)


def test_region_sweep_deterministic():
    """H: pass-0 must be bit-reproducible. The forward-backward BP sweep is sequential Python (no parallel
    reduction), so the same input must give a bit-identical belief AND identical emitted messages run-to-run —
    a prerequisite for any confidence claim about the solver. Uses the unstranded (κ=½), fully message-driven
    case, the one most sensitive to any ordering nondeterminism."""
    a, capa = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    b, capb = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    for nm in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg"):
        x, y = np.asarray(getattr(a, nm)), np.asarray(getattr(b, nm))
        assert np.array_equal(x, y, equal_nan=True), (nm, x, y)  # BIT-identical (not just close)
    # every unified imputation factor (the relay/combine output feeding the ψ solve) is bit-identical run-to-run
    for nm in ("mode_g", "prec_g", "mode_p", "prec_p", "mode_n", "prec_n"):
        x, y = np.asarray(capa[nm]), np.asarray(capb[nm])
        assert np.array_equal(x, y, equal_nan=True), (nm, x, y)


def test_float32_log_is_monotone_so_the_ambig_cube_may_hoist_it():
    """`_solve_ambig_logodds` computes `log f_pos` as `max(log f_grid, log floor)` rather than
    `log(max(f_grid, floor))` — the log on the (K,K_t) GRID instead of the (m,K,K_t) cube, ~140x fewer
    transcendentals for the same bits.

    That rewrite is EXACT iff numpy's float32 `log` is monotone on [0,1], the whole domain both arguments
    live in (a fraction, and `1/(n+1)`). It was verified exhaustively over all 1,065,353,217 float32 values
    there; this pins it against a numpy/platform change with a dense consecutive-value sweep per exponent
    band, plus the identity itself over the shape the solver actually forms."""
    with np.errstate(divide="ignore"):
        for e in range(-40, 1):  # one dense run of consecutive float32s per exponent band in [0,1]
            lo = np.float32(2.0**e).view(np.uint32)
            x = np.arange(lo, lo + 200_000, dtype=np.uint32).view(np.float32)
            assert (np.diff(np.log(x)) >= 0.0).all(), e

        rng = np.random.default_rng(4)
        grid = rng.random((60, 60)).astype(np.float32)  # the (K,K_t) fraction grid
        grid[rng.random(grid.shape) < 0.05] = 0.0  # the tau = +-1 boundaries: f_s exactly 0
        floor = (1.0 / (1.0 + np.exp(rng.uniform(0.0, 14.0, 500)))).astype(np.float32)[
            :, None, None
        ]
        direct = np.log(np.maximum(grid[None, :, :], floor))
        hoisted = np.maximum(np.log(grid)[None, :, :], np.log(floor))
    assert np.array_equal(direct.view(np.int32), hoisted.view(np.int32))
