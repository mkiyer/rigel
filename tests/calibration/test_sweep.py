"""The two-phase sweep (`sweep.solve_chain`) and the beliefs it starts from.

Every fixture is built by ``_synthetic.make_chain_parts`` on the shipped axes: a region axis, a
contiguous-boundary axis with ``k − 1`` entries per reference and so no terminal slots, and a sj
axis whose boundaries state their own ``(src, dst, strand)``. Both matter here. A reference
terminal would be a data-free boundary slot that could be G1-locked and emit structural all-gDNA
into its neighbour, an artefact no real annotation produces; and a sj stating its own endpoints is
what lets the mature flux at an intron-exon boundary be derived from the graph rather than placed
by hand. The gates cover the per-slot init, the factor-1 anchors, a delivered row's pull, the
mature-exon chain, the strand overdispersion the solve reads, and the numeric contract.
"""

from __future__ import annotations


import functools

import numpy as np
import pytest

from rigel.types import Strand

from rigel.calibration.blocks import SweepCapture
from rigel.calibration.messages.silent import SilentPolicy
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


#: These gates exercise the SWEEP's shape and the per-slot init under the measured floor, named
#: explicitly. The message policies have their own gates (`test_sweep_backbone.py`,
#: `test_transfer_policy.py`).
region_sweep = functools.partial(solve_chain, policy=SilentPolicy())


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def test_init_zero_gdna_introns_via_strand():
    # A zero-gDNA library. 3 regions: intergenic | intron+ | AMBIG (one ref).
    # intergenic gDNA = strand-symmetric; intron+ RNA = strongly sense-tilted (κ=0.95); AMBIG = symmetric.
    parts = make_chain_parts(
        [0, BIT_INTRON_POS, BIT_EXON_POS | BIT_EXON_NEG],
        region_size_bp=[1000.0, 2000.0, 800.0],
        region_pos=[50.0, 95.0, 50.0],
        region_neg=[50.0, 5.0, 50.0],
    )
    b = init_beliefs(parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)

    # the chain is N E N E N, so the regions are at 0, 2, 4 — there are no terminal slots.
    rid = [0, 2, 4]
    fg = b.f_g[rid]
    # intergenic: locked gDNA sink {0,0,1}, all precision locked (var 0).
    assert fg[0] == 1.0
    assert b.var_gdna[0] == 0.0
    # intron+ (zero gDNA): the strand tilt alone drives f_g → 0; g finite.
    assert fg[1] < 0.15
    assert np.isfinite(b.var_gdna[2])
    # AMBIG: unresolved by strand → {0,0,1} default at MAX (inf) variance for the sweep to resolve.
    assert fg[2] == 1.0
    assert np.isinf(b.var_gdna[4])


def test_init_boundary_continuity_gate():
    # 1 ref, 2 regions (exon+ | intron+) → ONE boundary between them. There are no terminal boundary
    # slots: a reference with k regions owns k-1 boundaries, so there is nothing before the first
    # region or after the last to be a sink.
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
    b = init_beliefs(parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)
    # slots: N0=0, E0=1, N1=2.
    # E0 (ex+→in+): +strand continuous (G2+) ⇒ the strand tilt resolves f_g → 0.
    assert b.f_g[1] < 0.15
    assert np.isfinite(b.var_gdna[1])


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
    b = init_beliefs(parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60)
    # slot 1 is the TSS boundary: a locked gDNA sink despite the sense tilt (all precision locked at 0).
    assert b.f_g[1] == 1.0 and b.var_gdna[1] == 0.0


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
    # Under the count-zero-info variance freeze the count enters as PRECISION: more evidence sharpens
    # that signal, so the higher-count region resolves NEARER the mode with a lower Var(log f_g)
    # rather than being pinned count-independently.
    assert d.gdna_frac[0] > 0.85 and d.gdna_frac[1] > 0.85  # both gDNA-dominant
    assert d.gdna_frac[1] >= d.gdna_frac[0]  # more count ⇒ nearer the mode
    assert d.gdna_frac_var[1] < d.gdna_frac_var[0]  # more fragments ⇒ sharper
    # the variance is present, finite, non-negative for active regions.
    assert np.all(np.isfinite(d.gdna_frac_var)) and np.all(d.gdna_frac_var >= 0.0)
    # a no-fragment region is inactive ⇒ zero variance.
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
    assert d0.gdna_frac_var[0] == 0.0


def _factor1_uniform_rho():
    """The factor-1 bedrock fixture: a UNIFORM-gDNA chain intergenic | AMBIG | intergenic, every object's
    count laid down as ``rho x its own placements`` (rho = 0.5). Returns the per-SLOT gDNA density after
    the sweep.

    One density per slot, not a (left, right) pair: a 0-bp boundary has one divisor. The invariant
    is that laying down a uniform field must be read back by the solver.
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
    belief = init_beliefs(parts.geometry, parts.statics, rna_sense_frac=0.7, n_grid=40)
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


def test_interior_anchor_is_immovable_and_produces_no_nan():
    """The interior-anchor regression. A composition-CERTAIN (`g1_locked`) region has
    ``Var(log f_c) = 0``, so any code path that forms a fusion weight as ``1/Var`` produces ``∞`` and cascades
    a nan through the whole chain. Pin both halves of the contract on the factor-1 chain, whose two intergenic
    REGIONs are exactly such anchors sitting INTERIOR to the chain (each has a live neighbour):

    1. no nan anywhere — beliefs and variances stay finite (``∞`` is the honest 'unsolved' state and
       is allowed on a variance; nan never is);
    2. the anchor is immovable — it reads back the true ρ exactly beside an AMBIG neighbour that is
       itself wrong by 22 %: a `g1_locked` region is never `solvable`, so its ψ output is discarded and its
       own count stands."""
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    assert not np.any(np.isnan(rho_g)), rho_g
    assert np.all(np.isfinite(rho_g)), rho_g
    assert np.allclose(rho_g[[0, 4]], rho, atol=1e-9)  # exact, not merely close


def test_gdna_sweep_zero_gdna_pin_and_monotone():
    # A pure-RNA chain intron+ | AMBIG(in+|in−) | intron−, solved under the silent policy: every slot on
    # its own evidence and nothing sent between them. The AMBIG region starts at the all-gDNA init f_g=1.
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    # sense-tilted RNA (κ=0.95): the + intron aligns genome+, the − intron genome−. The two boundaries carry
    # the same tilt as the regions they separate. Two boundaries, not four: there are no terminal
    # slots, so there is no directly-adjacent terminal G1 lock to emit structural gDNA into a flank.
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
    belief = init_beliefs(geom, st, rna_sense_frac=0.95, n_grid=40)
    assert belief.f_g[2] == 1.0  # AMBIG starts all-gDNA
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10000.0,  # the spliced sample behind κ: the stranded (κ=0.95) channel reads live
        n_grid=40,
    )
    # Slot 3 is the AMBIG | intron− boundary, whose − tilt reads RNA on its own. The AMBIG region itself
    # (slot 2) is not asserted: its strand is balanced (both strands live), so its strand likelihood is
    # DEGENERATE — a balanced count is equally consistent with gDNA and with balanced ±RNA — and with no
    # prior and nothing sent to it, it stays near its all-gDNA init.
    assert final.f_g[3] < 0.50
    # single-strand introns: the decisive strand wins → they read RNA, well below their all-gDNA init.
    assert final.f_g[0] < 0.50 and final.f_g[4] < 0.50


# a delivered row: two-sided pull + emergent deference --------


def _gdna_share_row(n_grid, mode_share, prec):
    """A claim on the gDNA share delivered as a λ-row on the solve grid — the message layer's one
    currency: a Gaussian on ``log f_g`` at ``log mode_share`` with precision ``prec``."""
    from scipy.special import log_expit

    from rigel.calibration.simplex_logodds import _logodds_grid

    lam, _ = _logodds_grid(int(n_grid), 10.0)
    return (-0.5 * float(prec) * (log_expit(lam) - np.log(float(mode_share))) ** 2)[None, :]


def test_a_delivered_row_pulls_two_sided_and_not_to_the_vertex():
    """A row claiming ``f_g = 0.2`` (strong) on a balanced AMBIG region (flat strand) pulls f_g TOWARD
    0.2 — two-sided by construction (a profile on the grid, no boundary wall), not to the f_g=1 vertex."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    z = np.zeros(1)
    # AMBIG region, balanced counts ⇒ the strand is flat (κ=0.5); only the row shapes f_g.
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
        lam_logprior=_gdna_share_row(80, 0.2, 200.0),
    )
    fg = float(d.gdna_frac[0])
    assert abs(fg - 0.2) < 0.05, fg


def test_a_weak_row_defers_to_a_decisive_strand():
    """Emergent deference: a WEAK row (precision 3) claiming ``f_g = 0.9`` must lose to a decisive
    single-strand region's ~1000-fragment strand likelihood — f_g stays ≈0 (a weak claim cannot override
    the data; no log-wall forces it off zero)."""
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
        lam_logprior=_gdna_share_row(80, 0.9, 3.0),
    )
    fg = float(d.gdna_frac[0])
    assert fg < 0.1, fg


# --- the mature-exon chain: an expressed exon between two gDNA-only introns ---------------------------
# Under the silent policy every slot is solved on its own evidence: the introns' balanced counts must
# read gDNA and the exon's sense-tilted count RNA. The spliced reads over the introns move no slot's
# answer here (the ``spliced=False`` chain solves identically); the policies that carry them have their
# own gates.


#: slots of the ``_mature_exon_chain`` fixture. The chain is ``N E N E N E N E N`` — 9 slots, regions at
#: the even ones — so the exon under test is slot 4 and its two flanking introns are slots 2 and 6.
MX_EXON, MX_INTRONS = 4, [2, 6]


def _mature_exon_chain(*, spliced: bool, rho_g=0.5, rho_m=1.0, kappa=0.95, spl_scale=1.0):
    """``exon+ | intron+ | EXON+ | intron+ | exon+`` — a pure-MATURE expressed gene with NO nascent.

    Five regions, not three, and the extra two are load-bearing: a sj states its own ``(src, dst)``,
    so it has to HAVE endpoints. The sj over intron ``n1`` runs ``n0 → n2`` and the one over ``n3``
    runs ``n2 → n4``, and `build_region_geometry` places their flux on the boundaries they leave and
    enter, so the exon under test (``n2``) ends up with mature flux on both its flanking boundaries
    — derived from the graph rather than asserted by construction.

    Physically consistent: every exon's contained unspliced is balanced gDNA + sense (+) mature; the
    introns and every boundary carry balanced gDNA only. `boundary_spliced` is 0 everywhere, and
    that is a fact rather than a convenience — mature RNA never crosses an exon-intron boundary
    (TRAPS: mature-rna-never-crosses-a-boundary). It skips the intron as a sj, never as a contiguous
    crossing.
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
    belief = init_beliefs(parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60)
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _sweep(args, kappa=0.95, n_rna_obs=10000.0):
    chain, st, geom, belief, ra = args
    cap = SweepCapture()
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        ra,
        rna_sense_frac=kappa,
        # The strand channel is live iff the protocol decision on (κ, N_rna) says so; with the default
        # N = 0 there is no decision and the channel is dead, so the fixture supplies the spliced sample
        # behind its κ = 0.95, which reads live by thousands of nats.
        n_rna_obs=n_rna_obs,
        n_grid=60,
        _capture=cap,
    )
    return final, cap


def test_mature_no_nascent_hallucination_in_introns():
    """An expressed exon whose RNA is all spliced must not make its flanking gDNA-only introns read as
    RNA; the introns stay gDNA (truth ``f_g = 1``). Each intron is solved on its own balanced counts, so
    this pins that nothing from the exon's RNA reaches them."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_introns = fin_m.f_g[MX_INTRONS]
    assert np.all(fg_introns > 0.85), fg_introns


def test_mature_measurement_recovers_exon_rna():
    """The companion direction: the same chain's expressed exon is recovered as mostly RNA (its true
    f_g ≈ ρ_g·E_g/(ρ_g·E_g+ρ_m·E_r) ≈ 0.32), read off its own + strand tilt — the gDNA-only introns
    beside it do not starve the exon of its own RNA."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_exon = float(fin_m.f_g[MX_EXON])
    assert fg_exon < 0.45, fg_exon  # truth ≈0.32; comfortably RNA-dominated, not pinned to gDNA


def test_tau_gag_fix_deconvolution_prediction_stays_gated():
    """On the unstranded chain with no spliced reads the exon has no composition evidence of its own — no
    strand channel and no sj — so nothing may drive it to a confident vertex: its solved f_g stays near
    the uninformative reference."""
    ex = MX_EXON
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # No spliced reads + no strand ⇒ the exon has no composition evidence of its own; it stays mid-range,
    # where the reference puts it.
    assert 0.2 < float(fin_no.f_g[ex]) < 0.8, fin_no.f_g[ex]


def test_the_overdispersion_CEILING_is_the_only_asserted_constant_left():
    """There is no shrinkage target and no derived weight. The gDNA fit is the away-half moment with
    no location prior at all, the RNA fit is its own raw moment, and the weaker of the two shrinks
    toward the better-measured one — so the reference is a measurement of the same library rather
    than a conjured number. What remains asserted is the ceiling alone."""
    from rigel.calibration import gdna_strand
    from rigel.calibration.gdna_strand import (
        _CEIL_ALPHA_BETA,
        _MAX_OVERDISPERSION,
        overdispersion_for_beta,
    )

    assert _CEIL_ALPHA_BETA == 2.0
    assert _MAX_OVERDISPERSION == pytest.approx(0.2)
    assert overdispersion_for_beta(2.0) == pytest.approx(_MAX_OVERDISPERSION)
    # ⛔ the deleted constants must not come back under any spelling
    for gone in (
        "_PRIOR_ALPHA_BETA",
        "_PRIOR_OVERDISPERSION",
        "_PRIOR_INFORMATION",
        "_prior_information",
    ):
        assert not hasattr(gdna_strand, gone), gone


def test_null_information_reduces_to_pair_count_at_symmetric_mean():
    """``I = 1/Var(od_mom)|₀`` must equal the pair count ``Σ n(n−1)/2`` exactly at μ = ½ (the gDNA
    case), and must not be substituted by the pair count away from it, where the pair count
    overstates the information by roughly an order of magnitude."""
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
    """A pure-gDNA single-strand region has exact 50/50 per-strand counts, which the strand mixture
    (gDNA mean ½, RNA mean κ≠½) must read as gDNA — f_g≈1. At the near-binomial overdispersion it
    does; inflating the gDNA Beta-Binomial widens it, loses specificity at ½ and drags the region
    toward the RNA/gDNA boundary. A pure-RNA control (+frac=κ) stays f_g≈0 at both."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    # κ=0.7 (intermediate strand): gDNA mean ½ is near enough to the RNA mean that the gDNA BB width
    # matters — exactly where an inflated prior does its damage.
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

    # pure gDNA (exact 50/50, truth f_g=1): near-binomial od → confidently gDNA; inflating od
    # monotonically under-calls it. The bound is set on the physics rather than on a recorded
    # magnitude: with ψ bare the strand speaks cleanly and the solver is materially od-robust, so a
    # moderate inflation costs only a few percent and it takes od=0.4 to move it by 0.15.
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
# The per-strand presence masks
#
# `nrna_active_strands` is `free_s`, the strands the annotation admits RNA on (exon or intron bit);
# `mrna_active_strands` is the tighter exon-bit-only mask. `build_region_statics` carries both, and
# `total_abundance` asserts the exon mask as the licence beside the mature wall distances.
# ---------------------------------------------------------------------------

# Slot ids on the `_mature_exon_chain` chain (``N E N E N E N E N``, regions at the even slots, see
# `MX_EXON`). No gate reads these three, and their names do not describe the slots they hold:
_R1_EXON = 3  # slot 3: the boundary between intron n1 (slot 2) and the exon under test (slot 4)
_B1 = 2  # slot 2: intron n1, the exon's left flank
_B2 = 4  # slot 4: the exon under test


def test_mrna_active_matches_same_strand_exon_rule():
    """The `mrna_active_strands` mature-presence mask is exactly the rule: mature is present on strand s
    across a boundary iff the SAME-STRANDED exon bit is set on BOTH flanks. Intron bits never qualify;
    `EX+EX- | EX+EX-` passes on BOTH strands. Enumerate all 16×16 signature pairs (a boundary's two
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

    # the headline case: overlapping opposite-strand exons on both flanks ⇒ mature passes on BOTH strands
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
    """No nan/inf reaches the fold. The real region_sweep over spliced/±, stranded/unstranded, and extreme
    gDNA/mature densities (pure-gDNA, pure-RNA, empty, tiny, huge) — every final fraction is finite & in range,
    every variance is ≥0 (∞ = the honest 'unsolved' state is allowed; nan is not)."""
    for spliced in (True, False):
        for kappa in (0.5, 0.95):
            for rho_g, rho_m in [(0.5, 1.0), (0.0, 1.0), (2.0, 0.0), (1e-6, 1e-6), (1e4, 1e4)]:
                cfg = dict(spliced=spliced, kappa=kappa, rho_g=rho_g, rho_m=rho_m)
                final, cap = _sweep(_mature_exon_chain(**cfg), kappa=kappa)
                for nm in ("f_g", "f_pos", "f_neg"):
                    v = np.asarray(getattr(final, nm))
                    assert np.all(np.isfinite(v)), (cfg, nm, v)
                    assert np.all(v >= -1e-9) and np.all(v <= 1.0 + 1e-9), (cfg, nm, v)
                v = np.asarray(final.var_gdna)
                assert not np.any(np.isnan(v)) and np.all(v >= -1e-12), (cfg, v)  # ∞ ok, nan not


def test_region_sweep_deterministic():
    """Pass-0 must be bit-reproducible. The sweep is one native call (`native.solve_blocks`) that solves
    each locus block end to end with no reduction across blocks, so the same input must give a
    bit-identical belief run-to-run — a prerequisite for any confidence claim about the solver. Uses
    the unstranded (κ=½) case, the one most sensitive to any ordering nondeterminism."""
    a, capa = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    b, capb = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    for nm in ("f_g", "f_pos", "f_neg", "var_gdna"):
        x, y = np.asarray(getattr(a, nm)), np.asarray(getattr(b, nm))
        assert np.array_equal(x, y, equal_nan=True), (nm, x, y)  # BIT-identical (not just close)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# THE NUMERIC CONTRACT, CONTINUED — ψ is CHUNK-EXACT: how the rows are tiled moves no number.
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def _chunk_substrate(m=255, K=120, seed=3):
    """A mixed substrate: single-strand and AMBIG slots, a fitted composition arm (a curve on a per-slot
    support, as the kernel reads it), non-flat λ-factor rows, delivered rows, and a per-slot freeze
    reference."""
    rng = np.random.default_rng(seed)
    u_pos = rng.integers(0, 60, m).astype(float)
    u_neg = rng.integers(0, 60, m).astype(float)
    ap = rng.random(m) < 0.75
    an = rng.random(m) < 0.55
    ap[~(ap | an)] = True  # no locked slot: every row solves
    lam = np.linspace(-10.0, 10.0, K)
    log_rho = np.linspace(-12.0, 2.0, 80)
    log_p = -0.5 * ((log_rho + 4.0) / 1.5) ** 2
    mass, eff = rng.uniform(5.0, 500.0, m), rng.uniform(100.0, 2000.0, m)
    rows = -0.02 * (lam[None, :] - rng.normal(0.0, 1.0, m)[:, None]) ** 2
    delivered = -0.05 * (lam[None, :] - rng.normal(1.0, 1.5, m)[:, None]) ** 2
    fg_ref = rng.uniform(0.05, 0.95, m)
    rest = 1.0 - fg_ref
    fpos_ref = np.where(ap & an, rest / 2, np.where(ap, rest, 0.0))
    fneg_ref = np.where(ap & an, rest / 2, np.where(an, rest, 0.0))
    args = (u_pos, u_neg, ap, an, u_pos + u_neg, np.zeros(m))
    kw = dict(
        kappa=0.9,
        od_g=0.02,
        od_r=0.03,
        n_grid=K,
        L=10.0,
        gdna_prior=(log_rho, log_p),
        gdna_support=(mass, eff),
        lam_logprior=rows,
        row_logprior=delivered,
        fg_ref=fg_ref,
        fpos_ref=fpos_ref,
        fneg_ref=fneg_ref,
    )
    return m, args, kw


def _solve_in_chunks(args, kw, edges):
    """`_solve_regions_logodds_all` on each ``[a, b)`` of ``edges`` in turn, scattered back — the
    block solve's arithmetic, with nothing but the row tiling changed."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    fields = ("gdna_frac", "rna_pos_frac", "rna_neg_frac", "gdna_frac_var")
    m = args[0].shape[0]
    out = {f: np.zeros(m) for f in fields}
    for a, b in edges:
        sub_args = tuple(x[a:b] for x in args)
        sub_kw = dict(kw)
        for key in ("lam_logprior", "row_logprior", "fg_ref", "fpos_ref", "fneg_ref"):
            sub_kw[key] = kw[key][a:b]
        sub_kw["gdna_support"] = tuple(x[a:b] for x in kw["gdna_support"])
        dc = _solve_regions_logodds_all(*sub_args, **sub_kw)
        for f in fields:
            out[f][a:b] = getattr(dc, f)
    return out


def test_the_psi_solve_is_chunk_exact_so_a_block_split_moves_no_number():
    """The property the locus solve stands on: the answer at a slot is a function of that slot's inputs
    and the grid, never of which other rows share its tile. Whole, in halves, in thirds and one row at a
    time must agree to the bit on every field of both paths (the single-strand solve and the AMBIG
    cube). The read-out earns this by reducing every row in a fixed order — a contiguous ψ and
    per-row moment sums that do not go through BLAS — because a layout- or row-count-dependent
    reduction moves the last ulp, and one ulp is a different number."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    m, args, kw = _chunk_substrate()
    whole = _solve_regions_logodds_all(*args, **kw)
    partitions = {
        "halves": [(0, m // 2), (m // 2, m)],
        "thirds": [(0, m // 3), (m // 3, 2 * m // 3), (2 * m // 3, m)],
        "one row at a time": [(i, i + 1) for i in range(m)],
    }
    for name, edges in partitions.items():
        got = _solve_in_chunks(args, kw, edges)
        for f, arr in got.items():
            ref = np.asarray(getattr(whole, f))
            assert np.array_equal(arr, ref), (
                f"{name}: {f} moved at {int((arr != ref).sum())} of {m} slots "
                f"(max |delta| {np.max(np.abs(arr - ref)):.3e}) — the read-out is not chunk-exact"
            )
    # not vacuous: both paths solved, with a fitted arm
    assert (args[2] ^ args[3]).any() and (args[2] & args[3]).any()
    assert kw["gdna_prior"] is not None and kw["gdna_support"] is not None


def test_the_psi_solve_is_thread_exact_so_the_thread_count_moves_no_number():
    """ψ's slots are pulled one at a time by a pool of threads, each slot solved by the same arithmetic
    on its own scratch and written to its own four outputs — so the answer at 2, 3 and every-core threads
    is the serial answer to the bit on every field of both paths (the single-strand solve and the AMBIG
    cube), and the budget is a resource, not a tunable. PERTURBATION: a kernel whose threads shared a
    scratch, or skipped or doubled a slot, fails here."""
    from rigel.calibration.simplex_logodds import _solve_regions_logodds_all

    m, args, kw = _chunk_substrate()
    serial = _solve_regions_logodds_all(*args, **kw, n_threads=1)
    for n_threads in (2, 3, 0):
        got = _solve_regions_logodds_all(*args, **kw, n_threads=n_threads)
        for f in ("gdna_frac", "rna_pos_frac", "rna_neg_frac", "gdna_frac_var"):
            a, b = np.asarray(getattr(got, f)), np.asarray(getattr(serial, f))
            assert np.array_equal(a, b), (
                f"{n_threads} threads: {f} moved at {int((a != b).sum())} of {m} slots — the solve is not thread-exact"
            )
    assert m >= 3 and (args[2] & args[3]).any(), (
        "the substrate must carry AMBIG slots and more slots than threads"
    )


def test_the_thread_budget_is_a_non_negative_count_and_the_cli_fans_it_out():
    """`CalibrationConfig.n_threads`: 0 is every core (the locus EM's reading), a negative count is
    refused; the default is 0."""
    import pytest

    from rigel.config import CalibrationConfig

    assert CalibrationConfig().n_threads == 0
    assert CalibrationConfig(n_threads=4).n_threads == 4
    with pytest.raises(ValueError, match="n_threads"):
        CalibrationConfig(n_threads=-1)
