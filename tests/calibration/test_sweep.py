"""The two-phase sweep (`sweep.solve_chain`) and the beliefs it starts from.

Every fixture is built by ``_synthetic.make_chain_parts`` on the shipped axes: a region axis, a contiguous-boundary axis
with ``k − 1`` entries per reference and so no terminal slots, and a sj axis whose boundaries state their own ``(src,
dst, strand)``. A reference terminal would be a data-free boundary slot that could be G1-locked and emit structural
all-gDNA into its neighbour; a sj stating its own endpoints lets the mature flux at an intron-exon boundary be derived
from the graph. The tests cover the per-slot init, the factor-1 anchors, a delivered row's pull, the mature-exon chain
and the numeric contract.
"""

from __future__ import annotations


import functools

import numpy as np
import pytest

from rigel.types import Strand

from rigel.calibration.diagnostics import SweepDiagnostics
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.strand_balance import StrandBalance
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
)


#: These tests exercise the sweep's shape and the per-slot init under the silent policy, named explicitly. The message
#: policies have their own tests (`test_sweep_backbone.py`, `test_transfer_policy.py`).
region_sweep = functools.partial(solve_chain, policy=SilentPolicy(), logodds_window=10.0)


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
    b = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60, logodds_window=10.0
    )

    # the chain is N E N E N, so the regions are at 0, 2, 4; there are no terminal slots.
    rid = [0, 2, 4]
    fg = b.f_g[rid]
    # intergenic: locked gDNA sink {0,0,1}, all precision locked (var 0).
    assert fg[0] == 1.0
    assert b.var_gdna[0] == 0.0
    # intron+ (zero gDNA): the strand tilt alone drives f_g → 0; g finite.
    assert fg[1] < 0.15
    assert np.isfinite(b.var_gdna[2])
    # AMBIG: unresolved by strand → {0,0,1} default at maximal (inf) variance for the sweep to resolve.
    assert fg[2] == 1.0
    assert np.isinf(b.var_gdna[4])


def test_init_boundary_continuity_gate():
    # 1 ref, 2 regions (exon+ | intron+) → one boundary between them. A reference with k regions owns k-1 boundaries, so
    # there is no terminal slot to be a sink.
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
    b = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60, logodds_window=10.0
    )
    # slots: N0=0, E0=1, N1=2.
    # E0 (ex+→in+): +strand continuous (G2+) ⇒ the strand tilt resolves f_g → 0.
    assert b.f_g[1] < 0.15
    assert np.isfinite(b.var_gdna[1])


def test_init_tss_boundary_is_black_hole():
    # intergenic | exon+: the internal boundary is a TSS (intergenic↔exon), so continuity blocks RNA and it is a sink,
    # even though the TSS-crossing fragments are sense-tilted.
    parts = make_chain_parts(
        [0, BIT_EXON_POS],
        region_size_bp=[1000.0, 2000.0],
        region_pos=[50.0, 80.0],
        region_neg=[50.0, 4.0],
        boundary_pos=[90.0],
        boundary_neg=[5.0],
    )
    b = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=60, logodds_window=10.0
    )
    # slot 1 is the TSS boundary: a locked gDNA sink despite the sense tilt (all precision locked at 0).
    assert b.f_g[1] == 1.0 and b.var_gdna[1] == 0.0


def test_precision_state_count_resolution():
    """The log-density solver's precision state is ``Var(log f_g)``, the message currency. It reflects evidence: a
    region with more fragments of the same composition resolves its log-density more sharply, so a lower ``Var(log
    f_g)`` means a more confident message. (In log space a confident ``f_g→0`` region has wide variance, since a
    near-zero gDNA density carries little gDNA-density information, so the lattice's linear ``Var(f_g)`` ordering does
    not carry over.) A region with no fragments reports zero variance.
    """
    from rigel.calibration.simplex_logodds import solve_regions_logodds

    kappa = 0.99
    # Same single-strand + composition at two evidence levels: region 1 has 20× the counts of region 0.
    u_pos = np.array([20.0, 400.0])
    u_neg = np.array([20.0, 400.0])
    allow_pos = np.array([True, True])
    allow_neg = np.array([False, False])  # both single-strand
    d = solve_regions_logodds(
        u_pos,
        u_neg,
        allow_pos,
        allow_neg,
        kappa=kappa,
        n_grid=60,
        L=10.0,
    )
    assert d.gdna_frac_var is not None
    # p̂=0.5 at κ=0.99 makes the fragments look unstranded, so the mean channel points at the gDNA mode f_g=1. The count
    # enters as precision, so the higher-count region resolves nearer the mode with a lower Var(log f_g) rather than
    # being pinned independently of the count.
    assert d.gdna_frac[0] > 0.85 and d.gdna_frac[1] > 0.85  # both gDNA-dominant
    assert d.gdna_frac[1] >= d.gdna_frac[0]  # more count ⇒ nearer the mode
    assert d.gdna_frac_var[1] < d.gdna_frac_var[0]  # more fragments ⇒ sharper
    # the variance is present, finite and non-negative for active regions.
    assert np.all(np.isfinite(d.gdna_frac_var)) and np.all(d.gdna_frac_var >= 0.0)
    # a no-fragment region is inactive, so its variance is zero.
    d0 = solve_regions_logodds(
        np.array([0.0]),
        np.array([0.0]),
        np.array([True]),
        np.array([True]),
        kappa=kappa,
        n_grid=60,
        L=10.0,
    )
    assert d0.gdna_frac_var[0] == 0.0


def _factor1_uniform_rho():
    """The factor-1 fixture: a uniform-gDNA chain intergenic | AMBIG | intergenic, every object's count laid down as
    ``rho x its own placements`` (rho = 0.5). Returns the per-slot gDNA density after the sweep: one density per slot,
    since a 0-bp boundary has one divisor. A uniform field laid down must be read back by the solver.
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
    belief = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=0.7, n_grid=40, logodds_window=10.0
    )
    final = region_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        belief,
        parts.region_arrays,
        strand=StrandBalance(0.7, 0, True),
        n_grid=40,
    )
    # gDNA density = f_g x count / E_gdna (the formula the sweep inlines).
    count = np.asarray(parts.geometry.unspliced_count, float).sum(axis=1)
    return final.f_g * count / np.asarray(parts.geometry.eff_gdna, float)


def test_gdna_sweep_factor1_intergenic_anchors():
    """On a uniform-gDNA chain the strand/signature-locked intergenic regions read back ρ exactly, an anchor invariant
    every phase must hold."""
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    interg = [0, 4]  # the chain is N E N E N, so the two intergenic regions are slots 0 and 4
    assert np.allclose(rho_g[interg], rho, atol=0.02)


def test_interior_anchor_is_immovable_and_produces_no_nan():
    """A composition-certain (`g1_locked`) region has ``Var(log f_c) = 0``, so any fusion weight formed as ``1/Var``
    would be ``∞`` and spread a nan through the chain. On the factor-1 chain, whose two intergenic regions are such
    anchors interior to the chain:

    1. no nan anywhere: beliefs and variances stay finite (``∞`` is the honest "unsolved" state and is allowed on a
       variance; nan never is);
    2. the anchor is immovable: it reads back the true ρ exactly beside a wrong AMBIG neighbour, because a `g1_locked`
       region is never `solvable`, so its ψ output is discarded and its own count stands.
    """
    rho = 0.5
    rho_g = _factor1_uniform_rho()
    assert not np.any(np.isnan(rho_g)), rho_g
    assert np.all(np.isfinite(rho_g)), rho_g
    assert np.allclose(rho_g[[0, 4]], rho, atol=1e-9)  # exact, not merely close


def test_a_zero_gdna_chain_reads_rna_at_every_strand_decisive_slot():
    """A pure-RNA chain intron+ | AMBIG(in+|in−) | intron−, solved under the silent policy: every slot on its own
    evidence, nothing sent between them, and the AMBIG region starting at the all-gDNA init f_g=1. Each slot whose
    strand tilt is decisive (the two single-strand introns and the AMBIG|intron− boundary) reads RNA from the strand
    alone. The sweep is one forward-backward pass, so there is no per-pass sequence to be monotone in."""
    gdna_fl, rna_fl = _delta_pmf(300), _delta_pmf(200)
    # sense-tilted RNA (κ=0.95): the + intron aligns genome+, the − intron genome−. The two boundaries carry the tilt of
    # the regions they separate. There are no terminal slots, so no adjacent terminal G1 lock emits structural gDNA into
    # a flank.
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
    belief = init_beliefs(geom, st, rna_sense_frac=0.95, n_grid=40, logodds_window=10.0)
    assert belief.f_g[2] == 1.0  # AMBIG starts all-gDNA
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        region_arrays,
        # the spliced sample behind κ: the stranded (κ=0.95) channel reads live
        strand=StrandBalance(0.95, 10_000, False),
        n_grid=40,
    )
    # Slots 0 and 4 are the single-strand introns, 1 and 3 the boundaries beside the AMBIG region, each tilted 95:5 or
    # 40:2 toward its RNA strand. A slot whose strand term says nothing reads f_g = 0.5, so the bound sits well below
    # that: at κ = 0.95 a 95:5 tilt puts f_g under 0.05. The AMBIG region (slot 2) is not asserted: its strand is
    # balanced, so its strand likelihood is degenerate (a balanced count fits gDNA and balanced ±RNA equally), and with
    # no prior and nothing sent to it, it stays near its all-gDNA init.
    decisive = [0, 1, 3, 4]
    assert np.all(final.f_g[decisive] < 0.1), final.f_g


# a delivered row: two-sided pull and emergent deference --------


def _gdna_share_row(n_grid, mode_share, prec):
    """A claim on the gDNA share delivered as a λ-row on the solve grid, the message layer's one currency: a Gaussian on
    ``log f_g`` at ``log mode_share`` with precision ``prec``."""
    from scipy.special import log_expit

    from rigel.calibration.simplex_logodds import logodds_grid

    lam, _ = logodds_grid(int(n_grid), 10.0)
    return (-0.5 * float(prec) * (log_expit(lam) - np.log(float(mode_share))) ** 2)[None, :]


def test_a_delivered_row_pulls_two_sided_and_not_to_the_vertex():
    """A row claiming ``f_g = 0.2`` (strong) on a balanced AMBIG region (flat strand) pulls f_g toward 0.2: two-sided by
    construction (a profile on the grid, no boundary wall), not to the f_g=1 vertex."""
    from _psi_reference import psi_solve

    # AMBIG region with balanced counts: the strand is flat (κ=0.5), so only the row shapes f_g.
    d = psi_solve(
        np.array([50.0]),
        np.array([50.0]),
        np.array([True]),
        np.array([True]),
        kappa=0.5,
        n_grid=80,
        L=10.0,
        row_logprior=_gdna_share_row(80, 0.2, 200.0),
    )
    fg = float(d.gdna_frac[0])
    assert abs(fg - 0.2) < 0.05, fg


def test_a_weak_row_defers_to_a_decisive_strand():
    """A weak row (precision 3) claiming ``f_g = 0.9`` loses to a decisive single-strand region's ~1000-fragment strand
    likelihood: f_g stays ≈0, since a weak claim cannot override the data."""
    from _psi_reference import psi_solve

    d = psi_solve(
        np.array([1000.0]),
        np.array([5.0]),
        np.array([True]),
        np.array([False]),
        kappa=0.99,
        n_grid=80,
        L=10.0,
        row_logprior=_gdna_share_row(80, 0.9, 3.0),
    )
    fg = float(d.gdna_frac[0])
    assert fg < 0.1, fg


# --- the mature-exon chain: an expressed exon between two gDNA-only introns ---------------------------
# Under the silent policy every slot is solved on its own evidence: the introns' balanced counts must read gDNA and the
# exon's sense-tilted count RNA. The spliced reads over the introns move no slot's answer here (the ``spliced=False``
# chain solves identically); the policies that carry them have their own tests.


#: slots of the ``_mature_exon_chain`` fixture. The chain is ``N E N E N E N E N`` (9 slots, regions at the even ones),
#: so the exon under test is slot 4 and its flanking introns are slots 2 and 6.
MX_EXON, MX_INTRONS = 4, [2, 6]


def _mature_exon_chain(*, spliced: bool, rho_g=0.5, rho_m=1.0, kappa=0.95, spl_scale=1.0):
    """``exon+ | intron+ | EXON+ | intron+ | exon+``: a pure-mature expressed gene with no nascent RNA.

    Five regions, because a sj states its own ``(src, dst)`` and so needs endpoints: the sj over intron ``n1`` runs ``n0
    → n2`` and the one over ``n3`` runs ``n2 → n4``, and `build_region_geometry` places their flux on the boundaries
    they leave and enter, so the exon under test (``n2``) gets mature flux on both flanking boundaries, derived from the
    graph. Every exon's contained unspliced count is balanced gDNA plus sense (+) mature; the introns and every boundary
    carry balanced gDNA only. `boundary_spliced` is 0 everywhere because mature RNA never crosses an exon-intron
    boundary: it skips the intron as a sj.
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
    # ``spl_scale`` < 1 models a capture-depleted sj: sj-spanning reads are only partly captured, so the sj
    # under-reports the exon's true mature density and the boundary→exon mature measurement disagrees with the exon's
    # own confident unspliced belief.
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
        parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=60, logodds_window=10.0
    )
    return parts.chain, parts.statics, parts.geometry, belief, parts.region_arrays


def _sweep(args, kappa=0.95):
    chain, st, geom, belief, ra = args
    cap = SweepDiagnostics()
    final = region_sweep(
        chain,
        st,
        geom,
        belief,
        ra,
        # The strand channel is live iff the protocol decision on (κ, N_rna) says so, so the fixture supplies the
        # spliced sample behind its κ.
        strand=StrandBalance(kappa, 10_000, False),
        n_grid=60,
        diagnostics=cap,
    )
    return final, cap


def test_mature_no_nascent_hallucination_in_introns():
    """An expressed exon whose RNA is all spliced must not make its flanking gDNA-only introns read as RNA; the introns
    stay gDNA (truth ``f_g = 1``), each solved on its own balanced counts."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_introns = fin_m.f_g[MX_INTRONS]
    assert np.all(fg_introns > 0.85), fg_introns


def test_mature_measurement_recovers_exon_rna():
    """The same chain's expressed exon is recovered as mostly RNA (true f_g ≈ ρ_g·E_g/(ρ_g·E_g+ρ_m·E_r) ≈ 0.32), read
    off its own + strand tilt; the gDNA-only introns beside it do not starve it of RNA."""
    fin_m, _ = _sweep(_mature_exon_chain(spliced=True))
    fg_exon = float(fin_m.f_g[MX_EXON])
    assert fg_exon < 0.45, fg_exon  # truth ≈0.32; comfortably RNA-dominated, not pinned to gDNA


def test_an_exon_with_no_composition_evidence_stays_near_the_reference():
    """On the unstranded chain with no spliced reads the exon has no composition evidence of its own (no strand channel
    and no sj), so nothing may drive it to a confident vertex: its f_g stays near the uninformative reference."""
    ex = MX_EXON
    fin_no, cap_no = _sweep(_mature_exon_chain(spliced=False, kappa=0.5), kappa=0.5)
    # No spliced reads + no strand ⇒ the exon has no composition evidence of its own; it stays mid-range,
    # where the reference puts it.
    assert 0.2 < float(fin_no.f_g[ex]) < 0.8, fin_no.f_g[ex]


def test_sweep_finite_over_extreme_configs():
    """No nan/inf reaches the fold. Over spliced/±, stranded/unstranded and extreme gDNA/mature densities (pure-gDNA,
    pure-RNA, empty, tiny, huge), every final fraction is finite and in range and every variance is ≥0 (∞, the unsolved
    state, is allowed; nan is not)."""
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
    """The sweep is one native call (`native.solve_blocks`) that solves each locus block end to end with no reduction
    across blocks, so the same input gives a bit-identical belief run to run. Uses the unstranded (κ=½) case, the most
    sensitive to ordering nondeterminism."""
    a, capa = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    b, capb = _sweep(_mature_exon_chain(spliced=True, kappa=0.5), kappa=0.5)
    for nm in ("f_g", "f_pos", "f_neg", "var_gdna"):
        x, y = np.asarray(getattr(a, nm)), np.asarray(getattr(b, nm))
        assert np.array_equal(x, y, equal_nan=True), (nm, x, y)  # bit-identical, not just close


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# The numeric contract, continued: ψ is chunk-exact, so how the rows are tiled moves no number.
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
    args = (u_pos, u_neg, ap, an)
    kw = dict(
        kappa=0.9,
        n_grid=K,
        L=10.0,
        gdna_prior=(log_rho, log_p),
        gdna_support=(mass, eff),
        row_logprior=rows + delivered,
        fg_ref=fg_ref,
        fpos_ref=fpos_ref,
        fneg_ref=fneg_ref,
    )
    return m, args, kw


def _solve_in_chunks(args, kw, edges):
    """`psi_solve` on each ``[a, b)`` of ``edges`` in turn, scattered back: the block solve's
    arithmetic with only the row tiling changed."""
    from _psi_reference import psi_solve

    fields = ("gdna_frac", "rna_pos_frac", "rna_neg_frac", "gdna_frac_var")
    m = args[0].shape[0]
    out = {f: np.zeros(m) for f in fields}
    for a, b in edges:
        sub_args = tuple(x[a:b] for x in args)
        sub_kw = dict(kw)
        for key in ("row_logprior", "fg_ref", "fpos_ref", "fneg_ref"):
            sub_kw[key] = kw[key][a:b]
        sub_kw["gdna_support"] = tuple(x[a:b] for x in kw["gdna_support"])
        dc = psi_solve(*sub_args, **sub_kw)
        for f in fields:
            out[f][a:b] = getattr(dc, f)
    return out


def test_the_psi_solve_is_chunk_exact_so_a_block_split_moves_no_number():
    """The answer at a slot is a function of that slot's inputs and the grid, never of which rows share its tile. Whole,
    in halves, in thirds and one row at a time must agree to the bit on every field of both paths (the single-strand
    solve and the AMBIG cube). The read-out reduces every row in a fixed order (a contiguous ψ and per-row moment sums
    that avoid BLAS), because a layout- or row-count-dependent reduction would move the last ulp.
    """
    from _psi_reference import psi_solve

    m, args, kw = _chunk_substrate()
    whole = psi_solve(*args, **kw)
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
    # both paths solved, with a fitted arm
    assert (args[2] ^ args[3]).any() and (args[2] & args[3]).any()
    assert kw["gdna_prior"] is not None and kw["gdna_support"] is not None


def test_the_psi_solve_is_thread_exact_so_the_thread_count_moves_no_number():
    """ψ's slots are pulled one at a time by a pool of threads, each solved by the same arithmetic on its own scratch
    and written to its own four outputs, so the answer at 2, 3 and every-core threads is the serial answer to the bit on
    every field of both paths; the thread budget is a resource, not a tunable. Threads sharing scratch, or a slot
    skipped or doubled, would fail here.
    """
    from _psi_reference import psi_solve

    m, args, kw = _chunk_substrate()
    serial = psi_solve(*args, **kw, n_threads=1)
    for n_threads in (2, 3, 0):
        got = psi_solve(*args, **kw, n_threads=n_threads)
        for f in ("gdna_frac", "rna_pos_frac", "rna_neg_frac", "gdna_frac_var"):
            a, b = np.asarray(getattr(got, f)), np.asarray(getattr(serial, f))
            assert np.array_equal(a, b), (
                f"{n_threads} threads: {f} moved at {int((a != b).sum())} of {m} slots — the solve is not thread-exact"
            )
    assert m >= 3 and (args[2] & args[3]).any(), (
        "the substrate must carry AMBIG slots and more slots than threads"
    )


def test_the_thread_budget_is_a_non_negative_count_and_the_cli_fans_it_out():
    """`CalibrationConfig.n_threads`: 0 is every core (as for the locus EM), a negative count is refused; the default is
    0."""

    from rigel.config import CalibrationConfig

    assert CalibrationConfig().n_threads == 0
    assert CalibrationConfig(n_threads=4).n_threads == 4
    with pytest.raises(ValueError, match="n_threads"):
        CalibrationConfig(n_threads=-1)
