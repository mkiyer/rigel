"""The Stage-3 message policy — three ABUNDANCES, recipient-decides, one static population table.

       Gate: ``tests/calibration/test_currency_policy.py``

⭐⭐⭐ **THE ARCHITECTURE (owner, 2026-08-19).** A message is the three ABUNDANCES ``{gDNA, RNA+,
RNA-}`` (counts/bp — ⛔ the word is ABUNDANCE; "level" is retired vocabulary) with per-component
precisions. **The sender sends blindly; the RECIPIENT decides** how to interpret, by asking a STATIC
table one question per hop: *is my transcript population the same as the source's?*

* a REGION's population is the transcripts covering it; a BOUNDARY's population is the transcripts
  that CROSS it — a transcript STARTING at the boundary is not in it. So
  ``pop(B) = pop(left flank) ∩ pop(right flank)``, and a hop ``R ↔ B`` is population-equal iff **R's
  flank gains nothing at B** — one bit per (boundary, side), shared by both directions of the pair,
  and it is :func:`~..region_geometry.terminus_flank_gain`'s own quantity. ⛔ The raw terminus BIT is
  not the licence: a boundary carrying another transcript's TSS is still population-equal with the
  flank that gains nothing (the owner's worked example, gated verbatim).
* populations EQUAL ⇒ the **COMPOSITION strategy**: rescale by the ratio of the two slots' MEASURED
  total abundances (no belief enters the rescale — owner ruling) with the §3.5e splice-flux
  arithmetic at sj boundaries; transfer variance from the binomial up/down-sampling derivation
  (concept C — not yet wired).
* populations DIFFERENT ⇒ the **ABUNDANCE strategy**: the abundances cross unscaled — they are claims
  about the components the source could see, and the destination's residual mass stays UNACCOUNTED
  (there is no mass rescale in this policy, structurally); transfer variance from the total-abundance
  disagreement (concept D — not yet wired).

**Sources are the slots' INITIAL BELIEFS** (``ctx.own`` — the g1 slots, the intron factory's
deconvolution at init, the strand model on stranded data); the policy invents no source. **Lifecycle**
(owner's walkthrough): a slot with NO belief RELAYS the running message unmodified; a slot WITH one
fuses precision-weighted, so a strong local belief dominates ("the message will likely dampen or die
there"); at solve time the two directional messages combine, and a beliefless slot's combined message
becomes its belief (ψ's fuse does that downstream).

⛔ A refused or inadmissible claim loses its VALUE and its PRECISION in one statement — there is no
separate precision for a later operator to hand back (``TRAPS: zero-the-precision-with-the-value``,
made structural). ⭐ A zero abundance is an ordinary VALUE at its own precision, so "there is none
here" crosses every hop intact (``TRAPS: a-ratio-cannot-carry-zero`` dissolves: the claim itself is
never a ratio's denominator).

⚠ **CONCEPT LADDER** (owner: test one new concept at a time): A (the population table) and B (the
own-sourced three-arm relay) are wired; C (composition rescale + §3.5e flux + binomial-sampling
variance) and D (disagreement damping) land next, each gated first. The §3.5e operators below are
already gated on the owner's worked numbers and wire at C.

⛔ ``RelayPolicy`` is NOT touched by this file's existence: the A/B is one config value
(``CalibrationConfig.message_policy``) and the shipped default installs the relay unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..region_geometry import terminus_flank_gain
from ..simplex_logodds import _log_fg, _logodds_grid
from ..splice_graph import (
    FLAG_ACCEPTOR_NEG,
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_NEG,
    FLAG_DONOR_POS,
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
)
from . import NeighbourState, PsiMessage, StepContext
from .variance import count_logvar

__all__ = [
    "CurrencyPolicy",
    "HopStructure",
    "composition_sampling_logvar",
    "hop_structure",
    "population_equal_from_left",
    "splice_in_densities",
    "splice_out_unspliced",
]

_EPS = 1.0e-9


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the static tables — structure bits and population equality, annotation only
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class HopStructure:
    """Per-slot, per-strand structure bits, read off the annotation alone in ``prepare()``.

    ⭐ Each strand reads ONLY its own bits — a boundary that is a ``+`` splice site and a ``−``
    terminus is both, independently. ⛔ Where a strand carries BOTH its term and its sj bit, that
    strand's hop is ruled a TERMINUS (owner: a terminus at a splice junction is still a terminus)."""

    term_pos: np.ndarray
    term_neg: np.ndarray
    sj_pos: np.ndarray
    sj_neg: np.ndarray


def hop_structure(kind_is_boundary, flags) -> HopStructure:
    """Classify every slot's structure bits, and ASSERT the classification covers every boundary.

    A BOUNDARY carrying neither a terminus nor a splice-site bit is a broken input, not a hop type
    (``TRAPS: an-object-class-does-not-see-a-terminus`` keys the table on these bits)."""
    is_b = np.asarray(kind_is_boundary, bool)
    fl = np.asarray(flags, np.uint16)
    term_p = (fl & (FLAG_TSS_POS | FLAG_TES_POS)) != 0
    term_n = (fl & (FLAG_TSS_NEG | FLAG_TES_NEG)) != 0
    sj_p = (fl & (FLAG_DONOR_POS | FLAG_ACCEPTOR_POS)) != 0
    sj_n = (fl & (FLAG_DONOR_NEG | FLAG_ACCEPTOR_NEG)) != 0
    uncovered = is_b & ~(term_p | term_n | sj_p | sj_n)
    assert not np.any(uncovered), (
        f"{int(uncovered.sum())} BOUNDARY slot(s) carry neither a terminus nor a splice-site flag; "
        "the hop table is keyed on the structure bits and cannot classify them "
        "(TRAPS: an-object-class-does-not-see-a-terminus)"
    )
    return HopStructure(
        term_pos=term_p & is_b, term_neg=term_n & is_b, sj_pos=sj_p & is_b, sj_neg=sj_n & is_b
    )


def population_equal_from_left(ctx: StepContext) -> np.ndarray:
    """Is the hop ``left[i] -> i`` population-equal? — the COMPOSITION licence, per pair.

    ``pop(B) = pop(left) ∩ pop(right)``, so the hop between a boundary and a flank is equal iff THAT
    FLANK gains no transcript at the boundary. For a BOUNDARY slot the left hop's flank is its LEFT
    one; for a REGION slot the left hop's flank is the region itself, which is its left boundary's
    RIGHT flank — one bit per (boundary, side), read twice (the relay's pair algebra).
    ⭐ Gated verbatim on the owner's worked example: a boundary carrying TB+'s TSS is still
    population-equal with its left flank ({TA+, TC−} on both), and unequal with its right (gains
    TB+)."""
    rgain, lgain = terminus_flank_gain(ctx.boundary_flags)
    is_b = np.asarray(ctx.is_boundary, bool)
    left = np.asarray(ctx.left, np.int64)
    sl = np.clip(left, 0, int(is_b.shape[0]) - 1)
    return np.where(is_b, ~np.asarray(lgain, bool), ~np.asarray(rgain, bool)[sl]) & (left >= 0)


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the §3.5e operators — pure functions, gated on the owner's worked numbers (wired at concept C)
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


def splice_out_unspliced(*, src, dst_unspliced, sj):
    """SPLICE OUT (`EQUATIONS.md` §3.5e(i)): the unspliced message an exon's abundances imply for its
    boundary. ``src`` and ``dst_unspliced`` are ``(gDNA, RNA+, RNA-)`` abundances; ``sj`` is the
    boundary's measured flux ``(SJ+, SJ-)`` at the boundary's own scale.

    The rescale is by the ratio of totals with the sj flux counted in the BOUNDARY's total; the sj is
    removed at the exon's scale; the rescale is undone. Worked: {10, 90, 0} against {2, 1, 0} + SJ+ 17
    → {2, 1, 0} exactly. Returns ``None`` where either total is 0 — a refusal, not a number."""
    src_tot = float(src[0]) + float(src[1]) + float(src[2])
    dst_tot = float(dst_unspliced[0]) + float(dst_unspliced[1]) + float(dst_unspliced[2])
    dst_tot_with_flux = dst_tot + float(sj[0]) + float(sj[1])
    if src_tot <= 0.0 or dst_tot_with_flux <= 0.0:
        return None
    r = src_tot / dst_tot_with_flux
    return (
        float(src[0]) / r,
        float(src[1]) / r - float(sj[0]),
        float(src[2]) / r - float(sj[1]),
    )


def splice_in_densities(*, src_with_flux, rho_tot_dst, rho_tot_src):
    """SPLICE IN (`EQUATIONS.md` §3.5e(ii)) — the inverse: the boundary sends its abundances INCLUDING
    the spliced-in flux (those fragments continue contiguously into the exon), and the exon rescales
    by ``rho_tot_dst / rho_tot_src``. Worked: {2, 18, 0} at totals 100/20 → {10, 90, 0} exactly.
    Returns ``None`` where either total is 0."""
    if float(rho_tot_src) <= 0.0 or float(rho_tot_dst) <= 0.0:
        return None
    r = float(rho_tot_dst) / float(rho_tot_src)
    return tuple(float(c) * r for c in src_with_flux)


def enrichment_ratio(ctx: StepContext, *, backward: bool) -> np.ndarray:
    """The per-hop ENRICHMENT RATIO ``T_dst / T_src``, from the MODEL-FREE total abundance.

    ⛔⛔ **NEVER ``mass / effective_length``.** That divisor is a function of the composition being
    solved for — gDNA and RNA carry different fragment-length distributions, so 100 counts in a 500 bp
    region reads 0.25 as pure gDNA and 0.33 as pure RNA — and an enrichment ratio built from it is
    circular. The accumulator already deposits the RECIPROCAL OPPORTUNITY per fragment (``1/(w−1)``
    crossing a BOUNDARY, ``1/(ell−w+1)`` contained in a REGION), whose expectation is the density
    EXACTLY, for any length distribution and any composition: ``ctx.inv_abundance``.
    ⚠ An equal-fragment-length panel cannot see the difference, which is why this had to be reasoned
    about rather than measured (owner, 2026-08-20).

    Returns 1.0 (perfect agreement — no evidence of enrichment) where either end reads nothing."""
    inv = np.asarray(ctx.inv_abundance, np.float64)
    nbr = np.asarray(ctx.right if backward else ctx.left, np.int64)
    src = np.clip(nbr, 0, inv.shape[0] - 1)
    t_src, t_dst = inv[src], inv
    ok = (nbr >= 0) & (t_src > 0.0) & (t_dst > 0.0)
    return np.where(ok, t_dst / np.where(ok, t_src, 1.0), 1.0)


def rescale_weight(*, log_ratio: float, var_ratio: float) -> float:
    """⭐⭐⭐ **THE KNOB — how much of the observed enrichment to believe. One number, no tunable.**

    The ABUNDANCE strategy and the COMPOSITION strategy are not two mechanisms to choose between: they
    are the two ENDS of one continuum (owner, 2026-08-20). ABUNDANCE says the enrichment is 1 —
    transport the abundances as they are; COMPOSITION says it is exactly what the totals report —
    rescale by the ratio. Both are point hypotheses about the same unknown, the log enrichment ``eta``.

    Fusing them by inverse variance — the observed ``log r`` with its sampling variance ``v``, against
    the no-enrichment premise whose error, if it is wrong, is ``log r`` itself — gives a SHRINKAGE
    estimator with no free constant::

        w  =  (log r)^2 / ( (log r)^2 + v )        eta_hat = w · log r        Var(eta_hat) = w · v

    so a disagreement that is indistinguishable from counting noise is shrunk to nothing (the
    capture-OFF end, ``w → 0``: transport the ABUNDANCE) and one that dwarfs the noise is believed in
    full (the capture-ON end, ``w → 1``: transport the COMPOSITION). ⭐ **The delivered claim is
    ``rho · r^w`` — a continuous interpolation in log space between the two strategies, chosen per hop
    BY THE DATA**, which is what makes this a derivation rather than a switch."""
    lr2 = float(log_ratio) * float(log_ratio)
    v = float(var_ratio)
    if lr2 <= 0.0:
        return 0.0
    return lr2 / (lr2 + v) if v > 0.0 else 1.0


def composition_rescale_factor(*, rho_g, rho_p, rho_n, E_g_dst, E_r_dst, M_dst):
    """The COMPOSITION strategy's rescale ``k``, EXACT under the premise and free of the DESTINATION's
    belief — which is what makes it admissible where the relay's mass rescale was not.

    Under the premise the two slots share a population and differ by ONE common enrichment factor
    ``k``: ``rho_c,dst = k · rho_c,src`` for every component. Laying the SOURCE's abundances down in
    the DESTINATION's geometry and requiring them to account for the destination's OBSERVED mass
    determines it::

        S  =  rho_g,src · E_g,dst  +  (rho_+,src + rho_-,src) · E_r,dst
        k  =  M_dst / S

    ⭐ **Only the source's BELIEF and the destination's OBSERVATIONS + CONSTANTS enter** — `StepContext`
    permits exactly that, and the destination's own belief never touches it
    (``TRAPS: a-message-from-the-destinations-belief``). ⛔ It is NOT the relay's mass rescale, which
    filled the components the source did not supply from the DESTINATION's belief and thereby turned a
    refused claim into "all your mass is gDNA".
    ⭐ It reproduces `EQUATIONS.md` §3.5e's worked numbers exactly in both directions, and unlike a
    ratio of ``M/E_g`` totals it is correct per arm when the two components' opportunities differ —
    which between a REGION and a BOUNDARY they always do.
    Returns ``None`` where the premise cannot be evaluated (no mass, or a source claiming nothing)."""
    S = float(rho_g) * float(E_g_dst) + (float(rho_p) + float(rho_n)) * float(E_r_dst)
    if S <= 0.0 or float(M_dst) <= 0.0:
        return None
    return float(M_dst) / S


def _sampling_logvar_array(n_c, n_tot):
    """:func:`composition_sampling_logvar` for arrays — ONE formula, two shapes (the scalar twin is
    what the sequential scan reads; this is what the vectorised combine reads)."""
    return count_logvar(np.asarray(n_c, np.float64)) - count_logvar(
        np.asarray(n_tot, np.float64) + 1.0
    )


def composition_sampling_logvar(*, n_c: float, n_tot: float) -> float:
    """the SAMPLING variance — the COMPOSITION transfer's own: binomial up/down-sampling (owner's derivation ask).

    The rescale transfers the source's SHARES, and a share estimated from the source's raw counts
    keeps the raw counts' precision however large the scale factor is ("the precision is still based
    on 2 fragments"). Under the Dirichlet-Jeffreys posterior over the three components
    (``alpha_c = 1/2`` each — AXIOM 0 fixes K = 3),
    ``Var(log share_c) = trigamma(n_c + 1/2) - trigamma(n_tot + 3/2)`` — every constant is Jeffreys,
    and a zero-count component gets a WIDE but finite variance, never a refusal."""
    return float(count_logvar(n_c) - count_logvar(n_tot + 1.0))


# ══════════════════════════════════════════════════════════════════════════════════════════════════════
# the policy
# ══════════════════════════════════════════════════════════════════════════════════════════════════════


class _CurrencyRelay:
    """One sweep's working object: the static tables, the initial beliefs, and the grid domain."""

    def __init__(self, ctx: StepContext):
        self.ctx = ctx
        own = ctx.own
        # ── the SOURCES: each slot's initial belief, exactly as region_init built it ─────────────────
        self._val0 = tuple(np.asarray(a, np.float64) for a in (own.rho_g, own.rho_pos, own.rho_neg))
        self._prec0 = tuple(
            np.asarray(a, np.float64) for a in (own.prec_g, own.prec_pos, own.prec_neg)
        )
        self._fp = np.asarray(ctx.free_pos, bool)
        self._fn = np.asarray(ctx.free_neg, bool)
        # ── the static tables ────────────────────────────────────────────────────────────────────────
        self.structure = hop_structure(ctx.is_boundary, ctx.boundary_flags)
        eq_l = population_equal_from_left(ctx)
        right = np.asarray(ctx.right, np.int64)
        sr = np.clip(right, 0, int(eq_l.shape[0]) - 1)
        self._eq_l, self._eq_r = eq_l, np.where(right >= 0, eq_l[sr], False)
        # ── the delivery frame (observations) and the grid domain ───────────────────────────────────
        self._M = np.asarray(ctx.mass, np.float64)
        self._E_g = np.asarray(ctx.eff_gdna_global, np.float64)
        self._E_r = np.asarray(ctx.eff_rna, np.float64)
        # ── the BELIEF-FREE face totals and flank flux (owner: no belief enters the rescale) ─────────
        # One total per FLANK (§3.6c's pairing): a face's total is the slot's unspliced abundance in
        # the gDNA frame plus the sj flux whose body lies on that side. ⚠ M/E_g mixes the RNA share
        # through the gDNA divisor — a deliberate belief-free approximation; the two sides of a hop
        # carry the same convention, so the frame error is common-mode in the RATIO.
        g = ctx.geometry

        def _flux(cnt, eff):
            c = np.asarray(cnt, np.float64)
            e = np.asarray(eff, np.float64)
            live = (c > 0.0) & (e > 0.0)
            return np.where(live, c / np.where(live, e, 1.0), 0.0)

        if g is not None and hasattr(g, "sj_count_lo"):
            self._flux_lo = (
                _flux(g.sj_count_lo[:, 0], g.eff_sj_lo[:, 0]),
                _flux(g.sj_count_lo[:, 1], g.eff_sj_lo[:, 1]),
            )
            self._flux_hi = (
                _flux(g.sj_count_hi[:, 0], g.eff_sj_hi[:, 0]),
                _flux(g.sj_count_hi[:, 1], g.eff_sj_hi[:, 1]),
            )
        else:  # a substrate with no sj banks (unit fixtures): no flux anywhere
            z = np.zeros(int(self._M.shape[0]))
            self._flux_lo = (z, z.copy())
            self._flux_hi = (z.copy(), z.copy())
        # ── the PER-HOP TABLES, both directions, indexed by the DESTINATION slot ─────────────────────
        # ⭐ ONE derivation read by BOTH twins (the sequential scan and the vectorised deliver) — the
        # relay's scalar/vector twin drift is what this precomputation exists to prevent.
        n = int(self._M.shape[0])
        is_b = np.asarray(ctx.is_boundary, bool)
        self._tables = {}
        for backward in (False, True):
            nbr = np.asarray(ctx.right if backward else ctx.left, np.int64)
            sn = np.clip(nbr, 0, n - 1)
            fdst = self._flux_hi if backward else self._flux_lo
            fsrc = self._flux_lo if backward else self._flux_hi
            # ⭐ the MODEL-FREE enrichment ratio and the variance of its log — the knob's two inputs.
            # The variance is the counting variance of the two abundance observations in the tree's own
            # convention (`count_logvar`, exact for the Gamma-Jeffreys posterior), so the knob carries
            # NO constant of its own: it is entirely trigamma of observed counts.
            r_enr = enrichment_ratio(ctx, backward=backward)
            n_obs = np.asarray(ctx.n_slot, np.float64)
            self._tables[backward] = {
                "log_r": np.log(np.maximum(r_enr, _EPS)),
                "v_r": count_logvar(n_obs) + count_logvar(n_obs[sn]),
                "eq": self._eq_r if backward else self._eq_l,
                "fx_p": np.where(is_b, fdst[0], fsrc[0][sn]),
                "fx_n": np.where(is_b, fdst[1], fsrc[1][sn]),
                "dst_is_b": is_b,
                "src": sn,
            }
        lam, _ = _logodds_grid(int(ctx.n_grid), float(ctx.logodds_window))
        dom = _log_fg(lam)
        self._dom_lo, self._dom_hi = float(dom[0]), float(dom[-1])
        cap = ctx.capture
        if cap is not None:  # the census hook — an instrument can read the tables this sweep used
            cap.setdefault("_currency_static", {}).update(
                pop_equal_from_left=eq_l.copy(),
                term_pos=self.structure.term_pos.copy(),
                term_neg=self.structure.term_neg.copy(),
                sj_pos=self.structure.sj_pos.copy(),
                sj_neg=self.structure.sj_neg.copy(),
                own_prec_g=self._prec0[0].copy(),
            )

    def scan(self, *, backward: bool):
        """One directional pass of the three-abundance message.

        Concept B: the value crosses every hop unscaled (C adds the composition rescale on
        population-equal hops; D adds the disagreement variance on unequal ones). At the destination,
        an inadmissible strand's claim is zeroed — value AND precision, one statement — and the claim
        fuses with the destination's own belief precision-weighted, so a beliefless slot relays
        unmodified and a strong local belief dominates."""
        fp_l, fn_l = self._fp.tolist(), self._fn.tolist()
        # the per-hop tables for THIS direction — forward: the source is left[i], the source presents
        # its HI face and the destination its LO (backward is the mirror); the flux on the hop is the
        # BOUNDARY end's face flux (§3.6c's role pairing). Shared with deliver — ONE derivation.
        tab = self._tables[backward]
        log_r_l = tab["log_r"].tolist()
        v_r_l = tab["v_r"].tolist()
        eq = tab["eq"].tolist()
        fx_p = tab["fx_p"].tolist()
        fx_n = tab["fx_n"].tolist()
        dst_is_b = tab["dst_is_b"].tolist()
        E_g_l, E_r_l, M_l = self._E_g.tolist(), self._E_r.tolist(), self._M.tolist()
        # the RUNNING state (mutated in place) and the OWN beliefs (read-only), as Python lists
        vg, vp, vn = (a.tolist() for a in self._val0)
        pg, pp, pn = (a.tolist() for a in self._prec0)
        og, op_, on = (a.tolist() for a in self._val0)
        qg, qp, qn = (a.tolist() for a in self._prec0)

        def _fuse(ov: float, oq: float, tv: float, tp: float) -> tuple[float, float]:
            p = oq + tp
            return ((oq * ov + tp * tv) / p, p) if p > _EPS else (ov, 0.0)

        def _cap(p: float, v: float, w: float = 1.0, v_r: float = 0.0) -> float:
            """The transported claim's precision. The SAMPLING variance is a CAP, not an accumulating
            penalty — the claim may never be more precise than its implied source counts support ("the
            precision is still based on 2 fragments"), and a cap is idempotent, so a claim relayed over
            many hops is not damped once per hop for the same original counts. On top of it sits
            ``w·v_r``: what the knob paid for believing the ratio it just applied."""
            if p <= 0.0:
                return 0.0
            if v > 0.0:
                p = min(p, 1.0 / v)
            add = float(w) * float(v_r)
            return p / (1.0 + p * add) if add > 0.0 else p

        def _damp_disagreement(p: float, lr: float, w: float) -> float:
            """The other half of the same statement (owner: *"the greater the disagreement, the less
            trustworthy the message … dampen the message precision"*). Whatever fraction of the observed
            enrichment the knob did NOT absorb is an unexplained mislift the claim still carries:
            ``(1 − w)·log r`` in log space, variance ``((1−w)·log r)²``. ⭐ At the COMPOSITION end
            (``w → 1``) it vanishes — the rescale accounted for it — and at the ABUNDANCE end
            (``w → 0``) it is the full squared disagreement. ONE expression spans the continuum, which
            is what makes the two strategies one mechanism rather than two."""
            resid = (1.0 - float(w)) * float(lr)
            v = resid * resid
            return p / (1.0 + p * v) if (p > 0.0 and v > 0.0) else p

        def step(s: int, i: int) -> None:
            tg, tpg = vg[s], pg[s]
            tp_, tpp = vp[s], pp[s]
            tn_, tpn = vn[s], pn[s]
            # ⭐⭐ THE COMPOSITION LICENCE IS TWO CONDITIONS, AND THE SECOND IS NOT OPTIONAL:
            #  (a) the populations are EQUAL (the static table), and
            #  (b) the source SUPPLIED every component the destination admits — a statement about
            #      PRECISION, not about a value.
            # Without (b) the mass identity has no right-hand side: the k-form would scale a PARTIAL
            # claim up to account for mass its missing components hold, which is precisely the
            # "all your mass is gDNA" pathology (measured here: it invented 73,728 gDNA fragments at a
            # zero-gDNA control before the supply test landed). Unsupplied ⇒ the ABUNDANCE strategy.
            supplied = tpg > 0.0 and (tpp > 0.0 or not fp_l[i]) and (tpn > 0.0 or not fn_l[i])
            # ⭐⭐⭐ THE KNOB: how much of the observed enrichment to believe, per hop, FROM THE DATA.
            # ``w → 0`` is the ABUNDANCE end (transport as-is), ``w → 1`` the COMPOSITION end (full
            # rescale), and every point between is reachable. ⛔ It is only ASKED where the composition
            # is licensed at all — populations equal AND every admissible component supplied; elsewhere
            # it is 0 by structure and the whole disagreement becomes a VARIANCE instead.
            lr, vr = log_r_l[i], v_r_l[i]
            w = rescale_weight(log_ratio=lr, var_ratio=vr) if (eq[i] and supplied) else 0.0
            if w > 0.0:
                # ── the COMPOSITION strategy (§3.5e) ────────────────────────────────────────────────
                # The population ENTERING the destination: a SPLICE IN carries the source boundary's
                # flux with it (those fragments continue contiguously into the exon); a SPLICE OUT
                # does not (they splice out and land elsewhere) and the flux is removed AFTER the
                # rescale, at the destination's own scale.
                in_p = tp_ + (0.0 if dst_is_b[i] else fx_p[i])
                in_n = tn_ + (0.0 if dst_is_b[i] else fx_n[i])
                # ⭐ the destination's OBSERVED mass for THIS hop's population: at a SPLICE OUT the
                # fragments that splice out here ARE observed at the boundary, as its sj flux, so the
                # face's mass is flux-inclusive (`EQUATIONS.md` §3.6c's per-flank total). At a SPLICE
                # IN the flux is already inside the message and the region's own mass is the total.
                m_face = M_l[i] + (fx_p[i] + fx_n[i] if dst_is_b[i] else 0.0)
                kk_full = composition_rescale_factor(
                    rho_g=tg,
                    rho_p=in_p,
                    rho_n=in_n,
                    E_g_dst=E_g_l[i],
                    E_r_dst=E_r_l[i],
                    M_dst=m_face,
                )
                # ⭐ the knob applied: ``k^w`` interpolates in LOG SPACE between no rescale (w = 0) and
                # the full mass-identity rescale (w = 1). The exponent is the data's answer, not a flag.
                kk = None if kk_full is None else kk_full**w
                if kk is None:
                    tpg = tpp = tpn = (
                        0.0  # the premise cannot be evaluated — no claim, one statement
                    )
                    tg = tp_ = tn_ = 0.0
                else:
                    tg = tg * kk
                    if dst_is_b[i]:  # SPLICE OUT: the sj flux leaves at the destination's scale
                        tp_ = max(in_p * kk - fx_p[i], 0.0)
                        tn_ = max(in_n * kk - fx_n[i], 0.0)
                    else:  # SPLICE IN: the flux is already in the message
                        tp_ = in_p * kk
                        tn_ = in_n * kk
                    # the sampling variance — the up/down-sampling cap: shares keep the SOURCE counts' precision
                    n_g = vg[s] * E_g_l[s]
                    n_p = vp[s] * E_r_l[s]
                    n_n = vn[s] * E_r_l[s]
                    n_tot = n_g + n_p + n_n
                    tpg = _cap(tpg, composition_sampling_logvar(n_c=n_g, n_tot=n_tot), w, vr)
                    tpp = _cap(tpp, composition_sampling_logvar(n_c=n_p, n_tot=n_tot), w, vr)
                    tpn = _cap(tpn, composition_sampling_logvar(n_c=n_n, n_tot=n_tot), w, vr)
            # populations differ ⇒ the ABUNDANCE strategy: the claim crosses unscaled (concept D adds
            # the disagreement variance here)
            # ⭐ whatever disagreement the knob did NOT absorb is a variance on EVERY hop — this is
            # the ABUNDANCE strategy's transfer variance, and it is the same statement as the knob.
            if lr != 0.0:
                tpg = _damp_disagreement(tpg, lr, w)
                tpp = _damp_disagreement(tpp, lr, w)
                tpn = _damp_disagreement(tpn, lr, w)
            # an inadmissible strand at the destination: no claim — value and precision together
            if not fp_l[i]:
                tp_, tpp = 0.0, 0.0
            if not fn_l[i]:
                tn_, tpn = 0.0, 0.0
            vg[i], pg[i] = _fuse(og[i], qg[i], tg, tpg)
            vp[i], pp[i] = _fuse(op_[i], qp[i], tp_, tpp)
            vn[i], pn[i] = _fuse(on[i], qn[i], tn_, tpn)

        def publish():
            return tuple(np.asarray(a, np.float64) for a in (vg, vp, vn, pg, pp, pn))

        return step, publish

    def _transport(self, nb: NeighbourState, *, backward: bool):
        """The VECTORISED twin of :meth:`scan`'s ``step`` transform — the same shared hop table, so the
        two cannot drift (the predecessor's scalar/vector twins did, deliberately, and it cost a
        session to prove which behaviour was intended).

        ⛔ ``nb.state`` arrives ALREADY INDEXED AT THE SOURCE, so this has no way to read a
        neighbour's relayed belief at the destination — assertion 1, made structural by the
        backbone."""
        vg, vp, vn, pg, pp, pn = (np.asarray(a, np.float64) for a in nb.state)
        tab = self._tables[backward]
        eq = np.asarray(tab["eq"], bool) & np.asarray(nb.valid, bool)
        # (b) the SUPPLY test — see the scan's note: without every admissible component the k-form
        # has no right-hand side and would inflate a partial claim to fill the destination's mass.
        eq = eq & (pg > 0.0) & ((pp > 0.0) | ~self._fp) & ((pn > 0.0) | ~self._fn)
        # ⭐ the KNOB, vectorised — the same shrinkage the scan applies, from the same shared table
        log_r = np.asarray(tab["log_r"], np.float64)
        v_r = np.asarray(tab["v_r"], np.float64)
        lr2 = log_r * log_r
        w = np.where(eq & (lr2 > 0.0), lr2 / np.maximum(lr2 + v_r, _EPS), 0.0)
        fx_p, fx_n = tab["fx_p"], tab["fx_n"]
        dst_is_b = np.asarray(tab["dst_is_b"], bool)
        # the population ENTERING the destination — direction-dependent (a SPLICE IN carries the flux)
        in_p = vp + np.where(dst_is_b, 0.0, fx_p)
        in_n = vn + np.where(dst_is_b, 0.0, fx_n)
        # the k-form: the source's claim laid down in the destination's geometry must account for the
        # destination's OBSERVED mass. Exact under the premise, free of the destination's BELIEF.
        S = vg * self._E_g + (in_p + in_n) * self._E_r
        m_face = self._M + np.where(dst_is_b, fx_p + fx_n, 0.0)
        ok = eq & (S > 0.0) & (m_face > 0.0) & (w > 0.0)
        k_full = np.where(ok, m_face / np.where(S > 0.0, S, 1.0), 1.0)
        k = np.where(ok, k_full**w, 1.0)  # ``k^w`` — the knob, in log space
        out_g = np.where(ok, vg * k, vg)
        out_p = np.where(ok, np.where(dst_is_b, np.maximum(in_p * k - fx_p, 0.0), in_p * k), vp)
        out_n = np.where(ok, np.where(dst_is_b, np.maximum(in_n * k - fx_n, 0.0), in_n * k), vn)
        # a population-equal hop whose premise cannot be evaluated makes NO claim — value AND precision
        dead = eq & ~ok
        # the SAMPLING CAP — up/down-sampling, from the SOURCE's implied counts
        src = tab["src"]
        n_g = vg * self._E_g[src]
        n_p = vp * self._E_r[src]
        n_n = vn * self._E_r[src]
        n_tot = n_g + n_p + n_n

        def _cap(p, n_c):
            v = _sampling_logvar_array(n_c, n_tot)
            lim = np.where(v > 0.0, 1.0 / np.where(v > 0.0, v, 1.0), np.inf)
            out = np.where(ok, np.minimum(p, lim), p)
            # the enrichment estimate's own variance (``w·v_r``) plus the disagreement the knob did
            # NOT absorb (``((1−w)·log r)²``) — one continuum, both ends covered
            resid = (1.0 - w) * log_r
            add = w * v_r + resid * resid
            out = np.where(add > 0.0, out / (1.0 + out * add), out)
            return np.where(dead, 0.0, out)

        out_g = np.where(dead, 0.0, out_g)
        out_p = np.where(dead, 0.0, out_p)
        out_n = np.where(dead, 0.0, out_n)
        return out_g, out_p, out_n, _cap(pg, n_g), _cap(pp, n_p), _cap(pn, n_n)

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        """Precision-weighted combine of the two directional messages, as ψ's log-share channels.

        The mode is the claim's share of the destination's OWN measured mass, clamped to the grid's
        domain — a zero abundance arrives as the lowest expressible share at its own precision, and an
        over-unit claim as the highest. A slot no message reaches gets precision EXACTLY 0."""
        ctx = self.ctx
        lg, lp, ln, lpg, lpp, lpn = self._transport(left, backward=False)
        rg, rp, rn, rpg, rpp, rpn = self._transport(right, backward=True)
        vl = np.asarray(left.valid, bool)
        vr = np.asarray(right.valid, bool)

        def _combine(av, ap, bv, bp):
            ap = np.where(vl, ap, 0.0)
            bp = np.where(vr, bp, 0.0)
            p = ap + bp
            v = np.where(p > _EPS, (ap * av + bp * bv) / np.maximum(p, _EPS), 0.0)
            return v, p

        cg, cpg = _combine(lg, lpg, rg, rpg)
        cp, cpp = _combine(lp, lpp, rp, rpp)
        cn, cpn = _combine(ln, lpn, rn, rpn)

        M, E_g, E_r = self._M, self._E_g, self._E_r
        lo_share = np.exp(self._dom_lo)

        def _mode(v, p, E):
            live = (p > 0.0) & (M > 0.0) & (E > 0.0)
            share = np.where(live, v * E / np.where(live, M, 1.0), 0.0)
            m = np.where(
                live,
                np.clip(np.log(np.maximum(share, lo_share)), self._dom_lo, self._dom_hi),
                0.0,
            )
            return m, np.where(live, p, 0.0)

        mo_g, pr_g = _mode(cg, cpg, E_g)
        mo_p, pr_p = _mode(cp, cpp, E_r)
        mo_n, pr_n = _mode(cn, cpn, E_r)
        # belt: the annotation's admissibility, applied to the delivered claim too — one statement
        mo_p = np.where(self._fp, mo_p, 0.0)
        pr_p = np.where(self._fp, pr_p, 0.0)
        mo_n = np.where(self._fn, mo_n, 0.0)
        pr_n = np.where(self._fn, pr_n, 0.0)
        cap = ctx.capture
        if cap is not None:
            cap.setdefault("_currency", []).append(
                {
                    "cg": cg.copy(),
                    "cp": cp.copy(),
                    "cn": cn.copy(),
                    "pg": pr_g.copy(),
                    "pp": pr_p.copy(),
                    "pn": pr_n.copy(),
                    "mo_g": mo_g.copy(),
                    "mo_p": mo_p.copy(),
                    "mo_n": mo_n.copy(),
                }
            )
        return PsiMessage(
            gdna_mode=mo_g, gdna_prec=pr_g, rna_mode=(mo_p, mo_n), rna_prec=(pr_p, pr_n)
        )


class CurrencyPolicy:
    """The Stage-3 policy under development. It replaces ``RelayPolicy`` when finished (owner);
    the name is a placeholder the owner has deprioritised.

    ⛔ **There is no strategy SWITCH, and that is the design.** A first version had one — the
    COMPOSITION strategy behind a flag — and the measurement said the flag was the wrong shape: the
    rescale helped under capture (0.61x) and hurt off it (8.7x stranded). Both readings are correct,
    because the two strategies are the two ENDS OF ONE CONTINUUM (owner, 2026-08-20), and the point on
    it is a property of the HOP, not of the run. :func:`rescale_weight` is that point, derived by
    shrinkage from the hop's own observed enrichment and its own counting noise."""

    name = "currency"

    def prepare(self, ctx: StepContext) -> _CurrencyRelay:
        return _CurrencyRelay(ctx)
