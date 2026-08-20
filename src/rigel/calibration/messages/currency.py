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

import math
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
    by_terminus = np.where(is_b, ~np.asarray(lgain, bool), ~np.asarray(rgain, bool)[sl])
    # ⭐⭐ THE SECOND POPULATION PREDICATE (owner, 2026-08-18): **a strand ADMISSIBLE on one side of a
    # hop and not the other is a different RNA population by construction.** The terminus flags cannot
    # see it — a gene edge beside an intergenic region need carry no TSS/TES bit at all — so a table
    # keyed on the flags alone licenses a composition across it. ⛔ What that costs is specific: a
    # structurally pure-gDNA slot's composition is a legitimate ``f_g = 1`` and an EMPTY one's
    # abundance is ~0, so transporting its COMPOSITION says "you are all gDNA" where transporting its
    # ABUNDANCE says "there is almost no gDNA here". At a zero-gDNA library the second is true and the
    # first is the whole error.
    fp = np.asarray(ctx.free_pos, bool)
    fn = np.asarray(ctx.free_neg, bool)
    same_strands = (fp == fp[sl]) & (fn == fn[sl])
    return by_terminus & same_strands & (left >= 0)


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
    lo = inv + np.asarray(ctx.inv_sj_lo, np.float64).sum(axis=1)
    hi = inv + np.asarray(ctx.inv_sj_hi, np.float64).sum(axis=1)
    # ⭐⭐ ONE TOTAL PER FACE, and the sj flux is what makes it a total (`EQUATIONS.md` §3.6c's flank
    # pair, in model-free units). Mature RNA cannot cross an exon|intron boundary contiguously, so the
    # exon's mature fragments appear at the boundary as FLUX rather than as crossings; comparing the
    # boundary's unspliced abundance against the exon's total reports a large "depletion" where there
    # is no probe at all. ⛔ Measured: omitting it cost 30.8x at `g50 ss0.50 capture_off`.
    # The pairing is the scan's: travelling low→high the destination presents its LOW face and the
    # source its HIGH one; backward is the mirror. A REGION's two faces are equal (it carries no flux).
    nbr = np.asarray(ctx.right if backward else ctx.left, np.int64)
    src = np.clip(nbr, 0, inv.shape[0] - 1)
    t_dst = hi if backward else lo
    t_src = (lo if backward else hi)[src]
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


#: ⛔⛔ **THE MASS-IDENTITY RESCALE WAS BUILT HERE AND IS DELETED — the record, so it is not rebuilt.**
#: ``k = M_dst / Sigma_c rho_c,src · E_c,dst`` is exact under the premise and free of the destination's
#: belief, and it reproduced §3.5e's worked numbers. ⛔ It is nevertheless WRONG AS A TRANSPORT, and the
#: reason is structural rather than a tuning failure: it rests on the SOURCE's belief being right, so a
#: source whose claim is small returns a huge ``k``. Measured on the test chromosome: a source holding
#: gDNA 3.9e-4 and RNA exactly 0 was rescaled by **235,800x** to account for a 23,889-fragment exon, and
#: since the only non-zero component was gDNA the exon was handed "all your mass is gDNA" — ``f_g``
#: 1.000 against a truth of 0.000 at a ZERO-gDNA library. That is the relay's own mass-rescale
#: pathology re-entering through a different door.
#: ⭐ The transport is the MEASURED enrichment ratio instead (:func:`enrichment_ratio`): it agrees with
#: ``k`` exactly where the source's belief is right — which is what the worked numbers check — and it
#: cannot amplify a claim the source does not have, because it never reads the claim.


def composition_rescale_factor(*, rho_g, rho_p, rho_n, E_g_dst, E_r_dst, M_dst):
    """⛔ RETAINED ONLY AS THE FALSIFICATION of the note above: the mass-identity rescale, which the
    policy no longer uses. `test_the_mass_identity_rescale_amplifies_a_weak_claim` is what keeps the
    reason measured rather than remembered."""
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
        # ⭐⭐⭐ A COMPOSITION CLAIM REQUIRES COMPOSITION EVIDENCE, and a structural lock is evidence
        # only where the STRUCTURE determines the composition — i.e. where no RNA strand is admissible
        # (``g1_locked``). A slot that admits an RNA strand and still reports certainty is reporting a
        # DEFAULT (``f_g = 1``, "all gDNA") rather than a measurement, and believed as a composition it
        # is catastrophic: transported into a 23,889-fragment exon at a ZERO-gDNA library it delivered a
        # gDNA share of 1.0000 and ψ solved ``f_g = 1.000`` against a truth of 0.000 (measured, and it
        # was 75,963 false-positive fragments on that condition alone).
        # ⛔ The RNA arms of such a slot therefore make NO CLAIM — value and precision zeroed in ONE
        # statement (``TRAPS: zero-the-precision-with-the-value``) — which also withdraws the SUPPLY the
        # composition licence needs, so the slot can still transport an ABUNDANCE and never a
        # composition. ⭐ Its gDNA claim SURVIVES: "no fragments over this much opportunity" is a real
        # observation and it is what the zero controls run on.
        _lock = np.asarray(getattr(own, "struct_lock", np.zeros(len(self._fp), bool)), bool)
        _mis_lock = _lock & (self._fp | self._fn)
        if np.any(_mis_lock):
            v_p, v_n = self._val0[1].copy(), self._val0[2].copy()
            p_p, p_n = self._prec0[1].copy(), self._prec0[2].copy()
            v_p[_mis_lock] = 0.0
            v_n[_mis_lock] = 0.0
            p_p[_mis_lock] = 0.0
            p_n[_mis_lock] = 0.0
            self._val0 = (self._val0[0], v_p, v_n)
            self._prec0 = (self._prec0[0], p_p, p_n)
        self._mis_lock = _mis_lock
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
        # the flux banks, per FACE and per TRANSCRIPT strand: the model-free abundance and the raw
        # count that gives it its precision
        self._inv_sj_lo = np.asarray(ctx.inv_sj_lo, np.float64)
        self._inv_sj_hi = np.asarray(ctx.inv_sj_hi, np.float64)
        _z2 = np.zeros_like(self._inv_sj_lo)
        self._sj_lo = (
            np.asarray(getattr(g, "sj_count_lo", _z2), np.float64) if g is not None else _z2
        )
        self._sj_hi = (
            np.asarray(getattr(g, "sj_count_hi", _z2), np.float64) if g is not None else _z2.copy()
        )

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
            # ⭐⭐⭐ THE SPLICE IN's CERTIFIED-RNA MEASUREMENT, per strand: the flux whose bodies lie on
            # the DESTINATION's side enters the exon as RNA that cannot be gDNA. Its abundance is the
            # model-free bank and its precision is its own COUNT's — it is a measurement, so it is
            # neither licensed nor rescaled.
            # ⚠ the SAME face pairing the ratio uses: travelling low→high the SOURCE presents its HIGH
            # face, so the flux entering the destination is the source's HIGH-end flux; backward is the
            # mirror. Getting this backwards makes the operator silently inert (measured: it did).
            f_in = self._inv_sj_lo if backward else self._inv_sj_hi
            c_in = self._sj_lo if backward else self._sj_hi
            self._tables[backward] = {
                "splice_in_rho": f_in[sn],
                "splice_in_prec": np.where(
                    c_in[sn] > 0.0, 1.0 / count_logvar(np.maximum(c_in[sn], 0.0)), 0.0
                ),
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
        si_rho = tab["splice_in_rho"].tolist()
        si_prec = tab["splice_in_prec"].tolist()
        fx_p = tab["fx_p"].tolist()
        fx_n = tab["fx_n"].tolist()
        dst_is_b = tab["dst_is_b"].tolist()
        E_g_l, E_r_l = self._E_g.tolist(), self._E_r.tolist()
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

        def _frame(p: float, lr: float) -> float:
            """The SPLICE-IN bound's frame variance — ``(log r)²``, the tree's own
            :func:`~.variance.splice_in_frame_logvar`, identically 0 where the two faces agree."""
            v = float(lr) * float(lr)
            return p / (1.0 + p * v) if (p > 0.0 and v > 0.0) else p

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
                # ⛔⛔ **THE FLUX IS NOT IN THE TRANSPORTED POPULATION, AND THIS IS A DELIBERATE
                # DEVIATION FROM `EQUATIONS.md` §3.5e(ii)'s WORKED ARITHMETIC — recorded, not slipped
                # in.** §3.5e adds the spliced-in flux to the message and rescales the whole thing; that
                # is self-consistent in one frame, which is the frame the worked example is in. Under
                # CAPTURE it is not: those fragments' BODIES lie in the destination exon, so they were
                # enriched by the EXON's probes, and multiplying them by a boundary→exon enrichment
                # ratio enriches them twice. ⭐ They enter once, below, as a MEASUREMENT of the exon's
                # RNA at their own count's precision (`DESIGN.md` §0c.1's four properties), which is
                # also what stops them being damped by a transport they never made.
                in_p = tp_
                in_n = tn_
                # ⭐⭐ the transport is the MEASURED enrichment ratio, raised to the knob: ``r^w``
                # interpolates in LOG SPACE between no rescale (w = 0) and believing the ratio in full
                # (w = 1). ⛔ NOT the mass-identity ``k`` — see the note on
                # :func:`composition_rescale_factor` for the 235,800x amplification that cost.
                kk = math.exp(w * lr)
                if True:
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
            # ⭐⭐⭐ THE SPLICE IN — a MEASUREMENT, not an imputation: a spliced fragment cannot be
            # gDNA, so the flux entering this exon is CERTIFIED RNA of its own strand, at its own
            # COUNT's precision. ⛔ It is added AFTER the transport AND after the disagreement damping,
            # because it did not travel as an imputation and must not inherit the hop's distrust — that
            # is what "it carries its own precision" means. ⛔ And BEFORE the admissibility refusal
            # below, so a refused arm's zeroing stays the last word
            # (``TRAPS: zero-the-precision-with-the-value``: an operator that GRANTS a precision must
            # run before the zeroing, never after).
            if not dst_is_b[i]:
                # ⛔⛔ **THE CERTIFIED FLUX IS A LOWER BOUND, NOT AN ESTIMATE, and treating it as an
                # equality is measurably worse.** It counts MATURE fragments only; the exon's RNA is
                # mature PLUS nascent, so an equality under-states every exon with unspliced RNA in it
                # and ψ books the shortfall as gDNA. Measured on the test chromosome: as an equality it
                # won the capture-ON zero controls (33,281 → 18,699) and LOST four other rows, worst
                # 666 → 4,172 at `g50 ss0.99 capture_off` — a condition with no probes at all.
                # ⭐ So it is applied ONE-SIDED: it may RAISE an RNA claim and may never lower one, and
                # a bound already satisfied says nothing and adds no precision — that is the whole
                # difference between a bound and a measurement.
                # ⭐⭐ And where it DOES bind it fuses by precision rather than replacing the claim: the
                # flux is measured in the sj's own frame and over-states an exon's contained density
                # (`EQUATIONS.md` §3.6b's two divisors of opposite sign), so jumping the claim to it
                # over-corrects — measured, replacement under-called gDNA by 20,752 fragments at
                # `g50 ss0.50 capture_on`. A fusion raises the claim as far as the two precisions
                # justify and no further.
                gp, gpp = si_rho[i][0], si_prec[i][0]
                gn, gpn = si_rho[i][1], si_prec[i][1]
                # ⚠ the bound is measured in the sj's OWN frame — its fragments span the junction and
                # sample a different position distribution from the exon's contained ones, which under
                # capture is a different enrichment. That frame step is charged as a variance, and it is
                # the tree's existing derivation for this very operator (`splice_in_frame_logvar`,
                # ``(log r)²``), so the bound is believed in proportion to how comparable the two
                # frames are rather than absolutely.
                if gpp > 0.0 and gp > tp_:
                    tp_, tpp = gp, tpp + _frame(gpp, lr)
                if gpn > 0.0 and gn > tn_:
                    tn_, tpn = gn, tpn + _frame(gpn, lr)
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
        in_p, in_n = vp, vn  # the flux enters as a MEASUREMENT below, never in the transport
        # the k-form: the source's claim laid down in the destination's geometry must account for the
        # destination's OBSERVED mass. Exact under the premise, free of the destination's BELIEF.
        ok = eq & (w > 0.0)
        k = np.where(ok, np.exp(w * log_r), 1.0)  # ``r^w`` — the MEASURED ratio, at the knob
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

        # ⭐ the SPLICE IN's certified-RNA measurement, vectorised — same table, same rule
        si_r = np.asarray(tab["splice_in_rho"], np.float64)
        si_p = np.asarray(tab["splice_in_prec"], np.float64)
        into_region = ~dst_is_b
        for arm in (0, 1):
            add_v = np.where(into_region, si_r[:, arm], 0.0)
            add_p = np.where(into_region, si_p[:, arm], 0.0)
            base_v, base_p = (out_p, pp) if arm == 0 else (out_n, pn)
            # ONE-SIDED: the bound may raise the claim, never lower it — see the scalar twin's note
            binds = (add_p > 0.0) & (add_v > base_v)
            fr = log_r * log_r  # the frame step between the sj's own frame and the exon's
            add_p = np.where(fr > 0.0, add_p / (1.0 + add_p * fr), add_p)
            v = np.where(binds, add_v, base_v)
            p = np.where(binds, base_p + add_p, base_p)
            if arm == 0:
                out_p, pp = v, p
            else:
                out_n, pn = v, p
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
