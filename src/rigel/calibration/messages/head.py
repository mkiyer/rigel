"""HeadPolicy — every operator the evolved solver carried, each behind a NAMED switch.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (per-array parity against the shipped answer)

⛔⛔ **THIS FILE IS A RESTRUCTURE, NOT A REWRITE, AND THE ARITHMETIC IS MOVED VERBATIM.** With every switch
ON it must reproduce the shipped panel BYTE-IDENTICALLY on all 1,872 scored fields of all 72 rows. That
identity is the whole value of the move: any difference is a bug, full stop. The alternative — a clean
rebuild — was tried, came out **+103 %**, and took two sessions to decompose into a correct derivation, a
UNIT ERROR and a STRUCTURAL disconnection, one of which no amount of derivation review would have caught.

⭐⭐ **WHY EVERY OPERATOR IS A SWITCH.** Removing them as a block is already measured: it is a net *harm* on
three of the four strata and a large regression on one. So the block figure says nothing about any member,
and TRAPS: all-small-singly-large-jointly says exactly what to do about it — *when every single ablation is small and the joint one
is large, go one stage upstream*. One switch per independently-ablatable operator is what lets
``ladder_arm_ab.py`` price them ONE AT A TIME, per stratum, which is the next step and not this one.

⚠ **The two twins are still twins.** ``_step`` (scalar, one slot at a time) and ``deliver`` (vectorised) run
the same transform in the same order, and the duplication is DELIBERATE: routing the sequential scan through
the numpy form was measured at **15.7×** per operation. They also differ in three edge cases ON PURPOSE, and
those differences are NOT bugs to unify — see :meth:`_HeadRelay.scan`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields

import numpy as np

from ..node_geometry import g1_locked, node_total_density, terminus_flank_gain
from ..node_init import own_composition_logvar
from . import NeighbourState, PsiMessage, StepContext
from .variance import (
    _fmax,
    composition_logvar,
    graft_frame_logvar,
    graft_frame_logvar_scalar,
    graft_premise_logvar,
    mismatch_deflate,
    mismatch_gap,
    peel_continue_share,
    peel_continue_share_scalar,
    peel_share_logvar,
    residual_level,
    residual_level_scalar,
    transfer_logvar,
)

__all__ = ["HeadPolicy", "HeadSwitches"]

_EPS = 1.0e-9


@dataclass(frozen=True, slots=True)
class HeadSwitches:
    """One boolean per independently-ablatable operator. **All default ON**, which is the shipped answer.

    ⛔ Turning one OFF is an EXPERIMENT, and the operator ledger predicts the sign for several of them: most
    are ~0 at the prior-free pass, and several are NEGATIVE on three of the four strata. Score PER STRATUM,
    never pooled — the panel total hides a sign flip between strata (TRAPS: never-pool-the-strata).
    """

    #: the composition reframe ``r = rho_tot(dst)/rho_tot(src)``. OFF ⇒ every message passes through at
    #: ``r = 1``. ⚠ Ledger: RE-DERIVE — keep the rule, re-derive the scope.
    reframe: bool = True
    #: the frame is a PAIR, one total per flank (``rho_lo``/``rho_hi``). OFF ⇒ one total on both sides,
    #: which inflates the lift on whichever side faces an intron.
    flank_pair: bool = True
    #: ⭐⭐ the gDNA SCALE RULE: a gDNA LEVEL crosses UNSCALED unless the source SUPPLIED both components
    #: of the pair. OFF ⇒ the gDNA arm is reframed like the RNA arms (``r_g = r``).
    gdna_level_scale: bool = True
    #: the POPULATION half of the composition licence — do the two objects measure the same RNA
    #: population? OFF ⇒ granted everywhere.
    terminus_population: bool = True
    #: the relay's mass pin ``Sigma_c rho_c E_c = M``, licensed in exactly two states. ⚠ Ledger: the
    #: ceiling says deleting it outright cost the panel **+0.0002** — it landed on the derivation.
    mass_pin: bool = True
    #: GRAFT: a flanking EDGE's measured junction density joins the RNA claim entering an EXON.
    graft: bool = True
    #: the-graft-frame-variance — the graft's FRAME-MISLIFT variance. Identically 0 at ``r = 1``.
    graft_frame_var: bool = True
    #: P1d — the graft's PREMISE variance, ONE library-level scalar fitted from flanking seam pairs.
    #: ⚠ Ledger: CUT from the backbone — re-derive per class when TSS/TES land.
    graft_premise_var: bool = True
    #: PEEL: an EXON's RNA leaving into an EDGE is scaled by the CONTINUING share ``w`` (the-continuing-share).
    peel: bool = True
    #: sigma^2_transfer = Var(log r) (the-reframe-scale-variance) — the reframe's own scale sampling, 0 on the matched graft.
    transfer_var: bool = True
    #: the DerSimonian-Laird composition-mismatch deflation (the-mismatch-deflation). Silent where the destination has no
    #: self-solve, i.e. never on unstranded data.
    mismatch_var: bool = True
    #: P1e — the conservation SURPRISE as a DL damping term. ⚠⚠ Its own comment: *"PARTLY A DEBT — THIS
    #: PRICES A BIAS AS A VARIANCE … a variance cannot move a mode toward truth."* Ledger: CUT.
    conservation_var: bool = True
    #: the lambda-EMISSION gate — a source carrying ONE component has no SPLIT to claim.
    lam_emission_gate: bool = True
    #: the four psi channels, individually. ⭐ ``rna_channel`` OFF is the certified-RNA measurement, which
    #: carries 75 % of the posterior precision on the confidently-wrong unstranded x capture-OFF exons.
    gdna_channel: bool = True
    rna_channel: bool = True
    lam_channel: bool = True
    theta_channel: bool = True

    def names(self) -> tuple[str, ...]:
        return tuple(f.name for f in fields(self))

    def off(self) -> tuple[str, ...]:
        """Which switches are down — the arm's own label, and its TRAPS: could-the-arm-have-fired check."""
        return tuple(f.name for f in fields(self) if not getattr(self, f.name))


class HeadPolicy:
    """The shipped message layer. ``HeadPolicy()`` — every switch ON — is the shipped answer."""

    name = "head"

    def __init__(self, switches: HeadSwitches | None = None):
        self.switches = switches if switches is not None else HeadSwitches()
        if self.switches.off():
            self.name = "head-" + "-".join(f"no_{s}" for s in self.switches.off())

    def prepare(self, ctx: StepContext) -> _HeadRelay:
        return _HeadRelay(ctx, self.switches)


class _HeadRelay:
    """One sweep's worth of prepared state. Everything derived here is a POLICY derivation."""

    def __init__(self, ctx: StepContext, sw: HeadSwitches):
        self.ctx = ctx
        self.sw = sw
        cap = ctx.capture
        n = ctx.n_slots

        # ── the arrays this policy works in, named as the shipped solver named them ───────────────────
        is_bnd_a = np.asarray(ctx.is_edge, bool)
        ex_a = np.asarray(ctx.is_exon_node, bool)
        fp_a, fn_a = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
        M = np.asarray(ctx.mass, np.float64)
        E_g = np.asarray(ctx.eff_gdna_global, np.float64)
        E_r = np.asarray(ctx.eff_rna, np.float64)
        SPL = np.asarray(ctx.junction_count, np.float64)
        ESP = np.asarray(ctx.eff_junction, np.float64)
        _n_node = np.asarray(ctx.n_slot, np.float64)
        ni = ctx.own
        self._is_bnd_a, self._ex_a, self._fp_a, self._fn_a = is_bnd_a, ex_a, fp_a, fn_a
        self._M, self._E_g, self._E_r, self._SPL, self._n_node = M, E_g, E_r, SPL, _n_node
        self._is_amb = fp_a & fn_a  # AMBIG: both strands live ⇒ the tilt is a free DOF

        # ── the-reframe-scale-variance sigma^2_transfer = Var(log r): the per-slot total-density log-variance ──────────────────
        # Var(log rho_tot) = trigamma-count + [(1/E_g - 1/E_r)/B]^2 * Var(f_g), with Var(f_g) =
        # (f_g(1-f_g))^2/tau_lam CAPPED at f_g(1-f_g) and 0 for a composition-certain slot. Evaluated at
        # the INPUT belief, consistent with the frame pair, which shares exactly this belief.
        # ⛔ TRAPS: a-variance-cap-asserts-certainty: the cap ``f_g(1-f_g)`` is 0 at the ``f_g = 1`` default of an evidence-free
        # slot, so the composition half vanishes there and sigma^2_transfer is the counting term alone.
        # That is a PROVEN defect carrying a STRICT xfail; it is preserved here because its fix is
        # panel-negative alone and this commit changes nothing.
        _fg0 = np.clip(np.asarray(ctx.belief_fg, np.float64), 0.0, 1.0)
        _fgfr = _fg0 * (1.0 - _fg0)
        _tau = np.asarray(ni.tau_lam, np.float64)
        _var_fg = np.where(
            np.asarray(ni.struct_lock, bool),
            0.0,
            np.where(_tau > _EPS, np.minimum(_fgfr * _fgfr / np.maximum(_tau, _EPS), _fgfr), _fgfr),
        )
        logvar_tot = np.asarray(composition_logvar(_fg0, E_g, E_r, _var_fg, _n_node), np.float64)
        if not sw.transfer_var:
            logvar_tot = np.zeros_like(logvar_tot)
        self._logvar_tot = logvar_tot

        # ── the per-TRANSCRIPT-STRAND junction DENSITY at each line ────────────────────────────────────
        def _mature_rho(strand: int) -> np.ndarray:
            c, e = SPL[:, strand], ESP[:, strand]
            live = (c > _EPS) & (e > _EPS)
            return np.where(live, c / np.where(live, e, 1.0), 0.0)

        spl_p = _mature_rho(0)  # + transcript junction density at this line (0 on NODE slots)
        spl_n = _mature_rho(1)
        self._spl_p, self._spl_n = spl_p, spl_n

        _li_a = np.asarray(ctx.left, np.int64)
        _ri_a = np.asarray(ctx.right, np.int64)
        _vl_a, _vr_a = _li_a >= 0, _ri_a >= 0
        _sl_a, _sr_a = np.clip(_li_a, 0, n - 1), np.clip(_ri_a, 0, n - 1)
        self._vl_a, self._vr_a, self._sl_a, self._sr_a = _vl_a, _vr_a, _sl_a, _sr_a

        # ── ⭐⭐⭐ THE POPULATION HALF OF THE COMPOSITION LICENCE ───────────────────────────────────────
        # An EDGE counts what spans it CONTIGUOUSLY, so T(EDGE) = T(left) ∩ T(right), and a transcript
        # TERMINUS at the EDGE makes one flank's RNA population strictly larger. Where two objects do not
        # measure the same population, a density discrepancy is not attributable to capture enrichment —
        # the two are indistinguishable — so the composition may not be imputed and the level crosses
        # unscaled. ⛔ TERMINI ONLY: a DONOR/ACCEPTOR EDGE also changes the population, but there the flux
        # is MEASURED and the graft and the peel exist to route it.
        # THE PAIR ALGEBRA: the chain strictly alternates, so exactly one slot of an adjacent pair is an
        # EDGE and the pair ``(i, left[i])`` IS the pair ``(left[i], right[left[i]])`` — so one array
        # answers it for the left-hand step of every slot and the other is that array read through
        # ``right``.
        if sw.terminus_population:
            _rgain, _lgain = terminus_flank_gain(ctx.edge_flags)
            _pop_l_a = np.where(is_bnd_a, ~_lgain, ~_rgain[_sl_a]) & _vl_a
            _pop_r_a = np.where(_vr_a, _pop_l_a[_sr_a], False)
        else:
            _pop_l_a, _pop_r_a = _vl_a.copy(), _vr_a.copy()
        self._pop_l_a, self._pop_r_a = _pop_l_a, _pop_r_a

        # ── own per-component densities + precisions — the message-free SELF-SOLVE ─────────────────────
        og, op, on = ni.rho_g, ni.rho_pos, ni.rho_neg
        pg_own, pp_own, pn_own = ni.prec_g, ni.prec_pos, ni.prec_neg
        self._og, self._op, self._on = og, op, on

        # ── THREE-STREAM precision seeds (the single-lambda combine) ───────────────────────────────────
        #   * tau (COMPOSITION): the Schur lambda-precision; 0 at anchors + unstranded non-factory slots.
        #   * mg/mp/mn (MEASUREMENT): the anchor gDNA count (own, struct_lock only) + the spliced RNA
        #     count (added at the graft). A density-level bound, NOT a composition vote.
        _struct = np.asarray(ni.struct_lock, bool)
        tau_own = np.asarray(ni.tau_lam, np.float64)
        mg_own = np.where(_struct, np.asarray(pg_own, np.float64), 0.0)
        mp_own = np.zeros_like(np.asarray(pp_own, np.float64))
        mn_own = np.zeros_like(np.asarray(pn_own, np.float64))

        v_own_g, v_own_r = own_composition_logvar(ni.f_g, tau_own, _struct)
        # the same three states on the lambda axis, for the single-DOF composition stream: Var(lambda) =
        # 1/tau_lam — the two per-component arms are perfectly anti-correlated, so their Jacobians sum to
        # 1 and the count cancels.
        v_own_lam = np.where(
            _struct, 0.0, np.where(tau_own > _EPS, 1.0 / np.maximum(tau_own, _EPS), np.inf)
        )
        self._v_own_g, self._v_own_r, self._v_own_lam = v_own_g, v_own_r, v_own_lam

        # ── v_mu uses the spliced COUNT, never a mass ──────────────────────────────────────────────────
        # The accumulator deposits fragments fractionally, so at a junction face the median count is 33
        # against a median mass of 11 and the mass would over-state v_mu ~3x. With one count per line the
        # rule is structural rather than a discipline: ``junction_count`` IS the junction fragment count.
        self._mu_s = (spl_p, spl_n)
        self._v_mu_s = tuple(
            np.where(c > 0.0, 1.0 / np.maximum(c, _EPS), np.inf) for c in (SPL[:, 0], SPL[:, 1])
        )

        # ── ⭐⭐⭐ THE REFRAME FRAME IS A PAIR, ONE TOTAL PER FLANK ──────────────────────────────────────
        # ``rho_lo[k]`` is the total slot ``k`` presents to its genomic-LOW neighbour and ``rho_hi[k]`` the
        # one it presents to its genomic-HIGH neighbour; they differ only at an EDGE, and only by which
        # junctions' flux belongs on which side.
        # ⭐⭐ THE PAIRING RULE: a hop always joins two ADJACENT slots ``(k, k+1)``, and whichever is the
        # source the pair is the SAME pair, so ``r = rho_lo[k+1] / rho_hi[k]``. ⛔ DIRECTION DOES NOT
        # ENTER. It is NOT expressible as one array per direction: within ONE forward pass an EDGE at a
        # junction's low end is the DESTINATION of a hop from its low flank (junction flux INCLUDED) and
        # the SOURCE of the next hop into its high flank (EXCLUDED). Two arrays indexed by ROLE is what
        # expresses that.
        # ⚠ TRAPS: a-message-from-the-destinations-belief DEBT, NAMED: this reads the belief at BOTH ends, so the frame at a hop is a function of the
        # destination's belief. The ledger prices it — a *solved* belief sets the frame at slots carrying
        # 57-77 % of library mass.
        rho_lo, rho_hi = node_total_density(ctx.geometry, np.asarray(ctx.belief_fg, np.float64))
        if not sw.flank_pair:
            _one = np.maximum(rho_lo, rho_hi)
            rho_lo, rho_hi = _one, _one
        self._rho_lo, self._rho_hi = rho_lo, rho_hi

        # ── the graft's PREMISE variance: ONE library-level scalar, fitted by method of moments ─────────
        vgp_prem, vgn_prem = self._seam_pair(rho_lo, rho_hi)
        if not sw.graft_premise_var:
            vgp_prem, vgn_prem = np.zeros_like(vgp_prem), np.zeros_like(vgn_prem)
        self._vgp_prem, self._vgn_prem = vgp_prem, vgn_prem

        # ── ⭐ THE SEQUENTIAL SCAN'S OPERANDS, AS PYTHON LISTS ──────────────────────────────────────────
        # The scalar step reads every one of these ONE ELEMENT AT A TIME — ~6 M edge-iterations at genome
        # scale, ~40 reads each. ``.tolist()`` is exact on float64/int64/bool (identical IEEE-754 doubles),
        # so this is BIT-IDENTICAL by construction, and it buys ~3x: ``lst[i]`` costs a third of
        # ``arr[i]`` and yields a PYTHON float, whose arithmetic is ~3x ``np.float64``'s because no
        # intermediate is boxed in a 0-d array.
        (
            self._og_l,
            self._op_l,
            self._on_l,
            self._pg_own_l,
            self._pp_own_l,
            self._pn_own_l,
            self._mg_own_l,
            self._mp_own_l,
            self._mn_own_l,
            self._tau_own_l,
            self._M_l,
            self._E_g_l,
            self._E_r_l,
            self._n_node_l,
            self._logvar_l,
        ) = (
            np.asarray(a, np.float64).tolist()
            for a in (
                og,
                op,
                on,
                pg_own,
                pp_own,
                pn_own,
                mg_own,
                mp_own,
                mn_own,
                tau_own,
                M,
                E_g,
                E_r,
                _n_node,
                logvar_tot,
            )
        )
        self._ex_l, self._bnd_l, self._fp_l, self._fn_l = (
            a.tolist() for a in (ex_a, is_bnd_a, fp_a, fn_a)
        )
        # the destination's own composition CERTAINTY — case (ii) of the mass pin's licence. Both axes: an
        # ``intergenic|exon`` EDGE is as structurally pure-gDNA as an intergenic NODE, and it is gated on
        # ``g1_locked`` rather than on ``node_init``'s node-only ``struct_lock`` because the object in
        # question is an EDGE.
        self._g1_l = g1_locked(fp_a, fn_a).tolist()
        self._spl_p_l, self._spl_n_l = spl_p.tolist(), spl_n.tolist()
        self._SP_l, self._SN_l = SPL[:, 0].tolist(), SPL[:, 1].tolist()
        self._mu_l = [c.tolist() for c in self._mu_s]
        self._v_mu_l = [c.tolist() for c in self._v_mu_s]
        self._vgp_l, self._vgn_l = vgp_prem.tolist(), vgn_prem.tolist()
        self._rho_lo_l, self._rho_hi_l = rho_lo.tolist(), rho_hi.tolist()
        self._pop_l_l, self._pop_r_l = _pop_l_a.tolist(), _pop_r_a.tolist()

        # ── the static half of the diagnostics capture ──────────────────────────────────────────────────
        if cap is not None:
            cap.setdefault("_uni_static", {}).update(
                M=M,
                E_g=E_g,
                E_r=E_r,
                spl_p=spl_p,
                spl_n=spl_n,
                og=og,
                op=op,
                on=on,
                pg_own=pg_own,
                pp_own=pp_own,
                pn_own=pn_own,
                rho_lo=rho_lo,
                rho_hi=rho_hi,
                logvar_tot=logvar_tot.copy(),
                mature_pos=SPL[:, 0].copy(),
                mature_neg=SPL[:, 1].copy(),
                tau_own=tau_own.copy(),
                mg_own=mg_own.copy(),
                struct_lock=_struct.copy(),
                is_bnd=is_bnd_a,
                is_exon=ex_a,
                fp=fp_a,
                fn=fn_a,
            )

    # ── the graft's premise variance ───────────────────────────────────────────────────────────────────
    def _flank_dom(self, rho_lo_a, rho_hi_a, spf):
        """Per slot: the flux each of its two flanking EDGES sends it, ALREADY lifted into this slot's
        frame, per strand.

        ⭐ The lift is the same reframe ``r`` the scan uses, so it takes the same FLANK PAIR and the same
        role pairing: for the LEFT flank the neighbour is the genomic-LOW slot, so this slot presents
        ``rho_lo`` and the neighbour ``rho_hi``; the RIGHT flank is the mirror. Reading one
        junction-inclusive total on both sides lifted each flux by a ratio inflated on exactly the side
        facing an intron."""
        is_bnd_a, vl, vr, sl, sr = self._is_bnd_a, self._vl_a, self._vr_a, self._sl_a, self._sr_a

        def _side(ok, nb, rho_here, rho_there):
            r = np.where(
                ok & (rho_there[nb] > _EPS) & (rho_here > _EPS),
                rho_here / np.maximum(rho_there[nb], _EPS),
                1.0,
            )
            return np.where(ok & is_bnd_a[nb], spf[nb] * r, 0.0)

        return (
            _side(vl, sl, rho_lo_a, rho_hi_a),
            _side(vr, sr, rho_hi_a, rho_lo_a),
        )

    def _seam_pair(self, rho_lo_a, rho_hi_a):
        """Per strand: the graft's premise log-variance — ONE library-level scalar, fitted by method of
        moments from the destination-frame disagreement of exons' flanking seam PAIRS, and applied to every
        graft edge.

        ⚠⚠ **A DEBT, not a model.** The one scalar stands in for a quantity that splits >=30x on whether
        the boundary carries a transcript TERMINUS — a bit the region map does not have.

        ⭐ The POOLED fit is applied to EVERY graft edge and the per-edge value is NOT used. ``d^2`` from
        ONE pair is a single draw of a scaled chi^2_1, so its own coefficient of variation is sqrt(2) — a
        per-edge "measurement" of a variance is mostly noise, and the UNDER-charging half does the damage
        because it REPLACES the population value on the ~48 % of edges where it fires. It also removes a
        real BP objection: with the per-edge form the message from the LEFT seam carried a variance
        computed from the RIGHT seam's counts, so a non-adjacent slot's data reached the destination
        twice."""
        cap = self.ctx.capture
        vl, vr, sl, sr = self._vl_a, self._vr_a, self._sl_a, self._sr_a
        logvar_tot = self._logvar_tot
        out = []
        for spf, vmu in ((self._spl_p, 0), (self._spl_n, 1)):
            fl, fr = self._flank_dom(rho_lo_a, rho_hi_a, spf)
            # each seam's own noise: its spliced COUNT (never the mass) ⊕ its lift's scale sampling (the-reframe-scale-variance's
            # source leg; the destination's leg is common to both lifts and cancels in ``d``).
            _lv = np.where(np.isfinite(logvar_tot), logvar_tot, 0.0)
            per, pooled = graft_premise_logvar(
                fl,
                fr,
                np.where(vl, self._v_mu_s[vmu][sl] + _lv[sl], np.inf),
                np.where(vr, self._v_mu_s[vmu][sr] + _lv[sr], np.inf),
            )
            if cap is not None:  # inert: the fitted scalar and the population it was fitted on
                _ok = (fl > _EPS) & (fr > _EPS)
                _d = np.log(np.maximum(fl, _EPS)) - np.log(np.maximum(fr, _EPS))
                _vv = np.where(vl, self._v_mu_s[vmu][sl] + _lv[sl], 0.0) + np.where(
                    vr, self._v_mu_s[vmu][sr] + _lv[sr], 0.0
                )
                cap.setdefault("_glv", []).append(
                    {
                        "strand": vmu,
                        "omega": pooled,
                        "n_pairs": int(_ok.sum()),
                        "ok": _ok.copy(),
                        "d": _d.copy(),
                        "noise": _vv.copy(),
                        "Ed2": float((_d[_ok] ** 2).mean()) if _ok.any() else 0.0,
                        "Enoise": float(_vv[_ok].mean()) if _ok.any() else 0.0,
                    }
                )
            out.append(np.full_like(per, pooled))
        return out[0], out[1]

    # ── the peel share, both twins ─────────────────────────────────────────────────────────────────────
    def _peel_share(self, tg, tpg, tp, tn):
        """The continuing share ``w`` and ``Var(log w)`` per strand, at every slot, for a message whose
        gDNA claim is ``(tg, tpg)`` and whose RNA claim is ``(tp, tn)``. Returns
        ``((w_p, vw_p), (w_n, vw_n))``; ``Var(log w) = +inf`` (⇒ zero precision, an inert message) only
        where the level does not exist.

        A boundary's unspliced crossing is ``gDNA + the RNA that CONTINUES``, so the honest operator is a
        SCALING by the continuing share ``w = rho_nu/(rho_nu+rho_mu)`` and not the subtraction it replaced.
        ``w`` is enrichment-free — capture multiplies the continuing and the splicing channels alike, so
        ``e`` cancels identically inside the ratio — and a scaling COMMUTES with the scale error the
        reframe carries where a difference AMPLIFIES it by ``u = 1/w``.

        ⚠ TWIN of :meth:`_peel_share_scalar` — mirror any change into both."""
        cap = self.ctx.capture
        M, _n_node, E_g, E_r = self._M, self._n_node, self._E_g, self._E_r
        _vg = np.where(np.asarray(tpg, np.float64) > 0.0, 1.0 / np.maximum(tpg, _EPS), np.inf)
        _nu_m, _vlog_m, _vl_m = residual_level(M, _n_node, tg, E_g, E_r, _vg)
        _A = np.asarray(tp, np.float64) + np.asarray(tn, np.float64)
        _a_p = np.where(_A > _EPS, np.asarray(tp, np.float64) / np.maximum(_A, _EPS), 0.0)
        out = []
        for _a, _mu, _vmu in (
            (_a_p, self._mu_s[0], self._v_mu_s[0]),
            (1.0 - _a_p, self._mu_s[1], self._v_mu_s[1]),
        ):
            # ── THE FUSE, in LINEAR density space ──────────────────────────────────────────────────────
            # Linear is the coordinate that lets a confident near-zero claim actually pull the level down;
            # a geometric fuse of positive modes cannot reach zero, which is why "the intron sets the
            # level" kept failing to. the-residual-level's level is a lambda-axis TOTAL, split onto this strand by the
            # message's own tilt. Mask the inf BEFORE the product — the standing ``0*inf = nan`` trap,
            # ``np.where`` evaluating both branches.
            _nu_ms = _nu_m * _a
            _vl_s = np.where(np.isfinite(_vl_m), _vl_m, 0.0) * _a * _a
            _pm = np.where(np.isfinite(_vl_m) & (_nu_ms > _EPS), 1.0 / np.maximum(_vl_s, _EPS), 0.0)
            _pt = _pm
            _live = _pt > _EPS
            _nu = np.where(_live, (_pm * _nu_ms) / np.maximum(_pt, _EPS), 0.0)
            # the-residual-level already returned this level's log-variance; use it. The round trip it replaces
            # (``k = rho^2/V`` then ``psi'(k)``) is EXACTLY ``psi'(1/v_log)``, and ``residual_level``'s
            # ``k >= 1`` floor is an exact limit of the TRUNCATION — but here ``k`` is a reciprocal
            # variance and not a count, so the same floor became a hard CEILING of ``psi'(1) = pi^2/6``,
            # over-stating confidence 6x on exactly the seams where the level is least determined.
            _v_nu = np.where(_live, _vlog_m, np.inf)
            _w = np.where(_live, peel_continue_share(_nu, _mu), 0.0)
            if cap is not None:
                cap.setdefault("_lvl", []).append(  # inert: the level's provenance, per slot
                    {
                        "nu": np.asarray(_nu).copy(),
                        "v_nu": np.asarray(_v_nu).copy(),
                        "pm": np.asarray(_pm).copy(),
                        "nu_m": np.asarray(_nu_ms).copy(),
                        "mu": np.asarray(_mu).copy(),
                        "w": np.asarray(_w).copy(),
                        "v_g": np.asarray(_vg).copy(),
                        "vl_m": np.where(np.isfinite(_vl_m), _vl_s, np.inf),
                        "phi": np.asarray(tg * E_g / np.maximum(M, _EPS)).copy(),
                    }
                )
            # a spliced DENSITY with no spliced COUNT cannot be priced ⇒ no claim at all.
            _ok = _live & (np.isfinite(_vmu) | ~(_mu > _EPS))
            out.append(
                (
                    _w,
                    np.where(
                        _ok,
                        peel_share_logvar(
                            1.0 - _w,
                            np.where(_live, _v_nu, 0.0),
                            np.where(np.isfinite(_vmu), _vmu, 0.0),
                        ),
                        np.inf,
                    ),
                )
            )
        return out

    def _peel_share_scalar(self, i, tg, tpg, tp, tn):
        """The SCALAR twin of :meth:`_peel_share`, for one slot ``i`` — see that docstring for the model.

        ⚠ TWIN: mirror any change into both. It exists because the scan is sequential and calls this once
        per slot per direction, so the array form runs ~50 numpy ops on 0-d arrays per call — 0.5-0.7 us
        each against ~0.02 us for the float expression. Same arithmetic in the same association order;
        every ``np.where`` becomes the branch it encodes, so the dead arms are no longer evaluated and
        discarded. Measured 25x on this path, which is the bulk of a genome-scale calibration."""
        _vg = 1.0 / _fmax(tpg, _EPS) if tpg > 0.0 else math.inf
        _nu_m, _vlog_m, _vl_m = residual_level_scalar(
            self._M_l[i], self._n_node_l[i], tg, self._E_g_l[i], self._E_r_l[i], _vg
        )
        _fin = math.isfinite(_vl_m)
        _A = tp + tn
        _a_p = tp / _A if _A > _EPS else 0.0
        out = []
        for _a, _mu, _vmu in (
            (_a_p, self._mu_l[0][i], self._v_mu_l[0][i]),
            (1.0 - _a_p, self._mu_l[1][i], self._v_mu_l[1][i]),
        ):
            _nu_ms = _nu_m * _a
            _pm = (
                1.0 / _fmax((_vl_m if _fin else 0.0) * _a * _a, _EPS)
                if (_fin and _nu_ms > _EPS)
                else 0.0
            )
            if not _pm > _EPS:  # no level ⇒ no claim: w = 0 at zero precision, an inert message
                out.append((0.0, math.inf))
                continue
            # ⚠ NOT ``_nu_ms``: ``(a*b)/a != b`` in floating point. This is the one-estimator remnant of a
            # three-estimator fuse (the two illegal arms were deleted); kept verbatim so the delta stays
            # attributable.
            _nu = (_pm * _nu_ms) / _pm
            _v_nu = _vlog_m  # the-residual-level's own log-variance — see the twin
            _w = peel_continue_share_scalar(_nu, _mu)
            _wm = 1.0 - _w
            _ok = math.isfinite(_vmu) or not _mu > _EPS
            out.append(
                (
                    _w,
                    _wm * _wm * (_v_nu + (_vmu if math.isfinite(_vmu) else 0.0))
                    if _ok
                    else math.inf,
                )
            )
        return out

    # ── THE SCAN: one hop, accumulated in place — the forward half of forward-backward ─────────────────
    def scan(self, *, backward: bool):
        """Return ``(step, publish)``. ``step(s, i)`` relays from source ``s`` into destination ``i``.

        ⭐ **The whole direction dependence is one swap.** The FORWARD pass reads its source at ``i-1``,
        the slot's genomic-LOW neighbour, so the destination presents its LOW-flank total and the source
        its HIGH-flank one; the backward pass is the mirror. Nothing else differs.

        ⚠ **This is the SCALAR twin of :meth:`deliver`, and the two differ deliberately in three edge
        cases which are NOT bugs to unify:** (1) the scan skips invalid edges with an early return while
        the combine masks them to ``r = 0``; (2) ``_damp`` here uses the raw ``p`` where the combine uses
        ``max(p, _EPS)``, which differ only for ``0 < p < _EPS``; (3) the scan short-circuits the graft
        block under ``if _gr`` while the combine evaluates the frame variance on every edge and masks.
        """
        sw = self.sw
        # ⚠ which NEIGHBOUR array the scan walks is the BACKBONE's business, not this policy's — all the
        # direction means here is which flank total each side of a hop presents.
        rho_dst = self._rho_hi_l if backward else self._rho_lo_l
        rho_src_a = self._rho_lo_l if backward else self._rho_hi_l
        pop = self._pop_r_l if backward else self._pop_l_l

        # every operand below is a Python float or bool — see the ``.tolist()`` block in __init__
        og_l, op_l, on_l = self._og_l, self._op_l, self._on_l
        pg_own_l, pp_own_l, pn_own_l = self._pg_own_l, self._pp_own_l, self._pn_own_l
        mg_own_l, mp_own_l, mn_own_l = self._mg_own_l, self._mp_own_l, self._mn_own_l
        tau_own_l = self._tau_own_l
        M_l, E_g_l, E_r_l, logvar_l = self._M_l, self._E_g_l, self._E_r_l, self._logvar_l
        ex_l, bnd_l, fp_l, fn_l, g1_l = self._ex_l, self._bnd_l, self._fp_l, self._fn_l, self._g1_l
        spl_p_l, spl_n_l, SP_l, SN_l = self._spl_p_l, self._spl_n_l, self._SP_l, self._SN_l
        vgp_l, vgn_l = self._vgp_l, self._vgn_l
        peel_share_scalar = self._peel_share_scalar

        # the running state: each slot's context belief IN ITS OWN FRAME
        rg, rp, rn = og_l.copy(), op_l.copy(), on_l.copy()
        pg, pp, pn = pg_own_l.copy(), pp_own_l.copy(), pn_own_l.copy()  # full → MODE fusion
        mg, mp, mn = mg_own_l.copy(), mp_own_l.copy(), mn_own_l.copy()  # MEASUREMENT
        tau = tau_own_l.copy()  # COMPOSITION (tau_lambda) → the lambda-message

        def _damp(p, s2t):
            """sigma^2_transfer per-hop damping: 1/p → 1/p + sigma^2_transfer."""
            return 1.0 / (1.0 / p + s2t) if p > 0.0 else 0.0

        def _damp_v(p, v):
            """Add a log-variance to a precision. 1/p → 1/p + v ; v = inf ⇒ 0 (no claim)."""
            return p / (1.0 + p * v) if (p > 0.0 and math.isfinite(v)) else 0.0

        def _fuse(a, pa, b, pb):  # scalar precision-weighted density fuse
            p = pa + pb
            return ((pa * a + pb * b) / p, p) if p > _EPS else (a, 0.0)

        def step(s, i):
            # ⭐ ``rho_dst`` is read at the destination and ``rho_src_a`` at the source, and the DIRECTION
            # pairs them — each slot presents the total belonging to the flank it is compared against.
            rho_src = rho_src_a[s]
            rho_dst_i = rho_dst[i]
            if sw.reframe:
                r = (rho_dst_i / rho_src) if (rho_src > _EPS and rho_dst_i > _EPS) else 1.0
            else:
                r = 1.0  # no frame ⇒ pass-through
            # GRAFT (EDGE → EXON): the EDGE's measured junction flux is a density AT THE SOURCE, so it
            # joins the source's RNA BEFORE the reframe; the peel is measured at the destination and so is
            # applied after. Only an EXON receives the graft — an intron carries no junction flux.
            _gr = sw.graft and ex_l[i] and bnd_l[s]
            # sigma^2_transfer = Var(log r) (the-reframe-scale-variance): 0 on the matched-set GRAFT (r is common-mode across
            # {g,R} and cancels in the composition — charging it there is a double-count), Var(log r)
            # elsewhere. The COMPOSITION-mismatch term is the combine's job: the scan has no destination
            # self-solve to measure a gap against, and its running belief is already fused with the
            # messages, so a gap here would be feedback and not evidence.
            s2t = 0.0 if _gr else (logvar_l[i] + logvar_l[s])
            gp = spl_p_l[s] if _gr else 0.0
            gn = spl_n_l[s] if _gr else 0.0
            # ⭐⭐ THE gDNA SCALE. ``_lend``: may this source lend a COMPOSITION? Two conditions, and they
            # are different questions about the same step — SUPPLY (did the source state both components
            # of the pair from its OWN crossing population? "supplied" is a statement about PRECISION,
            # not about a density's value) and POPULATION (is the source measuring the same RNA population
            # as the destination?). Where either fails the gDNA LEVEL crosses UNSCALED and the pin is off
            # with it.
            _lend = pop[i] and pg[s] > 0.0 and (pp[s] + pn[s]) > 0.0
            r_g = r if (_lend or not sw.gdna_level_scale) else 1.0
            tg, tp, tn = rg[s] * r_g, (rp[s] + gp) * r, (rn[s] + gn) * r
            tpg, tpp, tpn = _damp(pg[s], s2t), _damp(pp[s], s2t), _damp(pn[s], s2t)  # full (mode)
            tmg, tmp, tmn = _damp(mg[s], s2t), _damp(mp[s], s2t), _damp(mn[s], s2t)  # measurement
            ttau = _damp(tau[s], s2t)  # composition
            # The grafted junction flux is a MEASUREMENT (a COUNT), not an imputation, so it carries its
            # own precision and is NOT tau-gated — the source's PREDICTION precision is 0 on unstranded
            # data and would otherwise drop the graft on the floor. It enters BOTH the mode fusion and the
            # MEASUREMENT stream, never the composition tau: a count is not a composition vote.
            if _gr:
                # the-graft-frame-variance: the grafted spliced density is measured in the DESTINATION exon's frame, so it has no
                # matched gDNA partner to cancel ``r`` against and the-reframe-scale-variance's graft-zero does not cover it.
                # Charge the frame step it is implicitly mis-lifted by. Identically 0 at r = 1.
                _s2f = s2t + (graft_frame_logvar_scalar(r) if sw.graft_frame_var else 0.0)
                _sps = SP_l[s]
                _spc = _sps / (1.0 + _sps * _s2f) if _sps > _EPS else 0.0
                _sns = SN_l[s]
                _snc = _sns / (1.0 + _sns * _s2f) if _sns > _EPS else 0.0
                tpp += _spc
                tpn += _snc
                tmp += _spc
                tmn += _snc
                _vgp, _vgn = vgp_l[i], vgn_l[i]
                tpp, tmp = _damp_v(tpp, _vgp), _damp_v(tmp, _vgp)
                tpn, tmn = _damp_v(tpn, _vgn), _damp_v(tmn, _vgn)

            if sw.peel and bnd_l[i] and ex_l[s]:  # EXON → EDGE: PEEL by COMPOSITION
                (_wp, _vwp), (_wn, _vwn) = peel_share_scalar(i, tg, tpg, tp, tn)
                tp, tn = tp * _wp, tn * _wn
                tpp, tmp = _damp_v(tpp, _vwp), _damp_v(tmp, _vwp)
                tpn, tmn = _damp_v(tpn, _vwn), _damp_v(tmn, _vwn)
            if not fp_l[i]:
                tp, tpp, tmp = 0.0, 0.0, 0.0
            if not fn_l[i]:
                tn, tpn, tmn = 0.0, 0.0, 0.0
            # ── PIN THE CONTEXT TO THIS SLOT'S OBSERVED MASS, WHERE NO BELIEF ENTERS THE BUDGET ────────
            # ``Sigma_c rho_c E_c = M`` is an IDENTITY under the imputation premise, not an approximation.
            # The pin restores it with a common factor ``k = M/S``, and the budget ``S`` fills every
            # component the context does NOT supply from the destination's own density — which is what
            # keeps a partial claim partial, and also what can make the delivered value a function of the
            # destination's own BELIEF. So it is licensed in exactly the two states where no belief reaches
            # ``S``:
            #   (i)  ``_lend`` — the context SUPPLIED the composition, so the premise is granted and the
            #        identity is implied. This is the reframe's own predicate: one licence, one place.
            #   (ii) ``g1_l[i]`` — the destination is a structurally pure-gDNA object, so there is no
            #        unsupplied component to fill in. ``f_g = 1`` there is STRUCTURE, not a belief, and
            #        ``S = rho_g E_g`` makes the pin hand the object its OWN MEASURED density. ⭐ That is
            #        the operator the capture landscape travels on.
            # ⛔ ANYWHERE ELSE IT WAS TRAPS: a-message-from-the-destinations-belief AT FULL STRENGTH: the closed form ``k = 1/(phi_msg + R_own)`` is a
            # saturating map whose fixed point is ``(1-R_own) rho_tot`` with ``R_own`` the RNA share of the
            # destination's OWN self-solve — EXACTLY 1/2 at a slot with no composition evidence. ⚠ It hid
            # because the running product of the rescales TELESCOPES back to 1 at the far end of a gene,
            # so no aggregate, endpoint or conservation check could see it.
            if sw.mass_pin and (_lend or g1_l[i]):
                _sg = tg if tpg > 0.0 else og_l[i]
                _sp = tp if tpp > 0.0 else op_l[i]
                _sn = tn if tpn > 0.0 else on_l[i]
                _sv = _sg * E_g_l[i] + (_sp + _sn) * E_r_l[i]
                if _sv > _EPS and M_l[i] > _EPS:
                    _k = M_l[i] / _sv
                    tg, tp, tn = tg * _k, tp * _k, tn * _k
            rg[i], pg[i] = _fuse(og_l[i], pg_own_l[i], tg, tpg)
            rp[i], pp[i] = _fuse(op_l[i], pp_own_l[i], tp, tpp)
            rn[i], pn[i] = _fuse(on_l[i], pn_own_l[i], tn, tpn)
            # measurement + composition precisions are INDEPENDENT evidence → additive fuse
            mg[i] = mg_own_l[i] + tmg
            mp[i] = mp_own_l[i] + tmp
            mn[i] = mn_own_l[i] + tmn
            tau[i] = tau_own_l[i] + ttau

        def publish():
            # back to arrays for the vectorised combine — exact, they are the same doubles
            return tuple(
                np.asarray(a, np.float64) for a in (rg, rp, rn, pg, pp, pn, mg, mp, mn, tau)
            )

        return step, publish

    # ── THE COMBINE: both neighbours transported into this slot's frame, then the message packet ───────
    def _transport(self, nb: NeighbourState, rho_dst, rho_src_a, pop):
        """The VECTORISED twin of :meth:`scan`'s ``step`` — see the do-not-merge note there.

        ⛔ ``nb.state`` arrives ALREADY INDEXED AT THE SOURCE, so this function has no way to read a
        neighbour's relayed belief at the destination. That is TRAPS: a-message-from-the-destinations-belief made structural."""
        sw = self.sw
        cap = self.ctx.capture
        src, valid = nb.src, nb.valid
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = nb.state
        is_bnd_a, ex_a, fp_a, fn_a = self._is_bnd_a, self._ex_a, self._fp_a, self._fn_a
        M, E_g, E_r, SPL = self._M, self._E_g, self._E_r, self._SPL
        og, op, on = self._og, self._op, self._on
        _n_node, logvar_tot = self._n_node, self._logvar_tot
        spl_p, spl_n = self._spl_p, self._spl_n

        # A slot with no frame (no mass ⇒ no rho_tot) cannot reframe: the message passes through at r = 1.
        # Falling back to ``rho_src = 1.0`` instead made ``r`` the destination's ABSOLUTE density — a raw
        # scale masquerading as a ratio.
        framed = valid & (rho_src_a[src] > _EPS) & (rho_dst > _EPS)
        if sw.reframe:
            r = np.where(
                framed, rho_dst / np.maximum(rho_src_a[src], _EPS), np.where(valid, 1.0, 0.0)
            )
        else:
            r = np.where(valid, 1.0, 0.0)
        # GRAFT before the reframe (a density measured AT the source); PEEL after (measured at the dst).
        graft = (ex_a & is_bnd_a[src] & valid) if sw.graft else np.zeros_like(valid)
        gp = np.where(graft, spl_p[src], 0.0)
        gn = np.where(graft, spl_n[src], 0.0)
        # ⭐⭐ THE gDNA SCALE — ``lend`` asks two things of the step: the lambda-emission gate's predicate
        # of the SOURCE (it may lend a composition only if it SUPPLIED both components of the pair) and,
        # of the PAIR, whether the two objects measure the same RNA POPULATION. Where either fails the
        # reframe is a false premise and the gDNA LEVEL crosses UNSCALED.
        #
        # ⚠ **A GRAFT edge does NOT license it, and that is a deliberate divergence from the lambda gate**,
        # which counts the grafted junction precision as RNA supplied. TRAPS: mature-rna-never-crosses-a-seam: mature RNA does not
        # cross an intron<->exon EDGE contiguously, so that EDGE's OWN spanning population is gDNA and
        # unspliced RNA, and the junction flux is a measurement of RNA that lives in the DESTINATION — the
        # routing operator exists precisely because that component cannot cross by imputation. Using it to
        # license the imputation would be circular.
        #
        # ⚠ **The gDNA conjunct is INERT in practice and is kept for faithfulness, not for effect** —
        # recorded because a perturbation dropping it fires no gate, and that is a fact about the solver
        # rather than a hole. ``pg[src] == 0`` with RNA precision live requires the source's own gDNA
        # density to be 0, and then ``rg[src]*r_g`` is 0 whatever the scale.
        lend = pop & (pg > 0.0) & ((pp + pn) > 0.0)
        r_g = np.where(lend, r, np.where(valid, 1.0, 0.0)) if sw.gdna_level_scale else r
        tg, tp, tn = rg * r_g, (rp + gp) * r, (rn + gn) * r
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

        def _dv(p, s2=s2t):
            return np.where(valid & (p > 0.0), 1.0 / (1.0 / np.maximum(p, _EPS) + s2), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)  # full → mode fusion
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)  # measurement (anchor gDNA + spliced RNA)
        ttau = _dv(tau, s2t)  # composition (tau) → the lambda-message
        # the graft's MEASUREMENT precision — never tau-gated (see the twin). ``_sp`` > 0 only on a GRAFT
        # edge, where s2t is identically 0, so the inf→0 substitution below touches only already-masked
        # entries (a zero-count slot has logvar_tot = +inf ⇒ s2t = inf, and ``0*inf`` would nan the masked
        # branch ``np.where`` evaluates).
        _sp = np.where(graft, SPL[:, 0][src], 0.0)
        _sn = np.where(graft, SPL[:, 1][src], 0.0)
        _s2t_spl = np.where(np.isfinite(s2t), s2t, 0.0)
        if sw.graft_frame_var:
            _s2t_spl = _s2t_spl + np.where(graft, graft_frame_logvar(r), 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2t_spl), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2t_spl), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc  # into the mode-fusion precision …
        tmp, tmn = tmp + _spc, tmn + _snc  # … and the measurement stream (a count, never tau)
        # ⭐ P1d — the graft's PREMISE variance, applied to the WHOLE RNA claim after the spliced arm is
        # folded in, because the premise is about the SUM: measured FLAT in the spliced share, so charging
        # the spliced arm alone would reach only 10-93 % of the delivered confidence while the error
        # contaminates 63-95 % of the delivered density.
        _vgp = np.where(graft, self._vgp_prem, 0.0)
        _vgn = np.where(graft, self._vgn_prem, 0.0)
        tpp, tmp = tpp / (1.0 + tpp * _vgp), tmp / (1.0 + tmp * _vgp)
        tpn, tmn = tpn / (1.0 + tpn * _vgn), tmn / (1.0 + tmn * _vgn)

        peel = (is_bnd_a & ex_a[src] & valid) if sw.peel else np.zeros_like(valid)
        (_wp, _vwp), (_wn, _vwn) = self._peel_share(tg, tpg, tp, tn)
        tp = np.where(peel, tp * _wp, tp)
        tn = np.where(peel, tn * _wn, tn)

        def _dv_arr(pr, vv):
            _f = np.isfinite(vv)
            return np.where(peel, np.where(_f, pr / (1.0 + pr * np.where(_f, vv, 0.0)), 0.0), pr)

        tpp, tmp = _dv_arr(tpp, _vwp), _dv_arr(tmp, _vwp)
        tpn, tmn = _dv_arr(tpn, _vwn), _dv_arr(tmn, _vwn)
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        if cap is not None:  # inert: the PRE-PIN state
            cap.setdefault("_pin", []).append(
                {
                    "src": np.asarray(src).copy(),
                    "valid": np.asarray(valid).copy(),
                    "lend": np.asarray(lend).copy(),
                    "r": r.copy(),
                    "r_g": np.asarray(r_g).copy(),
                    "tg": tg.copy(),
                    "tp": tp.copy(),
                    "tn": tn.copy(),
                    "tpg": tpg.copy(),
                    "tpp": tpp.copy(),
                    "tpn": tpn.copy(),
                    "s2t": np.where(np.isfinite(s2t), s2t, 0.0),
                    "n_src": np.asarray(_n_node)[src].copy(),
                    "spl_p": _sp.copy(),
                    "spl_n": _sn.copy(),
                    "spl_prec": (_spc + _snc).copy(),
                    "graft": np.asarray(graft).copy(),
                }
            )
        # ── P1e: the conservation SURPRISE as a DerSimonian-Laird damping term ─────────────────────────
        # The claim asserts ``S = Sigma_c rho_c E_c`` fragments; the slot observed ``M``. That is an
        # IDENTITY, and the pin restores it by fiat and DISCARDS the residual. Price the residual instead.
        # ⚠⚠ **PARTLY A DEBT — THIS PRICES A BIAS AS A VARIANCE.** On a large share of its firing mass
        # ``delta`` is systematic, and a variance cannot move a mode toward truth. It is landed because it
        # is the only change measured to improve ACCURACY and honest PRECISION together, not because the
        # bias half is derived. **When the bias strata are diagnosed, this term must SHRINK.**
        _p3 = np.stack([tpg, tpp, tpn], axis=-1)
        _sup = _p3 > 0.0
        _mc = np.where(_sup, np.stack([tg, tp, tn], axis=-1), np.stack([og, op, on], axis=-1))
        _mc = _mc * np.stack([E_g, E_r, E_r], axis=-1)
        _S = _mc.sum(axis=-1)
        _okc = valid & (_S > _EPS) & (M > _EPS)
        _al = _mc / np.maximum(_S, _EPS)[..., None]
        _vc = np.where(_sup, 1.0 / np.maximum(_p3, _EPS), 0.0)
        _s2c = (np.where(np.isfinite(s2t), s2t, 0.0) + 1.0 / np.maximum(_n_node[src], _EPS))[
            ..., None
        ]
        _sv = np.where(_sup, _s2c + _al * np.maximum(_vc - _s2c, 0.0), 0.0)
        _aSa = np.sum(_al * _sv, axis=-1)
        _dlt = np.where(_okc, np.log(np.maximum(M, _EPS) / np.maximum(_S, _EPS)), 0.0)
        _den = _aSa + 1.0 / np.maximum(_n_node, _EPS)
        _b2 = np.maximum(_dlt * _dlt - _den, 0.0)
        # ── ⭐ THE SCOPE: only the OVER-claim direction is evidence against the MESSAGE ────────────────
        # ``S`` is a COMPLETE budget: the partial-claim semantics fill every component the message does not
        # supply from the slot's OWN density. So a shortfall can be the slot's own density being too low
        # just as easily as the message being wrong — it does not attribute. The OVER-claim direction does:
        # every density is non-negative, so nothing the unsupplied components could be would rescue a
        # budget that already exceeds ``M``.
        _b2 = np.where(_dlt < 0.0, _b2, 0.0)
        # ⭐ THE COMMON DIRECTION (the derived law). Because ``alpha^T 1 == 1`` — alpha is a share vector
        # over the same budget — adding the SCALAR to every supplied component's log-variance satisfies
        # the constraint exactly, and it leaves ``Var(lambda)`` **identically unchanged**: a common shift of
        # both arms cannot move the split. The rank-1 form borrows the CONDITIONAL MEAN's direction and
        # applies it as a variance inflation, which MC shows over-damps lambda 5x; it was implemented,
        # MC-refuted and DELETED — do not rebuild it.
        _dg = _dp = _dn = np.where(_sup.any(axis=-1), _b2, 0.0)
        if sw.conservation_var:
            tpg, tmg = tpg / (1.0 + tpg * _dg), tmg / (1.0 + tmg * _dg)
            tpp, tmp = tpp / (1.0 + tpp * _dp), tmp / (1.0 + tmp * _dp)
            tpn, tmn = tpn / (1.0 + tpn * _dn), tmn / (1.0 + tmn * _dn)
        # The mass frame — for the mismatch COMPARISON only. The delivered ``tg/tp/tn`` are left exactly as
        # the source measured and the reframe delivered them: a component's LEVEL is an absolute rate, and
        # re-normalising it against a budget built from the destination's own belief is what made the
        # message self-confirming. The other two claims are unaffected either way — the lambda and tilt
        # claims are scale-free, so the pin's common factor cancels from them identically.
        pin_g, pin_p, pin_n = self._pin_v(tg, tp, tn, tpg, tpp, tpn)
        # ── THE lambda-EMISSION GATE (structural, and PRIOR to any damping question) ───────────────────
        # A composition message is a claim about the SPLIT. A source carrying only ONE component has no
        # such claim to make — lambda is not "large" for it, it is UNDEFINED. The canonical case is exactly
        # ``intergenic | seam | EXON``: RNA cannot cross a gene boundary, so the seam is structurally
        # RNA-free while the exon it feeds has RNA. The combine builds the lambda mode as ``mo_g - mo_R``
        # with ``mo_R`` floored, so a zero-RNA message silently becomes the maximally confident assertion
        # "this slot is 100 % gDNA" — a numerical artifact of the floor — while its precision is real,
        # having been accumulated by slots that never said that. Gating at EMISSION is structural, the same
        # kind of presence test as the strand gates, and no threshold. ⭐ It is NOT a loss of information:
        # a pure-gDNA source's authority is a DENSITY LEVEL, and that travels on the measurement stream.
        # "Supplied" is a statement about PRECISION, not about the density's value: a component at zero
        # density with live precision IS supplied — it is the claim "there is none of this here", which is
        # exactly a composition claim.
        if sw.lam_emission_gate:
            ttau = np.where((tpg > 0.0) & ((tpp + tpn) > 0.0), ttau, 0.0)
        if sw.mismatch_var:
            # ── the COMPOSITION half of the cliff cost: the DL mismatch deflation, in the PINNED frame.
            # Both sides then account for the same total, so the common scale (the reframe residual) is
            # gone from the gap and only the share drift is left. Every stream is deflated — measured: the
            # pin recovers only when the composition tau-stream is damped alongside the measurement one.
            g_g, c_g = mismatch_gap(pin_g, og)
            g_p, c_p = mismatch_gap(pin_p, op)
            g_n, c_n = mismatch_gap(pin_n, on)
            tpg = mismatch_deflate(tpg, g_g, c_g, self._v_own_g)
            tpp = mismatch_deflate(tpp, g_p, c_p, self._v_own_r)
            tpn = mismatch_deflate(tpn, g_n, c_n, self._v_own_r)
            tmg = mismatch_deflate(tmg, g_g, c_g, self._v_own_g)
            tmp = mismatch_deflate(tmp, g_p, c_p, self._v_own_r)
            tmn = mismatch_deflate(tmn, g_n, c_n, self._v_own_r)
            # the composition stream is ONE DOF, so its gap is measured on the lambda axis — the message's
            # log(f_g/f_R) minus the own belief's, built from the SAME quantities the combine builds the
            # lambda mode from, so the gap is exactly the error of the claim psi receives. A contradiction
            # on EITHER arm contradicts the claim (a message with no RNA at all is asserting lambda = +inf).
            g_R, c_R = mismatch_gap(pin_p + pin_n, op + on)
            _tau_pre = ttau
            ttau = mismatch_deflate(ttau, g_g - g_R, c_g | c_R, self._v_own_lam)
            if cap is not None:  # inert: the per-message gaps + the tau-stream kill
                cap.setdefault("_dl", []).append(
                    {
                        "G_g": g_g.copy(),
                        "G_p": g_p.copy(),
                        "G_n": g_n.copy(),
                        "G_lam": (g_g - g_R),
                        "contra": (c_g | c_R).copy(),
                        "tau_pre": _tau_pre.copy(),
                        "tau_post": ttau.copy(),
                    }
                )

        # ── THE MESSAGE PACKET: a LEVEL claim per component, and SEPARATE SPLIT + TILT claims ──────────
        # A message must let the destination do three independent things: set each component's LEVEL, set
        # the SPLIT, and set the TILT. Those are three different claims at three different precisions — so
        # the message states the split and the tilt EXPLICITLY. Reading them back off the fused densities
        # instead delivers them weighted by the LEVEL precisions rather than their own, a mode/precision
        # mismatch: a message with almost no composition evidence but a large, well-counted density would
        # set the split, and one with real composition evidence but little mass would not.
        # Both are scale-free — the slot's mass cancels out of ``log(f_g/f_R)`` and out of the tilt — so a
        # message can state them without the destination reconstructing any density, and no rescaling of
        # the level claims can perturb them.
        _tR = tp + tn
        tlam = np.where(
            (tg > _EPS) & (_tR > _EPS),
            np.log(np.maximum(tg * E_g, _EPS)) - np.log(np.maximum(_tR * E_r, _EPS)),
            0.0,
        )
        # ⭐ the tilt is the ANGLE ``arcsin(tau)``, which is psi's own coordinate — TRAPS: off-grid-message-mode: a raw
        # log-odds delivered into this slot spans 2.9x the grid's whole domain and pins the tilt at the
        # boundary, and that was 74 % of a zero control's error.
        tth = np.arcsin(
            np.clip(np.where(_tR > _EPS, (tp - tn) / np.maximum(_tR, _EPS), 0.0), -1.0, 1.0)
        )
        return tg, tp, tn, tpg, tpp, tpn, tmg, tmp, tmn, ttau, tlam, tth

    def _pin_v(self, g, p, n, pg_, pp_, pn_):
        """The message's densities re-expressed in the destination's MASS FRAME.

        ⚠⚠ **THE RESULT IS A COMPARISON FRAME FOR THE MISMATCH TEST — IT IS NOT THE MESSAGE.** It used to
        replace the delivered densities, and that was a belief-propagation violation ON THE MODE: the
        budget fills every component the message does not supply from the DESTINATION'S OWN density, so the
        delivered claim became a function of the destination's own self-solve. Measured: with an RNA-only
        message and ``E_g = E_r`` the delivered RNA fraction was **exactly ``1/(1 + f_g_own)``** (verified
        to 2.1e-16) — the slot confirming itself. On unstranded data ``f_g_own`` is the uninformative 1/2,
        so the pin reserved **33.6 %** of the budget for gDNA the message never claimed and a zero-gDNA
        library read back **29.3 %** gDNA. The reservation WAS the false-positive rate.

        Feeding the DL mismatch test is what it is legitimately for: that test is a two-study random-effects
        comparison **against the destination's own self-solve**, so destination information is the point of
        it — and it sets a VARIANCE, which can mis-weight a message but cannot invent a location."""
        E_g, E_r, M = self._E_g, self._E_r, self._M
        sg = np.where(pg_ > 0.0, g, self._og)
        sp = np.where(pp_ > 0.0, p, self._op)
        sn = np.where(pn_ > 0.0, n, self._on)
        s_ = sg * E_g + (sp + sn) * E_r
        k = np.where((s_ > _EPS) & (M > _EPS), M / np.maximum(s_, _EPS), 1.0)
        return g * k, p * k, n * k

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        """The four ψ channels, from the two NEIGHBOUR states only.

        ⚠ The SAME role pairing as the scan: the left-hand message's source is the genomic-LOW neighbour,
        so the destination presents ``rho_lo`` and the source ``rho_hi``; the right-hand message is the
        mirror. The two calls differ ONLY in which neighbour they read and, with it, which flank total each
        side presents."""
        sw = self.sw
        cap = self.ctx.capture
        E_g, E_r, M = self._E_g, self._E_r, self._M
        rho_lo, rho_hi = self._rho_lo, self._rho_hi

        ag, ap, an, apg, app, apn, amg, amp, amn, atau, alam, ath = self._transport(
            left, rho_lo, rho_hi, self._pop_l_a
        )
        bg, bp, bn, bpg, bpp, bpn, bmg, bmp, bmn, btau, blam, bth = self._transport(
            right, rho_hi, rho_lo, self._pop_r_a
        )

        def _fuse_add(
            a, b
        ):  # additive (inverse-variance) fuse of two independent precision streams
            return np.asarray(a, np.float64) + np.asarray(b, np.float64)

        def _fuse_v(a, pa, b, pb):
            p = pa + pb
            return np.where(p > _EPS, (pa * a + pb * b) / np.maximum(p, _EPS), 0.0), p

        cg, cpg = _fuse_v(ag, apg, bg, bpg)  # density MODE (full precision-weighted)
        cp, cpp = _fuse_v(ap, app, bp, bpp)
        cn, cpn = _fuse_v(an, apn, bn, bpn)
        cm_g, cm_p, cm_n = _fuse_add(amg, bmg), _fuse_add(amp, bmp), _fuse_add(amn, bmn)
        c_tau = _fuse_add(atau, btau)
        mo_g = np.log(np.maximum(cg * E_g / np.maximum(M, _EPS), _EPS))
        mo_p = np.log(np.maximum(cp * E_r / np.maximum(M, _EPS), _EPS))
        mo_n = np.log(np.maximum(cn * E_r / np.maximum(M, _EPS), _EPS))
        # ── THE THREE-STREAM SINGLE-lambda COMBINE (the the-single-lambda-combine rank-1 fix) ─────────────────────────────────
        # (1) COMPOSITION → ONE lambda-message, precision ``c_tau`` (the fused Schur tau) — so psi counts
        #     the composition DOF ONCE, not twice. **Each claim is fused by ITS OWN precision**: the
        #     lambda mode is the tau-weighted mean of the two messages' lambda, not a ratio read back off
        #     the density fuse (where each component was averaged by its own MODE-FUSION precision). That
        #     mismatch delivered the split at a confidence it was never weighted by.
        cR = cp + cn
        lam_msg = np.where(c_tau > _EPS, (atau * alam + btau * blam) / np.maximum(c_tau, _EPS), 0.0)
        # A lambda message exists only where BOTH components of the pair reached this slot — the
        # structural presence test. The per-message gate cannot catch every case: a message may carry an
        # RNA DENSITY while contributing zero mode-fusion PRECISION, in which case ``cR`` collapses here.
        c_tau = np.where((cg > _EPS) & (cR > _EPS), c_tau, 0.0)
        # (2) ANCHOR gDNA MEASUREMENT → gdna_imp. (3) SPLICED RNA MEASUREMENT → rna_imp. INDEPENDENT of
        #     the composition, so fused separately (an RNA-only spliced measurement constrains f_g via
        #     f_R with NO gDNA info). (4) the TILT (AMBIG only) — fused by the MEASURED RNA that carries
        #     it, for the same reason as lambda.
        _tha, _thb = amp + amn, bmp + bmn
        th_msg = np.where(
            (_tha + _thb) > _EPS, (_tha * ath + _thb * bth) / np.maximum(_tha + _thb, _EPS), 0.0
        )
        th_prec = np.where(self._is_amb, cm_p + cm_n, 0.0)

        if cap is not None:  # inert diagnostic: the fused per-component densities + the frames
            cap.setdefault("_uni", []).append(
                {
                    "cg": cg.copy(),
                    "cp": cp.copy(),
                    "cn": cn.copy(),
                    "pg": cpg.copy(),
                    "pp": cpp.copy(),
                    "pn": cpn.copy(),
                    "ag": ag.copy(),
                    "ap": ap.copy(),
                    "an": an.copy(),
                    "bg": bg.copy(),
                    "bp": bp.copy(),
                    "bn": bn.copy(),
                    "apg": apg.copy(),
                    "app": app.copy(),
                    "bpg": bpg.copy(),
                    "bpp": bpp.copy(),
                    # ⭐ the FLANK PAIR, both published: an instrument reconstructing a hop's ``r`` must
                    # pair them by role — ``rho_lo[high slot] / rho_hi[low slot]`` — and one array cannot
                    # express that.
                    "rho_lo": rho_lo.copy(),
                    "rho_hi": rho_hi.copy(),
                    "mo_g": mo_g.copy(),
                    "mo_p": mo_p.copy(),
                    "mo_n": mo_n.copy(),
                    "lam_msg": lam_msg.copy(),
                    "c_tau": c_tau.copy(),
                    "cm_g": cm_g.copy(),
                    "cm_p": cm_p.copy(),
                    "cm_n": cm_n.copy(),
                    "amp": amp.copy(),
                    "bmp": bmp.copy(),
                    "amn": amn.copy(),
                    "bmn": bmn.copy(),
                    "cpg": cpg.copy(),
                }
            )
            # ⚠ The RAW per-slot relay state (``fwd_*``/``bwd_*``) is published by the BACKBONE, not
            # here — ``deliver`` is handed those arrays already gathered AT THE SOURCE, which is what makes
            # TRAPS: a-message-from-the-destinations-belief structural, so this policy genuinely cannot see the un-indexed form. That is the
            # assertion working, not a gap.
            # ⚠ ``prec_*`` here is the MODE-FUSION precision, NOT the precision ψ receives — ψ gets the
            # separate MEASUREMENT stream ``cm_*``, which is published beside it in ``_uni``. The two are
            # different quantities and an instrument that conflates them reads a channel's confidence off
            # the wrong array. Kept as the shipped solver published it, because the dissect loop's
            # replays are calibrated against exactly these keys.
            cap.update(
                mode_g=mo_g,
                prec_g=cpg,
                mode_p=mo_p,
                prec_p=cpp,
                mode_n=mo_n,
                prec_n=cpn,
            )
        # ⛔ The RNA MEASUREMENT psi factor was ablated behind a flag once. Do not re-try it as a deletion:
        # it carries 75 % of the posterior precision on the confidently-wrong unstranded x capture-OFF
        # exons, and it is the only thing that lets a zero-gDNA library say "my mass is all RNA".
        return PsiMessage(
            gdna_mode=mo_g if sw.gdna_channel else None,
            gdna_prec=cm_g if sw.gdna_channel else None,
            rna_mode=(mo_p, mo_n) if sw.rna_channel else None,
            rna_prec=(cm_p, cm_n) if sw.rna_channel else None,
            lam_mode=lam_msg if sw.lam_channel else None,
            lam_prec=c_tau if sw.lam_channel else None,
            theta_mode=th_msg if sw.theta_channel else None,
            theta_prec=th_prec if sw.theta_channel else None,
        )
