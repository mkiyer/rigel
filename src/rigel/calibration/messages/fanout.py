"""The PASS-0 FAN-OUT policy — stages 3 AND 4 of the first-pass redesign: solved single-stranded
introns speak to their flanking ``intron|exon`` boundaries, and the boundaries speak on to the
licensed exons. Nothing else says anything, and nothing is spoken back.

       Gate: ``tests/calibration/test_fanout_policy.py``

**Stage 3 (depth 1)** — the intron's OWN local solve on the λ channel: mode ``logit(f_g)`` clipped
into the grid domain (`TRAPS: off-grid-message-mode`), precision ``tau_lam`` (the strand +
intron-factory Fisher; the reference's location does not enter, §9c.2's ruling), NO drift allowance —
measured, not assumed: `hop_currency.py`'s excess-over-floor for this hop is 0.2–0.4 % on both
panels, capture ON and OFF (2026-08-23). The premise (one population on both sides) is certified by
the stage-0 matrix, 32/32; `intron_phi`'s refutation was the UNCONDITIONAL form with no fuse.

**Stage 4 (depth 2)** — the scan's ``step`` COMPOSES at each intron-flanking boundary: the intron's
claim fused with the boundary's OWN strand-λ evidence (precision-weighted; the boundary is the
SOURCE of the onward hop, so reading its own solve is what a message is). ``deliver`` then runs the
jointly-derived transfer (:func:`flank_to_exon_lambda`) per licensed flank, taking the sj FACE on the
exon's side of the flank — the population of a message is direction-dependent — and routes by the
completeness bit: a COMPLETE flank's claim is a two-sided λ estimate (two complete flanks fuse,
precisions adding — independent data); an INCOMPLETE flank's is the one-sided at-least-this-much-RNA
bound on the exon's live strand share (`PsiMessage.rna_one_sided`), at the Jacobian precision
``tau/f_g²``; two incomplete flanks keep the BINDING bound (the larger RNA floor — a precision fuse
of two one-sided bounds is not itself a bound, so it is not used).

⛔ Licences: sources are ``ss_intron_region``; stage-3 destinations are IMPLIED (a boundary flanking
an ss intron is always ``ss_intron_boundary`` — derivation in the stage-3 history, `762065bc`^);
stage-4 destinations are the licensed flanks (``exon_flank_left/right``). ⭐ **The destination mask on
the stage-3 channel is now REAL, exactly as the depth-1 docstring predicted**: with composed states
in the scan, an unmasked deliver would hand the boundary's composed claim BACK to the intron — its
own strand evidence returned at message precision — so boundaries alone may receive the raw-claim
fuse, and regions receive it nowhere. The no-echo gate pins this.
"""

from __future__ import annotations

import numpy as np

from . import NeighbourState, PsiMessage, StepContext
from .variance import mismatch_deflate

__all__ = ["FanOutPolicy", "flank_to_exon_lambda"]

_EPS = 1.0e-12


class _FanOutRelay:
    def __init__(self, ctx: StepContext):
        claims = ctx.claims
        self._L = L = float(ctx.logodds_window)
        f = np.clip(np.asarray(ctx.own.f_g, np.float64), _EPS, 1.0 - _EPS)
        lam_own = np.clip(np.log(f) - np.log1p(-f), -L, L)
        tau_own = np.asarray(ctx.own.tau_lam, np.float64)
        source = np.asarray(claims.ss_intron_region, bool)
        # the raw per-slot claim each SOURCE publishes; a non-source region publishes (0, 0) — no
        # claim, no precision, one statement (`TRAPS: zero-the-precision-with-the-value`).
        self._lam0 = np.where(source, lam_own, 0.0)
        self._tau0 = np.where(source, tau_own, 0.0)
        self._dst_b = np.asarray(claims.ss_intron_boundary, bool)
        # the boundary's OWN strand-λ evidence, for the stage-4 compose (gated to the boundaries the
        # fan-out passes through, so an unrelated slot's solve can never leak into a composed state).
        self._own_lam = lam_own
        self._own_tau = np.where(self._dst_b, tau_own, 0.0)
        self._flank = (
            np.asarray(claims.exon_flank_left, bool),
            np.asarray(claims.exon_flank_right, bool),
        )
        self._complete = (
            np.asarray(claims.exon_flank_left_complete, bool),
            np.asarray(claims.exon_flank_right_complete, bool),
        )
        # the destination's own λ-axis variance, for the stage-4 mismatch deflation below.
        # ⛔ No ``struct_lock`` branch: a licensed exon is single-stranded WITH mRNA, so it is
        #    ``solvable`` and therefore never locked unless it is EMPTY — and an empty exon carries no
        #    mass to misplace. A clause that provably cannot fire is the kind this policy has already
        #    deleted twice (the stage-3 destination mask, the ``tau_lam > 0`` guard).
        self._v_own = np.where(tau_own > 0.0, 1.0 / np.maximum(tau_own, _EPS), np.inf)
        self._ctx = ctx

    def scan(self, *, backward: bool):
        """Depth-TWO fan-out: regions carry their RAW claims (never rewritten), and the step composes
        the incoming intron claim with the boundary's own λ at each pass-through boundary — so the
        state a boundary publishes toward the exon already holds both. Each direction's scan carries
        its own side's intron, which is exactly the flank the destination exon reads."""
        lam = self._lam0.copy()
        tau = self._tau0.copy()
        dst_b = self._dst_b
        own_lam, own_tau = self._own_lam, self._own_tau

        def step(s: int, i: int) -> None:
            ts = tau[s]
            if ts > 0.0 and dst_b[i]:
                p = ts + own_tau[i]
                lam[i] = (ts * lam[s] + own_tau[i] * own_lam[i]) / p
                tau[i] = p

        def publish():
            return lam, tau

        return step, publish

    def _received(self, mode, prec, arrived):
        """The DESTINATION's reading of an arriving stage-4 claim: keep the value, discount the
        confidence by how far the claim sits from what this slot's own data says.

        The transfer that produced ``mode`` prices SAMPLING noise only — nothing in it states that the
        flank's composition IS the exon's, and under hybrid capture it measurably is not: every probe
        base lands in an EXON and introns carry exactly zero, so an intron|exon flank is a probe cliff
        and its gDNA route arrives ~``e^-1`` low. Worse, that transfer is the message layer's only
        precision-AMPLIFYING single-source transform (``tau_e = tau_b/a**2`` with ``a = 1 - f*w_mu ->
        0``), where every operator in the shipped relay only ever deflates.

        ⭐ The discount is the relay's OWN law, not a new one — the DerSimonian-Laird between-source
        mismatch variance (`variance.mismatch_deflate`, shipped ON as ``RelaySwitches.mismatch_var``),
        on the same lambda axis, at the OBSERVED gap ``G = lam_msg - lam_own``::

            b_hat^2 = max(0, G^2 - v_msg - v_own)      p = 1 / max(v_msg, G^2 - v_own)

        whose closed form states the safety property exactly: *an arriving claim may out-weigh this
        slot's own belief only when it agrees with it to within* ``sqrt(2)*sigma_own``. Probe placement
        is unknowable; its FOOTPRINT — this disagreement — is not.

        ⭐⭐ **AND THE REGIME THAT MAKES IT SAFE IS DERIVED, NOT GATED.** Where the destination has no
        own composition evidence (``tau_lam = 0``, which at ``kappa = 1/2`` is EVERY exon, because the
        strand lambda-term is identically zero there) ``v_own`` is infinite, ``b_hat^2`` is 0, and the
        claim passes UNTOUCHED. That is the whole unstranded half of the panel, both zero controls and
        the deferred stratum — the cells this policy WINS — protected by the algebra rather than by a
        condition someone has to remember.

        ⛔ ``arrived`` scopes this to the STAGE-4 destinations, and that scope is measured: stage 3's
        premise (an intron and its flanking boundary are ONE population) is certified by the stage-0
        matrix, while discounting the boundary claim too costs boundary wins.
        ⛔ ``contradicted`` is structurally False here — it marks a component asserted ABSENT, and a
        lambda claim clipped into ``[-L, +L]`` cannot assert infinite odds.
        """
        if not bool(np.any(arrived)):
            return mode, prec
        damped = mismatch_deflate(
            prec, mode - self._own_lam, np.zeros(prec.shape, bool), self._v_own
        )
        return mode, np.where(arrived, damped, prec)

    def deliver(self, left: NeighbourState, right: NeighbourState) -> PsiMessage:
        ll, lt = (np.asarray(a, np.float64) for a in left.state)
        rl, rt = (np.asarray(a, np.float64) for a in right.state)
        lt = np.where(left.valid, lt, 0.0)
        rt = np.where(right.valid, rt, 0.0)

        # ── stage 3: boundaries ← the RAW flanking-region claims (regions carry only those). Two
        #    licensed flanks (an intron|intron boundary) are two independent measurements of one
        #    shared composition: precision-weighted mean, precisions ADD. ──
        prec = np.where(self._dst_b, lt + rt, 0.0)
        mode = np.where(prec > 0.0, (lt * ll + rt * rl) / np.maximum(lt + rt, _EPS), 0.0)

        # ── stage 4: licensed exons ← the flanks' COMPOSED states through the joint transfer ──
        ctx = self._ctx
        g = ctx.geometry
        eg = np.asarray(ctx.eff_gdna, np.float64)
        er = np.asarray(ctx.eff_rna, np.float64)
        U = np.asarray(ctx.n_slot, np.float64)
        # the exon's single strand picks the sj column (a stage-0 licence guarantees ss)
        col = np.where(np.asarray(ctx.free_pos, bool), 0, 1)
        rows = np.arange(self._lam0.shape[0])
        sides = []
        for flank, complete, st_lam, st_tau, valid, src, S_all, E_all in (
            # an exon's LEFT flank delivers RIGHTWARD ⇒ the flux entering the exon is the flank's
            # HIGH face; the RIGHT flank mirrors. The faces are the direction-dependent population.
            (
                self._flank[0],
                self._complete[0],
                ll,
                lt,
                left.valid,
                left.src,
                np.asarray(g.sj_count_hi, np.float64),
                np.asarray(g.eff_sj_hi, np.float64),
            ),
            (
                self._flank[1],
                self._complete[1],
                rl,
                rt,
                right.valid,
                right.src,
                np.asarray(g.sj_count_lo, np.float64),
                np.asarray(g.eff_sj_lo, np.float64),
            ),
        ):
            lam_e, tau_e = flank_to_exon_lambda(
                st_lam,
                st_tau,
                U[src],
                S_all[src, col[rows]],
                E_all[src, col[rows]],
                eg[src],
                er[src],
                eg,
                er,
                self._L,
            )
            sides.append((lam_e, np.where(flank & valid, tau_e, 0.0), complete))

        (lam_l, tau_l, cm_l), (lam_r, tau_r, cm_r) = sides
        # two-sided: the COMPLETE sides fuse on λ (exon slots are disjoint from boundary slots)
        t2l = np.where(cm_l, tau_l, 0.0)
        t2r = np.where(cm_r, tau_r, 0.0)
        p2 = t2l + t2r
        mode = np.where(p2 > 0.0, (t2l * lam_l + t2r * lam_r) / np.maximum(p2, _EPS), mode)
        prec = np.where(p2 > 0.0, p2, prec)

        # ── RECEPTION: the destination decides what to believe of what arrived ────────────────────
        # ⭐⭐⭐ **THE PARADIGM, and it is why this is not a contract breach** (owner, 2026-08-23). A
        # SENDER publishes its claim UNCHANGED — it does not tailor, scale or hedge it for whoever is
        # listening. Deciding how much of an arriving claim to believe is the RECIPIENT's job, and the
        # recipient here is the exon. So an exon that receives a composition wildly at odds with what
        # its own data says may discount it, and doing so reads its OWN belief, never the sender's.
        # ⛔ That is the line `TRAPS: a-message-from-the-destinations-belief` actually draws: its nine
        # costumes all BUILT a claim's VALUE out of the destination's belief, which manufactures
        # agreement out of nothing. Reception cannot — :func:`_received` only ever LOWERS a precision
        # and never touches a mode, so it can discard information but never invent it.
        mode, prec = self._received(mode, prec, p2 > 0.0)
        # one-sided: the INCOMPLETE sides — keep the BINDING bound (smaller λ = larger RNA floor)
        t1l = np.where(~cm_l, tau_l, 0.0)
        t1r = np.where(~cm_r, tau_r, 0.0)
        left_binds = (t1l > 0.0) & ((t1r <= 0.0) | (lam_l <= lam_r))
        lam1 = np.where(left_binds, lam_l, lam_r)
        tau1 = np.where(left_binds, t1l, t1r)
        bound = tau1 > 0.0
        if not bool(bound.any()):
            return PsiMessage(lam_mode=mode, lam_prec=prec)
        f1 = 1.0 / (1.0 + np.exp(-lam1))
        rna_m = np.where(bound, np.log1p(-f1), 0.0)
        rna_p = np.where(bound, tau1 / np.maximum(f1 * f1, _EPS), 0.0)
        pos = col == 0
        return PsiMessage(
            lam_mode=mode,
            lam_prec=prec,
            rna_mode=(np.where(pos, rna_m, 0.0), np.where(~pos, rna_m, 0.0)),
            rna_prec=(np.where(pos & bound, rna_p, 0.0), np.where(~pos & bound, rna_p, 0.0)),
            rna_one_sided=bound,
        )


class FanOutPolicy:
    """Pass-0's message policy: intron → boundary (λ), boundary → licensed exon (λ two-sided, or the
    one-sided RNA bound where the flank's account is incomplete). Nothing else, in either direction."""

    name = "fanout"

    def prepare(self, ctx: StepContext) -> _FanOutRelay:
        return _FanOutRelay(ctx)


def flank_to_exon_lambda(
    lam_b,
    tau_b,
    U,
    S_face,
    eff_sj_face,
    eff_gx,
    eff_rx,
    eff_gc,
    eff_rc,
    logodds_window,
):
    """Stage 4's transfer — route-merge ∘ reframe, derived JOINTLY: the flank's composed unspliced
    claim ``(lam_b, tau_b)`` plus the spliced flux entering the exon, delivered as ``(lam, tau)`` in
    the EXON's frame.

    The routes merge as DENSITIES (each count over ITS OWN opportunity — the two footings differ by a
    measured 7–10×, so raw-count pooling is inadmissible), and the ratio is re-formed through the
    exon's opportunities, so the absolute level CANCELS — a uniform capture pull moves nothing, which
    is what the currency map certified about this hop::

        f     = sigma(lam_b)                       rho_g  = f·U/E_gx        rho_nu = (1−f)·U/E_rx
        rho_mu = S/E_sj                            rho_r  = rho_nu + rho_mu
        lam   = log(rho_g·E_gc) − log(rho_r·E_rc)          clipped into [−L, +L]

    ⛔ **ONE function, not route-merge THEN reframe through a (density, variance) interface — and that
    is `TRAPS: two-gaussians-one-latent` speaking**: ``rho_g`` and ``rho_nu`` share ``lam_b`` and
    ``U``, so a split interface loses their covariance and double-counts both. Carried jointly (the
    delta method on the three RAW statistics ``lam_b``, ``log U``, ``log S`` — each entering once,
    D.1's requirement)::

        Var(lam) = [(1−f) + f·rho_nu/rho_r]²/tau_b + (rho_mu/rho_r)²·(1/U + 1/S)

    whose limits are the falsification: at ``S = 0`` the transfer is a pure reparameterization —
    ``lam = lam_b + log(E_rx·E_gc/(E_gx·E_rc))`` at EXACTLY ``1/tau_b``, the count noise cancelling in
    the ratio — and at ``S ≫`` the gDNA side's Poisson noise carries. Monte-Carlo-gated.

    A claim requires ``tau_b > 0``, live opportunities, and live densities (``U = 0`` withdraws the claim through ``rho_g``) — anything else returns ``tau = 0``,
    value and precision withdrawn in one statement. ``S_face``/``eff_sj_face`` are the FACE of the sj
    flux on the exon's side (`RegionGeometry.sj_count_lo/hi`, ``eff_sj_lo/hi`` — the population of a
    message is direction-dependent).
    """
    lam_b = np.asarray(lam_b, np.float64)
    tau_b = np.asarray(tau_b, np.float64)
    U = np.asarray(U, np.float64)
    S = np.asarray(S_face, np.float64)
    E_sj = np.asarray(eff_sj_face, np.float64)
    gx = np.asarray(eff_gx, np.float64)
    rx = np.asarray(eff_rx, np.float64)
    gc = np.asarray(eff_gc, np.float64)
    rc = np.asarray(eff_rc, np.float64)
    L = float(logodds_window)

    # ⭐ no explicit U guard: with U = 0 the gDNA density is 0 and the ``rho_g > 0`` refinement below
    #   withdraws the claim — the perturbation sweep proved a separate conjunct did no work.
    live = (tau_b > 0.0) & (gx > 0.0) & (rx > 0.0) & (gc > 0.0) & (rc > 0.0)
    f = 1.0 / (1.0 + np.exp(-np.clip(lam_b, -L, L)))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_g = np.where(live, f * U / np.maximum(gx, _EPS), 0.0)
        rho_nu = np.where(live, (1.0 - f) * U / np.maximum(rx, _EPS), 0.0)
        rho_mu = np.where((S > 0.0) & (E_sj > 0.0), S / np.maximum(E_sj, _EPS), 0.0)
        rho_r = rho_nu + rho_mu
        live = live & (rho_g > 0.0) & (rho_r > 0.0)
        lam = np.where(
            live,
            np.clip(
                np.log(np.maximum(rho_g * gc, _EPS)) - np.log(np.maximum(rho_r * rc, _EPS)),
                -L,
                L,
            ),
            0.0,
        )
        w_mu = np.where(live, rho_mu / np.maximum(rho_r, _EPS), 0.0)
        a = np.where(live, (1.0 - f) + f * rho_nu / np.maximum(rho_r, _EPS), 0.0)
        var = a * a / np.maximum(tau_b, _EPS) + w_mu * w_mu * (
            1.0 / np.maximum(U, _EPS) + np.where(S > 0.0, 1.0 / np.maximum(S, _EPS), 0.0)
        )
        tau = np.where(live & (var > 0.0), 1.0 / np.maximum(var, _EPS), 0.0)
    return lam, tau
