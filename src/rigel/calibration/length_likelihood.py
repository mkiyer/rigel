"""The FRAGMENT-LENGTH composition likelihood — the fourth information source, on the ``λ`` axis.

     (P2) · Gate: `tests/calibration/test_length_likelihood.py`

⭐ **WHY THIS EXISTS.** ``node_init`` gives every slot its composition precision from three live sources:
the structural lock, the intron-factory density deconvolution, and the strand Beta-Binomial. Measured on
the chr22 pilot (P0), that leaves **13.3–40.1 % of library mass with no own evidence at all**
in *every* condition — the Schur gate silences strand on both-strand (AMBIG) nodes by construction — and
**100 %** on an unstranded or zero-gDNA library, because the strand Fisher information is
``I(f_g) ∝ (2κ−1)²`` and is identically 0 at ``κ = ½``.

Meanwhile the accumulator has measured, per object and per strand, the two length channels
``inv_length_sum`` (``Σ 1/weight``) and ``length_sum`` (``Σ L``) — and **nothing reads them.** This module
turns them into a log-likelihood on the same ``σ(λ)`` grid the strand term lives on.

⭐ **THE MODEL, IN FOUR LINES.** At one object let ``A(w)`` be the opportunity (the number of start
positions from which a length-``w`` fragment lands here) and ``u(w)`` the weight the accumulator deposited.
A component-``c`` fragment that LANDED here is drawn from the **opportunity-tilted** pmf::

    g_c(w) = f_c(w)·A(w) / E_c[A]          and  E_c[A] IS eff_gdna / eff_rna at this slot

Conditional on the integer count ``N``, each landed fragment is gDNA with probability ``pi`` — which is
exactly ``f_g``, the quantity the ``λ`` grid parameterises — so ``(D, S) = (Σu, Σw)`` is a sum of ``N``
i.i.d. draws from the mixture ``pi·g_g + (1−pi)·g_r``::

    mean  = N · ( pi·m1 + (1−pi)·m1' ,  pi·m2 + (1−pi)·m2' )
    cov   = N · [ mixture second moments − mixture first moments squared ]
    ll(pi) = −½ (x − mean)' cov⁻¹ (x − mean) − ½ log det cov

⭐ **The Gaussian is on a SUM, not on a ratio, and that is what makes it legitimate here.** The heavy tail
recorded in (realised sd 187.5 against a predicted 0.375) is the *ratio estimator*
``phi-hat``'s tail; its robust scale matched prediction to ~10 %. This is the same quasi-likelihood
`scripts/design/observable_efficiency.py` uses to produce every efficiency number in the design log — it
is being lifted into the solver, not invented.

⚠ **STRAND-AGNOSTIC.** Whether a read aligned ``+`` or ``−`` says nothing about whether its molecule was
gDNA or RNA, so both genome-strand columns are summed. The strand Beta-Binomial keeps its own columns and
is untouched by this module.

⭐ **WHY IT SPEAKS WHERE STRAND CANNOT.** The strand likelihood is rank-1 in the tilt ``theta``, so on an
AMBIG node — where ``theta`` is a free nuisance — its Schur complement on ``lambda`` is exactly zero. This
likelihood **does not depend on ``theta`` at all**, so that argument does not apply to it: it informs
``lambda`` directly on every node, stranded or not. That is the whole point of the step, and it is why
``tau_len`` must NOT be gated to single-strand nodes the way ``i_strand`` is.

⛔ **IT IS NOT A SUBSTITUTE FOR THE LENGTH MODELS.** ``m1``/``m2`` are functionals of the two FITTED pmfs,
so this channel inherits their error — it is a second *use* of the length model, not a way around it. The
model-free level channel is a different question (the ``1/A`` weight), and a
model-free channel provably carries zero composition information anyway.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .node_chain import EDGE, NODE, NodeChain

__all__ = [
    "LandedMoments",
    "contained_moments",
    "crossing_moments",
    "build_slot_moments",
    "length_loglik",
]


@dataclass(frozen=True, slots=True)
class LandedMoments:
    """The five moments of the opportunity-tilted length distribution ``g_c`` at each object.

    ``m1 = E[u]``, ``m2 = E[w]``, ``q1 = E[u²]``, ``q2 = E[w²]``, ``q12 = E[u·w]`` — everything the
    conditional mean and covariance of ``(Σu, Σw)`` need, and nothing else. ``eff`` is the tilt's own
    normaliser ``E_c[A]``, carried so a consumer can assert it against the divisor the solver used
    (they are the same quantity, and two implementations of one quantity is
    trap 27).

    All arrays are per object, or scalars broadcastable over objects.
    """

    m1: np.ndarray
    m2: np.ndarray
    q1: np.ndarray
    q2: np.ndarray
    q12: np.ndarray
    eff: np.ndarray


def _pmf_cumulants(fl_pmf: np.ndarray):
    """Cumulative sums of ``f(w)·w^k`` for ``k ∈ {−2,−1,0,1,2,3}`` — the whole of the node frame.

    ⚠ ``w = 0`` contributes 0 to every reciprocal sum: a zero-length fragment does not exist, and the
    pmf is 0 there in every real model. Guarding it here rather than trusting the input is what keeps a
    stray ``f(0) > 0`` from producing an infinity three call frames away.
    """
    p = np.asarray(fl_pmf, dtype=np.float64)
    total = float(p.sum())
    p = p / total if total > 0.0 else p
    w = np.arange(p.shape[0], dtype=np.float64)
    inv = np.zeros_like(w)
    np.divide(1.0, w, out=inv, where=w > 0.0)
    return (
        np.cumsum(p),  # F   = Σ f
        np.cumsum(p * inv),  # C1  = Σ f/w
        np.cumsum(p * inv * inv),  # C2  = Σ f/w²
        np.cumsum(p * w),  # S1  = Σ w f
        np.cumsum(p * w * w),  # S2  = Σ w² f
        np.cumsum(p * w * w * w),  # S3  = Σ w³ f
    )


def contained_moments(node_len_bp: np.ndarray, fl_pmf: np.ndarray) -> LandedMoments:
    """Moments of the tilted pmf for the CONTAINED population: ``A(w) = (ell − w + 1)+``, ``u(w) = 1/w``.

    ⭐ **Closed form, O(n_nodes).** Each raw moment is ``(ell+1)·<cumsum> − <next cumsum>``, exactly the
    shape `effective_length.contained_eff_length` uses for the denominator — so no
    ``n_nodes × max_len`` array is ever materialised (that would be 8 GB at human scale).

        E[A]     = (ell+1)·F  − S1        <- IS contained_eff_length
        E[A·u]   = (ell+1)·C1 − F
        E[A·w]   = (ell+1)·S1 − S2
        E[A·u²]  = (ell+1)·C2 − C1
        E[A·w²]  = (ell+1)·S2 − S3
        E[A·u·w] = (ell+1)·F  − S1        <- ⭐ identical to E[A], because u(w)·w = 1 at a node

    ⭐ That last identity means ``q12 ≡ 1`` at every node, for both components. It is not a shortcut: it
    says the two channels' cross-moment carries no composition information in the contained frame, and
    it falls out of the deposit rule rather than being imposed.
    """
    F, C1, C2, S1, S2, S3 = _pmf_cumulants(fl_pmf)
    n = F.shape[0]
    node = np.asarray(node_len_bp, dtype=np.float64)
    i = np.clip(np.floor(node).astype(np.int64), 0, n - 1)
    a = node + 1.0

    eff = np.maximum(a * F[i] - S1[i], 0.0)
    return _normalise(
        eff,
        a * C1[i] - F[i],
        a * S1[i] - S2[i],
        a * C2[i] - C1[i],
        a * S2[i] - S3[i],
        eff,  # E[A·u·w] == E[A] because u(w)·w == 1
    )


def crossing_moments(fl_pmf: np.ndarray) -> LandedMoments:
    """Moments for the CROSSING population at UNBOUNDED reach: ``A(w) = (w−1)+``, ``u(w) = 1/(w−1)``.

    Every entry is a scalar — under unbounded reach a line's opportunity does not depend on where it is,
    which is the "every edge has the same expectation" property stated in moments.

        E[A]     = mu − 1                    <- IS crossing_eff_length at UNBOUNDED_REACH
        E[A·u]   = P(w >= 2)
        E[A·w]   = E[w²] − mu
        E[A·u²]  = Σ f(w)/(w−1)
        E[A·w²]  = E[w³] − E[w²]
        E[A·u·w] = mu                        (u(w)·w = w/(w−1), so Σ f(w)(w−1)·w/(w−1) = mu)

    ⚠ **Unbounded reach only, matching `build_node_geometry`'s default.** With the TRAPS: prove-the-substrate taper switched on
    (``edge_rna_reach``) the opportunity becomes per-edge and these moments would have to as well. The
    taper was measured as a null (≤ 0.0002), so the default path is the one wired;
    a consumer that turns the taper on must extend this function rather than silently mismatch.
    """
    p = np.asarray(fl_pmf, dtype=np.float64)
    total = float(p.sum())
    p = p / total if total > 0.0 else p
    w = np.arange(p.shape[0], dtype=np.float64)
    ok = w >= 2.0  # a length-0 or length-1 fragment cannot cross a 0-bp line
    inv = np.zeros_like(w)
    np.divide(1.0, w - 1.0, out=inv, where=ok)

    eff = float((p * np.maximum(w - 1.0, 0.0)).sum())
    return _normalise(
        np.asarray(eff),
        np.asarray(float(p[ok].sum())),
        np.asarray(float((p * np.maximum(w - 1.0, 0.0) * w).sum())),
        np.asarray(float((p * inv)[ok].sum())),
        np.asarray(float((p * np.maximum(w - 1.0, 0.0) * w * w).sum())),
        np.asarray(float((p * w)[ok].sum())),
    )


def _normalise(eff, e_u, e_w, e_uu, e_ww, e_uw) -> LandedMoments:
    """Divide the raw ``E[A··]`` moments by ``E[A]`` to get the tilted-pmf moments.

    ⛔ Zero opportunity ⇒ every moment is 0, never a floored division.
    A slot with no opportunity for a component contributes nothing, and the caller's ``det > 0`` gate
    then makes the whole term inert there.
    """
    eff = np.asarray(eff, dtype=np.float64)
    live = eff > 0.0

    def d(x):
        x = np.broadcast_to(np.asarray(x, dtype=np.float64), eff.shape)
        return np.divide(x, eff, out=np.zeros(eff.shape, dtype=np.float64), where=live)

    return LandedMoments(m1=d(e_u), m2=d(e_w), q1=d(e_uu), q2=d(e_ww), q12=d(e_uw), eff=eff)


def build_slot_moments(chain: NodeChain, region_arrays, fl_pmf: np.ndarray) -> LandedMoments:
    """Scatter the two frames' moments onto the chain: contained at NODE slots, crossing at EDGE slots.

    The same slot layout `build_node_geometry` uses for ``eff_gdna``/``eff_rna``, so ``moments.eff`` is
    that array — asserted by ``test_eff_matches_the_solver_divisor``, because two implementations of one
    quantity is how a ½ went unnoticed for months.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_node, is_edge = kind == NODE, kind == EDGE
    n = int(chain.n_slots)

    node_len = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    node_m = contained_moments(node_len, fl_pmf) if node_len.shape[0] else None
    edge_m = crossing_moments(fl_pmf)

    fields = {}
    for name in ("m1", "m2", "q1", "q2", "q12", "eff"):
        out = np.zeros(n, dtype=np.float64)
        if node_m is not None:
            out[is_node] = getattr(node_m, name)[obj[is_node]]
        out[is_edge] = float(getattr(edge_m, name))
        fields[name] = out
    return LandedMoments(**fields)


def length_loglik(
    moments_g: LandedMoments,
    moments_r: LandedMoments,
    count: np.ndarray,
    inv_length_sum: np.ndarray,
    length_sum: np.ndarray,
    fg_grid: np.ndarray,
) -> np.ndarray:
    """The per-slot length log-likelihood on the ``σ(λ)`` grid → ``(n_slots, K)``.

    ``count`` / ``inv_length_sum`` / ``length_sum`` are per slot, already summed over the two genome-strand
    columns (this channel is strand-agnostic; see the module docstring).

    ⛔ **A slot contributes NOTHING — a flat row, not a small one — when it cannot speak**: no count, no
    opportunity for one of the two components, or a degenerate covariance. A flat row carries zero
    information through `density_deconv.density_factor_precision` (which tests ``ptp > 0``) and adds a
    constant to ψ, which normalises away. That is the honest encoding of "no evidence", and it matches
    `effective_length`'s own contract: never a floored value.

    ⚠ **The rows are NOT normalised.** ψ is a log-density up to an additive constant per node and is
    normalised once, at the end, by the caller's log-sum-exp — subtracting a per-row maximum here would
    be harmless but redundant, and doing it *twice* is how a term silently loses its scale.
    """
    fg = np.asarray(fg_grid, dtype=np.float64)[None, :]  # (1,K)
    n = np.asarray(count, dtype=np.float64)[:, None]  # (m,1)
    d_obs = np.asarray(inv_length_sum, dtype=np.float64)[:, None]
    s_obs = np.asarray(length_sum, dtype=np.float64)[:, None]

    def mix(g, r):
        return (
            fg * np.asarray(g, np.float64)[:, None]
            + (1.0 - fg) * np.asarray(r, np.float64)[:, None]
        )

    mu_d = mix(moments_g.m1, moments_r.m1)  # (m,K) per-fragment means
    mu_s = mix(moments_g.m2, moments_r.m2)
    # the mixture's second moments, then the central ones. Var over N i.i.d. draws is N times per-draw.
    v_dd = n * (mix(moments_g.q1, moments_r.q1) - mu_d * mu_d)
    v_ss = n * (mix(moments_g.q2, moments_r.q2) - mu_s * mu_s)
    v_ds = n * (mix(moments_g.q12, moments_r.q12) - mu_d * mu_s)

    r_d = d_obs - n * mu_d
    r_s = s_obs - n * mu_s
    det = v_dd * v_ss - v_ds * v_ds

    # live iff there is a count, BOTH components have opportunity here (else the mixture is not a
    # mixture and pi is unidentified), and the covariance is non-degenerate.
    # ⭐ **STRUCTURAL discrimination gate, and it must be EXACT.** If the two components' moments are
    # identical the row is constant in ``pi`` — but computed as a difference of large floats it comes out
    # flat only to ~1e-11, not to 0. That is not good enough: `density_factor_precision` tests
    # ``ptp > 1e-12``, so a 1e-11 row reads as LIVE and its near-uniform posterior returns
    # ``tau = 1/Var(uniform over the grid) = 0.029`` — **the grid's own width sold as evidence**, which is
    # the exact failure that function's docstring warns about. Measured: 689 slots, max tau 0.02902,
    # matching ``1/Var(lambda grid)`` to five digits. So the gate is on the MOMENTS, not on the output.
    discriminates = np.zeros(n.shape, dtype=bool)
    for a, b in (
        (moments_g.m1, moments_r.m1),
        (moments_g.m2, moments_r.m2),
        (moments_g.q1, moments_r.q1),
        (moments_g.q2, moments_r.q2),
        (moments_g.q12, moments_r.q12),
    ):
        discriminates |= (
            np.broadcast_to(np.asarray(a, np.float64), n.shape[:1])
            != np.broadcast_to(np.asarray(b, np.float64), n.shape[:1])
        )[:, None]

    live = (
        (n > 0.0)
        & discriminates
        & (np.asarray(moments_g.eff, np.float64)[:, None] > 0.0)
        & (np.asarray(moments_r.eff, np.float64)[:, None] > 0.0)
        & (det > 0.0)
        & (v_dd > 0.0)
        & (v_ss > 0.0)
    )
    safe_det = np.where(live, det, 1.0)
    quad = (v_ss * r_d * r_d - 2.0 * v_ds * r_d * r_s + v_dd * r_s * r_s) / safe_det
    ll = -0.5 * (quad + np.log(safe_det))
    ll = np.where(live, ll, 0.0)
    # ⚠ A row that is live at SOME grid points and not others would be discontinuous; require the whole
    # row to be live, so a slot either speaks everywhere on the grid or is inert everywhere.
    return np.where(live.all(axis=1, keepdims=True), ll, 0.0)
