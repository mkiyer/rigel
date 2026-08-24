"""The RNA-anchored evidence factor — the RNA side of the unspliced count, anchored on quantities
hybrid capture cannot mis-scale (owner design, 2026-08-24).

**The idea in one paragraph.** ψ's hardest slots are the stage-0 substrate's claimed populations —
`solvable_exon` REGIONs and `ss_intron_boundary` BOUNDARYs — where the strand channel is often dead
and, until 2026-08-24, a reference location decided the answer. The gDNA side of those slots cannot
be anchored (its expected rate under capture needs an enrichment estimate that varies per exon and
per panel), but the RNA side can: at a COMPLETE flank every molecule overlapping the exon passes
through the flank's splice junction, so the junction's spliced flux — certified RNA by construction,
gDNA cannot splice — predicts the exon's contained RNA; and mature RNA cannot cross an exon|intron
boundary unspliced (it spliced), so a boundary's unspliced crossing is gDNA + nascent only, and the
adjacent intron's excess over the intergenic background predicts the nascent. Each factor scores
``(1 − f_g) · unspliced count`` against its prediction, with a width carrying the prediction's
HONEST uncertainty. **No gDNA rate for an enriched slot appears anywhere**, which is
what makes the factor capture-robust: anchor and target share the same probe footprint, so
enrichment cancels pair-locally (measured: the RNA-frame boundary→exon ratio has the same
distribution capture-ON and capture-OFF, median 1.06–1.24, while the gDNA frame shifts 2.9–3.4×).

**The width is the load-bearing half, and it was paid for.** The first prototype carried only the
anchor's counting noise and failed both zero controls (a tight factor centred ~15 % low confidently
converts that residual into phantom gDNA). Three derived terms fix it, no constant introduced:

1. the PAIR-LEVEL dispersion of the transport. PRIMARY estimator: the TWO-FLANK DISAGREEMENT —
   an exon with two complete flanks has two independent junction estimates of the SAME RNA, so
   ``Var(log(rate_left / rate_right)) / 2``, MAD-robust and with the fluxes' own counting noise
   subtracted, measures the transport dispersion gDNA-FREE and capture-free at every condition
   (the same between-source-heterogeneity idea as the relay's DerSimonian–Laird law) — ⛔ but it is
   BLIND TO COMMON-MODE transport error (both flanks share the transport into the exon), and that
   under-read leaked 97k phantom gDNA at a zero-gDNA control when it ran alone. SECOND estimator:
   the NEGATIVE residuals of observed count vs prediction — an observation can exceed its RNA
   prediction because of gDNA but can only undershoot through dispersion — which sees the FULL
   spread but STARVES on gDNA-rich data (at `g98` no observation ever falls below its RNA
   prediction). The factor takes the MAX of whichever are estimable — conservative by
   construction, live at every condition one of them can see. 0.6745 is the normal-MAD constant,
   mathematics not tuning; the variance contribution is ``expm1(V)·mu²``;
2. the anchor's own counting noise — the Gamma-posterior variance of a Poisson rate under the same
   Jeffreys-½ convention `fit_gdna_background` uses; a zero-flux anchor therefore predicts the
   POSTERIOR MEAN ``½/opportunity``, never a hard zero;
3. at boundaries, the propagated noise of the (intron − background) subtraction — as the nascent
   excess falls inside its own counting noise the factor goes flat rather than asserting.

The factor is a COUNT-SCALE GAUSSIAN, ``−½·((1−f_g)·C − mu)² / var``, with
``var = mu + pair_var + anchor_var`` (Poisson + the two terms above). ⛔ The first build used a
NegBinomial with ``size = 1/CV²``, and its ``size → 0`` limit is a TRAP the gates caught before it
shipped: a NegBinomial with vanishing size is not "uninformative" — it concentrates all mass at
zero count, so exactly where the factor knew LEAST it asserted HARDEST (unknown dispersion read as
a confident gDNA call). The Gaussian's ``var → ∞`` limit is genuinely flat, and its confident-zero
limit (an empty prediction backed by a deep anchor) stays a tight deliver — the two limits the
NegBinomial parametrization could not separate.

Consumed by `calibrate` exactly as the intron factory's factor is: summed into the per-slot
λ-factor array, so it enters the local solve, the relay, and — via `density_factor_precision` —
the slot's own-evidence precision. Gate: ``tests/calibration/test_rna_anchor.py``. Priced on the
panel before landing (the panel-before-src rule): worst in-scope row (`g98` ss.99 capture-ON)
73k → 12k misplaced fragments, deferred stratum 1.78M → 166k, and the three message policies
converge once the slots carry their own evidence.

⚠ **The recorded residuals, so they are hunted deliberately rather than rediscovered**: (a) the
boundary factor leaks ~8k phantom gDNA at a zero-gDNA capture-ON library — the boundary's nascent
crossing straddles into the probed exon and is capture-enriched relative to the intron's contained
nascent, so the anchor under-predicts there (deliver is immune: zero × anything = zero); (b) exons
carry no nascent term, so at the panel's 20 % STRESS nascent share the unstranded capture-OFF exon
refute error regresses (6.4k → 56k; a stress reading under the NASCENT SCOPE RULING) — the
bias-corrected intron-nascent term is the measured-but-not-yet-clean fix (the naive form helped
there and regressed `g98`-ON exons through the positive bias of ``max(noisy − background, 0)``).
"""

from __future__ import annotations

import numpy as np
from scipy.special import ndtr

from .region_chain import BOUNDARY, REGION
from .signature import coarse_type_array
from .structural_claims import build_structural_claims

_EPS = 1.0e-12

#: Minimum anchored slots for the left-tail dispersion fit — below this the negative-residual MAD
#: is a quantile of a handful of points and the dispersion is unknowable, so the factor must be
#: WEAK rather than confidently mis-scaled (same population-minimum reasoning as `_MIN_TRAIN` in
#: `calibrate`; the gate asserts the weak limit). ⚠ A named threshold, flagged at review.
_MIN_PAIRS = 20

#: The normal distribution's MAD-to-standard-deviation factor — mathematics, not a tunable.
_NORMAL_MAD = 0.6745

#: Standard-normal quantiles at 25 % and 10 % — mathematics, not tunables. Two lower quantiles of
#: the pair residuals determine BOTH the transport's center and its spread (`left_fit`): gDNA can
#: only inflate a residual upward, so the lower quantiles are RNA-transport-only.
_Z25 = 0.6745
_Z10 = 1.2816


def flank_disagreement_log_variance(flux_a, opp_a, flux_b, opp_b) -> "float | None":
    """The transport dispersion from exons whose TWO complete flanks disagree about the same RNA.

    ``d = log(rate_a / rate_b)`` carries twice the per-flank transport dispersion plus both fluxes'
    counting noise; the MAD-robust variance of ``d`` halved, with the median counting contribution
    subtracted (floored at 0), estimates the per-flank dispersion — with no gDNA anywhere in the
    computation, so it cannot starve on contaminated data. ``None`` below the population minimum."""
    fa = np.asarray(flux_a, np.float64)
    fb = np.asarray(flux_b, np.float64)
    ok = (fa > 0.0) & (fb > 0.0)
    if int(ok.sum()) < _MIN_PAIRS:
        return None
    ra = (fa[ok] + 0.5) / np.maximum(np.asarray(opp_a, np.float64)[ok], _EPS)
    rb = (fb[ok] + 0.5) / np.maximum(np.asarray(opp_b, np.float64)[ok], _EPS)
    d = np.log(ra / rb)
    mad_var = (np.median(np.abs(d - np.median(d))) / _NORMAL_MAD) ** 2
    count_var = float(np.median(1.0 / (fa[ok] + 0.5) + 1.0 / (fb[ok] + 0.5)))
    return max(0.5 * (mad_var - count_var), 0.0)


def left_tail_log_variance(observed, predicted) -> "float | None":
    """The center-free WIDTH fallback: the MAD of the negative residuals (gDNA can only inflate a
    residual upward). Runs whenever the joint center fit refuses — a refused center must never
    take the width down with it (measured: losing this width under refusal over-tightened the
    factor at a zero control, 207k → 402k)."""
    obs = np.asarray(observed, np.float64)
    mu = np.asarray(predicted, np.float64)
    r = np.log(np.maximum(obs, _EPS) / np.maximum(mu, _EPS))
    neg = r[r < 0.0]
    if neg.size < _MIN_PAIRS:
        return None
    sd = float(np.median(np.abs(neg))) / _NORMAL_MAD
    return sd * sd


def left_fit_center_spread(observed, predicted) -> "tuple[float, float] | None":
    """The transport's CENTER and log-VARIANCE, fitted jointly from the two LOWER quantiles of the
    pair residuals ``log(observed / predicted)``.

    gDNA can only inflate a residual upward, so the lower quantiles are RNA-transport-only; under a
    log-normal transport, ``q25 = m − 0.6745·s`` and ``q10 = m − 1.2816·s`` determine both
    parameters. ⛔ The predecessor estimated the SPREAD alone and pinned the center at the raw
    prediction — and the transport's measured ~10–25 % systematic center (nascent invisible to the
    flux, frame effects) then converted into phantom gDNA at every RNA-rich exon, 15× a clean
    library's whole-library error. Returns ``(center_log, variance_log)``, or ``None`` when the
    left tail is too thin to support the 10 % quantile (⇒ the caller must go weak and center-free,
    never confident)."""
    obs = np.asarray(observed, np.float64)
    mu = np.asarray(predicted, np.float64)
    r = np.log(np.maximum(obs, _EPS) / np.maximum(mu, _EPS))
    if int((r < 0.0).sum()) < _MIN_PAIRS:
        return None
    q10, q25 = np.quantile(r, [0.10, 0.25])
    sd = (q25 - q10) / (_Z10 - _Z25)
    if not np.isfinite(sd) or sd <= 0.0:
        return None
    m = float(q25 + _Z25 * sd)
    # ⛔ SELF-CONSISTENCY GUARD — the fit's own validity condition, derived not tuned. The
    # two-quantile fit assumes the 10–25 % quantiles are RNA-transport-only; whether that held is
    # CHECKABLE: a log-normal (m, sd) predicts its own negative-residual fraction Φ(−m/sd), and
    # gDNA contamination of the fitting quantiles drives the observed fraction below it. Refuse
    # when they disagree beyond the binomial noise band (±2·√(p(1−p)/n)). Measured: without this
    # guard the center fit fires on `g98`'s gDNA-saturated tail with a spurious +m, inflating the
    # RNA prediction and giving back most of the high-gDNA win (release metric 1.20M → 1.68M);
    # with it, `g98` refuses (falling back to the flank-disagreement width, center-free) while the
    # clean-library fit stays accepted.
    p_pred = float(ndtr(-m / sd))
    p_obs = float(np.mean(r < 0.0))
    band = 2.0 * float(np.sqrt(max(p_pred * (1.0 - p_pred), _EPS) / r.size))
    if abs(p_obs - p_pred) > band:
        return None
    return m, float(sd * sd)


def _factor_rows(unspl, mu, var, fg_grid) -> np.ndarray:
    """``−½·((1 − f_g)·C − mu)² / var`` per slot over the grid, each row offset so its max is 0
    (an f-independent constant is irrelevant to ψ; the offset keeps the numbers scaled). A row
    whose variance is infinite (or not finite) is exactly zero — unknown means WEAK."""
    rna_share = np.clip(1.0 - np.asarray(fg_grid, np.float64), _EPS, 1.0)[None, :]
    C = np.asarray(unspl, np.float64)[:, None]
    mu = np.asarray(mu, np.float64)[:, None]
    var = np.asarray(var, np.float64)[:, None]
    g = rna_share * C
    with np.errstate(invalid="ignore", divide="ignore"):
        rows = np.where(np.isfinite(var) & (var > 0.0), -0.5 * (g - mu) ** 2 / var, 0.0)
    return rows - rows.max(axis=1, keepdims=True)


def exon_anchor_rows(
    *, unspl, flux, flux_opportunity, rna_opportunity, fg_grid, pair_log_variance=None
) -> np.ndarray:
    """The exon factor: expected contained RNA from the certified splice flux, per anchored exon.

    ``flux``/``flux_opportunity`` are pooled over the exon's COMPLETE flanks (each complete flank
    independently accounts every molecule, so pooling is averaging, never double counting). The
    prediction is the Jeffreys-Poisson posterior mean re-centered by the population's own fitted
    transport center (`left_fit_center_spread` — gDNA-free by the lower-quantile argument); the
    width is the larger of that fit's spread and the two-flank disagreement (each estimator is
    blind to what the other sees: the disagreement to common-mode transport error, the left fit to
    nothing but it starves on gDNA-rich data), plus the anchor's own counting noise."""
    flux = np.asarray(flux, np.float64)
    opp = np.asarray(flux_opportunity, np.float64)
    mu = ((flux + 0.5) / np.maximum(opp, _EPS)) * np.asarray(rna_opportunity, np.float64)
    fit = left_fit_center_spread(unspl, mu)
    if fit is not None:
        mu = mu * float(np.exp(fit[0]))
    width = fit[1] if fit is not None else left_tail_log_variance(unspl, mu)
    ests = [v for v in (pair_log_variance, width) if v is not None]
    V = max(ests) if ests else None
    pair_cv2 = float(np.expm1(V)) if V is not None else np.inf
    anchor_var = mu**2 / (flux + 0.5)  # Gamma-posterior variance of the flux rate, scaled
    var = mu + pair_cv2 * mu**2 + anchor_var  # Poisson + pair scatter + anchor noise
    return _factor_rows(unspl, mu, var, fg_grid)


def boundary_anchor_rows(
    *,
    unspl_boundary,
    unspl_intron,
    gdna_opportunity_intron,
    rna_opportunity_boundary,
    background_rate,
    fg_grid,
    pair_log_variance=None,
) -> np.ndarray:
    """The boundary factor: expected nascent from the adjacent intron's excess over background.

    Mature RNA cannot be in an unspliced exon|intron crossing (it spliced), so the crossing is
    gDNA + nascent and the intron's own measured excess predicts the nascent half. The width adds
    the per-slot propagated noise of the (intron − background) subtraction, which sends the factor
    flat as the excess falls inside its own counting noise."""
    C_b = np.asarray(unspl_boundary, np.float64)
    C_i = np.asarray(unspl_intron, np.float64)
    E_i = np.maximum(np.asarray(gdna_opportunity_intron, np.float64), _EPS)
    E_rb = np.asarray(rna_opportunity_boundary, np.float64)
    rho_i = (C_i + 0.5) / E_i  # Jeffreys-Poisson posterior mean
    rho_nascent = np.maximum(rho_i - float(background_rate), 0.0)
    mu = rho_nascent * E_rb
    fit = left_fit_center_spread(C_b[rho_nascent > 0], mu[rho_nascent > 0])
    if fit is not None:
        mu = mu * float(np.exp(fit[0]))
    width = (
        fit[1]
        if fit is not None
        else left_tail_log_variance(C_b[rho_nascent > 0], mu[rho_nascent > 0])
    )
    ests = [v for v in (pair_log_variance, width) if v is not None]
    V = max(ests) if ests else None
    pair_cv2 = float(np.expm1(V)) if V is not None else np.inf
    anchor_var = ((C_i + 0.5) / E_i**2) * E_rb**2  # propagated (intron − background) noise
    var = mu + pair_cv2 * mu**2 + anchor_var
    return _factor_rows(C_b, mu, var, fg_grid)


def _selection(chain, statics, geometry, region_arrays):
    """The two anchored populations and their per-slot inputs — ONE selection, shared by the
    builder and by :func:`eligible_slots` so the two cannot drift apart."""
    claims = build_structural_claims(chain, statics)
    kind = np.asarray(chain.kind)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    eff_r = np.asarray(geometry.eff_rna, np.float64)
    eff_g = np.asarray(geometry.eff_gdna, np.float64)
    eff_sj = np.asarray(geometry.eff_sj, np.float64).sum(axis=1)
    sj_cnt = np.asarray(geometry.sj_count, np.float64).sum(axis=1)
    unspl = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)

    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    slot_type = rtype[np.clip(np.asarray(chain.obj_idx, np.int64), 0, rtype.shape[0] - 1)]

    # exons: complete flanks with junction opportunity, flux pooled per exon; the per-side
    # rates of DOUBLE-complete exons feed the disagreement estimator
    comp_l = np.asarray(claims.exon_flank_left_complete, bool)
    comp_r = np.asarray(claims.exon_flank_right_complete, bool)
    exon_sel = (comp_l | comp_r) & (unspl > 0.0) & (eff_r > 0.0)
    e_idx = np.flatnonzero(exon_sel)
    flux = np.zeros(e_idx.size)
    opp = np.zeros(e_idx.size)
    pair_l = []
    pair_r = []
    for j, e in enumerate(e_idx):
        sides = []
        for cflag, side in ((comp_l[e], left[e]), (comp_r[e], right[e])):
            if cflag and side >= 0 and eff_sj[side] > 0.0:
                flux[j] += sj_cnt[side]
                opp[j] += eff_sj[side]
                sides.append(side)
        if len(sides) == 2:
            pair_l.append(sides[0])
            pair_r.append(sides[1])
    ok = opp > 0.0
    e_idx, flux, opp = e_idx[ok], flux[ok], opp[ok]
    pair_l = np.asarray(pair_l, np.int64)
    pair_r = np.asarray(pair_r, np.int64)

    # boundaries: the claimed ss-intron boundaries whose intron flank is measurable
    intergenic = (kind == REGION) & (slot_type == 0) & (eff_g > _EPS)
    background_rate = float(unspl[intergenic].sum() / max(eff_g[intergenic].sum(), _EPS))
    b_sel = np.asarray(claims.ss_intron_boundary, bool) & (kind == BOUNDARY) & (unspl > 0.0)
    b_idx = []
    b_intron = []
    for b in np.flatnonzero(b_sel):
        for side in (left[b], right[b]):
            if side >= 0 and kind[side] == REGION and slot_type[side] == 1 and eff_g[side] > 0.0:
                if eff_r[b] > 0.0:
                    b_idx.append(int(b))
                    b_intron.append(int(side))
                break
    b_idx = np.asarray(b_idx, np.int64)
    b_intron = np.asarray(b_intron, np.int64)
    return {
        "unspl": unspl,
        "eff_r": eff_r,
        "eff_g": eff_g,
        "exon_idx": e_idx,
        "exon_flux": flux,
        "exon_flux_opportunity": opp,
        "pair_flux_left": sj_cnt[pair_l],
        "pair_opp_left": eff_sj[pair_l],
        "pair_flux_right": sj_cnt[pair_r],
        "pair_opp_right": eff_sj[pair_r],
        "boundary_idx": b_idx,
        "boundary_intron": b_intron,
        "background_rate": background_rate,
    }


def eligible_slots(chain, statics, geometry, region_arrays) -> np.ndarray:
    """Boolean (n_slots,): exactly the slots the factor may touch — the leakage gate's contract."""
    sel = _selection(chain, statics, geometry, region_arrays)
    out = np.zeros(int(np.asarray(chain.kind).shape[0]), bool)
    out[sel["exon_idx"]] = True
    out[sel["boundary_idx"]] = True
    return out


def build_rna_anchor_factor(
    chain, statics, geometry, region_arrays, *, n_grid: int, logodds_window: float
) -> "np.ndarray | None":
    """The (n_slots, K) log-factor over the solve grid, zero at every unanchored slot; ``None``
    when nothing is anchored. Belief-free — counts, opportunities and structure only — so it is
    built once per grid and summed with the intron factory's factor."""
    from .simplex_logodds import _logodds_grid

    sel = _selection(chain, statics, geometry, region_arrays)
    if sel["exon_idx"].size == 0 and sel["boundary_idx"].size == 0:
        return None
    _, fg = _logodds_grid(int(n_grid), float(logodds_window))
    # ONE transport dispersion for both families — the two-flank disagreement, gDNA-free at every
    # condition; ``None`` (too few double-complete exons) falls back to each family's left tail.
    V = flank_disagreement_log_variance(
        sel["pair_flux_left"],
        sel["pair_opp_left"],
        sel["pair_flux_right"],
        sel["pair_opp_right"],
    )
    out = np.zeros((int(np.asarray(chain.kind).shape[0]), int(n_grid)), np.float64)
    if sel["exon_idx"].size:
        out[sel["exon_idx"]] += exon_anchor_rows(
            unspl=sel["unspl"][sel["exon_idx"]],
            flux=sel["exon_flux"],
            flux_opportunity=sel["exon_flux_opportunity"],
            rna_opportunity=sel["eff_r"][sel["exon_idx"]],
            fg_grid=fg,
            pair_log_variance=V,
        )
    if sel["boundary_idx"].size:
        out[sel["boundary_idx"]] += boundary_anchor_rows(
            unspl_boundary=sel["unspl"][sel["boundary_idx"]],
            unspl_intron=sel["unspl"][sel["boundary_intron"]],
            gdna_opportunity_intron=sel["eff_g"][sel["boundary_intron"]],
            rna_opportunity_boundary=sel["eff_r"][sel["boundary_idx"]],
            background_rate=sel["background_rate"],
            fg_grid=fg,
            pair_log_variance=V,
        )
    return out
