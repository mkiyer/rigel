"""The RNA-anchored evidence factor — the RNA side of the unspliced count, anchored on quantities
hybrid capture cannot mis-scale (owner design 2026-08-24; round 2 same day, after the adversarial
review and the two-round prototype pricing, 2026-08-24).

**The idea in one paragraph.** ψ's hardest slots are the stage-0 substrate's claimed populations —
`solvable_exon` REGIONs and `ss_intron_boundary` BOUNDARYs — where the strand channel is often
dead. The gDNA side of those slots cannot be anchored (an enriched slot's expected gDNA rate needs
a per-exon, per-panel enrichment estimate; capture GATING is owner-refuted — capture is a
spectrum), but the RNA side can: at a COMPLETE flank every molecule overlapping the exon passes
through the flank's splice junctions, so the junctions' spliced flux — certified RNA by
construction, gDNA cannot splice — predicts the exon's contained RNA; and mature RNA cannot cross
an exon|intron boundary unspliced (it spliced), so a boundary's unspliced crossing is gDNA +
nascent only, and the adjacent intron's excess over the intergenic background predicts the
nascent. Enrichment cancels pair-locally because anchor and target share each exon's own probe
footprint (measured: the RNA-frame boundary→exon ratio is capture-invariant while the gDNA frame
shifts 2.9–3.4×).

**ROUND 2 — the three review-driven changes, each priced separately in the round-2 review:**

1. ⛔ **THE ROUTE-SUM POOLING (a bug fix).** A flank's junctions are DISJOINT ROUTES — each
   molecule crosses exactly one — so the flank's molecule rate is the SUM of per-route rates
   ``Σ_J flux_J / A_J``. Round 1 pooled ratio-of-sums (the opportunity-weighted MEAN), under-
   predicting a k-route exon ~k× (measured: multi-route exons read truth/prediction median
   2.0–2.3; 0.99–1.15 after the fix) and manufacturing a fake heavy-tail "transport dispersion"
   out of route-count asymmetry (the deep-flux dispersion floor HALVED under the fix). Pooling
   within one shared junction stays ratio-of-sums — those genuinely are one measurement.
2. **THE MARGINAL SCORING.** The factor is the intron factory's own pattern — a NegBinomial, the
   Poisson marginal of the Gamma flux posterior (``size = flux + ½``, the Jeffreys convention
   `fit_gdna_background` uses) — averaged by a small median-preserving quantile quadrature over
   the multiplicative transport scatter. A 3× miss then costs a few polite nats instead of the
   count-scale Gaussian's tens, which is what stopped the anchor bullying an accurate relay at
   clean libraries (the g00 residual). ⛔ The round-1 Gaussian's tail spend, and before it the
   NB ``size → 0`` trap (a vanishing size ASSERTS mass-at-zero rather than going flat), are the
   two recorded wrong shapes.
3. **THE NASCENT TERM, MARGINALIZED.** Nascent RNA never splices, so the flux cannot see it; the
   adjacent intron measures it as excess over background. ⛔ BOTH point-estimate forms are
   positively biased at nascent-free truth — the plug-in ``max(rate − bg, 0)`` by ~0.40σ and the
   truncated-posterior MEAN by ~0.56σ (Jensen: the posterior mean of a convex clamp exceeds the
   clamp of the posterior mean; measured as the round-1 arm-B g98 regression, claimed B
   20.3k → 53.4k). The honest treatment PLUGS NOTHING: quantile nodes of the truncated-excess
   posterior enter the quadrature, so a clean intron's null carries its atom at exactly zero and
   a noisy intron contributes a SPREAD, not a phantom. The excess counts convert to a rate with
   the intron's RNA opportunity (nascent carries the RNA length distribution) — identical on the
   equal-length panel BY DESIGN, correct on real data; the fl-gap side panels falsify it.

⛔ **The quadrature scores EXONS only.** The boundary factor keeps the guarded-Gaussian family:
its prediction is near-zero wherever capture enriches the boundary crossing against the intron
the anchor reads (the cliff), and a quadrature there asserts that near-zero with counting-only
width — priced as a 734 → 27.5k zero-control leak. The Gaussian's guarded fit reads the cliff as
an enormous spread and goes honestly flat. Gated both ways: the capture-cliff flat contract and
the builder-dispatch sentinel (the correct family once shipped as DEAD CODE behind green gates).

**The width** is per-slot counting (the NB size) ⊕ the transport scatter, estimated per library
by the PER-ROUTE two-flank disagreement (two complete flanks of one exon predict the same RNA;
gDNA-free at every gDNA level; deliberately CONSERVATIVE — no counting-noise subtraction, because
the campaign's measured lesson is that over-confidence is the killer and over-width is cheap)
max-combined with the guarded lower-quantile fit's spread. When nothing is estimable the factor
is flat — the conservative unknown-dispersion limit (all-16 panel conditions estimate fine; the
limit exists for tiny substrates).

**Round-2 price** (`residual_proto3.py`, vs the round-1 shipped factor, wins-retention checked):
`g00`-OFF relay 148.8k → 55.2k whole-library region error; `g05`-OFF silent 124.7k → 31.8k
(claimed exons 80.7k → 5.5k); the stress-nascent row's claimed exons restored to the pre-anchor
base (22.9k → 6.6k vs base 6.4k); every `g98`/`g50` win within noise.

⚠ **Recorded residual candidates, deliberately not in this round** (round-2 review items): single-complete-flank exons inherit a center bias the pooled fit cannot see;
the sj strand columns are summed rather than matched to the exon's strand; short-exon flank pairs
share fragments (correlated disagreement); the capture-ON +12 % mature offset stays in the width.

**CITIZENSHIP (owner ruling, 2026-08-25): THE ANCHOR IS A MESSAGE.** The flank's spliced-fragment
observation is a one-hop imputation into the exon, so it lives in the message framework, not the
local solve: `prepare_flux_evidence` packages the sender-side observations once per library,
`RelayPolicy` (behind its `certified_flux` switch) delivers the claim as `PsiMessage.lam_rows`,
and the backbone sums the rows into the FINAL solve only — never phase-A, never the own-evidence
precision. The silent policy carries nothing, which restores it as the measured control.
`build_rna_anchor_factor` remains as the composed arithmetic REFERENCE the stream's parity gate
holds the policy to. Gate: ``tests/calibration/test_rna_anchor.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln, logsumexp, ndtr, ndtri
from scipy.stats import gamma as _gamma_dist

from .effective_length import crossing_eff_length
from .region_chain import BOUNDARY, REGION
from .signature import coarse_type_array
from .structural_claims import build_structural_claims

_EPS = 1.0e-12

#: Minimum population for the dispersion/center estimators — below this a quantile of a handful of
#: points is not an estimate, and the factor must be WEAK rather than confidently mis-scaled (the
#: same population-minimum reasoning as `_MIN_TRAIN` in `calibrate`). ⚠ A named threshold.
_MIN_PAIRS = 20

#: The normal distribution's MAD-to-standard-deviation factor — mathematics, not a tunable.
_NORMAL_MAD = 0.6745

#: Standard-normal quantiles at 25 % and 10 % — mathematics, not tunables (`left_fit_center_spread`).
_Z25 = 0.6745
_Z10 = 1.2816

#: Quadrature resolution: scatter nodes (midpoint quantiles of the scatter law) × nascent nodes
#: (quantiles of the truncated-excess posterior). Resolution knobs in the same sense as the solve
#: grid's K — accuracy parameters, not model parameters; the s = 0 identity and the gates hold at
#: any setting.
_N_SCATTER_NODES = 9
_N_NASCENT_NODES = 3


# ── estimators ───────────────────────────────────────────────────────────────────────────────────


def left_tail_log_variance(observed, predicted) -> "float | None":
    """The center-free WIDTH fallback: the MAD of the negative residuals (gDNA can only inflate a
    residual upward). Runs whenever the joint center fit refuses — a refused center must never
    take the width down with it (measured: 207k → 402k at a zero control when it did)."""
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
    pair residuals ``log(observed / predicted)`` (gDNA only inflates a residual upward, so the
    lower quantiles are RNA-transport-only; two log-normal quantiles determine both parameters).

    ⛔ SELF-CONSISTENCY GUARD — the fit predicts its own negative-residual fraction Φ(−m/s) and
    REFUSES when the observed fraction disagrees beyond the binomial noise band: a gDNA-saturated
    tail (`g98`) would otherwise fit a spurious positive center and give back the high-gDNA win
    (measured, 1.20M → 1.68M on the release metric)."""
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
    p_pred = float(ndtr(-m / sd))
    p_obs = float(np.mean(r < 0.0))
    band = 2.0 * float(np.sqrt(max(p_pred * (1.0 - p_pred), _EPS) / r.size))
    if abs(p_obs - p_pred) > band:
        return None
    return m, float(sd * sd)


def route_pair_log_variance(rate_left, rate_right) -> "float | None":
    """The transport dispersion from exons whose two complete flanks disagree — computed on
    PER-ROUTE-SUMMED flank rates (route-count asymmetry no longer reads as dispersion), MAD-robust,
    halved (each flank carries half the pair disagreement). ⛔ Deliberately CONSERVATIVE: the
    counting noise is NOT subtracted, so this over- rather than under-states the width — the
    campaign's measured lesson is that under-width (over-confidence) is what does damage.
    ``None`` below the population minimum (⇒ the factor goes weak, never confident)."""
    rl = np.asarray(rate_left, np.float64)
    rr = np.asarray(rate_right, np.float64)
    ok = (rl > 0.0) & (rr > 0.0)
    if int(ok.sum()) < _MIN_PAIRS:
        return None
    d = np.log(rl[ok] / rr[ok])
    mad = float(np.median(np.abs(d - np.median(d)))) / _NORMAL_MAD
    return max(0.5 * mad * mad, 0.0)


# ── the marginal scoring ─────────────────────────────────────────────────────────────────────────


def _log_negbinom_rows(g, mu, size) -> np.ndarray:
    """Vectorized ``log NegBinom(g; mean mu, size r)`` for continuous ``g`` — the shipped
    `density_deconv._log_negbinom` parametrization with per-slot arrays; the identity is GATED
    elementwise against the shipped function, so the two cannot drift."""
    g = np.asarray(g, np.float64)
    mu = np.maximum(np.asarray(mu, np.float64), _EPS)
    r = np.maximum(np.asarray(size, np.float64), _EPS)
    rpm = r + mu
    return (
        gammaln(g + r) - gammaln(r) - gammaln(g + 1.0) + r * np.log(r / rpm) + g * np.log(mu / rpm)
    )


def _quadrature_rows(unspl, mu, size, *, scatter_log_variance, nascent_count_nodes, fg_grid):
    """The factor rows: the NegBinomial marginal averaged over the transport scatter
    (stratified equal-mass midpoint-quantile nodes, symmetric ⇒ median-preserving — a
    mean-preserving −s²/2 offset is review-refuted under quantile calibration) and, when given,
    over the nascent-excess posterior nodes. ``scatter_log_variance = 0``
    reduces exactly to the pure NegBinomial (gated). ``None`` scatter ⇒ flat rows (the
    conservative unknown-dispersion limit)."""
    C = np.asarray(unspl, np.float64)
    mu = np.asarray(mu, np.float64)
    size = np.asarray(size, np.float64)
    rna_share = np.clip(1.0 - np.asarray(fg_grid, np.float64), _EPS, 1.0)
    g = rna_share[None, :] * C[:, None]
    if scatter_log_variance is None:
        return np.zeros_like(g)
    s = float(np.sqrt(max(float(scatter_log_variance), 0.0)))
    if s == 0.0:
        scatter = np.ones(1)
        log_w = np.zeros(1)
    else:
        # STRATIFIED EQUAL-MASS nodes (midpoint quantiles of the scatter law, equal weights) — the
        # configuration the panel PRICED. ⚠ Gauss–Hermite was tried and measured WORSE end to end:
        # its weights concentrate centrally, making the factor bulk-tighter, and the bulk is where
        # the clean-library damage lives (claimed exons 5.5k → 45k at `g05`-OFF under GH). The
        # stratified rule under-weights the deep tail and softens the bulk — the measured-winning
        # trade; a genuinely heavier-tailed scatter law is the recorded refinement. Symmetric
        # nodes keep the mixture median at MU (a mean-preserving −s²/2 offset is review-refuted).
        u = (np.arange(_N_SCATTER_NODES) + 0.5) / _N_SCATTER_NODES
        scatter = np.exp(s * ndtri(u))
        log_w = np.full(_N_SCATTER_NODES, -np.log(_N_SCATTER_NODES))
    if nascent_count_nodes is None:
        nas = np.zeros((C.shape[0], 1))
    else:
        nas = np.asarray(nascent_count_nodes, np.float64)
    comps = []
    comp_logw = []
    for j, nd in enumerate(scatter):
        for q in range(nas.shape[1]):
            mu_jq = np.maximum(mu * nd + nas[:, q], _EPS)
            comps.append(_log_negbinom_rows(g, mu_jq[:, None], size[:, None]))
            comp_logw.append(log_w[j] - np.log(nas.shape[1]))
    stacked = np.stack(comps) + np.asarray(comp_logw)[:, None, None]
    rows = logsumexp(stacked, axis=0)
    return rows - rows.max(axis=1, keepdims=True)


def _gaussian_rows(unspl, mu, var, fg_grid) -> np.ndarray:
    """The boundary family's scoring — the round-1 count-scale Gaussian, kept because every
    pricing round ran the boundaries in this form and the round-2 quadrature boundary was
    measured to LEAK at a zero control (claimed B 734 → 27.5k at `g00`-ON: its NB size from the
    intron count is confident, and the boundary's nascent is capture-enriched relative to the
    intron's — the recorded cross-cliff residual — so a tight factor converts the bias to phantom
    gDNA where the estimator-widened Gaussian stays honest). A row with unknown or infinite
    variance is exactly zero."""
    rna_share = np.clip(1.0 - np.asarray(fg_grid, np.float64), _EPS, 1.0)[None, :]
    C = np.asarray(unspl, np.float64)[:, None]
    mu = np.asarray(mu, np.float64)[:, None]
    var = np.asarray(var, np.float64)[:, None]
    g = rna_share * C
    with np.errstate(invalid="ignore", divide="ignore"):
        rows = np.where(np.isfinite(var) & (var > 0.0), -0.5 * (g - mu) ** 2 / var, 0.0)
    return rows - rows.max(axis=1, keepdims=True)


def boundary_anchor_rows(
    *,
    unspl_boundary,
    unspl_intron,
    gdna_opportunity_intron,
    rna_opportunity_intron,
    rna_opportunity_boundary,
    background_rate,
    fg_grid,
    pair_log_variance=None,
) -> np.ndarray:
    """The boundary factor: expected nascent from the adjacent intron's excess over background —
    mature RNA cannot cross an exon|intron boundary unspliced, so the crossing is gDNA + nascent
    and the intron's excess predicts the nascent half. Round-1 form (guarded center fit, MAD
    fallback, max-combined widths, per-slot subtraction noise) with the round-2 DIVISOR SWAP: the
    excess converts to a rate with the intron's RNA opportunity (nascent carries the RNA length
    distribution) — a no-op on the equal-length panel, correct under a length gap."""
    C_b = np.asarray(unspl_boundary, np.float64)
    C_i = np.asarray(unspl_intron, np.float64)
    E_gi = np.maximum(np.asarray(gdna_opportunity_intron, np.float64), _EPS)
    E_ri = np.maximum(np.asarray(rna_opportunity_intron, np.float64), _EPS)
    E_rb = np.asarray(rna_opportunity_boundary, np.float64)
    rho_i = (C_i + 0.5) / E_gi  # Jeffreys-Poisson posterior mean, gDNA frame
    rho_nascent = np.maximum(rho_i - float(background_rate), 0.0) * (E_gi / E_ri)
    mu = rho_nascent * E_rb
    fit = left_fit_center_spread(C_b, mu)
    if fit is not None:
        mu = mu * float(np.exp(fit[0]))
    mad = left_tail_log_variance(C_b, mu)
    ests = [
        v for v in (pair_log_variance, fit[1] if fit is not None else None, mad) if v is not None
    ]
    V = max(ests) if ests else None
    pair_cv2 = float(np.expm1(V)) if V is not None else np.inf
    anchor_var = ((C_i + 0.5) / E_gi**2) * (E_gi / E_ri) ** 2 * E_rb**2
    var = mu + pair_cv2 * mu**2 + anchor_var
    return _gaussian_rows(C_b, mu, var, fg_grid)


def nascent_rate_nodes(*, counts, gdna_opportunity, rna_opportunity, background_rate):
    """Quantile nodes of the intron's TRUNCATED-EXCESS posterior, as RNA-frame rates ``(m, Q)``.

    The intron's total-rate posterior is Gamma(counts + ½, gDNA opportunity) in the gDNA frame
    (its gDNA half and the intergenic background share that frame and cancel); the nascent excess
    is ``max(rate − background, 0)``, whose posterior QUANTILES — not its mean — enter the
    quadrature: the null keeps its atom at exactly zero (both point-estimate forms are positively
    biased there; Jensen, and the measured round-1 g98 regression). The excess then converts to
    the RNA frame by the opportunity ratio, because nascent fragments carry the RNA length
    distribution — a no-op on the equal-length panel, correct under a length gap."""
    C = np.asarray(counts, np.float64)
    Eg = np.maximum(np.asarray(gdna_opportunity, np.float64), _EPS)
    Er = np.maximum(np.asarray(rna_opportunity, np.float64), _EPS)
    u = (np.arange(_N_NASCENT_NODES) + 0.5) / _N_NASCENT_NODES
    q = _gamma_dist.ppf(u[None, :], (C + 0.5)[:, None], scale=(1.0 / Eg)[:, None])
    excess = np.maximum(q - float(background_rate), 0.0)
    return excess * (Eg / Er)[:, None]


# ── the route table and the selection ────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class RouteTable:
    """Per-junction flux and crossing opportunity, indexed by the region each route serves.

    ``into[r]`` / ``outof[r]`` list the sj indices whose destination / source region is ``r`` —
    the routes through r's left / right flank. One row per distinct genomic junction (the sj
    axis), which is one route; the routes at a flank are DISJOINT (each molecule crosses exactly
    one), which is why their rates SUM."""

    flux: np.ndarray
    opportunity: np.ndarray
    into: dict
    outof: dict


def build_route_table(sj_geometry, substrate, rna_fl_pmf) -> RouteTable:
    """The per-route resolution the per-slot summed arrays cannot give: each junction's certified
    flux (genome-strand columns summed) and its crossing opportunity under the RNA pmf."""
    flux = np.asarray(substrate.sj.count, np.float64).sum(axis=1)
    opp = np.asarray(
        crossing_eff_length(rna_fl_pmf, sj_geometry.reach_lo, sj_geometry.reach_hi), np.float64
    )
    into: dict = {}
    outof: dict = {}
    src = np.asarray(sj_geometry.src_region, np.int64)
    dst = np.asarray(sj_geometry.dst_region, np.int64)
    for j in range(src.shape[0]):
        if opp[j] > 0.0:
            into.setdefault(int(dst[j]), []).append(j)
            outof.setdefault(int(src[j]), []).append(j)
    return RouteTable(flux=flux, opportunity=opp, into=into, outof=outof)


def _flank_route_rate(flux, opportunity) -> tuple:
    """One flank's molecule rate: the SUM of per-route Jeffreys rates over its disjoint routes
    (the ½ split across the flank's routes so the flank carries one Jeffreys pseudo-count total).
    Returns ``(rate, total flux)``; ``(0, 0)`` for a route-less flank."""
    flux = np.asarray(flux, np.float64)
    opp = np.asarray(opportunity, np.float64)
    if flux.size == 0:
        return 0.0, 0.0
    rate = float(np.sum((flux + 0.5 / flux.size) / np.maximum(opp, _EPS)))
    return rate, float(np.sum(flux))


def _selection(chain, statics, geometry, region_arrays, routes: "RouteTable | None" = None):
    """The two anchored populations and their per-slot inputs — ONE selection, shared by the
    builder and :func:`eligible_slots` so the two cannot drift apart. With ``routes`` the exon
    side carries route-summed flank rates and the per-flank pair rates for the disagreement
    estimator; without (eligibility checks) only the membership is derived."""
    claims = build_structural_claims(chain, statics)
    kind = np.asarray(chain.kind)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    idx_c = np.asarray(chain.obj_idx, np.int64)
    eff_r = np.asarray(geometry.eff_rna, np.float64)
    eff_g = np.asarray(geometry.eff_gdna, np.float64)
    eff_sj = np.asarray(geometry.eff_sj, np.float64).sum(axis=1)
    unspl = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)

    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    slot_type = rtype[np.clip(idx_c, 0, rtype.shape[0] - 1)]

    comp_l = np.asarray(claims.exon_flank_left_complete, bool)
    comp_r = np.asarray(claims.exon_flank_right_complete, bool)
    exon_sel = (comp_l | comp_r) & (unspl > 0.0) & (eff_r > 0.0)
    e_idx = np.flatnonzero(exon_sel)

    exon_rate = np.zeros(e_idx.size)
    exon_flux = np.zeros(e_idx.size)
    pair_l: list = []
    pair_r: list = []
    if routes is not None:
        for i, e in enumerate(e_idx):
            r = int(idx_c[e])
            flank_rates = []
            for cflag, js in ((comp_l[e], routes.into.get(r)), (comp_r[e], routes.outof.get(r))):
                if cflag and js:
                    rate, fx = _flank_route_rate(routes.flux[js], routes.opportunity[js])
                    if rate > 0.0:
                        flank_rates.append(rate)
                        exon_flux[i] += fx
            if flank_rates:
                exon_rate[i] = float(np.mean(flank_rates))
            if len(flank_rates) == 2:
                pair_l.append(flank_rates[0])
                pair_r.append(flank_rates[1])
        keep = exon_rate > 0.0
        e_idx, exon_rate, exon_flux = e_idx[keep], exon_rate[keep], exon_flux[keep]
    else:
        # membership only — a complete flank with any junction opportunity anchors the exon
        has = np.zeros(e_idx.size, bool)
        for i, e in enumerate(e_idx):
            for cflag, side in ((comp_l[e], left[e]), (comp_r[e], right[e])):
                if cflag and side >= 0 and eff_sj[side] > 0.0:
                    has[i] = True
        e_idx = e_idx[has]

    # nascent inputs per anchored exon: pooled adjacent introns beyond the flank boundaries
    nas_c = np.zeros(e_idx.size)
    nas_eg = np.zeros(e_idx.size)
    nas_er = np.zeros(e_idx.size)

    def _intron_beyond(bslot):
        for s_ in (left[bslot], right[bslot]):
            if s_ >= 0 and kind[s_] == REGION and slot_type[s_] == 1 and eff_g[s_] > 0.0:
                return int(s_)
        return -1

    for i, e in enumerate(e_idx):
        for side in (left[e], right[e]):
            if side >= 0:
                isl = _intron_beyond(side)
                if isl >= 0:
                    nas_c[i] += unspl[isl]
                    nas_eg[i] += eff_g[isl]
                    nas_er[i] += eff_r[isl]

    # boundaries: the claimed ss-intron boundaries whose intron flank is measurable
    intergenic = (kind == REGION) & (slot_type == 0) & (eff_g > _EPS)
    background_rate = float(unspl[intergenic].sum() / max(eff_g[intergenic].sum(), _EPS))
    b_sel = np.asarray(claims.ss_intron_boundary, bool) & (kind == BOUNDARY) & (unspl > 0.0)
    b_idx: list = []
    b_intron: list = []
    for b in np.flatnonzero(b_sel):
        isl = _intron_beyond(b)
        if isl >= 0 and eff_r[b] > 0.0:
            b_idx.append(int(b))
            b_intron.append(isl)

    return {
        "unspl": unspl,
        "eff_r": eff_r,
        "eff_g": eff_g,
        "exon_idx": e_idx,
        "exon_rate": exon_rate,
        "exon_flux": exon_flux,
        "exon_nascent_counts": nas_c,
        "exon_nascent_gdna_opp": nas_eg,
        "exon_nascent_rna_opp": nas_er,
        "pair_rate_left": np.asarray(pair_l, np.float64),
        "pair_rate_right": np.asarray(pair_r, np.float64),
        "boundary_idx": np.asarray(b_idx, np.int64),
        "boundary_intron": np.asarray(b_intron, np.int64),
        "background_rate": background_rate,
    }


def eligible_slots(chain, statics, geometry, region_arrays) -> np.ndarray:
    """Boolean (n_slots,): exactly the slots the factor may touch — the leakage gate's contract."""
    sel = _selection(chain, statics, geometry, region_arrays, routes=None)
    out = np.zeros(int(np.asarray(chain.kind).shape[0]), bool)
    out[sel["exon_idx"]] = True
    out[sel["boundary_idx"]] = True
    return out


def prepare_flux_evidence(chain, statics, geometry, region_arrays, routes: RouteTable):
    """The SENDER-SIDE half of the certified-flux stream (owner ruling 2026-08-25: the anchor IS
    a message): every grid-independent observation the claim is built from — the route-summed
    flank rates, the eligibility selections, the nascent frame arrays, the intergenic background —
    prepared once per library and handed to the message policy at construction. Belief-free by
    construction: counts, opportunities and structure only. Returns ``None`` when nothing is
    anchored (a toy with no complete flank), which the policy treats as a silent stream."""
    sel = _selection(chain, statics, geometry, region_arrays, routes)
    if sel["exon_idx"].size == 0 and sel["boundary_idx"].size == 0:
        return None
    sel["n_slots"] = int(np.asarray(chain.kind).shape[0])
    return sel


def flux_rows(evidence, *, n_grid: int, logodds_window: float) -> "np.ndarray | None":
    """The RECIPIENT-SIDE arithmetic of the certified-flux stream: the delivered observation
    scored on the solve grid — route-sum + the NB-marginal quadrature at claimed exons, the
    guarded-Gaussian family at eligible boundaries — as an ``(n_slots, K)`` λ-factor row array,
    zero at every unanchored slot. ``None`` evidence ⇒ ``None`` (no claim)."""
    from .simplex_logodds import _logodds_grid

    if evidence is None:
        return None
    sel = evidence
    _, fg = _logodds_grid(int(n_grid), float(logodds_window))
    out = np.zeros((int(sel["n_slots"]), int(n_grid)), np.float64)
    rho_bg = sel["background_rate"]

    V_pair = route_pair_log_variance(sel["pair_rate_left"], sel["pair_rate_right"])

    e_idx = sel["exon_idx"]
    if e_idx.size:
        mu = sel["exon_rate"] * sel["eff_r"][e_idx]
        C = sel["unspl"][e_idx]
        fit = left_fit_center_spread(C, mu)
        if fit is not None:
            mu = mu * float(np.exp(fit[0]))
        # width = the MAX over ALL available estimates — an accepted fit's sd must never suppress
        # a larger MAD (measured: at capture-ON the flank disagreement is common-mode blind,
        # sd 0.11 where the certified scatter is ~0.44, and letting it rule leaked 27.5k at a
        # zero control through the exon factor's spillover onto boundaries)
        mad = left_tail_log_variance(C, mu)
        ests = [v for v in (V_pair, fit[1] if fit is not None else None, mad) if v is not None]
        V = max(ests) if ests else None
        has_i = sel["exon_nascent_gdna_opp"] > 0.0
        nas = np.zeros((e_idx.size, _N_NASCENT_NODES))
        if has_i.any():
            nodes = nascent_rate_nodes(
                counts=sel["exon_nascent_counts"][has_i],
                gdna_opportunity=sel["exon_nascent_gdna_opp"][has_i],
                rna_opportunity=sel["exon_nascent_rna_opp"][has_i],
                background_rate=rho_bg,
            )
            nas[has_i] = nodes * sel["eff_r"][e_idx][has_i, None]
        out[e_idx] += _quadrature_rows(
            C,
            mu,
            sel["exon_flux"] + 0.5,
            scatter_log_variance=V,
            nascent_count_nodes=nas,
            fg_grid=fg,
        )

    b_idx = sel["boundary_idx"]
    if b_idx.size:
        b_int = sel["boundary_intron"]
        # the round-1 Gaussian family, NOT the quadrature: the quadrature boundary asserts a
        # near-zero prediction with counting-only width, and at a capture-ON zero-gDNA control the
        # boundary's nascent crossing is capture-enriched relative to the intron the prediction
        # reads — the priced leak was 734 → 27.5k. The Gaussian's guarded fit reads that cliff as
        # an enormous spread and goes honestly flat there, while staying live where the anchor
        # tracks (measured: this exact swap is what the interim pricing selected)
        out[b_idx] += boundary_anchor_rows(
            unspl_boundary=sel["unspl"][b_idx],
            unspl_intron=sel["unspl"][b_int],
            gdna_opportunity_intron=sel["eff_g"][b_int],
            rna_opportunity_intron=sel["eff_r"][b_int],
            rna_opportunity_boundary=sel["eff_r"][b_idx],
            background_rate=rho_bg,
            fg_grid=fg,
            pair_log_variance=V_pair,
        )
    return out


def build_rna_anchor_factor(
    chain,
    statics,
    geometry,
    region_arrays,
    routes: RouteTable,
    *,
    n_grid: int,
    logodds_window: float,
) -> "np.ndarray | None":
    """The two halves composed — the arithmetic REFERENCE the stream's parity gate holds the
    policy to. Production goes through the halves (`prepare_flux_evidence` at `calibrate`,
    `flux_rows` inside the relay's certified-flux stream); this composition exists so one
    function still states the whole claim."""
    return flux_rows(
        prepare_flux_evidence(chain, statics, geometry, region_arrays, routes),
        n_grid=n_grid,
        logodds_window=logodds_window,
    )
