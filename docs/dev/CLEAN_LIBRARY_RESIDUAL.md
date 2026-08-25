# THE CLEAN-LIBRARY RESIDUAL — derivation, design, prototype plan (2026-08-24)

    ⚠ A DEV DOC, written for the owner to follow and contribute to. Plain terms first, math second.
    Every number here was measured this session from certified truth with no solver
    (`transport_decomposition.py` in the session scratchpad; re-derivable in minutes).

## 1. THE PROBLEM, IN ONE PARAGRAPH

The RNA anchor predicts "this exon's unspliced count should contain about MU fragments of RNA;
the rest is gDNA." On gDNA-rich libraries that call is right and huge (the g98 wins). On CLEAN
libraries the relay was already solving those exons almost perfectly, and the anchor — wrong by a
modest amount per slot but confident about it — FIGHTS the relay and loses us 130–280k fragments
per g00 row. The diagnostic signature: under silence the anchor barely changes g00 (194k → 178k);
under relay it wrecks it (1.3k → 50k). So the residual is not "the anchor is very wrong" — it is
"the anchor claims more confidence than it has, exactly where a better voice already exists."

## 2. WHAT WE MEASURED (certified truth, no solver, all 16 conditions)

Think of the anchor's error as three separate dials: the CENTER (is the prediction aimed right?),
the SCALE (how wide should its error bars be?), and the SHAPE (how fast should confidence die off
beyond the error bars?). We measured each:

**CENTER — small, and its parts are now identified.**
The prediction `flux-rate × exon RNA opportunity` against certified mature RNA is centered at
**1.05–1.08 capture-OFF** (nearly aimed right; the frames are consistent) and **1.17–1.19
capture-ON** (spliced fragments capture slightly differently than contained ones — a small
on-target-vs-on-target geometry term). Nascent RNA is INVISIBLE to the junction flux (it never
splices), and adding the adjacent intron's nascent rate closes most of the remaining capture-OFF
gap (median full-RNA ratio 1.21 → 1.03). ⚠ The median anchored exon has ZERO nascent even on the
stress panel — nascent concentrates in half the genes — so the nascent term must be a per-exon
measurement, never a global correction.

**SCALE — the honest error bar is ±40–65 % and it does NOT shrink with depth.**
Stratified by junction depth, the per-pair disagreement between prediction and certified truth
floors at |log ratio| ≈ 0.28 (OFF) / 0.44 (ON) for the deepest pairs, where counting noise alone
predicts 0.12. So beyond counting noise there is an intrinsic transport dispersion (isoform mix,
local geometry) of roughly sd 0.4–0.65 in log units. The anchor can never be a precise per-slot
oracle at RNA-rich slots; pretending ±20 % when the truth is ±50 % is exactly how it out-shouts
the relay.

**SHAPE — the tails are heavier than log-normal, at every depth.**
The p90/p50 ratio of the |log| residuals is ~3.6 where a log-normal predicts 2.4: a minority of
pairs genuinely disagree 3–4× even at high depth. Under today's Gaussian factor those tail slots
get punished with tens of nats — they are the slots that bully the relay. A correct shape spends
a few nats there and stops.

## 3. THE DESIGN — three changes, each tied to one measured dial

**① CENTER, derived not estimated:** `MU = flux-rate × RNA-opportunity + nascent-term`, where the
nascent term is the adjacent intron's excess over background — but as the POSTERIOR MEAN of the
truncated rate, `E[max(rho − background, 0)]` under the intron's Gamma posterior (a two-line
closed form via the Gamma CDF), not the plug-in `max(estimate − background, 0)`. The plug-in has a
positive bias at nascent-free genes (half the noise survives the clamp) — that bias is what made
the v3 nascent experiment regress g98 exons. The posterior mean is unbiased-in-the-right-way and
its variance flows into the width. The capture-ON +12 % geometry term is left in the width (it is
small against sd 0.5) — a later refinement can derive it from `capture_eff_length`.

**② SCALE, per-slot counting ⊕ per-library intrinsic:** width² = counting noise of the anchor
(known analytically: flux Gamma + nascent truncated-posterior variance) ⊕ the intrinsic transport
dispersion, estimated per library from the TWO-FLANK DISAGREEMENT (two complete flanks of one
exon predict the same RNA; their disagreement is gDNA-free at every gDNA level) with the existing
self-consistency-guarded left fit as the second view. This is the machinery we already have — the
measurement above says its target value is ~0.4–0.65 and stable across strata.

**③ SHAPE, the actual fix for the clean-library damage:** replace the count-scale Gaussian with
the MARGINAL the codebase already ships for the intron factory, generalized one step: the intron
factory scores counts with NegBinomial = Poisson marginalized over a Gamma rate posterior. Our
rate posterior has one extra ingredient — the multiplicative transport scatter — so the marginal
is a small QUADRATURE (5–9 fixed nodes) of that same NegBinomial over the scatter:

    log P(g) = logsumexp_j [ log w_j + log NB(g ; MU·exp(s·x_j − s²/2), size = flux + ½) ]

with the node set drawn heavy-tailed (a scale-mixture: the scatter's own scale varies per pair,
which is what the measured p90/p50 = 3.6 says; the tail weight is derived from that measurable
quantile ratio of the two-flank disagreement — per library, guarded, no constant). Limits, all
structural: s → 0 gives the intron factory's exact NB; s unknown gives a weak factor; flux ≈ 0
gives the confident deliver; a 3× outlier costs a few nats instead of forty.

**In one sentence for the record: aim the prediction with derived components, put the measured
±50 % in the error bar, and make disbelief beyond the error bar polite rather than fanatical.**

## 4. WHY THIS SHOULD FIX THE RESIDUAL WITHOUT TOUCHING THE WINS

The g98 deliver mechanism needs only `MU ≈ 0 ⇒ the count is gDNA` — unchanged under ①–③ (NB at
tiny mean is decisively concentrated whatever the tails). The clean-library damage lives in the
tail slots where the anchor's prediction misses by 2–4× and today's Gaussian burns tens of nats
against the relay's correct answer; ③ caps that spend at a few nats, ① removes the part of the
miss that was systematic, ② stops under-stating the bar. The zero controls and the stress-nascent
row are the falsification set, plus g98 must stay within noise of today's wins.

## 5. PROTOTYPE PLAN (monkeypatch first, panel-before-src)

Arms, priced on: both zero controls, `g00`/`g05`-OFF (the residual), `g98`/`g50` ss.99 ON (the
wins), `g50`-OFF unstranded (stress nascent), release metric + claimed populations:

    A  today (Gaussian + guarded center fit)          — the baseline
    B  ① + ② only (Gaussian kept)                     — prices the center/scale half
    C  ① + ② + ③ (NB-quadrature, log-normal nodes)    — prices the marginal
    D  C with heavy-tailed nodes                       — prices the shape's tail weight

One thing varied per step; every arm at shipped refits, silent + relay (the residual is a relay
interaction). The winning arm graduates to `src/` with the gate set: the existing rna_anchor
gates re-derived for the new form + a tail-cost gate (a 3× outlier's penalty is bounded) + the
NB-limit identity (s = 0 reproduces `_log_negbinom` exactly).

## 5b. THE ADVERSARIAL REVIEW (four independent lenses) — two refutations, both confirmed

The design above was put through four independent adversarial verifiers before any prototype
graduated. Verdicts that CHANGED the design:

**REFUTED ①: the ROUTE-POOLING BUG (frames lens; empirically confirmed, `route_sum_check.py`).**
The shipped anchor pools a flank's junctions as ratio-of-sums — the opportunity-weighted MEAN of
per-route rates — but the routes are DISJOINT (each molecule crosses exactly one junction at a
flank), so the flank's molecule rate is the SUM of per-route rates. At a k-route flank the
prediction is deterministically ~k× low. Measured: multi-route exons (~26 % of anchored) read
median truth/prediction **2.0–2.3 under the shipped pooling and 0.99–1.15 under the route sum**;
single-route exons are identical under both; the deep-flux dispersion floor **halves**
(|log r| p50 0.284 → 0.178, p90 1.04 → 0.62). So roughly half the "intrinsic transport
dispersion" and most of the heavy-tail signature were the estimator's own pooling convention —
the SUM fix shrinks all three dials at once and goes in FIRST, before any shape machinery.

**REFUTED ②: the nascent "posterior mean removes the plug-in bias" claim (two lenses, same Jensen
argument).** E[max(rho−b,0)] ≥ max(E[rho]−b,0) pointwise — the posterior mean of the clamp is
MORE positively biased at nascent-free truth (0.56σ vs the plug-in's 0.40σ), not less. The
truncated-moment closed forms are correct and keep their width role, but the CENTER needs mass at
"no excess": a spike-and-slab weight (the posterior probability of a real excess against the
background-only model — derivable, no constant), which is also the NASCENT SCOPE RULING's own
default stated as an estimator. My v3-regression causal story was wrong as written.

**Corrections adopted from WEAKENED verdicts:** the −s²/2 mean-preserving offset is wrong under
median calibration (and undefined under t tails) — nodes are median-preserving `MU·exp(s·x_j)`
(the prototype already did this; the doc was wrong); the nascent divisor must be the intron's RNA
opportunity, not its gDNA opportunity (identical on the equal-length panel BY DESIGN, wrong on
real data — falsify on the fl-gap side panels); when the scatter scale is unknown the factor must
NOT flatten — the MU≈0 deliver is robust to ANY scatter (no multiplicative scatter un-zeroes a
zero prediction), so unknown-s should keep the deliver and weaken only the positive-MU rows; the
two-flank estimator needs per-route rates before differencing, shared-fragment pairs (short
exons) excluded, and the strand column matched to the exon rather than summed.

**Survived attack:** the truncated-Gamma closed forms; the flux≈0 deliver (20–30 nats of margin,
robust to arbitrarily large s); C=0 slots are structurally flat (the factor can only re-label
fragments that exist); the g98 mature/flux=1.4 reading as zero-truncation selection (with a named
stratified check).

## 5c. ROUND-2 PRICING — the review-corrected arms, measured (residual_proto3.py)

Arms: A = shipped anchor; R = A + the route-sum fix ONLY; RQ = R + the NB-quadrature shape
(median-preserving); RQN = RQ + the nascent term MARGINALIZED into the quadrature (truncated-
posterior quantile nodes — no plug-in, no spike weight needed). Whole-library abs err
(region/boundary) and claimed populations, shipped refits:

| row | pol | A | RQN | reading |
|---|---|---|---|---|
| g00-OFF (worst residual) | relay | 148.8k / 198.1k | **55.2k / 87.7k** | 2.7×/2.3×; claimed E 50.2k → 11.0k |
| g05-OFF | silent | 124.7k / 38.5k | **31.8k / 21.4k** | 3.9×; claimed E 80.7k → **5.5k** (14×) |
| g00-ON (zero ctl) | silent | 66.3k / 59.8k | **34.8k / 26.5k** | ~2× |
| g50-OFF (stress nascent) | silent | claimed E 22.9k | **6.6k** | ≈ the pre-anchor base (6.4k) — the stress row is FIXED by the marginalized nascent |
| g98/g50 stranded-ON + g98-OFF (the wins) | both | — | **within noise of A** | every win retained |

Attribution: the ROUTE-SUM fix is the largest single lever (A → R covers most of the g00/g05
distance); the quadrature shape adds a little (as expected once the tails were shown to be ~half
pooling artifact); the MARGINALIZED NASCENT is the second-largest lever at capture-OFF and is what
restores the stress-nascent row — while the empirically-confirmed-biased plug-in version (round 1
arm B) had regressed g98-ON claimed B 20.3k → 53.4k, exactly as the Jensen refutation predicted.

Remaining gap at g00-OFF relay vs the pre-anchor base: 55.2k vs 17k — down from 8.7× to 3.2×,
not yet closure; the review's unaddressed items (single-flank center bias, strand-column match,
shared-fragment pairs) are the named candidates for the rest.

**Ship list for the src graduation (next step):** the route-sum pooling in `_selection` (a bug
fix), the median-preserving NB-quadrature scoring, the marginalized nascent (quantile nodes), the
per-route two-flank estimator, the nascent divisor swap (intron RNA opportunity — panel-invisible,
real-data-correct), and the unknown-scatter deliver preservation. Gates: the existing anchor gate
set re-derived, the s = 0 ⇒ NB identity, a bounded-tail-cost gate, and the round-2 table as the
falsification set (wins within noise, residual rows ≤ the RQN readings, stress row ≤ base).

## 5d. THE SRC GRADUATION — measured, ACCEPTED (src_price_final3, 2026-08-24)

The shipped implementation reproduces every RQN falsification row to ~4 significant figures:
g00-OFF relay **55,379 / 88,345** (claimed E 11,022); g05-OFF silent **31,871 / 21,510** (claimed
E 5,539); g00-ON zero control silent **34,874 / 26,619** — claimed B **718** against the anchor-OFF
base 710; stress-nascent claimed E **6,645** vs base 6,368; g98-ON silent 80.7k → 68.5k lib-region
with claimed E 5,698 (relay 182.2k → 145.9k), g98-OFF claimed B 40.9k → 21.5k. Suite
3,660 / 8 xfail, all anchor-gate perturbations firing.

⛔ **The incident that cost a pricing round, kept as its gate's reason-to-exist:** two departures
from the priced RQN configuration were caught ONLY by re-pricing — Gauss–Hermite scatter nodes
(bulk-tight: g05-OFF claimed E 5.5k → 45k) and the QUADRATURE BOUNDARY (asserts a near-zero
nascent prediction with counting-only width; at a capture-ON zero-gDNA control the boundary's
crossing is capture-ENRICHED relative to the intron the prediction reads — claimed B 734 → 27.5k).
The boundary revert was then written as a standalone function but the builder kept dispatching to
the quadrature — DEAD CODE behind 17 green gates, found because the re-price still read 27.5k.
Both are now gated in `test_rna_anchor.py`: the capture-cliff flat contract (with the quadrature's
assertion measured in-gate so the flat assertion cannot go vacuous) and the builder-dispatch
sentinel, each verified firing on the exact shipped bug. The quadrature stays correct for EXONS
(its prediction there is flux-anchored, never near-zero); the boundary keeps the round-1 Gaussian
family whose guarded fit reads the cliff as enormous spread and goes honestly flat at ON — the
boundary factor's ON-condition contribution is therefore ~nil BY DESIGN, and the g98-ON boundary
wins come through the exon factor plus messages.

## 6. OPEN QUESTIONS FOR THE OWNER

* The intrinsic dispersion is real physics (isoform mix at a complete flank still permits ±40 %
  per-exon flux-vs-contained disagreement — depth does not remove it). Worth a deeper look later:
  if the per-exon isoform structure predicts part of it, the opportunity model could absorb it.
* The capture-ON +12 % mature offset: leave in the width now, or derive from `capture_eff_length`
  in this pass?
* The nascent term reaches only exons with a measurable adjacent intron; single-intron-side exons
  inherit the flank we have. Acceptable, or pool gene-wide?
