# The variance ledger — every term in a message's precision, what it prices, and why none of them overlap

**Owner-directed 2026-07-26**, in response to: *"We need to be careful about introducing variance metrics that
overlap. First it leads to over-engineering and second it may harm performance by compounding variance
measures. Do these new variance sources replace existing ones?"*

This is the single place that answers that. **Update it whenever a variance term is added, removed or
re-scoped** — an unaudited addition is exactly the compounding risk the question is about.

---

## 1. The ledger

A message's precision is `p = 1/Var(log ρ_c^msg)`. Every term below is an *additive* contribution to that
log-variance except where noted.

| # | term | what it prices | where it applies | status |
|---|---|---|---|---|
| M1 | `Var(log f_c^src)` | the source's own composition uncertainty (the `τ_λ` Schur Jacobian) | every message | live |
| M1 | `1/n_src` | the source's Poisson count | every message | live |
| M5 | `σ²_transfer = Var(log r)` | the reframe's **scale** sampling | every message **except the graft** (there `r` is common-mode across the matched set and cancels) | live |
| M8 | `(log r)²` (`graft_frame_logvar`) | the graft's **frame mislift**: the measured spliced already sits in the destination exon's frame, so `r` is *not* common-mode for it | the grafted spliced component only; **identically 0 at `r = 1`** | live |
| M10 | `w_μ²(v_ν + v_μ)` (`peel_share_logvar`) | the **peel share**: what fraction of the exon's RNA continues past the seam | the peel (exon → boundary) only | live |
| M11 | `ψ'(k) + 1/(n·E[f_R]²)` (`residual_level`) | the **level** the peel share is formed from | inside M10's `v_ν` | live |
| M7 | `b̂²` (`mismatch_deflate`) | the imputation **premise** being false — DerSimonian–Laird against the destination's own self-solve | every stream, **but inert wherever `τ_own = 0`** (all AMBIG, all unstranded ⇒ 84 % of error mass) | live |
| — | `1/n_spl` | the spliced measurement's own count | the graft | live |
| **P1d** | **the graft's extrapolation** | a junction measured over ~100 bp used as the RNA density of a ~2,100 bp exon | the grafted spliced component only | **NOT BUILT** — §3 |
| P1e | the conservation **surprise** `δ²/(αᵀΣα)` | how implausible the claim is against the node's own observed mass | every message | not built |
| — | ~~NPMLE projection `var_proj + (Δμ_proj)²`~~ | a **density** disagreement between two nodes | — | **RETIRED**; `DensityNPMLE.project` is called nowhere. ⚠ `CalibrationConfig`'s docstring still describes it and is stale |

## 2. ⭐ Why they do not compound — the two rules

### 2.1 The additive terms are disjoint by construction

M1 / M5 / M8 / M10 / `1/n_spl` price **different physical events** and, where two could plausibly overlap, one
is switched off structurally rather than by judgement:

* **M5 vs M8** — both concern the reframe, and they are mutually exclusive by design: M5 is set to **0 on the
  graft** (`transfer_logvar`'s `graft` argument), because there `r` cancels; M8 exists *only* on the graft, to
  price the one component for which it does not. Neither fires where the other does.
* **M10 vs M8** — different edges. The peel is exon → boundary; the graft is boundary → exon.
* **M11 vs M1** — M11 is the level *inside* M10's share, not a second opinion about the same density.
* **`1/n_spl` vs P1d** — sampling versus premise. Charging a measurement its Poisson count says "I counted
  this many"; charging it an extrapolation says "and this many are not necessarily representative of the
  region I am applying them to". A count cannot see the second and never will, however large it is.

### 2.2 M7's DerSimonian–Laird **cannot** double-count with anything — this is provable

The deflated precision is

```
    b̂² = max(0, G² − v_msg − v_own)     ⇒     p_eff = 1 / max( v_msg , G² − v_own )
```

The `max` is the whole point. **Adding any new term to `v_msg` reduces `b̂²` one-for-one until `b̂²` hits its
floor at 0, after which `p_eff = 1/v_msg`.** So a newly-modelled variance never stacks on top of M7 — it
*replaces* the part of the gap M7 was previously absorbing as unexplained drift. M7 is the residual, by
construction, and it shrinks exactly as the model improves.

That is the direct answer to *"do these new variance sources replace existing ones?"*: **they do not replace
any named term, and they cannot compound with M7. The only place compounding is possible is among the additive
terms M1/M5/M8/M10, and §2.1 shows each pair is either structurally exclusive or physically disjoint.** The
discipline to keep is: **before adding a term, name the physical event it prices and show which existing term
would otherwise have absorbed it.**

## 3. P1d — the one gap the ledger exposes, and it is large

**The event:** the graft hands an exon a spliced density `ρ_μ = S/E_spl` measured over the junction's
~100 bp window and uses it as a claim about the exon's RNA density over ~2,100 bp — a **12–21×
extrapolation**. Nothing in the ledger prices that.

The one term that comes close is **M8**, and it does not cover it: M8 charges `(log r)²`, i.e. it assumes the
*only* reason a junction and an exon differ is the capture step between them. **Off-capture `r = 1` and M8
charges exactly zero — but a point and a region still differ**, because RNA density varies along an exon
(nascent gradients, alternative ends, isoform structure). That is a spatial effect, not a capture effect.

**Measured** (`scratchpad/x1_graft.py`; the graft share
`φ = [ρ_ν(bnd) + ρ_μ(bnd)] / ρ_R(exon)` in a common frame, using the ORACLE gDNA density ratio as the capture
step, since gDNA content is uniform along the genome):

| condition | edges | φ p25 | φ med | φ p75 | **sd(log φ)** |
|---|---|---|---|---|---|
| `gdna300 ss0.99 present capOFF` | 594 | 0.934 | 1.108 | 1.306 | **0.764** |
| `gdna300 ss0.50 present capOFF` | 599 | 0.920 | 1.126 | 1.354 | **0.758** |
| `gdna300 ss0.99 present capON` | 458 | 1.174 | 2.452 | 5.056 | 1.698 |
| `gdna300 ss0.50 none capON` | 449 | 0.944 | 2.454 | 5.822 | 1.873 |
| `gdna100 verystrong` | 237 | 3.392 | 5.545 | 9.111 | 1.261 |

Read the **capture-OFF rows** — that is where M8 is identically zero and the term would be doing all the work
alone. The bias is small (median φ = 1.11, i.e. +0.10 nats) but the **spread is 0.76 nats ⇒ a variance of
0.58**, against a count variance the graft is actually charged of as little as `1/2025 ≈ 0.0005`. **The graft
is roughly a thousand times more confident than it has earned.** Under capture the spread is larger still, but
there it is conflated with the frame step M8 already prices — which is exactly why any P1d term must be
**capture-independent**, or it double-counts M8.

## 4. What pricing it is worth — measured, with an oracle-free probe

`RIGEL_XVAR` adds a flat extra variance to the grafted spliced component (a **probe to size the prize**, not a
candidate: a flat constant is not landable). At 0.3:

| | suite error | **confidently-wrong** | z2\|Q1 all | **exon-single** | exon-AMBIG |
|---|---|---|---|---|---|
| landed (`pk2`) | 11,928,101 | 1,186,552 (9.9 %) | 15.5 | **8.4** | 92.1 |
| + extrapolation 0.3 | 12,216,265 | **689,639 (5.6 %)** | **8.2** | **1.2** | **62.6** |

**Single-strand exons go from 8.4× over-confident to essentially HONEST (z2 = 1.2), and 42 % of all
confidently-wrong mass disappears, for 2.4 % more error mass.** On mwae it costs (0.0855 → 0.0875), which is
the expected trade — damping a channel always costs mass when its mode is doing work — but mwae is not what
pass-0 owes Phase 2. This is the single largest trust improvement available anywhere on the list.

## 5. ⛔→✅ THE ESTIMATOR EXISTS. The earlier refutation was MY STATISTICAL ERROR.

### 5.1 The correction

An exon's two flanking junctions each claim its RNA density, so under the graft's premise they should agree,
and their disagreement is a free two-study estimate of the premise error. I first reported this as **refuted**
— "under-states by 20×, the two junctions are ~0.95 correlated". **That conclusion was wrong.** It compared
the *median* |gap| (0.14 — a robust statistic) against the full *sd* of log φ (0.76 — a non-robust one) on a
distribution whose variance is 7–9× concentrated in its tails. Apples to oranges.

With **matched statistics** (variance against variance):

| condition | `Var(log φ)` (truth, oracle) | `Var(gap)` | `Var(gap)/2` (the estimate) | ratio |
|---|---|---|---|---|
| `gdna300 ss0.99 present capOFF` | 0.583 | 1.287 | **0.644** | **0.91×** |
| `gdna300 ss0.50 present capOFF` | 0.574 | 1.319 | **0.659** | **0.87×** |

The two junctions are **essentially independent**, and `Var(gap)/2` recovers the target to within 9–13 %.
Subtracting the Poisson part gives the method-of-moments form the calibrator already uses twice
(`gdna_strand_overdispersion`, `rna_strand_overdispersion`):

```
    ω_graft = max( 0 , Var(gap) − E[ 1/n_left + 1/n_right ] ) / 2
            = 0.594  and  0.601        against a true Var(log φ) of 0.583 and 0.574   —  within 2 %
```

**Prior-free, no constant, and it lands on the magnitude the probe showed is worth 42 % of the
confidently-wrong mass.** (This does not resurrect the *left/right message gap* refuted in
`SESSION_2026_07_25_HANDOFF_10.md` §4 — that was a different statistic, comparing two fused BP messages that
share a reframe. These are two independent spliced COUNTS at two different loci, each already in the
destination exon's frame by M8's argument.)

### 5.2 ⭐ The open piece is now the SHAPE, not the magnitude

A single pooled ω would be a crude first cut, because the dispersion is strongly count-dependent — and **not
as `1/n`**:

| junction spliced count | `Var(log φ)` | `1/n` (Poisson) | excess |
|---|---|---|---|
| < 30 | 1.799 | 0.125 | **1.674** |
| 30–100 | 0.451 | 0.015 | 0.436 |
| 100–300 | 0.558 | 0.006 | 0.552 |
| 300–1000 | 0.145 | 0.002 | 0.143 |
| > 1000 | 0.042 | 0.0004 | **0.041** |

A **40× range** in the premise error across the count range, with the Poisson term explaining almost none of
it. A pooled ω = 0.59 would over-charge the well-counted junctions by 14× and under-charge the sparse ones by
3×.

The physical reading, and it is the thing to confirm against real data: **the count is a proxy for the
junction's SHARE of its exon's RNA.** A junction carrying most of the exon's transcripts is representative
(φ → 1, small variance); a minor-isoform junction carries little and represents little (φ ≪ 1, and `log φ` is
volatile). If that is the mechanism, the right per-junction observable is the junction's flux **relative to
the exon's total RNA flux** — for which the flux at the exon's *other* junction is the obvious prior-free
stand-in.

## 6. (superseded) the original "two candidates refuted" note


The magnitude needed is ~0.3–0.6. It must come from the data, not a constant.

* ⛔ **The two-junction disagreement.** An exon's two flanking junctions each claim its RNA density, so their
  gap looks like a free two-study estimate. Measured median `|log(ρ_μ,left/ρ_μ,right)|` = **0.14–0.32** ⇒ an
  implied `Var(log φ) ≈ 0.028`, against a true 0.58 — it **under-states by 20×**, because the two junctions of
  one exon share almost all of the premise error (implied correlation ≈ 0.95). This is the same
  correlated-estimator failure already on record for the left/right message gap
  (`SESSION_2026_07_25_HANDOFF_10.md` §4). Do not rebuild it.
* ⛔ **M7's `b̂²`.** It would price this in principle, but it is inert exactly where the problem is worst
  (`τ_own = 0` on all unstranded data), and §2.2 shows it would then merely absorb what a modelled term makes
  explicit.

**The leading candidate is a FITTED LIBRARY-LEVEL DISPERSION**, in the same family as the quantities the
calibrator already fits by method-of-moments from the data — `rna_sense_frac` (κ) and both strand
overdispersions. Those are *fitted parameters*, not tuned constants, and the "no magic numbers" directive has
always permitted them. What is still missing is the **observable** to fit against: the two-junction gap is
refuted, so a different statistic is needed — one whose correlation with the shared premise error is not ≈1.

Until that observable is found, P1d stays unbuilt. The probe (`RIGEL_XVAR`, default off, verified
bit-identical when unset) stays so the next attempt can size any candidate estimator against the same numbers.
