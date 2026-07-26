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
| **P1d** | `ω_graft` = `max(0, E[d²]−E[noise])/2` (`graft_premise_logvar`) — **ONE library-level scalar**, MoM. ⚠ **A DEBT, not a model — §3.5** | the graft's **premise**: a seam's RNA used as the exon's RNA, when it is only a LOWER BOUND on it | the WHOLE grafted RNA claim, both code paths | **LIVE** (2026-07-26) — §3, and ⚠ **its mechanism is NOT what §3 originally said**, see §3.0 |
| **P1e** | `b̂²_cons = max(0, δ²−αᵀΣα−ν)` on every LEVEL (`conservation surprise`) | the claim asserting more fragments than the node observed | every message, **scoped to `δ < 0`**; λ and θ untouched | **LIVE** (2026-07-26) — ⚠ **partly a DEBT, §6** |
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
  this many"; charging it a premise variance says "and what I counted does not necessarily speak for the
  region I am applying it to". A count cannot see the second and never will, however large it is.
* **M8 vs P1d** — both sit on the graft and they are disjoint in the way §2.1 requires: M8 prices the CAPTURE
  step (`(log r)²`, **identically 0 off-capture**), P1d prices the structural premise (measured
  capture-INVARIANT: residual 0.41–0.46 across gdna1…gdna300 and both strandednesses). M8 is the only term
  that could plausibly have absorbed P1d, and off-capture it charges exactly zero while the premise error is
  at its full 0.485 — which is how the gap was found.
* **M10 vs P1d** — the exact mirror pair, on opposite edges. M10 prices the share on the peel (exon →
  boundary); P1d prices the graft's (boundary → exon), which had **no** share term at all. Neither fires on
  the other's edge.

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

## 3.0 ⭐ P1d IS BUILT — and the mechanism below it is NOT extrapolation. Read this before §3–§5.

**Landed 2026-07-26** as `enrichment_frame.graft_premise_logvar`, ON by default, `RIGEL_GLV_OFF=1` ablates it
(verified bit-identical to the pre-P1d path, 32/32). §3–§5 below are kept as the evidence trail, but **three
of their load-bearing claims are now refuted by measurement**, and one of the refutations reverses the name of
the term.

**① "Extrapolation" is the wrong mechanism.** §3 attributes the error to a ~100 bp junction speaking for a
~2,100 bp exon, a 12–21× extrapolation. Binned by that very ratio the premise residual is **FLAT — 1.13× over
a 6.7× range** — and flat in exon length too. It is not a geometry effect.

**② The 40× count gradient of §5.2 is a PROXY, not a law.** Stratify by one structural bit — does the boundary
carry a **transcript terminus** — and the count gradient vanishes *inside* each stratum (junction-only excess
0.036 → 0.073 across a 200× count range), while the strata themselves differ by **≥30×**:

| stratum | edges | median φ | **ω̂** | share of squared deviation |
|---|---|---|---|---|
| boundary carries a **TERMINUS** | 20.8 % | 0.70 | **1.7–1.9** | **71.7 %** |
| **junction-only** boundary | 79.2 % | 1.10 | **0.04–0.06** | 28.3 % |

Verified three ways including variance components across independent simulated libraries (`Var_within` = the
counting noise = 0.1414; floor+step = 0.1406, a 0.6 % match). That check is also what caught the first pass
leaving the ORACLE capture step's own two gDNA counts inside the junction-only figure — 0.086 → **0.04**.
**This is `ROADMAP.md` §6's deferred structural item ("the region/boundary map has no TSS/TES", deferred as
"expected to be low-impact"), and it is 62–72 % of the graft's premise error.**

**③ The Poisson floor was under-subtracted, but it is not the story.** §3's `sd(log φ) = 0.76 ⇒ 0.58` is a
statistic built from FIVE counts with only one subtracted. With all five subtracted at exact trigamma, Poisson
explains **16–24 %** and the residual is **0.485** (MC null returns −0.017, so the law-of-total-variance
subtraction is sound). ⚠ §5.2's own "excess" column additionally used the spliced **mass**, which is ~2× the
count — the error `bp_solver.py`'s own comment warns about.

**The physical statement that replaced "extrapolation".** Every molecule counted at a seam is in the exon; the
exon may hold molecules that never touch that seam. So the graft knows an **inequality**,
`ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)`, and asserts it as an equality (`φ ≡ 1` to 2.2e-16 in both code paths). The
premise fails hardest exactly where RNA does not flow through the seam — a transcript terminus.

**What landed, and what did not.**

| variant | verdict |
|---|---|
| the **VARIANCE** as ONE library-level MoM scalar `ω_graft`, applied to every graft edge, with M5's lift noise subtracted | ✅ **LANDED** |
| the same term **per-edge**, using each exon's own two-seam gap | ⛔ **REVERTED same day** (owner's objection, and the measurement agrees). `d²` from ONE pair is a χ²₁ draw — CV = √2 — so a per-edge variance is mostly noise, and it *replaces* the population value on the ~48 % of edges where a pair exists, UNDER-charging as often as over. Pooled-only is better on every axis: confidently-wrong **870,245 → 762,000**, `z2` ALL **11.12 → 8.98**, exon-single **5.68 → 2.46**, and mwae **+0.0016 → +0.0010**. It also removed a real BP violation: the LEFT seam's message was carrying a variance built from the RIGHT seam's counts |
| the **MODE** fix (replace the claim with the dominant seam's flux, `max`) | ⛔ **not landed.** Right physics, wrong regime: capture-OFF **14/14 better, 0 worse**; capture-ON 8 worse. Under capture the graft ALREADY over-states (median φ = 2.45, an M8 frame problem), so no bound-tightening can help. Framing the two seams before comparing halves the damage but does not remove it. `min` is worse than either bound alone (0.514 vs 0.487) and `sum` over-states by +0.61 nats — both exactly as the algebra predicts, which is what makes the `max` result structural rather than a selection artifact |
| omitting M5's lift noise from the MoM subtraction | ⛔ inflates the fit in proportion to gDNA depth and over-charges precisely where the RNA claim is near-exact and the gDNA claim is the wrong one (P1b's population). Costs +0.0010 mwae to fix, and it is the honest form |

**Measured effect (32 conditions, `RIGEL_GLV_OFF=1` vs default):**

| | mwae r0 / r1 | **confidently-wrong** | z2 exon-single | z2 exon-AMBIG | **z2 ALL** |
|---|---|---|---|---|---|
| pre-P1d (`pk2`) | 0.0855 / 0.0671 | 1,200,951 | 8.60 | 92.10 | 15.55 |
| per-edge + pooled (reverted) | 0.0871 / 0.0680 | 870,245 (−27.5 %) | 5.68 | 65.28 | 11.12 |
| **P1d, pooled scalar only** | **0.0865 / 0.0676** | **762,000 (−36.6 %)** | **2.46 (−71 %)** | **64.39 (−30 %)** | **8.98 (−42 %)** |

9 better / 6 worse / 17 flat at refit=0; 8 / 4 / 20 at refit=1. **The trade is +1.2 % error mass for −36.6 %
confidently-wrong mass** — the trade `PASS0_FINISH_PLAN.md` §0b says decides pass-0. The un-landable flat
probe `RIGEL_XVAR=0.3` bought −42 % for +2.4 % and reaches only one of the two code paths; P1d gets most of
that for half the mass, prior-free and with no constant.

⚠ **The residual regression is localized and has a known mechanism**: unstranded × gDNA-rich × capture-OFF. That is P1b's population, where the
RNA claim is a near-exact measurement (36.6453 vs 36.6734) and the **gDNA** claim is 47× too big — so damping
the RNA arm is the wrong arm. Charging the premium to the arm whose premise actually failed is the follow-up.

⚠ The suite is **Poisson by construction**, so everything above is structural/annotation-driven dispersion
with no overdispersion component, and the terminus/junction split must be re-measured on a real annotation.
The suite's isoform structure is nested 5'/3' truncations plus exon skips of a shared chain — no alternative
TSS/TES *inside* an exon, no mutually exclusive exons, no retained introns.

## 3. P1d — the one gap the ledger exposes, and it is large (⚠ superseded in part by §3.0)

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

## 3.5 The fitted scalar `ω_graft` — what it is, and the four questions answered about it

> ### ⚠⚠ P1e IS PARTLY A DEBT — IT PRICES A BIAS AS A VARIANCE. DO NOT LET THIS BE FORGOTTEN
>
> P1e damps a message by the unexplained part of its conservation violation `δ = log(M/S)`. **On a large
> share of its firing mass `δ` is a systematic BIAS, not scatter** — measured `E[δ]` ≈ −0.5 to −1.5, with a
> bias share of **53–77 %** on graft × one-component messages and **98.9–99.2 %** at intergenic
> destinations. A variance prices random scatter; pricing a bias as a variance does not move the mode toward
> truth, it only weakens the message — **and it never shrinks**, which breaks the ledger §2.2 guarantee that
> a DerSimonian–Laird residual vanishes as the model improves.
>
> **It is landed anyway** because it is the only change measured to improve ACCURACY and honest PRECISION at
> the same time (suite 0.0870 → 0.0841, fit-substrate 0.0883 → **0.0845**, held-fixed `z2` ALL 8.53 → **1.74**,
> exon-single 3.57 → **0.88**), and because pass-0 is explicitly allowed to be weak-and-correctable. That is
> a pragmatic trade, not a derivation.
>
> **Four things that must not be overlooked:**
> 1. ⛔ **It was hoped the bias-dominated strata were inert** (intergenic nodes are `solvable = False`).
>    **Measured and REFUTED: 90–100 % of the damping mass lands on solvable destinations.** The bias is live.
> 2. **The magnitude is not what works.** Control arms: a flat pooled constant on the same firing set beats
>    the derived `b̂²` on 3 of 4 conditions, and `b̂² := δ²` with no null subtraction is identical. Permuting
>    `b̂²` within edge class FAILS and damping all grafts uniformly FAILS — so **`δ` identifies WHICH message
>    to distrust, and the calibrated magnitude adds nothing.** Exactly ω_graft's shape.
> 3. **`δ` is not the hard observable it was scoped as.** It is **M-free at every region node** (verified to
>    1.9e-16 — the message is built ∝ `M`, so `M` cancels). On a plain edge it is a composition disagreement
>    with the destination's own belief. The genuine content is in the ROUTING: **84.2 % of Σδ² is on peel and
>    graft edges.**
> 4. **The right fix for the bias half is a MODE fix, not a variance.** When the bias strata are diagnosed,
>    P1e must SHRINK — if it does not, the model has not improved.
>
> ### ⚠⚠ `ω_graft` IS A DEBT, NOT A MODEL — DO NOT LET THIS BE FORGOTTEN
>
> P1d's fitted `ω_graft` **partially compensates for a FAILURE IN OUR STRUCTURAL REPRESENTATION.** The
> region/boundary map has no TSS/TES, so the solver **cannot tell a splice junction from a transcript
> terminus** — and that distinction is the whole of the effect `ω` prices: measured `ω̂` = **1.7–1.9 at
> terminus boundaries vs 0.04–0.06 at junction-only ones, a ≥30× split**, with 20.8 % of edges carrying
> 71.7 % of the error. `ω` is a single library-wide average standing in for a bimodal quantity. It works
> because over-charging a variance is cheap and under-charging is expensive — **not because it is right.**
>
> **Three consequences that must not be overlooked:**
> 1. **It is expected to be fragile on real data.** It is fitted on ~200 exons with a 30–50 % standard
>    error, on a suite that is Poisson by construction and whose isoforms are only nested truncations plus
>    exon skips. Real annotations have alternative TSS/TES *inside* exons, mutually exclusive exons and
>    retained introns — none of which this suite contains.
> 2. **`ω` MUST be re-derived as a per-class quantity the moment TSS/TES enter the region map.** Same
>    estimating equation, two (or more) scalars keyed on the structural bit, each fitted over a
>    homogeneous population. The partial-pooling code already in `graft_premise_logvar` is the plug-in point.
> 3. **Until then, treat every `ω`-dependent result as provisional**, and never quote the fitted value as
>    if it were a measured constant of the library.


`ω = max(0, E[d²] − E[noise])/2` over exons with spliced flux on **both** flanking boundaries.
`d = log(ρ_μ(left)/ρ_μ(right))`, both lifted into the exon's frame; `noise` = the two seams' spliced counts
⊕ the two lifts' `logvar_tot`. **One scalar per strand per library.** Measured (relay fit, + strand):

| condition | pairs | E[d²] | E[noise] | **ω** | SE(ω) | sd(d) | med\|d\| |
|---|---|---|---|---|---|---|---|
| gdna300 ss0.99 capOFF | 192 | 0.692 | 0.099 | **0.297** | 0.111 | 0.827 | 0.102 |
| gdna300 ss0.50 capOFF | 194 | 0.676 | 0.209 | **0.233** | 0.108 | 0.819 | 0.115 |
| gdna300 ss0.99 capON | 153 | 1.155 | 0.312 | **0.421** | 0.159 | 1.074 | 0.331 |
| gdna100 verystrong | 105 | 0.995 | 0.177 | **0.409** | 0.172 | 0.997 | 0.230 |
| none ss0.50 capOFF | 197 | 0.620 | 0.220 | **0.200** | 0.093 | 0.786 | 0.046 |

Fit population: **~190–200 of 1,698 regions (≈11 %), 100 % exons** — an intron never qualifies, because its
flanking boundaries carry their spliced flux on the face pointing at the *exon*. Against the variance the
graft was charged before (0.008–0.09), `ω` dominates by **5–36×**: it says a seam's RNA claim about its exon
is good to about **±0.55 nats — a factor of 1.7 — before a single read is counted**.

**Q1 — is it belief-dependent?** Weakly, yes. The two fluxes are lifted by frames read off `ρ_tot(f_cur)`, so
it is refit 3× per sweep: 0.2968 / 0.2968 / **0.2937**, <2 % movement. The count and `logvar_tot` legs are
fixed data. The feedback can only WIDEN, so it cannot manufacture confidence.

**Q2 — should the average be confidence-weighted (DerSimonian–Laird)?** No, and for a reason, not a
preference: DL is the efficient estimator of a **common** effect, and there is no common effect
(`τ²_between` = 1.35–3.0 vs sampling 0.26–0.78). The quantity to charge a message is the expected squared
premise error of a randomly-drawn edge — `E[ω_i]`, the arithmetic mean — which the plain moment difference
estimates and DL does not. Measured they agree to ~1 % anyway (762,000 vs 755,727 confidently-wrong, DL
slightly worse on mass), so the argument is correctness, not score.

**Q3 — can it be per-node instead of library-wide?** ⛔ **Tested twice, lost twice.**

| variant | confidently-wrong | z2 ALL | z2 exon-single |
|---|---|---|---|
| raw per-edge (own `d²`) | 870,245 | 11.12 | 5.68 |
| **derived PARTIAL POOLING** (`B = τ²_b/(τ²_b+Var(ω̂_i))`, `Var(ω̂_i) = (2ω+v_i)²/2` from the χ²₁) | 765,281 | 9.23 | 2.68 |
| **pooled scalar** | **762,000** | **8.98** | **2.46** |

And the shrinkage weight comes out **B = 0.82–0.89**, i.e. the estimator itself says the between-node
variation is real and dominant — the heterogeneity is NOT a fitting artifact. It still loses, and the reason
is the **loss function's asymmetry**: over-charging only widens a message (cheap), under-charging leaves a
node confidently wrong (poisons the hyperprior). `ω_i` is heavily right-skewed — the top decile of pairs
carries **78–92 %** of `Σd²` — so any per-node point estimate, shrunk or not, charges the median node ≈0.
**The population mean is the right single number precisely because the loss is asymmetric.**

**Q4 — so what would actually fix it?** Not a better *estimate* of `ω_i` — a **classifier** for which nodes
are in the high-`ω` class. That is the terminus bit (P1g): `ω̂` 1.7–1.9 vs 0.04–0.06, a ≥30× split the
current map cannot see. `graft_premise_logvar` then becomes two scalars keyed on a structural bit, fitted by
the same equation, and the partial-pooling machinery already in the law is where it plugs in.

⚠ **Do not over-read the fitted value.** `E[d²]` is tail-dominated over ~200 pairs, so SE(ω) is **30–50 %
relative**, and the value is the mean of a bimodal population. That P1d works is robust to it: the flat probe
`RIGEL_XVAR=0.3` and DL's 0.11–0.27 all buy most of the same gain. **The fitted-ness matters for portability
to real data, not for the score here.**

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


## 6. ⚠ P1e — what is derived, what is landed, and what is owed

**Derived and MC-validated:** the COMMON direction. `αᵀ1 ≡ 1` (α is a share vector over one budget), so
adding the scalar `b̂²` to every supplied component's log-variance satisfies `αᵀΣ'α = αᵀΣα + b̂²` exactly and
leaves **`Var(λ)` identically unchanged** — a common shift of both arms cannot move the split. The rank-1
alternative (inflate along `Σα`) borrows the *conditional mean's* direction and applies it as a variance;
MC shows it over-damps λ **5×** on a pure scale error (λ z² 1.00 → 0.21) while still leaving the RNA arm
over-confident (2.88). **A scalar observation cannot identify a direction.**

**Derived:** the scope. `S` is a COMPLETE budget — `_pin_v` fills every unsupplied component from the node's
OWN density — so a shortfall (`δ > 0`) can be the node's own density being too low just as easily as the
message being wrong; it does not attribute. The over-claim direction does: all densities are non-negative,
so nothing the unsupplied components could be would rescue a budget already exceeding `M`. Measured, the
`δ > 0` half is **inert** (scoped 0.0841 vs unscoped 0.0841 to 4 dp), so scoping costs nothing and removes
the unlicensed half. A strictly harder test (`sup` — the SUPPLIED components alone over-claim) is available
behind `RIGEL_P1E_SCOPE=sup` and costs +0.0004.

**Landed (all 32 conditions, `RIGEL_P1E_OFF=1` ablates, verified bit-identical 32/32):**

| | HEAD | **P1e** |
|---|---|---|
| suite mwae r0 / r1 | 0.0870 / 0.0697 | **0.0841 / 0.0679** |
| **fit-substrate mwae** | 0.0883 | **0.0845** |
| suite ERR reads | 12,167,916 | **11,729,507** |
| capture OFF | 0.0457 | **0.0349** (4 better / **0 worse**) |
| unstranded × capON | 0.1572 | 0.1657 (1 / 4) ← the cost |
| held-fixed z2 ALL | 8.53 | **1.74** |
| z2 exon-single | 3.57 | **0.88** |
| z2 exon-AMBIG | 59.4 | **19.8** |
| z2 intron (the measurement class) | 1.55 | **1.54** — untouched |

**Owed — see the banner at the top of this file.** The bias half, the uncalibrated magnitude, and the fact
that `δ` is M-free at region nodes (so it is a composition disagreement with the destination's own belief on
a plain edge, and only the ROUTING — 84.2 % of Σδ² — carries hard-observable content).
