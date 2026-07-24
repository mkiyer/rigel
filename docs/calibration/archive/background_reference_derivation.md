# The aggregate DNA-background reference — derivation (2026-07-17)

**Purpose.** Derive, from first principles, how to give the calibration prior a correct notion of the *DNA
level* — the thing the total-density NPMLE is structurally blind to, and the missing piece behind both the
flagship under-call (gdna300 unstranded) and the zero-gDNA false positives (gdna_none unstranded). The released
tool already had this (a background DNA density from intron+intergenic regions, with a prior on top); the NPMLE
project superseded it and, in doing so, lost it. This derives how to put it back **elegantly**, unified with the
NPMLE, not bolted alongside it.

Companion evidence: `flagship_prior_asymmetry_diagnosis.md` (the 24-scenario A/B that rejected every *fixed*
prior); `scripts/debug/background_reference_probe.py` (the validation cited in §5).

---

## 1. The problem, stated precisely

The 24-scenario sweep showed: in an **unstranded** library the truth of `f_g` is *opposite* by DNA level —
gDNA-present wants `f_g→1`, gDNA-absent wants `f_g→0` — and the node's own composition (total density, strand)
is identical in the two cases. No **fixed** prior can be right for both; a down-bias (V0) crushes real gDNA, an
up-bias/barrier (V4) manufactures gDNA where there is none. Therefore `f_g` cannot be decided by a prior belief:
**the DNA level is a measured DATA quantity.** The question is how to measure it and fold it in.

The total-density NPMLE cannot supply it, for a structural reason: it is fit on `ρ_total = ρ_g + ρ_r`, so in a
DNA-free library it learns the RNA distribution and then mis-attributes it to DNA (a node at an RNA mode is
"rewarded" for `f_g=1`). It has **no channel that measures DNA separately**, hence no way to know DNA is absent.

---

## 2. The resolution principle — the background is a genome-wide SCALAR, not a per-region quantity

Let region `i` have effective length `E_i` and DNA count `g_i ~ Poisson(ρ_{g,i} E_i)`. The Fisher information a
single region carries about `log ρ_{g,i}` is `I_i = ρ_{g,i} E_i` — the *expected count*. When `ρ_{g,i} E_i < 1`
(a faint background over a short region), the region carries **less than one unit of information**: its count is
0 or 1, and a count of 0 is consistent with any `ρ ∈ [0, ~1/E_i]`. **A single region cannot measure DNA below
`~1/E_i`.**

A per-region model (the NPMLE, which represents `P(ρ)` allowing per-region variation) inherits this limit: at
the faint end it cannot *locate* the background — it can only say "somewhere below `1/E`". It spends the data on
the *shape* of a distribution it cannot resolve.

The escape is to change the estimand. Under a **uniformity assumption** (the DNA contamination rate is the same
genome-wide, up to local enrichment — see §7 for the tumor-CNV caveat), the background is **one number** `ρ_bg`.
Pool the pure-DNA regions into a single aggregate support:

```
    ρ_bg  =  (Σ_{i∈B} g_i) / (Σ_{i∈B} E_i) ,        Var(log ρ_bg) ≈ 1 / Σ_{i∈B} g_i .          (2.1)
```

Now the information is `ρ_bg · Σ_{i∈B} E_i` — the *sum* of the per-region informations. Even when every term
`ρ_bg E_i ≪ 1`, the sum over a genome-scale `Σ E_i` is large, so `ρ_bg` is resolved — and `ρ_bg = 0` (true zero)
is resolved *most* sharply of all (`Σg = 0` exactly). **This is the crux: the DNA background is estimable only as
a pooled scalar; it is not a per-region quantity, and no per-region prior (NPMLE or KDE) can recover it.**

`B` = the **pure-DNA regions**: intergenic + RNA-free introns. They carry DNA only (no RNA), so their pooled
density *is* `ρ_bg`. Under capture they are the depleted off-target floor; without capture they are the uniform
DNA level.

---

## 3. Factor the DNA rate into SCALE × SHAPE

Write the per-region DNA rate as

```
    ρ_{g,i}  =  ρ_bg · e_i ,                                                                    (3.1)
```

- `ρ_bg` — the **scale**: one aggregate scalar (§2), precise, the *absolute* DNA level of the library.
- `e_i` — the **shape**: the local enrichment (`e_i = 1` off-target / no-capture; `e_i > 1` on-target). It is a
  per-region quantity, resolvable *only where the DNA signal is strong* (on-target, high count). Because capture
  is nucleic-acid-agnostic (it enriches DNA and RNA alike), the *enrichment geometry* `e_i` is readable from the
  **total** density landscape — which is exactly the valid content of the current total-density NPMLE (the
  "enrichment profile", `sigma2_transfer_derivation_state`).

The NPMLE's error was to fit the **absolute** density — conflating scale and shape. The fix is to let the NPMLE
carry only the **shape** `e_i` (the enrichment, where it is resolvable) and let the **aggregate scalar** carry
the **scale** `ρ_bg` (where per-region resolution fails). This is precisely "background from intron+intergenic
+ a prior on top", read as `ρ_g = ρ_bg · e`.

---

## 4. The `f_g` prior anchor — a DATA-DRIVEN floor

For a node with observed total density `ρ_tot,i = M_i/E_i`, the prior DNA fraction is

```
    E[f_{g,i}]  =  ρ_{g,i} / ρ_tot,i  =  (ρ_bg · e_i) / ρ_tot,i ,                               (4.1)
```

with spread set by `Var(log ρ_bg)` (2.1) plus the enrichment uncertainty. The **un-enriched default** (`e_i=1`)
is the clean, load-bearing case:

```
    E[f_{g,i}] |_{default}  =  ρ_bg / ρ_tot,i .                                                 (4.2)
```

This single quantity resolves both corners that defeated every fixed prior:

- **gdna_none:** `ρ_bg = 0` (2.1, `Σg=0`) ⇒ `E[f_g] = 0`, sharply (the aggregate variance is tiny at true
  zero). The enrichment `e_i` is irrelevant — it multiplies a zero scale. **No false-positive DNA.**
- **gdna300:** `ρ_bg > 0` ⇒ a **floor** `f_g ≥ ρ_bg·e_i/ρ_tot` that the strand/messages cannot be crushed
  below. On-target (`e_i` large) the floor rises toward the truth; off-target it sits at the depleted level.

So the "barrier against `f_g→0`" that V4 tried to impose *uniformly* (and thereby broke gdna_none) is instead
**imposed by the data**: it is present exactly when and where DNA exists (`ρ_bg > 0`), and absent when DNA does
not. That is the whole difference between a fixed prior and a measured reference.

**Elegance check.** (4.2) is the released architecture — "DNA background from intron/intergenic, prior on top" —
recovered as the low-order term of (4.1), with the NPMLE now supplying the *enrichment* `e_i` (the "prior on
top") rather than the whole absolute density. Nothing new is invented; the two pieces are re-factored into
scale (aggregate) × shape (NPMLE), each used only where it is statistically valid.

---

## 5. Empirical validation (`background_reference_probe.py`, all 24)

**Discriminator (2.1).** The pooled pure-DNA density separates the DNA level cleanly and strand-independently:

| library | `ρ_bg` (capOFF / capON) | `bg f_g = ρ_bg/ρ_tot` |
|---|---|---|
| gdna_none | **0.00 / 0.00** (Σg = 0) | **0.000** |
| gdna100 | 0.19 / 0.03 | 0.06–0.55 |
| gdna300 | 0.57 / 0.09 | 0.10–0.79 |

`ss0.50 ≈ ss0.99` to three digits (gdna300 capOFF 0.567 both) — strand-independent. Under capture `ρ_bg` drops
to the depleted floor (0.57→0.09), exactly the "depleted floor" reading; the enrichment `e_i` restores the
on-target DNA above it.

**Resolution (§2).** In gdna_none, **100 %** of pure regions have count 0 — no single region sees the background
— yet the aggregate (`Σg = 0`) declares `ρ_bg = 0` with certainty. Where DNA is present, `Σg = 0.2–5.8 M` pins
`ρ_bg` to a Poisson relative error of 0.000–0.002. The total-density NPMLE's dominant low mode, by contrast,
sits at log10 ρ ≈ −3…−4 in gdna_none — a nonzero mode fit to sparse RNA, structurally unable to report "DNA
absent".

---

## 6. Incorporation into the current model

1. **Measure `ρ_bg` (+ `Var log ρ_bg`)** by pooling the count-observable **pure-DNA** regions (intergenic +
   RNA-free introns) — the same signature partition `density_model` already isolates (`region_count_observable`
   / `boundary_count_observable`). This resurrects the retired `bp_solver._floor_estimate` concept as a proper
   scalar-with-variance rather than a point floor. Under capture this is the depleted off-target floor; the
   uniformity assumption makes it a single number.
2. **Keep the NPMLE for the enrichment shape `e_i`**, but re-based on `ρ_bg`: the NPMLE models `e = ρ_g/ρ_bg`
   (dimensionless enrichment), so it is library-DNA-level-independent and its low end is *pinned* to the
   aggregate reference instead of being fit per-region where it cannot resolve.
3. **`ψ` gDNA arm** = `logP_g(ρ_g)` with `ρ_g = ρ_bg · e`, anchored per (4.1). At `ρ_bg → 0` its mass
   concentrates at `f_g → 0` (from the tight aggregate), giving the correct RNA default — replacing both the
   Jeffreys `½log f_g` barrier (which was level-blind) and the clamped left tail. The RNA arm **is** the
   reference measure `½·log(1−f_g)` (a fitted `logP_r` was tested and rejected 2026-07-19 — see
   `CALIBRATION_MASTER.md`).

This is a genuine simplification: one scalar (with variance) carries the DNA level; the NPMLE carries the
enrichment; the per-node solve reads `f_g` off `ρ_bg · e / ρ_tot`, refined by strand and messages.

---

## 7. Assumptions and where they break (for real data)

- **Uniformity** of the DNA contamination rate is what makes `ρ_bg` a single scalar. **Tumor copy-number
  changes violate it** — a locus with amplification/deletion has a locally different DNA rate. For real cancer
  samples the scalar becomes a slowly-varying field (a per-arm / per-CNV-segment `ρ_bg`), and `Var log ρ_bg`
  widens. We proceed with the single scalar now and flag that real data will need `ρ_bg` allowed to vary at the
  CNV scale (still an *aggregate* over each segment, never per-region).
- **Capture depletion.** Under capture the intergenic/intron pool is the *depleted off-target* rate, not the
  library-average DNA. That is correct for (4.2)'s floor; the on-target enrichment `e_i` (nucleic-acid-agnostic)
  lifts it. If capture is so deep that even the pooled off-target `Σg` is small, `Var log ρ_bg` reports the loss
  of precision honestly — the reference degrades gracefully to "faint, uncertain" rather than snapping to a
  wrong mode.
- **Nucleic-acid-agnostic enrichment** (`e_i` shared by DNA and RNA) is the one modelling assumption behind
  reading `e_i` from total density. It is the same assumption the σ²_transfer / enrichment-profile work already
  rests on.

---

## 8. ⚠ Robustness to strong capture — `ρ_bg` is a ONE-SIDED FLOOR, never a scale or denominator

**The concern (owner, 2026-07-17):** real libraries often carry LOW DNA (1–10 %), not the 100–300 % of our
scenarios. Under strong exon capture the off-target (intergenic/intron) DNA is depleted toward zero **even when
DNA is present** — capture is simply doing its job. So the depleted floor can be arbitrarily faint, and any
formula that **divides by it or uses it as a multiplicative scale is unstable** as capture strengthens. §3–4 as
first written are guilty: reconstructing on-target DNA via `ρ_g = ρ_bg·e` requires the enrichment
`e = ρ_g/ρ_bg` (or `ρ_tot/ρ_tot,bg`) — a division by a depleted quantity. **That is retracted.** Our validation
(§5) masked it: at 300 % DNA even the depleted floor (0.09) is count-rich; at 5 % DNA under strong capture it is
indistinguishable from the gdna_none zero.

**The inference is one-sided.** `ρ_bg` measures `ρ_DNA · c_off` (DNA × off-target capture efficiency), a LOWER
BOUND on the DNA level. `ρ_bg > 0` proves DNA is present (DNA ≥ floor); **`ρ_bg ≈ 0` does NOT prove absence** —
it cannot separate "no DNA" from "capture depleted the off-target". So the §4 mechanism "`ρ_bg=0 ⇒ f_g=0`" is
too strong: it would false-NEGATIVE the low-DNA/strong-capture case exactly as V4 false-POSITIVED gdna_none —
the same error mirrored.

**The corrected role.** `ρ_bg` enters ONLY as a **soft one-sided lower bound in log-density**, with variance:

```
    log ρ_{g,i}  ≳  log ρ_bg   (± σ_bg) ,      NEVER   f_g = (…)/ρ_bg   or   ρ_g = ρ_bg · e .
```

Three properties make this safe across all capture strengths:
1. **No division by the depleted quantity.** `ρ_bg` sits in a numerator / as a log lower-bound; the only
   denominator is the *observed* `ρ_tot` (never zero for a live node).
2. **Graceful vanishing.** As capture strengthens, `ρ_bg → 0`, `log ρ_bg → −∞`, and the bound recedes to "no
   lower bound" — correct, because off-target genuinely cannot lower-bound DNA then. It never inverts or blows
   up.
3. **One-sided ⇒ structurally safe.** It can lift `f_g` off a crush *only up to the demonstrable background*;
   it can never push `f_g` toward a target (no false positive) nor down (no crush). It cannot manufacture or
   erase DNA.

**The honest limit.** A one-sided floor does NOT solve the **enriched on-target** case — the flagship *and* the
realistic low-DNA/strong-capture sample. For an enriched node (`ρ_tot=34`, `ρ_bg=0.09`) the floor
`ρ_bg/ρ_tot ≈ 0.003` cannot lift a true-DNA exon to `f_g≈0.98`. Under capture the DNA level lives **on-target**,
in the junction-crossing unspliced signal entangled with RNA — i.e. the deconvolution itself. Therefore the
robust global DNA abundance is an **iteratively-refined OUTPUT** (the refit), anchored *below* by the one-sided
off-target floor and *lifted* by the on-target evidence — **not** a fixed scale extrapolated from the floor. The
floor is the safe anchor that cannot hurt; it is not the quantity everything is built on.

**Test-suite gap — NOW FILLED + VALIDATED (2026-07-17, `scripts/sim/configs/lowdna_capture_strength.yaml`,
`background_reference_probe.py`).** Generated low-DNA × capture-strength scenarios (gDNA ∈ {0, 1%, 5%, 100%} ×
capture ∈ {off, on, verystrong=binding60/off-target0.15}, unstranded + nascent). Result — `ρ_bg / Σg / relerr`:

| gDNA | off | on | verystrong |
|---|---|---|---|
| none | 0 / 0 / **1.0** | 0 / 0 / **1.0** | 0 / 0 / **1.0** |
| 1 % | 1.9e‑3 / 12867 / .009 | 2.9e‑4 / 1955 / .023 | **2.6e‑4 / 885 / .034** |
| 5 % | 9.5e‑3 / 64948 / .004 | 1.4e‑3 / 9791 / .010 | 9.5e‑4 / 4457 / .015 |
| 100 % | 0.19 / 1.3M / .001 | 0.028 / 202k / .002 | 0.011 / 89k / .003 |

Confirms the one-sided-floor design: (1) **true zero detected at every capture** (`Σg=0`, `relerr=1`); (2)
**graceful degradation** — `ρ_bg` shrinks smoothly but even at the worst realistic corner (1 % + verystrong)
still measures 885 DNA counts (`ρ_bg=2.6e‑4`), cleanly separated from true-zero's 0; (3) **`relerr` reports its
own reliability** (rises to .034 as evidence thins, →1 at true zero) — safe as a floor-with-variance that widens
to "no lower bound" rather than asserting a wrong value; (4) the background **fraction** collapses to .001–.002
under strong capture — empirical proof it is a floor, not a reconstructable scale; (5) the total-NPMLE low mode
is ≈ −4.5 for BOTH `1% verystrong` and `none verystrong` — the NPMLE cannot separate faint DNA from none; only
the aggregate `ρ_bg` can. Push capture/DNA further and `Σg→0`, `relerr→1`: the reference honestly stops claiming
to know and hands off to the on-target/iterative estimate — the derived graceful hand-off, observed.

## 8.5 The prior formulation — background ⊕ NPMLE without double-counting (derivation)

**The question (owner):** how do we combine the aggregate background with the NPMLE?
`(1) prior = bg + NPMLE`, or `(2) prior = bg + (NPMLE − bg)`? The released tool used (1) (background added to
the KDE) and we called it an approximation. Here is why (2) is correct, and its clean realization.

**Two estimators, two resolution ranges.** A region resolves its DNA rate only for `ρ ≳ 1/E_i` (§2). So:
- the **NPMLE** (a *variation* model, `P(ρ) = Σ w_j K(ρ; ρ_j)`) is valid on the **resolved** range `ρ ≳ 1/E`
  — the enrichment modes, where nonzero counts localize `ρ`;
- the **aggregate** `ρ_bg ± σ_bg` (a *constant* model, one pooled scalar) is valid on the **sub-resolution**
  range `ρ ≲ 1/E` — where per-region counts are 0/1 and only pooling locates the mass.

These ranges are **complementary**, exactly as intuited: the background measures the `[0, 1/E)` gap the NPMLE
cannot see.

**Why the NPMLE cannot own the background (the deeper reason).** Fit on faint pure-DNA regions (mostly count 0),
the NPMLE does not return a small uniform rate — by Kiefer–Wolfowitz **atomicity** it puts an atom at `ρ=0` plus
atoms at the few nonzero rates, because "most regions are exactly 0, a few are high" fits the zeros better than
"all regions share a small `ρ`". A *variation* model misrepresents a *uniform* faint background. The aggregate
(constant model) is the right — and efficient (all data → one number) — estimator for the uniform part.

**Why formula (1) double-counts.** The NPMLE fit on all regions already assigns some (mis-located, atomic) mass
`π_lo` to `[0, 1/E)` from the zero-count regions. The aggregate background *is those same regions*. So
`bg + NPMLE` counts the sub-resolution regions **twice** — the low end carries `≈ 2·π_lo`. That is precisely the
"approximation" in the released tool.

**Formula (2), made exact.** Partition the regions by resolution — *disjoint sets*:
```
    P(ρ)  =  π_bg · Bg(ρ; ρ_bg, σ_bg)     +     (1 − π_bg) · Enrich(ρ) ,                       (8.5)
```
- `Bg` = the aggregate background component (constant model), owning the sub-resolution mass at its correct
  location `ρ_bg` with spread `σ_bg`;
- `Enrich` = the NPMLE fit on the **resolved** regions only — the enrichment shape;
- `π_bg` = the fraction of DNA mass that is sub-resolution (the depleted-floor fraction).

`(NPMLE − bg)` in the owner's formula (2) is exactly "the NPMLE with its unreliable sub-resolution mass removed"
= `(1−π_bg)·Enrich`; adding `π_bg·Bg` gives (8.5). Disjoint supports ⇒ **no double-count.**

**The clean realization — a PINNED background component (no threshold, no literal subtraction).** Fit the NPMLE
mixture with **one extra component whose location is pinned at the aggregate** `ρ_bg` (kernel width `σ_bg`):
```
    P(ρ)  =  w_bg · N(log ρ; log ρ_bg, σ_bg²)   +   Σ_j w_j K(ρ; ρ_j) ,       {ρ_j : free enrichment kernels}
```
EM fits the weights `{w_bg, w_j}`; only the background component's **location is fixed** (by the aggregate,
which is precise where the per-region data is not). The sub-resolution mass then flows to the pinned component —
placed correctly by the aggregate rather than misplaced by atomicity — and the free kernels capture the resolved
enrichment. This is (8.5) realized at the fitting stage: elegant, one extra kernel, no explicit `1/E` threshold.

**How the model changes, precisely.**
- **Before:** `P(ρ_g) = DensityNPMLE.fit(total density)` — a single NPMLE mixture on total density, which both
  conflates DNA/RNA and (by atomicity) misrepresents the faint background.
- **After:** (i) measure `ρ_bg ± σ_bg` from the pooled pure-DNA (intergenic + RNA-free intron) regions; (ii)
  fit the mixture with the background component **pinned** at `ρ_bg`; (iii) the prior enters the `f_g` solve as a
  **one-sided log-floor** (§8), so the pinned background anchors the low end and recedes to `−∞` as `ρ_bg → 0`
  — never a denominator.

This keeps the released architecture's spirit ("background from intron/intergenic + a prior on top") but makes
the combination **exact** (pinned component = no double-count) and **resolution-honest** (constant model for the
uniform background, variation model for the enrichment).

## 9. Summary (revised)

The DNA background is a genome-wide aggregate `ρ_bg` from the pure-DNA (intergenic + intron) regions — resolvable
even at true zero because pooling turns per-region non-resolution into one precise scalar (§2). But it measures
`ρ_DNA · c_off`, a LOWER BOUND that capture drives toward zero, so it enters the prior ONLY as a **soft
one-sided log-floor with variance — never a scale, factor, or denominator** (§8). `ρ_bg > 0` demonstrates and
floors DNA; `ρ_bg ≈ 0` is inconclusive (strong capture ≡ no DNA off-target). The enriched on-target DNA level —
the flagship and the realistic low-DNA/strong-capture case — is carried by the on-target deconvolution and an
iteratively-refined global abundance, with `ρ_bg` as the robust lower anchor. This keeps the release tool's
"background from intron/intergenic + prior on top" while making it stable as capture → ∞.
