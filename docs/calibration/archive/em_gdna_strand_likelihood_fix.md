# The per-fragment gDNA strand-likelihood leak — diagnosis & fix design

**Status:** design proposal, not yet implemented. **Date:** 2026-06-14.
**Context:** the last-standing cause of the flagship gDNA→RNA leak after the displaced-mass
prior fix (`main@17d0209`) and after the effective-length was *exonerated* (the "eff-len too
long" finding was a smear-counting-oracle artifact — see
`memory/flagship_leak_em_not_calibration.md`). The `×0.5` eff-len knob that "fixed" the leak is a
**band-aid** for the bias described here.

---

## 1. TL;DR

The per-locus EM scorer gives the **gDNA component a flat strand likelihood** (`LOG_HALF` =
`log 0.5` for every fragment, `scoring.cpp:516`), while the RNA components (mature isoforms **and**
nascent spans) use the **stranded** model (`scoring.cpp:553–563`). gDNA is genuinely unstranded, so
`0.5` is the *correct per-fragment* likelihood — but it is a structural **handicap**: for the
protocol-expected (common) read orientation, a stranded RNA component scores `≈ log(0.99) = −0.01`
versus gDNA's `−0.69`, a **+0.68 log-odds advantage to RNA on every such read**. Because gDNA is
50/50, **half of its reads carry that orientation and leak to the stranded RNA twins** — preferentially
the **nascent** span (unspliced, gene-body-spanning, a geometric twin of gDNA), which in a
no-nascent library is a pure sink.

The fix cannot be a different per-fragment gDNA strand term — gDNA really is unstranded. The only
thing that distinguishes a 50/50 gDNA pile from a κ-stranded nascent pile is the **population strand
balance**, which **calibration already deconvolves** into `f_g`. The EM re-litigates it per-fragment
and overrides calibration. The fix is to stop that override.

---

## 2. The exact mechanism

### 2.1 The two strand terms (code)

- **gDNA hypothesis** (`scoring.cpp:512–520`), unspliced fragments only:
  `hit_log_ll = gdna_fl + gdna_log_sp + LOG_HALF + log_nm` → the strand term is **`log 0.5 = −0.693`,
  unconditionally**, for every fragment regardless of its alignment orientation.
- **RNA hypothesis** (mature + nascent transcript rows, `scoring.cpp:553–563`):
  `log_strand = same ? log_p_sense_ : log_p_antisense_`, where `same = (read orientation == transcript
  strand)`. At the measured `κ = rna_sense_frac ≈ 0.01` (reverse-stranded library):
  `log_p_sense_ = log(κ) ≈ −4.61` (read *matches* transcript strand — rare),
  `log_p_antisense_ = log(1−κ) ≈ −0.01` (read *opposite* — the common, expected orientation).

### 2.2 The per-fragment likelihood-ratio table (ss = 0.99)

| read orientation (rel. transcript) | gDNA | nascent/mature | winner |
|---|---|---|---|
| antisense-to-transcript (**common**, ~99% of true RNA) | −0.69 | **−0.01** | **RNA by +0.68** |
| sense-to-transcript (rare) | −0.69 | −4.61 | gDNA by +3.92 |

gDNA fragments are 50/50, so **~half are antisense-to-transcript and are out-competed by RNA by
0.68 log-units each.**

### 2.3 Why the EM actively prefers the leak (fixed-point argument)

Consider a locus whose unspliced gene-body reads are *truly all gDNA* (50/50 orientation, N reads).
Compare the strand likelihood of two assignments:

- **all-gDNA:** `0.5^N`.
- **split** (nascent claims the N/2 antisense reads, gDNA keeps the N/2 sense reads):
  `0.99^(N/2) · 0.5^(N/2)`.

Ratio split/all-gDNA = `(0.99/0.5)^(N/2) = 1.98^(N/2)` — astronomically `> 1`. **The EM's
strand likelihood is maximized by splitting unstranded gDNA into gDNA(sense) + nascent(antisense).**
This is the nascent component **overfitting the antisense half of unstranded gDNA**, and it is a model
**misspecification**: a κ-stranded component legitimately *expects* an all-antisense pile, which is
exactly what the antisense half of gDNA looks like. **Per fragment, the two are indistinguishable.**

Only the prior (`θ_gdna ≫ θ_nascent`) opposes this, and it is a soft pseudocount that the accumulated
per-fragment evidence overruns — hence the leak survives a roughly-correct prior, and cranking the
prior or shrinking the eff-len (both raise `θ_gdna`) mask it.

---

## 3. Why the obvious fixes do **not** work

1. **A different per-fragment gDNA strand term.** gDNA is genuinely unstranded; `0.5` is the correct
   per-fragment conditional. Any other value is wrong and would misbehave on real data.
2. **A per-*component* Binomial strand term.** `Binomial(s; n, μ) = C(n,s)·μ^s(1−μ)^(n−s)` equals the
   per-fragment product up to the assignment-independent `C(n,s)` — **no discrimination gain.**
3. The discriminator (50/50 vs κ) is a property of the **joint** read-orientation distribution. It is
   captured only by (a) **overdispersion coupling** (a Beta-Binomial ties the reads via a shared latent
   `p`) or (b) a **binding population constraint**. Calibration already computes the optimal population
   answer (`f_g`, via the strand Beta-Binomial deconvolution). **The strand is effectively used twice**
   — once correctly in calibration, once destructively per-fragment in the EM — and the EM use wins.

---

## 4. Fix options

### Option A — per-component Beta-Binomial strand likelihood *(gold standard, invasive)*

Score each component's assigned-read strand balance `(s_c, n_c)` against its expected mean via
`BetaBinomial(s_c; n_c, μ_c, od_c)`, using the calibration's **already-fitted** dispersions:
`μ_gdna = 0.5, od_g = gdna_strand_overdispersion`; `μ_rna = κ, od_r = rna_strand_overdispersion`.
The BB couples a component's reads through the shared latent `p`, so an **imbalanced** gDNA assignment
(all-sense, or all-antisense) is *low* likelihood under `BB(0.5, od_g)` → the EM is forced to keep gDNA
**balanced** → it cannot shed its antisense half to nascent. Correctly degrades: at ss = 0.5 every
`μ → 0.5`, the BB is strand-uninformative, and the count/FL/prior govern; with *real* nascent the pile
genuinely is κ-stranded and fits nascent.

- **Pros:** statistically correct; robust to prior/`f_g` error; one mechanism for all conditions.
- **Cons:** the E-step is no longer a simple per-fragment softmax — the BB couples fragments within a
  component, requiring per-component strand sufficient statistics and a mean-field / gradient update
  inside `em_solver.cpp`, plus rework of the pool-separated multimap pruning. Largest change; touches
  the hot loop and the EM's mathematical form.

### Option B — bind the calibration `f_g` as the gDNA-vs-(unspliced-RNA) split *(M-step regularizer, recommended)*

Calibration's `f_g` is the population-strand-optimal gDNA fraction of the unspliced mass. Make the
locus's **aggregate** `θ_gdna : θ_unspliced-RNA` ratio a **binding** constraint at `f_g : (1−f_g)`
(a hard constraint, or a Dirichlet pseudocount scaled to the locus depth so it cannot be overrun),
instead of the current soft pseudocount. The per-fragment strand may still reallocate *which* reads land
where, but it cannot change the *count* split — so the **abundance / net-flow leak is removed** even
though individual hard labels remain (acceptably) unrecoverable.

- **Pros:** principled ("use the strand once, in calibration; don't re-litigate per-fragment");
  surgical (the grouped prior update already computes `a_g/a_r`; this strengthens/binds it at the
  *correct* value rather than the band-aid `×5`); directly targets the abundance metric we care about.
- **Cons:** trusts calibration's `f_g`. Validated at high SS (the spatial check confirms the per-region
  gDNA is right); at ss = 0.5 the strand is uninformative and `f_g` falls back to the count module, so
  binding is only as good as the count module there — a *separate* calibration-quality question, and a
  regression to watch.

### Option C — gate the nascent-component prior on nascent evidence *(cheap complement)*

Tie each nascent span's prior to the **local unspliced strand imbalance** (the nascent strand signal):
a strand-*balanced* (gDNA) gene body yields nascent prior `≈ 0`, removing the specific sink. Complements
B. Must **not** fire when nascent is real (nrna_rnd): there the unspliced strand is genuinely imbalanced,
so the nascent prior is correctly `> 0`.

---

## 5. Recommendation

**Primary: Option B.** It is sufficient for the abundance/net-flow goal (it fixes the *count* split,
which is what the leak metric measures), it is the principled version of the eff-len/prior band-aid
(enforce the *correct* `f_g` rather than an arbitrary multiplier), and it is the smallest change. Pair
with **C** if a residual nascent sink remains.

**Fallback: Option A**, if B proves insufficient — e.g., if the within-RNA gDNA-vs-nascent split needs
finer per-fragment power than the aggregate constraint provides, or if we want the EM itself to be
robust to `f_g` errors (especially at ss = 0.5). A is the long-term "right" model but a much larger EM
rework; defer until B is measured.

**Decisive prediction / built-in check:** after B (or A), **re-measure the eff-len demand** with
`scripts/real_validation.py`. If the leak was the strand bias, the eff-len `×1` (the real IPR) should now
recover gDNA on its own — confirming the eff-len was never the problem and we can drop the band-aid for
good. If a genuine eff-len/IPR magnitude error remains, it will now show *cleanly* against an unbiased
likelihood.

---

## 6. Validation plan & regression guards

1. **Flagship** (`gdna300/ss0.99/capture-on`): gDNA → ~3.0M, nRNA → ~0, mRNA → ~1.0M (real pipeline,
   `summary.json` totals + net-flow per-locus).
2. **Full 20-condition net-flow benchmark** (`calibration-benchmark` skill). Watch especially:
   - **ss = 0.50** (strand uninformative): the fix must **degrade gracefully** — no gDNA *over*-call;
     B must not bind a bad `f_g`, A's BB must go uninformative.
   - **gdna_none / low-gDNA**: the fix must **not suppress real RNA** (no spurious gDNA from the
     constraint/BB).
   - **nrna_rnd (real nascent)**: the fix must **not kill genuine nascent** — B respects `f_g` (which
     detects real nascent via the strand imbalance); C must not zero a real nascent prior.
3. **Unit tests + goldens** regenerated; conservation (`gdna + rna = total`) preserved.
4. The per-fragment hard-label vs soft-mass divergence (workflow caveat) should shrink once the
   abundance is bound.

---

## 7. Relationship to the prior decision

This reopens the decision recorded in `memory/em_strand_per_fragment_penalty.md`, which chose *not* to
build a per-fragment EM strand term on the belief the leak was a hard-label artifact. The **net-flow
evidence contradicts that belief** — it is a real *abundance* bias, and it is entangled with the
eff-len (the band-aid). Note Option B is *not* the per-fragment penalty that was previously rejected; it
is a population-level M-step constraint. Option A is a per-*component* (not per-fragment) term. Both
sidestep the original objection (a per-fragment penalty cannot discriminate gDNA from nascent anyway —
§2.3).

---

## 8. Working the math (2026-06-14) — Option A is **refuted**; the leak is a competition/optimization problem, not a missing strand term

Before committing a strand term to the (pristine) EM, work the likelihood algebra. It overturns the
plan.

### 8.1 The decisive likelihood calculation

Take a locus whose unspliced gene-body reads are **truly all gDNA**: `n` reads, balanced (`n/2` sense,
`n/2` antisense). Components: gDNA (`q = 0.5`) and a nascent twin (`P(antisense) = 1−κ ≈ 0.99`). Compare
the **complete-data** log-likelihood (abundance term `log θ_c` + per-read strand term) of two
assignments:

- **all-gDNA** (`θ_g = 1, θ_n = 0`): every read scores `log θ_g + log P(orient|gDNA) = 0 + log 0.5`.
  Total `= n·log 0.5 = −0.69 n`.
- **split** (`θ_g = θ_n = 0.5`; gDNA keeps the sense half, nascent takes the antisense half):
  sense→gDNA `= log 0.5 + log 0.5 = −1.39` (×`n/2`); antisense→nascent `= log 0.5 + log 0.99 = −0.70`
  (×`n/2`). Total `= −1.04 n`.

**all-gDNA (`−0.69 n`) strictly beats the split (`−1.04 n`).** The per-fragment strand term *plus the
abundance term* already make "no leak" the optimum. **The strand likelihood is not the bug.** (The
`θ_g = 0.5` penalty the split pays on the *sense* half is what dominates — splitting wastes half of
gDNA's abundance.)

### 8.2 Why no per-component strand term (incl. the Beta-Binomial) can help — Option A refuted

The nascent component's strand model is `Bernoulli(q = κ)`. Its likelihood of an **all-antisense**
pile is `(1−κ)^n = 0.99^n` — the **maximum possible** for nascent (every read in its expected
orientation). So nascent is **rewarded, not penalized**, for cherry-picking the antisense half of gDNA.
*No* per-component strand model changes this:

- **fixed-`q`** → the above (rewards all-antisense).
- **Beta-Binomial** → for **distinguishable** reads (we observe *which* read is which orientation, not
  just the count), the marginal is a **Pólya urn**: `B(a+α, n−a+β)/B(α,β)`, which for any mean is
  **maximized at the imbalanced extremes** (`a = 0` or `a = n`) — it *rewards* imbalance, the **opposite**
  of "enforce balance." The balance-peaked shape only appears with the `C(n,a)` combinatorial factor,
  which belongs to the *indistinguishable-count* model — not ours.

So "each component must approximate its `ss`" does **not** constrain a κ-component: an all-antisense pile
*is* `ss ≈ 0.99`. The discriminator between "the antisense half of gDNA" and "a nascent component" lives
**only** in the mixture/population (the `+`-fraction identifies `f_g`), which calibration already
computes and the per-fragment EM already embodies. **Option A is geometrically/statistically void.**

### 8.3 The real mechanism — the gDNA-vs-nascent **abundance** competition

If all-gDNA is the optimum, why does the EM leak? Because the **split is also a fixed point**, and the
EM lands in (or near) its basin:

- The split holds iff `θ_n > (0.5/0.99)·θ_g ≈ 0.505·θ_g` (then antisense reads prefer nascent). With
  `θ_n = 0` the all-gDNA fixed point is **absorbing** (no read can seed nascent). So the leak requires
  something to **seed `θ_n > 0`** — the **prior's nascent allocation** (the unspliced-RNA pseudocount is
  split across RNA components incl. nascent) and the **warm start**.
- Once seeded, `θ_c = count_c / eff_len_c`: if the **nascent eff-len** (its span) is small **relative to
  the gDNA eff-len** (the IPR), nascent's per-base abundance is inflated and the split basin deepens.
  This is why **`eff_len×0.5` and `prior×5` both "fixed" the leak** — both raise `θ_g`, pushing the EM
  out of the split basin. They are band-aids on the *competition*, confirming the lever is **relative
  abundance (eff-len + nascent seeding)**, not the strand likelihood.
- It also explains the **monotonic** response to eff-len (a smooth competition shift) and the
  **hard-label vs soft-mass divergence** (the soft EM sits part-way into the split basin).

### 8.4 Corrected direction — and what must be verified before any code

The fix is **not** a strand term. It is to make the EM **reach the all-gDNA optimum** rather than the
split — by attacking the seeding/competition:

- **(i) Don't seed nascent without evidence.** A real nascent component has **intron coverage** (the
  spatial check showed nrna_none introns are *empty*); the nascent-vs-mature split of the unspliced-RNA
  prior should be driven by **position** (intron support), so a no-intron-signal locus seeds nascent
  `≈ 0` and the absorbing all-gDNA fixed point holds. *(This is principled prior allocation by geometry —
  not the band-aid "gate nascent" of Option C.)*
- **(ii) Make the gDNA-vs-nascent relative eff-len honest**, so the split basin doesn't deepen
  artificially.
- **(iii) Failing both, a deterministic-annealing / init-from-`f_g` escape** to the global optimum
  (mild, EM still free to move).

**Must verify first (the EM is pristine — do not design on an unconfirmed mechanism):**
1. **Local-optimum vs true-optimum:** initialize `θ_nascent = 0` (or from `f_g`) and check whether the
   leak vanishes (local optimum, seeding) or persists (the default eff-len makes the split the *true*
   optimum). The monotonic eff-len response hints "smooth/true-optimum," but confirm.
2. **The trigger:** measure, on the leaky loci, the nascent component eff-len vs the gDNA IPR and the
   nascent prior share — which one opens the basin.
3. Only then design the minimal, principled change and re-run the flagship + 20-condition net-flow with
   the ss=0.5 / low-gDNA / real-nascent guards.

**Net:** the strand likelihood is correct; Options A and B are both off the table (A void, B refuted by
the "don't fully trust calibration" principle); the leak is an abundance-competition / EM-convergence
problem seeded by the nascent prior and shaped by the relative eff-len. The next step is **empirical
mechanism confirmation**, not implementation.

---

## 9. Empirical mechanism confirmation (2026-06-15) — §8.4-(ii) is the dominant lever; the relative eff-len is dishonest *by construction*

The §8.4 experiments were run (diagnostics: `/tmp/em_basin_measure.py`, `/tmp/em_efflen_predict.py`,
`/tmp/gdna_haircut_factor.py`) on the flagship `gdna300/ss0.99/nrna_none/cap_on` (so **every** nascent
count is phantom leak). Result: **§8.4-(ii) [relative eff-len] is confirmed as the dominant basin-opener.**

### 9.1 The relative eff-len rigs the competition (measurement)

Per leaky locus, the nascent component's EM eff-len vs the gDNA component's IPR eff-len:

- **mass-weighted `nascent_eff / gDNA_IPR = 0.644`, on 100.0% of the phantom-nascent mass** (75 loci,
  182k phantom). The nascent component is systematically ~36% *shorter* than the gDNA component ⇒
  `θ_nascent = count/eff` is inflated ⇒ nascent out-competes gDNA for the shared unspliced gene-body reads.
- The user's *"nascent can never be shorter than mature"* invariant is **mostly satisfied** and is **not**
  the driver: `nascent_eff / mature_eff = 2.13` mass-weighted; only **2.7%** of phantom mass sits where
  nascent < mature. The leak is nascent-vs-**gDNA**, not nascent-vs-mature.

### 9.2 The decisive prediction (§5) — scaling the gDNA eff-len recovers truth

Re-running the per-locus EM with `gdna_eff_len × factor` (RNA eff-lens untouched):

| `gdna_eff ×` | gdna_em | nrna phantom | mrna |
|---|---|---|---|
| 1.000 | 2,751,681 | 88,207 | 1,109,732 |
| 0.800 | 2,818,497 | 57,053 | 1,074,070 |
| 0.644 | 2,877,704 | 32,536 | 1,039,380 |
| 0.500 | 2,939,488 | **15,037** | 995,095 |

Flagship truth (gdna300 = 3:1) ⇒ gDNA ≈ 2.96M, RNA ≈ 0.99M. At `×0.5` the EM lands almost exactly there
(gDNA 2.94M, mrna 0.995M, phantom 15k). **The leak is an eff-len competition** — the gDNA component
eff-len is too *long* relative to the RNA components, and correcting the ratio monotonically collapses the
leak and recovers truth. The historical `×0.5` "band-aid" was approximating this correction.

### 9.3 Root cause: gDNA and RNA eff-lens are built by *different formulas* (different units)

- **gDNA component** (`priors.assemble_priors`): `eff_g = IPR = (G+1)²/[Σ(g²/L)+…]`, capped at span — a
  raw concentration measure in **bp**, **no FL haircut**, **no further span-contraction** (the IPR *is* the
  contraction).
- **RNA components** (`capture_eff_length.transcript_capture_eff_lengths`): `eff_r = fl · [w·(IPR_r/span)+(1−w)]`
  — the **FL-marginal** length `fl ≈ span − E[frag]` (FL haircut) **times** the IPR/span contraction ratio.

So for a shared footprint `eff_nascent/eff_gDNA ≈ (fl/span)·(IPR_nascent/IPR_gDNA)`. Two inconsistencies:

1. **FL-haircut asymmetry** — RNA is FL-haircut, gDNA is not. *Quantified:* the principled fix (FL-haircut
   the gDNA IPR with the gDNA FL pmf, `region_eff_length(IPR, gdna_fl)/IPR`) is only **0.922×**
   mass-weighted — the IPR footprints are large (~6.6 kb), so subtracting the ~384 bp gDNA FL mean barely
   moves it. **FL-haircut alone is a real but minor (~15%) correction — NOT the fix.**
2. **IPR-construction asymmetry [DOMINANT]** — the nascent eff-len is contracted by the IPR over its
   **full span (exons + capture-depleted introns)**, which over-concentrates it toward the probed exons,
   while the gDNA component uses the locus-region IPR. This is most of the 0.644.

### 9.4 Corrected fix direction (supersedes §5's lean toward seeding/Option-C)

**The gDNA component must NOT be FL-haircut** (user correction, 2026-06-15). The FL-marginal haircut
`L − E[frag]` accounts for the **end-effect**: near a bounded region's two ends there are fewer valid
fragment start positions. A **transcript has real ends** (it terminates) → the RNA components correctly
get the haircut. **gDNA has no ends within a locus — it bleeds across transcript and locus boundaries**
(a gDNA fragment can start upstream of a region and extend in), so there is no end-effect and the raw IPR
is the correct gDNA measure. So the §9.3-(1) FL-haircut asymmetry (`0.922×`) is **correct physics, not a
bug** — RNA *should* be slightly shorter than gDNA by the end-effect. Do **not** pursue FL-haircutting the
gDNA (the earlier candidate is retracted).

That leaves the **dominant** §9.3-(2) term as the sole dishonest part: the RNA/nascent eff-len is
over-contracted because it is built as `fl·(IPR/span)` with the IPR taken over the component's **full
span including capture-depleted introns**, collapsing it toward the probed exons. The principled fix is a
**unified enrichment-based effective length** in which every component (mRNA, nRNA, gDNA) derives its
contraction from the **same per-region capture-enrichment ratio** (component density / global density),
with the FL/end-effect convention applied only where it physically applies (RNA ends: yes; gDNA: no). The
user is formulating this enrichment-based strategy; design + prototype to follow with their input. A
uniform `×0.5` is explicitly **unsafe** (a global down-scale manufactures phantom gDNA on zero-DNA
libraries). **This reopens the §5 decision** — the data shows the relative eff-len is the dominant,
measurable, parameter-free lever. **Validation gate for any fix: flagship + 20-condition net-flow with the
ss=0.5 / low-gDNA / nrna_rnd / zero-DNA guards.**
