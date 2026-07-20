# The node solver prior — the "is mature RNA present?" classifier

**Status:** design + implementation plan (branch `calib-ambig-init-wip`, off `v0.7.1`).
**Motivation:** the unstranded + hybrid-capture gDNA→RNA collapse (`ambig_dense_10mb`,
`gdna_*_ss_0.50_*_capture_on`), where enriched exons that are 70–99 % gDNA in truth calibrate to
`f_g ≈ 0` (~13 M-fragment leak). See `scripts/debug/calib_pool_benchmark.py`.
**Supersedes:** the gated `nascent0` grid-MAP + `_NASCENT_TIEBREAK` path committed at `880ac429`
(§7 — why it was the wrong mechanism).
**Reads on:** `CALIBRATION_ARCHITECTURE.md` (count-zero-information, **authoritative**) and
`calibration_initialization.md` (the precision currency). This doc adds the one missing piece those
two do not pin down: **which prior each node's solve uses, and why.**

---

## 0. The whole design in one question

> **"Is unspliced *mature* RNA present in this node's fragment mass — yes or no?"**

That single question selects the solver prior. It is answerable purely from the region annotation
(the 4-bit signature), before any counts are read. Everything below is the consequence.

- **Introns** — mature RNA is spliced out ⇒ unspliced mass is **nascent + gDNA only**.
- **Exons** — mature RNA is present ⇒ unspliced mass is **mature + nascent + gDNA**.
- **Intergenic** — no transcript ⇒ **gDNA only**.

---

## 1. The regular solver already does the strand-specific deconvolution

For a single-strand node the mixture plus-strand probability ([`simplex.py:33`](../../src/rigel/calibration/simplex.py#L33)) reduces to a straight line in `f_g`:

```
p = ½·f_g + κ·f_pos + (1−κ)·f_neg          # general
  = κ − (κ−½)·f_g                           # single-strand (+): f_neg=0, f_pos=1−f_g
```

| `f_g` | composition | model `p` |
|------|-------------|-----------|
| 0 | all (nascent) RNA | **κ** (stranded) |
| 1 | all gDNA | **½** (unstranded) |

The strand likelihood peaks where model `p` matches the observed `p̂ = u_pos/n`:

```
f_g* = (κ − p̂) / (κ − ½)          # "the peel"
```

This **is** the boundary self-solve — no separate estimator. The Fisher information the strand tilt
carries about `f_g` is `∝ (κ − ½)²·N`; it **vanishes at κ = ½**. So:

- **strand-specific data (κ ≫ ½)** → the likelihood dominates and deconvolves nascent-vs-gDNA, exactly
  as it deconvolves RNA-vs-gDNA for a region. Nothing new is needed here.
- **unstranded data (κ → ½)** → the likelihood is flat; the answer is set **entirely by the prior**.

**The prior is the whole game for unstranded data.** That is the one thing the current code gets
wrong at a nascent-only node, and the entire subject of this document.

---

## 2. The count-zero-information leak (why a naive prior swap is not enough)

The Beta-Binomial hands `f_g` information through **two** channels:

1. **The mean** — `p = κ − (κ−½)·f_g` — the legitimate strand tilt (Fisher `∝ (κ−½)²N`, → 0 at κ=½). ✅
2. **The variance normalizer** — `−½·log(var)` with
   `var = … + (n·f_g)²·¼·od_g + (n·f_pos)²·κ(1−κ)·od_r`
   ([`simplex.py:38-45`](../../src/rigel/calibration/simplex.py#L38)) — this term **also** depends on
   `f_g`, and at κ=½ (mean flat) it still tilts `f_g` toward the variance-minimizing value, **scaled by
   count²**. ❌

Channel 2 is a **count-zero-information violation**: the raw count magnitude silently steers the
composition, which `CALIBRATION_ARCHITECTURE.md` forbids. At the flagship's high-count enriched exons
it *overrides* a weak prior and drags `f_g` back toward ½. The solve must take `f_g`-information from
**the mean channel only** — the count enters as *precision*, never as a composition preference.

---

## 3. The node taxonomy — one classifier for regions and boundaries

Per node, per strand `s ∈ {+, −}`, two boolean flags — both **derived from the 4-bit signature**, which
is already computed and persisted at index build ([`signature.py`](../../src/rigel/calibration/signature.py)):

| flag | meaning | region (contained) | boundary (crossing) |
|---|---|---|---|
| `nrna_active_s` | **nascent** RNA on `s` (transcript present — exon *or* intron) | `sig & (EXON_s\|INTRON_s)` | `nrna_active_s(L) ∧ nrna_active_s(R)` |
| `mrna_active_s` | **mature** RNA on `s` (exon) | `sig & EXON_s` | `mrna_active_s(L) ∧ mrna_active_s(R)` |

The names map directly to your "nascent RNA active / mature RNA active": `nrna_active` = nascent (present
wherever a transcript is — introns *and* exons); `mrna_active` = mature (exons only). So
`mrna_active_s ⟹ nrna_active_s` (an exon has both; an intron has only nascent). `nrna_active_s` **is** the
existing transcript-continuity gate `free_s`
([`node_geometry.py:340-341`](../../src/rigel/calibration/node_geometry.py#L340)); `mrna_active_s` is the
**new** bit the code does not yet carry — the whole taxonomy hinges on it.

**Where the logic lives (code-review answer).** Do *not* add new stored index columns — the 4-bit
signature already **is** the build-time flag, and separate columns would be redundant state that can
drift. Put the derivation in one authoritative place: two pure helpers in `signature.py` (the
parameter-free encoding layer) — `nrna_active_strands(sig)` and `mrna_active_strands(sig)`, each
returning `(pos, neg)` masks. `_region_strand_stats` calls them on the region's own signature;
`_boundary_strand_stats` ANDs the two flank regions'. That is the whole wiring — the flag is a one-line
bit-mask, defined once, computed where used.

The three per-strand RNA states and their prior:

| state | condition | RNA in the unspliced mass | prior on strand `s` |
|---|---|---|---|
| **no RNA** | `¬nrna_active_s` | none | strand locked off (`var=0`) |
| **nascent-only** | `nrna_active_s ∧ ¬mrna_active_s` | nascent only (mature would be spliced away) | **nascent≈0 prior** (§4.2) |
| **mature-capable** | `mrna_active_s` | mature + nascent | **region-neutral prior** (§4.1, current) |

The node's overall class falls out of the two strand states:

| both strands | node class | solve |
|---|---|---|
| both `¬nrna_active` | **G1 gDNA sink** | none — locked `{0,0,1}` |
| exactly one `nrna_active` | **single-strand (1-D)** | `_solve_nodes_logodds`, prior per that strand's `mrna_active_s` |
| both `nrna_active` | **AMBIG (2-D)** | `_solve_ambig_logodds`, special init (§5) |

### The four boundary types (the user's taxonomy) are corollaries

1. **intergenic ↔ exon** — the intergenic flank has no transcript ⇒ neither strand `present` ⇒ G1
   gDNA sink. **No solve.**
2. **intron ↔ exon, single-strand** — one strand `present` (transcript on both flanks) but `¬mature`
   (the intron flank has no exon bit) ⇒ **nascent-only** ⇒ 1-D solve with the nascent≈0 prior.
3. **exon ↔ exon** — the shared strand is `mature` (exon both flanks) ⇒ **mature-capable** ⇒ 1-D solve
   with the region-neutral prior (current behavior — correct as-is).
4. **ambig ↔ ambig** — both strands `present` ⇒ AMBIG boundary ⇒ 2-D solve + special init (§5).

**Generalization (confirmed — O1 = yes).** The identical rule classifies *regions*: an **intron region**
is nascent-only, an **exon region** is mature-capable, an **intergenic region** is a gDNA sink. Today the
solver applies the region-neutral (Jeffreys) reference to *all* single-strand nodes including introns —
correct only when strand data resolves them, tippy toward ½ when unstranded. Introns **are included**:
the classifier is applied uniformly to regions and boundaries — this *is* the "general node solver".

---

## 4. The prior

### 4.1 mature-capable → region-neutral (unchanged)

Keep the current path: the Jeffreys Beta(½,½) reference + the change-of-variable Jacobian
([`simplex_logodds.py:174-180,213`](../../src/rigel/calibration/simplex_logodds.py#L174)), posterior
**median** read-out, plus the sweep's global gDNA prior. Mature RNA is *expected* in an exon; a neutral
RNA-vs-gDNA default is right and must not be disturbed.

### 4.2 nascent-only → the density-anchored nascent≈0 prior

The physical statement: *the unspliced mass is gDNA up to the gDNA background density; nascent RNA (any
excess) is a priori ≈ 0.* This is **not a new prior term** — it is the **global gDNA prior**
(count-zero-info's third source) applied *alone*, with the neutral Jeffreys reference removed.

Its behavior, in one closed anchor (`ρ_bg` = the expected gDNA density for the node's stratum, `E_gdna`
= the gDNA-FL effective length, `M` = observed unspliced mass):

```
f_g^anchor = clip( ρ_bg · E_gdna / M , 0, 1 )
```

- `ρ_obs = M/E_gdna ≲ ρ_bg` (density at or below background) → `f_g → 1`. **Nascent≈0 emerges** — it is
  a *consequence* of the density prior, not a hand-set default. ✓
- `ρ_obs ≫ ρ_bg` (density far exceeds background) → `f_g < 1`; the excess mass cannot be gDNA, so the
  prior **yields** and the data (strand tilt and/or density) dictates. ✓ — exactly the user's *"high
  density → let the prior disappear."*

**Strength is derived, not magic.** The prior is `log ρ_g ∼ 𝒩(log ρ_bg, σ²_bg)` with `ρ_bg` and the
spread `σ²_bg` fit from the **intergenic/intron density distribution** — the population baseline that
already anchors the global gDNA prior. No unit-information constant, no tunable ε. This *retires* the
`_NASCENT_TIEBREAK` question of §7 and the "magic strength" worry entirely: the intergenic distribution
sets the strength.

### 4.3 Why density is allowed but count is not (count-zero-info reconciliation)

`CALIBRATION_ARCHITECTURE.md` forbids the raw **count** from tilting composition; it explicitly permits
the **global gDNA prior**. Density = `M / E_gdna` is a *normalized* quantity compared against the gDNA
density *distribution*: 1000 fragments in a tiny region (high density) and 1000 in a huge region (low
density) imply different gDNA-ness precisely because gDNA has a characteristic background density. That
is the global prior doing its job — the third sanctioned source — not a count leak.

### 4.4 Capture: the prior must be BIMODAL, and available at init (the pass-0 KDE)

Under hybrid capture the gDNA background is **not** unimodal. Intergenic and intron nodes stay depleted
(the background prior fits them well), but **captured boundary nodes are enriched** — their density
greatly exceeds the depleted background. A single depleted `ρ_bg` then reads that excess as "not gDNA ⇒
nascent RNA" and **spawns nascent RNA at every enriched boundary** — the flagship collapse.

Two defenses, and when each applies:

1. **Strand + count (stranded data).** When the library is stranded, the strand balance says "no nascent
   RNA here" (unstranded-looking counts ⇒ gDNA), and the mean channel (§1) overrides the density-excess
   spawn. Stranded capture is largely self-correcting.
2. **A bimodal gDNA prior (the only defense for *unstranded* capture).** With no strand information the
   density-excess spawn is unchecked. The cure is a gDNA prior with **two modes** — a depleted mode
   (intergenic/introns) and an **enriched mode** (captured boundaries) — so an enriched boundary is
   scored against the enriched gDNA level and stays gDNA.

**We already build exactly this bimodal prior — the KDE** (`gdna_density_prior.py`): a fitted density map
that applies the prior in a density-weighted (bimodal) fashion against each node's observed density. The
gap is *timing*: the KDE is a **Phase-2 (sweep)** object today, and we need it **at init**. Call this the
**pass-0 gDNA KDE prior** effort (§8 Phase C). It is feasible without a prior solve because the
enriched-clean training sample is selectable from the **strand-free** boundary landscape (x = unspliced ≈
gDNA, y = spliced ≈ RNA → pick enriched, RNA-free intron↔exon boundaries → their adjacent exons are a
clean enriched gDNA sample), as validated in `memory: ambig_g1emit_stratified_kde_2026_07_12`. This is
the dominant lever for the flagship; the §3 classifier + the §2 count-zero-info fix are necessary but
not sufficient without it.

---

## 5. AMBIG boundary — the special initializer

At an AMBIG node both strands are live, so the mixture `p` depends only on the RNA **tilt** `τ` and is
**blind to `f_g`** (the magnitude is strand-invisible — `simplex.py` docstring). The strand likelihood
therefore contributes *nothing* to `f_g`; it is set **entirely by the prior**. So the AMBIG initializer
is simply §4.2 applied to the magnitude: `f_g = f_g^anchor` from the density prior (→ 1 at background,
yielding at excess), with the tilt marginalized and the honest count precision. The sweep then refines
it from neighbour messages. This is the "special initializer" — no strand shortcut exists.

---

## 6. How value and precision combine (robust, no tricks)

On the log-odds grid, per single-strand node:

```
ψ(λ) = mean-channel strand loglik(λ)         # channel 1 only (§2); count enters as precision
       + density prior on log ρ_g            # §4.2 for nascent-only; Jeffreys+global for mature-capable
f_g   = posterior MEDIAN over λ              # transform-invariant, robust to skew
Var(log ρ_c) = 1/N_c,  N_c = f_c·M           # honest per-component count precision (φ=0)
```

At κ=½ the likelihood is flat ⇒ the posterior = the density prior ⇒ `f_g → f_g^anchor` (→ 1 at
background). As κ moves off ½ the mean channel sharpens and overrides the prior — the crossover lands
where `(2κ−1)²N ∼ 1/σ²_bg`. This is the graceful shrinkage the user asked for ("nascent → 0 **only** for
unstranded"), and it is a plain Bayesian posterior — **no MAP, no `_NASCENT_TIEBREAK`.** The
confidence-weighted peel validated in `scripts/debug/boundary_ab.py` (flagship single-strand
`f_g → 0.994`) is the closed-form limit of this posterior; the grid computes it exactly.

---

## 7. What this supersedes (`880ac429`)

The committed gated `nascent0` path reads the **MAP** with a `−ε·(1−f_g)` tie-break to force the f_g=1
vertex ([`simplex_logodds.py:204-207,269-279`](../../src/rigel/calibration/simplex_logodds.py#L204)).
It is **not robust at real κ ≈ 0.499**: the near-flat likelihood + the max-of-sides `u_pos≠u_neg`
imbalance divided by the ~0 slope tips the MAP to an extreme (`scripts/debug/density_frame.py` used
*exactly* κ=0.5, which hid this). The fix is not a better tie-break — it is §4.2 (a real density prior)
+ §6 (posterior median, not MAP) + §2 (mean-channel only). The `nascent0` grid-MAP path and
`_NASCENT_TIEBREAK` are **to be removed**, not extended.

---

## 8. Implementation plan (phased; each phase independently testable)

**Guardrails.** Develop/validate on toys + `boundary_ab.py` before the big suite
(`memory: calibration_feature_dev_toy_scenarios`, `ambig_needs_robust_training_subset`). Any shared-code
change (the count-zero-info freeze) gets a **region A/B** to prove no regression. **Verify against the
oracle at every calibration step, not just end-to-end** (`scripts/debug/oracle.py`): the truth is known, so
each stage — post-init, post-sweep, post-projection — is compared to it, and improving initialization must
land the post-init state measurably *closer to truth* (else the change is not earning its keep). Any
behavior-changing phase regenerates the goldens (`--update-golden`) with the diffs reviewed for direction.

### Phase A ✅ DONE — carry the `mrna_active_s` classification (cheap, pure index arithmetic)
- Added two pure helpers to `signature.py`: `nrna_active_strands(sig)` and `mrna_active_strands(sig)` →
  `(pos, neg)` masks — the one authoritative definition (§3).
- `_region_strand_stats` / `_boundary_strand_stats`
  ([`node_geometry.py`](../../src/rigel/calibration/node_geometry.py)): compute `mrna_active_pos` /
  `mrna_active_neg` (region: own signature; boundary: AND of the two flanks) alongside the existing
  `free_pos` / `free_neg` (which are `nrna_active`); `_POS_BITS`/`_NEG_BITS` retired into the helpers.
- `NodeStatics` + `build_node_statics`: carry `mrna_active_pos` / `mrna_active_neg` onto the chain.
- Tests: the classifier + the 4-boundary taxonomy in `test_signature.py`; a `build_node_statics`
  integration test (`test_node_prior_classifier.py`) over every region + boundary type.
- **Landed behavior-neutral:** goldens bit-identical, full suite 1095 passed. `mrna_active` carried but
  not yet consumed by any solver.

### Phase B — split into B1 (the freeze, isolated) and B2 (the prior selection, with the density)

Phase B bundles three changes of different risk and coupling. They are split so the shared-code
correctness fix lands in isolation, and the prior selection lands *with* the density that makes it
magic-number-free.

#### Phase B1 — the count-zero-info freeze (isolated; **no new constant**)
- Add `f_g_ref, f_pos_ref, f_neg_ref` params to `_mixture_strand_loglik`
  ([`simplex.py:25`](../../src/rigel/calibration/simplex.py#L25); the 2 call sites are both in
  `simplex_logodds.py`): the model **mean stays live** (grid `f_g`), the **variance is frozen** at the
  reference (§2, §9). Thread the reference through `_local_loglik_logodds` / `_solve_ambig_logodds` /
  `_solve_nodes_logodds_all`.
- **Reference source:** `init_beliefs` passes a structural-neutral reference (`f_g=0.5` split among the
  *live* strands); `node_sweep._local_solve` passes the **incoming belief** `(f_g, f_pos, f_neg)`, so the
  variance — hence the message precision — is evaluated near the truth and stranded ripple stays small.
  (The reference is low-stakes for the *location* but sets the *precision*, so sourcing it from the belief
  rather than a flat ½ matters — `frozen_variance_sandbox.py` Test C.)
- **Delete** `_NASCENT_TIEBREAK` and every `nascent0` branch — superseded (§7), and referenced only in
  `simplex_logodds.py` (no tests, no production callers).
- No prior-selection change yet: nascent-only nodes keep Jeffreys; `mrna_active` (Phase A) stays
  carried-but-unconsumed.
- **Validation:** region A/B on stranded scenarios (expect ~neutral — the mean channel dominates);
  **post-init AND post-sweep compared to the oracle** (`scripts/debug/oracle.py`) — improving the leak must
  move the post-init state closer to truth; goldens regenerated, diffs reviewed for the leak-fix direction.
  NB the freeze is *not guaranteed small* — it shifts message precisions — which is exactly why the region
  A/B is mandatory, not a formality.

#### Phase B2 — the nascent≈0 prior selection (ships with Phase C)
- Consume `mrna_active`: `nascent_only_s = nrna_active_s ∧ ¬mrna_active_s`. Nascent-only nodes drop the
  Jeffreys reference and take the density-anchored prior (§4.2), reading the posterior **median** (§6).
- **`_type_belief` surgery** ([`node_geometry.py`](../../src/rigel/calibration/node_geometry.py)): the
  nascent-only case cuts across the current G2 (single-strand) class, which today applies the plain strand
  solve — it must instead apply the nascent≈0 prior. This is the risky init-structure change, kept out of
  B1.
- **Why merged with C, not standalone:** the nascent≈0 prior is magic-number-free *only* because its
  strength comes from the density (the intergenic distribution / the pass-0 bimodal KDE). Selecting the
  prior before the density anchor exists would force a temporary unit-information constant — the magic
  number we agreed to avoid. So B2 and Phase C ship together.
- **Persistent coupling:** init `f_g` feeds `_gdna_seed_estimate` / `_floor_estimate`
  ([`bp_solver.py:665,671`](../../src/rigel/calibration/bp_solver.py#L665)) → the global prior → the whole
  sweep, so B2 cannot be validated at init alone — compare post-init AND post-sweep to the oracle.

### Phase C — the pass-0 gDNA KDE prior (ships with B2; the flagship lever)
- Fit the **bimodal** gDNA KDE (`gdna_density_prior.py`) **before init**, selecting the enriched-clean
  training sample from the strand-free boundary landscape (no prior solve needed — §4.4). Apply it as the
  init density prior so enriched boundaries anchor to the enriched gDNA mode, not the depleted one.
- Distinct from Phases A/B (the classifier + solver mechanics). Builds on the stratified enriched-mode KDE
  work (`memory: ambig_g1emit_stratified_kde_2026_07_12`, currently uncommitted on `main`).
- **Gate (before committing Phase C):** extract the enriched-clean training-sample selection into a
  standalone diagnostic and plot the density histogram of the selected nodes across all 24 `ambig_dense`
  conditions (+ a real low-coverage run). The bimodal split must be visibly clean in practice, not just in
  theory. Where the enriched sample is thin the enriched mode must **degrade gracefully to the depleted
  background** (= current behavior, no worse) and gate off — it slightly *hurt* stranded (−0.9 %) in prior
  work, so it must not fire where the strand channel already resolves the enrichment.
- Test: `calib_pool_benchmark.py` on the flagship (`gdna_gdna300_ss_0.50_*_capture_on`) — target the
  ~13 M leak; **perfect self-gating** on `gdna_none` / `capture_off`; confirm stranded capture is not
  regressed.

### Phase D — AMBIG boundary initializer (§5)
- Route AMBIG nodes' `f_g` to `f_g^anchor` (density prior only; tilt marginalized) in
  `_solve_ambig_logodds`; verify the sweep refines from messages.
- Test: AMBIG-dense benchmark (`ambig_dense_10mb`, 305 AMBIG regions) with ample single-strand nodes
  present (`memory: ambig_needs_robust_training_subset`).

---

## 9. Decisions (resolved) + the count-zero-info tradeoff

- **O1 — generalize to intron regions? → YES.** The classifier is applied uniformly to regions and
  boundaries (§3); intron regions are nascent-only.
- **O2 — count-zero-info fix everywhere? → YES, with eyes open.** See the tradeoff below — the current
  design is *not* an accident, and we should preserve what it buys while removing the leak.
- **O3 — sparse capture → the pass-0 bimodal KDE (§4.4, Phase C).** The depleted background alone spawns
  nascent RNA at enriched boundaries; the bimodal KDE is the fix, needed at init. Remaining unknown: how
  the enriched mode's weight shrinks toward the depleted mode when the enriched-clean sample is thin (must
  be derived from the sample size, not tuned).

### The count-zero-info tradeoff (O2 — why the current design is not an accident)

`_mixture_strand_loglik` is a **correct heteroscedastic Gaussian approximation to the Beta-Binomial**:
mean `N·p`, variance = binomial + per-component overdispersion, and `−½·log(var)` is the proper Gaussian
normalizer. Making the variance depend on composition is textbook-correct — gDNA (`od_g`) and RNA
(`od_r`) have genuinely different overdispersions, so a given strand count is more/less surprising
depending on the mix. That is real statistical information; whoever wrote it wrote a *locally correct*
likelihood.

The pathology is narrow and bites only when the **mean is degenerate** (κ→½): there the mean channel
carries no `f_g` information, yet the variance normalizer still discriminates compositions — so the count
magnitude (via `n²`) leaks a composition preference toward the variance-minimizing `f_g`. Invisible until
unstranded data; hence "not an accident," just an untested regime.

**The fix is not "drop the variance" — it is "freeze the variance's composition argument at a reference."**
Evaluate `var` at a fixed reference composition `f_ref` so it is a function of the count + `f_ref`,
**not** of the `f_g` being solved. Critically, **only the variance freezes** — the model mean in the
numerator `(u_pos − N·p(f_g))²` stays **live** (that is the legitimate mean channel; at κ=½ that mean is
flat in `f_g` anyway, which is exactly why the location goes prior-driven there). This is standard
Fisher-scoring / expected-information practice, and it:
- **keeps** the heteroscedastic precision (the good part — the count still sets a composition-aware
  precision), and
- **removes** the `f_g`-tilt from *both* the `−½·log(var)` and the `(u−mean)²/var` terms (the leak).

**The reference is the prior anchor, not flat ½.** For a nascent-only node `f_ref` = the nascent≈0 anchor
(`f_g ≈ 1`, near the truth for pure gDNA); for a mature-capable node the neutral. This evaluates the
variance near the truth with no iteration and avoids the bootstrap failure of a flat-½ reference (which
over-states the overdispersion for near-vertex nodes and mis-sets their precision). The reference choice
is low-stakes for the *location* — at κ=½ the median is prior-driven regardless of `f_ref`
(`scripts/debug/frozen_variance_sandbox.py` Test C: spread 0), and at κ≫½ it matters only **modestly**
(Test C: ~0.04 median shift at `od=0.2`, from the wide overdispersed posterior — modest, *not* fully
second-order, which is exactly why the prior anchor rather than flat ½ is the right reference). An optional
one-step re-evaluation of `var` at the solved median is cheap insurance for high-overdispersion nodes.

**Acceptance test (the freeze is correct iff), validated in `scripts/debug/frozen_variance_sandbox.py`:**
at κ=½, sweeping the count `n` from 10 → 10 000 must **sharpen** the posterior (raise precision) but **not
move** the median `f_g`, and the count must no longer pull `f_g` toward the variance-minimum
`od_r/(od_g+od_r)`; at κ≫½ the median must still track the mean-channel peel `(κ−p̂)/(κ−½)` (no stranded
regression). Validate in the sandbox before touching the shared `_mixture_strand_loglik`, then run the
region A/B.
