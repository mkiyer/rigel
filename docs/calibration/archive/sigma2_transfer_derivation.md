# The message transfer-variance model — a derivation, proven belief-free on the oracle

**Status:** derivation + empirical proof, written 2026-07-16 (branch `calib-ambig-init-wip`), the session that
followed the mature-gate dismantle. This is the load-bearing §6.1 of
[`message_model_derivation.md`](message_model_derivation.md): what `σ²_transfer` **is**, how to estimate it
**without a solve** (the circularity firewall), and the one design question the owner posed — **source-only vs
pair** — settled on data. **Companions:** [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md)
(count-zero-information — authoritative), [`npmle_roadmap.md`](npmle_roadmap.md) (the NPMLE gDNA prior + the KEY
NEGATIVE RESULT this design must not repeat), [`archive/imputation_variance_model.md`](archive/imputation_variance_model.md)
(the Poisson/σ²_bio decomposition this specializes). **Evidence:** `scripts/debug/transfer_variance_diag.py`
(the oracle measurement) and `scripts/debug/sigma2_transfer_regime_proto.py` (the belief-free-stratifier proof);
both run on the 7-condition `_selfsolve_cache`.

> **The one-line result.** `σ²_transfer,g` — the missing third term of the message precision — is **not** the
> total-density adjacent disagreement (that is 12–28× too large off-capture, all RNA-expression noise), and it
> is **not** a source-only property (source-only blends the reliable same-mode edges with the crossing and gets
> both wrong). It is a **pair** quantity, **stratified on the belief-free total-density regime of both
> endpoints**, and fit **belief-free on RNA-free anchor pairs**. The NPMLE-on-total-density the owner already
> ships is exactly the regime classifier this needs; the transfer variances are read from adjacent RNA-free
> pairs, **not** from the NPMLE's marginal mode width (which over-states adjacent spread by an order of
> magnitude). This breaks the circularity: the regime label and `σ²_transfer` are both belief-free and **fixed
> across passes**; only the prior and the beliefs evolve.

---

## 1. The gap, precisely

The message from a source node about a destination component `c` carries precision
```
    w_c  =  1 / ( Var(log f_c^src)   +   1/M_src   +   σ²_transfer,c )
             ╰── composition ──╯    ╰ magnitude ╯   ╰─ transport ─╯
```
The first two terms are honest and shipped (the strand Fisher info folded into `Var(log f_c^src)`, and the
Poisson magnitude `1/M_src`; `bp_solver._scan`, `pr = sm/(sm·v + 1)` = `1/(v + 1/M)`). **The third term is
zero.** With `σ²_transfer = 0` the precision saturates at `1/Var(log f_c^src)` — for a dense, confident source
that is **thousands** of pseudo-observations (`message_model_experiments.py` E1: an exon message at precision
2453 drags a boundary that self-solved correctly to 0.963 down to 0.81). The message has **no per-hop decay**:
a dense node's reach is unbounded, so it overrules a sparse neighbour's correct self-solve. `1/σ²_transfer,c` is
exactly the missing **cap** — the maximum confidence a message may carry no matter how many fragments the source
counted. This is the count-zero-information principle at the message layer: *a high count buys statistical power,
never accuracy; the accuracy is capped by how much the world's true rate can drift across the edge.*

`σ²_transfer,c` is a property of the **world** — "how much does the true log-`ρ_c` differ between adjacent
nodes?" — not of our knowledge. It is **count-independent** by construction (it is the residual dispersion after
the Poisson `1/M` term is removed).

---

## 2. Why the shipped estimator is the wrong quantity (measured)

The production scalar `bp_solver.adjacent_disagreement_variance` measures the **total-density** adjacent
disagreement. `transfer_variance_diag.py` measures, against the oracle, what that scalar reports vs the **true
gDNA** transfer variance:

| condition | shipped σ²_imp (total) | TRUE gDNA σ²_transfer,g | ratio |
|---|---|---|---|
| ss0.50 none **capOFF** | 1.016 | **0.036** | **28×** |
| ss0.99 none capOFF | 1.019 | 0.040 | 25× |
| ss0.50 present capOFF | 0.865 | 0.064 | 12× |
| ss0.50 none **capON** | 3.92 | 6.70 | ~1× |

Off-capture the total-density estimator is **12–28× too large**: gDNA is genomically smooth
(`σ²_transfer,g ≈ 0.04`), but the total density carries `σ²_r ≈ 17` of pure **RNA-expression** variation, and
the total-disagreement scalar reads that RNA noise as if it were gDNA transport noise. **This is what gags every
message off-capture** — the "impotent message precision" the owner described. Under capture total ≈ gDNA (gDNA
dominates the enriched mode), so the scalar is accidentally ~right there, but a single scalar cannot express
what §3 shows the truth to be.

**Corollary (the fit set):** because total-density disagreement = gDNA + RNA and the RNA half is the larger and
regime-dependent, `σ²_transfer,g` **cannot** be read from total density directly. It must be measured where
RNA is absent, so that `ρ_total ≡ ρ_g` identically — the **RNA-free anchor**. On the anchor the raw and
belief-free numbers coincide with the oracle gDNA value (capOFF anchor 0.041 vs true 0.036) — the belief-free
estimator is validated.

---

## 3. The transfer variance is STRATIFIED — the central fact

Under capture the true `σ²_transfer,g` is not one number. Stratifying adjacent pairs by whether each endpoint is
in the **depleted** or **enriched** density regime (`transfer_variance_diag.py`, oracle gDNA regime):

| condition | dep–dep | enr–enr (bio) | **MIXED (crossing)** |
|---|---|---|---|
| ss0.50 none capON | **0.00** | 1.63 | **25.0** |
| ss0.99 none capON | 0.00 | 1.52 | 25.5 |
| ss0.50 present capON | 0.00 | 2.49 | 24.0 |
| ss0.50 none **capOFF** | 0.00 | 0.00 | 0.006 |

Read this as message caps `1/σ²`:

* **dep–dep: `σ²=0` ⇒ cap `=∞`.** The depleted gDNA background is *genuinely uniform* — adjacent depleted
  nodes differ only by counting noise. A depleted source may message a depleted neighbour at **full** source
  precision. This is what lets the intergenic/intron background propagate a reliable gDNA level.
* **enr–enr: `σ²≈1.6` ⇒ cap `≈0.6`.** Enriched nodes vary (probe efficiency differs), so they inform each other
  only weakly.
* **MIXED / crossing: `σ²≈25` ⇒ cap `≈0.04`.** At a probe boundary the depleted and enriched sides have
  genuinely different rates; the message is **essentially off**. This is correct — the two sides *should not*
  cross-contaminate. (It is also exactly the E1 pathology, cured: the dense enriched exon can no longer drag its
  depleted-boundary neighbour.)
* **capOFF: flat ≈0.** No modes, no crossings; one small number for every stratum. Messages flow, capped at
  `1/0.04 ≈ 25` (enough to prevent the "precision 2453" runaway, loose enough not to gag).

The dynamic range is **~25/0.006 ≈ 4000× within capON, and ~25× between capON-crossing and capOFF** — no scalar,
and no source-only law, can represent it (§5).

---

## 4. The circularity, and how the regime label breaks it

To stratify we need a regime label per node; the regime is a gDNA-density property; the gDNA density needs the
solve; the solve needs the messages; the messages need `σ²_transfer`. That is the circularity the owner named,
and the reason the *previous* honest attempt backfired (`npmle_roadmap.md`: `message_precision.
adjacent_imputation_variance`, fit on the **not-yet-peeled** belief — adjacent wrong nodes agree ⇒ `σ²` small ⇒
confident messages propagate the error; *"honesty measured against a wrong belief is not honesty."*).

**The break: the regime label is read from the BELIEF-FREE total density, and it is a good enough proxy for the
gDNA regime exactly where the stratification matters.** Total density `ρ_t = M/E` (asserting `f_g=1`) is the
naive, assumption-free landscape — the same substrate the shipped NPMLE (`GdnaRatePrior`) already fits at
pass-0. Its modes are the regimes. The concern is that `ρ_t = ρ_g + ρ_r`, so a depleted-gDNA but high-RNA exon
is mis-labelled "enriched". `sigma2_transfer_regime_proto.py` measures whether that mis-label costs anything, on
**all** edges (apply-time, RNA-bearing exons included):

| condition | regime agreement (all / exon) | gDNA σ² by **oracle** regime | gDNA σ² by **belief-free** regime |
|---|---|---|---|
| ss0.50 none **capON** | 0.956 / 0.893 | [0, 1.63, **25.0**] | [0, 1.19, **25.1**] |
| ss0.99 none capON | 0.960 / 0.906 | [0, 1.52, 25.5] | [0, 1.28, 25.2] |
| ss0.50 present capON | 0.816 / 0.841 | [0, 2.49, 24.0] | [0, 2.09, 25.1] |
| ss0.50 none **capOFF** | 0.841 / 0.700 | [0, 0, 0.006] | [0, 0, 0.003] |
| ss0.50 present capOFF | 0.617 / 0.768 | [0.004, 0, 0.032] | [0.012, 0.016, 0.038] |

Two regimes emerge, and both are safe:

* **Under capture** the belief-free regime agrees with the oracle gDNA regime **82–96 %** and the stratified
  `σ²` lands in the same place (ordering exact, magnitudes within ~20 %). Why: capture enriches gDNA at the same
  probes that enrich total density, so `ρ_t`-enrichment and `ρ_g`-enrichment **coincide**. The stratifier is
  accurate exactly where `σ²` is large and stratified.
* **Off capture** the agreement falls (RNA contaminates `ρ_t`), but every stratum is `< 0.04` — the mis-label is
  **benign** because `σ²_transfer,g` is flat there. The worst case (present, capOFF, 60 % agreement) still
  produces caps `≥ 25` everywhere.

So the belief-free total-density regime is a self-limiting classifier: it is accurate when it must be (capture,
bimodal landscape) and harmless when it is not (off-capture, unimodal). **The regime label — and therefore
`σ²_transfer` — is belief-free and fixed across all passes.** The circularity is cut.

---

## 5. Source-only vs pair — settled on data

The owner's principled instinct was source-only ("look up the source on the NPMLE; the message shouldn't depend
on the destination"). The data refutes it. `sigma2_transfer_regime_proto.py`, source-only blend (the variance
over all edges leaving a source of a given regime):

| condition | dep-source blend | enr-source blend | the pair truth it must replace |
|---|---|---|---|
| ss0.50 none capON | **20.8** | 3.71 | dep–dep **0**, enr–enr 1.6, crossing 25 |
| ss0.50 present capON | 20.4 | 4.20 | dep–dep 0, enr–enr 2.5, crossing 24 |

A depleted source's true reliability is **bimodal**: `σ²=0` to a depleted neighbour, `σ²=25` to an enriched
one. Source-only collapses that to one number (20.8) — which **gags the reliable dep–dep messages** (the very
ones that propagate the background) *and* is still too confident for some crossings. The blend is dominated by
whichever edge type is more numerous, which is an artifact of geometry, not physics. **Source-only is
structurally incapable of separating the same-mode edge from the crossing, and that separation is the entire
signal (0 vs 25).**

The pair model is the honest one, and it does **not** violate count-zero-information: the destination enters
only through its **belief-free regime label** (an observable that sets *reliability*), never as a composition
vote. It is "disagreement-aware", but via the **population-fit per-stratum variance**, not the single-edge
disagreement (which is the 1-sample-noisy, self-fulfilling estimator that caused the historical `base_var`
runaway). The destination's total density decides *which curve to read*, not *what the answer is*.

**One correction to the source-only proposal that survives:** the NPMLE's *marginal mode width* is **not**
`σ²_transfer`. The enriched mode may span ~2 decades of population spread while *adjacent* enriched nodes differ
by only `√1.6 ≈ 1.3` nats — adjacent nodes are spatially correlated. `σ²_transfer` must be measured on
**adjacent RNA-free pairs**; the NPMLE supplies the *mode membership* (the regime label), not the variance.

---

## 6. The model (recommended)

Per component `c ∈ {g, r₊, r₋}`, count-independent, fit once, belief-free, fixed across passes:

**6.1 gDNA leg.** Classify every node's regime from the belief-free total-density NPMLE (`GdnaRatePrior` modes):
`regime(node) ∈ {depleted, enriched}` = which mode `ρ_t = M/E` falls in (the antimode/valley is the boundary —
no magic median; a unimodal landscape ⇒ one regime ⇒ one small `σ²`). Fit
```
    σ²_transfer,g[ regime_src , regime_dst ]  =  Poisson-removed Var( log ρ_g^dst − log ρ_g^src )
```
over **RNA-free anchor** edges (`ρ_t ≡ ρ_g`): intergenic regions (structurally RNA-free) + RNA-free single-
strand SJ boundaries (the `dual_enriched` set of `npmle_fusion.py`, covering the enriched mode) + introns with
no nascent. Three strata (dep–dep, enr–enr, mixed); `mixed` is symmetric in src/dst. Measured values:
`≈ [0, 1.6, 25]` on-capture, `≈ 0` off-capture — **read from the data, not tuned**.

**6.2 RNA leg.** `σ²_transfer,r` is large everywhere (`rna_raw ≈ 8–17` across all conditions — expression is
spiky, as the architecture says). Fit it on deep constitutive exon segments (RNA-dominated, `ρ_t ≈ ρ_r`), or
default it conservatively high. Consequence: RNA messages are **weak** (cap `≈ 0.06–0.12`) — a nudge, not a
hammer — which is correct: nascent imputation should be gated **structurally** (`free_s`), and where it fires it
is a gentle discrepancy signal (the §6.2 nascent factor of `message_model_derivation.md`), never a confident
pin.

**6.3 Apply.** In `_scan`, replace `pr = sm/(sm·v + 1)` with `pr = 1/(v + 1/M_src + σ²_transfer,c[regime_src,
regime_dst])`. The message-free self-solve then becomes a **floor** the messages can only improve on (dep–dep
`σ²=0` still lets a certain source teach; the crossing `σ²=25` caps the message below the self-solve so it
cannot break it).

**6.4 The roadmap ordering (the circularity resolved).**
1. Fit the NPMLE on belief-free total density — **shipped** (`GdnaRatePrior` pass-0). Its modes are the regime
   classifier.
2. Self-solve all nodes (`init_beliefs`: strand + reference + NPMLE prior) — **shipped**.
3. Fit `σ²_transfer` belief-free on the RNA-free anchor (gDNA) + deep exons (RNA), stratified by the regime
   label from (1) — **NEW, this document's deliverable**.
4. Forward-backward sweep with the honest stratified `σ²_transfer` in the message precision.
5. Refit the NPMLE prior on the peeled belief (`calib_refit_iters`, shipped) — but **not** `σ²_transfer` or the
   regime label; those stay belief-free and fixed. Only the prior + beliefs evolve.

Nothing in 3–5 reads a solved belief into a message precision. The regime label (belief-free total density) and
`σ²_transfer` (anchor + deep-exon fit) are both fixed inputs. This is the owner's "workflow 2" — fit the
landscape first — with the landscape used as a **regime classifier for a belief-free variance fit**, not as a
prior that would need the solve to interpret.

---

## 7. Magic-number ledger

**No new tuned constants.** The `σ²_transfer` values are **measured** per stratum (a reliability, not an ε).
The regime boundary is the NPMLE **antimode** (derived, not a median). The number of strata (3) is the minimal
faithful discretization of "same-mode vs crossing"; a continuous alternative — `σ²_transfer,g` as a smooth
function of the NPMLE landscape distance between src and dst — is strictly cleaner and carries **zero** discrete
knobs (deferred; the 3-stratum form is the proven MVP). The anchor fit-set is structural (signature-defined),
not a threshold.

---

## 8. Established vs assumed

* The total-density estimator is 12–28× too large off-capture (pure RNA noise) — **measured** (§2).
* `σ²_transfer,g` is stratified `[0, 1.6, 25]` on-capture, flat off-capture — **measured** (§3).
* The belief-free total-density regime recovers the oracle gDNA regime where `σ²` is large, and its errors are
  benign where `σ²` is flat — **measured, all 7 conditions** (§4).
* Source-only conflates same-mode with crossing (0 and 25 → one blend) — **measured** (§5).
* The RNA-free anchor recovers the true gDNA transfer variance belief-free — **measured on capOFF** (§2);
  under capON the anchor brackets the truth (4.99 vs 6.70 raw) — **measured, weaker** (Poisson over-subtraction
  on fractional-mass counts; the *ordering/agreement* claims do not depend on the absolute Poisson removal).
* The belief-free **structural** anchor (intergenic + RNA-free SJ boundaries + no-nascent introns) covers both
  regimes with enough pairs to fit all three strata — **ASSUMED**; the oracle-RNA-free proof used the oracle
  mask for the measurement set. **Verify next**: does the structural anchor reproduce §3's strata? (`dual_enriched`
  already identifies the enriched RNA-free set, so this is expected, not a blocker.)
* Plugging honest stratified `σ²_transfer` into the sweep improves accuracy without breaking high-nascent cases
  — **NOT yet run** (§6.3 is the wiring; this document is derivation + belief-free-estimator proof only).

---

## 9. Verification update (2026-07-16) — the structural-anchor correction (`anchor_verify.py`)

Running the §8 "verify" item surfaced a design correction, not a confirmation. Three findings:

1. **The structural RNA-free anchor only covers the DEPLETED stratum.** Requiring *both* endpoints to be
   belief-free RNA-free (intergenic region OR RNA-free SJ boundary) yields **zero** adjacent pairs in the
   enr-enr and MIXED strata (all 7 conditions, `nB=0`). RNA-free nodes are adjacent to each other only in
   intergenic stretches; an enriched RNA-free SJ boundary is surrounded by RNA-bearing exon/intron regions, so
   it never forms an RNA-free *edge*. **The direct belief-free response `Δ log ρ_g` exists only where RNA is
   absent, and RNA is absent only in the depleted background.** This is fundamental, not a mask bug: off the
   depleted background there is no `f_g`-independent gDNA-density difference to regress on (the same wall
   `honest_precision_unified_design.md` §3.4 hit for the RNA channel).

2. **Total density reproduces the gDNA strata under capture, and only under capture.** `C_total` (the
   total-density disagreement on *all* edges) matches the oracle gDNA strata on-capture (`[0, 1.3, 25]` vs
   `[0, 1.4, 25]`) because gDNA dominates total there (`total_over_gdna ≈ 1`); off-capture it is RNA-
   contaminated (enr-enr 1.37 vs true 0), because there total = RNA-expression. So the enriched/crossing
   `σ²_transfer,g` is inferrable from total density **iff gDNA dominates total** — a capture property.

3. **Coverage-depth refutes the sampling-bias worry but demands a continuous form.** Splitting each stratum by
   coverage: capON MIXED `σ²` is **0.98** on the low-coverage half vs **25.7** on the high-coverage half. The
   high `σ²` is driven by *well-covered* crossings (real density jumps at probe edges), **not** low-count noise —
   so the strata are structural physics, not Poisson artifact. But the 26× within-stratum range shows the binary
   regime is coarse: the true driver is the *magnitude* of the density gap, which a **continuous** `σ²_transfer,g`
   as a function of the belief-free total-density gap captures and a 3-step lookup does not.

**The corrected belief-free estimator (two pieces):**
* **Depleted baseline** `σ²_transfer,g ≈ 0.04`: measured directly on the intergenic RNA-free anchor (both
  endpoints truly RNA-free; the only place a belief-free gDNA response exists). RNA-immune, always available.
* **Enriched / crossing**: inferred from the total-density disagreement, **gated on a belief-free capture-
  enrichment detector** = median(on-target total density) / median(intergenic total density) (both belief-free;
  the intergenic level is structural — this is the separation the NPMLE already renders). Large ⇒ capture ⇒
  gDNA dominates total ⇒ the total-density strata are valid. ≈1 ⇒ no capture ⇒ fall back to the flat depleted
  baseline everywhere (the off-capture truth is flat ≈0.04 regardless).

This replaces §6.1's "fit on the RNA-free anchor stratified by regime" (which assumed the anchor covered all
strata — it does not). The regime label is still the belief-free total-density landscape; the *response* is the
intergenic anchor (depleted) + capture-gated total density (enriched). Both remain belief-free and fixed across
passes.
