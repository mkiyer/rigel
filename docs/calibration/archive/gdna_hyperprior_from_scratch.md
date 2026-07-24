# Building the gDNA Hyperprior from scratch — clean-slate design (Role B only)

**Purpose.** A ground-up rebuild of the **gDNA composition hyperprior** (Role B), starting from the simplest
possible model and adding features only as the poles demand. Organized around one principle — **gravity** —
and answering the eight open questions + the third-party reviewer's four fixes.

**Hard scope boundary — two NPMLEs, kept separate.**
- **Role A — transfer-variance / enrichment NPMLE** (`DensityNPMLE.project` → `mu_proj`, `var_proj`; fit on
  every node's **total** density before the pass). **LEAVE IT ALONE.** This document does not touch it.
- **Role B — the gDNA composition hyperprior** (`DensityNPMLE.logprior`; fit AFTER pass-0 on the deconvolved
  gDNA). **This is the only thing we rebuild.** Role B may *read* Role A's `mu_proj` (enrichment) as an input,
  but the two fits never merge (reviewer Fix 3).

**Starting point.** The pass-0 solved node belief state: per node `(f_g, Var(log f_g), M, E)`. That is the raw
material; the hyperprior is fit on it and then re-consulted in the refit.

---

## 1. The one principle: GRAVITY

> A node should be pulled toward the **nearest** mode of the landscape, and the pull is **stronger the closer
> it is**. A far mode — however tall — must exert only weak pull. And the pull must be able to go **UP**, not
> only down.

Everything below is in service of this. The current system violates it in three concrete ways, each proven:

1. **The projection is a global search, not a gravitational one.** The δ-pin slices the landscape along the
   ray `ρ_g = f_g·M/E` and the read-out picks the *dominant* mode on that ray. A tall depleted mode at
   log10 ρ_g ≈ −2 out-votes a short enriched mode at +1 **by their height ratio (~50:1)**, regardless of which
   one the node is actually near. That is anti-gravity: the far tall mode wins. (Reviewer point 4.)
2. **The up-pull was deleted.** With a fitted prior, `_gdna_arm` **replaces** the `½·log f_g` Jeffreys
   reference (confirmed: `simplex_logodds.py:161-163`). That reference is an **up-pull** — it rises with `f_g`
   (−∞ at f_g→0, 0 at f_g→1) and is the only thing holding an unstranded node off the f_g→0 vertex. The RNA
   side keeps its `½·log(1−f_g)` (a down-pull). Net: the reference is **asymmetric — down-pull only**.
   (Reviewer point 3.)
3. **The floor doesn't cover where the node is.** The interior ε-floor is bounded to the fitted support
   `[q0.5, q99.5]` (confirmed: `npmle.py:287-296`); the enriched tail where node 1055 sits is **outside** the
   support and keeps the full Gaussian penalty (−7 nats). So the one node that needs gravity gets none.

**Evidence** (all reproducible, no full pipeline): `scripts/debug/crush_dissection.py` (swap the prior shape →
the solver is correct: an enriched-mode prior lands node 1055 at 0.863), `scripts/debug/proximity_projection_test.py`
(global δ-pin 0.0001 vs proximity 0.25 on the same prior), `docs/figures/pass0_npmle_vs_oracle.png` (the fit
*has* an enriched shoulder — shifted to the under-called density — but the global search ignores it; stranded
nodes fit the oracle exactly). Full write-up: `gdna_crush_dissection_node1055.md`.

---

## 2. The two projection philosophies (the central decision)

| | **GLOBAL δ-pin** (current + v0.7.1 KDE) | **GRAVITY / proximity** |
|---|---|---|
| what it asks | "over all f_g, where is the landscape tallest along the ray?" | "which mode did *this node's* density come from?" |
| respects distance? | **No** — far tall mode wins (50:1) | **Yes** — near mode dominates |
| anchor | none (unstranded ⇒ posterior = the raw landscape) | the node's observed density (+ precision) |
| node 1055 | 0.0001 (crush) | 0.25 (honest, ≈ pass-0) |

**Answering Q6 — "if global worked with the KDE, it should work here, right?"** Global did **not** fully work
with the KDE. It crushed the *same* enriched-unstranded nodes — the `ss0.50 + capture_on` collapse was never
solved by any release. It *appeared* to work because global is correct for the two easy majorities: **stranded**
nodes (the strand likelihood anchors the search — gravity comes for free) and **depleted** nodes (the global
mode *is* the right answer). It fails only for the **unstranded-enriched minority**, which has no anchor. The
KDE softened the failure with two crutches the current code dropped — an **up-pull** (`−log(1−f_g)`, its RNA
scale-free Jeffreys) and the **mixture bridge** — but even with them it did not reach the truth at node 1055.

**Conclusion.** Global projection is an acceptable *starting frame* only if we give it back gravity through the
reference terms; true per-node gravity needs enrichment-conditioning. We do **not** need to rewrite the
projection into a responsibility integral on day one — the reference floor (§4) buys most of the gravity inside
the global frame. That is the simplify-first path.

---

## 3. The build ladder (simplest first — do NOT skip ahead)

Each stage is shippable and independently testable against the two poles (§6). Stop as soon as the poles pass.

**STAGE 0 — the minimal viable hyperprior + restored gravity floor.**
- Additive KDE (Q1), unit mass per node (Q2), v0.7.1 substrate (Q3), fixed bandwidth (Q5), **global** projection (Q6).
- **The one fix that matters here:** put the Jeffreys reference **back as a smooth floor under the fitted
  prior** (Q4 / reviewer Fix 2), covering the *whole* f_g range (fixes Q7). This restores the up-pull and stops
  the below-pass-0 crush. *This is expected to be ~80% of the win for near-zero added complexity.*
- Success = node 1055 no longer collapses below its pass-0 0.375; `none_*` unchanged.

**STAGE 1 — counting-precision bandwidth (reviewer Fix 4).**
- Convolve the landscape along the log-rate axis by each node's counting width (`~1/√M` ⊕ `Var(log f_g)`).
  A tall sharp depleted peak flattens far more than a broad enriched hump → the 50:1 ratio shrinks → the
  enriched mode survives the comparison. This is gravity-via-smoothing.
- Success = the enriched mode is no longer out-voted purely on height.

**STAGE 2 — enrichment-conditioning (reviewer Fix 1) = true per-node gravity.**
- Re-parameterize the prior over the **intrinsic** rate `d_g = ρ_g / e(x)`, with `e(x)` the node's enrichment
  from Role A's `mu_proj` (READ only — fits stay separate). At an enriched node `e(x)` re-centers the prior
  expectation **up**, placing mass where the node actually sits — without trusting pass-0's under-called f_g.
- Success = node 1055 pulled *up* toward the truth, not merely held at pass-0.

Proximity-as-a-responsibility-integral (the literal "which mode did it come from") is the same idea as Stage 2;
enrichment-conditioning is the **count-zero-information-safe** way to implement it (the anchor is total density
— the enrichment axis — never a gDNA/RNA composition vote).

---

## 4. The eight questions — recommendations

**Q1 — Additive model?** **Yes, adopt (settled).** No EM competition ⇒ a minority (capture-enriched) mode is
never competed to zero weight. Each training node deposits its own kernel. Retire the competitive-EM path for
Role B (keep it only if Role A needs it — separate object).

**Q2 — Unit density?** **Yes, start here (simplest).** Every training node contributes exactly one unit of
mass. Precision, if used at all, is a **width** (Stage 1), never a weight (weighting down-weights the whole
low-precision *floor* — the exact opposite of what we want).

**Q3 — Training subset? (open)** **Align to v0.7.1:** single-strand + structural REGION nodes, **including
zero-count intergenic/intronic** (they *are* the depleted anchor — no separate floor/tower, no double-count),
**AMBIG excluded** (the two-root ambiguity it must resolve — circularity), boundaries per the release. Keep
Role A / Role B training strictly separate (Fix 3). Revisit stratification-by-signature only if a pole
demands it — not on day one.

**Q4 — Jeffreys pull / gravity / the sign flip?** **The key Stage-0 fix.** Two distinct terms, do not conflate:
- The **gDNA reference `½·log f_g`** (an up-pull / f_g→0 barrier) was **deleted** when the fitted prior
  replaced it → asymmetric down-only reference → the crush. **Restore it as a smooth floor** (below).
- The **RNA scale-free `−log(1−f_g)`** (v0.7.1's extra up-pull) was also removed, deliberately, as an improper
  grid-dependent ramp. **Do not restore that one** — the symmetric reference floor is the principled version.

**Q5 — Bandwidth? Adaptive?** **Start fixed** (`npmle_bandwidth`; the release used Silverman/LSCV floored at
the weighted-median noise). **Adaptive = Stage 1** (Fix 4: widen by counting precision). Don't couple the
day-one build to a bandwidth estimator.

**Q6 — Proximity or global?** See §2. **Global for Stage 0**, but only *with* the reference floor restored
(which supplies crude up-gravity). Global alone will never solve the enriched-unstranded pole — it didn't for
the KDE either. **True gravity = Stage 2** (enrichment-conditioning).

**Q7 — Is the ε-floor big enough / added correctly?** The problem is **coverage, not magnitude**. The floor is
bounded to the fitted support `[q0.5, q99.5]`; the enriched tail (node 1055's location) is *outside* it and
keeps the full penalty. So it never helps the node that needs it. **Replace/augment it with the Jeffreys
reference floor over the whole f_g range** (Q4) — a floor at the *uninformative* level, everywhere. That is
the correctly-designed version of the same idea.

**Q8 — Must pull UP, not only down.** This is the acceptance criterion. Stage 0's floor gives crude up-gravity
(the prior's unreliable tails can no longer penalize *below* uninformative, so high-f_g becomes acceptable);
Stage 2's enrichment-conditioning gives *targeted* up-gravity (re-centers per node). Both are on the ladder.

### The Stage-0 floor, concretely (and smoothly)

The reviewer's `np.maximum(global_logprior, ref)` is correct in spirit but is a **hard cliff** — against our
"no clamps/cliffs" rule. Use a **smooth** floor so the fitted prior sits *on top of* the uninformative
reference and can never penalize below it:

```python
def _gdna_arm(lam, global_logprior):
    ref = _JEFFREYS_REF * _log_fg(lam)[None, :]          # ½·log f_g — the up-pull / f_g→0 barrier
    if global_logprior is None:
        return ref
    # smooth floor: the fitted prior ADDS mass above the reference but never subtracts below it.
    # (soft-max via logaddexp; or an ε-mixture P=(1−ε)·e^{logP}+ε·e^{ref}. NOT np.maximum — that is a cliff.)
    return np.logaddexp(global_logprior, ref)
```

Effect: where the fitted prior has data (the depleted bulk) it dominates the reference and is unchanged; where
it is a deep extrapolated tail (the enriched region) the reference floor lifts it to *uninformative* — so an
enriched node is no longer actively pushed to 0. Apply identically in **both solvers**: the 1-DOF
`_solve_nodes_logodds` and the 2-DOF `_solve_ambig_logodds` (both consume `global_logprior` the same way).

---

## 5. The reviewer's four fixes — assessment

| fix | verdict | stage | note |
|---|---|---|---|
| **1. Enrichment-condition** `d_g = ρ_g/e(x)` | **Agree** — the real gravity | 2 | read `mu_proj` from Role A; fits stay separate |
| **2. Retain Jeffreys floor** | **Agree** — the simplest, biggest win | 0 | but **smooth**, not `np.maximum` (cliff) |
| **3. Separate Role A/B, stratify/additive** | **Agree** — already the intent | 0 | enforce separation; additive + v0.7.1 substrate |
| **4. Counting-precision bandwidth** | **Agree** — gravity-via-smoothing | 1 | flattens the 50:1 height ratio |

All four are compatible and complementary; the ladder sequences them simplest-first.

---

## 6. Acceptance test — the two poles + the invariant

Run on the 32 cached `ambig_dense_10mb` scenarios (`OMP_NUM_THREADS=1`, no full pipeline); refit-vs-oracle by
node type.

- **Pole 1 — crush** (node 1055, `gdna100 ss0.50 capON`, true f_g 0.902): must **not** fall below pass-0
  (Stage 0); should be pulled **up** toward truth (Stage 2).
- **Pole 2 — zero-gDNA** (`none_*`, true f_g 0): must **not** over-call. The floor must lift the enriched tail
  to *uninformative*, not to *enriched* — verify the floor level doesn't manufacture false positives.
- **Invariant — stranded** (`ss0.99`): already matches the oracle (see the figure). **Do not regress it.**
- **Count-zero-information**: composition must never read the raw count. Stage 2's `e(x)` uses *total density*
  as the **enrichment** axis (allowed — it is composition-neutral), never as a gDNA vote.

---

## 7. What to STRIP (the deconstruction)

The current object accreted; simplify it down to:
- **One** density estimator for Role B — the **additive KDE** (retire the competitive-EM Role-B path).
- **One** floor — the **Jeffreys reference floor over the full range** (retire the support-bounded ε-floor and
  the `n_regions` aggregate-cell / background-tower machinery; the structural training nodes are the floor).
- **One** projection frame — global + reference floor now (Stage 0), enrichment-conditioned later (Stage 2).
- Keep Role A (`project`/`mu_proj`/`var_proj`) entirely untouched and clearly separated.

Target: the Role-B path readable in one sitting — fit (additive KDE on the v0.7.1 substrate) → arm (fitted
prior `logaddexp` the Jeffreys reference) → read-out (posterior median, both solvers).

---

## 8. Code touch-points

- **Fit:** `calibration/calibrate.py:_fit_gdna_hyperprior` + `calibration/npmle.py:DensityNPMLE.fit` (additive
  path) — substrate (Q3), unit mass (Q2), bandwidth (Q5), retire the support-bounded ε-floor (Q7).
- **Arm / floor:** `calibration/simplex_logodds.py:_gdna_arm` (:155) — the smooth Jeffreys floor (Q4/Fix 2).
- **Read-out:** `_solve_nodes_logodds` (:357) and `_solve_ambig_logodds` (:500) — **both** solvers.
- **Build:** `calibration/bp_solver.py:327` — where `logprior` is projected onto the grid.
- **Enrichment (Stage 2):** read `mu_proj` from the Role-A `DensityNPMLE.project` already computed in the pass.
- **DO NOT TOUCH:** the Role-A enrichment NPMLE / `project` / σ²_transfer path.

---

## Appendix — evidence index

| artifact | shows |
|---|---|
| `scripts/debug/crush_dissection.py` | swap the prior shape → enriched-mode prior lands node 1055 at 0.863 (no solver bug) |
| `scripts/debug/proximity_projection_test.py` | global δ-pin 0.0001 vs proximity 0.25, same prior |
| `docs/figures/pass0_npmle_vs_oracle.png` | the fit *has* an enriched shoulder (shifted to the under-call); stranded fits the oracle exactly |
| `docs/calibration/gdna_crush_dissection_node1055.md` | the 5-hop code path + the full dissection |
| `docs/calibration/gdna_hyperprior_clean_slate.md` | the earlier strategic study (KDE vs NPMLE, precision weight-vs-width) |
