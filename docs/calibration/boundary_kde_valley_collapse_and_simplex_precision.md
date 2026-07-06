# Boundary KDE-valley collapse & interdependent-simplex precision

> **UPDATE 2026-07-06 — implementation & validation (read this first; it supersedes the framing below).**
>
> **Fix 1 (mixture bridge) is IMPLEMENTED, VALIDATED, and recommended to SHIP.** **Fix 2
> (joint/interdependent-simplex precision) is IMPLEMENTED but EMPIRICALLY NET-NEGATIVE — do NOT ship it as
> implemented.** Both are env-gated (default off ⇒ baseline is bit-exact).
>
> **Corrected framing (the KDE input).** The KDE scores the node's **current-belief gDNA density**
> ρ_g = f_g·d (the candidate the solve evaluates, refined by strand + spliced mass + boundary crossing +
> messages), **never** total density d (RNA's >10⁴× range makes total density a useless gDNA observation).
> The valley collapse was the KDE refusing to let a mixture-density belief be interpreted as a mixture. The
> earlier "total density isn't diagnostic → identifiability wall" note (in §2.2/§2.3 below) was wrong and
> needlessly pessimistic.
>
> **Fix 1 = the mixture bridge** (`bp_solver._kde_logprior`, env `RIGEL_MIX_BRIDGE=ε`). The KDE, estimated
> from clean (unimodal) *region* nodes, has a deep valley between the depleted and enriched modes; a node
> whose gDNA density is a spatial MIXTURE (boundary, sparse-probe region — this is NOT boundary-specific; the
> collapse is a whole density BAND d∈~[0.1,8]) lands in that valley by construction and collapses to f_g≈0,
> emitting a pathologic RNA message. Mixing the KDE with a uniform "any-mixture" bridge over the observed
> gDNA-density support floors the valley (no collapse) while leaving the KDE's real tails OUTSIDE the support
> intact (FP suppression unchanged). **The level is not a finicky magic number**: the peak/valley gap is ~138
> nats (~60 orders of magnitude), so any small ε defeats the cliff; **ε=0.01** is the operating point (a
> peak-relative "the mixture floor is a few e-folds below the tallest mode"; the "1st-percentile-of-heights"
> idea FAILS because the empty inter-mode gap dominates the support by area, so height percentiles sit *in*
> the valley). Continuity-from-flanks (the §2.3 framing) is NOT needed for the prior — messages carry that
> spatially and enter the *solve*, not the prior.
>
> **Fix 1 validation** (in-process, `total_threads=1`, deterministic; the net-flow hard-label benchmark is
> insensitive to soft prior shifts — `benchmark_ab_methodology_cautions` — so the per-node calibration Σ|err|
> is the sensitive instrument):
> * Region 1085 (the #1 error node): f_g 0.624 → **0.927** (true 0.920); boundary 2170 collapse
>   0.003 → 0.993; it leaves the top-error set entirely — **Fix 1 alone resolves the 1085 class** (the
>   trojan-horse message is never emitted once the source boundary is correct).
> * Calibration Σ|err| across all 16 conditions: capture-on failure mode **−11.0%** (ss0.99/none) and
>   **−5.3%** (ss0.99/rnd); **inert** on capture-off and clean gDNA-free libraries (no manufactured FP);
>   small (<2%) regressions confined to unstranded (ss0.50, +1.0%) and no-gDNA-nascent capture-on cases.
>
> **Fix 2 = joint (conditional-variance) message disagreement** (`bp_solver` + `simplex_logodds`, env
> `RIGEL_JOINT_DISAGREE=1`; the log-fraction cross-covariances from the grid → the conditional variance as
> the disagreement anchor, restricted to AMBIG destinations). **Empirically net-negative** (4-arm A/B on the
> target: baseline Σ|err| 34,875; Fix 1 31,020; **Fix 2-only 38,569 (+10.6% WORSE)**; both 32,906 > Fix 1).
> **Why**: for a strand-null AMBIG node the neighbour messages ARE the load-bearing signal; Fix 2 makes a
> confident AMBIG node skeptical of *all* disagreeing messages, but most disagreeing messages are the
> *correct* cross-node continuity correction, not trojan horses — and you can't tell good from bad at the
> destination. The trojan-horse defense therefore belongs at the SOURCE (don't emit the bad message = Fix 1,
> which already removes the collapsed-boundary source), not at the destination. **Fix 1 subsumes Fix 2's
> goal.** Recommendation: ship Fix 1; keep Fix 2 env-gated-off or remove it (the covariance machinery is kept
> only as the record of this negative result).

---

**Status:** diagnosis + design theory (2026-07-05); implemented + validated 2026-07-06 (see the UPDATE above).
**Scenario:** `gdna300 / ss0.99 / capture-on` (the stranded, capture-enriched, moderate-gDNA cell).
**Instrumentation:** `scripts/debug/pass_trace.py` (pass-by-pass per-node trace → self-contained
`pass_trace.npz`) + an adversarial-verification workflow (C1/C2/C5 CONFIRMED, C3/C4 PARTIAL). Every
number below was recomputed in-NPZ (no BAM re-scan); pin `total_threads=1` to reproduce (the ~2.6%
cross-process scanner nondeterminism, `calibrate_cross_process_nondeterminism.md`, does not affect the
in-NPZ decomposition).

Related: `CALIBRATION_ARCHITECTURE.md` (count-zero-info), `precision_state_design.md` +
`node_state_representation.md` (the current 3-variance belief), `dispersion_aware_message_precision.md`
(the per-edge σ²_edge that fails here), `KDE_boundary_prior_review.md` +
`ambig_boundary_spliced_deconvolution.md` (the KDE + the AMBIG boundary), `PHASE2_gdna_mixture_fit_design.md`
(the KDE fit). **This document supersedes the "pass-2 KDE self-reinforces the under-call" claim in
`exon_gdna_peak_smoothing_diagnosis.md` — that mechanism is refuted below (the KDE HELPS; the pass-2
message hurts).**

---

## 0. TL;DR

The worst calibration error node (region 1085, an AMBIG exon, true gDNA fraction 0.92, solved 0.62) is
produced by a **two-pronged** failure, and each prong needs its own fix:

1. **A neighbouring boundary node collapses to spurious ~0 gDNA** because the region-population gDNA-density
   KDE is applied to a node whose density is, *by construction*, a mixture of an enriched exon and a
   depleted intron — landing in the KDE's between-modes valley. This is the **error source**.
2. **The collapsed boundary's confident RNA message crushes the exon**, and the exon cannot defend itself
   because node precision is stored as **three independent per-component variances**. A message that
   nominally targets the exon's *uncertain* RNA− component moves its *confident* gDNA component through the
   simplex constraint — a trojan horse. This is the **propagation mechanism**.

The fixes are, respectively: (1) **gDNA is continuity-pinned, RNA is the residual** — give boundary/mixture
nodes a gDNA prior from genomic continuity of the flanks, not the region-population KDE; (2) **represent the
node belief as a joint distribution on the 2-simplex** (a precision *matrix* in log-ratio coordinates), so
disagreement and message combination respect that the three components are interdependent slices of one pie.

---

## 1. The empirical case

### 1.1 Region 1085 traced end-to-end (all stages CONFIRMED)

Region 1085 is an AMBIG exon (`free_pos & free_neg`; balanced strand, sense 0.535), contained unspliced
mass `raw = 5536`, true gDNA `5094` ⇒ true `f_g = 0.920`. The per-node gDNA fraction at each solver stage:

| stage | f_g | what changed |
|---|---|---|
| `init` | 1.000 | signature-binary AMBIG init `{0,0,1}` |
| `strand` | 0.766 | strand-only local belief — pulls the init **down** (strand is null for AMBIG) |
| `+global floor` | 0.700 | weak floor prior, capped at 1 pseudo-obs |
| `+messages (pass 1)` | 0.700 | inert on this node |
| `+KDE (pass 2 local)` | **0.901** | the trained KDE lifts the local belief **UP toward truth 0.920** |
| `+messages (pass 2)` | **0.624** | pass-2 messages **crush** it back down — the −1637 error |

The single dominant error stage is the **pass-2 message**, and it *overrides an otherwise-correct KDE lift*.
This corrects the prior recorded belief that the KDE reinforces the under-call: the 0.70→0.62 is the
pass1-final → pass2-final *net*; its cause is the message, not the KDE.

### 1.2 The two error classes (top-60 error nodes)

All 60 top-|error| nodes are exons (45 single-strand, 15 AMBIG), in two mechanistically distinct classes:

* **Class A — strand-envelope residual (43/60):** single-strand exons whose strand solve already lands
  within 0.01–0.04 of truth and which nothing downstream moves. They rank high purely by **mass
  amplification** (raw 14k–70k), are sign-balanced, and are a *strand-solve precision ceiling* — **not** the
  1085 pathology. Do not conflate.
* **Class B — AMBIG high-gDNA crush (the 1085 class):** of 15 AMBIG nodes, 12 have true `f_g > 0.85`; the KDE
  correctly lifts them (12/15) and a pass-2 neighbour message crushes them (9/15). Seven share 1085's exact
  profile (AMBIG · KDE-up · msg-down): **1085, 238, 229, 525, 1029, 926, 242**.

**Scope of the boundary collapse:** 39–40 of 132 AMBIG boundaries flip `b1_fg > 0.5 → b2_fg < 0.1`, every one
adjacent to an AMBIG exon; **55 of 146 AMBIG exons (38%)** neighbour ≥1 collapsed boundary. Region 1085 is the
extreme member of a broad class, not a one-off.

**Pool vs per-node:** strand under-calls the AMBIG-init exons by −47.6k net, the KDE corrects +75.4k, so the
**pool net is only −1.9k** while the per-node **|residual| is 34.3k**. Any net-flow benchmark reports this
run as nearly perfect and is blind to the crush — **fixes must be scored per-node.**

---

## 2. Prong 1 — the mixture-node KDE-valley collapse

### 2.1 What happens

The Phase-2 gDNA prior is a KDE of `log ρ_g` (gDNA density) trained on solved **region** nodes. It is bimodal:
a tall **depleted** mode at `ρ_g ≈ 0.012` (log −4.39, from 462 intergenic+intron training nodes) and a broad
**enriched** mode at `ρ_g ≈ 16.1` (log +2.78, from captured exon nodes), with a deep valley between. Per-node
the solve maximises `KDE(log ρ_g) + Jeffreys(ρ_rna)` where `ρ_g = f_g·d`, `ρ_rna = (1−f_g)·d`, `d = M/E` the
node's total unspliced density, and the RNA Jeffreys term is `−log(1−f_g)`.

Boundary node 2170 (the left flank of region 1085) has `d = 5.17`. Sweeping f_g:

```
 f_g     ρ_g     log ρ_g   KDE      Jeffreys   sum
 0.003   0.016   −4.17     −0.26    +0.00     −0.26   ← the only peak (depleted mode)
 0.20    1.03    +0.03    −101.6    +0.22    −101.4
 0.50    2.58    +0.95     −45.8    +0.69     −45.1
 0.80    4.13    +1.42     −25.8    +1.61     −24.2
 0.99    5.12    +1.63     −18.7    +4.61     −14.1
```

The objective's global maximum is at `f_g ≈ 0.003` — the boundary collapses to **near-zero gDNA**, splitting
its mass ~50/50 across the two RNA strands (`fp 0.501 / fn 0.496`). This is grossly wrong: both flanks are
almost pure gDNA (intron 1084 true density 0.02, exon 1085 true density 17.4), so the boundary's crossing
fragments are overwhelmingly gDNA.

### 2.2 Why the likelihood does *not* point to the enriched mode

The user's intuition is correct in spirit — a **gDNA density of 5.17 is ~485 nats more likely under the
enriched mode than the depleted mode** (6.0 vs 31.7 bandwidths away). If we were asking *"which mode did this
gDNA density come from?"* the enriched mode wins overwhelmingly, and the collapse to the depleted mode is
indefensible. So why doesn't the solver see this?

Two reasons, and neither is a fixable "which mode is more likely" rule on its own:

1. **Conservation caps `ρ_g ≤ d`.** The solver never *observes* `ρ_g = 5.17` and asks which mode generated it;
   it *chooses* f_g, so `ρ_g` ranges only over `[0, d] = [0, 5.17]`. The enriched mode (`ρ_g = 16`) is not a
   candidate — not because of an arbitrary `if (d < enriched) then depleted` cliff, but because a node whose
   *total* density is 5.17 **cannot** have a gDNA density of 16. The only KDE peak inside `[0, 5.17]` is the
   depleted-mode spike at `f_g ≈ 0.003`; everywhere else is the valley. So the argmax is the depleted spike.
2. **The node's true gDNA density is genuinely in the valley.** Even at the correct high f_g, `ρ_g ≈ 0.8·5.17
   = 4.1` (log 1.4) — *between* the two population modes. A boundary is, by definition, a **mixture** of two
   regions; when it straddles an enriched exon (ρ_g ~17) and a depleted intron (ρ_g ~0.02), its crossing gDNA
   density is intermediate. **No clean region lives at intermediate density, so the region-population KDE has
   a valley exactly where boundary/mixture nodes live.** Applying the region prior to a boundary is a category
   error — the correct answer isn't at a mode either.

The **free-RNA escape hatch** seals it: because RNA density is priored by a scale-free Jeffreys
(`p(ρ_rna) ∝ 1/ρ_rna`, chosen because RNA spans >10⁴× — `unstranded_spliced_derived_rho`), the model pays
*no* cost to explain `d = 5.17` as `ρ_g = 0.016` (tall depleted mode) `+ ρ_rna = 5.15` (free). gDNA is treated
as the residual *below* a free RNA; the "is 5.17 a plausible gDNA density?" question is dodged entirely.

### 2.3 The fix: gDNA is continuity-pinned, RNA is the residual

The physically-correct signal the KDE is standing in for is **spatial**: gDNA density varies *smoothly* along
the genome (whole-genome shotgun background × a smooth capture-enrichment profile set by probe design +
accessibility — *not* by expression), whereas RNA density is *structured* by transcripts (exons high, introns
low, sharp boundaries). This is the same asymmetry as `rna_imputation_transcript_structure.md` (gDNA flows
genomically; RNA is structured).

So the fix inverts the current logic **for mixture/boundary nodes**:

> **gDNA density is the floor set by genomic continuity; RNA is what sits on top.** A boundary's `ρ_g` is
> pinned by its flanks' gDNA densities (gDNA is continuous across the seam); `ρ_rna = d − ρ_g` is the
> residual. This replaces the region-population KDE — the wrong prior for a mixture node — with a
> continuity prior.

Concretely (candidates, to prototype and A/B):

* **(a) Continuity prior mean.** Set the boundary's gDNA-density prior mean to a crossing-geometry-weighted
  blend of its flanks' solved `ρ_g` (the enriched flank, solved well at the KDE's enriched mode, pins it up),
  at high precision (gDNA continuity is physically strong). The intron flank contributes little gDNA *mass*,
  so weight by gDNA count, not by density, to avoid the depleted side dragging the blend to ~0.
* **(b) Continuity lower bound.** `f_g` of a boundary ≥ the f_g implied by the min sustainable flank gDNA
  density — a boundary cannot carry *less* gDNA than continuity requires. Cheapest guard against collapse.
* **(c) Don't apply the region KDE to boundary nodes at all** — boundaries get gDNA from continuity + strand
  only; the KDE stays a *region* prior.

This same principle disambiguates the general "high total density = captured gDNA vs expressed RNA" question
that AMBIG nodes cannot resolve from strand: the gDNA part is fixed by the smooth continuity profile, and the
rest is RNA — a pure-RNA expressed exon inherits low gDNA from its depleted flanks (correct), a captured exon
inherits high gDNA from the enrichment profile (correct).

**Concerns / regressions to watch:**
* Boundaries legitimately own one-sided spliced/mature RNA + nascent crossing (`spliced_mature_nascent_message.md`);
  the fix must pin gDNA *without* zeroing that RNA — RNA stays the residual `d − ρ_g`, not forced to 0.
* Ordering: continuity needs the flanks solved first. The enriched-exon flank is solved by the pass-2 KDE, so
  the continuity pin belongs in/after pass 2, and must handle the AMBIG-flank case (a flank that is itself
  strand-blind) without a circular pull.
* Do **not** touch the KDE's *region/exon* behaviour — it is correctly lifting 12/15 AMBIG exons to their
  enriched mode. The bug is only its application to *mixture* densities.
* Guard the no-gDNA false-positive cell: continuity from depleted flanks must yield ~0 gDNA on a genuinely
  gDNA-free library (the flanks are all depleted → the pin is ~0), so this should *not* regress FP; verify.

---

## 3. Prong 2 — interdependent-simplex precision (the trojan horse)

### 3.1 Why three independent variances is wrong

A node has ONE total unspliced density `d`, partitioned into a pie `(f_pos, f_neg, f_g)` with
`f_pos + f_neg + f_g = 1`. **The three slices are interdependent: you cannot change f_g without changing
f_pos and/or f_neg.** The uncertainty of a composition is therefore a 2-D distribution *on the simplex* — it
has a **covariance**, not three independent marginal variances.

The current belief (`precision_state_design.md`, `node_state_representation.md`) stores exactly three
independent variances `(var_pos, var_neg, var_gdna)` and messages carry three independent per-component
precisions. This throws away the off-diagonal coupling — which is precisely the information needed to defend a
confident component.

**Worked example (the user's).** A node with `d = 10`, signature EX+|IN−, strand ss = 0.8. The strand tilt
pins the + component: `ρ_pos ≳ 6` (precise). The remaining ~4 is shared between RNA− and gDNA, which both feed
the − reads and are **anticorrelated and jointly uncertain**: `(ρ_neg, ρ_g) = (0, 4)` and `(4, 0)` are both
consistent. The honest belief is: `Var(f_pos)` small, `Var(f_neg)` large, `Var(f_g)` large, and **`Cov(f_neg,
f_g)` strongly negative**. Three marginal variances record the two large variances but *lose the
anticorrelation* — the one fact that says "these two trade off against each other."

### 3.2 The trojan-horse mechanism (verified on region 1085)

Region 1085 in pass 2 has a **confident gDNA** belief (local `f_g = 0.901`, `var = 0.0167` ⇒ precision 59.8)
but an **uncertain RNA− split** (`f_neg = 0.080`, `var = 3.36` ⇒ precision 0.30). The collapsed boundary sends
an RNA− message: implied `f_neg = 0.291` at precision **22.6**.

The disagreement-aware precision (`dispersion_aware_message_precision.md`) is *supposed* to silence a
surprising message. It computes surprise **per component, on the f_neg axis**:

```
residual  = log(f_neg_msg) − log(f_neg_local) = −1.234 − (−2.527) = 1.293
residual² = 1.673
vn_loc    = 3.362                     ← 1085's own f_neg uncertainty
s²_edge   = max(residual² − (base_var + vn_loc), 0) = 0     ← "not surprising" ⇒ FULL precision (22.6)
```

Because 1085 is *uncertain about f_neg*, a message pushing f_neg to 0.29 is "not surprising" **on that axis**,
so it keeps full precision. But **accepting `f_neg = 0.29` forces `f_g` down from 0.90** (the pie sums to 1),
and `f_g` is *confident* (precision 59.8). The per-component surprise metric **never sees the induced f_g
move** — it measures surprise in the wrong space. The least-precise component (f_neg) is a weak link through
which a message reaches and moves the most-precise component (f_g). The most-precise component should have
*stabilised* the solve; instead the least-precise one drove it.

### 3.3 The fix: a joint distribution on the simplex

Represent the node belief as a Gaussian on the 2-simplex in **log-ratio coordinates** (Aitchison geometry),
with a full **2×2 precision matrix**, and combine messages by **precision-matrix addition** — not scalar
per-component precision.

* **Coordinates.** ALR with gDNA as reference is intuitive here (gDNA is the "defended", genomically-continuous
  component): `y = (log(f_pos/f_g), log(f_neg/f_g)) ∈ ℝ²`, with `f_g = 1/(1 + e^{y₁} + e^{y₂})`. (ILR/CLR are
  the symmetric alternative; note the numerical care as `f_g → 1 ⇒ y → −∞`.)
* **The strand likelihood** constrains the + fraction `s = 0.5·f_g + p·f_pos + q·f_neg` — a **1-D** constraint
  ⇒ a rank-1 precision (tight in the strand-balance direction, ~0 orthogonal). For balanced-strand AMBIG it is
  near-zero: the honest statement that AMBIG has a strand null space.
* **The gDNA prior (KDE/continuity) + Jeffreys** act on the "size toward the gDNA vertex" direction ⇒
  precision along that direction. For 1085 this is the confident (precision-59.8) direction.
* **Messages** map a neighbour's densities to a mean + precision in `y`, combined `Λ_post = Σᵢ Λᵢ`,
  `μ_post = Λ_post⁻¹ Σᵢ Λᵢ μᵢ`.
* **Disagreement in the joint space.** A message's surprise is the **Mahalanobis distance** of its implied
  composition from the node's belief under the *full* local precision:
  `surprise = (y_msg − μ_local)ᵀ Λ_local (y_msg − μ_local)`. A message that would move a high-precision
  direction (drop a confident f_g) is now correctly surprising and down-weighted — even if it looks benign on
  the f_neg axis alone. This is the joint-space generalisation of the current per-edge σ²_edge.

Under this representation the 1085 crush self-corrects: the RNA− message implies a composition with
`f_g ≈ 0.62`, which is many σ from the confident local `f_g = 0.90` along the high-precision gDNA direction ⇒
large Mahalanobis surprise ⇒ the message is heavily down-weighted ⇒ f_g stays near the KDE-correct 0.90.

**Concerns / regressions to watch:**
* This is a solver-representation change (belief object, message construction, the combine, the disagreement
  metric). Scope it behind an env gate and A/B against the current 3-variance path.
* RNA imputation is load-bearing for nascent-in-intron and legitimately RNA-continuous nodes; the joint metric
  must not blunt a message that agrees with the *full* belief — verify the beneficial pass-2 up-moves
  (e.g. regions 1083, 226, where pass-2 correctly raises gDNA) survive.
* The strand rank-1 precision must be built honestly (near-zero for balanced strand) or it will re-introduce a
  false strand confidence on AMBIG.
* Numerical stability near the simplex vertices (`f_g → 1`) in the chosen log-ratio coordinate.

---

## 4. How the two prongs interact

They are cause and propagation of one failure, and each fix helps the other:

* Prong 1 (continuity-pinned boundary gDNA) removes the **source**: the boundary no longer collapses, so the
  crushing RNA− message is never emitted (its `ρ_neg 3.86 → ~0`).
* Prong 2 (joint precision) removes the **vulnerability**: even if some node emits a bad RNA message, a
  confident-gDNA destination defends itself.
* Either fix alone helps 1085; both are needed for robustness (Prong 1 for the 40 collapsing boundaries,
  Prong 2 for any residual bad message and for the general AMBIG-defence). **Prong 1 is higher priority** — it
  is the smaller change and it dries up the error at the source.

---

## 5. Recommendations & validation plan

**Order:** Prong 1 first (boundary continuity prior / lower bound), then Prong 2 (joint-simplex precision) if
residual Class-B error remains.

**Non-targets (confirmed):** the KDE's region/exon behaviour (correct), the strand solve (its AMBIG drift is
structural and the KDE repairs it), the Class-A single-strand residual (a strand-precision ceiling, separate
and lower value).

**Validation:**
1. **Metric = per-node residual**, not net-flow (§1.2). Track `Σ|obs_g − true_g|` over exon nodes and the
   Class-B profile-mates (1085, 238, 229, 525, 1029, 926, 242) specifically; watch the boundary flip count
   (39–40 → ?).
2. **A/B in-process, `total_threads=1`** (nondeterminism, `calibrate_cross_process_nondeterminism.md`).
   Use `scripts/debug/pass_trace.py` to re-trace 1085 pass-by-pass after each change (expect: boundary 2170
   no longer collapses; 1085 stays near 0.90).
3. **No-regression cells:** no-gDNA (false-positive gDNA), ss0.50 (the mechanism depends on strand being null
   — re-run the taxonomy there before committing), capture-off.
4. **Downstream:** trace one Class-B locus through `assemble_priors` → EM to confirm the calibration fix moves
   the *deliverable* (abundance / gDNA-vs-RNA split), not just the calibration intermediate.

**Open questions:**
* Is the boundary collapse better fixed at the prior (continuity) or by making the KDE mixture-aware (a third,
  broad "mixture" component)? Continuity is more principled and physically grounded; the KDE-component route
  keeps a single mechanism but risks a tuned mode.
* The nascent-siphon link (`ambig_boundary_spliced_deconvolution.md`): at ss0.99 the net exon gDNA error is
  ~0, so calibration is not the primary source of the +72k nascent hallucination there (that is downstream EM,
  consistent with prior notes) — but the RNA-message-into-exon bias is the calibration-side seed. Quantify via
  the §5.4 EM trace.

---

## Appendix — verified numbers (region 1085 / boundary 2170, in-NPZ, `total_threads=1`)

| quantity | value |
|---|---|
| 1085 true f_g / final f_g / error | 0.920 / 0.624 / −1637 frags |
| 1085 stage f_g | init 1.000 → strand 0.766 → +floor 0.700 → msg1 0.700 → **+KDE 0.901** → **msg2 0.624** |
| 1085 pass-2 messages | gDNA `f_g 0.205 @ prec 0.616`; RNA+ `0.144 @ 1.77`; **RNA− `0.291 @ 22.6`** |
| 1085 local precisions (pass 2) | f_g 59.8 · f_pos 0.28 · f_neg 0.30 |
| RNA− disagreement (per-axis) | residual² 1.673 < vn_loc 3.362 ⇒ s²_edge = 0 ⇒ full precision (not silenced) |
| boundary 2170 flip | b1_fg 0.766 → b2_fg 0.003 (local `fg_loc 0.003`, so the collapse is local, not messages) |
| boundary 2170 density / flanks | d = 5.17; flanks ρ_g: intron 1084 = 0.02, exon 1085 = 17.4 (mixture) |
| KDE modes | depleted `ρ_g 0.012` (log −4.39) · enriched `ρ_g 16.1` (log +2.78) · bw 0.19 |
| ρ_g = 5.17 vs modes | 6.0 bw from enriched, 31.7 bw from depleted (≈485 nats favouring enriched) |
| scope | 39–40/132 AMBIG boundaries flip; 55/146 AMBIG exons neighbour a collapsed boundary |
| pool vs per-node | strand net −47.6k, KDE net +75.4k, pool net −1.9k, per-node |residual| 34.3k |
