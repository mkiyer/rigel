# Enrichment-Sensitivity Recovery — Detailed Work Log, Discoveries & the Pass-0 Pivot

Definitive record of the DNA-prior enrichment-sensitivity work (the workflow trajectory, the concept discovered,
the elegant current best, and the pivot to the solver). **Companion to the quick-start
`dna_prior_session_resume.md`.** When we return from solver work, resume the enrichment-sensitivity thread from
§5 (the elegant core) and §7 (what's left).

---

## 0. The problem, and the pivot

**The problem.** Pass-0 (the prior-free belief-propagation solve) **systematically UNDER-estimates gDNA density
at enriched, unstranded nodes.** The DNA prior we then fit inherits this: its enriched mode is under-massed, so
the *projection* (observed density → anchor `μ*`) cannot pull enriched nodes up — `enr_recovery ≈ 0`.

**Two ways to respond, and the owner's call.** We *found* a projection fix (an **asymmetric-upward** read-out
that de-biases the under-call — §5), and it works. **But it is a Band-Aid**: it corrects a symptom (the biased
observation) downstream, rather than the cause (the biased pass-0 solve). The owner's framing — **the systematic
under-call is a clue that pass-0's assumptions fail on unstranded data; fixing it at the source would cascade to a
better landscape AND better enrichment sensitivity.** So the plan is: **lock in this projection/landscape work
(this document), then dissect pass-0 to recover the gDNA at the source** (§6).

---

## 1. The measured under-call (hard evidence)

Single-strand REGION nodes with real gDNA (`G>0`), pass-0 `obs=log10(ĝ/E)` vs oracle `tru=log10(G/E)`; `(obs−tru)<0`
= under-call. `single`, NOT `ambig` (the AMBIG nodes are what the prior is meant to predict — excluded).
(`/tmp/undercall_diag.py` — regenerate anytime from the cache.)

| cap / ss | mean under-call (all single-strand) | mean under-call (ENRICHED single-strand) |
|---|---|---|
| OFF / **0.50** (unstranded) | +0.415 (depleted OVER-call) | **−0.135** |
| OFF / 0.99 (stranded) | −0.026 | ~0 |
| ON / **0.50** (unstranded) | +0.097 | **−0.169** |
| ON / 0.99 (stranded) | −0.048 | **−0.012** (recovered) |

**Findings:** (1) the enriched gDNA under-call is a **stranding failure** — stranded nodes recover gDNA (−0.01),
unstranded lose it (−0.13…−0.17, up to −0.30). (2) It is a **compression**, not a uniform bias: unstranded
DEPLETED nodes *over*-call (+0.1…+0.4) while ENRICHED under-call — everything drifts toward the moderate default
(the count-zero-information identifiability). (3) It is **worse with nascent RNA present** (nascent competes for
the unspliced mass) and **at low DNA**. Worst scenarios (both suites agree): `gdna300_ss0.50_capON_nrna`
(n_enr≈217, −0.19), `gdna1_ss0.50` (−0.28…−0.30), `none_ss0.50_capON` (−0.22). → the **dissection targets**.

---

## 2. The endpoint & metric (the projection)

The DNA prior is the **third line of defense** (after strand, message-propagation): a gentle anchor, not a
solver. Its endpoint is the **projection**: a node's observed DNA density `d` → anchor `μ*` onto the fitted
landscape (a *sampling likelihood*, "what DNA level was this drawn from?"). Metric (`L.enriched_sensitivity` /
`scripts/debug/gdna_projection.py::score_enriched`): `enr_recovery = mean(μ*−obs)` over truly-enriched nodes
(want **positive**), `enr_abs_err = mean|μ*−oracle|`, `fabrication` = μ* on zero-DNA capOFF exons (specificity
canary). **Projection mechanism validated**: projecting the truth onto the *oracle* landscape returns it to
±0.06; projecting the under-called observation onto the oracle landscape *pulls it up* — so a real enriched mode
is what the projection needs, and the projection itself is sound.

---

## 3. The three workflows — what they did and discovered

All ran as **pure numpy on a cached substrate** (pass-0 once, pickled; no calibrate re-runs) via
`scripts/debug/gdna_explore_lib.py`, over `ambig_dense_10mb` (32) + `quick_3to1_5mb` (16). Metric = EMD-to-oracle
(Goal 1) then the projection (Goal 2).

### 3.1 Landscape exploration — `wf_5df36e4c` (Goal 1)
5 parallel families + synthesis. **Discoveries:** additive **Poisson** (zero-native, `e^{−ρE}`) ≫ point-estimate
KDE (which floors at `1/E` and can't represent sub-resolution DNA) — everywhere. **Zero-count structural** nodes
(intergenic/intron with mass≈0; `Poisson(0|ρE)` = the decay-toward-zero) are the **zero-DNA representation**; the
`live=mass>0` filter was silently dropping them (60–70% of the total error). A **confidence filter** (train on
low-variance nodes) removes the over-call pollution. Winner (synthesis): reliability-mass weighting + boundaries +
zero-native Poisson, **mean EMD-to-oracle 0.267** cross-suite (from a 0.87 baseline). **Boundaries help** (halve
worst-case EMD; they self-solve via spliced+intron so they distinguish real enrichment by *evidence*).

### 3.2 Landscape simplification — `wf_f5ebad7b`
4 families (ablate / unify / minimize-constants / robustness) + synthesis. **Discoveries:** the variance ceiling
that "confidence" cuts against is the **count-zero-info max-entropy variance ~2.83 nat²**, a solver *constant*.
The struct_zero anchor and the boundary info-gate both **EMERGE from one weight** by tempering the belief variance
by the count's add-one share: `w = S0/(v·g/(g+1) + S0)` (`g=0` → `w=1` anchor emerges; high-`g`+high-`v` →
self-suppresses). The **adaptive `tau_prior` gate was proven INERT** — a fixed bandwidth reproduces it byte-
identically (over-engineering, removed). Result: the **unified recipe** — 1 uniform rule, **0 special cases, 1
constant** (`S0=(0.15·ln10)²`), EMD still 0.267, generalizes (identical on both suites — no overfit).
(`scripts/debug/gdna_landscape_recipe.py`.)

### 3.3 Enriched-sensitivity — `wf_12e0fdd0` (Goal 2)
3 families (read-out / enriched-substrate / weighting-bandwidth) + synthesis, against the projection metric.
**The unifying discovery (§4).** Each family closed the ~+0.22 `enr_recovery` gap a different way — read-out via
an asymmetric projection (+0.25), substrate via a no-oracle midpoint relocation (+0.20), weighting via
gentler-weight + density-uplift + balloon smoother (+0.10) — and the synthesis combined them (+0.28) but became a
**15-constant monstrosity that even re-added the inert `tau_prior` gate**. **Hand-distillation showed the elegant
core needs only the read-out** (§5).

---

## 4. The KEY CONCEPT — the directional, bracketed under-call

Pass-0 does not under-call *randomly*; it under-calls **directionally and boundedly**. For any node the true gDNA
density is **bracketed**:

```
   observed  log10(ĝ/E)   ≤   TRUE gDNA density   ≤   total  log10(mass/E)
   (deconvolution under-calls)                      (RNA-inflated over-count)
```

Validated: enriched exons obs +0.79, total +1.09, oracle +0.96 (≈ the midpoint); boundaries obs +0.15, total
+0.44, oracle +0.44. This single fact is why every sensitivity lever worked, and it can be applied either in the
**landscape** (relocate carriers to the bracket midpoint) or — more elegantly — in the **projection** (the
observation is a *lower bound*, so look UP). **It is also the fingerprint of the pass-0 problem**: the solve lands
between the two brackets but too close to the lower one at unstranded enriched nodes.

---

## 5. The elegant core (current best — the Band-Aid that works)

Two pieces, **4 constants total, 0 special cases**:
- **Landscape** — `scripts/debug/gdna_landscape_recipe.py`: `w = S0/(v·g/(g+1) + S0)`, zero-native additive
  Poisson. Robust (EMD 0.267) but enriched-blind on its own.
- **Projection** — `scripts/debug/gdna_projection.py::project_asym(P, d_obs, hup=0.70, hdn=0.02, cap=1.0)`: the
  observation is a lower bracket → trust landscape mass ABOVE (wide `hup`), refuse BELOW (`hdn→0`), bounded to
  `cap` decades; mean read-out.

**Result** (`python scripts/debug/gdna_projection.py`): `enr_recovery` −0.05 → **+0.25** (ambig +0.27 / quick
+0.23), `enr_abs_err` 0.226 → ~0.20 — **matching the 15-constant synthesis at 3 projection constants, no landscape
change.** `hup` larger recovers more but overshoots (`0.9→+0.31/err0.24`; `1.2→+0.33/0.26`); `0.70` is the sweet
spot where recovery AND error both beat baseline.

**Why it's a Band-Aid:** `hup`/`cap` encode "how much pass-0 under-called" — a downstream correction of an
upstream bias. `hdn→0` is principled (obs is a floor); `hup`/`cap` are the open magic (§7).

---

## 6. The pivot — recover gDNA in pass-0 (the source)

**Thesis (owner):** the systematic enriched under-call at unstranded single-strand nodes means pass-0's
assumptions fail there; recovering the gDNA at the source cascades to a better landscape and better sensitivity —
and may make the projection Band-Aid unnecessary.

**Dissection plan (next session):** pick a worst scenario (`gdna300_ss0.50_capON_nrna` — many enriched nodes, clear
−0.19 under-call, realistic; or the `gdna1_ss0.50` extremes). Among its **single-strand** region nodes (NOT
ambiguous), find those that under-call gDNA and ask **why** — per node, look at: the observed strand ratio (should
be ~50/50 in an unstranded library → uninformative), the **spliced** evidence (mature RNA — is the RNA level
over-attributed?), the **neighbour messages** (intron/intergenic gDNA background — is it under-propagating under
capture's non-uniform coverage?), and the **nascent** competition (worst when nrna present).

**Root-cause hypotheses to test:**
- **H-a Identifiability drift.** Without strand, the reference/messages pull the node toward the moderate default
  → enriched under-call, depleted over-call (the compression in §1). *Is the default too strong / the enriched
  signal too weak?*
- **H-b Nascent competition.** The unspliced mass = gDNA + nascent RNA; pass-0 over-attributes to nascent →
  under-calls gDNA (worst when nrna present, §1). *Can spliced/junction evidence bound the nascent share?*
- **H-c Message under-propagation under capture.** The gDNA background from neighbouring introns/intergenic should
  raise the enriched node's gDNA, but capture breaks the uniform-coverage assumption the messages rely on
  (`docs/calibration/CALIBRATION_ARCHITECTURE.md`) → the message is too weak. *Is the enrichment-crossing
  transfer-variance over-damping the gDNA message?*
- **H-d Spliced-derived gDNA unused.** The spliced fragments measure mature RNA directly; the residual unspliced
  minus a mature-derived RNA estimate is a confident gDNA floor that isn't being exploited at single-strand exons.

Whatever pass-0 recovers, re-measure §1's under-call and §5's `enr_recovery` — improvement at the source should
lift both.

---

## 7. What's left (both tracks)

**Enrichment-sensitivity track (paused here — resume from §5):**
1. **Derive `hup`/`cap` per-node** — the last projection magic; scale with the node's deconvolution/belief
   uncertainty (more uncertain → under-called more → wider `hup`). Target 0 tunables.
2. **Specificity round** — the fabrication canary drifted up (zero-DNA μ* −0.83 → −0.4…−0.6); fix by
   **evidence-gating** (boundaries/spliced distinguish real enrichment from over-call) — likely buys specificity
   back "for free."
3. **Wire into the solver** — unified landscape + `project_asym` as the gentle Gaussian anchor in ψ, both solvers,
   A/B via `gdna_hyperprior_eval.py`; remove the abandoned δ-pin gates (resume guide §8).
4. **Bigger battery** — cache `gdna_benchmark_5mb`, `cfrna` (REAL), `complex_25mb`; confirm generalization.

**Solver (pass-0) track (the pivot — next up):** §6 — dissect the worst scenario's single-strand under-callers,
test H-a…H-d, recover gDNA at the source, re-measure. **This is the priority per the owner** (fix the cause, not
the symptom).

---

## 8. DISSECTION RESULTS — the pass-0 source diagnosis (2026-07-21)

Ran the REAL production `calibrate()` (refit=0) on `gdna_gdna300_ss_0.50_nrna_present_capture_on`, dissected the
**single-strand region nodes** that under-call gDNA via the `_debug`/`_capture` hook (strand-only → +prior →
final f_g, received message, neighbours). Tool: `scripts/debug/flagship_interrogate.py` (class table) +
`scratchpad/pass0_ss_dissect.py` (the H-a…H-d probes). **The hypotheses resolved cleanly to ONE root cause.**

**H-a — CONFIRMED, this is the root.** Every single-strand node has `fg_strand ∈ [0.49, 0.51]` (100% within 0.05
of 0.5). On unstranded data (κ=½) the strand Fisher info is EXACTLY 0, so the message-free solve pins every exon
at f_g=0.5 — the count-zero-info identifiability floor, by design. The final only reaches mean 0.564 vs oracle
0.909 over enriched nodes (messages recover ~13% of the 0.5→0.909 gap). This is a **compression toward 0.5**:
enriched nodes under-call, depleted over-call (§1).

**H-c — CONFIRMED as the failed rescue, and the mechanism is CAPTURE STARVATION.** Messages are the only escape
from 0.5, but the enriched exons (mass 50–80k fragments) are flanked by boundary nodes carrying **2–17
fragments**. *Capture concentrates gDNA mass CONTAINED inside exons and empties the CROSSINGS.* Since the BP chain
routes all gDNA info through boundaries, the enriched exons are effectively **disconnected**. So the received
message is (a) low-precision (63% gagged, prec<0.05 — the tiny source count `sm`), (b) **wrong-signed** (mode −6
to −7: a mass-2 boundary honestly reports "≈0 gDNA density here," correct for it but wrong for the exon), and (c)
σ²_transfer-damped across the ~3-decade capture cliff. Net: f_g barely lifts off 0.5. Only 15% are true islands
(no oracle-enriched neighbour) — the info is mostly nearby but **cannot propagate** through the empty crossings.

**H-b — consistent** (nascent competes for the unspliced mass; worse with nrna present, §1) but subsumed by H-a:
with no strand, gDNA and nascent are locally indistinguishable regardless.

**H-d — CONFIRMED as an unused lever, but weak.** The motif-stranded spliced reads are a true-strand RNA
measurement, but they measure MATURE, not the nascent that shares the exon's unspliced mass — an indirect bound.
(They live on boundary nodes, not region faces, and never enter ψ by design.)

### The verdict — this IS the "honest ceiling," and the source-fix is the fragment-length channel

The single-strand under-callers are precisely the failure `CALIBRATION_MASTER.md` §9 already names — *"a
fully-unstranded (κ=½) node that is an enriched, spliced-free island with depleted neighbours is at the
identifiability ceiling — strand is silent, messages have no enriched-level neighbour to borrow from, and the
gDNA hyperprior can only supply the background level."* The dissection re-derived it from the pass-0 under-call.

**The one signal that breaks it: the FRAGMENT-LENGTH mixture.** gDNA and RNA fragment lengths are near-disjoint in
this suite (gDNA 111±32 vs RNA 200±4.5; total-variation 0.886; **per-fragment Fisher info for f_g = 3.18** ⇒
se(f_g)=0.0025 at a 50k-fragment exon) — a strand-INDEPENDENT, per-fragment composition signal that calibration's
accumulator collects (the FL pools) but the per-node ψ marginalizes away into eff-lengths and never uses. It is
NOT a count-zero-info violation: the count still enters only as precision; the fragment *lengths* vote, exactly
as the strand split does (§1.1's "only intrinsic signal" is why FL was never wired — but the MASTER §9 and a prior
29-agent adversarial workflow both name FL as **the** resolver for this ridge, "already in the scorer, absent from
calibration"). It is the doc's own named escape — currently gated behind "only-if-we-plateau," which the owner's
pivot (fix the source, not the Band-Aid) reopens.

**CAVEAT — sim-inflation.** This suite generates fragments at delta-function lengths (gDNA=100, RNA=200,
frag_std=0), so the FL separation is a near-oracle toy artifact. Real gDNA/RNA fragment lengths overlap heavily,
and hybrid CAPTURE biases FL (probe hybridization length preference) so the library-average gDNA FL model may not
match the per-node captured FL. The FL channel is principled and buildable, but the toy will over-state its
real-data payoff — exactly the "toy may not extrapolate" risk. **Decision point for the owner (see below), not a
unilateral build.**

**Infra:** `scratchpad/pass0_ss_dissect.py` (the H-a…H-d dissection, re-runnable), `flagship_interrogate.py`
(production-faithful node trace, `--condition`/`--node`). The FL Fisher measurement is a ~10-line check on
`inp["gdna_fl_pmf"]`/`inp["rna_fl_pmf"]`.

### 8b. The full per-node DOSSIER — the RNA competitor is NASCENT (single-exon-like) [owner probe, 2026-07-21]

Full pass-0 dossier (init → strand counts → every directional message with its source → outcome + transcript
structure) on the worst SOLVABLE single-strand under-callers: `scratchpad/pass0_ss_dossier.py` (captures the
production `init_beliefs` + `build_node_statics` via a thin wrapper — no production edits). Findings:

1. **The single worst node (1781) IS a single-exon nascent span** (`G0186.1`, n_exons=1, is_nrna=True; oracle
   f_g=0.947, solved 0.243; mass_spliced=0). Genuine floor — no junction, no strand.
2. **But single-exon is the TAIL, not the bulk.** Enriched single-strand under-callers: **single-exon 11 nodes /
   6% mass, mean under-call −0.47 (worst per-node); multi-exon 178 / 94% mass, −0.21.** The bulk is exons of real
   multi-exon genes (`G0200`/`G0106`/… 3–6 exons).
3. **The competitor is NASCENT everywhere.** Every multi-exon under-caller has `mass_spliced=0` AND overlaps a
   `RIGEL_NRNA_…` single-exon nascent span. The RNA sharing the exon's UNSPLICED mass is nascent pre-mRNA —
   junction-blind (spliced/mature is a disjoint pool at the boundaries) and strand-blind at κ=½. So the owner's
   "single-exon is undetectable" GENERALIZES: the RNA that must be separated from gDNA in the unspliced mass is
   the undetectable (nascent, single-exon-like) kind. This is the identifiability floor restated concretely.
4. **The messages actively HURT.** e.g. node 1909: strand→0.51, but the gDNA message (mode −6.4 = "no gDNA", from
   a mass-5 capture-starved boundary) drags f_g DOWN to 0.394 (oracle 0.984). Starved-boundary messages make the
   under-call WORSE than the bare 0.5 floor.
5. **The one informative neighbour is gagged.** Node 1055's left boundary 1054 (mass 206, oracle gDNA 0.94,
   spliced 96) emits NOTHING — the τ-evidence gate zeroes its precision on unstranded data (I_strand=0 ⇒ τ=0).

**Two concrete levers (from this dossier):**
- **(A) Withhold single-exon transcript regions from hyperprior TRAINING** (owner's fix). Correct + easy — they're
  undetectable and worst per-node (−0.47), so they pollute the landscape most. Addresses the ~6% tail, not the
  94% bulk. Lives on the hyperprior/landscape track.
- **(B) The fragment-length channel** breaks the gDNA-vs-NASCENT degeneracy at the source (nascent = RNA length
  ~200 ≠ gDNA ~111; junction-blind and strand-blind, but FL-separable). This is the bulk fix. See §8.

### 8c. The INTERGENIC↔EXON (TSS/TES) boundary — a pure-gDNA signal we corrupt into the wrong sign [owner lead, 2026-07-21]

Owner's insight: on unstranded data the ONLY signal sources are (1) intergenic regions, (2) intronic regions, (3)
boundaries with spliced mass. A single-exon transcript has two boundaries — TSS + TES — with NO introns ⇒ no
spliced ⇒ no RNA detection. BUT the TSS/TES boundary is **intergenic↔exon**, and NO RNA crosses a transcript
start/end (transcription begins AT the TSS), so its crossing mass is **structurally pure gDNA**. That pure-gDNA
composition should pull the adjacent exon toward gDNA (imprecisely — the crossing is thin + capture-depleted).

**What we currently do there (`scratchpad/tss_boundary_probe.py`, verified):**
- 430 intergenic↔exon (TSS/TES) seams; **oracle f_g = 1.00 EXACTLY** (p10–p90 all 1.00) — pure gDNA confirmed. We
  DO lock them to f_g=1 (G1 init) and they DO emit inward. Contrast: exon↔exon seams oracle 0.74 (RNA-contaminated
  — the legit phantom concern); intron↔exon 0.79.
- **BUT the inward message is the DENSITY mode** `log(seam_gDNA_density / exon_total_density)`. Capture depletes
  the seam crossing (mass ~5) while the exon is enriched (mass ~50k), so the ratio is tiny and the message reads
  as f_g≈0: e.g. seam 2 (oracle 1.00) → exon 3 sends implied f_g=0.007. **38% of the 430 seam messages read as
  f_g<0.5 — a pure-gDNA source telling the exon "no gDNA here."** This ACTIVELY pulls the exon below the 0.5
  strand floor (node 1909: strand 0.51 → 0.39). The composition certainty (f_g=1) is destroyed in density space.
- **Root:** `use_shift = (not exon@src) and (not exon@dst)` drops the composition mode because the DST is an exon,
  gagging even a trustworthy pure-gDNA SOURCE. And `struct_lock` (composition-certain) requires `is_region`, so
  ALL boundary seams are lumped + excluded to avoid the exon↔exon phantom — punishing the safe intergenic↔exon
  case for the unsafe exon↔exon one.

**THE FIX (proposed, surgical):** a structurally pure-gDNA source (intergenic region OR intergenic↔exon seam)
sends a **COMPOSITION** message `f_g→1` at honest LOW precision (count-limited by the thin depleted crossing) —
a gentle inward pull toward gDNA, not the current density message that pulls away. Corrects the sign for all 430
seams; keeps precision low (the phantom cascade concern was HIGH-precision seam emitters in gdna_none — separable
from the mode). **Must re-check the zero-gDNA FP guard.** This is the pass-0 SOURCE fix, using existing structure
(no new signal) — distinct from and complementary to the FL channel (§8). Prototype next.

---

## Appendix — provenance
Workflows: `wf_5df36e4c` (landscape-explore), `wf_f5ebad7b` (landscape-simplify), `wf_12e0fdd0`
(enriched-sensitivity) — scripts under `…/workflows/scripts/gdna-*.js`, full results in the task `.output` files.
Infra: `scripts/debug/{gdna_explore_lib,gdna_cache_build,gdna_landscape_recipe,gdna_projection}.py`; caches in the
scratchpad. Quick-start: `dna_prior_session_resume.md`. Rejected dead-ends: `π_r` RNA prior, enrichment-
conditioning log-shift, Stage-0 logaddexp floor.
