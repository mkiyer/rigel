# PLAN — THE CALIBRATION FIRST PASS, REDESIGNED. Stages 0–2.

    ⚠ **A DEV DOC.** Nothing may cite it and it is not the state. When a stage settles, its RULING
    moves to `DESIGN.md`, its DERIVATION to `EQUATIONS.md`, its LESSON to `TRAPS.md` as a named rule,
    and its current number to `ROADMAP.md` §0 — and is deleted from here in the same edit.

⭐⭐⭐ **THE OWNER'S PLAN, 2026-08-22.** A stepwise redesign of calibration's first pass, presented three
stages at a time; the later stages are being written while these are built. **Stage 0** derives a
SUBSTRATE from structure alone. **Stage 1** fixes the intergenic slots and models the gDNA background
from them. **Stage 2** deconvolves single-stranded introns from two independent channels, message-free.

⛔ **JUDGE ONLY WHAT IS ACTIVELY SOLVED** (owner): while these stages are in flight, every score is
restricted to the slots the stage claims. The rest of the chain is not yet this design's business, and a
pooled number over all slots would hide both the win and the damage.

---

## §A THE ONE THING THIS DESIGN RESTS ON, STATED FIRST

⭐⭐⭐ **THE SUBSTRATE IS A gDNA ANCHOR SET ONLY BECAUSE NASCENT RNA IS SPARSE.** Every slot the plan
calls solvable is one where **mature** RNA is structurally absent — that is what the annotation can
prove. What remains in an unspliced count there is `gDNA + nascent`. So:

* under the NASCENT SCOPE RULING (`DESIGN.md` §0b) nascent is absent by default and the slot is a gDNA
  anchor, which is what makes stages 1–2 work at all;
* where nascent IS present the slot is a mixture, and stage 2's second channel is precisely what has to
  notice. ⛔ So "sparse" is not a convenience here — it is the load-bearing assumption, and the design
  must fail LOUDLY rather than silently when it is violated.

⚠ This is why the sparse-nascent panel rebuild had to come first: on the retired uniform panel every
intron carried nascent at the local mature level, so no intron was ever a gDNA anchor and this substrate
could not have been validated.

---

## §B STAGE 0 — THE SUBSTRATE, DERIVED FROM STRUCTURE

### B.1 What is in it

| class | predicate, in words | why it is solvable |
|---|---|---|
| intergenic REGION | no transcript covers it | no RNA of any kind is admissible: composition is certain |
| intergenic\|exon BOUNDARY | a gene edge | same, one dimension down |
| single-stranded intron REGION | exactly one strand admissible, intronic | no mature RNA is continuous here; unspliced mass is `gDNA + nascent` |
| single-stranded intron\|exon BOUNDARY | as above, at a boundary | same |
| ⭐ **SOLVABLE EXON REGION** | single-stranded, and **at least one flanking boundary is a splice site with no contiguous exon crossing** | see B.2 — the flanking boundary is an anchor and the sj flux is certified RNA |

### B.2 ⭐ WHY THE EXON RULE WORKS — the mechanism, because a substrate without one is a guess

Take an exon REGION `R` whose flanking boundary `b` carries a splice junction and across which **no
mature transcript is continuous**. Three facts hold at once, and each is structural:

1. **`b`'s unspliced crossing population has no mature term.** Mature RNA reaching `b` splices; it does
   not cross contiguously. So `b`'s unspliced count measures `gDNA + nascent` — sparse ⇒ gDNA.
2. **`b`'s spliced flux is certified RNA.** A spliced fragment cannot be gDNA. So the sj count at `b`
   measures mature RNA entering `R`, with no deconvolution required.
3. **`R`'s own mass is `gDNA + nascent + mature`**, and `R` shares its gDNA field with `b` — they are
   adjacent, so the gDNA density is locally the same quantity.

⭐ So `R`'s gDNA is reachable WITHOUT imputing anything from far away: its neighbour supplies a gDNA
level whose mature term is structurally zero, and the sj flux accounts for the mature RNA. ⛔ **This is
not the message layer.** It is one hop along a chain whose endpoint is structurally determined — but that
distinction is exactly the one `TRAPS: one-hop-lifted-out-is-still-the-relay` refuses to accept, so the
design must either use the message framework's own operator for the hop or state why this is a different
thing. **OPEN DECISION `the-one-hop`, in §G.**

⛔ **The AMBIG exclusion is not a preference.** On unstranded data the strand channel is identically
zero, so an AMBIG slot has neither strand nor a structurally-clean neighbour; it is the deferred
stratum's blindness restated. Excluding it from the TRAINING substrate is what keeps the prior honest;
it does not mean those slots go unanswered later.

### B.3 The predicate, in the arrays that already exist

    single_stranded(s)  :=  free_pos[s] XOR free_neg[s]
    solvable_exon(R)    :=  single_stranded(R) ∧ exonic(R)
                            ∧ ∃ b ∈ flanks(R) : is_splice_site(b) ∧ ¬mrna_active(b)

⚠ `¬mrna_active(b)` IS the owner's "no exon|exon boundary component": `mrna_active` means an annotated
mature transcript is continuous across the position. `is_splice_site(b)` distinguishes the splice-site
boundary from a terminus or a gene edge, which also satisfy `¬mrna_active` but carry no sj.

### B.4 Where it is computed — OPEN DECISION `substrate-home`

The owner's plan says index-creation time. Two homes, and the choice is not free:

* **at index time, persisted** — matches the precedent of `mature_walls` / `boundary_reach`, which are
  index-derived and persisted because the C++ accumulator needs them. Costs an index schema field and a
  hash concern; ⚠ `reach` is already covered by no existing hash, which is a recorded trap.
* **at calibration init, derived** — a pure function of statics already in memory, O(n_slots) array
  math, cannot go stale, no schema change, no cache invalidation.

⭐ **RECOMMENDATION: derive at calibration init**, in one named module with one gated definition, and
persist later only if profiling shows it matters. The substrate is consumed only by calibration; nothing
in the scan needs it. This is reversible in a way an index schema change is not.

---

## §C STAGE 1 — INTERGENIC INITIALIZATION

**What it asserts:** intergenic REGIONs and gene-edge BOUNDARIEs are `f_g = 1` at variance 0, fixed.
**What it then measures:** the gDNA background distribution, from intergenic REGIONs ONLY.

⭐ Most of this exists. What the redesign changes is that the LOCK becomes STRUCTURAL rather than
derived from a runtime solvability predicate — which is the standing defect that the emission lock is
scoped `~solvable ∧ REGION` instead of `g1_locked ∧ REGION`.

⛔⛔ **AND THAT IS NOT A FREE WIN — IT IS HALF OF A RECORDED CANCELLING PAIR.** Rescoping the lock alone
was priced and refused (`ROADMAP.md` §4.1; `TRAPS: a-cancelling-defect-pair`): it improves its target and
costs the zero controls, because the mis-scoped mask is load-bearing for an unrelated population. So
stage 1 must be priced WITH whatever replaces that load — not landed on the grounds that it is more
principled.

⚠ Boundaries are excluded from the background fit by the owner's own ruling, and the code already pools
REGIONs only. The reason is worth writing down where the fit lives: a boundary is a crossing count with
no genomic measure, so it estimates the same rate through a different opportunity model.

---

## §D STAGE 2 — INTRON DECONVOLUTION, TWO CHANNELS, MESSAGE-FREE

### D.1 The two channels are genuinely independent, and here is why that matters

The plan calls for an **honest precision-weighted combine** of an abundance/density answer and a strand
answer. That is only legitimate if the two are independent — and `TRAPS: two-gaussians-one-latent` records
what happens when they are not: adding the Fisher information of two channels built from one latent is
exactly 2× over-confident.

⭐ **They are independent here, and it is the standard Poisson factorisation.** A slot's fragments arrive
as a count `N`, and each is assigned to `+`/`−`. Condition on the slot: `N` is Poisson and the strand
split is multinomial given `N`, so **the total and the split are independent statistics**. The density
channel reads `N` against an opportunity; the strand channel reads the split. Their informations add
legitimately.

⛔ **What is NOT independent, and must be checked in the implementation:** both channels are functions of
the same underlying densities, so if either is re-parameterised in terms of the other's estimate the
independence is destroyed. The combine must take each channel's own raw sufficient statistic.

### D.2 The prior

The plan's candidate: `f_g_prior = ρ_bg / ρ_total(intron)`, clamped to 1 when `ρ_bg ≥ ρ_total`.

⭐ It is derivable rather than chosen, which is what makes it admissible under `TRAPS: no-magic-numbers`:
the background density is measured in stage 1 and the intron's total density is a composition-free
observable. ⚠ Two things to settle:

* **the clamp's semantics** — `ρ_bg > ρ_total` is a *measurement collision*, not a certainty. Clamping
  the PRIOR's location to the vertex is not the same as asserting `f_g = 1` with variance 0, and
  `ROADMAP.md` §4.2 records that a near-vertex location is a strong claim in nats whatever the strength
  budget. **OPEN DECISION `clamp-semantics`: is the clamp a location clamp, or does it also cap the strength?**
* **the prior's STRENGTH** — the plan does not state one, and location is strength on the λ scale
  (`strength = logit(m)`). The shipped structural reference carries one pseudo-observation; the same
  budget is the obvious default and needs no new constant.

### D.3 The combine — OPEN DECISION `combine-form`

Two forms, and the existing solver already implements one of them:

* **A — the lattice product (what exists).** Each channel contributes a log-factor on the λ grid and the
  posterior is their product. No Gaussian approximation, precisions fall out as curvature, and a channel
  that is uninformative contributes a flat factor and therefore nothing.
* **B — two point estimates plus precision weights (what the plan's words describe).** Each channel
  produces `(x_i, p_i)` and the combine is `Σ p_i x_i / Σ p_i`. Simpler to audit per channel; loses the
  shape; needs the independence argument of D.1 to be exactly right.

⭐ **RECOMMENDATION: keep the lattice product as the SOLVE (it is the honest combine, and at κ = ½ the
strand factor is exactly flat so the degenerate case is automatic), and PUBLISH each channel's own point
estimate and precision as a diagnostic.** That gives the plan's auditability — you can see what each
channel said and what the combine did with it — without approximating the posterior. If the two
disagree, the diagnostic is where that shows up.

---

## §B.5 ⭐⭐⭐ STAGE 0 STATUS — BUILT AND VALIDATED, 2026-08-26

* **The module**: `rigel.calibration.structural_claims` (layer 3, declared in `_layers.py`), derived at
  calibration init per the `substrate-home` ruling. Four disjoint classes + the licensing-flank arrays.
* **The unit gate**: `tests/calibration/test_structural_claims.py` — a designed six-locus GTF where
  every clause has a positive and a negative case, hand-enumerated membership, and a brute-force
  per-slot reference. Written first and verified failing; then NINE perturbations of the fixed code,
  each caught by its designated gate. ⭐ Two perturbations were SILENT on the first fixture (no AMBIG
  intron without exon shielding; no minus-strand locus, so `~free_pos` masqueraded as `g1_locked`) —
  both holes closed and the sweep re-run 9/9 FIRED.
* **The confusion matrix**: `scripts/design/structural_claims_audit.py` (self-test 8/8, including
  the judge-only-claimed-slots falsification). ⭐⭐⭐ **RESULT: every claim holds on every condition of
  BOTH panels — 32/32, zero violating fragments in any class.** Coverage of the unspliced library on
  claimed slots: 27.6 % (capture-ON) to 94 % (`g98` capture-OFF). Suite baseline moved 3,658 → 3,676,
  re-derived per the table.

## §C.1 ⭐⭐⭐ STAGE 1 STATUS — the assertion already ships; the rescope RE-PRICED AND REFUSED (2026-08-26)

Stage 1's two assertions are production code already (`g1_locked` pins both axes at Var 0;
`fit_intron_background` pools intergenic REGIONs only). The one proposed change — the emission lock
rescoped `~solvable ∧ REGION` → `g1_locked ∧ REGION` — was priced exactly as §C demands: TOGETHER
with the redesign's replacement for the mis-scoped mask's load (the stage-2 measured reference), as
`ladder_arm_ab --arm stage1_pair` and `--arm stage1_pair_onesided`, messages on, 16 ladder conditions.

⛔ **REFUSED, both compositions, against the stage-2-alone baseline**: the pair's marginal effect is
+1.9 % / +1.6 % / −1.9 % on the three in-scope strata, +19.0 % on the deferred one, and its wins are
confined to `g00` (−15.6 % further) — one-sided evidence by ROADMAP rank 2's own rule. Adding the
one-sided certified-RNA bound deepens `g00` to −64.4 % from base and improves unstranded × OFF to
−12.0 %, but flips stranded × ON to +1.4 % and takes the deferred stratum to +29.9 %. ⭐ The pass-0 vs
deliverable inversion at `g00` (pass-0 +141 % while the deliverable improves) is §1.1's
"prior-fidelity anti-correlated with deliverable" shape reproduced on the current panel. The five
pair-gated xfails stay xfail; the rank-11 brief carries the update.

## §D.4 ⭐⭐⭐ STAGE 2 STATUS — PRICED ON BOTH PANELS AND LANDED, default OFF (2026-08-26)

**The scope ruling came from the A/B, not from preference.** `stage2_intron_ref` (regions +
boundaries) vs `stage2_intron_ref_r` (regions only), test chromosome then ladder, messages on:

* **FULL arm — REFUSED.** Regressed stranded × capture-ON on both substrates (+38.7 % test chr,
  +4.2 % ladder with the ladder's boundary axis +31.7 %), and the messages-OFF run attributed the
  damage to the LOCAL solve on the BOUNDARY axis alone (+52.3 % boundaries, +0.1 % regions): under
  capture the intergenic rate is the wrong rate for a probe-adjacent crossing.
* **REGIONS-ONLY — ACCEPTED, every stratum improved on the ladder**: deliverable −11.7 % (stranded ×
  OFF), −7.3 % (unstranded × OFF), **−1.2 % (stranded × ON — the stratum every previous per-object
  reference regressed)**, −3.5 % deferred, **zero-gDNA control −19.6 %**; worst single cell +0.7 %.
  The test chromosome's unstranded × OFF +10.3 % did NOT transfer to the ladder (it improved) —
  `TRAPS: a-toy-and-a-panel-can-disagree-in-rank`, met again, in rank order.

**Landed**: `density_deconv.measured_reference_location` (layer 5), gate
`tests/calibration/test_measured_reference.py` (7 cases written first and verified failing; 7
perturbations each caught), plumbed `calibrate → solve_chain(reference_location=…)`, behind
`CalibrationConfig.measured_intron_reference` — **default OFF, bit-identical** (suite 3,685/9 xfail,
zero goldens moved). ⭐ The flag-ON location array is BYTE-IDENTICAL to the priced arm's on a real
ladder condition (8,640 slots, 4,967 collisions at `g50 ss0.50 ON`). ⚠ The location deliberately uses
the bare pooled ratio, not the factory's Jeffreys posterior — the 0/0-keeps-base branch needs an
exact zero on a zero-gDNA library; recorded in the function docstring. ⛔ The default flip is the
owner's call (`preflight.py --full` after).

**The full §E bar is met (2026-08-26):** DELIVER/REFUTE at claimed intron regions, four `g50`
conditions, never pooled — **DELIVER −21.8 % to −45.0 %, REFUTE within +5.1 %** (the FIRM trade,
quantified); and BOTH zero controls via `zero_controls.py`'s own machinery, flag OFF vs ON — the
zero-gDNA control exact (+0.0001) in both configs, and the flag FIXES two of the three standing
zero-RNA failures (silent −0.0363 → −0.0018, spliced_exons −0.0591 → −0.0029; the third is an EXON
slot outside the intron scope, a pre-existing defect this stage does not claim).

## §E VALIDATION — what has to be true before a stage is called done

Per stage, in this order, and nothing is quoted before it: ① a falsification test written FIRST and
verified failing, then each gate fired by perturbing the fixed code; ② the substrate itself checked
against the certified per-slot truth on both panels — for stage 0 that means "every slot this predicate
admits really does have the property it claims", scored as a CONFUSION MATRIX against truth, not as an
aggregate; ③ per-stratum scoring restricted to the claimed slots, both zero controls included; ④ the
test chromosome first, then the 16-condition ladder (`TRAPS: panel-before-src`).

⭐ **Stage 0's own test is the cheapest and the most decisive**: it needs no solver at all. Take every
slot the predicate admits and read its certified `true_f_g` — a solvable exon whose truth says otherwise
falsifies the predicate immediately.

---

## ⭐⭐⭐ VERIFIED AGAINST THE CODE — the facts that change the design

*Read by a grounding sweep, 2026-08-22 — all five mappers reported. Their full reports are not in the
repo; the load-bearing subset is here.*

**① STAGE 2's STRAND CHANNEL ALREADY PRODUCES A STANDALONE PER-SLOT ESTIMATE, in production.**
`region_geometry._type_belief` keeps a strand-only solve wherever
`g2_active = (free_pos ^ free_neg) & (mass_unspl > 0)` (`region_geometry.py:631-644`) — its docstring
says "the STRAND DECONVOLUTION alone resolves the pie". So the plan's channel-one estimate is not new
code; it is an existing quantity that needs surfacing. Measured accuracy at N = 1000, κ = 0.99:
`f_g` 0.0323 / 0.2312 / 0.7688 / 0.9677 against truth 0 / 0.25 / 0.75 / 1.0 — the residual is the
Beta(½,½) reference's vertex repulsion, not noise.

**② A PRECISION-WEIGHTED COMBINE ALSO ALREADY EXISTS, and its hazard is already written down.**
`tau_lam = I_strand + I_density` (`region_init.py:296-303`) is additive on the precision axis, and the
relay fuses two claims exactly as the plan asks — `c_tau = _fuse_add(atau, btau)` with the λ mode as a
τ-weighted mean `lam_msg = (atau·alam + btau·blam)/c_tau` (`relay.py:1101, 1116`). ⛔ Its comment
(`relay.py:1110-1114`) states the trap the plan must adopt: **each claim must be fused by ITS OWN
precision, never by a ratio read back off the other fuse** — "that mismatch delivered the split at a
confidence it was never weighted by". ⭐ This resolves `combine-form`: the idiom exists, is used, and
carries its own recorded failure mode.

**③ AT κ = ½ THE STRAND TERM IS EXACTLY ZERO BY IDENTITY — AND THE DEADBAND IS THE PART THAT MATTERS.**
`disc = 4·max(0, (κ−½)² − sig2_d)` (`region_init.py:194`) is bit-exactly 0 at κ = ½ for any depth and any
overdispersion. ⭐ But κ is FITTED, so on a genuinely unstranded library κ̂ ≠ ½, and the deadband
subtracting κ̂'s own sampling noise is what keeps the channel silent. A design that reasons only from the
algebraic identity will be surprised by a fitted κ̂.

**④ ⛔ A REAL DEFECT THE SWEEP FOUND, and stage 2 sidesteps it by restriction rather than by fixing it.**
The two precisions the code carries for the strand channel disagree about κ = ½ and only one is honest;
at an AMBIG slot `f_g` still moves by ~2.5 nats on strand data, **toward gDNA**, while
`has_own_composition_evidence` reports no own evidence there. Stage 2's single-stranded restriction
avoids it. ⚠ Any later stage that admits AMBIG slots inherits it, and it needs its own precision rather
than a zero.

**⑤ ⭐⭐ THE SUBSTRATE-HOME DECISION HAS A DIRECT PRECEDENT, AND IT CUTS THE OTHER WAY THAN I ASSUMED.**
`splice_graph.build_region_partition_arrays` already produces `region_types` — a per-region structural
CLASS array that crosses the scanner ABI and is inside `partition_hash` (`splice_graph.py:910-953`). So
an index-time structural class array is an established pattern with a hash that already covers it.
⛔ **But the per-BOUNDARY structural flags were DROPPED from the index schema at
`INDEX_FORMAT_VERSION 7`** and survive only in `edges.feather` (`region_arrays.py:139-143`, which records
the deletion), so a per-boundary substrate array would mean re-adding schema that was deliberately
removed. ⭐ That asymmetry — regions have a home, boundaries do not — is the real content of
`substrate-home`, and it was not visible before the sweep.

**⑥ ⛔⛔ `mrna_active` MEANS TWO DIFFERENT THINGS ON THE TWO AXES, and Stage 0 depends on the boundary
one.** At a REGION it is a MEMBERSHIP test (this region carries an exon bit); at a BOUNDARY it is the
per-strand AND over the two flanks — *mature RNA may cross here contiguously*
(`region_geometry.py:743-755`). ⭐ The boundary reading is exactly the owner's "no exon|exon boundary
component", so §B.3's predicate is right — but it is right for a reason that does not transfer to the
region axis, and anyone reusing the bit must know which axis they are on.

**⑦ A REGION CANNOT SEE ITS OWN FLANKING BOUNDARIES' FLAGS.** `boundary_flags` is 0 at every REGION slot
(`region_geometry.py:756-758`), so `solvable_exon` must GATHER through `chain.left` / `chain.right`
(`region_chain.py:51-56`). Mechanical, but it decides where the predicate can live: it needs the chain,
not just the region arrays.

**⑧ ⚠ THE CLOSEST PRIOR ART IS A FILE I DELETED EARLIER IN THIS SESSION.**
`scripts/design/exon_solvability.py` — retired because its docstring verdict was not reproducible from
the file — is, per the sweep, "the closest existing thing to the owner's solvable-exon predicate". ⛔ Read
it from git (`git show 3176f24d:scripts/design/exon_solvability.py`) before writing the predicate: it
already had truth-free predictors including hops-to-nearest-pure-gDNA-slot and the flank's sj/terminus
structure, and its measured finding was that the STRAND BIT dominates — which is why the owner's rule
leads with single-strandedness.

⭐ **Stage 1 is further along than the plan assumed**: `g1_locked = ¬free_pos ∧ ¬free_neg`
(`region_geometry.py:531-551`) already pins `{0,0,1}` at `Var(log f_g) = 0` on BOTH axes, so an
`intergenic|exon` BOUNDARY is already fixed at 100 % gDNA with zero variance in the belief. What is
missing is not the pin but the EMISSION licence's scoping — see §C.

---

## §F THE TWO RUNNING LOGS (owner-requested, updated as the design lands)

### F.1 ⛔ OBSOLETE — superseded by this design

* ✅ **EXECUTED 2026-08-23: the four pricing arms deleted from `ladder_arm_ab.py`**
  (`stage2_intron_ref{,_r}`, `stage1_pair{,_onesided}`) — their mechanism shipped as
  `measured_intron_reference` (default ON) or was refused, and under the shipped default they could
  only raise NEVER-FIRED. Verdicts survive in the config docstring, §C.1/§D.4 here, and ROADMAP
  rank 11; implementations in git one commit back.
* **`structural_reference_location`'s constant 0.75 — NOT deletable, and here is the boundary of the
  supersession**: the measured location overrides it ONLY at ss-intron REGIONs; the constant is still
  the reference for every exon and every intron BOUNDARY (where the boundary-inclusive measured form
  is measured-refused). The function survives as the measured form's ``base``.
* **`region_init.strand_evidence`'s `locked = ~solvable` feeder for `struct_lock`** — ⛔ **NOT YET
  obsolete: the re-pricing WITH the stage-2 prior as the replacement load was run (2026-08-26) and
  REFUSED** — see §C.1. The mis-scoped mask stays load-bearing; whatever eventually replaces it needs a
  better level channel than the one-pseudo-fragment reference (rank 11 stays open with the updated
  brief). ⭐ Stage 0's contribution stands: gene-edge BOUNDARIES carry EXACTLY ZERO RNA fragments on
  32/32 conditions, so the load is about zero-count slots emitting certainty, not contaminated counts.

### F.2 ⭐ SIMPLIFICATION — smaller, clearer, fewer branches

* **Instrument-side population predicates now have ONE gated home.** `structural_claims` (layer 3) is
  the single definition of intergenic / ss-intron / solvable-exon membership — the thing the retired
  `exon_solvability.py` re-derived privately, and the population the background fit and the stage-2
  prior both select on. An instrument that needs the class reads the module instead of restating bits.
* **The stage-2 prior needs NO fit machinery for its location** — `rho_bg` is the ratio of sums over
  the stage-0 intergenic REGION class, the same estimand `fit_intron_background` measures (its docstring
  says the pools are measured identical), so the location costs one pooled division rather than a
  second fitting path.
* ✅ **EXECUTED 2026-08-23 (owner's go): ψ's location plumbing collapsed to ONE construction site.**
  `solve_chain` lost the `structural_reference` bool and its fallback branch; `calibrate` constructs
  the location (structural base, measured override) and passes `reference_location`. Gated on
  BIT-IDENTITY: all 52 array/scalar fields of pass-0 and final hashed equal on a real ladder condition
  under the shipped config, suite 3,684 / 9 xfail unmoved. ⚠ The two instruments that tally the
  shipped location path (`vertex_ceiling`, `quant_accuracy`) were repointed at `calibrate`'s binding in
  the same batch — patching `sweep`'s would have counted nothing (`TRAPS: an-ablation-that-never-ran`).
* ✅ **EXECUTED 2026-08-23: the two wrong-fact `message_propagation` comment blocks corrected**
  (`calibrate.py`, `zero_controls.py`) — both now state the shipped default (True, RelayPolicy, since
  2026-08-18) and keep the 2026-08-07 OFF measurement as dated history.

---

## §H STAGES 3–4 — THE AGREED DESIGN (owner's plan + review, settled in conversation 2026-08-23)

⭐⭐⭐ **SCOPE (owner): this is PASS-0 ONLY.** It solves the stage-0 substrate and nothing else; its
product is the training set for the gDNA landscape that can see enriched-vs-depleted, and the full
solve comes AFTER the refit with that prior. Unsolved slots stay unsolved — that is a feature.

**Stage 3 — intron → flanking `intron|exon` boundaries.** The intron's ψ posterior is delivered as a
factor on the boundary's λ grid (the EXACT form of the Beta-Binomial predictive — the ``N²·Var(p)``
epistemic term falls out automatically; no resampling machinery, no ``(α₀,β₀)`` chosen), widened by a
per-hop-type DRIFT variance **measured by `hop_currency.py`**, never chosen — the honest cost of
imputation the owner requires. Combined at the boundary with its own strand likelihood + reference by
the existing lattice product. ⭐ The premise (the intron's contained fragments and the boundary's
unspliced crossings are ONE population) is CERTIFIED by the stage-0 matrix, 32/32. ⛔ The ghost:
`intron_phi` (2026-08-04) refuted the UNCONDITIONAL form, split BY STRAND — priced per stratum from
day one, unstranded rows the watch-list; the fuse, the honest precision, and the repaired intron
beliefs are what changed.

**Stage 4 — boundary → solvable exon.** AXIOM 0 governs the wording: one RNA, arriving by two ROUTES
(unspliced contiguous crossing, spliced junction crossing). At the boundary the routes MERGE into one
RNA quantity (owner's design), under two guards, both measured rather than stylistic:

* **Guard 1 — merge in DENSITY space, never raw counts.** Each route's count over ITS OWN crossing
  opportunity (`crossing_eff_length` for contiguous, `sj_opportunity` for spliced) — the two footings
  differ by a MEASURED 7–10×, enough to swing a transferred composition several-fold. In density
  units the routes are genuinely indistinguishable, which is the owner's point made exact. Variances
  ADD across the summed routes (a sum of components, not a combine of witnesses — the owner's
  "precisions are not additive" intuition, formalized). The delivered claim is the DENSITY RATIO,
  re-formed as a composition through the exon's own opportunities — the absolute level cancels, which
  is what makes the transfer capture-robust.
* **Guard 2 — sidedness is STRUCTURAL.** Two-sided iff every RNA route into the exon is accounted
  (each flank is an eligible splice site or a terminus); ONE-SIDED ("at least this much RNA") when
  an unaccounted flank exists — §1.1's inequality, the only form the zero-gDNA control has ever
  endorsed on every row. Decidable from stage-0 structure, no constant.

The exon's solve is then the lattice product of the flank message(s), its own strand likelihood, and
the reference. **After stage 4: fit the gDNA landscape on SOLVED slots only** (rank 4's ask — kills
the fit-on-your-own-belief circularity).

### §H.1 STAGE 3 STATUS — BUILT AND GATED (2026-08-23); pricing in flight

* **Landed**: `messages/fanout.FanOutPolicy` behind `message_policy = "fanout"` (shipped default
  untouched), `StepContext` gained the `claims` constant (built once in `solve_chain`; policies that
  never read it are unaffected — suite bit-identical), gate `tests/calibration/test_fanout_policy.py`
  (7 cases written first and verified failing; integration through the REAL backbone asserts ONLY the
  claimed boundaries move). Suite 3,696 / 9 xfail, +12 re-derived exactly.
* **The perturbation sweep simplified the policy**: the destination mask and the `tau_lam > 0` guard
  were PROVEN to do no work (the destination is implied by source-licence + adjacency — derivation in
  the policy docstring — and a zero-precision claim is no claim) and were DELETED; the `valid` mask
  stays as the `NeighbourState` contract. All four surviving clauses fire.
* **First pricing (test chromosome, 16 conditions, vs the relay base)**: pooled deliverable
  stranded × OFF −54.9 %, stranded × ON −26.7 %, deferred −6.6 %, unstranded × OFF +5.0 % on a small
  base. ⛔ **The pooled `g00` row (+25,146 %) is an ARTIFACT of the comparison, not a stage-3 defect**:
  outside its claimed slots the fan-out is silence BY DESIGN (its g00 deliverable is bit-equal to the
  messages-OFF run), while the relay base sweeps the whole chain — so the pooled instrument charges
  pass-0 for every slot it deliberately leaves unsolved. The honest verdict is §E's rule — claimed
  slots only.
* ⭐⭐⭐ **THE CLAIMED-SLOT VERDICT (ladder, misplaced gDNA fragments at the 17,518
  `ss_intron_boundary` slots, silent / relay / fanout)**: `g50` capture-OFF both strand settings —
  **fanout best, beating BOTH silent and the relay** (31,121 / 31,529 / **24,793** at ss0.50; 25,561 /
  28,266 / **20,944** at ss0.99); `g00` ss0.50 OFF — **fanout best** (100,297 / 8,716 / **7,860**),
  so the pooled g00 artifact is fully explained: at its OWN slots the fan-out is the best policy at
  the zero-gDNA control. `g50 ss0.99` capture-ON — silent best, fanout close (+7.7 %, 38,007 vs
  40,933), relay worst (66,637): under capture the intron's message is slightly net-harmful at
  boundaries vs silence — small (0.26 % of the crossing mass) but real; the stage-4 drift question's
  kin. Deferred stratum (`g50 ss0.50` ON) — the relay dominates (215,475 vs fanout's 828,266), as its
  recorded value-concentration predicts; fanout still beats silence (1,082,469).
* ⭐⭐⭐ **THE FULL 16-CONDITION CLAIMED-SLOT TABLE (`stage3_claimed` script, 2026-08-23)**: fanout is
  the BEST of the three policies on **all eight rungs of the two capture-OFF strata** (including the
  stranded `g00` row at 482 vs silent 4,323 / relay 2,360, and `g98` at −37 %/−51 %); on stranded ×
  capture-ON it is +3.9 % vs silent summed over the stratum (best at `g00`, slightly behind at
  `g50`/`g98`) with the relay far worst (−39 % fanout vs relay); the deferred stratum stays the
  relay's, with fanout still halving silence there.
* ⭐⭐ **POOLED LADDER (context, `--arm policy_fanout` vs the relay base): every IN-SCOPE stratum
  improves on the deliverable, 18/18 rows** — stranded × ON −39.1 %, stranded × OFF −39.9 %,
  unstranded × OFF −36.8 %. The deferred stratum (+199.6 %) and pooled `g00` (+611 %) regress because
  the fan-out deliberately does NOT do the relay's whole-chain imputation — in the redesign's flow
  that job belongs to the post-refit FULL solve with the landscape prior, not to pass-0. ⛔ So
  `"fanout"` is NOT a candidate shipped default in today's flow; it is pass-0's policy, and the
  default question arrives only with the full redesigned flow.

### §H.2 ⭐⭐⭐ STAGE 4 — THE DESIGN (2026-08-23, design-only turn; no code yet)

**One sentence:** the eligible flank tells the exon its composition, as a density ratio formed from
everything the flank sees — and a structural bit says whether that account is COMPLETE (a two-sided
estimate) or PARTIAL (a one-sided at-least-this-much-RNA bound).

#### The completeness theorem — the sidedness rule, derived rather than chosen

A molecule overlapping the exon either **contains the delivering flank** (crosses it contiguously, or
splices through it — both measured at the flank) or **originates exactly at it** — because the exon's
interior holds no boundary, a molecule that does not contain the flank must have its end AT the flank,
which is a transcript-end bit the graph already carries (`region_geometry._RNA_LOW_END` for a left
flank, `_RNA_HIGH_END` for a right). Therefore:

    a flank's account of the exon's RNA is COMPLETE  ⇔  the flank carries no transcript-end bit
                                                         facing the exon (on the exon's strand)

Complete ⇒ the transfer is a genuine TWO-SIDED estimate. Incomplete ⇒ RNA invisible to the flank may
exist, and the honest claim is ONE-SIDED — §1.1's inequality, the only form the zero-gDNA control has
ever endorsed on every row. ⭐ No constant, index-derived, and falsifiable with NO solver: at complete
flanks certified truth must satisfy the equality (flank-implied RNA ≈ exon RNA, within the noise
floor); at incomplete flanks the inequality's DIRECTION must hold. Worked check: TA's FIRST exon
(TSS at its left flank) delivering through its right donor is COMPLETE — every molecule covering it
continues through the donor — so the TSS costs nothing when the claim comes from the other side.

⭐ **A derived simplification (2026-08-23): the both-strand end-bit masks suffice.** At a LICENSED
flank, an other-strand transcript-end in the facing set would have to extend INTO the exon, which
makes the exon AMBIG and therefore unlicensed — so at licensed flanks the strand-specific check and
the plain ``_RNA_LOW_END`` / ``_RNA_HIGH_END`` masks are exactly equivalent, and no strand plumbing
is needed. The derivation lives in the code's docstring; no test can distinguish the two forms at a
licensed flank, which is the proof the simpler form is the right one.

#### The owner's scope ruling (2026-08-23), folded in

* **Imputation comes ONLY through eligible flanks** — one eligible flank ⇒ one claim; two ⇒ the
  precision-weighted fuse of two independent claims (the established pattern). An INELIGIBLE flank
  contributes NOTHING in pass-0: at a terminus there is no spliced measurement of RNA and therefore no
  composition measure at all, and reaching such slots is a MULTI-HOP problem that belongs to the full
  forward-backward relay run later WITH the landscape prior — never to pass-0.
* **Pass-0's product is the most CONFIDENT substrate derivable, not the largest**: it exists to train
  the gDNA landscape that must see enriched capture against the depleted floor. Widening the exon set
  would trade substrate confidence for coverage — refused.
* ⭐ **The completeness bit therefore carries a SECOND duty**: complete-flank exons enter the landscape
  TRAINING set as estimate-grade; incomplete-flank exons are bound-grade — naturally down-weighted by
  their one-sided-wide posteriors, or excluded outright (one line, priced at 4g when the training-set
  choice is made).
* **The three primitives (fuse, route-merge, reframe) are deliberately REUSABLE** — validated here,
  they are the candidate vocabulary for the rebuilt full relay downstream.

#### The hop, as three primitives — stage 3 and stage 4 are ONE operator with different currencies

1. **fuse** — precision-weighted λ combine (exists; stage 3's deliver, the relay's recorded idiom).
2. **route-merge** — at the flank, every route's count over ITS OWN opportunity (unspliced crossing
   over `crossing_eff_length`; spliced flux over `sj_opportunity`, taking the FACE that enters the
   exon — `inv_sj_lo/hi` by direction, hop_currency's direction-dependence); densities ADD, variances
   ADD (a sum of components, not a combine of witnesses — the owner's non-additive-precisions
   intuition, formalized). One RNA class throughout; AXIOM 0 untouched.
3. **reframe** — the density RATIO re-formed as the exon's composition through the EXON's own
   opportunities (destination CONSTANTS ✓). The absolute level cancels, which is what makes the
   transfer capture-robust — certified by the currency map (COMPOSITION beats LEVEL 20–30× under
   capture at this hop).

Stage 3 = fuse alone (identity currency: same population, same frame). Stage 4 = fuse at the boundary
(intron claim ⊗ the boundary's own strand λ — the scan's `step()` becomes real, ~15 lines), then
route-merge ∘ reframe at delivery, under the completeness sidedness. Depth stays structurally 2:
claims die at exons because exons publish nothing. ⭐ One policy, one sweep, no calibrate flow change.
⚠ Deliberately NOT transported: the boundary's reference term — a message carries DATA, and composing
own-λ (data curvature only) rather than the ψ posterior keeps the reference from being counted twice.

#### The drift (the imputation cost) — measured 5.1 %/10.2 % excess capture-ON at this hop, ~0 OFF

* **v1 (recommended): NO drift term** — the reframe through opportunities is a sharper transport than
  the currency map's plain-composition rule, so part of the measured excess may already be the frame
  error the reframe removes. Measure first.
* Fallback: a capture-gated variance derived from the measured excess (a measured constant; owner
  sanction required under no-magic-numbers).
* ⛔ REFUTED in advance: a flank-disagreement-fitted drift — both flanks biased the same way under
  capture is exactly `TRAPS: a-variance-cannot-fix-a-bias`'s recorded case.

#### What the codebase SIMPLIFIES (the §F ledger's stage-4 entries, prospective)

* `ONE_SIDED_RNA` (a process-global toggle) → a per-message sidedness (one optional `PsiMessage`
  field, default None ⇒ byte-identical); the global dies with the relay.
* The fan-out at depth 2 is the honest rebuild of the relay's splice-in + certified-RNA + reframe
  machinery. When the full flow lands and is priced, the deletion candidates are: `RelayPolicy`'s
  operator set, the mass rescale, the reframe and its licences — and **rank 11's cancelling pair
  DISSOLVES rather than being repaired**, because the fan-out never consults `struct_lock` for
  emission. ⛔ All gated on the deferred stratum (the relay's one stronghold) and the owner's flow
  ruling; nothing is deleted on promise.

#### The pieces, in order, each with its falsification

* **4a — the completeness bits** in `structural_claims` (structure only): designed-GTF unit cases
  (first/last exon, retained intron, the 1 % both-bits flank), then the confusion-matrix extension on
  both panels (equality at complete flanks, direction at incomplete) — stage 0's discipline replayed.
* **4b — route-merge as a pure function** (counts, opportunities, direction faces → densities +
  variance): brute-force unit tests; the face selection tested explicitly.
* **4c — reframe as a pure function** (+ variance propagation): hand-derived cases, and the
  capture-invariance gate — scaling both components by one factor must leave λ unchanged exactly.
* **4d — ψ's per-message sidedness**: the optional field, byte-identical when absent; the one-sided
  behaviour inherits the existing `ONE_SIDED_RNA` test machinery.
* **4e — wire into `FanOutPolicy`**: `step()` composes at boundaries; `deliver()` runs 4b∘4c at exon
  destinations under 4a's licence; integration + perturbation sweep as stage 3's.
* **4f — pricing**: PROMOTE the scratch claimed-slot scorer into `scripts/design/` (one instrument
  covering the boundary AND exon claimed populations, silent/relay/fanout, DELIVER/REFUTE split at
  truth-pure slots); test chromosome → ladder; both zero controls.
* **4g — the flow's ending** (separate stage): the landscape fit on SOLVED slots only.

**Owner sign-offs wanted before code:** ① the per-message sidedness field (it touches ψ);
② drift v1 = none, fallback measured-capture-gated; ③ the two-sided channel is λ, the one-sided
claim rides the RNA-share channel with the one-sided residual (the certified-RNA semantics it
generalizes).

### §H.3 ⭐⭐⭐ STAGE 4 STATUS — BUILT, GATED, PRICED (2026-08-23); one measured residual, one lever

All five pieces landed test-first with full perturbation sweeps (`762065bc` 4a, `17bf23a0` 4b+4c,
`5e5f232e` 4d, `d46f285e` 4e, plus `pass0_claimed_ab.py`, the 4f instrument — self-test 6/6, indexed).

**Pooled ladder (context):** vs the shipped relay, the complete pass-0 policy reads −37…−40 % on all
three in-scope strata, 18/18 rows; stage 4's own increment over stage 3 cost ±0.2 % in scope while
improving the deferred stratum −63.5 % and pooled `g00` −56.5 % (8/8) — the exon fan-out earns its
keep exactly where exons are blind. The relay's deferred-stratum edge shrank +200 % → +9.3 %.

**Claimed slots (the §E verdict; B = intron boundaries, E = solvable exons; silent/relay/fanout):**
* **B: the stage-3 verdicts stand** — fanout best on both capture-OFF strata and at `g00`-OFF; the
  stranded × ON +7.5 % residual unchanged; and on the deferred stratum's boundaries the refit now
  carries fanout to near-relay (251,503 vs 215,475, from 828k at depth 1).
* **E, the wins:** deferred-stratum exons **169,659 vs silent 697,694 and relay 258,207** — the
  design's target population, decisively taken; `g00` exons 9× better than silence (21,836 vs
  192,486), though the relay's whole-chain zero-spreading still leads there (1,548).
* **E, the named residuals:** ⛔ **stranded × capture-ON exons REGRESS: 40,528 vs silent 18,042**
  (2.25× at `g50 ss0.99 ON`); ⚠ unstranded × capture-OFF exons +21 % on a small base (7,755 vs
  6,415). Neither is pooled away, and neither gets a mechanism before its dissection speaks.

#### ⭐⭐⭐ THE OWNER'S ANALYSIS FRAME FOR THE STRANDED × CAPTURE-ON RESIDUAL (2026-08-23) — binding

* **Break the stratum BY gDNA RUNG (`g05`/`g50`/`g98` separately), never pooled**, and SURVEY the
  regions and boundaries that regress under the fan-out — the dissection is its own piece of work.
* **The central question: how accurate is the STRAND MODEL ALONE at single-stranded exons?** It is
  understood to be very accurate there — so the imputed composition arriving from the boundaries is
  EXPECTED to be less accurate than the exon's internal strand solve. That is anticipated, not a
  surprise; the question is why the message WINS the fuse anyway.
* **The problem is two-fold, and both halves are open:** ① ACCURACY (bias) — the imputation is good
  (it improves many scenarios, especially unstranded) but less accurate than the strand model here,
  and the brand-new implementation may simply carry bugs; ② PRECISION (too high) — the message's
  delivered precision is evidently high enough to harm the strand model's solve. **The elegant fix is
  honest precision on BOTH sides**: the strand model's exon precision must be honestly HIGH, and the
  composed-and-transferred boundary precision must be honestly (lower) — the computation and
  propagation of the message precision needs careful evaluation both THEORETICALLY and at the
  regressed slots.
* ⚠ **"The drift", defined once so it stops being jargon**: the premise-variance term — the measured
  excess error that remains when the TRUE source value is transported over this hop
  (`hop_currency.py`'s excess-over-floor: 5.1 % ladder / 10.2 % test-chromosome under capture-ON,
  ~0 OFF). It prices "the flank's composition is only locally-approximately the exon's," as a
  variance ADDED to the message — i.e., one principled candidate INSIDE the owner's precision branch,
  not a separate concept and not presumed to be the fix.
* ⛔ **The owner holds a suspicion and will NOT act on it until the dissection data over the
  regressed regions and boundaries is on the table.** No mechanism lands before then.

⛔ Next work, in order: ① the stranded × capture-ON dissection (per rung, per slot, bias-vs-precision
columns); ② the unstranded × capture-OFF exon cell's dissection; ③ then 4g — the landscape fit on
solved slots only, where the estimate/bound grades choose the training set. ⭐ ① and ② are DONE and
their data is §H.4 below; ③ is the next piece of work.

**Implementation shape:** the new message POLICY on the existing backbone (the `the-one-hop` ruling) —
emission licensed only FROM solved-substrate slots TOWARD adjacent unsolved ones, `SilentPolicy`
behaviour elsewhere; the forward scan carries the rightward fan-out, the backward the leftward, so the
two existing scans cover the depth-2 pattern with no backbone change. Stage 3 first; falsification
test first; test chromosome before the ladder; both zero controls on every rung.

### §H.4 ⭐⭐⭐ THE RESIDUAL DISSECTION (2026-08-23) — THE DATA. NO MECHANISM PROPOSED

⛔ This section is the deliverable the owner's frame (§H.3) asked for, and it stops at the data. Nothing
in `src/` was changed; every arm below patches a DELIVERED message or reads a recorded one.

⭐⭐ **THE SURVEY IS RE-DERIVABLE, NOT QUOTED**: it is now `pass0_claimed_ab.py --dissect`, which
prints readings ①–④ below off one cached payload, under two gates — the recorded `deliver` must
belong to the sweep the debug capture describes, and the per-slot `f_g·count` must BE
`mass_gdna_region` rather than merely agree with it. `--self-test` 15/15, and its perturbation sweep
closed two fixture holes it found on the first pass (a `gx == rx` fixture cannot see a dropped SOURCE
frame; a direction-only spliced check cannot see `U` at all, because `U` cancels between the two
unspliced routes — both are now pinned to hand-derived values). ⭐ `--dissect --sweep` adds readings
**⑤** and **④b**, which price the candidate DIRECTIONS rather than describing the tool, and cost ~15
solves — so they are opt-in, and everything in this section is re-derivable from one command.

#### ① WHERE IT IS — per rung, never pooled (stranded × capture-ON, claimed exons, silent → fanout)

| rung | silent | fanout | | DELIVER | REFUTE |
|---|---|---|---|---|---|
| `g00` | 1,749 | **1,344** | −23 % — a WIN at the zero control | 0 → 0 | 1,749 → 1,344 |
| `g05` | 10,344 | 12,419 | +20 % | 325 → 361 | 10,020 → 12,058 |
| `g50` | 18,042 | 40,528 | +125 % | 345 → **201** | 14,003 → 35,068 |
| `g98` | 18,032 | 28,370 | +57 % | 22 → 27 | 18,010 → 28,344 |

⭐ **The damage is REFUTE-ONLY at every rung** — slots carrying real RNA. DELIVER (certified pure gDNA)
is flat or IMPROVED. The regression is not monotone in gDNA: it peaks at `g50`.

#### ② WHAT DOES IT — stage 4's λ transfer, and demonstrably nothing else

* **`fanout_no_exon_lam` is EXACTLY EQUAL to silent** at the exon population, on every rung and at both
  refit settings. Stage 3's boundary message is provably inert here.
* The one-sided branch is worth ±8 fragments of 22,486. The two-sided λ transfer is the whole thing.
* The refit is not the amplifier: the gap is +24,860 at `calib_refit_iters = 0` and +22,486 at 3.
* All **+22,506** of the added error sits on the **4,548** exons that receive a λ message; the other
  5,614 claimed exons move by **−20** in total.

#### ③ THE TWO HALVES, SEPARATELY MEASURED — one table, one basis (`--dissect`, `calib_refit_iters = 0`)

Everything below is mass-weighted over the exons that HEAR a λ claim. ⭐ `median |z|` is the only
column that judges the precision against ANYTHING: `z = (lam_mode − lam_true)·sqrt(lam_prec)` at
interior-truth exons, and an honestly-calibrated claim reads **0.674**.

| condition | own `\|f−true\|` | message `\|f−true\|` | signed (msg) | `lam_prec/tau_lam` med | median `\|z\|` ⇒ over-claim | gDNA route |
|---|---|---|---|---|---|---|
| `g05` ss.99 ON | 0.0095 | 0.0234 | −0.018 | 7.8 | 2.19 ⇒ **11×** | −0.98 / −0.80 |
| `g50` ss.99 ON | 0.0131 | **0.1186** | **−0.110** | 14.9 | 4.50 ⇒ **45×** | −0.96 / −0.98 |
| `g98` ss.99 ON | 0.0128 | 0.0963 | −0.096 | 9.7 | 2.58 ⇒ **15×** | −1.01 / −1.01 |
| `g50` ss.99 OFF | 0.0127 | 0.0171 | +0.001 | **47.5** | 1.76 ⇒ 6.8× | **−0.07 / −0.09** |
| `g50` ss.50 OFF | **0.4640** | 0.0183 | +0.004 | — (own τ ≡ 0) | 1.48 ⇒ 4.8× | **−0.01 / −0.06** |

⭐⭐⭐ **READ THE LAST TWO ROWS AGAINST THE FIRST THREE — they are the whole verdict.** The stranded
capture-OFF row carries the HIGHEST precision ratio on the panel (47.5, the message outweighing the own
side at 94.7 % of slots) and the fan-out **WINS** that cell at `refits = 0` (6,671 against silence's
7,473). ⛔ **So a dominant message is not itself the fault; a dominant WRONG message is.** The bias
column, not the precision column, is what tracks the sign of the outcome.

⚠ The own side's effective precision was also measured directly rather than assumed — from the
observed posterior shift, `t_eff = p·(m − λ₁)/(λ₁ − λ₀)`, which agrees with `region_init.tau_lam` to a
median 1.9×, so `tau_lam` is a sound and conservative proxy. On that footing the message wins the fuse
at 77 % of slots at `g50` ON — but **at the top 2 % of slots by error the ratio is 0.83**, i.e. the
message LOSES the fuse there and still does the damage, because a claim off by −0.11 in `f_g` carried
at even half weight moves a 10⁴-fragment exon by hundreds. ⛔ A precision repair alone cannot reach the
mass.

#### ④ WHY THE CLAIM IS WRONG — the reframe does NOT cancel capture, and that is measured

Push the flank's **CERTIFIED** composition through the transfer instead of its solved one, and ask how
much of the delivered bias survives a PERFECT source. `delivered f_g − exon truth`, per side:

| condition | SOLVED source | CERTIFIED source | share of the bias the FRAME owns |
|---|---|---|---|
| `g05` ss.99 ON | −0.021 / −0.020 | −0.014 / −0.013 | **66 % / 62 %** |
| `g50` ss.99 ON | −0.110 / −0.122 | −0.063 / −0.073 | **57 % / 60 %** |
| `g98` ss.99 ON | −0.107 / −0.107 | −0.020 / −0.022 | 19 % / 20 % |
| `g50` ss.99 OFF | −0.001 / −0.000 | +0.003 / +0.003 | both ≈ 0 |

⭐⭐ **At the two rungs where the claim does most damage, a PERFECT source composition removes only
about 40 % of the delivered bias — the majority belongs to the frame.** At `g98` the ordering reverses
and the source's own error dominates, so the frame is not the whole story at every rung.

⛔ **A CORRECTION, recorded because it was nearly written into this document.** A first pass measured
this on the slots where a given SIDE's transfer was live rather than on the slots that actually HEAR a
claim, and on that population the certified source read WORSE than the solved one — which would have
supported "the two errors cancel, so fixing the boundary solve first makes the exon delivery worse."
**That is false**; it does not survive on the population the verdict is about, and `--dissect` now
prints both transports side by side so the comparison cannot be made on a private slot set again.
(`TRAPS: the-intermediate-is-not-the-deliverable`, in its population costume.)

The `gDNA route` column in ③ is `log(transported ρ_g / the exon's own certified ρ_g)` per side.
⭐⭐⭐ **It sits at `e^{−1.0}` — 2.7× too LOW — at EVERY capture-ON rung regardless of gDNA level, and
at `e^{−0.07}` capture-OFF.** The merged RNA route moves far less (−0.12 to −0.34), so the net odds
run ≈ 2.2× short of gDNA. **The geometric frame is sound; the deficit is capture and only capture.**

⛔ **The reason is structural: `eff_gdna` and `eff_rna` carry NO capture term.** `region_geometry`'s
`divisor` is `contained_eff_length` at a REGION and `crossing_eff_length` at a BOUNDARY, both purely
geometric — `capture_eff_length.py` is a TRANSCRIPT-level EM length for `assemble_priors`, not the
per-slot calibration opportunity. So §H.2's "the absolute level cancels, which is what makes the
transfer capture-robust" holds only under a UNIFORM capture pull, and a probe-covered exon interior
against an intron-spanning flank crossing is not uniform — the two routes are pulled differently
because the sj route's molecules live on the spliced (fully exonic, fully probed) template while the
crossing route's straddle into unprobed intron.

⚠ **`hop_currency.py`'s recorded 5.1 % / 10.2 % excess badly understates this hop for this
population** — measured here it is ~1.0 nat on the gDNA route. Do not quote the recorded figure as this
hop's drift.

#### ④b A NAMED, SEPARABLE SLICE OF THAT BIAS — the reframe constant, built from FITTED fl PMFs

`K = log(E_rx·E_gc) − log(E_gx·E_rc)` enters `lam_e` with coefficient **exactly 1** and carries **no
variance term at all**, and all four effective lengths come from FITTED fragment-length pmfs. Rebuild
the four from the simulator's own post-capture pmfs and deliver the difference as a pure mode shift:

| condition | fitted gDNA / RNA pmf mean | truth | median `dK` | share of the regression |
|---|---|---|---|---|
| `g05` ON | 254.04 / 221.53 | 240.53 / 229.36 | **+0.495** | 13.9 % |
| `g50` ON | 254.81 / 221.37 | 240.40 / 229.23 | **+0.467** | 21.1 % |
| `g98` ON | 254.85 / 221.61 | 240.43 / 229.47 | **+0.504** | 24.1 % |
| `g50` OFF | 216.58 / 204.09 | 216.67 / 212.66 | **+0.132** | n/a — the fan-out WINS this cell |

⭐ Under capture the fitted gDNA length model runs **6.0 % LONG** and the RNA one **3.4 % SHORT**;
capture-OFF the gDNA fit is exact to 0.04 %. So the mis-fit is itself capture-specific, and `K`
amplifies it into a stable **~0.5-nat shove toward RNA** at every capture-ON rung — a fifth of the
regression, in a term that claims no uncertainty whatsoever. ⛔ This is a DIAGNOSTIC: a truth pmf does
not exist at run time. It prices the frame constant's error; it does not propose reading truth.

#### ④c THE DELTA METHOD IS RUN ~175× OUTSIDE THE DOMAIN ITS GATE VALIDATES

`test_flank_to_exon_variance_matches_monte_carlo` pins `tau_b = 30, U = 400, S = 250`. Production runs
at **`tau_b` median 0.171, `U` median 51, `S` median 4**. Re-running that same brute force at the
production regime:

| regime | delta λ | MC λ | `1/tau` (delta) | MC var | MC/delta |
|---|---|---|---|---|---|
| the GATE's | −0.791 | −0.796 | 0.01369 | 0.01375 | **1.00** |
| `tau_b=0.17 U=51 S=4` | 1.563 | 1.134 | 1.055 | 1.793 | **1.70** |
| `tau_b=0.17 U=51 S=333` | −2.440 | −2.735 | 0.0686 | 0.4876 | **7.11** |
| `tau_b=0.41 U=103 S=333` | −1.657 | −1.688 | 0.01313 | 0.02052 | 1.56 |
| the tail, `tau_b=4.2` | −2.927 | −2.928 | 0.002597 | 0.002601 | 1.00 |

⛔ **The linearisation biases the MODE as well as the variance** (1.563 against 1.134; −2.440 against
−2.735), so this is not only a precision fault. ⭐ The gate is honest about what it tests and blind to
where the code runs — `TRAPS: a-falsification-test-needs-perturbation` applied to the FIXTURE's
operating point rather than to its assertions.

⚠ **④c accounts for only 1.6–7× of ③'s 44× over-claim; the remainder is the frame error of ④/④b,
which no sampling model contains.** The own side, by the same reading, is honest to slightly
conservative (median |z| ≈ 0.64).

#### ④e WHERE THE DAMAGE SITS STRUCTURALLY

**84.1 % of the added error (+18,916 over 3,351 exons) is on SINGLE-complete-flank exons**; the
two-sided fuse carries +3,590 over 1,197. ⛔ Any account resting on the two-flank fuse is capped at
16 %. ⚠ And the message-receiving exons are **not** RNA-dominated as a population: certified `f_g` is
**0.515 count-weighted, 0.474 error-weighted** at `g50` ON — capture puts probes on exons, so at 50 %
gDNA the exon is where the gDNA is. The same exons sit at median 0.038 capture-OFF, which is the other
half of why the identical machinery wins there.

#### ⑤ NEITHER HALF ALONE CLOSES IT (sweep at `calib_refit_iters = 0`; the silent floor is 17,359)

| `lam_prec ×` | 1 | 0.3 | 0.1 | 0.03 | 0.01 | 0 |
|---|---|---|---|---|---|---|
| claimed-exon error | 42,219 | 27,595 | 20,641 | 17,758 | **17,108** | 17,388 |

| `lam_mode +` (nats) | 0 | 0.5 | **0.8** | 1.0 | 1.5 |
|---|---|---|---|---|---|
| claimed-exon error | 42,219 | 30,543 | **27,810** | 28,883 | 37,235 |

* **Precision alone**: a **100×** cut is needed merely to reach parity, and the best attainable is
  **−1.4 %** against silence. Honest precision makes this message INAUDIBLE, not useful.
* **A scalar bias correction alone**: the optimum is **+0.8 nats**, which independently reproduces the
  measured median odds deficit — but it is still **+60 %** over silence, because the deficit varies
  slot to slot and a scalar cannot fix a scattered error.
* ⛔ `TRAPS: a-variance-cannot-fix-a-bias`, measured at this hop, in both directions.

#### ⑥ THE IMPLEMENTATION IS CORRECT — four plausible bug sites re-derived, none defective

* **the sj FACE.** A LEFT flank reads `sj_count_hi`, and that is right: `jc_hi` is filed at the
  ACCEPTOR (`region_geometry` fills `jc_lo` at the donor, `jc_hi` at the acceptor), and a molecule
  splicing INTO an exon arrives at its left boundary as an acceptor. The RIGHT flank mirrors.
* **`S_all[src, col[rows]]`.** `rows = arange(n)`, so this is `S_all[src, col]` — the SOURCE
  boundary's flux in the DESTINATION exon's strand column. That is the intended pairing; a licensed
  exon is single-stranded so `col` is well defined, and every other slot is masked out.
* **the delta-method Jacobians.** `∂λ_e/∂λ_b = (1−f) + f·ρ_ν/ρ_r`, `∂λ_e/∂log U = w_μ`,
  `∂λ_e/∂log S = −w_μ` — all three as coded.
* **the one-sided Jacobian.** `x = log(1−f)`, `dx/dλ = −f`, so `precision(x) = τ_λ/f_g²` — as coded.
* **the precision compose does not double-count.** A boundary's own `tau_lam` enters its composed
  state once; the stage-3 branch reads RAW region claims at boundary destinations (regions are never
  written by the scan), and an exon's two flanks are distinct slots.

⭐⭐ **So the defect is in the model's PREMISE, not in the new code.**

#### ⑦ THE UNSTRANDED × CAPTURE-OFF CELL IS A DIFFERENT MECHANISM, and smaller than it looked

Per rung, claimed exons, silent → fanout: `g00` 192,486 → **21,836** (8.8× better); `g05` 1,922 →
2,024 (+5 %); `g50` 6,415 → 7,755 (+21 %); `g98` 7,983 → **6,193** (−22 %). ⭐ **The recorded "+21 %" is
ONE RUNG**; over the three contaminated rungs the stratum is 16,320 → 15,972, a small net WIN.

⛔ **At κ = ½ `strand_evidence`'s `disc = 4·max(0, (κ−½)² − σ²_d)` is IDENTICALLY ZERO, so `tau_lam` is
EXACTLY 0 at all 5,063 message-receiving exons** — the exon has no own composition evidence on λ at
all, and "the message's precision is too high" is not the complaint here. **Here the message is the
ACCURATE voice**: mass-weighted `|f_g − true|` is **0.0183 for the message against 0.4640 for the own
solve** (signed +0.433 — the own side floods these exons with phantom gDNA). The message is **25×
more accurate** than what it displaces.

**The competitor is the REFIT.** At `calib_refit_iters = 0` silent reads 268,498 and the fan-out
**20,809 — 12.9× better**; at 3 refits silent falls to 6,415 and the fan-out only to 7,755. The
refitted gDNA prior imputes these exons better than the flank message does, and the message DISPLACES
it. ⭐ The one-sided bound is a small WIN here (20,809 against 22,059 with it muted) — the opposite of
its stranded irrelevance. ⭐⭐ **The same pattern appears at STRANDED × capture-OFF `g50`**: fan-out
6,671 against silence's 7,473 at `refits = 0`, and 5,935 against 5,460 at three. So a message and a
refitted prior are **two imputations competing at the same slot with nothing arbitrating them** —
a design question in its own right, independent of the capture bias.

The drift is the same mechanism at a tenth the size: gDNA route median −0.34 nats capture-OFF.

#### ⑧ WHAT IS COMMON TO BOTH RESIDUALS, in one sentence

`flank_to_exon_lambda`'s variance prices SAMPLING NOISE only — it is a statement about the source's
counts and never about whether the flank's composition is the exon's — so the delivered precision is
unrelated to the error the claim actually carries. In the stranded stratum that beats an accurate
strand model; in the unstranded one it beats an accurate refitted prior. **One missing term, seen
through two different competitors.**

⚠ **And the frame's share is rung-dependent, which a single ruling has to accommodate**: at `g05` and
`g50` capture-ON the frame owns ~60 % of the delivered bias and a perfect source would not rescue the
cell, while at `g98` the source's own error owns ~80 % of it. So neither "fix the frame" nor "fix the
source" is sufficient alone, and ④'s table is the thing to re-run after any change to either.

### §H.5 ⭐⭐⭐ THE PRECISION CHAIN, THE REFERENCE PRIOR, AND WHAT A DAMPENER WOULD DO (2026-08-23)

The owner's follow-up to §H.4, answered against the code. ⛔ Still no mechanism in `src/`.

#### ⓪ THE PROBE MECHANISM, MEASURED ON THE PANEL — and it is NOT mainly the splice-junction case

⭐⭐⭐ **Every probe base on this panel lands in an EXON; introns and intergenic carry EXACTLY ZERO.**
Measured directly off `reference/capture_panel.bed` intersected with the index's region signatures:

| region class | regions | bp | probe bp | coverage | regions with any probe |
|---|---|---|---|---|---|
| intergenic + the three intronic signatures (0, 4, 8, 12) | 11,117 | 91,945,594 | **0** | **0.00 %** | **0.0 %** |
| the twelve exonic signatures | 24,018 | 5,665,613 | 1,568,781 | 23–50 % | 33–59 % |

⛔ **So an intron|exon boundary is a PROBE CLIFF by construction.** An unspliced crossing fragment
straddles into a region with no probe at all, while the exon interior it speaks to may be a quarter
probe-covered at a binding weight of ~1200:1 against off-target. The flank's gDNA density is therefore
systematically below the exon interior's — which IS the measured `e^{−1.0}` gDNA-route deficit, and it
does not require a splice-junction-spanning probe at all.

⚠ **The panel DOES express the owner's case 3** — probes are tiled in TRANSCRIPT coordinates, so
**3,348 of 13,824 (24.2 %, and 24.2 % of probe bp) are multi-block in genomic space**, and the sampler
gives such a probe a `gdna_split_penalty` on the gDNA axis while mature RNA keeps full binding. But
matched on local probe fraction the capture-attributable deficit is **larger where case 3 is ABSENT**
(single-block flanks ≈ −0.76 nats against SJ-spanning ≈ −0.41). ⭐ **So case 3 is real and adds to the
bias, but case 1's converse — probes inside exons, none in introns — is the dominant term on this
substrate.** ⛔ That matters for scope: a fix aimed only at splice-junction probes would address the
smaller half.

#### ① HOW THE DELIVERED PRECISION IS BUILT — three steps, and each one's assumption

**STEP 1 — the intron's own λ precision** (`region_init.build_region_init`, `strand_evidence`):

    tau_lam = [single-strand] I_strand  +  I_density
    I_strand = N_eff · disc · [f_g(1−f_g)]² / (4 p(1−p)),   p = kappa + f_g(½−kappa)
    N_eff    = N / (1 + (N−1)·od_r)        disc = 4·max(0, (kappa−½)² − sigma²_d)

*Assumes*: the strand counts are the only composition evidence, and the overdispersed effective count
is the honest sample size. **This step is sound** — §H.4 ③ measures the own side at median |z| ≈ 0.64,
honest to slightly conservative.

**STEP 2 — the compose at the flank** (`messages/fanout.py`, `_FanOutRelay.step`):

    tau_b = tau_intron + tau_own(boundary)        lam_b = precision-weighted mean

*Assumes*: the intron's contained fragments and the boundary's unspliced crossings are two
**independent measurements of ONE composition**. The fragment sets are disjoint, so independence holds;
what is asserted beyond that is that the two slots have the SAME composition. Precisions ADD, so the
composed claim is **more confident than either input**.

⭐⭐ **AND THIS STEP IS HONEST — measured, not argued.** Scoring the composed `(lam_b, tau_b)` against
the BOUNDARY's own certified `true_f_g` gives median |z| ≈ **0.80** (L/R), i.e. an over-claim of only
**1.4–1.5×** capture-ON. ⛔ **None of the 44× exists before the transfer**, which localises the entire
defect to STEP 3 and retires "the compose adds precisions it should not" as a lead.

**STEP 3 — the transfer** (`flank_to_exon_lambda`), the delta method on three raw statistics:

    Var(lam_e) = a²/tau_b  +  w_mu²·(1/U + 1/S)
    a    = (1−f) + f·rho_nu/rho_r  = d(lam_e)/d(lam_b)          w_mu = rho_mu/rho_r

*Assumes*: the only uncertainty is SAMPLING NOISE in `lam_b`, `U` and `S`.

#### ①b ⭐⭐⭐ THE INVARIANT THE SHIPPED RELAY OBEYS AND THE FAN-OUT BREAKS

**Every precision transform in `RelayPolicy` moves precision in exactly one direction — DOWN.**
`_damp(p, s²) = 1/(1/p + s²)` (`relay.py:651`), `_damp_v(p, v) = p/(1+p·v)` (`relay.py:655`), `_dv`
(`relay.py:875`), `_dv_arr` (`relay.py:914`), and `mismatch_deflate` (`variance.py:531`) are all
`p ↦ p/(1 + p·v) ≤ p`. The only places precision RISES are additive fusions of genuinely INDEPENDENT
evidence — the two neighbours at `relay.py:1118`, the sj count at `relay.py:892`.

⛔⛔ **`flank_to_exon_lambda` is the message layer's ONLY precision-AMPLIFYING single-source transform**:
`tau_e = tau_b/a²` with `a → 0` gives `tau_e/tau_b` a p90 of ~210 capture-ON. ⭐ **So the fan-out did not
merely omit a deflation — it introduced an amplifier the shipped relay never contains, and then omitted
every damper the relay uses to keep such a thing safe.** Stated as a design invariant the backbone could
enforce: *a single-source transform may only REDUCE a precision; only independent evidence may add.*

#### ② WHY IT IS TOO HIGH — the failing assumptions, ranked

1. ⭐⭐⭐ **THERE IS NO PREMISE-VARIANCE TERM AT ALL.** Nothing in the chain prices "the flank's
   composition is only approximately the exon's", which §H.4 ④ measures at ~1 nat on the gDNA route
   capture-ON. ⛔ **`fanout.py` imports NOTHING from `messages/variance.py`** — not `transfer_logvar`,
   not `mismatch_deflate`, not `composition_logvar` — while `relay.py` and `currency.py` both do. It is
   the only SPEAKING policy with no premise-variance law. That is the defect stated structurally, and
   it is most of the over-claim. ⚠ **And `transfer_logvar` is NOT the missing law**: measured at this
   hop it is ~0.033 nats² (sd 0.18 nats) against a delivered composition error of median ~1.07 nats —
   **34× too small to touch it**, because it prices the SAMPLING uncertainty of the frame step and
   never its SIZE.
2. ⭐⭐ **THE DELTA METHOD IS EVALUATED AT A POINT THE SOURCE IS NOWHERE NEAR CONFIDENT ABOUT.**
   Since `rho_nu/rho_r = 1 − w_mu`, `a = (1−f) + f(1 − w_mu)`, so **`a → 0` iff `f → 1` AND `w_mu → 1`
   JOINTLY** — a gDNA-saturated flank whose RNA route is all sj flux (measured at `a < 0.1`: median
   `f` = 0.991, median `w_mu` = 0.998, both conditions always together). Then `tau_e = tau_b/a²`
   explodes: `a` reads p5 **0.0017**, p25 0.082, p50 0.461; `tau_e/tau_b` reads p50 3.25, p75 35.5,
   p90 **214**. Meanwhile `sd(lam_b) = 1/√tau_b` is 2–4 nats. §H.4 ④c prices the linearisation error at
   1.6–7× and shows it biases the MODE too.
3. ⭐⭐ **NOTHING CAPS THE CLAIM BY THE INFORMATION BEHIND IT.** The source flank holds a median 51
   unspliced crossings and 4 sj fragments; the destination exons carrying the error hold 10⁴. No term
   in the chain prevents ~55 fragments from out-weighing 10⁴.
4. ⚠ **A COORDINATE ASYMMETRY, and it is real rather than cosmetic.** The exon's own λ precision
   carries the Jacobian `[f_e(1−f_e)]²`, which crushes it as `f_e → 0` or `1`; the transferred `tau_e`
   carries no equivalent crush. So the same physical evidence buys the message far more weight than it
   buys the exon, purely because of where each sits in the λ coordinate.
5. ⚠ **THE FRAME CONSTANT `K` CARRIES COEFFICIENT 1 AND ZERO VARIANCE** (§H.4 ④b) — its four effective
   lengths are estimated, and their error is a bias that never confesses.

⛔ **Assumptions 1 and 5 are missing INFORMATION; 2 and 3 are missing CEILINGS; 4 is the coordinate.**
A dampener addresses 1 (and, incidentally, 5, because it reads the OUTCOME rather than the cause). It
does **not** address 3, and it interacts badly with 4 — see ④ below.

#### ③ THE REFERENCE PRIOR AT STAGE-3 AND STAGE-4 SLOTS — and 0.75 is DERIVED, not a magic number

`simplex_logodds.structural_reference_location` is four lines:

    a = b = _JEFFREYS_REF                                # ½
    m = min((a + 1)/(a + b + 1), sigmoid(L))             # = 0.75
    mature = mrna_active_pos | mrna_active_neg
    return where(mature, _NEUTRAL_LOCATION, m)           # ½ where mature, 0.75 elsewhere

⭐ **`0.75` is the posterior mean after ONE pseudo-observation of gDNA on the reference's own Jeffreys
exponents** — `(a+1)/(a+b+1)` at `a = b = ½` — so it is a composition mean derived from the reference's
own parameters, in the same units, introducing no new number, and it moves automatically if
`_JEFFREYS_REF` does. Its STRENGTH on the log-odds scale is `log(m/(1−m)) = log 3 = 1.0986` nats (3:1),
and the strength never scales with data. The swept alternative `m = σ(L)` (9.31 nats) was measured
**worst of all** at the REFUTE obligation, and the ladder cannot rank strength because it scores
DELIVER alone.

| stage-0 class | reference LOCATION | strength |
|---|---|---|
| `intergenic` REGION | 0.75 (¬mature) | 1.0986 nats |
| `ss_intron_region` | ⭐ the **MEASURED** background `ρ_bg·E_g / M`, clipped to the lattice cap (stage 2, `measured_reference_location`, default ON) | one pseudo-fragment |
| ⭐ `ss_intron_boundary` — **STAGE 3's destinations** | **0.75** — `measured_reference_location` is REGIONS ONLY (a measured refusal, 2026-08-26: the boundary-inclusive form regressed stranded × capture-ON on both substrates) | 1.0986 nats |
| ⭐ `solvable_exon` — **STAGE 4's destinations** | **NONE.** An exon is `mature`, so it gets `_NEUTRAL_LOCATION = ½` and **the reference term vanishes identically** | zero |

⭐⭐ **So stage 4's destinations receive no reference at all**, which is worth stating plainly because it
changes how the unstranded cell must be read.

**THE OWNER'S "LUCKY GUESS" HYPOTHESIS for unstranded × capture-OFF `g50` — REFUTED as stated, but the
intuition is right about the mechanism.** It cannot be the reference: exons get none. At κ = ½ the
strand λ-term is identically zero and the intron factory does not reach exons, so a claimed exon's own
λ evidence is **exactly 0** (§H.4 ⑦, measured at all 5,063 message-receiving exons) — and what solves
them is the **REFITTED gDNA prior**. §H.4 ⑦'s numbers are the proof: at `calib_refit_iters = 0` the
silent arm reads 268,498 and three refits take it to 6,415, a **42×** rescue, while the fan-out moves
only 2.7×. ⭐ So the owner is right that something other than the exon's own data is solving that cell
accurately and that the message displaces it — it is the refit, not the reference constant.

#### ④ WHAT THE SHIPPED DAMPENER WOULD DO — measured, as a DIAGNOSTIC arm

`messages.variance.mismatch_deflate` is the DerSimonian–Laird between-source mismatch variance, applied
by the relay on the λ axis at `relay.py:1029`:

    G   = lam_msg − lam_own                     (the observed composition gap)
    b̂²  = max(0, G² − v_msg − v_own)
    p   = 1/(v_msg + b̂²) = 1/max(v_msg, G² − v_own)

⭐ **Derived, not tuned — no magic number** — and its closed form states the safety property exactly:
*a message may out-weigh the destination's own belief only when it agrees with it to within
`√2·σ_own`.* Its three regimes fall out of `v_own` alone with no gate: `v_own = ∞` (no own evidence)
⇒ `b̂² = 0` ⇒ the message passes UNTOUCHED; `v_own` finite ⇒ agreement barely moves it and conflict
kills it; `v_own = 0` ⇒ the full `G²` is charged.

Deflating the STAGE-4 exon claim only, claimed-exon misplaced gDNA fragments:

| condition | silent | fanout | + deflate (exon only) | vs silent |
|---|---|---|---|---|
| `g00` ss.99 ON | 1,749 | 1,344 | **1,344 — untouched** | −23.1 % |
| `g05` ss.99 ON | 10,344 | 12,419 | **9,239** | **−10.7 % — a WIN** |
| `g50` ss.99 ON | 18,042 | 40,528 | **19,214** | **+6.5 %** (was +125 %) |
| `g98` ss.99 ON | 18,032 | 28,370 | 24,117 | +33.7 % (was +57 %) |
| `g50` ss.99 OFF | 5,460 | 5,935 | 5,632 | +3.2 % (was +8.7 %) |
| `g50` ss.50 OFF | 6,415 | 7,755 | **7,755 — untouched** | +20.9 % |
| `g00` ss.50 OFF | 192,486 | 21,836 | **21,836 — untouched** | −88.7 % |

⭐⭐⭐ **Both zero controls and the entire unstranded stratum are BYTE-IDENTICAL**, for the structural
reason the docstring gives: `τ_λ ≡ 0` at κ = ½ ⇒ `v_own = ∞` ⇒ `b̂² = 0`. The cells the fan-out wins
cannot be harmed by this term.

⛔ **Deflate the STAGE-4 claim ONLY.** Applying it to the stage-3 boundary claim as well costs boundary
wins (`g50` OFF boundaries go 20,886 → 23,480 against silence's 25,561; `g98` ON 52,127 → 59,187), and
that is principled rather than empirical: stage 3's premise is CERTIFIED by the stage-0 matrix 32/32,
stage 4's is the one that fails under capture.

⛔⛔ **IT IS NOT A COMPLETE FIX, AND THE LIMIT IS STRUCTURAL.** `g98` keeps most of its residual because
the DL test is RELATIVE: it fires on `G² − v_own`, and at `g98` the exon's own `f_g → 1` makes
`I_strand ∝ [f_g(1−f_g)]²` collapse, so `v_own → ∞` and `b̂² = 0` — **the dampener is disabled exactly
where failing assumption ④ (the coordinate crush) bites hardest.** The two failure modes interact, and
closing `g98` needs the missing INFORMATION CEILING (assumption 3), not a better disagreement test.

### §H.6 ⭐⭐⭐ THE REFERENCE, SWEPT ON THE CURRENT PANEL (2026-08-23) — a CONSTANT is the wrong shape

The owner's objection: `m = 0.75` at an intron|exon BOUNDARY is indefensible, because a non-exonic slot
has no annotated RNA and the only data-free prior there is 100 % gDNA — on the lattice, the highest
representable point `σ(L)`. ⛔ That form was swept and REFUTED once, but on a substrate holding
`nrna = 0`, and `simplex_logodds.py`'s own docstring names the falsification: *re-measure on a panel row
with `nrna > 0`.* **The sparse-nascent panel is that row, so the sweep is now runnable — and it was run.**

⚠ **The first run of this sweep was VOID and the reason is worth carrying**: `rigel.calibration.__init__`
rebinds the name `calibrate` to the FUNCTION, so `import rigel.calibration.calibrate as CAL` yields a
function and patching it is a silent no-op; and `calibrate.py` binds
`structural_reference_location` into its own namespace at import, so patching `simplex_logodds` alone
does nothing either. The arm reported "no effect at any m" — a perfect
`TRAPS: an-ablation-that-never-ran` signature. The fixed harness asserts its patch target.

**Claimed `ss_intron_boundary` slots, `Σ|est−true|` in fragments, `SilentPolicy` (so the reference is
isolated from the message), reference location swept at the ¬mature slots:**

| condition | certified true `f_g` | none (½) | 0.75 | 0.90 | 0.99 | grid `σ(L)` | argmin |
|---|---|---|---|---|---|---|---|
| `g00` ss.99 ON | 0.000 | **714** | 975 | 1,348 | 2,495 | 7,733 | **none** |
| `g00` ss.50 OFF | 0.000 | **48,363** | 100,297 | 136,537 | 155,919 | 157,662 | **none** |
| `g05` ss.99 ON | 0.662 | 27,529 | 22,714 | **22,363** | 26,282 | 38,557 | **0.90** |
| `g50` ss.50 OFF | 0.713 | 33,680 | 31,121 | **29,256** | 57,609 | 78,445 | **0.90** |
| `g50` ss.99 ON | 0.974 | 42,557 | 38,007 | 32,385 | **26,026** | 28,799 | **0.99** |
| `g98` ss.99 ON | 0.999 | 54,866 | 49,374 | 40,117 | 14,624 | **1,219** | **grid** |

⛔⛔⛔ **THE TABLE ABOVE IS CENSORED ON THE ENTIRE LOWER HALF OF THE LATTICE, AND THAT CHANGES THE
CONCLUSION.** Its floor is `m = ½` (zero nats), and at the zero controls the required location is far
BELOW it. Adding `σ(−L)` — the lattice's other end, "presume pure RNA" — to the same sweep:

| condition (silent, refits = 3) | `σ(−L)` | 0.11 | ½ | 0.75 | `σ(+L)` |
|---|---|---|---|---|---|
| `g00` ss.50 OFF | **5** | 4,294 | 48,363 | 100,297 | 157,655 |
| `g00` ss.99 ON | **2** | 345 | 714 | 975 | 4,526 |

⭐⭐⭐ **`σ(−L)` beats "assert nothing" by 9,673× and 357×.** So "no prior wins the zero controls" was an
artifact of the grid floor, and the corrected reading is far more interesting:

⭐⭐⭐ **THE OWNER'S FORM IS RIGHT AND ONLY ITS SIGN IS DATA-DEPENDENT.** A reference that sits at a
LATTICE EXTREME and waits to be overturned is the correct shape — it wins by three to four orders of
magnitude where it points the right way. What is wrong is that **which extreme** is right varies: pure
RNA at `g00`, pure gDNA at `g98`. ⛔ **And the shipped `structural_reference_location` can express
NEITHER end — it emits only ½ or `min(0.75, σ(L))`.** That is the defect, stated exactly: not a badly
chosen constant, but a function whose range excludes both of the answers the panel actually wants.

⚠ **TWO CORRECTIONS to how this section first read, both from re-derivation on all 16 conditions.**
① *"The argmin tracks certified truth row for row"* is TOO STRONG: argmin and the mass-weighted truth
disagree by 0.6–5.6 nats on the interior rungs and **flip sign at `g05` capture-OFF** (mass-median
0.0312 = −3.43 nats against an argmin of 0.90 = +2.20 nats). They corroborate in coarse direction only,
and `_location_term`'s PROPERTY 4 licenses "m is the prior's median at ONE slot", never a
fragment-weighted population median. ② Pooled over all 16 conditions the best single constant is
`σ(L)` at 3.22× the per-condition optimum against the shipped 0.75's **17.57×** — so on a pooled score
the owner's proposal is ~5.5× better than what ships. ⛔ Pooling across strata is exactly what this repo
forbids, so that is a caveat and not a verdict — but it does mean "0.75 is safer" is false in every
reading tried.

⭐ **The premise needs one refinement too.** *"A non-exonic boundary has no annotated RNA"* is true, but
it does not imply the slot is mostly gDNA: the unspliced population there is gDNA **+ nascent**, and at a
5 %-gDNA library nascent is the majority of it. Certified mass-weighted `true_f_g` at
`ss_intron_boundary`, `g05` capture-OFF: **0.1145**, with only **5.1 % of the mass** exactly pure gDNA.
So at low library gDNA the honest data-free statement at an intron boundary is nearly the opposite of
"100 % gDNA" — which is precisely why the location must be derived rather than asserted.

**A constant of any value is the wrong SHAPE, and the range must span BOTH ends.**

⭐ The derivation already exists and already ships — `density_deconv.measured_reference_location`,
`m_i = ρ_bg·E_g,i / M_i` from the intergenic pool — and it is **REGIONS ONLY** by a measured refusal
(2026-08-26): the boundary-inclusive form regressed stranded × capture-ON on both substrates, because
*"under capture the intergenic rate is the wrong rate for a probe-adjacent crossing."* ⭐⭐ **§H.5 ⓪ now
explains that refusal exactly**: every probe base is in an exon and introns carry zero, so the
intergenic pool is the UNPROBED rate and a probe-adjacent boundary is not on it. The refusal was right
and its cause is now named.

#### THE SECOND READING — a stronger boundary reference PARTLY CANCELS the transfer's bias

The same sweep under `FanOutPolicy`, at the claimed EXONS (whose own reference is identically zero, so
this is entirely the message moving):

| condition | none | 0.75 | 0.90 | 0.99 | grid | silence |
|---|---|---|---|---|---|---|
| `g05` ss.99 ON | 13,544 | 12,419 | 11,496 | **11,207** | 11,310 | 10,344 |
| `g50` ss.99 ON | 41,218 | 40,528 | 39,415 | 35,555 | **34,061** | 18,042 |
| `g98` ss.99 ON | 28,851 | 28,370 | 27,505 | 25,459 | **24,708** | 18,032 |

⭐ **Pushing the boundary reference toward pure gDNA monotonically improves the fan-out's exon delivery
under capture** — it pushes the source back against the probe cliff's RNA-ward pull. ⛔ **But it never
approaches silence** (34,061 against 18,042), so the reference is not a substitute for the dampener; the
two address different halves. ⚠ And it does nothing unstranded, structurally: at κ = ½ a boundary's own
`tau_lam` is 0, so its reference never reaches the composed claim.

#### THE STAGE-4 EXON REFERENCE — the owner's diagnosis, measured exactly

Certified gDNA densities, and what the intergenic-transported reference would assert at each exon:

| condition | exon / intergenic gDNA density | reference asserts | certified true | error |
|---|---|---|---|---|
| `g50` ss.50 OFF | **1.003×** | m = 0.0672 | 0.0691 | 5,998 fragments on 579,207 |
| `g50` ss.99 OFF | 1.013× | m = 0.0672 | 0.0696 | 6,295 on 580,653 |
| `g05` ss.99 ON | 586× | m = 0.0001 | 0.0593 | 69,676 |
| `g50` ss.99 ON | 516× | m = 0.0010 | 0.5455 | 698,490 |
| `g98` ss.99 ON | 497× | m = 0.0019 | 0.9834 | 1,369,418 |

⭐⭐ **The owner's reasoning is confirmed to the digit.** Capture-OFF the intergenic background
transported through the exon's own opportunity is a near-exact reference — a few percent relative on the
composition, with no solver anywhere in it. Capture-ON it is destroyed, asserting ~0 gDNA where truth is
0.55–0.98, because the intergenic pool is the one place the probes never reach.

**SOLVED, not just predicted.** Widening the SHIPPED estimator's selector from `ss_intron_region` to
`solvable_exon` (`dataclasses.replace`, so no arithmetic is re-implemented), scored at the claimed exons
with the stage-4 λ claim muted so the reference is the only thing varied:

* ⭐ **wins 10 of 16 conditions** — **all four `g00` rungs to EXACTLY 0 fragments** (182,112 / 61,040 /
  4,797 / 1,735 → 0, because `ρ_bg = 0` at a zero-gDNA library makes the transport trivially exact),
  and **all six contaminated capture-OFF rungs by 2.9–20.7 %**;
* ⛔ **loses all six contaminated capture-ON rungs**, by +0.03 % to +571 %. The damage is
  REFUTE-dominated (74.5 % of the added fragments) rather than concentrated on the DELIVER side, and it
  is not a lattice-clipping artefact.

⭐⭐⭐ **So the owner's proposal (A) is measured CORRECT wherever capture is off — including both zero
controls, where it is exact — and needs a capture gate for the six ON rungs.** That is a much better
starting point than the current NOTHING, and the gate is the only open question. ⚠ The enrichment
detector is trusted as a BOOLEAN only (`ROADMAP` rank 3), which is exactly the amount of information
this gate needs.

### §H.7 ⭐⭐⭐ STAGE 5 — THE MISMATCH DEFLATION, LANDED (2026-08-23, owner-approved)

**What changed, in one sentence:** `FanOutPolicy` now deflates its stage-4 exon λ claim by the shipped
DerSimonian–Laird mismatch variance against the destination's own λ solve — the relay's own law, on the
same axis, applied to the one hop that lacked any premise variance at all.

⭐ **The whole change is one guarded block in `messages/fanout.py` plus one import.** No new file, no
config flag, no constant: `message_policy` still defaults to `"relay"`, so **no shipped default moves.**

**Why it is admissible, restated as the three things that make it not a knob:**

1. **It is DERIVED** — DerSimonian–Laird between-source heterogeneity, `b̂² = max(0, G² − v_msg − v_own)`
   at the OBSERVED gap `G = lam_msg − lam_own`. No constant is introduced anywhere.
2. **It is the SHIPPED law** — `mismatch_var` defaults True in `RelaySwitches`, so this is production
   0.7.1 code applied at a second hop, not a new mechanism.
3. **Its safety is ALGEBRAIC, not gated** — where the exon has no own composition evidence
   (`tau_lam = 0`, which at κ = ½ is EVERY exon) `v_own = ∞` ⇒ `b̂² = 0` ⇒ the claim passes untouched.

**MEASURED, ladder, claimed slots, shipped refits (`pass0_claimed_ab.py`):**

| condition | silent | fanout BEFORE | fanout AFTER | |
|---|---|---|---|---|
| `g05` ss.99 ON | 10,344 | 12,419 | **9,239** | +20 % → **−10.7 %, a WIN** |
| `g50` ss.99 ON | 18,042 | 40,528 | **19,214** | +125 % → **+6.5 %** |
| `g98` ss.99 ON | 18,032 | 28,370 | **24,117** | +57 % → +33.7 % |
| `g50` ss.99 OFF | 5,460 | 5,935 | **5,632** | +8.7 % → +3.2 % |
| `g00` ss.99 ON | 1,749 | 1,344 | **1,344** | *byte-identical* |
| `g00` ss.50 OFF | 192,486 | 21,836 | **21,836** | *byte-identical* |
| `g05` ss.50 OFF | 1,922 | 2,024 | **2,024** | *byte-identical* |
| `g50` ss.50 OFF | 6,415 | 7,755 | **7,755** | *byte-identical* |
| `g98` ss.50 OFF | 7,983 | 6,193 | **6,193** | *byte-identical* |
| ⭐ `g50` ss.50 ON (deferred) | 697,694 | 169,659 | **169,659** | *byte-identical* |

⭐⭐⭐ **Every unstranded condition and both zero controls are byte-identical — including the DEFERRED
stratum, this policy's design target and largest win, unmoved to the fragment.** Boundaries move only
through the refit (`g50` ON 40,864 → 40,962), and the stranded capture-OFF boundary win is intact
(20,857 → 20,886 against silence's 25,561).

⚠ **Test chromosome** (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank` — naming the substrate): the
fan-out beats silence at every condition checked and on both populations, `g50 ss.99 ON` exons
8,226 → 1,712, the deferred stratum 12,799 → 5,334, the zero control 17,055 → 7,404.

**The gates** — six cases added to `tests/calibration/test_fanout_policy.py`, TWO of them verified
failing before the code existed and four constraining it afterwards; suite 3,720 → **3,726**. The
perturbation sweep caught **5 of 7** injected defects, and ⭐ **both misses are provably SEMANTICALLY
INERT rather than coverage gaps**: `mismatch_deflate` forms `b̂²` from `g*g`, so flipping the gap's sign
is the same function; and `known = isfinite(vo)` makes `v_own = 1/_EPS` and `v_own = ∞` identical
outcomes for any `|G| ≤ 2L`, which the λ clip guarantees. Recorded rather than papered over.

⛔ **WHAT THIS DOES NOT FIX, and it is structural.** `g98` keeps most of its residual: the DL test is
RELATIVE, so it cannot fire where the destination's own λ precision has itself collapsed
(`I_strand ∝ [f_g(1−f_g)]²` vanishes as `f_g → 1`). Closing that needs the information ceiling of
§H.5 ②-3, whose blunt form (`tau_e ≤ tau_b`) is measured-DISQUALIFIED at the deferred stratum (§H.5 ⑤).

### §H.8 ⭐⭐⭐ THE DERIVED REFERENCE — THE PLAN (2026-08-23). It is ROADMAP rank 3's rung ④, not a new track

⛔ **Read this first: do NOT build a second reference mechanism.** ROADMAP rank 3 already states rung
④'s door as *"ψ's location term widening from one scalar to `(m_lo, m_hi, w)`"* behind a
`composition_reference` flag. §H.6 is the REQUIREMENTS SPEC for that widening, measured. The work below
is rung ④ scoped and staged; nothing here needs new machinery invented.

**THE ESTIMAND, in one line.** ψ's location wants the slot's own expected composition *before any
solve* — "the gDNA I should expect HERE, over what I actually observe HERE" — asserted at one
pseudo-fragment so any real evidence overturns it:

    m_i  =  rho_bg(i) · E_g,i / M_i        clipped into the lattice, strength a + b = 1

⭐ **That form already ships** as `density_deconv.measured_reference_location` at ss-intron REGIONs. It
is not a new idea; the only open question is what `rho_bg(i)` may be at the OTHER claimed classes.

**WHAT IS MEASURED ABOUT IT (§H.6), which is what makes this a plan rather than a proposal:**

| where | with `rho_bg` = the intergenic pool | verdict |
|---|---|---|
| ss-intron REGION | shipped, default ON | ✅ landed |
| `solvable_exon`, capture-OFF | wins **all six** contaminated rungs by 2.9–20.7 % | ✅ take it |
| `solvable_exon`, **both zero controls** | **exactly 0 fragments** (`rho_bg = 0` ⇒ the transport is trivially exact) | ✅ take it |
| `solvable_exon`, capture-ON | loses **all six** contaminated rungs, +0.03 % … +571 % | ⛔ needs an on-rate anchor |
| `ss_intron_boundary` | 2.8–11.5× better than the constant capture-OFF, exact at `g00`; 2.5–4.3× WORSE capture-ON | ⛔ same |

⭐⭐⭐ **The single cause of every capture-ON failure is named and structural: `rho_bg` is pooled over
INTERGENIC, and intergenic carries EXACTLY ZERO probe bases (§H.5 ⓪). It is the UNPROBED rate.** A
destination that is probe-covered is simply not on it — measured, exonic gDNA density runs ~500× the
intergenic rate under capture and 1.003× off it. That is also the exact cause of the 2026-08-26
regions-only refusal, which said so in words (*"the intergenic rate is the wrong rate for a
probe-adjacent crossing"*) without yet knowing why.

**THE SIX REQUIREMENTS §H.6 IMPOSES ON RUNG ④** — each one is a measurement, not a preference:

1. ⭐⭐ **THE RANGE MUST SPAN BOTH LATTICE ENDS.** `structural_reference_location`'s range is
   `{½, min(0.75, σ(L))}`; the panel wants `σ(−L)` at `g00` (which beats asserting nothing by
   9,673×/357×) and `σ(+L)` at `g98`. **A function that cannot express the answer cannot be tuned into
   it** — this is the defect, more than any particular constant.
2. **`rho_bg` MUST BE ON-RATE FOR THE DESTINATION'S ENRICHMENT CLASS.** One library-wide rate cannot
   serve a probed and an unprobed slot; that is what `(rho_0, w)` and the enrichment DETECTOR are for.
3. **EXACT AT `g00`.** The transported form already is, structurally (`rho_bg = 0`). Any replacement
   must keep that — it is the cheapest possible falsification and it is free.
4. **STRENGTH STAYS ONE PSEUDO-FRAGMENT.** `_location_term`'s PROPERTY 2 — the tails are `e^(−|λ|/2)`
   for every `m`, only the LOCATION moves — is what makes widening the range safe. ⚠ But see §H.5 ③'s
   recorded risk: the FIRM clip already puts a large minority of intron REGIONs at ≈20.7 nats of
   effective range, ~19× the structural reference, and that has never been priced.
5. **REFERENCE AND MESSAGE INTERACT — PRICE THEM JOINTLY.** A stronger boundary reference monotonically
   improves the fan-out's exon delivery under capture (§H.6), because it pushes the source back against
   the probe cliff. Pricing either alone will mis-attribute the other.
6. **BOTH ZERO CONTROLS ON EVERY ARM; DELIVER/REFUTE NEVER POOLED.** The whole reference question is a
   DELIVER-versus-REFUTE trade and a pooled score picks the wrong end of it — that is precisely how
   `m = σ(L)` was once selected and once rejected.

**THE STAGING, cheapest-first, each step its own experiment:**

* **④a — the EXON selector, capture-gated.** Widen the shipped estimator from `ss_intron_region` to
  `solvable_exon`, gated OFF under capture. ⭐ Measured to win 10 of 16 conditions already; needs only
  the gate, and the enrichment detector is trusted as exactly the BOOLEAN this requires. **This is the
  cheapest real win on the board and it needs no new estimator.** ⛔ Its falsification is requirement 3:
  both zero controls must read exactly 0.
* **④b — the ON-RATE ANCHOR.** The open research question: what plays `rho_bg`'s role for a probed
  destination when probe coordinates are unavailable. This is rank 3's `(rho_0, w)` on the
  `AbundanceLandscape`, and §H.5 ⓪'s probe-footprint result is the physical statement it must reproduce.
  ⚠ Until it exists, capture-ON keeps the structural constant at boundaries — which is worse than
  nothing at `g98` and better than nothing at `g05`, so the gate matters more than the constant.
* **④c — THE RANGE WIDENING** (requirement 1), which is a one-function change to
  `structural_reference_location` but must land WITH a licence saying which end, or it is just a new
  constant with a bigger blast radius.

⛔ **What NOT to do, measured-refused:** replace 0.75 with any other CONSTANT — `σ(L)` is 3.3×/10.8×
worse than nothing at the zero controls, and 0.75 is optimal at no condition on the panel. The choice is
not between two constants; it is between a constant and a derivation.

### §H.9 ⛔ ④a — BUILT, PRICED ON ALL 16, AND DELETED (2026-08-24). The MEASUREMENT stands; the CODE is gone

⛔⛔ **READ THIS HEADER BEFORE THE SECTION.** `exon_measured_reference` was implemented, gated (6 cases,
5 verified failing first), priced on all 16 scenarios — and then **removed from `src/`, `config.py`
and the test suite**, because §H.10 and the external review settled that it is an instance of the
circular empirical Bayes the whole reference concept is refuted for: *fitting a prior on the data being
deconvolved*. ⭐ The table below is kept because it is evidence about the PROBLEM, not about a
mechanism — it is the cleanest demonstration that a data-derived reference wins wherever the background
rate is on-rate and fails catastrophically wherever it is not. ⛔ **Do not resurrect it as a flag**
(`no legacy, converge and delete`).

**Original section follows, unedited.**

#### ④a as it was built — `exon_measured_reference`, default OFF, and the 16-scenario ladder

**What changed:** `density_deconv.measured_reference_location` gained `include_exons`, widening its
SELECTOR from `ss_intron_region` to `ss_intron_region | solvable_exon`. Same `rho_bg`, same ratio, same
FIRM clip — **only the selector moves.** Wired to `CalibrationConfig.exon_measured_reference`,
**default `False` ⇒ bit-identical**.

**Gates:** five cases in `test_measured_reference.py`, all five verified failing first, plus ⭐ **a
WIRING gate in `test_calibrate.py`** that the perturbation sweep demanded: with `calibrate` hard-wiring
`include_exons=False`, every estimator-level gate still passed — a flag that nothing forwards is a flag
that does nothing. Sweep now **6/6**. Suite 3,726 → **3,732**.

⚠ **The import trap bit for the THIRD time this session and is now written into the gate itself**:
`rigel.calibration.__init__` rebinds the name `calibrate` to the FUNCTION, so
`import rigel.calibration.calibrate as m` hands back a function and patching it is a silent no-op; and
`calibrate.py` binds the estimator into its own namespace at import, so patching `density_deconv`
misses too. Reach the MODULE through `sys.modules`.

#### THE 16-SCENARIO LADDER — claimed EXONS, misplaced gDNA fragments, never pooled

| scenario | certified gDNA | silent | relay | fanout | fanout+④a | silent+④a |
|---|---|---|---|---|---|---|
| `g00` ss.50 OFF | 0 | 192,486 | 1,548 | 21,836 | **0** | **0** |
| `g00` ss.50 ON | 0 | 62,180 | 5,052 | 39,971 | **0** | **0** |
| `g00` ss.99 OFF | 0 | 4,935 | 1,040 | 1,234 | **0** | **0** |
| `g00` ss.99 ON | 0 | 1,749 | 1,181 | 1,344 | **0** | **0** |
| `g05` ss.50 OFF | 4,077 | 1,922 | 2,072 | 2,024 | 1,891 | **1,816** |
| `g05` ss.50 ON | 69,671 | 48,722 | 28,517 | **31,266** | 66,783 | 69,614 |
| `g05` ss.99 OFF | 4,038 | 1,761 | 1,912 | 1,809 | 1,737 | **1,711** |
| `g05` ss.99 ON | 69,751 | 10,344 | 11,109 | **9,239** | 50,525 | 69,686 |
| `g50` ss.50 OFF | 39,997 | 6,415 | 7,239 | 7,755 | 7,284 | **5,930** |
| `g50` ss.50 ON | 699,284 | 697,694 | 258,207 | **169,659** | 538,773 | 698,160 |
| `g50` ss.99 OFF | 40,409 | 5,460 | 5,941 | 5,632 | 5,362 | **5,232** |
| `g50` ss.99 ON | 699,689 | **18,042** | 24,631 | 19,214 | 20,845 | 20,982 |
| `g98` ss.50 OFF | 79,297 | 7,983 | 8,131 | 6,193 | **5,295** | 6,335 |
| `g98` ss.50 ON | 1,371,478 | 1,368,599 | 370,604 | **247,154** | 967,928 | 1,368,978 |
| `g98` ss.99 OFF | 79,519 | 6,346 | 7,357 | 6,115 | **5,281** | 5,371 |
| `g98` ss.99 ON | 1,372,029 | **18,032** | 35,706 | 24,117 | 25,479 | 22,292 |

⭐⭐⭐ **④a wins 10 of 16 and the split is exactly the predicted one**: **all 8 capture-OFF rungs**, and
**all 4 `g00` rows to EXACTLY ZERO** — including the capture-ON ones, because that is the `rho_bg = 0`
branch and it does not care about capture. It **loses all 6 contaminated capture-ON rungs**.

⛔⛔ **AND THE LOSS IS NOT UNIFORM — it is concentrated where the exon cannot argue back.** At
UNSTRANDED capture-ON (`g05` 48,722 → 69,614; `g50` fanout 169,659 → 538,773; `g98` fanout 247,154 →
967,928) the damage is catastrophic, because `τ_λ ≡ 0` at κ = ½ leaves the exon no own evidence to
overturn a reference asserting near-zero gDNA. At STRANDED capture-ON it is modest (`g50` 18,042 →
20,982; `g98` 18,032 → 22,292) — the strand solve overturns it. ⚠ `g05` ss.99 ON is the in-scope
exception and it is severe (10,344 → 69,686), because at 5 % gDNA the reference's error is the whole
available mass.

⭐ Boundaries are essentially untouched (④a only selects exons; the small movements are the refit).

⛔ **So default OFF stands, and a CAPTURE GATE would convert this into a clean win** — 8/8 capture-OFF
plus the zero controls, inert elsewhere. ⚠ **A fully-DERIVED subset is available with no detector at
all**: gate the widening on `rho_bg == 0`, which is exactly the four `g00` rows and is measurable
without any solve. ⛔ But `ROADMAP` rank 2 rules that **`g00` is ONE-SIDED — a win confined to it is
not a win** — so that subset is a safety property, not a result, and must not be sold as one.

#### 0.8.0's OWN METRIC, all 16, ranked, scope-tagged (`calibration_vs_oracle.py`, SHIPPED tree)

Worst IN-SCOPE is **`g98 ss.99 capture_on` at Σ|Δ| = 436,569 fragments**, then `g50 ss.99 ON` at
388,449; the two worst rows overall are the DEFERRED stratum (3,191,561 and 2,190,807) as always.
⭐⭐ **The residual is overwhelmingly UNDER-calling gDNA at the contaminated rows** — `g98 ss.99 ON`
reads over 8,359 against under 428,210.

⭐⭐⭐ **AND THE AVAILABILITY CENSUS IS WHAT MAKES THAT ACTIONABLE** (no solver, straight off certified
truth): **the four `g98` rows hold 64.06 % of the panel's 71.3 M certified gDNA fragments, `g98` + `g50`
together hold 96.73 %, and all eight `g00`/`g05` rows hold 3.27 %** — the four `g00` rows exactly zero.
⛔ So an error at `g05` cannot matter much however large its ratio looks, and **effort belongs at `g98`,
where both the worst in-scope error and two thirds of the panel's mass sit.**

### §H.10 ⭐⭐⭐ "THE REFERENCE IS TOO STRONG" — the DIAGNOSIS confirmed, the PRESCRIPTION refuted (2026-08-24)

**The owner's position (2026-08-24):** ψ's reference is not a prior at all — a prior is what you know
before the data, and pass-0's whole *job* is to learn it, so setting it FROM the data is cheating by
definition. All we truly know a priori is STRUCTURE. The term exists because the log-odds simplex needs
a numerical stabilizer. And its nominal "one pseudo-count" understates its real pull. **Hypothesis:
weaken it to a pure stabilizer.**

#### ① THE DIAGNOSIS IS RIGHT, AND THE CODEBASE ALREADY MEASURED IT

`simplex_logodds.py:660-663` records it verbatim and it was never acted on: *"the real overturn depth
is **996–1,993 fragments** — the prior is strongest precisely where the likelihood is flattest.
⛔⛔ On UNSTRANDED data (κ = ½) **it is never overturned at any depth** (measured 0.9998 against a
pure-RNA truth at 10,000 fragments)."*

Re-measured directly on ψ (one single-strand slot, TRUE `f_g = 0`, reference `m = 0.75`, noiseless
strand split — so the answer is the reference against the likelihood and nothing else):

| κ | N = 10 | N = 100 | N = 1,000 | N = 10,000 | N = 10⁶ | overturned at |
|---|---|---|---|---|---|---|
| 0.99 | 0.1203 | 0.0298 | 0.0092 | 0.0030 | 0.0004 | **N = 3** |
| 0.90 | 0.1836 | 0.0393 | 0.0118 | 0.0038 | 0.0005 | **N = 3** |
| 0.75 | 0.4414 | 0.0725 | 0.0201 | 0.0064 | 0.0008 | **N = 10** |
| **0.50** | **0.7471** | **0.7471** | **0.7471** | **0.7471** | **0.7471** | ⛔ **NEVER** |

⭐⭐⭐ **Read the last row: at κ = ½ the posterior is 0.7471 at every depth from 10 to 1,000,000 — the
answer does not depend on the data AT ALL.** And scaling the location's weight `c` moves it and nothing
else: `c = 1 → 0.7471`, `0.3 → 0.5806`, `0.1 → 0.5271`, `0 → 0.5000`, **identical at N = 100, 1,000 and
10,000.**

⭐ **So "one pseudo-fragment" is honest as a PSEUDO-COUNT and false as LEVERAGE.** Where the strand
channel is alive the term is worth ~one fragment exactly as advertised (κ = 0.99: `c = 1` gives 0.0298
against `c = 0`'s 0.0284 at N = 100 — a 5 % difference, and none at all by N = 1,000). Where the strand
channel is DEAD its leverage is unbounded, because leverage is the reference's nats against the data's
own Fisher information on λ, and that information is `N_eff·disc·[f_g(1−f_g)]²` — **identically zero at
κ = ½**. Divide by zero. ⭐⭐ **This is the SAME coordinate asymmetry that made the fan-out's message
over-strong (§H.5 ② ④): one mechanism, two symptoms.**
⭐ And `c = 0` at κ = ½ returns exactly **0.5000** — proper, finite, no pole. **The Jeffreys Beta(½,½)
alone IS a sufficient stabilizer**, so the location term is not needed for the numerical job.

#### ② THE PRESCRIPTION IS REFUTED ON THE PANEL — weakening makes every contaminated rung WORSE

Location weight swept on all 16 conditions, both policies, claimed populations (`c = 0` is
"Jeffreys stabilizer only"). Selected rows, and the pattern is uniform:

| condition | pol | B c=1 | B c=0.3 | B c=0.1 | B c=0 |
|---|---|---|---|---|---|
| `g00` ss.50 OFF | silent | 100,297 | 66,840 | 54,478 | **48,346** |
| `g00` ss.50 ON | silent | 22,738 | 9,333 | 7,115 | **6,228** |
| `g05` ss.99 ON | silent | **22,714** | 25,289 | 26,454 | 27,131 |
| `g50` ss.99 OFF | fanout | **20,886** | 24,504 | 25,245 | 25,505 |
| `g50` ss.50 ON | fanout | **251,503** | 530,581 | 656,091 | 712,437 |
| `g98` ss.50 OFF | fanout | **20,020** | 34,544 | 38,057 | 39,574 |
| `g98` ss.99 ON | fanout | **52,127** | 68,361 | 72,139 | 73,994 |

⛔⛔ **At the DEFERRED stratum — the fan-out's biggest win and exactly where the owner expected the
reference to be fighting — removing it is 2.0× worse on exons (169,659 → 346,513) and 2.8× worse on
boundaries (251,503 → 712,437).** The reference is HELPING there, not fighting.

⭐ Weakening helps only at the ZERO CONTROLS (2.1×–3.7×), which is one-sided by `ROADMAP` rank 2 and
therefore not a win.

#### ③ THE SYNTHESIS — and it changes what to build, not whether

Both halves are true and they are not in tension:

* Where the strand channel is dead, **the reference is not a weak prior being out-voted — it IS the
  entire answer, at any depth.** That is exactly as unprincipled as the owner says.
* But on this panel it happens to point the RIGHT WAY at intron boundaries (certified truth there runs
  0.66–0.999), so deleting it removes the only voice and the error doubles.

⛔ **So the fix is not to weaken the reference. It is to give those slots real INFORMATION so the
reference stops being the only voice** — which is precisely what the fan-out message does, and why it
wins the deferred stratum 4× (697,694 → 169,659) where no reference change comes close.

⭐⭐ **And this re-reads ④a's failure correctly.** ④a was catastrophic at unstranded capture-ON not
because "a reference is bad there" but because THAT reference **pointed the wrong way**: the
intergenic-transported location asserts ~0.001 gDNA where certified truth is 0.65–0.98. Same slot, same
strength, opposite sign, opposite outcome. ⭐ Which makes the on-rate anchor (§H.8 ④b) the thing that
matters, and weakening a knob the thing that does not.

## §G OPEN DECISIONS — ALL FOUR RULED (owner, 2026-08-26)

* **`the-one-hop` — RULED: it IS a message.** The exon's hop uses the message framework's operator when
  the exon-solve stage lands (a constraint recorded on the rebuilt `RelayPolicy`; structural termination
  becomes a licence state). Stages 0–2 write NO exon-transport code: stage 0 only ADMITS solvable exons
  to the substrate, whose claims are checkable against certified truth with no solve.
* **`substrate-home` — RULED: derived at calibration init.** One named layer-3 module, one gated
  definition, O(n_slots) array math. Persisting to the index is refused for now: per-BOUNDARY structural
  flags were deliberately dropped at `INDEX_FORMAT_VERSION 7`, and a persisted calibration-facing
  artifact would need a reach digest no existing hash provides.
* **`clamp-semantics` — RULED: FIRM.** At a collision (`ρ_bg ≥ ρ_total`) the prior's LOCATION clips to
  the lattice cap (the valid CAP use of `EQUATIONS.md` §9c.1, never a chooser) and the STRENGTH stays the
  standing one-pseudo-fragment budget — real evidence overturns it. Grounds: the collision population is
  dominated by truth-exactly-1 slots under sparsity, and firm beat soft on every stratum when priced
  (§9c.1). The stage-2 score must split DELIVER/REFUTE and never pool them; `g00` is the structural
  falsification (no collision may occur there).
* **`combine-form` — RULED: A, the lattice product** (owner's integrate-don't-invent rule): the existing
  grid combine is the solve, and each channel's own point estimate + precision is published as a
  diagnostic. D.1's independence check becomes a test (each channel consumes its own raw sufficient
  statistic).
