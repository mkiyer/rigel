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

**Implementation shape:** the new message POLICY on the existing backbone (the `the-one-hop` ruling) —
emission licensed only FROM solved-substrate slots TOWARD adjacent unsolved ones, `SilentPolicy`
behaviour elsewhere; the forward scan carries the rightward fan-out, the backward the leftward, so the
two existing scans cover the depth-2 pattern with no backbone change. Stage 3 first; falsification
test first; test chromosome before the ladder; both zero controls on every rung.

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
