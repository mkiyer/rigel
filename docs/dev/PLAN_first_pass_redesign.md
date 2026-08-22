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

## §F THE TWO RUNNING LOGS (owner-requested, updated as the design lands)

### F.1 ⛔ OBSOLETE — superseded by this design

*(populated as each stage lands; each entry names the file, the symbol, what supersedes it, the evidence,
and what has to be re-priced before it can go)*

### F.2 ⭐ SIMPLIFICATION — smaller, clearer, fewer branches

*(same, for code that survives but gets simpler)*

---

## §G OPEN DECISIONS FOR THE OWNER

* **`the-one-hop`** — the exon's single hop: use the message framework's operator, or state why a structurally
  terminated hop is a different thing? (`TRAPS: one-hop-lifted-out-is-still-the-relay`)
* **`substrate-home`** — at index time and persisted, or derived at calibration init? (recommendation:
  derived)
* **`clamp-semantics`** — location only, or strength too?
* **`combine-form`** — lattice product with per-channel diagnostics, or two point estimates?
