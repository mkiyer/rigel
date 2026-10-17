# Spliced-fragment artifacts: the derivation (2026-10-14, after the external review)

> **Review status, 2026-10-15:** this is the earlier derivation, not an implementation contract. The owner's
> complete-pair and observable-opportunity requirements now govern. `INTERGENIC_GDNA_REVIEW.md` demonstrates an
> existing orphan-placement bug and length-contrast collapse. Its normalized derivation and the revised
> `SPLICE_ARTIFACTS_IMPLEMENTATION.md` supersede this note where they conflict: unchanged opportunity in §4,
> automatic independence after deposit/move in §7, and unconditional allele identifiability/capacity/prior claims
> in §9 are **not established**. The second-pass strand omission described below was fixed in `41899ae4`.

This follows `SPLICE_ARTIFACTS_DESIGN.md` §12 and the review (`SPLICE_ARTIFACTS_REVIEW.md`). The review's verdict:
build the competing-placement architecture, but not as first specified. Here the model is written out term by term
before any code. Three mechanisms are covered:

- **Local contiguity**, derived fully in §1–§8.
- **The sample's own junctions** (deletion alleles), sketched in §9.
- **Remote origin** (retrocopies, paralogs), sketched in §10.

Owner, 2026-10-14: 0.8.0 waits for the mechanism and its validation on independent artifact-containing evidence.

## 1. The likelihood the EM already uses

A fragment's candidates are pairs: a component `k` (a transcript `t`, or the locus's gDNA `D`) and a physical path `P`
(the genomic intervals its sequenced bases came from). The EM's per-fragment weight is

```
w(k, P) = (θ_k / ℓ_k) · f_k(L_k(P)) · s_k(o) · e(x | P)
```

| term | meaning |
|---|---|
| `θ_k / ℓ_k` | abundance per unit of opportunity (unchanged) |
| `f_k(L)` | the RNA or gDNA fragment-length law, at the candidate's molecule length |
| `s_k(o)` | the strand likelihood of the observed orientation: κ or 1 − κ under RNA, ½ under gDNA |
| `e(x \| P)` | the read bases' likelihood along `P` |

The fragment's likelihood is the SUM over candidates: alternative origins add, they never multiply.

Today `e` is a common factor for all of a hit's candidates: the NM penalty, `log(mismatch_alpha)` per mismatch.
Every candidate shares the reported path, so the factor cancels. A candidate on a different path needs its own
`e`. That means the common factor plus a derived relative correction (§3), counted once.

## 2. The hypothesis space: the aligner proposes, it does not testify

The reported alignment is a deterministic function of the read's bases and the aligner's settings. It is a
**proposal of paths**, not evidence: STAR's score is never turned into a probability.

For a reported spliced alignment, the candidate paths are:

- **the reported spliced path,** read by the RNA candidates as today;
- **the contiguous paths,** one per anchor block, read by gDNA and also by RNA (an unspliced molecule is admissible
  RNA, Axiom 0). Every other block of the mate is laid contiguous with the anchor:
  - for `[s−a, s) | N | [e, e+b)` that gives right-anchored `[e−a, e+b)` and left-anchored `[s−a, s+b)`;
  - every sequenced base is placed; nothing is clipped;
  - a multi-junction mate lays all its other blocks contiguous with the anchor, so the gDNA path splices nowhere;
- **partial paths** (RNA only): one junction kept and another laid contiguous. This is the "one junction genuine, one
  artifact" fragment.

**Both mates** describe one molecule. A gDNA candidate is the combination of the two mates' contiguous paths that is
one interval in a valid orientation, with its length from the outer extent. Identical physical paths are deduplicated.

The aligner proposes no path it did not report. A remote origin, a retrocopy unspliced elsewhere, is therefore not
in this space. That is §10, a declared coverage gap until it is built.

## 3. The sequence emission

Under independent symmetric substitutions, a read base `x` against a reference base `r` has likelihood `1 − ε` if it
matches and `ε/3` if not. Two paths that place the same bases differently disagree only at the **discriminating
positions**, where their reference bases differ. Let:

- `m_s` = discriminating positions where the read matches the spliced reference;
- `m_c` = those where it matches the contiguous reference;
- bases that match neither contribute 1.

Then

```
e(x | P_contiguous) / e(x | P_spliced) = [ ε / (3(1 − ε)) ] ^ (m_s − m_c)
```

- **At h = 0** (no discriminating position) the ratio is exactly 1, whatever the read's errors. Nothing but the
  reference is needed.
- **At h ≥ 1** the read's own bases are needed: SEQ at the discriminating positions. NM cannot tell `m_s` from `m_c`
  (review, finding 7). Base qualities can replace ε per base later, as a separate A/B.
- **ε** is the aligned substitution divergence, Σ(NM − I − D) / Σ(M, =, X), over uniquely mapped unspliced fragments.
  It is read from the library, and it includes the sample's variants, which is acceptable for a ratio.

**Measured on VCaP** (`seq_prong.py`; single-junction units, short side ≤ 40 bp):

| | units | h = 0 (tie) | reads vote contiguous (`m_c > m_s`) | reads vote spliced |
|---|---|---|---|---|
| DNA (artifact), annotated | 2,626 | 2,124 | **0** | 502 |
| DNA (artifact), unannotated | 1,080 | 197 | 81 | 791 |
| RNA, annotated | 6,989,530 | 56,352 | 2 | 6,933,098 |

The annotated h ≥ 1 artifacts carry the SPLICED reference's bases. They are not mismatch-tolerant stitches: their
sequence really is exon-joined, which is the signature of a retrocopy (§10). Read-sequence evidence is required for
correctness, but it adds little on this library. False contiguous votes on the RNA side are 2 of 6.99 M.

## 4. Length and opportunity

- The contiguous gDNA candidate's length is its molecule's genomic extent, both mates included. gDNA's law is read
  there.
- Opportunity is unchanged: `ℓ_D` is the locus's gDNA effective length, and a path is one placement in it.
- A contiguous path that leaves the locus's span is out of scope for this mechanism: it lies in another locus.

## 5. Pruning

The existing rule applies to the full likelihood, `f · s · e` and the overhang term, inside the scorer. No candidate
is screened on a partial product: the review's counterexample (h = 2 times ½/κ on the wrong strand gives 4.4e-4,
above the 1e-4 bound) shows partial screens change the admitted set.

## 6. What the mechanism cannot do

Take a singleton at a likelihood tie: h = 0, the right strand, comparable lengths. It is decided only by the
calibrated share `q` and the effective-length ratio `B`.

- The MAP fixed point is `z = B(z+q) / (B(z+q) + 2 − z − q)`.
- With every prior unit on gDNA (`q = 1`) and `B = 0.05`, the read stays 0.947 on RNA.

No constant and no iteration creates information the observations lack. The fixtures report this matrix before and
after; the mechanism claims nothing at ties.

## 7. Calibration: the same candidates, its own decision

Calibration keeps its hard, seeded draw (ruling) over the same candidate paths and emissions. Two contract fixes come
first, and both are bugs today:

- **Unannotated junctions disagree between the stages.** Calibration banks an all-unannotated-intron fragment as
  unspliced; the EM keeps it RNA-only.
- **The second pass skips its strand term** when a motif was sequenced. A spliced flag alone would not restore the
  competition.

**When to decide.** The scan cannot know κ or the length laws, and `½/κ` has no bound as κ → 0. So no scan-time hold
rule is provably equal to the full rule, and holding every spliced fragment costs about 1–2 GB on VCaP. Two designs:

- **(a) Hold every spliced fragment** that has a contiguous candidate. Simple, but memory-heavy.
- **(b) Deposit as today, then move.** Record a compact side row per spliced fragment: its path ids, anchor, `m_s`,
  `m_c`, h, orientation and the contiguous extent. After pass one, draw each one's reading with the full likelihood.
  Reverse the spliced deposit and deposit the contiguous one for those drawn contiguous. Deposits are deterministic
  integers, so the reversal is exact and conservation is checkable.

**Recommended: (b).** It costs about 30 bytes per spliced fragment and adds an `undeposit` gate in the accumulator.
A held junction's support then comes from independent pass-one evidence, so no deposit of its own is read against it.

## 8. The training sets

The strand model (κ) and the RNA length law train on each spliced fragment weighted by its probability of being RNA
under a valid RNA geometry, not merely `1 − P(gDNA)`. The weighted fit statistics stay separate from calibration's
integer deposits.

## 9. The sample's own junctions (derivation draft, 2026-10-15)

**The physics.** A DNA read across a deletion or rearrangement in the sample's genome looks spliced and fits no
contiguous path in the reference. VCaP's unannotated h ≥ 1 artifacts fit nowhere 79 % of the time. Such a read is
still a gDNA molecule, drawn from the same gDNA component as its neighbours. What differs is only that the sample
carries an allele with this junction.

**The model.** Each unannotated junction `j` has an allele fraction `α_j ∈ [0, 1]`: the share of gDNA molecules at
`j`'s position that carry the junction. A fragment at `j` then gets a second gDNA placement, gDNA at the reported
spliced path, weighted

```
w(D_j) = (θ_D / ℓ_D) · f_D(L) · ½ · α_j
```

against its RNA readings as before. `α_j = 0` adds nothing. No new population appears: the allele is an observation
mechanism of gDNA (Axiom 0).

**Why `α_j` is identifiable without strand.**

- **Capacity.** gDNA can produce at most `α_j · ρ_D · W_j` junction reads. `ρ_D` is calibration's gDNA density at
  `j`'s ends, and `W_j` is the number of fragment starts whose fragment crosses the breakpoint with both anchors
  placed, an opportunity of the same kind calibration already computes.
- **The other allele.** The reference allele's reads run contiguously through `j`'s two ends. Together with the
  junction reads they fix `α_j`, as a variant caller's allele fraction would.
- **What follows.** A novel junction in an expressed gene carries far more reads than any `α_j ≤ 1` allows, so it is
  RNA. A singleton in a silent, gDNA-rich region fits `α_j > 0`. Between the two, strand decides on a stranded
  library, and on an unstranded one the posterior honestly stays split.

**The prior that `α_j > 0`** comes from the junction's class, as a likelihood ratio derived from the annotation, not
from any library:

- P(class | allele) is the chance that a random breakpoint lands on a known splice site, from the annotation's
  site density;
- P(class | RNA) is how often genuine junctions reuse known sites, from the annotation's own alternative splicing:
  how often an annotated junction shares an end with another;
- the base rate (how many junctions are alleles) is fitted jointly with the EM, so prevalence is counted once
  (review, finding 1).

**Fitting.** `α_j` is one more parameter of the locus EM. Its M-step is the expected allele junction reads over the
expected gDNA molecules spanning `j`'s ends, junction or reference. That is a binomial with the class prior as its
Beta-spike prior. No constant is introduced: the opportunity `W_j`, the density `ρ_D` and the class likelihood ratio
are all derived.

**Open:**

- `W_j` for a breakpoint whose ends lie in different regions, and for multi-junction reads;
- how a tumour's subclonal `α` (continuous) and a germline `α ∈ {½, 1}` share one prior;
- the annotation-withholding controls (review, finding 6): a withheld genuine junction must stay RNA, because its
  read count exceeds gDNA's capacity.

## 10. Remote origin (sketch): a junction-retrocopy table from the reference alone

The annotated h ≥ 1 artifacts carry the spliced sequence exactly (§3): genomic copies of the joined exons. If the
copy is identical across the read, the aligner meets a tie between the parent's spliced path and the copy's
contiguous one, and STAR's annotated-junction bonus picks the parent, so the read is reported unique.

**A reference-only prong.**

- At index time, search each annotated junction's spliced sequence contiguously in the genome. That is the last and
  first `L` bases of the joined exons, up to the read length; an exact or near-exact k-mer index suffices.
- Record the remote loci per junction.
- A spliced fragment at a listed junction then gets a remote gDNA candidate at the copy, with its emission against
  the copy and gDNA's density there.

This needs no `alignable`. It would cover the blacklist's unique catches (the review's 1,155 records) by evidence
rather than by a cut.

**Measured, 2026-10-14** (minimap2 `-x sr`, unspliced):

- **The artifact reads.** Of VCaP's DNA artifact reads at annotated junctions with h ≥ 1, 71 % fit the genome
  contiguously and near-perfectly (≥ 95 % covered, NM ≤ 2), and 996 of those 1,056 fit at ANOTHER locus. The h = 0
  artifacts fit locally (97 %). Unannotated h ≥ 1 artifacts fit almost nowhere (21 %): rearrangements, §9. The RNA
  control fits remotely 3 % of the time, from genes that have copies.
- **The table itself.** Each annotated junction's joined 100-mer was aligned to the genome. The alignment took 20 s
  and found a remote near-perfect copy for 2,658 of 404,168 junctions (0.66 %, median one copy).
- **Recall.** 533 of VCaP's 558 retrocopy-artifact fragments (95.5 %) sit on a listed junction.
- **Exposure.** 484,357 RNA units (3.2 %) sit on listed junctions. They get the copy as a candidate, decided by
  abundance; they are not cut.

**Index time against live realignment.** Realigning every spliced read unspliced at scan time would also find copies
of unannotated junctions: 75 of VCaP's 2,067 unannotated h ≥ 1 artifacts. It costs a genome index in memory at scan
time (minimap2's peaked at 15 GB) and a realignment of every spliced read. Every read across an annotated junction is
a window of that junction's joined sequence, so the index-time table answers the same question once per junction.
minimap2 is then needed only by `rigel index`. Open: the search's cost at index time, near-exact copies, and the multimapper intergenic hits
(`ISSUES: multimapper-intergenic-alignments`).

## 11. Mechanism fixtures, before any code

Each fixture is a falsification test, verified failing first, then the fixed code broken to watch it fire.

- an exact tie (h = 0), both anchoring directions, both mates;
- discriminating bases: the read matching the spliced reference, the contiguous one, or neither; equal NM with the
  mismatch moved on and off the discriminating position; an indel; a soft clip;
- an interior artifact junction, a two-junction mate, genuine plus artifact junctions (the partial RNA path kept);
- contradictory mates, missing XS, missing SEQ;
- a junction whose fragments are all held, so pass-one support is absent rather than zero;
- an omitted true annotation: a genuine novel isoform must not flow into gDNA;
- the singleton matrix (§6), under MAP and VBEM, from both starts;
- conservation: the deposits after the move equal a direct deposit of the drawn reading.

## 12. The order of work

Each step is A/B'd alone against its parent on the VCaP mix (both source libraries apart, both errors, the transcript
table beside them), plus harm checks on the ladder and both fragment-length gap panels.

1. The fixtures (§11).
2. The stage-contract fixes (§7), each its own A/B.
3. The contiguous candidates and the h = 0 emission, in the EM.
4. Calibration's move (§7 b).
5. SEQ at the discriminating positions, for h ≥ 1.
6. The training weights.
7. Then §10's table and §9's rate, each derived first.

The gate is the review's: independent artifact-containing evidence before release. That is the cluster.
