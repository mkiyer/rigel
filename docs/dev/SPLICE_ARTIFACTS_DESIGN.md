# Spliced-fragment artifacts: the design and its evidence (review packet, 2026-10-14)

This packet is for an external reviewer and the project owner. It covers `ISSUES: splicing-artifacts`, which blocks
the 0.8.0 release.

**Where things stand.** The owner agrees with the proposal in §5, subject to this review. Nothing has been built. All
numbers come from today's tree (`origin/main` 7805ede7) and the measurement scripts listed in §11. This file is
sandbox: nothing in it is authoritative until it moves to `DESIGN.md`.

## 1. Rigel in one paragraph

Rigel is a Bayesian RNA-seq transcript quantifier that separates RNA from genomic-DNA (gDNA) contamination.

- **The scan.** A single-pass C++ BAM scanner resolves each fragment against the annotation and tallies it.
- **Calibration.** It deconvolves every genomic region and boundary into gDNA, RNA+ and RNA−, the three populations
  and no fourth, from strand, density and fragment length. It then hands each locus a gDNA prior.
- **The locus EM.** It assigns fragments to transcripts and to one gDNA component per locus. Each fragment's
  likelihood under each candidate combines:
  - the fragment-length law, at the candidate's own length;
  - strand: κ or 1 − κ under RNA, ½ under gDNA;
  - overhang and mismatch penalties.

**Why spliced fragments matter.** A fragment whose alignment crosses a splice junction is **certified RNA**: gDNA
cannot splice. So spliced fragments get no gDNA reading. Calibration treats them as RNA without deconvolution, and
they train the strand model (κ) and the RNA fragment-length law. They are the tool's most specific signal, so an
artifact among them is costly.

## 2. The defect

An aligner sometimes writes a junction on a molecule that was not spliced there. Rigel then certifies the fragment as
RNA, and trains κ and the RNA length law on it.

The defect is **two-sided**. Today's guard, a splice-junction blacklist (§3), also rejects genuine junction reads.
Those become unspliced, gDNA-eligible fragments whose footprint still spans the intron. The cost lands on the gDNA
and nascent pools of exactly the listed genes. A repair is judged on both errors together.

## 3. What the tree does today (read from the code)

- **The blacklist** comes from the companion tool `alignable`. It records the junctions an aligner wrote on simulated
  genomic reads, kept at count ≥ 2, with the longest anchor seen on each side. It is a yes/no cut per alignment
  record:
  - the test is `anchor_left <= max_left || anchor_right <= max_right` (OR), on exact (ref, start, end), strand-blind;
  - anchors are cumulative across earlier junctions;
  - the count is stored but never loaded;
  - the cut runs per record before the mates are joined, so a junction rejected on one mate but kept on the other
    survives;
  - a fragment is labelled ARTIFACT only when **every** junction is rejected.
- **An ARTIFACT fragment** is held out of calibration entirely. In the EM, gDNA may explain it at a footprint that
  still spans the rejected intron.
- **A certified spliced fragment** has no gDNA candidate. `gdna_can_explain` admits only UNSPLICED, ARTIFACT and
  IMPLICIT. A unique, single-candidate SPLICED_ANNOT fragment skips the EM entirely, on the deterministic path.
- **An unannotated junction.**
  - A fragment whose junctions are all unannotated is SPLICED_UNANNOT. A fragment that also has an annotated junction
    is SPLICED_ANNOT.
  - Its EM candidates are the transcripts it overlaps by exon or span. Each is scored at its own length, which keeps
    the novel gap inside it, with overhang penalties.
  - It gets no gDNA candidate.
  - With no transcript overlap it is "intergenic spliced", which is counted as RNA.
  - **Corrected after the review:** calibration does NOT deposit it as spliced. Its deposit is spliced only when it
    uses an annotated junction (`accumulator.cpp`, `deposit_survivor`: `spliced = !sj_ids.empty()`). A fragment whose
    observed introns are all unannotated goes on the unspliced banks, gDNA-eligible, with its gaps still cut. The two
    stages disagree, as `ISSUES: splicing-artifacts` already lists.
- **Strand.** κ comes from a three-class fit over junctions: genuine at κ, artifact at ½, reversed at 1 − κ. The fit
  labels whole junctions and is used only to clean κ. No within-junction artifact share is estimated.
- **Implicit splices.** These are unspliced pairs whose mate gap could hold an annotated intron. They already carry
  both readings in both stages. Calibration holds the fragment, scores each hypothesis and takes one hard draw
  (owner's ruling). The EM scores each reading in its own likelihood. The two stages decide separately; the EM never
  sees the draw.
- **The evidence the scanner reads** is NM only: never MD, the read sequence, the qualities or the reference. The
  index already reads the FASTA.

## 4. The evidence: the VCaP mix, measured on today's units

**The substrate.** The mix blends two real VCaP libraries in silico. The read name's instrument field gives each
fragment's source:

- the DNA library's instrument prefix is exome **DNA** (150-bp reads), so every junction in it is an artifact;
- the RNA library's instrument prefix is the **transcriptome** (125-bp reads), mostly genuine.

The local index carries **no blacklist**, so today's tree runs as the no-blacklist arm.

**Units:** uniquely mapped, primary records with a motif strand (XS), one per (fragment, junction). That gives 6,517
artifact units and 14,926,264 RNA units. The wrong-strand rate is 0.503 for the artifacts and 0.00026 for the RNA
units. The library's genuine κ is 1.26e-4 (ss 0.99987).

### 4.1 The strand prong, as the owner first framed it

The owner's estimator: each junction's artifact share a_j is fitted from its strand split against κ. A read falls on
the wrong strand at `a/2 + (1 − a)κ`; 10 wrong and 90 right at ss 0.99 give â ≈ 0.18. Each fragment's gDNA reading is
then weighted by its posterior.

Two variants were screened. One learns the prior over a_j from all junctions at once (a grid NPMLE, no constant, with
a reversed branch). The other uses one library-wide share.

| | artifact mass caught | RNA mass lost | DNA units, posterior > 0.5 | RNA units, posterior > 0.5 |
|---|---|---|---|---|
| per-junction share, learned prior | 2,045 of 6,517 (31 %) | 1,961 | 2,031 | 460 |
| library-wide share alone, corrected after the review: the fragment-weighted share, π = 0.000702 | 2,414 (37 %) | 8,066 | 3,278 | 3,842 |
| *as first reported: the NPMLE's junction-weighted mean, 0.00355, wrongly used as a fragment prior* | *3,067 (47 %)* | *30,122* | *3,278* | *3,842* |

The NPMLE stopped at its 500-step cap, so its convergence is unchecked. 32,019 of the 82,001 RNA units at h = 0 come
from multi-junction records, so a local tie is not a whole-fragment tie. Both points come from the review, finding 8.

- **Strand sees only the wrong-strand half.** With the learned prior, the right-strand artifact mass caught is 516 of
  3,239.
- **There is little to pool.** The junctions where DNA is the majority hold 1.5 artifact units each, and 1,741 of
  them hold one. An n-read artifact junction shows no wrong-strand read with probability 2⁻ⁿ.
- **Unstranded libraries get nothing** from this prong.

### 4.2 The reference prong, h

Take a junction read with blocks `[s−a, s) | N | [e, e+b)`. Its unspliced alternative lays the short block contiguous
with the long one. h is the number of reference mismatches between the short block's spliced and contiguous
placements. At h = 0 the read's bases fit both readings, whatever they are, so the aligner's spliced call is no
evidence of splicing.

| short side (local block) | DNA units, h = 0 | RNA units, h = 0 |
|---|---|---|
| annotated, 1–5 bp | 2,005 of 2,154 (93 %) | 79,074 of 770,100 (10 %) |
| annotated, 6–10 | 75 of 210 | 1,309 of 1.21 M |
| annotated, ≥ 11 | 58 of 1,900 (3 %) | 749 of 12.8 M |
| unannotated, any | 222 of 2,253 (10 %) | 869 of 115 k |
| all | 2,360 of 6,517 (36 %) | 82,001 of 14.9 M (0.55 %) |

- On the DNA library at annotated junctions with a 1–10 bp short side, 81 % of the short blocks are on the acceptor
  side, where STAR's annotated-junction bonus breaks the tie toward splicing (2026-09-26 session, 0.7.1 outputs).

### 4.3 The harmful class, and unannotated junctions by their ends

The harmful artifacts sit on junctions where DNA is the majority: a spliced call in a gene the RNA does not express.
There are 2,944 such units.

- **Annotated junctions: 828**, of which 692 have h = 0.
- **Unannotated junctions: 2,116**, 72 % of the class. Most have long anchors and h ≥ 2, and most junctions hold one
  read.

The owner's hypothesis is that **genuine novel junctions mostly reuse known splice sites** (alternative isoforms). To
test it, each unannotated junction's two ends were checked against every annotated intron start and end.

| unannotated units | both ends known | one end known | neither end known |
|---|---|---|---|
| DNA (artifact) | 0 | 65 | **2,188 (97 %)** |
| RNA | 27,141 | 66,821 | 20,557 |
| RNA wrong-strand rate | 0.00004 | 0.002 | **0.142** |
| the class's artifact share from its pooled wrong-strand excess | ≈ 0 | 0.004 | **0.35** |

- 82 % of the RNA library's novel-junction units reuse a known site.
- The RNA library's own "neither known" units run 14 % wrong-strand, so about 28 % of them are artifacts too. The
  DNA labels give only a lower bound on artifacts (0.096 for the class).

### 4.4 What the production blacklist catches

These figures come from 0.7.1's annotated BAM, made on the cluster with the production blacklist (5.8 M junctions).
They cover single-junction records.

| DNA library (all artifacts) | records | rejected |
|---|---|---|
| annotated, h = 0 | 5,160 | 4,709 (91 %) |
| annotated, h ≥ 1 (mostly long anchors) | 1,872 | 1,155 (62 %) |
| unannotated, uniquely mapped | 2,385 | 187 (8 %) |
| unannotated, multimapped | 1,693 | 185 (11 %) |

- **Its cost on the RNA half:** 598,010 records rejected, where a pure-gDNA library of that size predicts about
  20,000. That is 3.4 % of the uniquely mapped single-junction records and 8.2 % of the multi-junction ones.
- **Its unique contribution** is the annotated long-anchor class that h cannot see (1,155 records). Paralogs and
  retrocopies are the hypothesis: `alignable` lists a parent's junction from every alignment of a simulated read,
  secondary ones included.

## 5. The proposal

### 5.1 One rule: gDNA may explain a spliced alignment, at a weight

**The rule.** A gDNA molecule produces the observed spliced alignment with probability `r`. The EM gives gDNA a
candidate on a spliced fragment at `log r` wherever that reading can compete, by the rule that already prunes every
candidate: within `-log(pruning_min_posterior)` of the best RNA reading. The EM then decides fragment by fragment, with what it
already weighs:

- abundance: the locus's calibrated gDNA against the transcripts' RNA;
- strand: ½ under gDNA, κ under RNA. This is the strand prong per fragment, pooled through the abundances, with no
  separate per-junction estimator and the strand evidence never counted twice;
- fragment length.

The deterministic bypass is turned off for a fragment that carries the candidate.

**`r` is a product of mechanisms, each derived, none tuned:**

- **The aligner's tie, at the short side.** Under gDNA the short block's bases came from the contiguous genome. The
  ratio of the bases' likelihoods, gDNA against spliced RNA, is therefore about `(ε/3)^h`:
  - ε is the library's per-base error rate, from NM on uniquely mapped unspliced fragments, not `mismatch_alpha`;
  - at h = 0 the ratio is 1, and the aligner's call carries no information. This is the owner's "edit the alignment",
    made continuous;
  - the read's own mismatches at the short side would refine it, but they need MD, which is not read today (a
    question below).
- **The junction's class.** A library-level artifact share for each class: annotated; unannotated with both, one or
  neither end known (§4.3).
  - On a stranded library it is estimable from the class's pooled wrong-strand excess, as in §4.3.
  - An unstranded library needs another source (a question below).
- **Remote origin and genome structure** stay with the cluster plan: realignment, `alignable`'s per-record origins,
  and the multimapper intergenic hits.

### 5.2 The gDNA reading's geometry

The owner's choice is the largest contiguous aligned block.

- At h = 0 the short block is laid contiguous with the long one: the same molecule, longer by the short side.
- At h ≥ 1, or for a long-anchor novel junction, the other block is clipped.

The footprint sets gDNA's fragment length in the EM and calibration's unspliced placement.

### 5.3 Calibration reads the same weight

A spliced fragment whose gDNA reading competes is offered to calibration's second pass as two placements, the
spliced reading and the gDNA reading, and one is drawn with `r` in the draw. Two changes are needed in the second pass:

- its scorer decides "spliced" from implied introns only, so it needs a per-reading spliced flag;
- its prior must not read the junction's own pass-one deposits, or the spliced reading becomes a decisive zero
  whenever all of a junction's fragments are held.

The volume is bounded by the same pruning rule, applied to `r` alone before any abundance is known, so no new
threshold is introduced. At ε ≈ 10⁻³, `(ε/3)^h` stays above the bound only at h ≤ 1 (an approximation: the other
likelihood terms are taken as comparable). Add the unannotated classes with a non-negligible rate, and on VCaP that is
about 2 % of spliced fragments (§4.2: 82 k units at h = 0 and 179 k at h = 1 of 14.9 M).

### 5.4 The training sets

The strand model (κ) and the RNA length law train on each spliced fragment weighted by its probability of being
genuine.

## 6. The owner's five strategies for unannotated junctions, against the data

| strategy | what §4.3 says |
|---|---|
| 1. discard all unannotated junctions | It removes 2,253 artifact units with 114,519 RNA units, 94 k of which reuse a known site and are mostly genuine. Too costly. |
| 2. deposit them as unspliced, at the longest block | The same cost: every genuine novel isoform becomes a gDNA-eligible unspliced fragment at the wrong length. |
| 3. both readings, spliced and longest block | This is §5's rule for this class. It needs `r`. |
| 4. penalize the RNA reading for unannotated junctions | This is `r` set by class, derived from the class's artifact share rather than chosen. |
| 5. split by whether the ends are known splice sites | The strongest signal in the data: 97 % of DNA artifacts use neither known end, 82 % of RNA units use at least one. |

**Recommended: 3 with 5,** with 4's weight derived per class. On VCaP the classes give: neither known 0.35, one known
0.004, both known ≈ 0.

## 7. Shipping without `alignable`

Neither mechanism in §5 needs `alignable`:

- h needs the FASTA, which the index already reads;
- the end classes need the GTF;
- the class rates and κ come from the library itself.

Dropping it would remove a 26 GB dependency, which makes it a candidate to ship without. The cost is the annotated
long-anchor class the blacklist catches and h cannot (§4.4): 1,155 of the DNA library's 7,032 annotated single-junction
records. Against that, the blacklist rejects about 30× a pure-gDNA library's expectation in genuine reads. Its
mappability is not used at all today.

## 8. Phases

Each phase is A/B'd alone, against the previous tree, on the VCaP mix (both errors by class, with the transcript table
read beside them). The ladder cannot see any of this: its oracle BAM holds no artifacts. So no phase ships before the
aligned ladder, or a second truth-labelled library, has run.

1. **h at index time, and ε from NM.** Per annotated junction, a prefix mismatch count on each side, up to the read
   length.
2. **The EM's gDNA candidate on spliced fragments,** at `log r` with the gDNA reading's length, and the deterministic
   bypass off. Annotated junctions only.
3. **The unannotated classes,** with their rates.
4. **Calibration's second placement,** with the two scorer changes.
5. **The training sets.**
6. **The blacklist's fate**, decided on §4.4's numbers re-measured against phases 2–3. Then the cluster plan for
   remote origin.

## 9. Known limits and open problems

- **One library.** The local truth is one highly stranded library (ss 0.99987) with exome DNA. Unstranded libraries,
  other aligners and other settings are untested.
- **Unstranded libraries.** The class rates (§5.1) are not estimable on them from strand.
- **The two stages decide separately.** Calibration draws and the EM weighs, as they do today for implicit splices.
- **NM only.** `(ε/3)^h` assumes the read matches the spliced placement at the short side; the read's own mismatches
  there need MD.
- **Multi-junction fragments.** `r` is per terminal block. A fragment with two artifact junctions, or an interior
  artifact junction, is not covered by the derivation as written.
- **Thin per-junction strand data.** A per-junction share from strand alone is noise at one or two reads, which is
  why the strand evidence lives in the EM instead.

## 10. Questions for the reviewer

1. Is §5.1's rule sound? A gDNA candidate on spliced fragments at `log r`, decided by the EM, with the strand evidence
   inside the EM rather than in a per-junction estimate. Is anything counted twice?
2. Is `(ε/3)^h` the right form for the aligner's tie? What should ε be, given only NM? Is MD worth reading? How should
   the aligner's own scoring, such as STAR's annotated-junction bonus and genomic-length penalty, enter, if at all,
   for a design meant to be robust to aligner settings?
3. **Unannotated junctions.** Is the end-class split (§4.3, §6) the right partition? Other features to consider
   include intron length, the splice motif, distance to the nearest known site, and whether the junction is supported
   in both mate orientations. How should a class rate be set on an unstranded library?
4. **Geometry (§5.2).** Is the largest-block reading right? How should it treat the clipped bases, and multi-junction
   fragments?
5. **Calibration (§5.3).** Is a hard draw with `r` in the draw consistent with the EM's soft weight, or should the two
   stages share one decision?
6. **The blacklist (§7).** Retire it, keep it as one more factor in `r` (the aligner's simulated rate), or keep the cut
   for the annotated long-anchor class only?
7. **The A/B (§8).** What evidence should gate shipping, given one truth-labelled library and a ladder that cannot see
   artifacts?
8. **An isolated read.** Take one h = 0 spliced read at a junction of a gene with no other RNA evidence: the harmful
   class's typical member. The EM can always build a transcript's abundance from that read alone. Does the
   calibrated gDNA prior decide it, under MAP-EM and under VBEM with their RNA and gDNA pseudocounts? If not, what
   would, without a constant?

## 11. Materials

- **Scripts,** in `$RIGEL_RUNS/prototypes/2026-10-14_splice_design/` on the owner's machine:
  - `extract_spliced.py`: every spliced record of the mix;
  - `strand_prong.py`;
  - `h_prong.py`;
  - `unannot_ends.py`.
- **The 2026-09-26 session's tables** (0.7.1 production blacklist): `$RIGEL_RUNS/prototypes/2026-09-26_splice_artifacts/`.
- **The code,** at `src/rigel/native/`:
  - `resolve_context.h` (classification), `bam_scanner.cpp` (the blacklist cut and calibration's adapter),
    `scoring.cpp` (`gdna_can_explain`);
  - `src/rigel/junction_fit.py`, `src/rigel/second_pass.py`.
- **The earlier plan and the owner's 2026-08 investigation notes,** in `docs/dev/`:
  - `GDNA_SPLICE_ARTIFACTS_PLAN.md`;
  - `GDNA_SPLICE_ARTIFACTS_OVERVIEW.md`, `_RIGEL.md`, `_ALIGNABLE.md`, `_STAR.md`.

## 12. After the external review (2026-10-14)

The review (`docs/dev/SPLICE_ARTIFACTS_REVIEW.md`): **build the competing-placement architecture after revision; do
not build §5 as specified.** Each finding was checked here before acceptance.

| finding | verdict here | checked how |
|---|---|---|
| 1 blocker: `r` multiplies an emission, a likelihood ratio and a posterior share; the EM's abundance then reuses the prevalence | **accepted** | derivation: a class share is P(gDNA \| class), so it already carries prevalence; only the ratio P(c\|D)/P(c\|R) = odds(a_c)/odds(π) may enter |
| 2 blocker: a local edit is not a whole-fragment origin; clipping changes the emission model | **accepted** | the gDNA reading must place every sequenced base, both mates, the whole path |
| 3 blocker: pruning on `r` alone is not the existing rule | **accepted** | arithmetic reproduced: h = 2 gives 1.113e-7, times ½/κ on the wrong strand gives 4.418e-4, above the 1e-4 bound |
| 4: the stage contract | **accepted; §3 corrected** | `accumulator.cpp`: `spliced = !sj_ids.empty()` |
| 5: the calibrated prior does not decide a singleton | **accepted** (answers question 8) | the reviewer's audit reproduces: q = 1, B = 0.05 leaves 0.947 on RNA |
| 6: alignment validity is not molecular origin; a reversed genuine fraction confounds the strand excess | **accepted** | the 0.35 for "neither known" is model-dependent |
| 7: NM does not determine the sequence ratio | **accepted** | the ratio is [ε/(3(1−ε))]^h over bases matching the RNA reference, and needs SEQ |
| 8: the screen's measurement contract | **accepted; §4.1 corrected** | 7,120 wrong of 14,932,781; π = 0.000702; 32,019 of 82,001 RNA h = 0 units are multi-junction. The ranking holds: pooling still wins the two errors |

**The revised direction,** to be derived before anything is built:

1. **Component likelihoods.** Each candidate's weight is `(θ/ℓ)·f(L)·s(o)·e(x)`, with every term counted once. The
   gDNA candidate on a spliced fragment is a physical placement with its own length, strand ½ and a sequence emission
   relative to the RNA reading. Class evidence enters only as a likelihood ratio from a coherent fit, and is neutral
   on an unstranded library.
2. **Placements.**
   - Enumerate the contiguous placements anchored on each side, both mates together, the whole path for a
     multi-junction fragment.
   - The contiguous geometry is also an admissible RNA reading (Axiom 0).
   - A fragment with one artifact junction and one genuine junction keeps its partial RNA path.
   - Remote origin is a separate mechanism, a declared coverage gap until the cluster work.
3. **The sequence emission** from the read's bases at the discriminating positions: new scanner infrastructure (SEQ
   at those positions, reference access). ε is the aligned substitution divergence, Σ(NM − I − D)/Σ(M,=,X), over
   uniquely mapped unspliced fragments.
4. **Pruning** by the existing rule on the full likelihood. Any calibration pre-filter needs a proven bound.
5. **One stage contract.** Calibration and the EM share candidates and emissions, and each keeps its own decision.
   - Fix the unannotated-junction disagreement, the second pass's observed-motif strand bypass, and a held junction's
     support.
   - Training weights are P(RNA under a valid RNA geometry).
6. **Limits stated, not hidden.** The singleton matrix (finding 5) before and after each mechanism. Ties stay ties.
7. **The blacklist** is retired only after a matched two-error A/B. Simulation evidence is kept only as an
   exposure-normalized likelihood.
8. **The gates** are the review's §6, including independent artifact-containing evidence before release.
