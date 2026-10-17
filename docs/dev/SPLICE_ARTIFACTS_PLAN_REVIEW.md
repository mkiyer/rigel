# Splice artifacts: the full plan, reviewed (2026-10-15)

> **Status, later on 2026-10-15:** superseded where it disagrees with the reviewer's revision of
> `SPLICE_ARTIFACTS_IMPLEMENTATION.md` and `INTERGENIC_GDNA_REVIEW.md` — in particular §3.2's N-free-insert rule
> (wrong: two mapped ends can span an unsequenced gap) and §6's hold rule. Threads are parked in `PARKED_WORK.md`.

For the owner, who asked on 2026-10-15 to step back and review the whole plan before anything more is implemented:
"The plan must be completely robust and correct. Then we can implement."

This review consolidates the design now spread over five sandbox files:

- `SPLICE_ARTIFACTS_DESIGN.md`: the first review packet, and in its §12 the verdicts on the first review;
- `SPLICE_ARTIFACTS_REVIEW.md`: the first external review;
- `SPLICE_ARTIFACTS_DERIVATION.md`: the derivation;
- `SPLICE_ARTIFACTS_IMPLEMENTATION.md`: the one-design summary;
- `INTERGENIC_GDNA_DESIGN.md`: the intergenic packet, now with a reviewer.

It reviews that plan critically, foundations first, and ends with the design work still to do, in order. Nothing in it
is built. Sandbox: once approved, it becomes the one design document and the others are retired.

## 1. Why the plan needs this review

- **On 2026-10-14 the design reached a reviewed shape.** Every artifact becomes a competing gDNA placement, both stages
  score it with the full likelihood, the read's own bases settle sequence ties, and both mates and the whole path are
  placed.
- **I then prototyped two pieces before the design was complete:** the local tie as a second hit, and a
  multimapper's intergenic hit. I measured them on a minimap2-aligned panel.
- **The panel exposed defects beneath the design.** The scanner builds placements from a single mate, and the index
  counts unmappable bases as opportunity. Every gDNA term in those measurements read a fragment-length law corrupted
  by both, and the intergenic prototype itself fed one-mate placements to the EM.
- **The design never checked its own foundations:** what a placement is, and what a length counts. This review starts
  there.

## 2. The rulings the design obeys

**Settled before today:**

- Axiom 0: three populations; an artifact is an observation mechanism of gDNA, never a fourth population.
- No magic numbers; a yes/no cut on a continuous quantity is a constant.
- The aligner proposes paths; it does not testify. Its score is never turned into a probability.
- Calibration keeps its hard, seeded draw, and the EM its soft allocation; the two share candidates and emissions.
- Class evidence enters only as a likelihood ratio, and is neutral where it cannot be identified (unstranded).
- Pruning acts on the full likelihood only; no partial-product screen without a proven bound.
- The release waits for validation on independent artifact-containing evidence.
- The genome is stored in the index; `mappy` runs at index time; unannotated copies are an accepted gap; source C is
  derived (all 2026-10-14).
- Calibration holds the fragments whose own bases do not favour the spliced reading (2026-10-14).

**New (owner, 2026-10-15):**

- **A placement holds both mates.** A paired-end fragment's alignment hypothesis places both mates. Mates placed
  apart are a chimera, a real but much less likely explanation when a valid paired placement exists. There is no
  single-mate fragment.
- **Abundance is per mappable base.** Unmappable sequence never enters a length that turns counts into a density.
- **The length defects are release-blocking.**
- **Design before implementation.**
- **STAR runs on the cluster,** once implementation and testing are done.

## 3. The foundations, and what is broken in them

### 3.1 A placement holds both mates

**Today.** `bam_scanner.cpp`, `group_records_by_hit` and `pair_multimapper_reads`, unchanged since the first commits
(2026-03-07), make a one-mate hit three ways:

1. A primary whose mate is unmapped is the fragment's only hit.
2. A secondary flagged mate-unmapped is a "singleton hit".
3. Without HI tags, each secondary record is a location for its own mate. Locations pair by shared transcripts, then
   by the closest distance on one chromosome, then by every remaining R1 × R2 combination. Whatever is left becomes a
   one-mate "SINGLETON" hit.

The resolver's chimera rule, under the existing maximum fragment length, then rejects pairs it cannot read as one
molecule, but nothing ever rejects a one-mate hit. `collect_placement` deposits it, the EM scores it, and
`num_hits = max(NH, hits)` counts it, chimeric pairs included.

**Scope** (`prototypes/2026-10-15_mate_pairing/count_one_mate.py`):

| input | fragments with an orphan secondary | multimappers only through one-mate hits | primary with mate unmapped |
|---|---|---|---|
| minimap2 panel, `g50 ss.99 OFF` | 771 k (7.7 %; 3.9 M orphan records) | 501 k (5.0 %) | 105 k (1.05 %) |
| minimap2 panel, `g00` | 510 k (5.1 %) | 346 k of 1.05 M "multimappers" | 25 k |
| VCaP, STAR (HI tags) | — | 4 | 318 one-mate fragments in 18.7 M |

STAR pairs its alignments and by default requires both mates, which is why no real library and no oracle panel
showed this.

**Consequences.**

- A one-mate hit's footprint is one read, so its length is the read's span (100 bp here). It enters calibration's
  length pools, the EM's length term, and the multimapper count.
- It is one of the two causes of the gDNA law's collapse on aligned input (§3.3).

**The rule, the ruling made operational:**

1. A placement is a pair of mate alignments that the resolver accepts as one molecule (not chimeric, under the
   existing maximum fragment length).
2. Without HI tags, every R1 location is tried with every R2 location, and the pairs that resolve as one molecule are
   the placements. This replaces the three pairing heuristics and the singletons. A mate record that no placement uses
   contributes nothing.
3. With HI tags, an HI group holding one mate only is no placement.
4. NH is the number of the fragment's distinct placements.
5. A fragment with no placement is an observation for neither stage. It is "mate unplaced" when a mate has no
   alignment, and chimeric when both mates are placed but never as one molecule. Both are counted in the summary.

**Open:**

- An aligner that writes NH above the placements it reports (STAR with `--outSAMmultNmax`) leaves placements
  unseen. Is the reported set the hypothesis space, with the deficit only recorded?
- A mate the aligner splits (a supplementary record) is merged into that mate today and judged by the chimera rule.
  Keep?

**Verification.** Fixtures, each verified failing first: a primary with its mate unmapped; an orphan secondary beside
a valid pair; two R2 secondaries against one R1; an HI group with one mate. Then each fixed gate is broken to watch
it fire. Then the minimap2 panel (the length laws, both errors) and VCaP, where it should be nearly a no-op.

### 3.2 Abundance per mappable base

**Today.** The index's region lengths include N bases.

| index | intergenic length | of which N | genic length that is N |
|---|---|---|---|
| the panel genome (chr21, chr22, ERCC) | 51.0 Mb | 18.2 Mb (36 %) | 0.1 % |
| GRCh38, VCaP's index | 1,318 Mb | 150 Mb (11.4 %) | 0.05 % |

**Consequences on real and aligned data:**

- a region that contains N reads its density low by its N share, and most of the N sits in a few large gap regions
  that read as empty;
- the length module's intergenic purity `a_0 = ρ·E_0/n_0` is overstated: `E_0` counts the gaps, while `n_0` has no
  reads there;
- an empty N gap reads as a depleted object to the landscape and to the one-sided density.

On the oracle ladder the simulator places gDNA in the gaps, so none of this showed.

**The rule.** An opportunity counts a fragment start only if the whole fragment lies in placeable sequence. Taking N
runs as holes, a region's opportunity at length L is the sum, over its N-free segments, of max(0, segment − L + 1).
That is today's geometry applied per segment, with no new constant; the genome stored in the index gives the N runs.
Sequence where reads can only multimap still counts ("unique or multimapping regions").

**Consequence for the oracle ladder.** Its N-gap reads would land where the opportunity is zero. They must leave the
oracle BAMs (they are identifiable: all-N sequence), or the ladder must be re-simulated. That is a decision (§10).

**Verification.** A brute-force enumeration of the opportunity with N runs, to the deposit rule's standard; then the
calibration instruments and the A/B.

### 3.3 The gDNA length law's conditioning

**Today** (`calibration/fl.py`, `_deconvolved_gdna_counts`). The law is a two-pool contrast:

```
g = [ (1 − a_1)·f_0 − (1 − a_0)·f_1 ] / (a_0 − a_1)
```

It assumes the two contained pools' contaminants share one length distribution. Any difference between the pools that
is not the gDNA/RNA mixture is amplified by `1/(a_0 − a_1)`. The resolution weight `sep²/(sep² + SE²)` guards only
against count noise in the purities, not against a shape the pools do not share. On the aligned panel the read-length
spike differed fourfold between the pools (8.9 % against 2.0 %), and the law put 76 % of its mass on it. VCaP's STAR
BAM is unaffected (its gDNA law has mean 276 bp).

**What to do.** First re-measure after the paired-placement rule and mappable opportunity, since both feed the spike and `a_0`. If the law recovers, the length law's conditioning stays a
robustness item, pinned by a stress fixture: a spike present in one pool must not dominate the law. If it does not
recover, the estimator needs a derivation that tests its own premise. the length law's conditioning is not a splice item, but every gDNA
candidate reads this law.

## 4. The placement model, consolidated

**The observation** is a paired-end fragment. **A hypothesis** is a component (a transcript, or gDNA) with a
placement: a two-mate path accepted as one molecule (§3.1).

**The weight** of component k at placement P:

```
w(k, P) = (θ_k / ℓ_k) · c(P) · f_k(L_k(P)) · s_k(o) · e(x | P)
```

| term | meaning | status |
|---|---|---|
| `θ_k / ℓ_k` | abundance per mappable, capture-weighted position | settled; ℓ changes with mappable opportunity |
| `c(P)` | capture efficiency of the objects P sits on; cancels within a placement, not across placements | derived (intergenic packet §4); not built |
| `f_k(L)` | the length law at the candidate's own molecule length | settled; the gDNA law is broken on aligned input (the length law's conditioning) |
| `s_k(o)` | κ or 1 − κ under RNA, ½ under gDNA | settled |
| `e(x \| P)` | the read's bases along P: across placements only the discriminating positions differ, `[ε/(3(1−ε))]^(m_s − m_c)`, with ε the aligned substitution divergence | derived; needs SEQ and the stored genome |

**The fragment's likelihood sums over placements for each component:** `P(u, k) = Σ_P w(k, P)`.

- Today each transcript keeps its best placement, and gDNA averages its placements (the mean).
- The sum is the derived form. On inputs carrying the three foundation defects it removed the mean's dilution under capture but made the
  capture-OFF tables slightly worse. It gets re-measured on the repaired foundations before any decision.

**An out-of-locus placement** (an intergenic copy) has no `θ_g/ℓ_g`. Its rate is calibration's gDNA density per
mappable base in its own region, held fixed inside the locus's EM, and its mass is booked as intergenic gDNA (the
intergenic packet's arm 2, now under review).

**Where placements come from:**

| source | the placement | found by | when |
|---|---|---|---|
| **D**, the aligner | each two-mate alignment it reports | the aligner | scan |
| **A**, the local contiguous origin | each spliced mate laid contiguous, anchored on either side, joined with its mate into one molecule; partial paths for RNA | the tool, from the stored genome | scan (index-time tie lengths speed h = 0) |
| **B**, a copy elsewhere | the fragment's whole transcript-coordinate span on a copy of the transcript | an index-time copy table (`mappy`) | index, then scan |
| **C**, the sample's allele | gDNA at the reported spliced path | the tool | scan; the allele fraction fitted in the EM |

- **Duplicates.** One physical placement counts once, for example an aligner secondary that is also a listed copy.
- **Pruning** is the existing full-likelihood rule.

## 5. The sources, reviewed

### A — the local contiguous origin

- **Settled.** The geometry (derivation §2) and the emission (§3). At h = 0 only the reference is needed, at h ≥ 1
  the read's bases. On VCaP's DNA library, 2,124 of the 2,626 annotated artifact units are ties, and none of the rest
  vote contiguous: they carry the spliced bases, which is source B's signature.
- **The prototype (step 1)** handled ties only, as a second EM hit, with calibration unchanged. It copies both mates.
  - On VCaP it moved 2,044 of 2,363 DNA tie artifacts to gDNA, at a cost of at most 441 RNA-library fragments (about
    180 genuine, by the strand split).
  - Side effect: isoform redistribution (ADAR 3,545 → 134). A tie's contiguous reading fits inside another isoform's
    longer exon. The physics is plausible: STAR's annotated-junction bonus splices such reads. But it is not verified
    without isoform truth.
- **Open:**
  - both mates spliced, at different junctions;
  - a multi-junction mate: one gDNA path per anchor, and the partial RNA paths;
  - the alternative must itself pass §3.1's rule;
  - the isoform side effect, judged on the minimap2 panel's truth once the three foundation defects are repaired, with a fixture: a read
    whose contiguous reading fits isoform B's longer exon.

### B — a copy elsewhere

- **Settled.**
  - Of the annotated h ≥ 1 DNA artifacts that fit the reference contiguously at all (1,056), 996 fit another locus.
  - The index-time table lists 2,658 of 404,168 junctions; it recalls 533 of 558 (95.5 %) of VCaP's retrocopy
    artifacts and exposes 3.2 % of RNA units.
  - The facts are per position, so read length is not a parameter.
  - Scoring is an out-of-locus placement (§4).
- **Open:**
  - the table's exact content: each transcript's copies as covered transcript intervals with their differing
    positions, so any fragment's span, mate included, can be tested;
  - a copy inside another gene, which brings its own RNA readings;
  - admitting a near-exact copy: the emission prices the differences, so the admission rule must be derived, not
    tuned to `mappy`'s score;
  - duplicates with aligner secondaries.

### C — the sample's own allele

- **Settled (draft, derivation §9).**
  - Each junction has an allele fraction `α_j`.
  - It is identifiable by gDNA's capacity, and by the reference allele's reads.
  - A class likelihood ratio comes from the annotation.
- **A new gap: non-reference retrocopies.** A processed-pseudogene insertion present in the sample but not in the
  reference gives gDNA reads spliced exactly at annotated junctions, at every junction of the gene, with no contiguous
  origin in the reference.
  - Only source C can explain them, and its class prior assumes a random breakpoint, which a retrocopy is not.
  - The derivation must cover a whole-transcript allele: one α shared by every junction of the transcript (one
    insertion). It is identifiable by capacity at every junction, and by strand where the library is stranded.
  - Candidates on VCaP: about 29 % of the annotated h ≥ 1 DNA artifacts (roughly 430) fit nowhere in the reference
    (derivation §10).
- **Open (from the draft):**
  - the opportunity `W_j` when the ends lie in different regions;
  - how a subclonal α and a germline α share one prior;
  - the annotation-withholding controls.

### D — the aligner's own secondaries

- **After the paired-placement rule:** only two-mate placements count.
- **Open:**
  - the sum against the mean;
  - each placement's own capture efficiency;
  - the out-of-locus density (the intergenic review);
  - alignment quality beyond NM. Soft clips are unpriced and `mismatch_alpha` is underived. One derived emission over
    every sequenced base would serve A, B and D alike.

## 6. The stages

**The EM.**

- Every placement is scored with every component, under full-likelihood pruning.
- The deterministic bypass is off for a fragment with alternatives.
- Out-of-locus placements enter at a fixed rate, and their mass is booked as intergenic gDNA.

**Calibration.**

- **The adopted hold rule is a screen on the emission alone,** and it has an expected bound, not a per-fragment one.
  - A fragment it does not hold is disfavoured by at least one discriminating base: a factor of at most
    `ε/(3(1−ε))`, about 3e-4.
  - But on the wrong strand of a stranded library `½/κ` is about 5,000, and in a gene with no RNA the abundance ratio
    has no bound.
  - Such fragments are gDNA molecules whose read carries an error matching the spliced reference: about ε/3 of the
    local artifacts, a handful on VCaP.
  - The EM still scores them with the full likelihood; only calibration's deposit carries the approximation.
  - **Recommended:** keep the rule; pin it with a fixture (a wrong-strand, h = 1, base-favoured read in a gene with
    no RNA); and measure on VCaP how many unheld spliced fragments the EM's full likelihood calls gDNA. If that number
    is not negligible, adopt the exact deposit-then-move (derivation §7 b).
- **The second pass** draws among a held fragment's placements with the same factors:
  - a held junction's support is absent, not zero;
  - the motif bypass is fixed (41899ae4);
  - source C's placement removes the disagreement over unannotated junctions, since both stages then see both
    readings;
  - only two-mate placements are offered (the paired-placement rule).

**The training sets.** κ and the RNA length law weight each spliced fragment by its probability of being RNA under a
valid RNA geometry. Weighted fit statistics stay apart from integer deposits.

**The blacklist** is retired after the matched two-error A/B, and `alignable` with it.

## 7. Robustness

**Across aligners:**

| behaviour | STAR (defaults; VCaP) | minimap2 `-x splice:sr` (the panel) | what the design needs |
|---|---|---|---|
| NH / HI | written | not written | the paired-placement rule must not depend on either |
| secondaries | paired | per mate, at least 80 % of the primary's score | only two-mate placements |
| one mate placed | not reported by default | reported: 1 % of fragments at `g50` | not an observation |
| annotated-junction preference | the `sjdbScore` bonus | the `-j` junction bonus | the tool's own alternatives (A–C) |
| soft clips | rare | frequent on secondaries | an emission over every sequenced base |
| junction strand tag | XS | `ts`, flipped on reverse-strand records | as today |

HISAT2 is not measured; it reports one-mate alignments unless `--no-mixed` is set, so the paired-placement rule matters there too.

**Across strata.**

- Stranded libraries carry strong strand evidence; on unstranded ones strand is neutral, and emission, length,
  density and the class likelihood ratio decide (neutral when unidentified).
- Capture: each placement's own efficiency.
- `g00` must stay free of false gDNA, by fixtures and the `g00` rows.

**Library and annotation.**

- Overlapping mates and adapter read-through: the emission must count each sequenced base once (review: shared
  variation across overlapping mates needs consistent treatment).
- Annotation completeness: the annotation-withholding controls (review, finding 6).

## 8. What the prototypes taught, and what must be redone

- **Source A step 1** is physically sound, places both mates, and was measured on VCaP and on the panel. Its isoform
  side effect is real and unjudged; it is re-measured once the three foundation defects are repaired.
- **The intergenic repair (arm 1) and the sum arm** were built on one-mate placements. At capture-ON, 118 k of the
  275 k mixed gDNA fragments had only a one-mate intergenic placement. Their numbers are not evidence for the design.
  What survives holds on the paired subset: the mean dilutes gDNA's term (33.4 % of the genic-origin fragments have a
  worse paired intergenic placement).
- **The minimap2 panel** is a valid local instrument once the three foundation defects are fixed, and its baseline is re-recorded then.

## 9. Validation design

- **The specification is the fixtures:** derivation §11, plus the paired-placement rule's and mappable opportunity's.
- **Harm:** the oracle ladder and both fragment-length gap panels, with the N-gap reads handled per §10.
- **Effect, locally:** the minimap2 panel (exact labels) and VCaP (source labels; the RNA side read as an upper bound
  with the strand-implied split).
- **The gate:** the STAR-aligned ladder on the cluster after implementation and testing (owner), with the review's
  aligner-setting controls.

## 10. Decisions and open questions

1. **the paired-placement rule.** The rule in §3.1 follows the ruling. Is there anything to add, such as reporting "mate unplaced" fragments
   in `summary.json`?
2. **mappable opportunity.** N only, for now. The oracle ladder's N-gap reads: filter them out of the BAMs, or re-simulate the ladder?
3. **The hold rule.** Keep it with its stated approximation and the VCaP measurement, or move to the exact
   deposit-then-move?
4. **Source C.** Should it cover non-reference retrocopies at annotated junctions?
5. **The multimapper sum** is decided on the re-measurement after the three foundation defects.
6. **External review** of the complete design before any implementation?

## 11. The design work, in order

Nothing is implemented until the design is complete.

1. **the paired-placement rule and mappable opportunity:** each written up as a one-page rule plus the fixtures that specify it. They are bugs under the
   rulings, not design choices, so the owner may let them proceed to prototype and A/B ahead of the splice design.
2. **Re-measure** the minimap2 panel on the repaired foundations; decide the length law's conditioning from it.
3. **The placement model (§4)** as one derivation: the sum, each placement's efficiency, the out-of-locus density,
   duplicates, and an emission over every sequenced base (soft clips, overlapping mates).
4. **Each source's derivation completed:** A's multi-junction and both-mates cases, B's table, C with retrocopy alleles.
5. **The stage contract (§6)** written as fixtures.
6. **The validation design** with its stop criteria.
7. **An external review** of the complete design.
8. **Then implementation,** phase by phase, each A/B'd.
