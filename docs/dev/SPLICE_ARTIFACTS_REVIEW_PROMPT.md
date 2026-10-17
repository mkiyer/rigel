# Reviewer prompt: spliced-fragment artifacts (copy everything below the line)

---

You are an independent technical reviewer. Your background is short-read alignment (STAR in particular),
splice-junction biology, and mixture models and EM for RNA-seq quantification.

You are reviewing a **design** for **Rigel**, a Bayesian RNA-seq transcript quantifier that separates RNA from
genomic-DNA (gDNA) contamination. The design is for its handling of **spliced-fragment artifacts**: alignments that
cross a splice junction on a molecule that was never spliced there. Rigel treats a spliced fragment as certified RNA,
so each such artifact is a gDNA fragment counted as RNA. Today's guard, a splice-junction blacklist, also rejects
large numbers of genuine junction reads.

**The decision state.** The project owner agrees with the proposal in the packet's §5, subject to your review. Nothing
has been built. The problem blocks the release (0.8.0).

Your job:

1. Find the proposal's weaknesses.
2. Improve it.
3. Rank the owner's five strategies for unannotated junctions.
4. Say what evidence should gate shipping.

## Materials

- **The review packet:** `docs/dev/SPLICE_ARTIFACTS_DESIGN.md`. **Read it first and completely.** It covers:
  - Rigel in brief;
  - today's behaviour, read from the code;
  - the evidence measured on a truth-labelled library;
  - the proposal;
  - the owner's five strategies for unannotated junctions;
  - shipping without the companion tool `alignable`;
  - the phases, the known limits and eight questions.
- **The evidence bundle:** `rigel_splice_review.zip`, described by its `HANDOFF.md`. It holds:
  - the packet and this prompt;
  - the project's rules and rulings (`CLAUDE.md`, `docs/DESIGN.md`, `docs/ISSUES.md`);
  - the code the packet reads;
  - the measurement scripts with their captured outputs;
  - the owner's earlier investigation notes.
- **The truth-labelled library is the VCaP mix:** two real libraries blended in silico, exome DNA and transcriptome. The
  read name gives each fragment's source, so every spliced fragment from the DNA half is an artifact. The BAM itself
  is private and not in the bundle.

## Rules your proposals must respect

These are the owner's standing rulings. You may argue that a ruling costs accuracy, but only with a measured or
derivable cost. Do not quietly design around one.

1. **Axiom 0.** Every genomic object holds at most three populations: gDNA, RNA on the + strand and RNA on the −
   strand. RNA inside an intron is RNA that has not spliced there; "mature" and "nascent" are not populations.
2. **No new empirical constant.** Every number must be derived, read from the library, or already defined in the code
   (for example, the EM's pruning bound). If a proposal needs a constant, name it as a discussion item; do not hide it
   in a functional form. A yes/no cut on a continuous quantity counts as a constant.
3. **Real data are test inputs, never design inputs.** Design from the mechanism, not from VCaP's numbers. VCaP is one
   highly stranded library; do not tune to it.
4. **Judge both errors together:** artifacts certified as RNA, and genuine reads rejected. Never read one alone.
5. **Robust to aligner settings.** A mechanism that depends on one aligner's scoring constants must say so, and say
   what happens under other settings.
6. **Rulings need benchmarks.** Every proposed change states how it would be A/B'd: one mechanism at a time, against
   the current tree, on the same inputs, with both errors and the transcript table read per library.
7. **Settled machinery, not re-litigated here:**
   - calibration's choice among a held fragment's readings is a hard, seeded draw (whole counts);
   - the EM is a per-locus MAP-EM/VBEM with SQUAREM and a calibrated gDNA prior;
   - implicit splices (an unspliced pair whose mate gap could hold an annotated intron) carry both readings in both
     stages today.

## What to deliver

A markdown report of about 3–6 pages, with these sections.

1. **Verdict** on the proposal (packet §5), in at most 10 lines. Is it sound to build as described? List every
   **blocker**: a defect that must be fixed before building, with the evidence that it is one.
2. **Findings**, numbered. Each one has:
   - a severity: **blocker**, **fix before building**, or **consider later**;
   - the packet section or table row it concerns;
   - the reasoning, labelled as *derived*, *measured* (cite the number) or *conjecture*;
   - a concrete, falsifiable test: what to run, and what result would confirm or refute it.
3. **Improvements**, ranked by expected benefit against risk and effort. Address at least:
   - **The rule's statistics.** Does adding a gDNA candidate at `log r` count any evidence twice? Is the strand
     evidence better placed inside the EM than in a per-junction estimate fitted beforehand (§4.1)? Does the
     identifiability problem in question 8 defeat the rule for the harmful class?
   - **The aligner's tie.** Is `(ε/3)^h` right? What should ε be, given only NM? Is reading MD worth it, given that
     today the scanner reads neither MD, the sequence nor the reference? How should the aligner's own scoring enter,
     if at all?
   - **The geometry** of the gDNA reading (§5.2), including the clipped bases and multi-junction fragments.
   - **Consistency between the stages** (§5.3): calibration's hard draw against the EM's soft weight.
4. **Unannotated junctions** (packet §4.3 and §6).
   - Rank the owner's five strategies.
   - Is the end-class partition (both, one or neither end a known splice site) the right one, or should other features
     join it? Candidates: intron length, splice motif, distance to the nearest known site, support in both mate
     orientations, read depth.
   - How should a class's artifact share be set on an **unstranded** library, where the strand excess cannot be read?
5. **The blacklist and `alignable`** (packet §4.4 and §7). Retire it, keep it as one more factor in `r` (the
   aligner's simulated rate), or keep the cut only for the annotated long-anchor class h cannot see? Weigh the 1,155
   catches it alone makes against the ~30× excess of genuine reads it rejects.
6. **Validation.** What evidence should gate shipping, given:
   - one truth-labelled library;
   - a benchmark panel whose simulated alignments hold no artifacts at all;
   - an aligned panel and further truth mixes only available on a cluster, later?

   Name the minimum set, and what result on each would stop the release.
7. **Answers to the packet's §10 questions,** one short paragraph each.

Be specific and terse. A finding without a test is a conjecture: say so. Where you are uncertain, give the experiment
that would settle it rather than a judgement.
