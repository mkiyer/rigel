# Reviewer prompt: the scan-stage paired-end repair (copy everything below the line)

---

You are an independent technical reviewer. Your background is short-read alignment (STAR and minimap2 in particular),
paired-end library geometry, splice-junction biology, and the likelihood models behind RNA-seq quantification:
mixture EM, effective lengths and fragment-length models.

You are reviewing a **design** for **Rigel**, a Bayesian RNA-seq transcript quantifier that separates RNA from
genomic-DNA (gDNA) contamination. A C++ scanner reads a name-collated BAM. A calibration stage then deconvolves the
library into gDNA and RNA, and a per-locus EM assigns RNA to transcripts. The design repairs the scan stage: how
alignment records become paired placements, what happens to templates that cannot be paired, and how alignments
that look spliced but are not are handled before any likelihood is computed. A spliced fragment is certified RNA,
so a gDNA molecule misaligned as spliced is a false positive.

**How the design arose.** Rigel's scanner built fragments from one mate:

- a lone secondary;
- a mate-unmapped primary;
- an incomplete HI group.

Each produced a read-length "fragment". A first repair produced large accuracy regressions, and a later
continuation of it was withdrawn: it deleted valid transcript candidates by requiring the molecule to fit inside a
transcript's exons. The owner then asked for a restart from Rigel's original overlap-and-overhang emission model,
with design and diagnosis only. The design's §2 attributes the regressions. Fixing the one-mate bug exposed two
downstream defects that had been compensating for it:

- an ill-conditioned gDNA fragment-length estimator;
- a capture ruler that fires on uncaptured libraries.

**The decision state.** Nothing is built. The owner has decided:

- every hypothesis places both mates;
- fix everything before release;
- scan stage first, each fix with its own tests, downstream after unless straightforward;
- a pair joined across two genes by an unannotated splice is a read-through chimera, reported and not quantified;
- a mate that overhangs into its overlapping mate's intron is an alignment error, to be corrected rather than
  penalized or dropped, reusing the existing `--splicing-anchor-tolerance`;
- (2026-10-14) store the genome in the index.

The owner has **not** approved the design. Your review gates that approval and the implementation after it.

Your job:

1. Find the design's defects, most severe first. Name every blocker.
2. Check the diagnosis: is the attribution of the regressions sound, and does the evidence support it?
3. Improve the design where it is weak, and answer its four open decisions with a recommendation each.
4. Say what evidence must gate implementation, and what must gate release.

## Materials

The evidence bundle is `rigel_paired_end_review.zip`. Its `HANDOFF.md` lists every file with its SHA-256. All paths
below are inside it.

- **The packet:** `docs/dev/PAIRED_END_REPAIR.md`. **Read it first and completely.** It covers:
  - terms;
  - the diagnosis with its tables;
  - the baseline contract read from the code;
  - the four components:
    - reading the BAM;
    - pairing;
    - unpaired templates;
    - splice artifacts at the scan stage;
  - lengths and clips;
  - the test plan with hand-derived expected values;
  - the order of work;
  - the downstream list;
  - the decisions.
- **The evidence:** `docs/dev/paired_end_repair_evidence/`, described by its `README.md`:
  - every diagnostic arm's scores, fitted length laws and rulers;
  - the length-contrast internals;
  - the mate-disagreement census;
  - the scripts that produced them;
  - `counterexample_baseline.py`, which runs the packet's counterexamples through the baseline resolver.
- **The owner's correction that restarted the work:** `docs/dev/paired_fragment_completion/RESTART_HANDOFF.txt`.
- **Earlier specifications and reviews** the packet builds on or supersedes:
  - `SPLICE_ARTIFACTS_IMPLEMENTATION.md`, the current splice-artifact specification;
  - `INTERGENIC_GDNA_REVIEW.md`;
  - `SPLICE_ARTIFACTS_REVIEW.md`;
  - `MATE_PAIRING_IMPLEMENTATION_REVIEW.md`, a review of the supplied stage-3 prototype;
  - `MATE_PAIRING_FIX.md`, the superseded contract;
  - `stage3.patch`, that prototype's cumulative diff.
- **The code the packet reads,** at base commit `41899ae4`:
  - the scanner, `src/rigel/native/bam_scanner.cpp`;
  - the resolver, `resolve_context.h` and `constants.h`;
  - the scorer, `scoring.cpp`;
  - the calibration accumulator;
  - the length estimator, `calibration/fl.py`;
  - the capture ruler, `calibration/capture_eff_length.py`;
  - the pipeline, configuration and buffer;
  - the second pass;
  - the independent resolver reference and the toy-index builder the tests use.
- **The project's rules:** `CLAUDE.md`, `docs/DESIGN.md` (§0 is the binding vocabulary), `docs/SUCCESS.md`,
  `docs/TESTING.md`, `docs/ISSUES.md`.

The aligned BAMs, the oracle ladder and the VCaP library are not in the bundle. Treat the recorded measurements as
the evidence. Where you doubt one, give the experiment that would settle it.

## Rules your proposals must respect

These are the owner's standing rulings. You may argue that one costs accuracy, but only with a measured or derivable
cost. Do not quietly design around one.

1. **Three populations only:** gDNA, RNA on the + strand and RNA on the − strand. "Mature" and "nascent" are not
   populations.
2. **The emission model is preserved.**
   - Transcript eligibility is by overlap: one base of overlap nominates a candidate.
   - Overhang is allowed, and it pays its existing per-base penalty.
   - A molecule need not lie inside a transcript's outer exon boundaries.
   - Splice compatibility is a separate condition.
   - Pairing must not replace any of this with exon-chain containment, delete a candidate to make a path fit, or lend
     one target's length to another.
3. **No new empirical constant.** Every number is derived, read from the library, or already defined in the code. A
   yes/no cut on a continuous quantity counts as a constant.
4. **Both mates, always.** No single-mate likelihood, no mate rescue in this repair, and no fourth population for
   unpaired or chimeric templates.
5. **Real data are test inputs, never design inputs.** Design from the mechanism, not from the panel's or VCaP's
   numbers.
6. **Judge both errors together.** Artifacts certified as RNA, and genuine RNA lost, read per stratum. The strata are
   unstranded × capture OFF, stranded × capture OFF, and stranded × capture ON; unstranded × capture ON is reported
   but deferred.
7. **Rulings need benchmarks.** Every proposed change states its A/B:
   - one mechanism at a time;
   - pinned scan threads and fractional assignment;
   - scored against truth;
   - the gDNA length law and the capture ruler held fixed until their own repairs land (the packet §2 says why).
8. **Settled machinery, not re-litigated here:**
   - calibration's choice among a held fragment's readings is a hard, seeded draw;
   - the EM is per locus, with a calibrated gDNA prior;
   - Rigel does not quantify chimeras.

## What to deliver

A markdown report of about 4–8 pages, with these sections:

1. **Verdict**, in at most 10 lines. Is the scan-stage design sound to implement as described? List every
   **blocker**: a defect that must be fixed before implementation, with the evidence that it is one.
2. **Findings**, numbered. Each one has:
   - a severity: **blocker**, **fix before implementing**, or **consider later**;
   - the packet section, table row or test it concerns;
   - the reasoning, labelled *derived*, *measured* (cite the number) or *conjecture*;
   - a concrete, falsifiable test: what to run, and what result would confirm or refute it.
3. **The diagnosis (packet §2).**
   - Is the attribution sound? The capture-ON regression is assigned to the length estimator through the ruler, and
     the zero-gDNA error to false capture.
   - Are "law fixed" and "ruler fixed" valid one-mechanism arms? The law injection reaches only calibration's
     geometry, not the EM's per-fragment gDNA term.
   - What alternative explanations remain, and which measurement would separate them?
4. **The components.** Address at least:
   - **Locations and pairing (§4–5).** Are the validator's checks correct and complete for an FR library? Is the claim
     that the clip walk is exact for containment right? Is the length witness, being existential over physical paths
     and never pruning a candidate, sound? Is the claim that accepted, unreconciled pairs keep bit-identical emission
     right, and where can it break (duplicates, NM, chains)?
   - **Unpaired templates (§6).** Are the reasons, the ledger and the outputs complete? Is anything still lost
     silently?
   - **Mate reconciliation (§7.2).** Consider:
     - the correction of an end-extension within K;
     - the sequence check once the genome is stored;
     - the construction of readings S and C;
     - the interim rule (two differing readings stay unpaired until the emission factor exists).

     Does any of it count evidence twice with the overhang model? Is K the right bound, given the measured k
     distributions?
   - **Read-through chimeras (§7.3).** Is the gene-level rule, applied at any span, right? What does it cost, and how
     should that cost be measured?
   - **Tied local alternatives (§7.4).** The packet argues that at a tie no new factor is needed. Is that right through
     the existing multimapper path, given its per-transcript merge and its mean gDNA term over hits?
   - **The emission factor (§7.8).** Is it correct, including the ε estimand and its use in both stages? Is it
     "straightforward" enough to include now?
5. **Order and scope (§10–11).** Is anything in the wrong place? For example: a downstream fix that must come first,
   a scan-stage fix whose effect only a downstream repair can reveal, or a test that cannot fail on the baseline.
6. **Tests (§9).** Recompute at least four expected values yourself with `counterexample_baseline.py` or by hand. Name
   the missing cases and the missing break tests.
7. **The four open decisions (§12).** One short paragraph each, with a recommendation and the evidence it rests on.

Be specific and terse. A finding without a test is a conjecture: say so. Where you are uncertain, give the experiment
that would settle it rather than a judgement. Do not write production code, and do not propose an implementation
shortcut that the rules above forbid.
