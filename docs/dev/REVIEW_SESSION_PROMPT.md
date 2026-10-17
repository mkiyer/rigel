# The prompt for the review-investigation session

*Paste the block below as the first message of a new session, with the reviewer's critique appended where
marked. Sandbox file.*

---

Read, in this order and before anything else: `CLAUDE.md`, `docs/SUCCESS.md`, the memory index and the two
memory notes on the ladder re-score of 2026-10-11 and the capture-weight expectation campaign, then
`docs/dev/CAPTURE_WEIGHTS_AND_DEFERRED_STRATUM_REVIEW_BRIEF.md` (the brief an external reviewer was given) and
`docs/dev/RNA_CAPTURE_WEIGHT_EXPECTATION_PLAN.md` (the two refuted readouts, with their A/B tables).

The state: `main` is converged and pushed (the reader's speed campaign, the per-node share lattice and the
multimapper second pass are landed; the suite's standing count is in `CLAUDE.md`). The ladder, both gap
panels and the test chromosome were re-scored on 2026-10-11 and their arms are current; earlier arms are
preserved beside them. Two problems are open and attributed: (1) the transcript-level regression on the
stranded strata, owned by the capture weights' noise on a capture-OFF field (g00 stranded OFF 13.61 → 1.93 %
with every weight forced to one); (2) the deferred unstranded × capture-ON stratum, where captured gDNA at
exons is booked as nascent RNA since the intron constraint's removal, and which the gap panels show was
workable with a length gap before that change.

Below is the external reviewer's critique of the brief. Treat it as a hypothesis set, not a ruling: for each
point, say what it claims, what measurement would confirm or refute it, and run that measurement on the
cheapest adequate instrument (the test chromosome first, the ladder and the gap panels to decide, VCaP for
real data), pinned and fractional, read per stratum, one mechanism at a time, scored against truth. Where a
point proposes a mechanism, derive it on paper in a sandbox note before any prototype, prototype outside the
main tree, and A/B it against what ships on every in-scope stratum and both gap panels before judging it.
The constraints in the brief's §6 are fixed: no detector, no library-level reference, no new constant, no
tuning to the panels. Do not commit; the owner drives commits. Budget every run and say what each will cost
before launching it; keep the machine idle for timings.

Deliver, at the end: a verdict per critique point with its numbers; the single mechanism you would take next
with its derivation and falsifiers; and what you would refuse and why. Update the memory notes as you go.

REVIEWER'S CRITIQUE:
<paste here>
