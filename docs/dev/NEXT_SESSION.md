# NEXT SESSION — the state after the 2026-08-31 RNA-length-law refutation

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md`. MOVE anything that settles into the permanent docs and DELETE this file.

## Where the tool is

**The fragment-length thread is CLOSED, both halves.** The gDNA side landed last session (two estimands,
routed and coupled — `ROADMAP.md`'s fragment-length line). The RNA side closed this session as a
**refutation**: `rna_pmf` is sound as shipped, the planned fix is not licensed, and the full record with
its stamped measurements is `docs/ISSUES.md` (CLOSED, `the-rna-length-law-fix`). The strand
overdispersion thread was already closed. ⭐ **Fragment lengths and strand estimation are no longer open
work items** — tool development moves on to `ROADMAP.md`'s NEXT section (the dissection loop, `ISSUES:
message-value-for-blind-slots`, and `ISSUES: per-transcript-prior-lane`).

## What the refutation established, in three lines

* The plan's −4 bp sj-observability bias was measured in the **undrained** frame; production fits the
  **drained** payload, where the residual is −0.2…−1.3 bp. The drain repairs the selection. Step 0
  ("re-measure before building") is what caught it — the whole plan would otherwise have been built.
* The −18k transcript "ceiling" was a **compensation lever**: shifting the shipped pmf +10 bp — far past
  truth — wins −171k at `g05 ss.99 ON` (and −43k of 542k at `g50 ss.99 OFF`) by nudging the EM's standing
  within-gene isoform misassignment at short-exon genes through the per-locus prior. On the 0.8.0 metric
  (`calibration_vs_oracle.py`) a perfect RNA pmf is worth ±0.3–2.9 %, mixed sign — nothing.
* ⛔ **Do not rank any calibration-input change on the transcript thermometer** — it rewards this
  compensation at up to 7 % of transcript error per condition. The 0.8.0 metric is the judge; the
  thermometer stays a thermometer. (Last session's "score on the metric its consumer reads" caution,
  now with the sharper mechanism: the lever's dissection is in the closed ISSUES entry.)

## Observations for the owner, found in passing (nothing changed)

* **The calibration instruments score a different substrate than production calibrates.**
  `calibration_vs_oracle.py`, `calibration_oracle.py` and `calibration_walk.py` run `calibrate()` on the
  UNDRAINED scan cache (`scan_cache.calibration_inputs`), while `pipeline.py` drains the side buffer
  BEFORE `build_fl_models`/`calibrate`. The cache stores pass one by design (so the drain can re-run at
  any seed) — the instruments could drain after loading (`pipeline._drain_side_buffer`, one call). The
  drained-vs-undrained gap is exactly what produced the RNA-law misdiagnosis, so the gap in the
  release-metric instrument deserves a decision: price it, then either drain there or record why not.
* **A stale comment**: `pipeline.py:416-419` says the drain's fl models are "the SAME de-tilted RNA pool
  the calibrator will read" — false since the calibrator reads the drained tally (line 382 has it right).
* **The dose–response ammunition for `ISSUES: per-transcript-prior-lane`**: the −171k/+10 bp lever is a measured lower bound on what
  the EM's isoform allocation at short-exon multi-exon genes is worth; it concentrates on expressed,
  multi-exon, median-exon ≤150 bp transcripts and is capture-independent. If that issue's sparsity mechanism
  needs a target population, this is it.
* `r_hat` (the two-pool contrast's discarded contaminant law) is accurate off capture (±0.2 bp vs
  nascent truth) and broken under it — kept as a fact in the closed ISSUES entry, licensing nothing.

## The other live thread

⭐⭐ **The message policy — the standing charter, still at rung 0.** `MessagePolicy` remains byte-identical
to `SilentPolicy` on all 30 test-chromosome conditions. The bar is unchanged: **win on unstranded, minimal
harm on stranded, never pooled.** Before building any mechanism, answer `ISSUES: message-value-for-blind-slots`
with no solver — is the information a blind unstranded or AMBIG slot actually present in its neighbours?
The anchored twin block was designed for exactly that. ⚠ Re-baseline first — the policy numbers in
`ROADMAP.md` predate the strand-estimator and fragment-length changes.

## Standing risks carried forward (unchanged from last session)

* **The gDNA contrast survives capture by a degeneracy, not by its premise** — a probe panel that put RNA
  back into intergenic space would break it silently — `ISSUES: capture-degeneracy-standing-risk`.
* **Geometry is knife-edge sensitive to the pmf tail** — any change to an `effective_lengths_em` input
  needs the reseed floor beside it.
* **Nothing has been run on real data** (`real-data-is-a-test-input`) — no fragment-length claim has met
  a real library.
