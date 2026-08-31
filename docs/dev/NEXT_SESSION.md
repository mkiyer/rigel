# NEXT SESSION — the state after the 2026-08-31 frame migration

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md`. MOVE anything that settles into the permanent docs and DELETE this file.

## What landed today (three commits; the changelog is git)

1. **The RNA length-law plan was refuted** — `rna_pmf` is sound as shipped; the −18k "fix" was a
   compensation lever on the EM's isoform misassignment. `docs/ISSUES.md` (CLOSED,
   `the-rna-length-law-fix`).
2. **The docs were reorganized** (owner decision): `ROADMAP.md` is the short ranked view,
   `docs/ISSUES.md` the issue log, git the changelog.
3. **The frame ruling landed, wave 1** (owner: instruments monitor the DRAINED frame).
   `scan_cache.calibration_inputs` drains at the production seed, runs the production fl models
   (two-pool contrast included — a third divergence closed), and publishes the lift box;
   `OracleTruth.from_cached_parts`/`lift_drain_parts` put cached oracle partitions in the same frame;
   the exact-zeros gate is frame-aware (pass-one refuses, drained RECORDS `gdna_spliced_leak`).
   All ten `calibration_inputs` consumers migrated; **`slot_truth` re-certified on both panels in the
   drained frame** (all gates green, `rna-strands-close` exact). Gates: `test_oracle.py`,
   `test_scan_cache.py`, each perturbed.

## What the next session should know before measuring anything

* ✅ **The in-scope baselines are RE-DERIVED (2026-08-31), drained frame, same-session floors** —
  `docs/dev/BASELINE_2026-08-31_drained_frame.md` is the snapshot; run the instrument for a current
  number. The re-derived facts to know: worst IN-SCOPE scenario is `g98 ss.99 capture-ON` on the
  0.8.0 metric (gDNA UNDER-called at near-pure objects — the pattern holds on every stratum) and
  `g50 ss.50 capture-OFF` on solvability's Σ|err|; the relay hurts the solvable set on 9/12
  contaminated conditions; MessagePolicy ≡ silent on all 46 in the new frame; the reseed floor is
  ~30–2,700 transcript fragments per condition. ⛔ Nothing measured before 2026-08-31 is comparable.
* ⛔ **Only wave-3 bank-readers still read pass one** (list in
  `ISSUES: instruments-calibrate-undrained-cache`): the census family, `structural_claims_audit`,
  `calibration_truth_ab`, `fl_pool_purity` and friends — migrate as touched. Everything reached
  through `calibration_inputs`, `pass0_vs_oracle.measure_condition` or `prior_vs_oracle` is drained.
* ⭐ **The owner's priority: the drain's certified-RNA leak** (`ISSUES: drain-contaminates-certified-rna`,
  elevated 2026-08-31 — "avoid all sources of false positives"). Up to 1.9 % of the certified-RNA
  channel at `g98 ss.99 ON` is drain-assigned gDNA. The issue carries three candidate repairs and the
  two-sided scoring requirement (false positives at zero-RNA AND the spliced-pool repair the drain
  provides — a fix must not resurrect the −4 bp bias the drain currently prevents). DERIVE, don't
  tune; the provenance-split candidate needs the owner (it touches "determinacy, not provenance").
* ⚠ `ISSUES: rename-the-drain` is filed (owner: "might consider") — ~530 sites, needs its own staged
  census pass; candidate verbs and their collisions are in the entry. Do not rename ad hoc.

## The other live thread

⭐⭐ **The message policy — the standing charter, still at rung 0.** Bar unchanged: win on unstranded,
minimal harm on stranded, never pooled. First answer `ISSUES: message-value-for-blind-slots` with no
solver; the anchored twin block was built for it. ⚠ Re-baseline first (see above).
