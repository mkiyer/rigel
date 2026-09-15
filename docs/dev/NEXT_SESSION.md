# NEXT SESSION — start here (2026-09-14, after the test-chromosome session)

This file is only how to begin. The port is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F; the θ thread's
derivations are `docs/dev/THETA_QUADRATURE.md`.

## The agreed order of the sessions (owner, 2026-09-14)

1. ~~Code review and cleanup~~ — DONE 2026-09-14.
2. ~~The test chromosome's new structures, both panels remeasured~~ — DONE 2026-09-14, the six commits after
   `1f885a7e` (one per structure, then the remeasured docs).
3. **THIS ONE — the ruler at zero gDNA** (`ISSUES: g00-shrinkage-upstream-repair`) and **the flux price's
   witness in whole-strand units** (`ISSUES: flux-price-witness-units`). One mechanism at a time, each its
   own commit; the defects the new structures exposed are listed below and in `SESSION3_DEFECTS.md`.
4. The perform## What this session left

Six commits on `main` after `1f885a7e`: the five structures, one each (`enc`, `inexon`, `shared`, `inintron`,
`div`), then the remeasured docs. Each structure's commit message carries its ladder mirror, its geometry,
its budget step and its census on certified truth. Nothing in `src/` changed: `rename_identity.py --check`
was BIT-IDENTICAL on all three `review_identity_*` references at the start and after every structure
(18/18), and the suite stayed at 3,374 passed / 3 xfail / 3,377 collected at every stage (a content-only
change moves the count by zero). The per-structure snapshots that made the commits stay in the session
scratchpad's `commits/` as the record.

nt-only change moves the count by
zero).

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, at the six test-chromosome commits after `1f885a7e`, or later). Read `CLAUDE.md`,
> then `docs/dev/NEXT_SESSION.md`. Run `scripts/design/preflight.py --full` and the suite before touching
> anything; `CLAUDE.md`'s baseline line is the count to reproduce (3,374 passed / 3 xfail / 3,377 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code; converge
> and delete — no legacy, no compatibility shims, no speculative code; the suite's count re-derived, never
> adjusted; each item its own commit, committed on my go; never `ruff format scripts/`.
>
> THIS SESSION IS THE RULER AT ZERO gDNA (`ISSUES: g00-shrinkage-upstream-repair`) and THE FLUX PRICE'S
> WITNESS IN WHOLE-STRAND UNITS (`ISSUES: flux-price-witness-units`), one mechanism at a time, each its own
> commit: DERIVE → DESIGN → PLAN → PROTOTYPE outside `src/` → A/B against the shipped policy on the test
> chromosome (`policy_prototype.py`, `calibration_vs_oracle.py --set`) → confirm on the ladder → only then
> `src/`. Judge every mechanism per stratum on `calibration_vs_oracle.py` and `policy_benchmark.py --by-class`
> on BOTH panels against the standing numbers in `DESIGN.md` §7, the two zero controls and the deferred
> stratum reported apart; a falsification test first, verified failing, then break the fixed code and watch
> each gate fire. `rename_identity.py --check` against the three `review_identity_*` references gates every
> `src/` change; re-freeze them only with a stated reason. The defects the new structures exposed are listed
> in this file — take the one your mechanism touches, never all of them. End with the handoff for session 4.

## The new standing numbers (both panels, this tree)

`DESIGN.md` §7 carries the table (the tree at `1f885a7e` on the twelve-block chromosome; the ladder is
unchanged and reproduces the location-floor landing's zero controls, 282 / 194 / 265 / 172 = 913, exactly).
In one line each, `policy_benchmark.py` silent → transfer and `calibration_vs_oracle.py` P/O gDNA:

* ladder: unstranded OFF 358,551 → 307,281 (0.9938); stranded OFF 292,673 → 248,975 (0.9949); stranded ON
  598,645 → 427,003 (0.9954); deferred 18.8 M → 3.58 M (0.8546); g00 367 / 258 / 353 / 233; ruler at g00 0.154.
* test chromosome (265 genes): unstranded OFF 54,436 → 51,230 (1.0088); stranded OFF 47,114 → 43,710 (1.0071);
  stranded ON 64,200 → 41,866 (1.0028); deferred 1.93 M → 194,992 (1.0220); ss 0.70 rows 150,048 → 113,062;
  g00 8,882 / 28 / 8,812 / 96 / 9,040 / 31, of which 8,873 / 8,799 / 9,032 are the shadow floor on
  `test_blank` (the annotated chromosome's own zero rows read 8–13 OFF, 28–96 ON); ruler at g00 0.141.
* where `transfer`'s error sits (`--by-class`): ladder unstranded OFF — introns 45 %, `exon|intron` 14 %,
  `exon|exon [term]` 13 %; ladder stranded ON — `exon|exon [term]` 27 %, `exon|exon` 18 %, walled exons 18 %,
  `exon|intron` 15 %; test chromosome stranded ON — licensed-face exons 49 %, walled exons 23 %; test
  chromosome unstranded OFF — the shadow floor 39 %, licensed-face exons 24 %, introns 22 %.

The full outputs: the session scratchpad's `final/` (`cvo_test.txt`, `cvo_ladder.txt`, `pb_test.txt`,
`pb_ladder.txt`, `cvo_test.json`) and `gen/pb_strata.py`, which sums a `--by-class` output per stratum.
All seven test-chromosome configs were re-simulated, re-cached and re-certified on the final chromosome
(the fl-gap and od05 arms' caches are current; their FIELD stamps are per row as before).

## What the new structures exposed — for session 3

Nothing the five structures exposed is a calibration defect that blocks 0.8.0; the certified truth and
the worst-object reads are in `SESSION3_DEFECTS.md` in the scratchpad. In brief:

* **The FIELD gate flakes on boundaries at λ ≈ 7.** Two of the five benign-panel builds (and the odg05
  rebuild) certified one row at COMPOSITION only: `gdna-field-uniformity` flags a slot at |z| > 4 and fails
  a row at a flag rate ≥ 1e-3, i.e. three flags among ~2,800 z-scored slots, and every flagged slot was a
  boundary with n ≈ 18 against λ ≈ 7 — once the terminus cluster's ten boundaries at 4,991,853–4,991,874,
  which share one crossing population within 21 bp and so flag together. Under Poisson ~0.3 such boundaries
  are expected per build against 3 and 11 seen: either the bar is too tight for a 4σ per-slot flag or the
  unspliced crossing count at a low-λ boundary is overdispersed against `rho × E_g`. Nothing the metric reads
  needs FIELD. An owner's call on the certifier, not a mechanism for session 3.
* **A balanced both-strand exon with no single-strand piece reads gDNA locally** — the shared exon 0.990 and
  the head-to-head overlap 0.983 at the eq loci — and the RNA level lanes bring both to within 0.03 of the
  truth (0.130 / 0.121; 0.123 / 0.148). The stress passes where a junction pins each strand; the in-exon
  block's slot (no junction anywhere) is where the tilt atom's accepted limit now has a panel number, read
  it off `policy_benchmark.py --by-class` on the walled-exon class of the stranded rows.
* **Unstranded rows: the new single-strand exon pieces read the ½ default plus a message** (0.35 against 0.23,
  0.27 against 0.23, 0.30 against 0.20) — the in-scope unstranded residual the ruler and the flux witness
  are meant to move; the encompassing block's intron∩antisense regions read within 0.02 on 1,000–1,600
  fragments.
* **Two instrument hazards met on the way:** `panel.py cache` for two panels with the same condition names
  cannot run side by side (the origin split's work dir is keyed by condition under `$RIGEL_SCRATCH`, default
  `/tmp`; the rule is now in `TESTING.md` §0a's recipe), and `panel.py status` ends with `next: panel.py
  score`, which is the quant-accuracy arms and not part of the calibration loop.

## The session scratchpad (persists; nothing in the tree cites it)

`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/23defdd7-4a51-4398-89bc-ac9b16e623dc/scratchpad/`:
`commits/` (the six snapshots, messages and `commit_series.sh`), `gen/` (`add_block.py` — the block
generator, one structure per invocation, `rebuild.sh`, `slot_census.py` — the new slots by node class off
`slot_truth.npz`, `--field-flags` for the FIELD gate's flagged slots; `final_measure2.sh`), `census/` (the
ladder-annotation census behind every "mirrored from" claim: `ladder_census.txt`, `expressed.py`),
`identity/` (18 logs), `rebuild/` (every build's logs), `worst/` (worst_objects.py per build), `final/`
(the standing-number runs and the six configs' logs), `handoff/SESSION3_DEFECTS.md` (the long form of the
list above). The pre-session derived set is parked at `~/Downloads/rigel_runs/test_reference_STALE_205genes_2026-09-14/`
(29 GB; delete once the new panels are trusted); the four intermediate builds were deleted, their logs kept
under `rebuild/`.

## The two pricing questions inside the port's unit — what the owner decides (carried from 2026-09-14)

Both live in `messages/transfer_rows.py` and `messages/lanes.py`, the port's first step, so a fix after the
port is written twice; the decision is fix-before or accept-and-record.

* **`ISSUES: flux-price-witness-units`.** A junction's certified flux is a lower bound on the exon's RNA
  level; its price compares the junction's route rate (whole-strand units) with the exon's own count on the
  read column, which holds only `(1 − κ')` of the strand's RNA plus half the gDNA, so every flux ceiling pays
  `log(1 − κ')²` nats² of spurious disagreement — 0.48 on unstranded data, an in-scope stratum. The naive fix
  was REFUSED; the owed form is the strand's RNA count in whole-strand units. Session 3 derives and A/Bs it.
* **`ISSUES: the-lower-bound-noise-ratchet`** (the encompassing flank's xfail; now also the `enc` block's
  single-strand flanks on the panel, 375–388 fragments each). Recommendation: accept and record; the port
  carves the level rule as it is.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere (the tilt count is derived); the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot.
* A few high-quality instruments, kept current; no suite gate polices instruments (owner, 2026-09-14).
* The source cites no doc (owner, 2026-09-14).
* `EMConfig.warm_start`'s `prior` and `uniform` arms stay; the `gdna_none_` condition names stay.
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's; parallelism
  waits for the port.
* CI runs on demand only (`.github/workflows/ci.yml`, owner, 2026-09-14); the identity references are
  `~/Downloads/rigel_runs/arms/review_identity_*` (frozen on `84923136`); `ISSUES: hygiene-ledger` holds the
  review's leftovers.
