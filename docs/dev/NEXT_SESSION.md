# NEXT SESSION — PHASE 1 OF THE BOTH-STRANDED LOCUS IS APPROVED; BUILD IT (handoff, 2026-09-08)

⭐⭐⭐ **THE REFERENCE IS `docs/dev/AMBIG_DESIGN.md`** — the design, the census that sized it, phase 0's
findings, the owner's rulings (2026-09-06/07/08) and the phase table. Read it whole before touching
anything. This file is the state; that file is the plan.

## WHAT STANDS (2026-09-08)

1. **The level lane is landed** (`DESIGN.md` §6b.12; uncommitted): the gDNA level as an absolute
   profile, the default rule of every face without a composition rule, lower-only, bounds intersect;
   ladder unstranded 8/8 at or below the pre-lane policy, stranded 7/8 (`g50 ss.99 ON` 0.948×,
   `g98 ss.99 ON` 0.908×); suite 3,765 / 8 xfail.
2. **Two blocks joined the test chromosome** (uncommitted): the TERMINUS-CLUSTER block (2026-09-05, the
   empty pieces) and the BOTH-STRANDED block (2026-09-07: `asin`, `asinrev`, `span`, the owner's `conv`;
   24 two-gene loci, host-only `cap` twins). 193 genes, 5.986 Mb, budget 960 k; all seven panels
   re-simulated, cached and certified 30/30; preflight green; the superseded derived sets are under
   `~/Downloads/rigel_runs/test_reference_STALE_*`. ⛔ `git checkout` of the YAML restores the COMMITTED
   chromosome and drops every uncommitted block — never do it; regenerate from the generators if needed.
3. **Phase 0's two findings** (`AMBIG_DESIGN.md` §4a): the both-stranded stretches are unreached at pass
   zero because a held composition is never re-issued as a level (the fix is that levels ALWAYS travel,
   a composition alongside where a map exists); and the spanning host's junctions inside the antisense's
   exon need a TWO-SIDED RNA level from their own intron (rung 1's own law: one shared unspliced
   population). Both are steps of phase 1.
4. **The owner's rulings for phase 1 (2026-09-08)**: the three levels travel TOGETHER in one message
   (components present only where measured; empties forwarded); the certified flux joins phase 1 as an
   RNA source; ONE representation everywhere — profiles on the solve grid, a held level evaluated at the
   density each cell implies, no Gaussian summary anywhere in the transfer policy; the bar is about one
   percent of a row.

## PHASE 1 — the order, each step a falsification gate watched firing, then an A/B against the landed policy, halves apart, pass zero beside the pipeline

1. **Levels always travel.** The gDNA lane emits on every face (a composition alongside where a map
   exists); the solve reads the composition from a side that sent both, the level otherwise. Gate: a
   stretch's entrance forwards; a node holding both from one side reads the composition. Judged on the
   ladder's sixteen rows and the three panels — it changes the gDNA lane's reach everywhere.
2. **The RNA lanes' faces and plumbing, no sources yet.** Per-strand faces from the flag bits (a junction
   of + stops RNA+ and not RNA−); two-sided only between an intron and its own boundary. Byte-identical.
3. **The sources.** (a) own strand profiles read as RNA levels of the live strand, priced by the strand's
   own counts; (b) the certified flux at an exon's junctions as that strand's RNA level at the exon. Gates:
   the coordinate round trip; the hop price by strand counts; no echo; the flux level is a lower bound at
   the route rate's Poisson width, one hop for the spliced claim itself.
4. **Delivery at AMBIG nodes** as per-strand rows evaluated on the ``(λ, θ)`` cube (the analogue of
   `lam_rows`; the relay's Gaussian RNA channels are not used). Gate: THE BRACKET THEOREM on a hand-built
   node — three lower bounds and the strand equation give a two-sided gDNA share, any one removed opens a
   side. Judged first at the spanning loci's boundaries, then everything.

Then phases 2–5 as the design's table has them (single-strand recipients, the remaining structures,
the tilt ruling). Prototype outside `src/` first (`policy_prototype.py --module`), A/B, only then `src/`.

## WHERE THE POLICY STANDS (the level lane UNCOMMITTED in the working tree — the owner drives commits)

Ten composition messages plus THE LEVEL LANE. Every directed face of the chain now carries a
composition rule or the lane, or leads into structural pure gDNA or off the chain — the completion
contract's "no skipped boundary or region" is a GATE, not a hope. What is still owed before the default
flips: the sj+terminus composition (D), the AMBIG complex's RNA levels and tilt (H), the terminus-cluster
block that lets the test chromosome see empty pieces, and the ship protocol. The unstranded FIRST PASS
(the plateau's median at one-sided nodes; the −26 % a two-sided lane would give) is the enrichment
witness's and the landscape's, after the architecture (`ISSUES: two-sided-exon-row`).

## THE NEXT CASES — `MESSAGE_PLAN.md` §6 (owner rulings 2026-09-06: the enrichment witness IS the landscape prior; bound-only nodes do not train it; the rules finish first)

1. **H. THE AMBIG COMPLEX** — ⭐ designed 2026-09-06, `docs/dev/AMBIG_DESIGN.md` (the RNA level lanes, the bracket theorem, the `asin`/`span` blocks, five phases, five owner decisions); the RNA levels per strand and the tilt: 9,912 AMBIG nodes carry 38 % of
   the remaining error on `g50 ss.99 ON` and 50 % on `g98 ss.99 ON` (`ambig_census.py`). Needs a
   both-stranded block on the test chromosome mirroring a real ladder locus (the owner authors).
2. **D. sj+terminus** — 386 boundaries, 0.5–0.7 % of any row; `outside_flank` with a junction present,
   item 5's map with the leaving flux.
3. **The ship protocol** — default flip, `silent`/`relay` retired; the g00 rows are the landscape's.
4. **Then the landscape** — the estimator at bound-only nodes and the training population.

## HOW TO START

* Read `MESSAGE_PLAN.md` first (the ten rules' status, the level rule, the order), then
  `TWO_PHASE_BACKBONE.md` §0 (the words), §3 (the skeleton as landed), §6d–§6f, §8–§9.
* The prototypes and logs are in the 2026-09-04 session scratchpad
  (`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/d7adfcc1-…/scratchpad/`): `bp_proto.py`,
  `flux_proto.py`, `flux_stage0.py`, `bp_identity.py`, `pass0_score.py`, `halves.py` and every `.out`.
  ⚠ Scratchpads are per session: `flux_stage0.py`, `bp_identity.py` and `pass0_score.py` are worth
  promoting into `scripts/design/` (each adds +4 collected cases and needs a docstring the gates accept).
* Every A/B: `policy_prototype.py --panel test|test_junction|test_sparse --all --by-class`, the halves
  apart, pass zero beside the full pipeline, then the ladder. ⛔ Never quote one panel.
* The suite baseline after the landing is in `CLAUDE.md`; re-derive, never adjust.

## OWNER DECISIONS OUTSTANDING FROM BEFORE (recorded, not re-litigated)

Item 2's `g05 ss.99 ON` residue; item 5's +46 at `g05 OFF`; item 6's +1.2 % inside exons at `g50 ON`;
item 7's low-gDNA capture-ON rows on the benign and junction panels and `g50 ss.99 ON` benign/sparse.
The pre-walled panels are parked at `~/Downloads/rigel_runs/test_reference/pre_walled/` (14 GB).

## THE LESSONS THIS SESSION PAID FOR

* **A level carried between locales is refuted by probe placement alone, and the tool never sees the
  probe panel** — `DESIGN.md` §6b.9's founding refusal, re-measured at 4–17× on the adversarial panels
  after a clean stage 0 and a 9/10 win on the exon-probed one. Run all three panels before believing
  any mechanism that reads a rate across a face.
* **A count-based row has no resolution at a three-fragment exon, and the refit prior trains on its
  false modes** (`ISSUES: gdna-landscape-trains-on-false-positives`): the ladder's zero control turned
  a +13 % test-chromosome loss into 3.3×.
* **Today's `prepare` ORDER was an implicit propagation rule** (the walled block's terminus chains):
  a pass replaces it, and the per-slot identity gate is the instrument that finds such rules.
* **The relay's pooled transport-centre fit refuses on small substrates and misfires where it accepts**;
  it is not a drop-in for the transfer policy.
