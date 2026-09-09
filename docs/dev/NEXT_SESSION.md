# NEXT SESSION — PHASES 1 AND 2 OF THE BOTH-STRANDED LOCUS AND THE sj+terminus RULE ARE LANDED; THE SHIP PROTOCOL IS NEXT (handoff, 2026-09-08)

⭐⭐⭐ **THE REFERENCE IS `docs/dev/AMBIG_DESIGN.md`** — §3 the design; §4a phase 0; §4b phase 1 as built and
measured; §4c the owner's two corrections measured; §4d phase 2. The rulings of record are `DESIGN.md`
§6b.13. Read them whole. This file is the state.

## WHAT STANDS (2026-09-08, end of session; the working tree is UNCOMMITTED — the owner drives commits)

1. **The base is committed**: `8a2c68e2` on `message-layer`. Everything below is in the working tree.
2. **LANDED in `src/` today, in this order, each gated and A/B'd before it went in:**
   * THE RNA LEVEL LANES (phase 1, steps 2–4): per-strand faces from the flag bits with the intron test
     per strand (`StepContext.exon_pos`/`exon_neg`); two-sided only between an intron and its own
     boundary (counting-only price); the sources; delivery at both-stranded nodes on ψ's cube
     (`PsiMessage.cube_rows`). The ladder: every non-zero row of both halves won.
   * THE FLUX AS A PRICED ESTIMATE (the owner's correction): the junction's route rate is an estimate of
     the exon's RNA abundance priced by the junction–exon pair's disagreement in the STRAND's abundance
     (`count_price`, the count read on the column κ points to, `read_column`), kept lower-sided
     (`flux_level(..., v)`). The ladder: stranded worst 1.000×, unstranded worst 1.011×.
   * THE CEILING AT SINGLE-STRAND NODES (phase 2): an RNA level of the live strand read as a ceiling on
     the gDNA share (`rna_row_of_level`, `_ceilings`), ONLY from a face that sent no composition — a
     licensed face's splice-in map already carries the flux as its cap. Test chromosome: stranded 14/16
     on every panel (worst 1.005×), unstranded 6–7/8 (worst in scope 1.000×), the unstranded zero control
     −15 % (−54 % at pass zero). Ladder: see §4d's line.
   * THE sj+terminus RULE (`DESIGN.md` §6b.14): the orientation helper reads the terminus alone, the
     junction's flux is placed at its exon's flank (`junction_exon_side`), rule 8's mode prediction keeps
     the crossing's scaling. Neutral in scope (the case is 0.6 % of a row at the counting floor of 10–40
     crossings), −14 % locally on the unstranded zero row, every ladder row within 0.3 %. ⚠ The first
     census over-sized it by counting the outside intron's error, which is the intron class's own.
   * THE sj+terminus BLOCK on the test chromosome (`sjterm` · `capsjterm`, 12 genes, RUNX1's and LARGE1's
     patterns on both strands; 205 genes, 6.401 Mb, budget 1,030 k); all seven panels re-simulated,
     cached and certified 30/30 (one transient oracle-validation failure on the sparse panel when three
     panels built at once; a resumable retry passed). The superseded derived set is at
     `~/Downloads/rigel_runs/test_reference_STALE_193genes_2026-09-08/`.
   **Suite 3,781 / 8 xfail, +16 gates re-derived (`CLAUDE.md`).** The landed source reproduces each
   prototype arm to the fragment on the conditions checked.
3. **REFUSED today, each recorded with its numbers** (`ISSUES.md`): "levels always travel" for the gDNA
   lane (twice: alone, and re-judged on top of the ceiling — 1.59× on `g25 ss.50 OFF`); the two-sided flux
   estimate (over-claims at a probe cliff); the total abundance as the flux's witness (charges the other
   transcript's RNA at an overlap); a both-stranded node emitting its determined gDNA level (ladder
   `g05` rows 1.4–1.7×); the naive ceiling that reads every level and every flux (counts the flux twice,
   the weak-κ zero control 36×).
4. **The prototypes and their instruments** are in this session's scratchpad and copied to
   `~/Downloads/rigel_runs/prototypes/2026-09-08_phase1/` (`lanes_proto.py`, `p1b_proto.py`, `p2_proto.py`,
   `lat2_proto.py`, `pass0b.py` — pass zero beside the pipeline — `halves_pass0.py`, `join_halves.py`, the
   gate scripts, the chain dumps). ⛔ They subclass `TransferPolicy`; compare `src` against `src` from
   here. Worth promoting to `scripts/design/`: `pass0b.py` (+4 collected cases, a docstring the gates
   accept).
5. **Before committing**: `preflight.py` green; `ruff` clean; the goldens unmoved (the shipped default is
   still the relay); no upward import.

## WHAT THE MEASUREMENTS EXPOSED AND LEFT OPEN

* **A probed junction over-reads the exon body 20–40×, and at a gDNA-rich exon no local witness sees it**
  (`ISSUES: flux-floor-dispersion`). Phase 1's flux at both-stranded exons pays this on the junction
  panel (`capspan_eq_H`'s host exon 0.988 → 0.640 at `g98 ss.99 ON`, ~175 fragments; the row still wins).
  The exposure of the single-strand ceiling is the deferred stratum (capture-OFF has no over-read; the
  strand channel defends stranded single-strand exons). What would see it: a witness of the junction's
  own enrichment — the owner's owed transport-dispersion decomposition (`transport_dispersion.py`).
* **The rung-2 face map already carries the flux as a cap** (`face_map_lambda`'s ``s``). Any RNA-side
  message at a licensed exon must not add the flux again. This is why phase 2 reaches only faces without
  a composition, and why "levels always travel" cannot be rescued by the ceiling: the floors it adds land
  on licensed exons, where the ceiling may not go. The remaining unstranded plateau at licensed exons is
  `ISSUES: two-sided-exon-row`'s (the intron row's missing upper side), not a lane's.
* **The flux estimate's scatter beyond counting where the pair agrees** (stage 0: 5–9 % at depth) is
  unpriced; a lucky over-read is a sharp floor a few points too high (`capspan_eq_H`'s host exon 0.586
  against 0.645). Same instrument as above.

## WHAT IS NEXT — the owner's call, in the tracker's order (`MESSAGE_RUNGS.md`)

1. **The ship protocol**: the default flip, `silent`/`relay` retired, the obsolescence pass (the relay's
   certified-flux stream is now inside the transfer policy as the flux level; `rna_anchor.py`'s pooled
   fits retire with it), goldens, `MANUAL.md`, `CLAUDE.md`.
2. **The landscape**: the estimator at bound-only nodes and the training population
   (`ISSUES: gdna-landscape-trains-on-false-positives`); phase 2 already cleaned the zero rows' pass-zero
   error by half, which is that population.
3. **The remaining both-stranded structures** (`div`; the antisense's nascent variant on `asin`) and the
   tilt ruling (phase 5: the `theta` channel retired from the contract; today it is unused by construction).

## THE LESSONS THIS SESSION PAID FOR

* **Check the closest refusal before extending a mechanism.** The flux-as-level refusal named the junction
  panel; measuring there first turned up both the 20–40× over-read and the face map's cap, and the
  no-double-count rule came from the second.
* **A witness must be read in the frame the claim is in.** The junction rate is in transcript-strand terms;
  a genome-strand count compared with it must be the column that strand's RNA reads on (κ ≈ 0.01 here).
* **The identity hop charges counting alone; the total is the wrong witness for a strand's lane at an
  overlap; a one-sided bound at a channel-free node is a tilt.** Three prices, each measured before it was
  believed.
* **Never `sed -i` a log a running process writes; never `sleep` in a foreground command; a waiter's own
  command line must not contain the pattern it waits on (it matched itself for half an hour).**
* **Judge a case at its destination class, not its neighbourhood.** The sj+terminus census summed the
  flanks' error and read the intron's ordinary error as the case's; at the boundary and the inside exon
  the case was already at the counting floor.
