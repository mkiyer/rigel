# THE MESSAGE RUNGS — the tracked list of messages and boundary cases still to finish (owner, 2026-09-02)

    ⚠ A DEV DOC and a TRACKER. It says what is done, what is next and in which order; the
    derivations and measurements live in `COMPOSITION_TRANSFER_STAGE01.md`, the substrate in
    `test_chr.yaml`'s header. Move a finished item's verdict to its permanent home and mark it here.

**The paradigm (owner rulings 2026-09-01/02):** one node type, one message, one boundary case at a
time; DERIVE → DESIGN → PLAN → PROTOTYPE (outside `src/`) → A/B on the test chromosome against
per-object truth, with the perturbation watched firing → confirm on the ladder, two halves apart,
adversarial probe panels in the loop → only then `src/`. A moved number must have ONE cause, so the
substrate grows ONE structure per step and a step is a whole session or more, validation included.
⛔ Nothing below is small; do not bundle two items to save a rebuild.

## Where the rungs stand

| rung | substrate | what it is | state |
|---|---|---|---|
| 1 | twin block | multi-exonic, single-isoform, single-stranded: the intron\|exon BOUNDARY | ✅ COMPLETE 2026-09-02 — intron → boundary (rung 1) and exon → boundary (item 1) both ship |
| 2 | twin block | the same substrate: the EXON region | ⚠ PARTIAL — the boundary → exon face map ships; both faces sum; the exon's reverse messages are owed (see rung 1) |
| 3 | mono block | single-exon transcripts: the intergenic\|exon EDGE | ✅ the sign-certified lower bound ships; the ceiling REFUSED (accepted error) |
| 4 | isoform block | multi-exonic, MULTI-isoform, single-stranded: exon\|exon boundaries | ⏳ substrate holds ONE structure (`altstart`); nothing ships; a jumped-ahead prototype is recorded and set aside |
| 5 | — | strand-change faces, both-stranded loci, the AMBIG tilt channel | ⛔ LAST, owner ruling |

## The ordered list

Finishing rungs 1 and 2 (owner: high priority, first):

| # | message / case | what must be derived | state |
|---|---|---|---|
| 1 | **exon → intron\|exon boundary** | DERIVED 2026-09-02 (owner's rescale/subtract/rescale ⇒ `f_b = f_E·(U_b+S_b)/U_b`, the enrichment ratio cancels; = the shipped splice-in map read backwards; checked unbiased on certified truth). Owed: the honest width (delta method through the map, diverging at `U_b → 0`), the imputation-cost premise priced on the adversarial panels, prototype + falsifiers, A/B | ⏳ PROTOTYPED 2026-09-02: wins every stranded/part-stranded capture-ON row (0.83–0.99×), exact silence elsewhere, both falsifiers fire; ONE benign cost (`g25 ss.99 ON` 1.014×, nascent-bearing probed boundaries) and sparse-panel harm (1.124×) = the premise cost, measured. Transfer variance DERIVED fresh: counting (1/S+1/U) + premise (log a, fitted by the two-witness estimator in log-ρ units — reads 0 off capture; under capture a BIAS a≈1.3 benign / 2.2 junction with ~no spread, recorded not corrected); width applied as a MARGINAL over log ρ; opportunities must be capture-blind for both components (geometric form). `formb_g` wins every stranded capture-ON row on all three panels. LADDER 2026-09-02: unstranded byte-identical, stranded capture-ON 0.987–0.995×, capture-OFF +8…+33 fragments, reversal fires. ✅ SHIPPED 2026-09-02 (`messages/transfer.py`, `simplex_logodds.strand_row_logodds`; ruling `DESIGN.md` §6b.4; 3 new gates, 3 perturbations watched firing) |
| 2 | **boundary → intron** | whether an intron with its own factory takes anything from its boundaries; possibly "nothing", measured | ☐ |
| 3 | **the exon solve with every face speaking** | both faces' incoming rows plus the exon's own evidence, once item 1 exists on both sides; the two-witness sum re-priced | ☐ |
| 4 | the factory on a region carrying BOTH exon and intron bits (raised 2026-09-02: today the factory runs only where no exon bit is set) | whether the density-against-background measurement is valid there and under capture | ☐ when first needed |

Rung 4, one structure per step (each step = one YAML change, one licence or message):

| # | boundary case | structure | what must be derived | state |
|---|---|---|---|---|
| 5 | **exon\|exon boundary WITH a terminus, solved from the OUTSIDE flank** | `altstart` (present) | the orientation table (TSS+/TES− body right, TES+/TSS− body left; mixed → no side), the opportunity shift, the solve; verified per object, then with orientation reversed | ☐ after items 1–3 |
| 6 | the region INSIDE the terminus | `altstart` | what, if anything, crosses into it (a level bound was derived and REFUTED on the sparse panel — record, do not re-litigate without a new measurement) | ☐ |
| 7 | exon\|exon boundary at an alternative SPLICE SITE | add `altss` | the crossing loses the transcript that splices out (join vs leave direction); the composed transport was prototyped and measured neutral on the ladder | ☐ |
| 8 | a region WALLED by two termini | add `nest` | what reaches it; measured 2026-09-02: the refit prior already serves it (0.666 vs 0.630) | ☐ |
| 9 | an alternative FIRST exon inside an intron (`exon\|intron[term]`) | add `instart` | the intron is always the outside flank; the inside exon's other face | ☐ |
| 10 | terminus faces with sj+term flags, chains of termini | (ladder only so far) | after 5–9 | ☐ |

Rung 5 (last): strand-change faces (246 exon|exon + 888 exon|intron on the ladder), the both-stranded
locus, the AMBIG tilt channel.

## Recorded and set aside (do not rebuild)

* The rung-4 prototype (`rung4_proto.py`, session scratchpad 2026-09-02): composed join-only transport
  NEUTRAL on the ladder (worst 1.006×/1.008×, −6 % at `g00 ss.50 OFF`); the inside bound REFUTED
  (in-scope +1.3 %, sparse probes +8 %); flips fire. Its pieces re-enter only as items 5–8 earn them.
* The census instruments (`exon_exon_census.py`, `terminus_pair_gap.py`, `directed_reach_census.py`)
  live in the same scratchpad; promotion into `scripts/design/` with `--self-test` is owed when item 5
  starts.
* `panel.py cache` lacks the g00 pre-warm / `_main` copy / certify steps (recipe in `TESTING.md` §0a).
