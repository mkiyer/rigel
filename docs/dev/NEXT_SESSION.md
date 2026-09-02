# NEXT SESSION — FINISH the `transfer` policy (then, and only then, the default flip)

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand and how to start, not what
    is settled — rulings are `DESIGN.md`, the ranked list is `ROADMAP.md`, the open problems are
    `ISSUES.md`. MOVE anything that settles into those and DELETE this file.

## Where the thread stands (2026-09-02, branch `message-layer`)

Three rungs of the ground-up message rebuild are SHIPPED as `message_policy = "transfer"`
(`bbdf067f`, `3c73f0a6`, `71c3b668`): intron→boundary composition transfer, the face-composed
exon transfer (monotone map, flux ceiling, derived ingredient width), and the intergenic|exon
edge's sign-certified lower bound (zero edges VACUOUS BY LAW). Ladder, both bars, never pooled:
**unstranded 7/8 wins (worst 1.01×), stranded 7/8 (worst 1.00×)** — vs the shipped relay's 5/8
(1.49×) and 2/8 (1.65×) — and best-or-near-best on BOTH adversarial probe panels, where the
relay's anchor is catastrophic. Suite baseline: `CLAUDE.md`. Thread record with every stage,
refutation and accepted error: `docs/dev/COMPOSITION_TRANSFER_STAGE01.md`.

⭐ **The owner's ruling: `transfer` will easily become the default, but an INCOMPLETE policy does
not ship. Finish it first — audit the remaining holes and tackle them one by one.**

## THE HOLE AUDIT (2026-09-02, ladder chain + `g00 ss.50` truth for the mass column)

| hole | population | crossing/slot mass | state |
|---|---|---|---|
| ⭐⭐⭐ **exon\|exon boundaries** | **12,811 slots** | **~2.6 M** (mature RNA crosses these contiguously — 14× the exon\|intron crossing mass) | NO rung touches them: unsolved as destinations, and 8,779 exons (~1.6 M mass) have ONLY such faces — fully unreached |
| ⭐⭐ terminus-refused intron\|exon faces (internal TSS/TES) | ~3,800 exons | ~600 k | rung 2's licence refuses (an UNMEASURED population change); nothing delivered through that face |
| strand-change faces (AMBIG ↔ single-strand) | unmeasured — census it | — | licence refuses |
| the AMBIG tilt channel | ladder AMBIG slots | — | `transfer` delivers λ rows only; θ never imputed |
| spliced lanes at refused faces | — | — | flux info enters only through rung 2's map (`s`); refused faces lose it. The relay's certified-flux anchor is DEAD (cross-locale assumption refuted on the probe panels) — any flux use must be face-local |
| chain-end exons ("no faces") | 92 exons | ~39 k | small; note, don't chase |

⛔ **THE SUBSTRATE GATES EVERYTHING: the twin block has ZERO exon|exon boundaries and ZERO
internal termini.** The owner authors the indicting structures first (`TESTING.md` §0a — the
paradigm: grow one structure at a time). The natural addition: a multi-isoform block (two
isoforms sharing exons with offset boundaries → exon|exon REGIONs and boundaries; an isoform
with an internal TSS/TES → terminus-flagged intron|exon faces). Mind the recorded rebuild
recipe: GTF+abundances+probes edits → `build_test_reference.py` → index → all SEVEN panels
simulate+cache → g00 prewarm + `_main` copy → `calibration_oracle.py` certify (the one-command
chain from this session lives in the scratchpad's `rebuild_all_test_panels.sh` pattern —
rewrite it, it died with the session).

## Derivation notes banked for the exon|exon rung (verify, don't trust)

* An exon|exon boundary's crossing population = {gDNA, RNA± INCLUDING mature} — the SAME set as
  its flanking exons when no terminus intervenes and strand sets match. The shared-population
  licence may therefore extend to exon↔exon|exon-boundary hops directly (composition transfer,
  no new currency). The open questions: neither flank has an intron factory (what carries the
  split evidence INTO the pair?), and chains of exon|exon hops raise the one-hop-vs-multi-hop
  question the foundation spec reserves.
* For internal-terminus faces: the population changes by an UNMEASURED amount (the licence is
  right); what survives is the gDNA side — the edge rung's sign-certified profile-bound pattern
  (`edge_bound_row`) may generalize (gDNA is continuous across a terminus), with the same
  vacuity law.

## Standing cautions

* DERIVE → DESIGN → PLAN → PROTOTYPE → A/B → only then `src/`; ONE mechanism per rung;
  fail-first gates and WATCH each perturbation fire (three gates this thread only earned their
  teeth after a watched non-firing).
* Never pool the two halves; develop on the test chromosome, believe nothing before the ladder;
  run the ADVERSARIAL PROBE PANELS in every rung's loop (`scenarios_probes_sparse` /
  `_junction` — they killed the anchor; they keep the next mechanism honest).
* Prototype via the `message_policy="message"` seam (patch `calibrate.MessagePolicy`
  in-process); compare `src`-vs-`src` across a landing, never harness-vs-`src` (a measured
  1-ULP/flat-posterior divergence family produced fake, flip-insensitive wins twice).
* A claim below its own evidence is SILENCE, never a near-zero row (the vacuity law; near-zero
  rows perturb flat posteriors through refit amplification).
* Accepted errors — do not re-litigate without new evidence: the zero-gDNA edge residual (the
  enrichment-ceiling upper side is owner-REFUSED as over-engineering), the deferred
  `g05 ss.50 ON` row.
* `ISSUES: gdna-landscape-trains-on-false-positives` is the exposed systemic issue (100 %
  fiction training at `g00 ss.50`; the naive exclusion REFUTED — it starves the bootstrap that
  generalizes message-delivered truth). Related, separate; the owner ranks it.

## After the policy is finished (the recorded order)

1. The 0.8.0-METRIC pricing under `transfer` — `calibration_vs_oracle.py` per stratum and
   `solvability_audit.py` (its `calib` column answers whether the rebuilt policy's declared
   precision is EARNED — the thread's original indictment).
2. The DEFAULT FLIP decision with the full trade table (the relay's remaining `g00 ss.50` lead
   is real and rests on the probe-fragile level transport), then the flip protocol:
   `preflight --full`, instruments-not-only-suite.
