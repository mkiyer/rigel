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

## 2026-09-02, later: THE SUBSTRATE FOR HOLES ① AND ② IS BUILT (uncommitted — the owner drives commits)

* ⭐ **The test chromosome now has ONE source, `scripts/sim/test_reference/test_chr.yaml`** (the
  `rigel sim` scenario schema + `probed` + `shadow_genes`; owner ruling: leverage the existing YAML
  architecture). `build_test_reference.py` renders the GTFs, abundances and three probe panels
  (versioned beside the YAML) and the FASTA; `--check` and `tests/test_test_reference_renders.py`
  refuse drift; `preflight.py` checks both. The conversion was falsified: from a YAML encoding the old
  chromosome, every render and the FASTA came back byte-identical (abundances as a set).
* ⭐ **The ISOFORM BLOCK** (8 types × 5 blocks, 40 genes) and **the strand rule** (every gene an explicit
  strand, chromosome balanced 42/43, every type on both strands; − geometry mirrored) — design and
  built-index verification in `docs/dev/ISOFORM_BLOCK.md`. Depth raised to 480 k in all seven configs.
* ⭐⭐ **THE LADDER CENSUS RE-FRAMED HOLE ①** (`ISOFORM_BLOCK.md` §1): the dominant exon|exon class is
  the internal TERMINUS (7,604 slots / 950 k), not the splice site (4,838 / 458 k); 10,259 exons
  (1.05 M) are unreachable by any licensed chain, walled by terminus faces. Holes ① and ② share one
  mechanism.
* ⭐⭐ **STAGE 0 IS DONE ON LADDER TRUTH** (§2 there): at a terminus boundary the OUTSIDE flank (the
  side the terminating transcript does not cover — direction read from the flag) shares the crossing's
  composition to within counting noise on every condition; the INSIDE flank does not. The directed
  licence is certified; the inside side is at most the gDNA-continuity bound (sign certified for
  `exon|intron[term]` by the edge's tool-scope assumption, NOT for `exon|exon[term]` under capture).
* The seven-panel rebuild ran from the scratchpad driver (recipe now in `TESTING.md` §0a, including
  the g00 pre-warm / `_main` copy / certify steps `panel.py cache` still lacks — a `panel.py` gap owed).
  ⛔ The new baseline (`policy_benchmark.py --panel test --policies silent transfer`) must be recorded
  BEFORE any rung-4 code; every earlier test-chromosome number is on the 45-transcript substrate.

## THE OWNER'S RESET (2026-09-02, late): back out, one message at a time — `docs/dev/MESSAGE_RUNGS.md` IS THE LIST

The four-structure block and the three-mechanism prototype jumped ahead. Owner rulings: (1) rungs 1 and 2
are NOT finished — the exon → boundary message was nullified and every exon-side message is owed, and
they come FIRST; (2) rung 4 adds ONE structure per step; the substrate was TRIMMED to `altstart` the same
day (55 genes, 27 + / 28 −, 30 probed; rebuild relaunched — `scratchpad/rebuild2.log`, the old 125-tx
derived set at `test_reference_STALE_125tx_2026-09-02`); (3) the tracker `MESSAGE_RUNGS.md` carries the
ordered list, and each item is at least part of a session. The trimmed benign panel is certified and its
silent/transfer baseline is RECORDED (thread record, last section: unstranded 14/20 worst 1.16×, stranded
8/10 worst 1.00×; `g25 ss.50 OFF` at 1.06× is the in-scope row to watch). Item 1 (form B, the owner's arithmetic; the enrichment ratio cancels) is DERIVED, checked on certified truth
and PROTOTYPED (`scratchpad/item1_proto.py`; thread record, the two Item 1 sections): wins every
stranded/part-stranded capture-ON row, exact silence elsewhere, falsifiers fire. The transfer variance was derived FRESH (thread record: counting + premise; the premise measured as a
BIAS a≈1.3 under benign capture with ~no spread, recorded not corrected; the width applied as a MARGINAL over
log ρ; opportunities capture-blind for BOTH components — the geometric form, which is what fixed the sparse
panel). Ladder: unstranded byte-identical, stranded ON 0.987–0.995×, OFF ≤ +33 fragments, reversal fires.
⭐ ITEM 1 SHIPPED 2026-09-02 (owner: "land this"): `transfer.py` gains `splice_out_row` + the exon block on a
shared licensed-face helper (rung 2 rewritten on it), `simplex_logodds.strand_row_logodds` is the public
strand row, `calibrate` passes the fitted strand triple; ruling in `DESIGN.md` §6b.4; 3 gates, 3
perturbations watched firing. NEXT: item 2 (boundary → intron) — the owner may want a cleanup session first. ⭐ START THERE: item 1, the exon → intron|exon
boundary message (the splice-out direction), derived and taught one piece at a time — the owner wants to
understand every piece. ⚠ Correction recorded: in an exon-covered locus the introns and sj still exist;
what is true is only that the intron FACTORY runs on regions without an exon bit (tracker item 4).

## RUNG 4 — where it stood before the reset (2026-09-02)

* DERIVE done (thread record, RUNG 4): the DIRECTED LICENCE (outside flank ↔ terminus boundary by
  composition; inside flank at most a gDNA lower bound), the COMPOSED TRANSPORT (one intron factory row
  through a chain of licensed maps, JOIN steps only — LEAVE steps refused, splice-out priced at 1,741
  ladder exons / 177 k), and the INSIDE BOUND (the row's running maximum through the level map with the
  crossing count's one-sided Poisson tail; `edge_bound_row` is its special case).
* PROTOTYPE built (`rung4_proto.py`, scratchpad) and measured on the rebuilt benign panel: identity arm
  byte-identical to `transfer`; `multi` −9 % on the deferred blind row (`g50 ss.50 ON` 27,892 → 25,445)
  with ≤ 3 fragments of stranded movement, +1–3 % on in-scope unstranded OFF rows (noise-level in
  absolute terms — dissected: a 10-molar gene's gDNA counting noise); the inside bound INERT (every
  inside flank on this chromosome has a licensed witness or is reached by the refit prior — the
  `nest` walls read 0.666 vs truth 0.630 under the SHIPPED policy against silence's 0.001); both
  falsifications fire (`flip` +6 k at `capinstart`, 0.460 vs 0.630 at the walls).
* Adversarial panels (thread record): `multi` robust (sparse −2 %, junction −32 % on the blind row,
  stranded ±0.2 %, zero controls untouched); **the inside bound HARMS on sparse probes (+8 %)** — the
  uncertified enrichment sign at an exon|exon terminus, as derived — so it is refuted for shipping and
  the `src/` candidate is the composed join-only transport ALONE.
* ⭐ THE LADDER RUN IS DONE (thread record, RUNG 4 — THE LADDER): `multi` is NEUTRAL (5/8 better in
  each half, worst 1.006× unstranded / 1.008× stranded, −6 % at the in-scope `g00 ss.50 OFF` zero
  control); `bound` −5 % on the two big DEFERRED rows but +1.3 % in scope and +8 % on sparse probes
  (REFUTED); `flip` catastrophic (falsification fires). Holes ① and ② are CLOSED STRUCTURALLY by the
  composed transport at ~zero cost and ~zero gain — the refit prior already served the walled exons.
  ⛔ OWNER DECISION OWED: ship `multi` as the policy's completion, or record it and call the holes
  closed by the prior's reach. Either way the remaining audit items are strand-change faces and the
  AMBIG tilt channel (the both-stranded step, last).
* (superseded) THE DECISION WAS THE LADDER RUN, launched at the end of the session:
  `rung4_proto.py ladder transfer,multi,bound,flip all` → `scratchpad/rung4_ladder.log` (~2 h). Read
  the two halves apart. If `multi` wins the unstranded half without stranded harm and the bound stays
  inert, the candidate for `src/` is the composed transport ALONE (one mechanism), the inside bound
  recorded as derived-but-valueless-on-measurement; then the adversarial probe panels, then promotion
  with fail-first gates (identity to rung 3 with MULTI off; the join/leave predicate on a − strand
  structure; the flip).
* ⚠ Two benign-panel rows (`g25 ss.70/.99 OFF`) certify at COMPOSITION level only (gdna-field-uniformity
  flags 1–2 of 981 slots; shared OFF simulation, so one event on every probe panel). Not chased.
* `panel.py cache` still lacks the g00 pre-warm / `_main` copy / certify steps (recipe in `TESTING.md`
  §0a; the driver died with the session) — folding them in is owed.

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
