# NEXT SESSION — THE MESSAGE LAYER IS MERGED TO `main` (2026-09-09); START FROM `ROADMAP.md`'s RANKING (handoff)

⭐⭐⭐ **THE REFERENCES:** `CLAUDE.md` (the scope, the message layer's shipped design in one place, the
instrument table, the standing baseline), `ROADMAP.md` (the ranking, audited on the merge day),
`ISSUES.md` (every open case with a priority; every refusal with its number), `DESIGN.md` §6b.4–§6b.14
(every ruling of the shipped policy with its measurement). This file is the state.

## WHAT STANDS (2026-09-09; `main`)

1. **`message-layer` is merged to `main` (a fast-forward).** The `transfer` policy is the shipped default;
   the relay and its baggage are gone (git carries them); the cleanup on the one shipped policy is done;
   the empty-piece flux source is landed; `vertex_ceiling.py` measures again. The owner's ruling: this
   closes the development of the new message policy — future work proceeds from `main`.
2. **THE SUITE: 0 failed / 3,593 passed / 2 xfail, 3,595 collected** (`CLAUDE.md`'s baseline line carries
   the accounting). The 2 xfails: the toy harness's intron-independence gate (`ISSUES: two-sided-exon-row`)
   and the antisense t2 prior-assembly casualty (rank 2's territory). `preflight.py --full`: everything
   present, 16/16 instrument self-tests. `module_census.py`: no upward import, no dead public surface.
3. **THE IDENTITY REFERENCES for the next cleanup** are the post-landing pair
   `~/Downloads/rigel_runs/arms/fluxsrc_identity_{gdna_g05_ss_0.99_nrna_mid_capture_on,gdna_g05_ss_0.50_nrna_mid_capture_off}.json`
   (`rename_identity.py --check --condition … --reference …`). The `retire_identity_*` pair is the
   pre-landing tree and is history.
4. **THE MEASURED STATE OF THE MESSAGE LAYER on the merge day (re-derive with the instrument named,
   never quote):** `calibration_vs_oracle.py` (mass-weighted |Δ gDNA share| per object, region / boundary
   axis) — stranded × OFF 0.0067 / 0.0104, stranded × ON 0.0104 / 0.0149, unstranded × OFF 0.0088 / 0.0124,
   the `g00` zero controls 0.0081 / 0.0077, the DEFERRED unstranded × ON 0.0742 / 0.1166.
   `policy_benchmark.py --panel ladder --by-class` — transfer beats silence 7/8 on each half (worst
   1.01× / 1.00×); the misplaced gDNA as a share of every non-intergenic fragment is 0.6–2.4 % on every
   in-scope contaminated row; by class, the intron 1–2.3 % of its own fragments off capture (the largest
   class by mass, 43–45 % of the absolute error), exon|intron boundaries 4–8 % (about one fragment each,
   the counting floor), every exon and exon|exon class ≤ 2.5 % except the `g98` capture-OFF rows at 4–9 %
   (near-pure gDNA under-called: the vertex atom). ⭐ The owner's reading (2026-09-09): messages are at
   diminishing returns in scope (the bar was 3 %); move on to the landscape prior.
5. **THE VERTEX CEILING under the shipped policy** (`vertex_ceiling.py`'s docstring carries the
   whole-ladder table; `ROADMAP.md` rank 4): at most 1 % on every stranded in-scope row, 2–7 % of the
   unstranded capture-OFF rows rising with gDNA, 25–35 % of the deferred stratum; the zero rows are total
   by construction. The in-scope value sits on silent genes' regions and nascent-free introns — the
   landscape prior's training population — and arrives through the refit prior, not at pass zero.

## WHAT IS NEXT — `ROADMAP.md`'s ranking, each with its first concrete step

1. **THE gDNA LANDSCAPE PRIOR** (`ISSUES: gdna-landscape-trains-on-false-positives`; the owner's ruling of
   2026-09-06 in its entry: the enrichment witness IS the landscape prior, and nodes whose only evidence is
   a bound do not train it). First step: census the training population per node class on the zero rows
   (which slots train the prior, with what evidence), then the estimator at bound-only nodes —
   `abundance_landscape_census.py`, `landscape_head_to_head.py`, `calibration_walk.py` (the refits rung
   alone). `docs/dev/PLAN_measured_prior.md` is the thread's sandbox record.
2. **THE POST-CALIBRATION, PRE-EM SETUP** (`priors.py` / `result.py` / `derive.py` against
   `prior_vs_oracle.py` and `calibration_vs_oracle.py`'s ruler column). First step: re-run
   `prior_vs_oracle.py` so the assembler's own error is re-recorded under the shipped policy before
   anything moves.
3. **THE MESSAGE POLICY, only where a row is above the bar** — one prototype arm at a time through
   `policy_prototype.py --module`, three panels and the ladder, halves apart, pass zero beside the
   pipeline: `ISSUES: flux-price-witness-units` first (the open defect the landed flux level shares: the
   column count as the witness at a gDNA-rich node), then `two-sided-exon-row`, `flux-floor-dispersion`,
   `ambig-node-as-a-gdna-source`, `message-layer-open-cases`. ⛔ Prototypes subclass the shipped
   `TransferPolicy` (`_LevelLane`, `prepared.lanes[...]`) and are compared `src` against `src` across a
   landing; the ship-day prototypes under `~/Downloads/rigel_runs/prototypes/2026-09-09_ship/` still
   name the pre-cleanup `_Lane`/`_RnaLane` and need a one-line rename before they import.
4. **THE VERTEX ATOM** — priced (above); a mechanism for it is the landscape's, not a message's.

Then `ISSUES: performance-memory-bounded-solve` (owner: mandatory before 0.8.0).

## FOUND, NOT CHANGED — proposals for rank 3 (each would move a number, so each is a measurement)

* **An edge boundary's own strand profile is overwritten by the edge level's marker.** `_claims` gives a
  single-strand boundary with a live channel its strand profile; `_edge_level` then sets `own[b] = zeros`
  at every intergenic|exon edge as the level claim's marker, so the RNA lanes never read an edge
  boundary's strand profile as an RNA source (the gDNA lane is unaffected: a gene edge's level is its
  Poisson count). Preserved exactly.
* **Two liveness predicates.** The policy's test of a node's strand channel is `tau_lam > 0` (`_Chain`);
  the instruments' one home is `has_own_composition_evidence` = `tau_lam > 1e-9`. They agree on every
  value observed and the one-home gate asserts it on a vector; a ladder census of `0 < tau_lam ≤ 1e-9`
  would close it for good.
* `PsiMessage.lam_rows`' docstring still narrates the retired relay's certified-flux stream as the
  channel's origin; harmless history, one paragraph.
* A few rung/item labels remain in `DESIGN.md` §6b.10–§6b.13's bodies and in `ISSUES.md`'s older
  entries; the owner's decision covered §6b.4–§6b.9.

## THE LESSONS THE LAST SESSION PAID FOR

* **A census before a prototype.** The empty-piece flux source's own example (the sj+terminus piece)
  cannot be reached by the lane at all — a strand's RNA level stops at that strand's own terminus by the
  face rule — and the test chromosome had no substrate for it; the census said where it could act (the
  AMBIG exon|exon boundaries of the overlapping loci, ~70 pieces per strand on the ladder) and which panel
  could judge it (the ladder only) before a line of mechanism was written.
* **An instrument's population is part of its measurement.** The vertex ceiling classified by the
  REALIZED vertex and priced chance: tiny boundaries whose few crossings were all gDNA by luck, pinned as
  certain and propagated, made the region axis 23 % worse. The PARAMETER vertex — a vertex by
  construction of the scenario — is the population a ceiling means.
* **A pin must enter where the solve reads.** Under the two-phase design a node's neighbours receive its
  OWN CLAIM and the node answers through ψ; rewriting a belief the backbone reads only into its
  diagnostics capture is inert, and the instrument's own comparator would have called it "did not fire".
* **A price at a source with no witness follows the hop's own rule.** The zero count's counting
  (`hop_price` at n = 0) kept the empty-piece source inside the bar on every in-scope row; the sharper
  price (the junction's counting alone) failed it by 1.8 % on a stranded capture-ON row. Nothing new was
  invented: the existing rule already said what a zero witness pays.
* **Two references, one per half, make a refactor provable one stage at a time.** Five cleanup stages,
  ten identity checks, zero bits moved; the one fold that could have moved a bit (`count_price` into
  `hop_price`) was proved safe by the gDNA lane's counts being positive, not by the suite.
* **Delete a thread record only after grepping its numbers in the permanent docs.** Three of the eight
  sandbox records deleted had refusals whose numbers lived nowhere else; the MOVE RULE turned those into
  ISSUES entries instead of losses.
