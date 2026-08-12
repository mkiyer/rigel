# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that separates **RNA from genomic-DNA contamination**. A
single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the library into gDNA vs
RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package `rigel-rnaseq`; the import and CLI
are `rigel`.

## ⛔⛔⛔ AXIOM 0 — RNA IS RNA. READ THIS BEFORE ANYTHING ELSE

**There are THREE populations and there is no fourth: `gDNA`, `RNA+`, `RNA−`.**

⛔ **"Mature" and "nascent" are NOT populations, NOT species, and NOT a degree of freedom.** RNA inside an
intron is RNA that has not spliced *at that position*. The only distinction that exists is whether a
fragment **is spliced** — certified RNA, gDNA cannot splice, needs no deconvolution — or **is not**, which
is the entire deconvolution problem.

⚠ **This paragraph replaced one that said Rigel "jointly models mRNA, nascent RNA and genomic DNA", and
that sentence caused the error it is written to prevent** (2026-08-06): a derivation opened with the
population set `{gDNA, nascent+, nascent−, mature+, mature−}` and every table built on it came out wrong.
The axiom was already in `DESIGN.md` §0 and in memory; the *first sentence of this file* contradicted it,
and the first sentence won.

⭐ **The set is a function of TWO BITS, which is what makes this structural rather than something to
remember:**

```
T(slot) = {gDNA}                                  # always — gDNA is genomically continuous
        ∪ {RNA+ if statics.free_pos[slot]}        # iff the annotation admits that strand here
        ∪ {RNA− if statics.free_neg[slot]}
```

so `|T| ∈ {1,2,3}`, always. **No expression in this codebase may produce a fourth population.**

⛔ **The tell that you have just violated this:** you have written a population set with more than three
members, or the words "mature"/"nascent"/"subspecies" while describing a *solver*, *composition* or
*population* question, or you have concluded that something needs "a third component". Stop, and re-ask the
question as **"what is this channel's OPPORTUNITY for RNA at this object?"** — the answer is a geometry, it
is derivable from the index, and it dissolves the objection every time.

⚠ The words survive in exactly two places, as **simulator inputs and never as model concepts**: the
simulator's `nrna_abundance` knob and the toy harness's `--nrna` arm, both of which exist to put RNA inside
introns so the solver can be tested on it.

## ⭐ The docs — read them in this order

There are six. They do not overlap, and none of them is a changelog.

| doc | what it is |
|---|---|
| ⭐⭐ **`docs/SUCCESS.md`** | **START HERE.** How performance is measured in this phase: **Stage A** (the accumulator — faithful, unbiased, sufficient?) and **Stage B** (calibration init + pass-0, against an oracle). Says what "done" means for each, and why the end-to-end number is a thermometer rather than a target |
| ⭐⭐ **`docs/ROADMAP.md`** | Where the tool is, in six measured numbers, and **what to do next in priority order** |
| ⭐ **`docs/TRAPS.md`** | **Mistakes this project has already made.** Read before designing anything. Lessons, not measurements — several of these were made twice |
| **`docs/EQUATIONS.md`** | The derivations the code depends on. Deposit rule, opportunity functions, the 2×2, strand, overdispersion, damping |
| **`docs/DESIGN.md`** | What is built and the rulings behind it. Settled — do not re-litigate. ⭐⭐ **§0 is the binding VOCABULARY** — NODE, EDGE, step, the mass pin, structurally pure-gDNA object — and it lists the banned synonyms. Read it before writing a comment or a docstring |
| **`docs/TESTING.md`** | The simulated panels, how to build them, the simulator's own gates, and what the suite can and cannot judge. ⭐⭐ **§0b is the TOY HARNESS** — how to define a test transcript and calibrate it in under a second; read it before hand-rolling any scenario |

Reference rather than design: `docs/MANUAL.md` (user-facing), `docs/PUBLISHING.md`.

⭐⭐ **`docs/dev/` IS THE SANDBOX, AND IT IS ENCOURAGED.** Working notes, a feature write-up, a handoff, a
half-finished argument — put them there. ⛔ **Nothing in `docs/dev/` is authoritative and nothing may cite
it**: not the source, not a test, and not one of the eight permanent docs. It is how the owner talks to
collaborators mid-flight, so it is expected to be provisional, contradictory and occasionally wrong.

⛔⛔ **THE ONE RULE THAT MATTERS: A DEV DOC MUST NEVER BECOME THE STATE.** Two of them once grew to
**1,181 lines — larger than DESIGN + ROADMAP + SUCCESS combined** — and a new session was being pointed at
one as "THE STATE". ⭐ The failure was not that they existed; it was that nothing ever moved out of them, so
the temporary copy became the real one. **When a finding is settled, MOVE it — a current number to
`ROADMAP.md`, a lesson to `TRAPS.md`, a ruling to `DESIGN.md`, a derivation to `EQUATIONS.md` — and delete
it from the dev doc in the same edit.** Copying is what creates two homes; moving does not.

⚠ A dev doc that has not been touched in a while is not a problem. A dev doc that a permanent doc *depends
on* is.

⛔ **THE SOURCE DOES NOT REFERENCE THE DOCS, AND MUST NOT START.** Docs evolve and rot — 73 % of the
citations that used to be in the source pointed at documents that had already been deleted. A docstring may
cite a **test** or the executable specification (`tests/native/_accumulator_reference.py`), because code
cannot rot silently; it may not cite a document.

⛔ **`tests/native/_accumulator_reference.py` is the executable specification** for the accumulator. The
C++ is gated on byte-identity to it; where it and a document disagree, it wins.

## ⭐⭐⭐ WHERE DOES A CHANGE GO? — the calibration layering

`src/rigel/calibration/` is 39 modules. It is **not a knot** — measured from the AST there are no import
cycles, and 18 of them have exactly one importer. It was a **FLAT PILE of peers**, which is the one shape
that cannot tell you where to add anything. `rigel/calibration/_layers.py` names the layers that were
already in the edges, and `tests/calibration/test_layering.py` enforces them.

**THE ONE RULE: an import may point DOWN a layer or SIDEWAYS within one. Never UP.** A module that needs
something from a higher layer is telling you the thing belongs lower — a TYPE almost always does.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally MEANS | **1 · the payload view** — `splice_graph` `substrate` `region_arrays` |
| how many places a fragment COULD have sat | **2 · opportunity** — `effective_length` `capture_eff_length` `junction_opportunity` `gdna_opportunity` `fl` |
| one slot's own numbers, and ψ | **3 · geometry + the per-slot solve** — `node_geometry` `simplex_logodds` |
| which strand a fragment came from | **4 · strand** — `gdna_strand` `strand_deconv` `strand_balance` `strand_summary`, and `strand_likelihood` (an executable REFERENCE, gated) |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `density_deconv` `npmle` `gdna_landscape` `background_reference` `run_fill` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `messages/` (the policy) + `node_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `diagnostics` `track` |

⭐ **Run `python scripts/design/module_census.py` rather than trusting this table** — it re-derives the
graph, prints every layer, and flags any import that points up. ⚠ It also flags **docstrings that name a
sibling with no import edge**: 14 were measured on 2026-08-07 and **6 were genuinely stale**, including one
claiming four consumers that had one. A layer is *not* a claim that its modules are the right SIZE — layer 4
is five modules for one concept — and that question is deliberately still open.

## ⛔⛔⛔ MESSAGE PROPAGATION IS OFF BY DEFAULT

`CalibrationConfig.message_propagation = False` (owner, 2026-08-07) installs `messages.silent.SilentPolicy`
— ψ carries each slot's OWN evidence alone. ⭐ A measurement put it there: muting the relay is a net
improvement on **three of the four strata** (−58.3 % / −43.7 % / −32.1 %, and 16/16 on the two stranded
ones). ⛔ **The price is large and concentrated**: +154.8 % on unstranded × capture-ON, which carries 73 %
of the panel's error, and on zero-gDNA golden scenarios the false-positive gDNA mass goes up ~1,900–3,100×.

⚠ **So this is a STUDY configuration, not a shipping one**, and it stays down until the tool is optimised
end to end across all scenarios (owner, 2026-08-10). Turning it back on is one config flag, and every
operator inside `HeadPolicy` is still behind its own named switch.

⛔⛔ **AND THE PLANNED EXIT HAS BEEN MEASURED AND IS NOT ONE — DO NOT REACH FOR IT.** This paragraph used
to say the exit was to give an AMBIG slot its own composition evidence via `length_likelihood`, "the only
θ-independent channel". That channel was A/B'd for the first time on 2026-08-10: on the flgap pair it is
better on 7 of 8 conditions and appears to resolve the blindness (0.0324 → 0.5222 against a truth of
0.507) — and on the `g00` ZERO-gDNA control it reports **54–57 % gDNA in a library containing none**. It
returns a near-CONSTANT ~0.5, and the flgap panel is all `g50`
(`TRAPS: a-single-level-panel-cannot-see-a-constant`). ⭐ `ROADMAP.md` §1.4 is the table; `EQUATIONS.md`
§3d–§3e are the derivations and what would have to change first.

## ⛔⛔ CITE A RULE BY ITS NAME. NUMBERED LABELS ARE BANNED

`TRAPS.md`'s rules used to be `A16`, `D4j`, `C0b`. They are **named**: cite them as
`TRAPS: off-grid-message-mode`, `TRAPS: a-cancelling-defect-pair`,
`TRAPS: frame-free-is-not-assumption-free`. ⭐ The name IS the identifier, so a citation says what it means
without a lookup and is still one greppable string with one home.
⛔ `tests/test_no_jargon_labels.py` enforces it and its allowlist is scoped, never blanket.

⚠ **The numbers were not merely opaque, they were AMBIGUOUS.** `A1` was a validation trap here *and*
`SUCCESS.md`'s FIDELITY criterion; `C1` was a pool trap *and* a moment variable in `length_likelihood.py`;
and **`G1` meant "no magic numbers" AND "a structurally pure-gDNA object" — 201 occurrences, and a reader
could not tell which.** That is this project's own `TRAPS: two-masks-one-name`, committed by the labelling
scheme. The rename rewrote **980 citations across 114 files** and correctly left 66 alone.

## Working rules

- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- ⭐⭐ **Measure the CEILING before building the correction** — *when a ceiling is the right instrument.*
  `calibration_truth_ab.py --ceiling`. It has re-ranked the project twice and is how you learn a phase is
  finished. ⛔ **But a ceiling prices something that may be unreachable, so it is NOT the default** (owner,
  2026-08-05). The loop that produced the 39 % win is: run the panel → take the worst scenario → dissect it
  to the single highest-error object → find the cause → fix → REPEAT. Reach for that first.
- **One thing varied per experiment**, and a baseline re-recorded from the current tree in the same
  session.
- **Score against TRUTH, not against the previous run.** The simulator writes per-fragment ground truth
  into the oracle BAM's read names.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. Code kept "for
  comparison with the old version" is a defect. ⛔ No version suffixes in file names — it is
  `accumulator.py`, never `accumulator_v5.py`.
- **No Greek letters in identifiers** (fine in maths write-ups).
- ⛔ **Real data is a TEST input, NEVER a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
  spectrum, not a sample of it. Sweep the plausible space, report the worst case, and bring the owner the
  domain call.
- **Profile on real cfRNA, never a small synthetic suite** — a toy ranks hotspots backwards.
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

> **Every build, test and lint command must run inside the activated `rigel` conda environment** — it holds
> htslib and the compilers, and the C++ build finds htslib via `$CONDA_PREFIX`.

⛔ **Run the suite as `python -m pytest`, never bare `pytest`.** Bare `pytest` does not put the repo root on
`sys.path`, so an intra-repo import raises and the suite reads one extra failure.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q                     # baseline: 0 fail / 3298 pass / 0 skip / 10 xfail
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

⭐ **THE STANDING BASELINE IS GREEN: 0 failures, 3,298 passing, 0 SKIPPED, 10 xfail** (measured
2026-08-12). **Any failure at all is a regression** — a stronger and
cheaper rule than counting expected ones.

⚠ **It was 3,235 / 7 xfail earlier the same day; the delta is ACCOUNTED, not adjusted.** `+62` is the
new per-script import gate (`test_scripts_index.py`, one case per file across `scripts/design/` and
`scripts/sim/`) and `+1 xfail` is `prior_units_check.py`, which is broken on a deleted API and carries a
DECISION rather than a repair.

⛔⛔ **AND THAT GATE EXISTS BECAUSE A GREEN SUITE HID FIVE DEAD INSTRUMENTS.** Nothing checked that a
script still IMPORTS — only that it was indexed and had a docstring, both true of a file that raises on
line 1. Two commits took five out: `94d283c0` deleted the fixed-point layer (`INV_LENGTH_SCALE`,
`inv_length_quantum`) and killed three, `0d9d422b` deleted `enrichment_frame` and killed one, and one
died on `_component_node_arrays`. ⭐ A `src/` deletion is the mechanism, so run the suite after one.

⛔ **RE-DERIVE THE COUNT, NEVER ADJUST IT** (`TRAPS: re-record-the-baseline`). Several suite tests are
PARAMETRISED OVER DOC, SOURCE AND SCRIPT FILES — roughly 2 per src module, 3 per script, 1 per doc — so
adding or retiring a file moves the total by a few. Account for the delta; a hand-carried number has
drifted every time it has been carried. ⚠ Four successive versions of this paragraph were stale when
read, which is why it no longer records what they said.

⭐⭐ **7 xfails, and they are NOT one kind of thing** — interrogated 2026-08-11, so do not treat the list
as uniform:

* **2 were GREEN and one config flag broke them.** `test_ambig_scenario` and `test_toy_harness`'s intron
  gate both go green the instant `message_propagation = True` (verified: ambig `0.4577 → 0.000247`). They
  are the *recorded price of the switch*, not defects. ⛔ Flipping it costs **+154.8 %** on the stratum
  carrying 73 % of panel error — the owner's ruling, not a test decision.
* **5 were WRITTEN as xfail and were never green** — executable records of two PROVEN defects found by
  the dissect loop: `struct_lock` is `~solvable & NODE` rather than `g1_locked & NODE`, and `Var(f_g)` is
  capped at `f_g(1−f_g)`, which is 0 at the `f_g = 1` default of an evidence-free slot. ⛔ For these
  "fix the test" is a CATEGORY ERROR: the test is right and the code is wrong. Rescoping `struct_lock`
  alone measures **−1.2 %** on its target stratum, **+0.4 %** worse at `g98` and **+3,207 %** on the
  zero-gDNA control — it is half of a cancelling pair (`TRAPS: a-cancelling-defect-pair`) and must be
  priced WITH the seam-composition arm or not at all.

⛔ **Two were closed on 2026-08-11 and the method is the precedent:** the R1-sense gap by BUILDING the
missing capability, and the paralog collapse by asserting the collapse STRUCTURALLY (`min == 0`,
`max == total`) — which still fires the day the tie is broken AND catches a partial break a strict xfail
could not see. ⚠ Neither was closed by widening a bound.

Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Tooling under `scripts/`

`docs/SUCCESS.md` lists these in the order you would run them. ⚠ Everything not listed here was deleted as
unrunnable; if a script is not in this table, treat it as not existing.
⛔ **That claim is currently FALSE for eight files** and the drift is recorded rather than silently
tolerated: `bam_spans.py`, `implicit_splice_census.py`, `inv_L_limits.py`, `length_sieve.py`,
`partition_test.py`, `reference_on_real_data.py`, `sigma_inv_L.py`, `spanning.py` are on disk and unlisted.
They pre-date 2026-08-05, none was run this campaign, and none should be trusted without re-running it —
either promote a row or delete the file, but do not assume the table is complete until that is done.

| | |
|---|---|
| **⭐⭐⭐ the validation campaign — `ROADMAP.md` §1** | |
| `design/prior_vs_oracle.py` | ⭐⭐⭐ **CALIBRATION'S ENDPOINT AGAINST TRUTH** — `LocusPriors` (`gdna_prior_count`, `rna_prior_count`, `gdna_eff_len`) is what the EM actually reads, and until 2026-08-07 nothing compared it to anything. Arms: **P** shipped, **O** the same assembler fed the origin-split truth masses, **S** plus each component's own true share, ⭐ **Fo** the EM's OWN candidate count per origin, **F** the first-base count. `P−O` is calibration's error and `O−Fo` the ASSEMBLER's — different repairs in different files. ⛔⛔ **`F` IS NOT THE EM's TARGET and every `O−F` printed before 2026-08-08 was scored against a quantity nothing consumes** (`TRAPS: score-the-consumers-own-count`); it is retained and priced on table ⑪. ⭐ Reports the count, the composition claim `phi = a_g/(a_g+a_r)` and the SCALE separately, per stratum. ⛔ Undrained on every arm, priced not waved |
| `design/prior_yardstick.py` | ⭐⭐⭐ **IS THE PRIOR ASSEMBLER RIGHT? — the flgap pair, DRAINED, against the EM's own count.** The one place the assembler can be scored the way production runs it: `prior_vs_oracle.py` must stay undrained (a drained oracle is inadmissible on the ladder), and this reads the flgap study cache, which holds both. ⭐⭐ **Verdict: with perfect masses AND perfect per-component shares the assembler reproduces `Fo` to `rel` 2.8e-5 … 2.0e-3, and the pooled share is 82–99 % of the whole residual — the "72 % open residual" was the yardstick.** ⛔ Prints the drained arm's contamination FIRST (flgap_short leaks 1 spliced gDNA record, flgap_long 8,641 — short is the verdict, long the cross-check). ⭐ Table ④ scores the composition claim against the UNSPLICED pool — the population `a_g:a_r` describes, because a spliced unit never gets a gDNA candidate — and it is exact there (`Δphi` ≤ 5e-4). ⛔ Scoring it against ALL RNA units reads a phantom +0.07…+0.10 tilt; that mistake was made and retracted on 2026-08-08. ⭐ It also reports the prior's STRENGTH, which is 1.000 pseudo-fragment per real one by construction — the posterior is a 50/50 blend of calibration and the EM's evidence, with no knob |
| `design/quant_accuracy.py` | ⭐⭐⭐ **THE TOOL END TO END, AND THE PRIOR-INJECTION CEILING ABOVE IT** — `--arm base` is the transcript-level accuracy table, `--arm oracle` re-quantifies with the perfect prior injected, and the three `oracle_<field>` arms say WHICH of the prior's three numbers carries the value. ⛔ Scores **count against count** against the condition's OWN `truth_abundances.tsv` (the realised observed fragment count); the suite-level table is the *pre-capture molar* abundance and is a different quantity. ⛔ Pins `EMConfig.seed` — the shipped default is `None` (`TRAPS: the-deliverable-is-not-reproducible-by-default`) — and `base_reseed` prints the noise floor beside the effect |
| `design/mass_prior_ab.py` | ⭐⭐⭐ **THE PRIOR AS A CONSERVED FRAGMENT COUNT — plan 5.5 against the oracle arms.** `edge_unspliced_mass` sums to ONE per fragment, so the assembler can stop manufacturing a count from a density. Subsamples a real condition by qname hash (so the whole and all three origin partitions stay consistent), then re-assembles `O` and scores it against `F`. ⭐ **Verdict: the capture-ON over-call goes `O/F` 1.149 → 1.019 on three gDNA levels; capture-OFF `Σ|Δ|` falls 61 %.** ⛔ Prints the SUBSTRATE check first — the subsample must reproduce the defect before its disappearance means anything — plus a `share≡1` ablation and the plan-8-Q4 ceiling. ⚠ It MEASURES the shipped bank and never re-derives it; the rule's own gates are `tests/native/test_conserved_mass.py` |
| `design/flgap_study_cache.py` | ⭐⭐⭐ **SCAN ONCE, INTERROGATE MANY TIMES** — the flgap pair's four conditions cached whole: scored `multi_loci`, calibration, the DRAINED payloads with their lifted origin partitions, the oracle masses and per-component shares, the raw per-region start counts, and ⭐ the per-unit `(true origin, is_spliced)` pair the `Fo` arm joins on, for both scoring stages. Building is ~5 min/condition; loading is ~1 s. ⭐ **`priors.py` is deliberately OUTSIDE the cache key**, so testing an `assemble_priors` change is a one-second loop — ⛔ which is only sound because **nothing `priors.py` produces is stored**: every arm (P, O, S, F) is derived on load. It used to store `p_arm` and `f_gdna` and served them stale beside a fresh O (`TRAPS: a-hash-that-misses-its-artifact`, second form). ⛔ Keyed on the scan manifest + a content hash of every producing source file, the builder itself included; a mismatch REBUILDS rather than warning. ⚠ Stores `gdna_spliced_leak` and `lift_ambiguous` — **1** on flgap_short vs **1,491 / 8,641** on flgap_long (⭐ BOTH spliced banks summed; the "1,010 / 5,827" once recorded here was `edge_spliced_count` alone, because `_validate` raises on the FIRST nonzero bank and its message was quoted as a total), which is why flgap_short is the admissible panel for a drained measurement and flgap_long a cross-check |
| `design/rescan_panels.py` | ⭐⭐⭐ **RE-SCANNING IS THE FIRST IRREVERSIBLE STEP, SO IT CARRIES ITS OWN FALSIFICATION.** After a deposit-rule or schema change every cache is refused and must be rebuilt — and the old ones are then gone. This rebuilds them and GATES each condition on byte-identity: every bank the change was not supposed to touch must come back identical, and the symmetric difference must be the SAME delta on every condition. ⛔ Reads the stale `.npz` BEFORE the write, because `write_scan_cache` overwrites in place and would otherwise destroy its own baseline; a condition that fails is NOT written. ⛔ No expected delta is hard-coded — it is derived from the first condition and required of the rest, so a side effect appearing only under capture or only at `g98` is caught. ⚠ Counts an ALL-ZERO baseline as VACUOUS rather than a pass (`g00`'s gdna, every `nrna_none`'s nrna). `--self-test` perturbs the comparator with no I/O; `--jobs N` shards conditions |
| **⭐⭐⭐ where to develop** | |
| `design/module_census.py` | ⭐⭐⭐ **WHERE DOES A CHANGE GO?** — the calibration package's real shape, re-derived from the AST: the layering with every upward import, the graph with each module's importers inside and outside, ⭐ **docstrings that name a sibling with no import edge** (14 measured, 6 genuinely stale), and dead public surface. ⛔ It REPORTS, it does not judge — an entry point looks dead, and a gated reference implementation looks duplicated |
| **⭐⭐⭐ the backbone** | |
| `design/arm_identity.py` | ⭐⭐⭐ **IS THIS ARM BYTE-IDENTICAL TO THAT ONE? — the A5 gate.** `arm_score.py` AGGREGATES, so a difference that cancels between two fields is invisible there; this compares **every scored field of every row** and exits nonzero on any difference. ⛔ Both of A5's recorded lies are gated: the row-key sets must be EQUAL (an arm with ZERO rows once scored "32/32 IDENTICAL") and both files' mtimes are printed (a stale baseline is what A5 warns about — re-record it, B8). ⭐ Falsified by perturbation: it resolves a **1-ULP** nudge and refuses an empty file |
| `design/backbone_parity.py` | ⭐⭐⭐ **WHAT DOES ONE MESSAGE OPERATOR DO, PER SLOT?** — two policies, one real 70,176-slot chain, one process, every output array compared element by element, plus every diagnostic `_capture` key and ⭐ the five backbone assertions with their ELIGIBLE sets (A14). Strictly stronger per condition than the panel and strictly weaker across conditions, so run it FIRST. `--arm-a head --arm-b no_<switch>`. ⛔ Its verdict on the restructure: **421,056 output elements and 18,245,830 diagnostic elements, zero differences** |
| **⭐⭐ the toy harness** | |
| `design/toy_panel.py` | ⭐⭐ **one toy spec × EVERY cached condition × an RNA-density ladder, scored per object** — what you run once a structure is the *target* rather than a probe. Prior-free pass-0 by default; the RNA density is a multiple of each donor's OWN gDNA density. ⛔ Reports which object carries the error, whether the messages helped or hurt it, and which objects are CONFIDENTLY wrong. ⚠ 13 s per condition; shard with `--conditions` |
| `design/verify_toy_substrate.py` | ⭐⭐⭐ **IS THE INPUT CORRECT? — no solver runs.** Every accumulator bank re-derived from the simulator's per-fragment TRUTH by an independent implementation, plus the splice combinatorics and the length marginal. ⛔ Run this on any new toy spec BEFORE reading a solver number off it |
| `design/verify_capture.py` | ⭐⭐ **hybrid capture, probes ON vs OFF on identical geometry** — the gDNA landscape, the length selection and the junction depletion, each gated on the direction the knobs predict |
| `design/toy_trace_error.py` | ⭐⭐⭐ **the full error trace for ONE scenario** — every node and edge ranked by error MASS, the four init sources, the relay hop by hop with the licence state, then a ψ channel ABLATION (fidelity-checked to 1e-6) and ⭐ a RELAY-level arm that withholds the composition licence in both twins at once. ⛔ This is what "dissect the error" means |
| `design/zero_controls.py` | ⭐⭐⭐ **THE TWO ZERO CONTROLS — and the owner requires them on every experiment.** ZERO RNA (silent genes, truth `f_g = 1.000` everywhere) and ZERO gDNA (the `g00` donor, truth `0.000`). The truth is a CONSTANT, so every deviation is a false positive with nothing to cancel it. ⭐ The zero-RNA arm is the **biologically dominant** case: most annotated transcripts are OFF. ⛔ Prints the three-rung ladder (strand → self-solve → final), what the MESSAGES delivered as an implied `f_g`, and flags any EMPTY object, because a zero arm can be degenerate and then it tests nothing (`TRAPS.md` A14) |
| `design/reframe_walk.py` | ⭐⭐⭐ **THE REFRAME, WALKED — every count, both flank totals, and every hop in BOTH directions, on one two-exon transcript.** Not a summary statistic: per hop it prints `r` as used, `r` as the predecessor's single total would have given it, what TRUTH says the same ratio is, and the true gDNA-density ratio beside it — so the include/exclude decision at each `(EDGE, side)` is a number. ⭐ Three gDNA arms (zero / high / very high), capture-OFF × unstranded, HIGH RNA so nothing is sparse. ⛔ The geometry is rebuilt and GATED against the solver's published frames to 1e-12 |
| `design/toy_dissect.py` | ⭐⭐ **one scenario opened to every slot and every channel** — the solver's own ladder (strand-only → self-solve → final) beside `cm_g` / `c_tau` / `cg`, so a dead channel is visible as a zero rather than inferred. ⛔ Step 3 of the debug loop, for a toy |
| `design/certified_rna_audit.py` | ⭐⭐⭐ **IS THE CERTIFIED-RNA CHANNEL WIRED?** — `edge_spliced` is a spliced fragment, so it cannot be gDNA and needs no deconvolution. Audits (a) is the bank populated, (b) does it have a divisor, (c) **is a precision emitted** — any one failing looks identical. Runs a TA × TB abundance grid on `tes_readthrough`. ⛔ Scores the **mass**, never `f_g`: the solver's `f_g` is the UNSPLICED fraction and the oracle's is spliced-inclusive, and confusing them reads a 1.5 % defect as a half-unit error. Its docstring carries the 24/24 verdict |
| `design/certified_q_census.py` | ⭐⭐ **CAN A CERTIFIED COUNT SPEAK ABOUT THE UNSPLICED SPLIT? — the answer is NO, and this is why.** Measures the splice-visibility `q = S/(S+C_R)` directly off the origin-split oracle on all 36 conditions; no solver runs, so the whole ladder is seconds. ⛔ Its verdict is a NEGATIVE and the docstring says so first: `q`'s mass-weighted median is 0.19–0.71, so the term §2d drops is the same size as the one it keeps, and the raw-count λ term is **worse than the uninformative reference on 12 of 36 conditions** (worst +0.4578). ⭐ Its two extreme rungs ARE the two zero controls, which is what makes the diagnosis certain |
| `design/vertex_ceiling.py` | ⭐⭐ **the re-solve ceiling on the REAL ladder** — pins oracle truth at a chosen object class in `build_node_init` and re-solves, with a `noop` arm that MUST be byte-identical (`TRAPS.md` A5). ⛔ Its own verdict is a NEGATIVE and the docstring says so first: the 24.4 % it measures is the value of MISSING INFORMATION, not headroom. ⭐ Reusable override plumbing for any `build_node_init` prototype |
| `design/toy_ceiling.py` | ⭐⭐ **the RE-SOLVE ceiling** — hand one object class a different own belief and re-solve the whole chain, six arms sharing one simulation. ⛔ This is what `TRAPS.md` B17 demands instead of a substitution, and `--arms base noop` is its own falsification (must be byte-identical). Its docstring carries what it measured |
| `design/psi_channel_ablation.py` | ⭐⭐⭐ **WHICH ψ CHANNEL IS DOING THE WORK?** — one level below `worst_objects.py`: that says which OBJECTS carry the error, this says which CHANNEL put it there. Records every argument of the final ψ combine, then re-solves with one `*_imp` channel nulled at a time, so an attribution is a re-solve of the REAL call. ⛔ A5 — the `as-is` arm must reproduce the run bit-identically, and that means reproducing the WRITE-BACK too (the first version read `max\|Δ\| = 1.0` for exactly that reason). ⭐ Read the per-slot table, not the totals: D4h says all-small-singly + large-jointly means the channels share an upstream quantity |
| `design/arm_score.py` · `design/arm_sweep.py` | ⭐⭐ score two `ladder_arm_ab` arms against each other, and sweep a family, **per stratum**. ⛔ Never pooled — the panel total hides a sign flip between strata, and both print `abs_err_all_final` (the deliverable) beside `abs_err_all` (pass-0) on every row, per `TRAPS.md` A15. `$RIGEL_ARMS` points at the `--out` directory |
| `design/ladder_arm_ab.py` | ⭐⭐ **the same override on the REAL 36-condition ladder**, scored by `solvability_audit.py`. ⛔⛔ **Run this before writing a mechanism into `src/`** — FOUR times now a toy- or single-condition-positive change has been panel-negative (`TRAPS.md` B18). ⭐⭐ **`--jobs 8` runs the whole 36-condition panel in ~2.2 min** (was ~20 min): the two unused `c_input_*` arms are no longer built and the MAIN scan payload is cached beside the oracle cache. ⭐ Carries the eight `zc_*` arms plus ⭐⭐ `msgfree_p0` / `msgfree_all` / `msgscale_<k>` / `onesided_rna`, which between them price the whole message layer. ⛔⛔ **`--messages {off,on}` IS PART OF THE ARM AND IS STAMPED INTO EVERY ROW** — under the shipped `off`, **22 of the 26 arms cannot move a number** and 16 of those score byte-identical to `base` while their fire counter reads healthy (`TRAPS: an-ablation-that-never-ran`, third form). An arm the policy cannot express is now REFUSED up front. ⭐⭐ **`--self-test` is the harness's own falsification**: every arm under both policies, gated INERT-under-one / MOVED-under-the-other — 25/25 in ~6 min. `--arm zc_noop` must still be BYTE-IDENTICAL to `base`. ⚠ Every arm file predating 2026-08-11 lacks the `messages` field, so `arm_identity.py` correctly refuses to compare one against a fresh run |
| `design/length_ceiling.py` | ⭐ **what a PERFECT length model is worth, on the ladder, ONE pmf at a time.** ⛔ Pricing the pair together hid a 14× split between them (`TRAPS.md` B21). Reports the pass-0 solvable yardstick beside the mass-weighted one, because they disagree |
| `design/toy_harness.py` | ⭐⭐ **a mini chromosome you define, calibrated in 0.1–5 s** (`docs/TESTING.md` §0b), with every object's answer beside per-object truth. The library-level priors a toy cannot fit are harvested from a real cached condition and INJECTED; the gDNA depth is DERIVED to match that donor, never set by hand. `--list` for the spec ladder. ⭐ Reach for this FIRST when a mechanism needs isolating — it found C1's mechanism candidate in one sweep |
| **the substrate** | |
| `design/simulator_gates.py` | ⛔ the simulator's own gates G-S1…G-S6, scored on per-fragment truth. Run before trusting the panel |
| `design/suite_resolves.py` | ⛔ can the suite resolve the axis you are changing? Run before quoting any suite number |
| `design/index_census.py` | re-derive an index's census — never quote a stored table, run this |
| `design/verify_index_rebuild.py` | nodes byte-identical, edges only in contiguous reach |
| **Stage A — the accumulator** | |
| `design/fl_anchor_gap.py` | the anchor and both length models vs truth. `--drain` measures before and after the second pass |
| `design/gdna_pool_census.py` | ⭐ the four gDNA pools, each against its own opportunity and against truth |
| `design/second_pass_accuracy.py` | the second pass scored PER FRAGMENT against the oracle BAM's read names |
| `design/observable_efficiency.py` | what fraction of the length information a storage choice keeps |
| `design/node_density_derivation.py` | the reciprocal-opportunity theorem, T0–T6, each perturbed |
| **the scoping number** | |
| `design/calibration_truth_ab.py` | ⭐⭐ the deliverable against truth, and `--ceiling` — what perfecting each length model is *worth* |
| **Stage B — calibration** | |
| `sim/configs/gdna_ladder.yaml` | ⭐ **the Stage-B panel** — 36 conditions, EQUAL fragment lengths, gDNA 0 → 98 % at a fixed 10 M total. An ORACLE CACHE ships beside it at `ladder/oracle_cache/` and stays valid across every calibration change. `docs/TESTING.md` §0 explains both panels |
| `design/solvability_audit.py` | ⭐⭐ **START HERE FOR PASS-0.** Which objects are SOLVABLE, which of those are solved wrong, and — the one that matters — which are **confidently** wrong, since a tight-variance wrong value outvotes correct neighbours and anchors the prior. ⛔ Excludes the undetermined population: in pass-0 an object with no own evidence reporting `f_g ≈ ½` at zero precision is *correct*. Carries the ablation ladder (strand → local → final) and a calibration curve that needs no threshold. `--suite` alone gives the whole panel |
| `design/pass0_vs_oracle.py` | T (the origin-split payload) vs P (`calib_refit_iters=0`) vs two levered ceilings, per object and per class. ⛔⛔ **ITS HEADLINE MASS-WEIGHTED ERROR IS THE WRONG YARDSTICK FOR PASS-0** — it scores every object with mass, so honest ignorance reads as error and buried a 0.0456 answer inside a 0.3150 one. Use it for the T/C/P decomposition and the oracle plumbing; use `solvability_audit.py` to judge pass-0. ⛔ undrained on every arm. `--oracle-cache` makes repeat runs cheap: the oracle depends on the accumulator and index, never on calibration, so one cache serves a whole debugging campaign |
| `design/worst_objects.py` | ⭐⭐ **step 3 of the debug loop** — one condition dissected to individual nodes/edges, ranked by error **MASS**. Read the CONCENTRATION curve first (concentrated ⇒ a mechanism exists; diffuse ⇒ a systematic bias). `fg_loc` vs `pred_fg` separates a wrong local solve from wrong messages |
| **diagnostics** | |
| `design/anchor_opportunity_census.py` | ⭐⭐ **IS A ZERO-COUNT ANCHOR'S DENSITY CLAIM TRUE OF ITS NEIGHBOURHOOD? — no solver runs.** Every slot's TRUE gDNA density re-derived from the origin-split oracle through the SHIPPED `build_node_geometry`, then the empty structurally-pure-gDNA population against what its own chain neighbours measure. ⛔ Its verdict is that the claim is false by **346×** under capture and true off it — which is NECESSARY for the "empty means no probe here" hypothesis and, as the panel arms then showed, not sufficient. ⛔⛔ **THAT VERDICT IS ABOUT 1,312 ANCHORS, NOT ABOUT `struct_lock`** — the docstring AND the mask both claimed "exactly `strand_evidence`'s `struct_lock`" until 2026-08-11 and the solver's mask is **15–23× larger** (30,423 vs 1,312 at `g00`), because it also holds every zero-count NODE. This instrument's population is `g1_locked ∧ NODE` — what the standing xfail wants `struct_lock` rescoped TO. ⭐ Both sizes now print on every row, and the measured-redundant `∧ intergenic` term is ASSERTED rather than deleted |
| `design/composition_evidence_census.py` | how much library mass reaches the solver with NO composition evidence. `--inject-kappa 0.5` is its falsification handle |
| `design/held_flux_census.py` | how often a held candidate has ZERO flux evidence, by cause |
| `design/prior_units_check.py` | the EM prior in fragment units vs the old incidence sum |
| `design/edge_q_population.py` | ⭐⭐⭐ **IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND? — DRAINED, no solver runs.** `assemble_priors` divides a line's crossing INCIDENCE by `q = mass/count` to get a fragment count — but `q` is measured on the POOLED population and applied to gDNA and RNA separately. Feeds a PERFECT `f_g` from the origin-split oracle so the `q` conversion is isolated from calibration's own error, and scores `Σ count_c·q_pooled` against the conserved truth `Σ mass_c`, EDGE-only and on the TOTAL prior the EM actually consumes. ⛔ **Its verdict: the defect is real, is NOT driven by fragment length (the equal-length null shows `q_g` 0.633 vs `q_r` 0.523, a LARGER error than a 102 bp gap), and is bounded at ≤ 0.6 pp of composition — so record the bound, build nothing.** ⭐ The driver is PLACEMENT: gDNA crosses lines in long intergenic nodes where `q → 1`, RNA in short exon nodes. ⭐ Gated on the pooled `q` reproducing the shipped `mass_per_crossing` bit-for-bit and on the partition summing to the whole; ⚠ its Ⓖ4 records that ladder capture-ON is NOT an equal-length null, because capture manufactures a +15.33 bp gap |
| **plumbing** | |
| `design/build_scan_cache.py` | scan once, calibrate many times. ⚠ re-run after any accumulator change — the payload schema digest invalidates every cache, by design |
| `design/native_parity_on_real_data.py` | the native-parity gate on real cfRNA at full scale |
| `design/scan_profile.py` | ns/fragment, regressed over several BAMs |
| **⭐⭐⭐ the workflow** | |
| `sim/panel.py` | ⭐⭐⭐ **THE SIMULATION + BENCHMARKING WORKFLOW, ONE COMMAND PER STAGE** — `status` / `build` / `simulate` / `cache` / `score` / `report`, every path derived from one panel YAML. ⛔ It adds NO measurement code: each stage shells out to the instrument that already owns it, because duplicating a scorer is how a baseline and a ceiling drift apart. ⭐⭐ **`status` FIRST** — every stage is expensive and resumable, and it names the next one. ⛔⛔ **`cache` builds BOTH caches**, and the oracle one is why this exists: it is the origin-split truth every scoring instrument reads, and until 2026-08-11 it was a SIDE EFFECT of `pass0_vs_oracle.py` with no step of its own, so the documented recipe produced a panel every scorer refused. ⚠ A cached oracle condition needs all FOUR parts (`gdna`/`mrna`/`nrna`/`_main`). Gated by `tests/test_panel_workflow.py`, falsified by perturbation |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | build the panel. ⚠ `panel.py build` drives the last two; the reference carve needs the SOURCE genome/GTF, which a panel config does not name, so it stays manual |
| `rigel.sim.net_flow` (a MODULE, not a script) | ⭐⭐ **WHERE DID EACH MISASSIGNED FRAGMENT GO?** — the one question `quant_accuracy.py` does not answer. It measures the MAGNITUDE of transcript error; this decomposes its DIRECTION, per transcript, into gDNA-sourced and RNA-isoform-sourced flow. ⛔ **All that survived `sim/analysis.py`, retired 2026-08-11** (owner): that was a 1,589-line SECOND SCORER rendering its own accuracy tables against its own truth, and two scorers is how a baseline and a ceiling drift apart. 618 lines kept with their tests, 971 deleted along with `scripts/sim/evaluate_suite.py`. Gate: `tests/test_net_flow.py` |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
