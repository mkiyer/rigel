# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that separates **RNA from genomic-DNA contamination**. A
single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the library into gDNA vs
RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package `rigel-rnaseq`; the import and CLI
are `rigel`.

⭐⭐⭐ **The release target is 0.8.0 and CALIBRATION is the focus.** The SCOPE section below is the frame
for every change in this phase — read it before choosing what to work on.

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

## ⛔⛔⛔ THE 0.8.0 SCOPE — READ THIS BEFORE CHOOSING WHAT TO WORK ON

Owner ruling, 2026-08-14. The version on disk is **0.7.1** (`pyproject.toml`); the target is **0.8.0**.
This is the frame every doc, every experiment and every ranked list hangs off.

⭐⭐⭐ **THE FOCUS IS CALIBRATION, AND THE METRIC IS THE CALIBRATION RESULT SCORED AGAINST ORACLE
CALIBRATION — not the end-to-end transcript number.** The transcript number stays a thermometer
(`SUCCESS.md` says why); a change is judged by `solvability_audit.py` and `prior_vs_oracle.py`.

**THE FOUR STRATA — THREE ARE THE TARGET AND ONE IS NOT:**

| stratum | 0.8.0 |
|---|---|
| unstranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-ON | ⭐ **IN SCOPE** |
| unstranded × capture-ON | ⛔ **DEFERRED** |

⛔⛔ **DEFERRED IS NOT DROPPED, AND THAT DISTINCTION IS THE WHOLE RULING.** unstranded × capture-ON
**remains in every benchmark and every measurement, and must keep being reported**. It is simply not a
development target until the other three are fully optimised. If it improves as a side effect of other
work, that is a free win.

⚠ **It is also where the error is, which is exactly why this had to be a RULING rather than a
follow-the-error loop.** On the rebuilt 16-condition ladder it carries **64.5 % of transcript error and
90 % of gene-level error** (measured 2026-08-13/14, noise floor 0.996–1.013); the three in-scope strata
carry the rest. ⛔ On that stratum the tool emits a near-zero gDNA fraction **regardless of truth** — exon
`f_g` **0.040 / 0.0016 / 0.0021** at `g05` / `g50` / `g98` against truths **0.054 / 0.518 / 0.982** — so
it looks acceptable at low gDNA BY COINCIDENCE. ⛔ The debug loop will point here every time; take the
worst **IN-SCOPE** scenario instead.

⛔⛔⛔ **THE LENGTH CHANNEL IS RETIRED UNTIL AFTER 0.8.0. DO NOT PROPOSE IT.** The fragment-length /
"length likelihood" channel **as a CALIBRATION composition channel** is future work: 0.8.0 ships without
it, it is not a candidate, and it may not be listed or ranked. ⚠ It does **not** exist in `src/` — it was
A/B'd once, on 2026-08-10, and the prototype was purged the same day (`b7ed7a0b`), so this is a scope
ruling and not a deletion. The measurement itself is under MESSAGE PROPAGATION below.

⭐ **Three things called "length" are NOT affected and must not be touched in its name:** layer 2's `fl` /
`effective_length` / `capture_eff_length` — the fragment-length **OPPORTUNITY** model, which calibration
and the EM both depend on; `length_likelihood` in `src/rigel/second_pass.py` — per-fragment second-pass
assignment, a different thing that happens to share the word; and the fl PMFs priced by
`length_ceiling.py`.

⭐⭐ **WHY THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS — and the real reason is stronger than the
one that used to be written down.** The old wording was "the length channel is neutralised, so residual
error is attributable to density and strand". ⛔ The actual reason: **the EM ALREADY USES THE
FRAGMENT-LENGTH DISTRIBUTION.** A large gDNA-vs-RNA length difference lets the EM assign fragments on
**LENGTH ALONE**, bypassing calibration entirely and **MASKING bugs in it**. Equal lengths **FORCE** the
calibration phase to be exercised. What calibration then has is strand and density — plus belief
propagation across objects, ON since 2026-08-18 (see MESSAGE PROPAGATION IS ON, below).

⭐ **SCENARIOS ARE CACHED, and that is a scope requirement rather than a convenience**: calibration must
re-run in seconds off a cached scan, so `panel.py cache` / `build_scan_cache.py` come before any
experiment. ⭐ **The ranked list of what to do next is `ROADMAP.md`** — this section says what counts, not
what is next.

## ⭐ The docs — read them in this order

There are six. They do not overlap, and none of them is a changelog.

| doc | what it is |
|---|---|
| ⭐⭐ **`docs/SUCCESS.md`** | **START HERE.** How performance is measured in this phase: **Stage A** (the accumulator — faithful, unbiased, sufficient?) and **Stage B** (calibration init + pass-0, against an oracle). Says what "done" means for each, and why the end-to-end number is a thermometer rather than a target. ⭐ **0.8.0's metric is Stage B's** — calibration against oracle calibration |
| ⭐⭐ **`docs/ROADMAP.md`** | Where the tool is, in measured numbers, and **what to do next in priority order** — the 0.8.0 ranked list. ⛔⛔ **IT IS NOT A CHANGELOG, AND A COMPLETED ITEM IS DELETED RATHER THAN STAMPED** (owner, 2026-08-14) — see below |
| ⭐ **`docs/TRAPS.md`** | **Mistakes this project has already made.** Read before designing anything. Lessons, not measurements — several of these were made twice |
| **`docs/EQUATIONS.md`** | The derivations the code depends on. Deposit rule, opportunity functions, the 2×2, strand, overdispersion, damping |
| **`docs/DESIGN.md`** | What is built and the rulings behind it. Settled — do not re-litigate. ⭐⭐ **§0 is the binding VOCABULARY** — REGION, BOUNDARY, step, the mass pin, structurally pure-gDNA object — and it lists the banned synonyms. Read it before writing a comment or a docstring |
| **`docs/TESTING.md`** | The simulated panels, how to build them, the simulator's own gates, and what the suite can and cannot judge. ⭐⭐ **§0b is the TOY HARNESS** — how to define a test transcript and calibrate it in under a second; read it before hand-rolling any scenario |

⛔⛔ **AND `ROADMAP.md` IS COMPACTED RATHER THAN APPENDED TO — OWNER RULING, 2026-08-14.** It exists to
give a **clear set of directions for the next steps**, not a history of what was implemented. **When an
item is DONE, DELETE IT.** ⚠ The rule it replaces — *"a measurement is never deleted, only STAMPED"* —
grew the file to **787 lines**, where the completed work competed with the ranked list for attention.

⭐ **The safeguard is the MOVE rule, which is the same one that governs `docs/dev/`:** a finished item is
sometimes the only home for something durable, so before deleting one, put its **lesson** in `TRAPS.md`
(as a NAMED rule), its **ruling** in `DESIGN.md`, its **derivation** in `EQUATIONS.md`. A number that
describes where the tool **is** stays in `ROADMAP.md` §0. ⛔ **MOVE it — do not copy.** Copying is what
creates two homes that then diverge.

⚠ **What is NOT "completed work" and must survive**: §0's state (including a ✅ row whose job is to say
*do not work here*), and the **GRAVEYARD** of refused mechanisms — that one points forward (*do not
rebuild these*) and is a guardrail, not history.

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

`src/rigel/calibration/` is **38 modules, 13,216 lines** (re-derived 2026-08-17 — ⛔ do not carry this number, run `module_census.py`; the `39` it replaces was one of five counts in three homes that had all drifted). It is **not a knot** — measured from the AST there are no import
cycles, and 18 of them had exactly one importer when this was measured on 2026-08-07. It was a **FLAT PILE of peers**, which is the one shape
that cannot tell you where to add anything. `rigel/calibration/_layers.py` names the layers that were
already in the boundaries, and `tests/calibration/test_layering.py` enforces them.

**THE ONE RULE: an import may point DOWN a layer or SIDEWAYS within one. Never UP.** A module that needs
something from a higher layer is telling you the thing belongs lower — a TYPE almost always does.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally MEANS | **1 · the payload view** — `splice_graph` `substrate` `region_arrays` |
| how many places a fragment COULD have sat | **2 · opportunity** — `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `fl` |
| one slot's own numbers, and ψ | **3 · geometry + the per-slot solve** — `node_geometry` `simplex_logodds` |
| which strand a fragment came from | **4 · strand** — `gdna_strand` `strand_deconv` `strand_balance` `strand_summary`, and `strand_likelihood` (an executable REFERENCE, gated) |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `density_deconv` `npmle` `gdna_landscape` `background_reference` `run_fill` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `messages/` (the policy) + `node_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `diagnostics` `track` |

⭐ **Run `python scripts/design/module_census.py` rather than trusting this table** — it re-derives the
graph, prints every layer, and flags any import that points up. ⚠ It also flags **docstrings that name a
sibling with no import boundary**: 14 were measured on 2026-08-07 and **6 were genuinely stale**, including one
claiming four consumers that had one. A layer is *not* a claim that its modules are the right SIZE — layer 4
is five modules for one concept — and that question is deliberately still open.

## ⛔⛔⛔ THE THIRD POLICY IS BUILT — AND THE PANEL SAYS IT IS THE WORST ARM ON EVERY IN-SCOPE STRATUM

⛔⛔ **`src/rigel/calibration/messages/currency.py` is a DEVELOPMENT BASELINE WITH A MEASURED DEFICIT,
not an improvement, and the first thing to know about it is that its toy result did not transfer.** On
the **test chromosome** it beat both shipped policies on all three IN-SCOPE strata; on the **16-condition
ladder** it is the WORST of the three on all three, and `SilentPolicy` wins 8 of 16 rows
(`ROADMAP.md` §0 carries both tables). ⛔ **`TRAPS: panel-before-src`, for the fifth recorded time — a
claim about this policy MUST name its substrate.**

⭐ It sits beside `SilentPolicy` and `RelayPolicy` on the same gated backbone; the A/B is one config
value, `CalibrationConfig.message_policy`, and the default is `relay`, so nothing that works is affected.
⭐ **Where it IS best on the ladder**: both capture-ON zero controls and the dense unstranded capture-ON
rows — where the local solve is blind. That is the shape of what it does: it moves accuracy from
evidence-BEARING objects to evidence-free ones.

⭐⭐ **THE MECHANISMS BELOW ARE SOUND AND ARE WORTH KEEPING EVEN IF THE POLICY IS NOT** — each was found
by measurement, and two of them are about the tool as a whole rather than about this policy. The
handoff is `docs/dev/NEXT_SESSION.md`.

⭐⭐ **THE FOUR RULINGS BEHIND IT, so they are not re-litigated** — `DESIGN.md` §0c.0b–d and
`EQUATIONS.md` §3.5f–h: ① the ABUNDANCE and COMPOSITION strategies are **one continuum**, and the point
on it is a FITTED shrinkage (`w = (log r)²/((log r)² + v)`), never a switch; ② the **recipient** decides,
from a static population table keyed per (boundary, side) — a boundary carrying another transcript's TSS
is still population-equal with the flank that gains nothing; ③ **an imputation is a weak predictor by
definition** and pays a fitted PREMISE on every hop, because a message layer whose every variance shrinks
with counts transports for free and tramples a measurement; ④ an enrichment ratio is **model-free** — it
comes from the accumulator's reciprocal-opportunity banks, never from `mass / effective_length`.

⛔⛔ **AND THE MEASUREMENT DISCIPLINE IS PART OF THE RULING: split any message-layer error by whether the
destination HAD ITS OWN COMPOSITION EVIDENCE.** "The messages help" and "the messages trample a
measurement" are different findings and a pooled number cannot tell them apart — it is how an 8.09×
degradation of the measured half of a chain was found inside a modest total.

## ⚠ HISTORY — the rebuild plan that produced the policy above. Read as PROVENANCE, not as a queue.

⛔⛔⛔ **READ THIS BEFORE FIXING ANYTHING IN `messages/relay.py`. THE OWNER RULED (2026-08-18) THAT THE
RELAY IS REBUILT, NOT REPAIRED BUG BY BUG** — *"the whole method from the top down is not simple, and it
is not elegant … it can be simple and elegant"*. The per-bug loop below is HISTORY: it converted several
mysteries into mechanisms, and every one of those mechanisms is now a design INPUT to the rebuild rather
than a ticket against `RelayPolicy`.

⭐ **THE PLAN, and it is the frame for all message-propagation work**:
`docs/dev/PLAN_message_propagation_rebuild.md`. Build a **THIRD policy** beside `SilentPolicy` and
`RelayPolicy` on the same gated backbone, driven by a **currency-per-hop-type table** — a hop carries a
LEVEL or a COMPOSITION, and the two have exactly complementary invariances (a gDNA LEVEL survives a
POPULATION change and dies at an ENRICHMENT change; a COMPOSITION is the reverse). ⛔ `RelayPolicy` is NOT
touched, so nothing that works can break and the A/B is a config value rather than a diff.

⭐⭐ **WHERE IT IS: stages 0 (VOCABULARY) and 1 (UNBLOCK) are DONE. The next step is STAGE 2, "THE MAP" —
a MEASUREMENT with no `src/` change**: re-run the hop census and the currency oracle on all 32 conditions
and let the data name the currency for each of the six hop types, BEFORE a line of policy code is written.
⚠ Every number in the plan's §2 was measured on the PRE-FIX caches and is stale by construction — the
chimera repair moved boundary crossings by 2.4 % — so Stage 2 is a re-measurement, not a review.

⛔ **DO NOT** start a per-bug dissection of `RelayPolicy`'s zero-control gap. It is a real, named, open
defect (below) and it is deliberately NOT the next thing: the policy that carries it is being replaced.
Its mechanism is preserved as a CONSTRAINT the new policy must satisfy.

### The history below is the measured standing of the relay being replaced — read it as INPUT, not as a queue

`CalibrationConfig.message_propagation = True` (owner, 2026-08-18), after ~11 days muted. ⛔ **The
2026-08-07 mute and its −58/−44/−32/+155 % table are HISTORY, not the current state** — measured with a
confirmed licence bug live, on a panel retired 2026-08-13. ⚠ **The LAST MEASURED standing — 2026-08-18, on the PRE-chimera-fix caches.** The repair moved 2.4 % of
every boundary crossing, which is exactly what the relay hops on, so every number below is due a re-price
on the rebuilt panel (`relay_pool_ab.py`):

* ⭐ the relay is a **large win** on unstranded × capture-ON (0.31× vs muted) and at every zero-gDNA
  control (0.04–0.45×), and it is what solves a slot with NO own evidence at all;
* ⛔ it still **costs the three IN-SCOPE contaminated strata ~1.5–1.65× vs muted** — down from ~2.0×
  before the per-strand licence landed. ⛔ "Closing bug by bug" WAS the plan and is not any more — see the
  ruling above.

⭐⭐ **THE DISSECTION LOOP IS THE METHOD AND IT HAS CONVERTED THREE MYSTERIES INTO MECHANISMS IN ONE
DAY**: the sliver-dominated background wall (fixed — the Gamma(Σg+½, ΣE) posterior), and the reframe's
composition smuggling at termini (fixed — the PER-STRAND population licence, `rna_level_scale`:
population intact on both strands ⇒ reframe by `r`; this strand intact, other changed ⇒ `r = 1`, a
DENSITY transfer; this strand changed ⇒ NO CLAIM, value and precisions zeroed). Owner's ruling: a
boundary where the RNA populations differ — TSS/TES on either strand, strand-activity changes — may
transfer a DENSITY, never a composition.

⭐ **THE ZERO-CONTROL GAP, ITS MECHANISM NAMED — now a CONSTRAINT on the new policy, NOT the open queue**
(2026-08-18, session 2): the per-strand licence's `g00` regression (154 k → 431 k) walked to ONE shape — a zero RNA claim WITH
LIVE PRECISION, relayed one hop, where the mass pin's licence reads both arms as supplied and rescales
the whole budget onto gDNA (`k = M/(tg·E_g)`, 467,000×) — "all your mass is gDNA" into a 26 k-fragment
exon. Two feeders: ① the SPLICE IN added the sj count's precision back onto an arm the licence had just
zeroed — **FIXED, both twins, gated** (`TRAPS: zero-the-precision-with-the-value`): in-scope
**0.931 / 0.999 / 1.000**, `g00` **431 k → 324 k**, all 8 zero-control rows improve, deferred 1.000;
② the mis-scoped `struct_lock` (`~solvable`, the strict xfail) gives an EMPTY exon "RNA = 0 @ 0.2026"
on every arm — **OPEN**. ⛔ The xfail's fix (`g1_locked ∧ REGION`) and a zero-mass-source pin guard were
BOTH priced on all 32 and both fail the written bar with NAMED losses (`ROADMAP.md` §1 rank 1) — the
sharp predicate must refuse a source's OWN zero-count artefact and keep a RELAYED composition passing
through an empty slot (capture's conduits). ⚠ `g00` is ONE-SIDED: any RNA-favoring message reads as
"right" there, and the empty exons' "gDNA = 0 @ 0.2026" is such a coincidence — the graveyard's pattern.
⭐⭐ **THAT MECHANISM IS NOW A CONSTRAINT ON THE NEW POLICY, NOT A TICKET AGAINST THE OLD ONE**: a refused
claim must lose its precision with its value, and no operator may hand it back
(`TRAPS: zero-the-precision-with-the-value`). Every rung of the rebuild passes BOTH zero controls before
the next begins, which is how the new policy is prevented from re-acquiring it.

⭐⭐ **THE VOCABULARY RENAME IS RUNNING, AND IT IS THE OWNER'S HIGHEST PRIORITY (2026-08-18).** Landed:
`graft` → **SPLICE IN** (`splice_in`, tree-wide — verified unambiguous) and the MESSAGE-OPERATOR sense of
`p‑e‑e‑l` → **SPLICE OUT** (`splice_out`, **SENSE-SCOPED**). ⛔ That token also named the DECONVOLUTION verb
("strip the gDNA off, RNA is the residual") in `density_deconv` / `calibrate` / `object_composition` and in
the rule now called `TRAPS: the-deconvolution-is-as-good-as-the-density-it-is-handed`; a global replace
would have corrupted it (`TRAPS: two-masks-one-name` — run `rename_census.py --sense <token>` BEFORE
renaming anything).
⭐ **ALSO LANDED 2026-08-19**: that second sense → **`deconvolve`** (owner: *"peeling is for
fruits"*), carrying the renamed rule `TRAPS: the-deconvolution-is-as-good-as-the-density-it-is-handed`;
`mass_pin` → **`mass_rescale`** and the capture key `_pin` → `_rescale` — SCOPED, because `pin` also means
the ordinary verb (a pinned seed, a pinned thread count); `HeadPolicy`/`HeadSwitches`/`head.py` →
**`RelayPolicy`/`RelaySwitches`/`relay.py`** — scoped by IDENTIFIER, so `_bam_header` and git HEAD are
untouched; `lend` → **`may_share_composition`**; `s2t` → **`hop_logvar`**.
⛔⛔ **STILL TO RENAME — `arm` ALONE, DEFERRED FOR A MEASURED REASON: IT HAS THREE SENSES.** An EXPERIMENT
arm (`--arm`, `arm_score.py`, `ladder_arm_ab.py`); a COMPONENT arm ({gDNA, RNA+, RNA−} — the sense the
owner ruled becomes `component`); and ⛔ **`__ARM_NEON` in the C++ scanner**. 1,358 sites. It needs its own
`--sense` pass, never a tail-end sweep.
⚠ **AND A RENAME REWRITES THE LIST OF THINGS TO RENAME** — this paragraph was itself rewritten by the pass
it describes (`mass_pin` → `mass_rescale` inside the sentence saying to rename it). Write such a list with
the tokens in a form the pass cannot match, or re-read it afterwards. ⛔ The land/sea analogy is an EXPLANATION and must never enter code or docs (owner, 2026-08-18) —
it is at zero occurrences and must stay there. `docs/dev/PLAN_message_propagation_rebuild.md` §5 is the
review; `rename_identity.py` is the gate.
## ⛔⛔ CITE A RULE BY ITS NAME. NUMBERED LABELS ARE BANNED

`TRAPS.md`'s rules used to be `A16`, `D4j`, `C0b`. They are **named**: cite them as
`TRAPS: off-grid-message-mode`, `TRAPS: a-cancelling-defect-pair`,
`TRAPS: frame-free-is-not-assumption-free`. ⭐ The name IS the identifier, so a citation says what it means
without a lookup and is still one greppable string with one home.
⛔ `tests/test_no_jargon_labels.py` enforces it and its allowlist is scoped, never blanket.

⚠ **The numbers were not merely opaque, they were AMBIGUOUS.** `A1` was a validation trap here *and*
`SUCCESS.md`'s FIDELITY criterion; `C1` was a pool trap *and* a moment variable in
`calibration/length_likelihood.py` — ⚠ a module PURGED since (`b7ed7a0b`): the moments were mis-filed, not
dead, and live in `effective_length.py`, so no fragment-length composition channel exists in `src/`;
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
  2026-08-05). ⛔⛔ **AND CHECK WHICH RULER THE CEILING LEFT INSTALLED.** Every measurement arm to date
  patches `assemble_priors`, and `pipeline.py` builds `effective_lengths_em` **before** it calls
  `assemble_priors` — so the effective-length shrinkage has never been inside any ceiling, and every
  ceiling number measured so far had the shipped (wrong) ruler still in place (2026-08-14). The loop that
  produced the 39 % win is: run the panel → take the worst **IN-SCOPE** scenario (0.8.0 defers
  unstranded × capture-ON, which is otherwise the worst one every time) → dissect it to the single
  highest-error object → find the cause → fix → REPEAT. Reach for that first.
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
python -m pytest tests/ -q                     # baseline: 0 fail / 3554 pass / 0 skip / 9 xfail
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

⭐ **THE STANDING BASELINE IS GREEN: 0 failures, 3,652 passing, 0 SKIPPED, 9 xfail, 3,661 collected**
(re-derived 2026-08-21, after **THE NPMLE RETIREMENT** — converge-and-delete, gated on a BIT-IDENTICAL
deliverable and verified: `rename_identity.py --check` reads
`✅ BIT-IDENTICAL — quant digest and every array's content unchanged`). ⚠ **The `−19` from 3,680 is
ACCOUNTED as a DELETION**: `−16` `tests/calibration/test_npmle.py`'s own cases with the file, `−2` its
meta pair (jargon + docs-boundary — a `tests/calibration/` file is `+2`, never `+3`), `−3` for
`src/rigel/calibration/npmle.py` (jargon, docs-boundary, **and layering**, which it only had because it
was declared in `_layers.py` — removed with it), `+1` in `test_report.py` (the one `from_prior` gate
became two `from_abundance_landscape` gates) and `+1` in `test_abundance_landscape.py` (the refusal gate
became a skip gate plus an asymmetry gate). ⛔⛔ **AND THE DEFAULT FLIP COULD NOT LAND AS
`PLAN_measured_prior.md` §3d WROTE IT — measure before flipping.** `abundance_landscape = True` with the
old REFUSAL broke **65 callers** (24 failed + 41 errors), all one cause: unit and toy fixtures have no
wall arrays and never wanted a QC panel. The refusal's own reason was *"rather than fitting on unmasked
totals"*, and the alternative was never to fit unmasked — it is not to fit; so the missing-walls case now
SKIPS with a WARNING and yields `None` (the report omits the panel, its pre-existing contract).
⛔ `background_abundance` KEEPS its refusal — that pair feeds ψ — and a gate now asserts the asymmetry so
the two cannot be conflated. ⚠ The previous line read 3,670 passing / 3,679 collected. ⚠ **The `+9` over 3,670 collected is ACCOUNTED as `8 own + 2 gates − 1 deleted doc`**: the two new
`scripts/design/` instruments (`landscape_head_to_head.py`, `transfer_variance_audit.py`) are **+4 each**
— the measured `+4` rule for a `scripts/design/` file, confirmed a fifth and sixth time (imports,
says-what-it-is-for, jargon, docs-boundary) — plus **+2** own cases appended to the EXISTING
`tests/calibration/test_abundance_landscape.py` (the `knn_scale` pass-through and the exported
`split_basins` rule, so no per-file meta drift), minus **−1** for `docs/dev/NEXT_SESSION.md`, deleted per
its own instruction once its findings had MOVED. ⛔ Derive it with
`pytest --collect-only -q | grep -E 'landscape_head_to_head|transfer_variance_audit'`, which prints
exactly eight lines. ⚠ **The `test_scripts_index` gate fired once during this work and it was RIGHT**: both
new instruments were on disk and absent from the table below, which is the "a script nobody indexed is a
script the next session rebuilds" failure; rows were added rather than the gate exempted. ⚠ The previous
line read 3,661 passing / 3,670 collected. ⚠ The
`+25` over 3,644 is `16 own + 9 meta`, re-derived with exact bracket matching:
`tests/calibration/test_abundance_landscape.py` carries its 16 own cases plus **+2** (jargon,
docs-boundary); `src/rigel/calibration/abundance_landscape.py` is **+3** (jargon, docs-boundary, and
the layering case that exists only because it is declared at layer 5 in `_layers.py`);
`scripts/design/abundance_landscape_census.py` is **+4** (imports, says-what-it-is-for, jargon,
docs-boundary). ⚠ The retirement census also found and fixed a LIVE breakage the suite was green over:
`toy_harness.describe()` still read the `background` field deleted earlier that day
(`TRAPS: a-green-suite-hid-five-dead-instruments`, the instrument form — found by a subagent census,
not by any gate).
⚠ It read **3,635 passing / 3,644 collected** before that (`background_reference` converge-and-deleted). ⚠ The `−19` from 3,654 is
a DELETION and is accounted: `−13` `test_background_reference.py`'s own cases with the file, `−5`
`test_npmle.py` fixtures that exercised the dead `background=` parameter, `−1` my own
`test_measure_background_takes_the_PAIR_too`, `−2` meta (the deleted module and its deleted test file
each carried a jargon + docs-boundary case), `+2` because the src module and test file removals also
drop their layering/parametrised entries — net `−19`, confirmed with `--collect-only`. ⭐ The one gate
that had to MOVE rather than go: `test_the_shipped_default_is_BIT_IDENTICAL_and_the_flag_is_NOT_INERT`
asserted on the deleted aggregate `background`; it now asserts on `intron_background`, which is the
background that actually reaches ψ — where it should have pointed from the start.
⚠ It read **3,654 passing / 3,663 collected** before that (the two background consumers swapped). ⚠ The `+7` over 3,647 is `6 own + 1 meta`: six swap gates
appended to the EXISTING `tests/calibration/test_total_abundance.py` (so no per-file meta drift) plus
one jargon case for one new `docs/dev/` file. ⚠ **The docs-boundary gate fired once during that work
and it was RIGHT**: `ROADMAP.md` had been given a path into `docs/dev/`, which is the one thing a
permanent doc may never do; the citation was rewritten to carry the finding instead of the path.
⚠ It read **3,647 passing / 3,656 collected** before that (the MEASURED TOTAL's validation rung). ⚠ **The `+32` over 3,615 is
ACCOUNTED as `23 own + 9 meta`, and the meta split confirms all three documented per-file rules at
once**: `tests/calibration/test_total_abundance.py` contributes its 23 own cases plus **+2** (jargon,
docs-boundary — a `tests/calibration/` file, never `+3`); the new src module
`src/rigel/calibration/total_abundance.py` is **+3** (jargon, docs-boundary, and ⭐ **layering**, which
exists only because it is declared in `_layers.py`); the new instrument
`scripts/design/total_abundance_audit.py` is **+4** (imports, says-what-it-is-for, jargon,
docs-boundary). ⛔ Derive it with `pytest --collect-only -q | grep total_abundance` and EXACT bracket
matching — the src module and the test file share a stem prefix. ⚠ **The jargon gate fired once during
this work and it was RIGHT**: the new test file's prose labelled two boundaries `B0`/`B1`, which is the
numbered-label ambiguity the rule exists to stop; the wording was fixed rather than exempted.
⚠ It read **3,615 passing / 3,624 collected** before that (the START/END/SPAN accumulator banks). ⚠ The
`+6` over 3,618 was the six new gates in `tests/native/test_accumulator_spec.py` (an existing file, so
no meta drift), each verified failing first; the C++ was additionally perturbed once (END deposited at
the start region) and the byte-identity parity gate fired, proving the new banks are inside it.
⚠ It read **3,609 passing / 3,618 collected** before that (end of the abundance/cleanup day). ⚠ The `+10` over 3,608 is ACCOUNTED:
`+4` meta for `exon_solvability.py` becoming tracked and `+1` for the weighted-premise gate (both
landed in `f295a313`, whose commit did not carry the number back here); `+2` dev docs (the dissection
session record and the spike-and-slab primer) `−1` (`NEXT_SESSION.md` deleted per its own
instruction); `+1` the truncation gate in `test_fragment_length_proof.py` (fragments OUTSIDE the
support read `rho·P(w<=ell)`, never `rho`); `+1` the `inv_abundance` fill gate in
`test_region_geometry.py`; `+2` more dev docs (the model-free-abundance plan and the spike-and-slab
plan). ⛔ Commits `a2b81b34` and `1ada4a3c` each RECORDED "3,618 collected" while their trees collect
**3,617** — an arithmetic slip (3,608+9), corrected here rather than by rewriting history
(`TRAPS: re-record-the-baseline`).
⚠ It read **3,599 passing / 3,608 collected** before that (the third policy's landing). ⚠ **The `+45` over 3,563 collected is ACCOUNTED, not
adjusted: `39 own + 6 meta`, from the THREE files this work added.** The 39 are
`tests/calibration/test_currency_policy.py`'s own cases — the policy is gated concept by concept and
every gate was fired by perturbation. The 6 are **jargon ×3** (the test file, `messages/currency.py`,
and the dev doc), **docs-boundary ×2** (the test file and the src file — a doc gets no docs-boundary
case), and ⭐ **layering ×1**, which is the one a reader would miss:
`test_no_import_points_UP_a_layer[messages/currency]` exists only because the new module was declared in
`_layers.py`. ⚠ Derive this with EXACT bracket matching — a bare `grep currency.py` also matches
`hop_currency.py` and over-counts by four. ⛔ **The DOC sweep that moved this session's findings into
`DESIGN.md`/`EQUATIONS.md`/`TRAPS.md`/`ROADMAP.md` moved the collected total by ZERO, and that is the
check** — content-only edits to existing files add no cases. Re-derive with
`pytest --collect-only -q | grep test_currency_policy`; never adjust.
⚠ **It read 3,546 / 3,555 collected before that** (the method-development test reference).
⚠ Before that: the per-arm hop map, the per-strand oracle partitions and `build_test_reference.py`
(`+4`, one new `scripts/sim/` file); and the index-driven simulator at 3,538 / 3,547.
⚠ The index-driven simulator's own `+2` was two new gates in `tests/test_sim_capture.py`. ⛔ **That commit MOVED FOUR GOLDEN SCENARIOS** — `nrna_moderate_ss90`,
`nrna_heavy_ss90`, `combo_moderate`, `combo_extreme` — regenerated with the diffs read first
(`effective_length` ≤ 0.7 %, `nrna_em_count` +6 to +20 %, `gdna_em_count` moves on the nascent-only
scenarios): the nascent model changed by design, so the simulated reads changed.
⚠ It read **3,536 / 3,545 collected** before that (Stage 2's `hop_currency.py`, the measured `+4` rule
for one new `scripts/design/` file, no `src/` change and no golden moved).
**Any failure at all is a regression** — a stronger and cheaper rule than counting expected ones.

⚠ **It read 3,532 / 9 xfail / 3,541 collected before that** (the SPLICE IN keeps its precision WITH its
value at a refused hop): the `+2` over 3,530 is the two SPLICE IN gates (VECTORISED + SCALAR) in
`test_gdna_scale_rule.py`, an EXISTING file so no meta-test drift; that fix moved NO golden.

⚠ **It read 3,521 / 9 xfail / 3,530 collected before that** (the PER-STRAND population licence): the `+3`
over 3,527 was the three-case licence contract's unit gates in `test_gdna_scale_rule.py`; the licence
work's golden moves were 1e-6–6e-4 relative, read before regenerating.

⚠ **It read 3,501 / 10 xfail / 3,511 collected before that, and every move is ACCOUNTED, not
adjusted.** `+16` collected = `+12` for three new `scripts/design/` instruments
(`calibration_oracle.py`, `calibration_walk.py`, `relay_pool_ab.py` — the measured `+4` rule, three
times) `+4` for `test_density_deconv.py`'s new falsification gates (sliver-invariance,
smooth-through-zero, empty-pool-honesty, populated-limit). ⭐⭐ **xfail 10 → 9 is `−2 +1` and BOTH
directions are the precedent kinds**: the two recorded-price-of-the-mute xfails
(`test_ambig_scenario`, `test_toy_harness`'s intron gate) went GREEN the day
`message_propagation = True` landed — the exact route their reasons named — and their measurements
moved into their docstrings; `test_nrna_multiexon_t2_low_ss` became a strict xfail recording the
prior-assembly `alpha = 0` nascent rule (owner diagnosis 2026-08-18), a casualty of a net-helpful
change whose repair belongs in the assembler, not calibration.

⚠ **It read 3,496 / 11 xfail that morning and BOTH moves are ACCOUNTED, not adjusted. The COLLECTED
total is unchanged at 3,508, which is the real check.** ⚠ It then moved to **3,511 collected** as the
instrument work landed: the gate was EXTENDED to `scripts/profiling/` (a third tree, previously covered by
nothing) while two dark `scripts/design/` instruments and one executed dev doc were deleted. Re-derive
with `--collect-only`; never adjust.
* **`+1` passing** is ONE `test_no_jargon_labels::test_no_numbered_rule_labels` case for one new
  `docs/dev/` file, re-derived with `pytest --collect-only -q | grep <stem>`, which printed exactly one
  line. The dev-doc `+1` rule is confirmed a third time.
* ⭐⭐ **`+1` passing / `−1` xfail is AN XFAIL CLOSED BY REPAIRING THE THING, which is the precedent this
  project keeps**: `test_scripts_index.py` calls `pytest.xfail` only while a script is BOTH listed in
  `BROKEN_ON_IMPORT` and genuinely broken. `prior_units_check.py` was repaired, so the gate stopped
  xfailing and — correctly — **FAILED on the stale exemption** until the entry was removed. That dict is
  now EMPTY and should stay so.
⛔ **The doc and instrument cleanup moved the collected total by ZERO and that is the check**: edits to
comments, docstrings and prose add no cases, so a content-only sweep of ~50 files must leave it where it
was. A moved total after a content-only sweep means a file was added or removed.

⛔⛔ **AND DERIVE THE FAILURE SET, NEVER EYEBALL THE TAIL** (`TRAPS: read-the-whole-failure-list`). pytest
prints the last screen; a long summary scrolls the rest away. 21 failures were once reported here as
"21 golden failures" from the visible eleven — **two were real**, and one of them was the most informative
result of that session. Use::

    python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c

⚠ **It read 3,476 / 11 xfail before that, and the delta is `+20 = 17 own + 2 meta + 1 dev doc`, all of it
`tests/calibration/test_composition_closes.py`** — re-derived with
`pytest --collect-only -q | grep test_composition_closes`. ⛔ That commit also regenerated **21**
`tests/golden/` scenarios; the diffs were read BEFORE regenerating (`em_effective_length` up to 19.4 %,
`count` ≤ 1.1e-2) rather than waved through.

⚠ **It read 3,459 / 10 xfail before that, and the delta is `+17 passed / +1 xfail = +18 COLLECTED`, all
of it `tests/calibration/test_structural_reference.py`** — 16 of its own cases (15 passing plus the strict
xfail that records the relay defect the flag EXPOSES rather than causes) plus **2** meta ones, re-derived
with `pytest --collect-only -q | grep structural_reference`, which prints exactly
`test_no_numbered_rule_labels` and `test_nothing_outside_the_sandbox_cites_into_it`. ⚠ The `+2` rule for a
`tests/calibration/` file is now confirmed twice.

⛔⛔ **AND FLIPPING THAT DEFAULT MOVED `tests/golden/` — 6 scenarios, regenerated with
`--update-golden` and the magnitudes recorded rather than waved through**: `em_effective_length` on the two
gDNA scenarios (the ruler is downstream of composition) and `count` on the four others, max relative
difference **0.00026 … 0.049**. ⚠ A golden update is where a regression is laundered into "intended", so
the diff was read first and the truth-scored instruments (panel, thermometer) checked before regenerating.

⚠ **It read 3,420 before that work; the `+39` is ACCOUNTED in two parts.** `+38` is
`tests/calibration/test_reference_location.py` — **36** of its own cases plus **2** meta ones (jargon and
docs-boundary); `+1` is the handoff dev doc, which the jargon gate covers. ⚠ **A `tests/calibration/`
file is `+2` meta, not the `+3` a top-level `tests/` file gets** — `test_scripts_index` does not
parametrise over it.

⚠ **It read 3,420 before that; the `+38` is `tests/calibration/test_reference_location.py` — 36 of its
own cases plus 2 meta ones (jargon and docs-boundary), re-derived with
`pytest --collect-only -q | grep reference_location`.** ⚠ A `tests/calibration/` file is `+2` meta, not
the `+3` a top-level `tests/` file gets: `test_scripts_index` does not parametrise over it.

⚠ **It read 3,413 before that work and the `+7` is ACCOUNTED in two parts, each re-derived with
`pytest --collect-only -q | grep object_composition`.** `+4` is `scripts/design/object_composition.py`
(`test_every_instrument_still_imports`, `test_every_instrument_says_what_it_is_for`,
`test_no_numbered_rule_labels`, `test_nothing_outside_the_sandbox_cites_into_it`) — the `+4` rule below
confirmed a third time; `+3` is `tests/test_object_composition_self_test.py` (its own case plus the
jargon and docs-boundary ones). ⚠ The `3,413` it replaces is itself `3,412` plus the **dev doc**
`docs/dev/rna_arm_design.md`, which the previous session added and did not carry back here — a dev doc
is `+1` (the jargon gate covers `docs/dev/`; the docs-boundary gate does not).

⚠ **It read 3,404 until the ruler work landed on 2026-08-14; that `+8` is ACCOUNTED in three parts,
each re-derived with `pytest --collect-only -q | grep <stem>` rather than reasoned about.** `+4` is
`scripts/design/calibration_vs_oracle.py` (`test_every_instrument_still_imports`,
`test_every_instrument_says_what_it_is_for`, `test_no_numbered_rule_labels`,
`test_nothing_outside_the_sandbox_cites_into_it`); `+3` is
`tests/test_calibration_vs_oracle_self_test.py` (its own case plus the jargon and docs-boundary ones);
`+1` is the dev doc, which the jargon gate covers and the docs-boundary gate does not.
⛔ **A new `scripts/design/` file is worth `+4`, not the `+3` this file used to estimate** — measured
2026-08-14 against four existing instruments, all four.

⚠ **It read 3,397 until earlier on 2026-08-14 and that number was stale by `+7`, ACCOUNTED in two parts.**
`+5` was already on disk: the panel rebuild (`59a341e0`) measured **3,402** in its own commit message
and did not carry the number back here. `+2` is the doc rewrite, which added
`docs/dev/2026-08-13_panel_preconditions_report.md` and
`docs/dev/session_2026-08-14_calibration_handoff.md` — one `test_no_jargon_labels` case each, confirmed
with `pytest --collect-only -q | grep <stem>`. ⚠ **A commit that measures the suite must update this
line**, or the next session reads a green run as a regression.

⚠ **It was 3,298 / 10 xfail on 2026-08-12; the `+99` to 3,397 was ACCOUNTED, not adjusted, in two steps.** `+56`
is the per-transcript RNA prior foundation (`test_grouped_prior_update` 24, `test_transcript_path` 8,
`test_warm_start` 6, 7 composition gates in `test_result_schema`, and the file-parametrised meta-tests
gaining one case per new file). `+43` is this session: `+6` for the conserved sj mass
(5 in `test_result_schema`, 1 in `test_calibrate`, and **no** meta drift because it added no file),
`+31` for `test_transcript_weights.py`, and `+6` meta cases from the two new files — two each in
`test_docs_boundary`, `test_no_jargon_labels` and `test_scripts_index`.
⭐ **Re-derive it the same way**: `pytest --collect-only | grep <new file stem>` prints exactly which
meta-tests gained a case, which is faster and more honest than reasoning about the multipliers.

⛔⛔ **AND THAT GATE EXISTS BECAUSE A GREEN SUITE HID FIVE DEAD INSTRUMENTS — the rule and its two
recorded mechanisms are `TRAPS: a-green-suite-hid-five-dead-instruments`, which now has a BODY.** ⚠ This
paragraph used to carry that rule's text; it was the two-homes failure the MOVE rule exists to stop, so it
has been moved rather than copied. ⭐ What you need here is only the trigger list: run the suite **and the
instruments** after a `src/` DELETION *and* after a **CONFIG DEFAULT FLIP** — the second was found on
2026-08-17, when six instruments were dead under the shipped `message_propagation = False` while the suite
stayed green, because the TEST readers install `RelayPolicy` themselves.

⛔ **RE-DERIVE THE COUNT, NEVER ADJUST IT** (`TRAPS: re-record-the-baseline`). Several suite tests are
PARAMETRISED OVER DOC, SOURCE AND SCRIPT FILES — roughly 2 per src module, 3 per script, 1 per doc — so
adding or retiring a file moves the total by a few. Account for the delta; a hand-carried number has
drifted every time it has been carried. ⚠ Four successive versions of this paragraph were stale when
read, which is why it no longer records what they said.

⭐⭐ **7 xfails, and they are NOT one kind of thing** — interrogated 2026-08-11, so do not treat the list
as uniform:

* **2 were GREEN and one config flag broke them.** `test_ambig_scenario` and `test_toy_harness`'s intron
  gate both go green the instant `message_propagation = True` (verified: ambig `0.4577 → 0.000247`). They
  are the *recorded price of the switch*, not defects. ⛔ Flipping it costs **+154.8 %** on
  unstranded × capture-ON — 0.8.0's DEFERRED stratum, which carried 73 % of panel error on the retired
  36-condition ladder — the owner's ruling, not a test decision.
* **5 were WRITTEN as xfail and were never green** — executable records of two PROVEN defects found by
  the dissect loop: `struct_lock` is `~solvable & REGION` rather than `g1_locked & REGION`, and `Var(f_g)` is
  capped at `f_g(1−f_g)`, which is 0 at the `f_g = 1` default of an evidence-free slot. ⛔ For these
  "fix the test" is a CATEGORY ERROR: the test is right and the code is wrong. Rescoping `struct_lock`
  alone measures **−1.2 %** on its target stratum, **+0.4 %** worse at `g98` and **+3,207 %** on the
  zero-gDNA control — it is half of a cancelling pair (`TRAPS: a-cancelling-defect-pair`) and must be
  priced WITH the boundary-composition arm or not at all.

⛔ **Two were closed on 2026-08-11 and the method is the precedent:** the R1-sense gap by BUILDING the
missing capability, and the paralog collapse by asserting the collapse STRUCTURALLY (`min == 0`,
`max == total`) — which still fires the day the tie is broken AND catches a partial break a strict xfail
could not see. ⚠ Neither was closed by widening a bound.

Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Tooling under `scripts/`

`docs/SUCCESS.md` lists these in the order you would run them.

⛔⛔ **THE TABLE IS AN INDEX OF `scripts/design/` PLUS FOUR `sim/` ROWS AND NOTHING ELSE — DO NOT READ AN
ABSENCE AS A DELETION.** It used to say *"everything not listed here was deleted as unrunnable"*, which was
false in three trees. Re-derived 2026-08-17: **73 `.py` under `scripts/` = 61 `design` + 8 `sim` +
3 `profiling` + `__init__`**.
* **8 `design/` files are on disk and unlisted** — `bam_spans.py`, `implicit_splice_census.py`,
  `inv_L_limits.py`, `length_sieve.py`, `partition_test.py`, `reference_on_real_data.py`,
  `sigma_inv_L.py`, `spanning.py`. ⭐ All eight were re-run this campaign and all eight COMPLETE;
  `bam_spans.py` carries a PROMOTE recommendation.
* **4 `sim/` files are indexed nowhere**: `build_toy_2exon_reference.py`, `locus_sweep.py`,
  `simulate_suite.py`, `snapshot_suite.py`.
* ⛔⛔ **THE WHOLE OF `scripts/profiling/` IS COVERED BY NO GATE** — `tests/test_scripts_index.py`'s
  `ALL_SCRIPTS` is `design + sim` — and **2 of its 3 files were dead** when that was measured
  (`pyspy_driver.py` on 2026-08-11, `scan_profile.py` on 2026-08-17). That tree has now rotted twice with
  the defect class the existing gate already catches, only out of its reach. **Extend the gate to a third
  tree, or delete the directory** — an owner decision, recorded in `scripts/README.md`.
⚠ **And a basename collides across trees**: `design/scan_profile.py` (indexed) vs
`profiling/scan_profile.py` (unindexed). Both are live; renaming one is an owner call.

⭐⭐ **The GROUPS are ordered by 0.8.0 priority** — calibration against oracle calibration first, then the
panel and its caches, then the thermometer above them. A row's own ⭐s rate the INSTRUMENT; the group says
whether this phase is pointed at it.

| | |
|---|---|
| **⭐⭐⭐ START A SESSION HERE** | |
| `design/preflight.py` | ⭐⭐⭐ **CAN THIS SESSION RUN AND REGENERATE EVERYTHING? — one command, one verdict, before anything else.** Owner, 2026-08-19: *"I want the next session to verify that they have access to all of the necessary tools and can rerun and regenerate everything."* Checks the TOOLCHAIN (the `rigel` env, the compiled native extension, the CLI), both REFERENCES (the panel's and the method-development test chromosome's, hand-edited sources vs DERIVED artifacts), both PANELS (scan caches, oracle caches with all FIVE partitions — the three `ORIGINS` plus `rna_pos`/`rna_neg` — and the certified `slot_truth`), and every `scripts/design/` instrument (imports, then `--self-test`; `--fast` skips the ~2 min sweep). ⛔ **It changes nothing and measures nothing** — every check is a read or an import, and a ✘ prints the exact command that regenerates the thing, because a missing DERIVED artifact is a command you have not run rather than damage. ⭐ Its first real run caught the 8 test-chromosome scenarios sitting uncertified. `--self-test` 8/8, the reporter perturbed |
| **⭐⭐⭐ 0.8.0'S METRIC — calibration against ORACLE CALIBRATION** | |
| `design/calibration_vs_oracle.py` | ⭐⭐⭐ **THE CALIBRATION RESULT ITSELF AGAINST AN ORACLE CALIBRATION — and THE ONLY INSTRUMENT THAT REACHES THE EFFECTIVE-LENGTH SHRINKAGE.** `P = calibrate(...)` against `O = dataclasses.replace(P, **override_masses(ra))`, same payload, only the six deconvolved arrays swapped. ⛔⛔ **A `CalibrationResult` HAS TWO CONSUMERS and every other arm in the tree patches the second one**: `transcript_capture_eff_lengths` (the EM's RULER) is built by `_setup_geometry_and_estimator` **before** `assemble_priors` runs, so the shrinkage had never been inside any measurement arm. ⭐ **Its verdict, 2026-08-14: at the `g00` zero control the ruler reads factor 0.0951 where the truth is exactly 1.000** — 1.06 **billion** bp of opportunity, on a library containing no gDNA — while on the three contaminated in-scope strata `P/O` is **1.008 / 0.997 / 1.022**. ⛔⛔ **DO NOT READ THAT `P/O` AS "THE RULER IS FINE THERE" — IT IS THE MISTAKE THIS ROW EXISTS TO PREVENT.** The aggregate barely moves while nearly every transcript is redistributed, and substituting the oracle ruler on those strata is **2.5–2.7× WORSE** end to end (`ROADMAP.md` §0). Read `ruler_n_moved`, never the aggregate. ⭐⭐ Carries a fourth column no other instrument has: **`U`, the NO-ENRICHMENT NULL** — O's own gDNA total laid down at exactly uniform density, which at capture-OFF is the physically correct field and must read 1.000. It reads **0.9963–0.9967** while BOTH the shipped ruler (0.932–0.939) and the ORACLE one (0.919–0.924) sit far below — ⛔ **so at capture-OFF the shipped ruler is CLOSER to the truth than the oracle one, and `oracle_ruler` is an A/B between two wrong rulers rather than a ceiling.** The `U` ruler is the ceiling arm that still needs building. ⛔ Reports BOTH aggregates because they disagree by 3.6× at `g00` (`Σeff/Σfl` **0.095** vs the unweighted mean **0.345** — the contraction falls hardest on LONG transcripts); rank on `factor`. ⛔ Reuses `pass0_vs_oracle.score_axis` and `prior_vs_oracle`'s strata rather than adding a second scorer, and does NOT re-score `LocusPriors` — that is `prior_vs_oracle.py`'s. ⭐ **No solver, no EM, no BAM re-scan: ~5-12 s/condition, the whole ladder in ~2 min**, which is the iteration loop. ⭐ It is the FIRST instrument to stamp the 0.8.0 IN-SCOPE/DEFERRED ruling onto every row. ⚠ Reads the four `g00` full payloads from `scan_cache/` (the oracle cache holds no `_main` for a zero-gDNA row) and PRINTS which source it used; sum-to-full is re-run either way. `--self-test` perturbs every comparator with no I/O — 21/21, gated by `tests/test_calibration_vs_oracle_self_test.py` |
| `design/object_composition.py` | ⭐⭐⭐ **DOES ψ's Beta REFERENCE HAVE TO BE ONE LIBRARY-WIDE NUMBER? — measured NO, and a per-object one is worth 25×.** ψ solves one object at a time and the gDNA arm's fitted term is ALREADY per slot `(n_slots, K)`; the reference is the only scalar left. Written per object as `m_i = ρ_g,i·E_g,i / (ρ_g,i·E_g,i + ρ_r,i·E_r,i)`, every term but the two DENSITIES is per-object and exactly known. ⛔ **No solver, no EM, no re-scan.** ⭐⭐ Scored as `Σ|m_i − f_g,i|·M_i` — fragments the prior misplaces if believed outright — against the shipped constant ½, per stratum: **0.040 / 0.326 / 0.040 / 0.326** and **exactly 0.000** at both zero controls; `g05`, which regressed 1.43× under every library-wide mean tried, reads **0.009**. ⛔⛔ **The class-pooled TRUTH arm is NOT a ceiling** — the prior-free form beats it everywhere, because it hands an on-target RNA density to objects mature RNA cannot occupy. ⛔⛔ **AND THE WIN IS NOT PER-OBJECT RESOLUTION AT EXONS**: within `R exon` `m_i` has sd **0.0021** against a true `f_g` sd of **0.4441**, and within `B exon|exon` it is exactly constant. The four ANNOTATION-DETERMINED strata go to 0.000; exons/`exon|exon`/AMBIG are the population the gDNA LANDSCAPE serves — the two terms PARTITION the object universe rather than competing. ⛔⛔ **THE BOUNDARY AXIS SPLITS ON WHETHER MATURE RNA CAN CROSS, NOT ON WHETHER A sj ATTACHES** (owner, 2026-08-15): an earlier 'has a sj' stratum measured true `f_g` **0.0000 over 955,428 fragments** at the zero control because it lumps `exon|exon` in. The same `mrna_active` predicate then gates the RNA density — `R intron` 0.498 → **0.000**, `B exon|intron` 0.485 → **0.000**. ⭐⭐ **Hybrid capture needs no detection step**: the ratio of the off-target anchor to the in-gene `exon|intron` one IS the enrichment — **0.98** without probes, **113–114** with, no threshold and no flag; the off-target anchor is 0.91–1.00× of truth on all 16. ⛔ The in-gene anchor is a DETECTOR, not a calibrated level — it under-reads on-target gDNA **2.6–3.6×** under capture (it sits at the EDGE of the probe footprint). ⭐⭐ Table ⑧ is the RUNAWAY BOUND: the post-solve update iterated with each object's LIKELIHOOD REMOVED — a strict upper bound on the feedback — converges from four starts spanning three decades to the same fixed point to six decimals on all 12 contaminated conditions and to exactly 0.000 at the control. ⛔ `boundary_spliced` is a SEPARATE bank from `boundary_unspliced`, so S SUBTRACTS rather than bounding `f_g`; a first draft wrote `f_g ≤ 1 − S/M` and the truth violated it by **302**, and subtracting it inside `m_i` was then a net LOSS (0.158 → 0.072 without it). ⚠ `R intron` / `B exon|intron` reading exactly 1.0000 is `nrna = 0` restated — only `R intergenic` is nascent-free. Seven strata, partition ASSERTED; `--self-test` 25/25, gated by `tests/test_object_composition_self_test.py` |
| `design/abundance_landscape_census.py` | ⭐⭐⭐ **WHAT DOES THE TOTAL-DENSITY FIELD LOOK LIKE, PER CONDITION — the `AbundanceLandscape` census, and the instrument rung 4's inputs are read from.** Fits `calibration.abundance_landscape` on the cached payload's wall-exact measured totals (no solver, no EM) and prints every mode with its basin mass, `rho_0`, the span `R`, the pooled-anchor consistency gap in nats, and the per-region enrichment responsibility `w` BY CLASS (exposure-weighted, never thresholded). ⭐⭐ **First full run (2026-08-21, 40 conditions over 4 panels): every contaminated ladder row anchor-CONSISTENT** — gaps 0.003–0.063 nats capture-OFF (`rho_0` 0.0511 vs a pooled anchor of 0.051 at `g50 ss0.99`) and 0.25–0.65 under capture; the enriched basin is CLAIMED BY EXONS (w̄ₑ 0.35 vs intron 0.002 at `g50` capture-ON). ⛔⛔ **The first real fit corrected the plan's own gate wording**: capture-OFF totals are NOT unimodal — expressed exons sit decades above the gDNA background in any TOTAL — and the invariant that holds EXACTLY is the non-exon one, `w[intergenic] = w[intron] = 0.0000` on every contaminated capture-OFF row; at `g98` capture-OFF (2 % RNA) the field does read unimodal. ⭐ So the landscape alone cannot distinguish expression-bimodality from capture-bimodality, and rung 4 must pair `(rho_0, w)` with the measured enrichment detector (the anchor ratio) — settled in the plan before any prior code exists. ⚠ `g00` rows have an EMPTY anchor pool (zero intergenic counts ⇒ gap `nan`, fallback stamped); ⚠ a tiny mode AT the grid's resolution wall is the zero-count population's own kernels, never selected by the anchor rule. `--json` dumps census + curve + rug for rendering elsewhere (a renderer would duplicate no scorer). `--self-test` 13/13, perturbed, no I/O |
| `design/landscape_head_to_head.py` | ⭐⭐⭐ **HOW GOOD IS THE TOTAL-DENSITY LANDSCAPE, AND WHAT IS ITS BANDWIDTH REALLY?** ⛔ **Its NPMLE arm was REMOVED with the NPMLE itself (2026-08-21) — you cannot A/B a deleted thing — and the verdict was MOVED to `ROADMAP.md` §1 rank 3 + `DESIGN.md` §3.1a-iii, not lost.** What it measured, for the record: No solver, no EM, nothing patched in `src/` except a per-arm `_N_GRID` (restored in a `finally`, and that restoration is gated). ⭐⭐ **The structural finding comes first because it decides what a score can mean**: the NPMLE is fit on `(unspliced_count, eff_gdna)` — a total over ONE component's OPPORTUNITY MODEL — while the landscape is fit on the START/END banks over the region's own LENGTH, and since the origin partitions' START banks sum to the full payload's, **the landscape's inputs ARE the truth up to Poisson noise while the NPMLE's are the truth pushed through a model**. Three scorers, none needing a chosen constant: ⓐ the AXIS OFFSET `log(M/E_gdna) − log(count/exposure)`, whose SPREAD is irreducible (IQR **0.12 nats capture-OFF, 1.66 UNDER CAPTURE** — comparable to the mode separation itself); ⓑ **HELD-OUT PREDICTIVE LOG-LIKELIHOOD** on index-parity halves, both arms scored on the SAME held-out pairs (⛔ EMD is monotone in smoothing and may not select anything — a predictive likelihood is not, which is what makes the sweeps a SELECTION); ⓒ the depleted level against the `gdna` partition's pooled intergenic rate, ⛔ STAMPED VACUOUS at `g00` where the truth is exactly 0. ⭐ **Verdict: the landscape wins the estimand — 0.0056 / 0.0264 nats vs the NPMLE's 0.1191 / 0.1266 on the two capture-OFF strata (4.8–21×) — and ⓑ TIES off capture (−5.36 vs −5.35) while winning by 0.47 nats under it.** Read both: the tie says the NPMLE is not a bad density estimate, the nats column says it is the wrong quantity. ⛔⛔ **AND ITS BANDWIDTH FINDING RETIRES THE OBVIOUS KNOB: `_KNN_SCALE` IS VERY NEARLY INERT AT PANEL SCALE** — `knn_widths` floors every kernel at one grid step and **98.9 % of kernels sit at that floor** at the shipped 0.5 (92 % at scale 2.0), so the bandwidth in force is the GRID STEP, `span/(_N_GRID−1)` ≈ 0.034 dec. `--grid-sweep` prices THAT: the likelihood knee, the reproducibility maximum (0.976, FALLING to 0.889 at 1040) and the anchor-gap collapse (0.126 → 0.003 nats) **all land on the shipped 260**, which is the selection evidence that constant never had. ⛔ The mode COUNT is a resolution artefact (3.7 → 17.7 across the sweep, never converging, buying 0.006 nats) — **no consumer may depend on it**. ⛔⛔ **`span_R` IS GRID-FRAGILE and rung 4 was going to consume it**: 58 → 95.6 → **1.9** on `g50 ss0.99` capture-OFF, because `split_basins` picks the enriched mode by BASIN MASS and over-resolution fragments the bulk into competitors. `rho_0` and the anchor verdict are grid-robust; `R` is not. `--self-test` **42/42** (49 before the NPMLE arm went), six comparators perturbed |
| `design/transfer_variance_audit.py` | ⭐⭐⭐ **WHERE DOES THE MESSAGE TRANSFER VARIANCE COME FROM, AND COULD THE `AbundanceLandscape` SUPPLY A BETTER ONE? — an AUDIT: it proposes nothing and changes nothing.** ⭐ Gate ① verifies the NPMLE's retirement **from the AST, not from the comment that claims it** — and since 2026-08-21 it checks the strongest available fact, that the MODULE IS GONE and nothing under `src/` calls `.project()` (the formula survives in one docstring as provenance, which is correct — a grep would false-positive on it). ⚠ Its own falsification had to move with the deletion: it used to prove the scanner could still SEE a call by scanning `test_npmle.py`, which no longer exists, so it now proves that on a synthetic source AND proves a prose mention does not count (`TRAPS: an-ablation-that-never-ran`). ⭐⭐ **THE CENTRAL FINDING: the retired role is NOT vacant in the tool, it is vacant in the SHIPPED policy.** The replacement `transfer_logvar` = `Var(log rho_tot^dst) + Var(log rho_tot^src)` is a counting term plus a composition term that vanishes at `Var(f_g) = 0`, so ⛔ **every term shrinks as either slot is counted more deeply, where the retired between-mode term did not** — measured, **25.6 % of real hops carry the whole damping term below 0.1 nats² off capture**, i.e. a deeply-counted transport arrives essentially undamped. `currency.premise_logvar` fills that hole by a better route, so the gap is `RelayPolicy`'s alone. ⛔⛔ **AND THE ANSWER TO THE OWNER'S QUESTION IS NO, FAILING IN BOTH DIRECTIONS AT ONCE**: the landscape's per-slot grid posterior is **0.92× the counting variance capture-OFF** — TIGHTER, the wrong direction, because a fitted population prior shrinks a slot toward a mode and substituting it would make every message MORE confident; while its population spread **over-states a fitted premise by 7.4–8.0×** because it contains the real enrichment/expression structure the message is trying to transport (signal, not uncertainty). ⭐ What survives: under capture the posterior LOOSENS (1.03×), and that increment IS mode-membership ambiguity — the retired term reappearing legitimately, on the right substrate and in the right direction; it is the honest seed of the per-hop-type premise `premise_logvar`'s docstring already names as the obvious refinement. ⛔ The composition-vacuity limit is unmoved: a TOTAL licenses a LEVEL hop's variance and the population premise, never a COMPOSITION hop's. ⚠ σ² is read at pass-0 init where `Var(f_g) = 0` makes `composition_logvar` reduce EXACTLY to the counting term — asserted per condition, never assumed; adjacency and the opportunity pair are the relay's own. `--self-test` **20/20**, perturbed |
| `design/calibration_oracle.py` | ⭐⭐⭐ **THE CERTIFIED PER-OBJECT GROUND TRUTH — run this before debugging calibration against anything.** Every REGION and BOUNDARY's accumulator count + realized `n_gdna`/`n_nrna`/`n_mrna` + `true_f_g`, REFUSED unless its named gates pass (sum-to-full, partition-projects-exactly, gdna-field-uniformity, exact-zeros, nascent-in-annotation) — a merely *plausible* oracle is how a calibration bug gets debugged against a truth bug and both survive. ⭐ Two certification LEVELS: COMPOSITION (`true_f_g`, no opportunity model anywhere in it) and FIELD (densities too). ⭐ **Its first real run (2026-08-18) caught a ~2.5 % boundary-crossing deficit — CLOSED 2026-08-19**: every REGION class sits at ratio 1.000 against `rho·E` while every BOUNDARY class reads 2–3 % LOW, two refs, three independent sims. ⭐ The cause was a transcript-set chimera predicate silently dropping 4,087 gDNA fragments per condition (`DESIGN.md` §3.1b-ii; `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`); the companion "R exon" half was a SCORER artefact and is RETRACTED (`TRAPS: a-ratio-needs-a-population-that-can-supply-its-numerator`). The capture-OFF field gate went **8/16 → 16/16**. Writes `slot_truth.npz` beside each oracle cache; `--self-test` 11/11, falsified by perturbation |
| `design/total_abundance_audit.py` | ⭐⭐⭐ **IS THE MEASURED TOTAL A TRUE TOTAL? — the VALIDATION rung of the measured prior, on all 32 conditions, no solver and no re-scan.** Five arms against the origin partitions, and the one to read first is ⭐⭐ **ⓔ START/END AGREEMENT, the only FIELD-FREE arm and the decisive test of the WALL RULE**: `S/ell` and `E/ell` estimate the same rate with the same opportunity in the same field, so both-exact regions must AGREE and a bound side must read LOW — a DIRECTION, which a random mask cannot fake and a sides-swapped one inverts. ⭐ **Measured on the ladder, every stratum and both zero controls: both-exact 1.0018–1.0057, start-exact/end-BOUND 0.712–0.879 (<1), end-exact/start-BOUND 1.242–1.428 (>1)** — 16/16, capture-ON included. ⭐ ⓐ the TRUNCATION LAW (`Sum Cinv / Sum (S/ell)·P(w<=ell)`, 1.000 is the claim) confirms the START bank is an untruncated TOTAL where the contained bank is a SHAPE: **1.00025 / 1.00045 at capture-OFF over 14.4 M fragments**, and the gDNA rows are the DECISIVE ones because gDNA cannot splice so the law is exactly formable there. ⛔⛔ **IT ASSUMES LOCAL FIELD UNIFORMITY AND CAPTURE VIOLATES THAT — a FINDING, not a failure**: the two banks weight a region's positions differently, so under a probe-shaped field they are different weighted averages of it (up to 30× on intron/intergenic), which is why a consumer swap under capture changes the QUANTITY rather than removing a bias. ⭐ ⓑ the absolute arm is **gDNA ONLY** — RNA density varies per transcript BY DESIGN, so a pooled RNA rate scores the panel's expression profile and not the estimator — and it is STAMPED VACUOUS under capture. ⭐ ⓒ re-derives the mask's mature distances with a per-transcript Python loop sharing no helper with the shipped kernel and demands EXACT agreement; **its first run caught a 57-distance divergence — the two implementations disagreed on SYNTHETIC spans**, and the population rule is now stated in both rather than inherited from a table's contents. ⛔ Gates refuse rather than warn: the ledger twice over on all four payloads, both new banks summing to the whole (max\|delta\| = 0 on 32/32), and the mask's two implementations agreeing — **194 gates, 0 failures across both panels**. `--self-test` 15/15, every gate perturbed |
| `design/calibration_walk.py` | ⭐⭐⭐ **WHICH STAGE OF CALIBRATION INTRODUCES THE ERROR?** The solve as a ladder — init → strand → message-free local (refit 0) → +messages → +landscape refits → SHIPPED — each rung scored per stratum against `calibration_oracle.py`'s certified table, which it REFUSES to run without. Every arm asserts what ran (`_uni` only under `RelayPolicy`; a muted arm must reproduce `fg_loc` bit for bit). ⭐ Its first run isolated the zero-gDNA nascent catastrophe to the LOCAL SOLVE at introns (½ → 0.80 with truth 0) and exonerated messages and refits — one function to open instead of a tool-wide hunt |
| `design/relay_pool_ab.py` | ⭐⭐ **MESSAGE PROPAGATION OFF vs ON, per condition, per pool, NEVER collapsed** — gDNA/RNA estimates against origin-split truth with `net` (signed) and `abs` (misplaced-mass) errors in FRAGMENTS, three truth pools printed so a gDNA over-call names which RNA it ate. ⛔ Rows whose nascent truth is exactly 0 are STAMPED not-a-measurement (`nrna = 0` restated is not a nascent verdict). Both arms in ONE process off one cached payload; the relay arm asserts it ran. `--self-test` 11/11, falsified by perturbation |
| `design/benchmark_report.py` | ⭐⭐ **THE BENCHMARK AS ONE HTML PAGE — every scenario, in COUNTS, never pooled until the end.** Owner, 2026-08-19: *"report the results across every scenario … DON'T POOL SCENARIOS … ONLY POOL AT THE END"*, and *"a script to produce an HTML with the results"*. One row per scenario per arm with `gdna_true`/`gdna_est`/`net`/`abs` and the RNA truth split beside them, the pooled total LAST and marked as a summary of the rows above it, plus a per-TRANSCRIPT observed-fragment table read from each scenario's own `truth_abundances.tsv`. ⛔ **It SCORES NOTHING** — it renders `relay_pool_ab.py --out`, because duplicating a scorer is how a baseline and a ceiling drift apart. `--self-test` 10/10 (the pooled row must come last and equal the rows above it; a transcript missing from a scenario renders 0, not blank) |
| `design/hop_currency.py` | ⭐⭐⭐ **STAGE 2 OF THE MESSAGE-PROPAGATION REBUILD — WHICH CURRENCY EACH HOP TYPE CARRIES, MEASURED ON ALL 32.** Every ordered adjacent pair `dst <- src` classified, the source's TRUE value transported BOTH ways (LEVEL `rho_g`, COMPOSITION `f_g` of the population that ENTERS the destination — direction-dependent, the SPLICE-IN population includes the spliced fragments and the sj flux on the exon's side) and scored at the destination in FRAGMENTS beside a Monte-Carlo NOISE FLOOR, plus the depth census. ⛔⛔ **ITS FIRST FINDING IS THAT THE SEVEN OBJECT CLASSES DO NOT NAME THE HOP TYPES**: a `B exon\|intron` is a SPLICE SITE or a TSS/TES inside another transcript's intron, and the whole of the composition error into exons sat on the second (209 k of 221 k on one condition); the hop-type key is `object class x {sj, term, sj+term}` from the boundary flags. ⭐ The map: LEVEL out of an exon and into one from a terminus; COMPOSITION into an exon from a splice site and from an intron into its own boundary; ⛔ NEITHER at a terminus into an exon UNDER CAPTURE (0.57–0.65 M fragments per hop type at `g50`, the whole residual census — a third mechanism). ⛔ The floor must reproduce the scorer's selection (zero-truncated at ~1 crossing per boundary) or `g05` reads a 10 % rule error that is not there. ⚠ Carries a PANEL caveat: the simulator enriches nascent RNA only by its own transcript's probes, so `nrna_mid x capture-ON` intron→boundary rows read 5–37 % that is not a rule failure. `--self-test` 36/36, falsified by perturbation; ~18 s for the ladder |
| `design/exon_solvability.py` | ⭐⭐⭐ **WHICH EXONS ARE SOLVED ACCURATELY, AND WHAT — WITHOUT TRUTH — TELLS YOU WHICH? The TRAINING-SUBSTRATE question for the gDNA landscape prior.** Calibration is circular: a per-object prior can only be trained on the exons pass-0 already solves, so the deliverable is WHICH OBSERVABLE CLASS is safe to train on. ⛔ **NO THRESHOLD ANYWHERE** — "solved well" is the mass-weighted CURVE (fraction of exon mass with `|f_g − true_f_g| ≤ x`) plus mass-weighted quantiles plus `Σ|f_g − true|·M` in FRAGMENTS; nothing branches on an error value. Five truth-free predictors: own composition evidence (the solver's OWN `has_own_composition_evidence`, imported not restated), AMBIG vs single-strand, the exon's own depth by mass-weighted QUARTILE (data-derived edges, equal-mass buckets), hops to the nearest structurally pure-gDNA slot (reuses `hop_currency.depth_to_measured`), and the flank's sj/terminus structure. ⭐⭐ **THE ANSWER IS THE STRAND BIT** — pooled per stratum on the ladder under the shipped relay, single-strand exons read **8.9 / 16.1 / 15.8 per 1k** misplaced (stranded×OFF, stranded×ON, unstranded×OFF) against AMBIG's **37.5 / 63.2 / 39.4**, and AMBIG is 18–20 % of exon mass carrying ~50 % of exon misplacement. ⭐⭐ On every stranded contaminated condition `has_own_composition_evidence` at an exon is the SAME MASK as single-strand, to the slot (the strand λ-term is credited to single-strand slots only) — reported, never asserted. ⭐ The second predictor is the exon's own DEPTH, the only one monotone on all three in-scope strata (**4.6/7.7/14.5/29.7** per 1k, deepest→shallowest quartile, stranded×OFF). ⛔⛔ **THE FLANK STRUCTURE REVERSES BETWEEN STRATA** — `sj only` is the best class at the zero control (3.9/1k) and the worst on stranded×capture-ON (49.6/1k) — so it is not a training criterion. ⛔ Names a class no other predictor does: the **47 ERCC references with a live exon have no pure-gDNA slot at all** (39,354 fragments, `unreachable`/`reference end`), unreachable by any message. ⛔⛔ **Its verdict on the Stage-3 CurrencyPolicy at exons: a regression on all three IN-SCOPE strata (1.08–2.21×) and a win only on the DEFERRED one (0.80×)** — it MOVES accuracy from evidence-bearing exons to evidence-free ones (AMBIG 39,401 → 24,532 while single-strand 41,321 → 74,668). ⚠ The test chromosome says the OPPOSITE (a 4× currency win at `g50` stranded×OFF) — `TRAPS: panel-before-src`. `--refits 0` is the PASS-0 reading, `--self-test` 39/39 falsified by perturbation |
| `design/solvability_audit.py` | ⭐⭐⭐ **START HERE FOR PASS-0 — and 0.8.0 is judged here, not on the transcript table.** Which objects are SOLVABLE, which of those are solved wrong, and — the one that matters — which are **confidently** wrong, since a tight-variance wrong value outvotes correct neighbours and anchors the prior. ⛔ Excludes the undetermined population: in pass-0 an object with no own evidence reporting `f_g ≈ ½` at zero precision is *correct*. Carries the ablation ladder (strand → local → final) and a calibration curve that needs no threshold. `--suite` alone gives the whole panel |
| `design/prior_vs_oracle.py` | ⭐⭐⭐ **CALIBRATION'S ENDPOINT AGAINST TRUTH** — `LocusPriors` (`gdna_prior_count`, `rna_prior_count`, `gdna_eff_len`) is what the EM actually reads, and until 2026-08-07 nothing compared it to anything. Arms: **P** shipped, **O** the same assembler fed the origin-split truth masses, **S** plus each component's own true share, ⭐ **Fo** the EM's OWN candidate count per origin, **F** the first-base count. `P−O` is calibration's error and `O−Fo` the ASSEMBLER's — different repairs in different files. ⛔⛔ **`F` IS NOT THE EM's TARGET and every `O−F` printed before 2026-08-08 was scored against a quantity nothing consumes** (`TRAPS: score-the-consumers-own-count`); it is retained and priced on table ⑪. ⭐ Reports the count, the composition claim `phi = a_g/(a_g+a_r)` and the SCALE separately, per stratum. ⛔ Undrained on every arm, priced not waved |
| `design/pass0_vs_oracle.py` | T (the origin-split payload) vs P (`calib_refit_iters=0`) vs two levered ceilings, per object and per class. ⛔⛔ **ITS HEADLINE MASS-WEIGHTED ERROR IS THE WRONG YARDSTICK FOR PASS-0** — it scores every object with mass, so honest ignorance reads as error and buried a 0.0456 answer inside a 0.3150 one. Use it for the T/C/P decomposition and the oracle plumbing; use `solvability_audit.py` to judge pass-0. ⛔ undrained on every arm. `--oracle-cache` makes repeat runs cheap: the oracle depends on the accumulator and index, never on calibration, so one cache serves a whole debugging campaign |
| `design/worst_objects.py` | ⭐⭐ **step 3 of the debug loop** — one condition dissected to individual regions/boundaries, ranked by error **MASS**. Read the CONCENTRATION curve first (concentrated ⇒ a mechanism exists; diffuse ⇒ a systematic bias). `fg_loc` vs `pred_fg` separates a wrong local solve from wrong messages |
| `design/calibration_truth_ab.py` | ⭐⭐ the deliverable against truth, and `--ceiling` — what perfecting each length model is *worth*. ⚠ "length model" here is the fl PMF inside the OPPORTUNITY model, not the composition channel deferred past 0.8.0 — the two share a word and nothing else. ⛔ Read the ceiling caution in Working rules before quoting a number: the effective-length shrinkage sits OUTSIDE every arm's patch point |
| **⭐⭐⭐ the panel, and the caches that make calibration a seconds-long loop** | |
| `sim/configs/flgap_rna_long.yaml` · `flgap_rna_short.yaml` | ⭐⭐ **THE fl-GAP SIDE PANEL — TWO ARMS, OPPOSITE SIGNS, and it exists to exert the ONE mechanism the ladder nulls by design.** gDNA 75 / RNA 250 and the reverse, each × ss `0.50/0.99` × capture off/on at `g50` — 4 conditions per arm. ⛔⛔ **A SIDE panel, NEVER a ladder rung**: the ladder equalises fragment lengths because the EM already reads the fl distribution, so a gap lets it split the origins on LENGTH ALONE and mask calibration bugs. ⛔ Therefore the TRANSCRIPT-level number here is not a calibration result; ⭐ but everything that stops before the EM is valid, which includes `total_abundance_audit.py` (banks only) and `calibration_vs_oracle.py` (it runs `calibrate`, never the EM). ⭐ **Why it exists**: the shipped `count / E_contained` and the pmf-free START/END pair differ by exactly `rho_r·(E_r/E_g − 1)` — the RNA share of a MIXED pool times a factor the ladder pins at 1 (+0.02…+1.6 % at its −12.4 bp gap) and a real library does not (−0.2…−78 % at +87.6 bp). ⚠ **The std MOVES with the mean and that is forced, not chosen**: the sampler truncates to `[frag_min, frag_max]`, so cfg (75, 98) REALISES at 137.65 bp — measure the realised gap off the partitions and quote that, never the configured number (the audit prints it on every row). ⛔ Both arms are required: RNA is NOT reliably longer than gDNA, and `E_r/E_g − 1` FLIPS SIGN between them, so a real repair must move the two in opposite directions |
| `sim/configs/gdna_ladder.yaml` | ⭐ **the Stage-B panel, and the only panel the tool is RANKED on** (⚠ it was "the ONLY panel on disk" until the fl-gap side panel landed 2026-08-21) — **16 conditions** (`g00/g05/g50/g98` × ss `0.50/0.99` × capture off/on), ⭐ **every one of them carrying NASCENT RNA at a 20 % fragment share** (owner, 2026-08-19: the nascent ON/OFF contrast is a TOY experiment, not half a panel — *"we don't need to spend 16/32 scenarios just to add nascent RNA"*), so no row of this panel is `nrna = 0` restated, EQUAL fragment lengths, gDNA 0 → 98 % at a fixed 10 M total. ⛔⛔ **WHY EQUAL LENGTHS — and the reason is stronger than the one that used to be written down: THE EM ALREADY USES THE FL DISTRIBUTION.** A large gDNA-vs-RNA length gap lets the EM assign fragments on LENGTH ALONE, bypassing calibration entirely and MASKING its bugs; equal lengths FORCE calibration to be exercised, leaving it strand and density (owner, 2026-08-14). ⛔ **REBUILT FROM SCRATCH 2026-08-13**, when `pilot`, `flgap_short` and `flgap_long` were deleted and the ladder went 36 → 16; `g05` survives because `suite_resolves.py`'s requirement (f) tests the RATE and demands one rung in `0 < rate ≤ 0.10`. ⚠ The ORACLE CACHE beside it at `ladder/oracle_cache/` holds the four origin partitions for all 16, but its `_main` full payload was written for only **12 of 16** BY DESIGN — `pass0_vs_oracle.py` holds every zero-gDNA condition out, so `panel.py status` prints ✘ on a complete panel. ⛔ **That "12" is no longer stable and a reader must not treat a `_main` as evidence of anything**: `pass0_vs_oracle.measure_condition` WRITES `<oracle_cache>/<condition>/_main` whenever `--oracle-cache` is passed, so any arm run over the `g00` rows creates the missing four (`TRAPS: shard-an-arm-sweep-by-condition`). `docs/TESTING.md` §0 |
| `sim/panel.py` | ⭐⭐⭐ **THE SIMULATION + BENCHMARKING WORKFLOW, ONE COMMAND PER STAGE** — `status` / `build` / `simulate` / `cache` / `score` / `report`, every path derived from one panel YAML. ⛔ It adds NO measurement code: each stage shells out to the instrument that already owns it, because duplicating a scorer is how a baseline and a ceiling drift apart. ⭐⭐ **`status` FIRST** — every stage is expensive and resumable, and it names the next one. ⛔⛔ **`cache` builds BOTH caches**, and the oracle one is why this exists: it is the origin-split truth every scoring instrument reads, and until 2026-08-11 it was a SIDE EFFECT of `pass0_vs_oracle.py` with no step of its own, so the documented recipe produced a panel every scorer refused. ⚠ A cached oracle condition needs all FOUR parts (`gdna`/`mrna`/`nrna`/`_main`). Gated by `tests/test_panel_workflow.py`, falsified by perturbation |
| `design/build_scan_cache.py` | scan once, calibrate many times. ⛔⛔ **THE CACHE KEY COVERS `accumulator.cpp`'s DEPOSIT AND NOT `resolve.cpp`'s FRAGMENT CONSTRUCTION** (rule restated 2026-08-20 — the old *"the schema digest does not see a deposit-rule change"* is stale). `scan_cache.py:220-224` hashes `deposit_digest()` into `payload_schema_digest()`: the NATIVE accumulator runs over a fixed 7-fragment fixture and every bank's raw bytes are hashed, gated by `tests/test_scan_cache.py`. ⚠ What the key still does NOT see is a change to which fragments are OFFERED — the 2026-08-19 `cached / skip` on all 32 was the chimera repair, a `resolve.cpp` change, and the digest was right not to move while the caches were wrong to be reused. For THAT class use `--force` (wired through `panel.py`) or delete the caches by hand. ⛔ `rescan_panels.py` is NOT the tool for it: it gates shared banks on BYTE-IDENTITY, right for a rename or a bank ADDITION, and stops on every condition when the counts legitimately move |
| `design/rescan_panels.py` | ⭐⭐⭐ **RE-SCANNING IS THE FIRST IRREVERSIBLE STEP, SO IT CARRIES ITS OWN FALSIFICATION.** After a deposit-rule or schema change every cache is refused and must be rebuilt — and the old ones are then gone. This rebuilds them and GATES each condition on byte-identity: every bank the change was not supposed to touch must come back identical, and the symmetric difference must be the SAME delta on every condition. ⛔ Reads the stale `.npz` BEFORE the write, because `write_scan_cache` overwrites in place and would otherwise destroy its own baseline; a condition that fails is NOT written. ⛔ No expected delta is hard-coded — it is derived from the first condition and required of the rest, so a side effect appearing only under capture or only at `g98` is caught. ⚠ Counts an ALL-ZERO baseline as VACUOUS rather than a pass (`g00`'s gdna, every `nrna_none`'s nrna). `--self-test` perturbs the comparator with no I/O — 14/14; `--jobs N` shards conditions. ⛔⛔ **"Byte-identity" is now EXACT FOR THE INTEGER BANKS AND A DERIVED BUDGET FOR THE SIX FLOAT ONES, and the difference is not a relaxation — it is a repair (2026-08-13).** The comparator demanded exactness because *"these are integer tallies"*; the 2026-08-10 one-numeric-convention ruling made six of them float64, float addition is not associative across worker threads, and the gate has been **UNSATISFIABLE ever since** — two scans of the same BAM by the same binary differ on exactly those six, by up to 3.5e-14 relative. ⭐ The budget is `deposits · eps · value` per element, with `deposits` read from the COUNT bank at the same object and the ×2 on the mass banks derived from the deposit rule — no tuned constant. ⚠ It cost most of a session, as `TRAPS: a-stale-gate-accuses-the-newest-change`: the failure named five banks the change had not touched, and the reading "my change had a side effect" was wrong |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | build the panel. ⚠ `panel.py build` drives the last two; the reference carve needs the SOURCE genome/GTF, which a panel config does not name, so it stays manual |
| **⭐⭐ the prior assembler, and the end-to-end thermometer above it** | |
| `design/quant_accuracy.py` | ⭐⭐⭐ **THE TOOL END TO END, AND THE PRIOR-INJECTION CEILING ABOVE IT** — `--arm base` is the transcript-level accuracy table, `--arm oracle` re-quantifies with the perfect prior injected, and the three `oracle_<field>` arms say WHICH of the prior's three numbers carries the value. ⛔ Scores **count against count** against the condition's OWN `truth_abundances.tsv` (the realised observed fragment count); the suite-level table is the *pre-capture molar* abundance and is a different quantity. ⛔ Pins `EMConfig.seed` — the shipped default is `None` (`TRAPS: the-deliverable-is-not-reproducible-by-default`) — and `base_reseed` prints the noise floor beside the effect. ⚠ **0.8.0 is judged on CALIBRATION against oracle calibration, so this table is the THERMOMETER above it**: optimise against `solvability_audit.py` and read this to see whether it moved |
| `design/mass_prior_ab.py` | ⭐⭐⭐ **THE PRIOR AS A CONSERVED FRAGMENT COUNT — plan 5.5 against the oracle arms.** `edge_unspliced_mass` sums to ONE per fragment, so the assembler can stop manufacturing a count from a density. Subsamples a real condition by qname hash (so the whole and all three origin partitions stay consistent), then re-assembles `O` and scores it against `F`. ⭐ **Verdict: the capture-ON over-call goes `O/F` 1.149 → 1.019 on three gDNA levels; capture-OFF `Σ|Δ|` falls 61 %.** ⛔ Prints the SUBSTRATE check first — the subsample must reproduce the defect before its disappearance means anything — plus a `share≡1` ablation and the plan-8-Q4 ceiling. ⚠ It MEASURES the shipped bank and never re-derives it; the rule's own gates are `tests/native/test_conserved_mass.py` |
| `design/transcript_truth.py` | ⭐⭐⭐ **THE TRUE PER-TRANSCRIPT COUNT, SPLIT BY SPLICEDNESS** — what a per-transcript prior is scored against, and the one quantity no other instrument produced (`prior_vs_oracle.py` is per multi-locus; `quant_accuracy.py` scores the tool's OUTPUT). One pass over the oracle BAM, **read names only, 41 s per 10 M-fragment condition**. ⛔⛔ **Splicedness comes from the read name's interval in SPLICED TRANSCRIPT coordinates against the transcript's interior cumulative-exon boundaries — NEVER from the CIGAR**, which misses every sj falling in the unsequenced inner gap: measured **1,233,428 of 10 M (12.3 %)** truly-spliced fragments carry no `N`, and **0** in the reverse direction. ⭐ Seven NAMED gates; the free falsification is that a CIGAR `N` where truth says unspliced is structurally impossible, and it reads **0**. ⚠ `t_df` encodes strand as `{1: POS, 2: NEG}`, not `±1` — reading it as a sign silently disables the minus-strand exon reversal. ⛔ The `indexed` gate currently FAILS by design: **329 fragments from 2 transcripts the index does not contain** (`ENST00000397293.6`, `ENST00000432560.6`), reported rather than absorbed as "unspliced" |
| `design/transcript_weights.py` | ⛔⛔ **THE SOFT-MIN PER-TRANSCRIPT WEIGHT — BUILT, PRICED AND REFUSED (2026-08-13). Read `ROADMAP.md` §4's row before touching it; this is a RECORD, not a proposal.** The owner's theorem — a transcript's density is bounded by the MINIMUM along its path — implemented as an opportunity-weighted power mean `M_p`, one dial from the pooled `Σmass/Σopportunity` (`p=1`, the control, no deconvolution at all) through the owner's harmonic (`p=-1`) to the theorem verbatim (`p→-∞`), times a transcript opportunity (`total` / `full` / `seen`). ⛔ **16 arms, both strata: worse than `base` at transcript level EVERYWHERE (1.26–1.60×), and the `g00` gene-level win of 0.395–0.527× collapses to 1.006–1.041× on the blind stratum.** ⭐ **What survives and is worth reusing:** the gated density identity `density(r) = Σ_{t∋r} A_t/E_t` (the object's own opportunity cancels — this is what makes the `min` a min of something); the measured multiplier ruling (TOTAL abundance beats the true unspliced count, 0.163× vs 0.202×); the opportunity weighting, which is a **deliberate deviation from the spec as written** because `TRAPS: a-mean-of-ratios-inherits-the-partition` was measured *while designing this prior*; and the monotone dial (`min` < `harmonic` < `geometric` < `arithmetic` on every axis), which says the theorem is not what failed. ⛔ What failed is `TRAPS: an-upper-bound-is-not-an-estimate`. ⚠ Standalone it runs ONLY its own falsification — the re-partition gate, where unweighted forms drift 1.67–3.03× on data that did not change and weighted ones do not move. Scoring belongs to `quant_accuracy.py`'s `alloc_*`/`allocg_*` arms; a second scorer here is how a baseline and a ceiling drift apart |
| `rigel.sim.net_flow` (a MODULE, not a script) | ⭐⭐ **WHERE DID EACH MISASSIGNED FRAGMENT GO?** — the one question `quant_accuracy.py` does not answer. It measures the MAGNITUDE of transcript error; this decomposes its DIRECTION, per transcript, into gDNA-sourced and RNA-isoform-sourced flow. ⛔ **All that survived `sim/analysis.py`, retired 2026-08-11** (owner): that was a 1,589-line SECOND SCORER rendering its own accuracy tables against its own truth, and two scorers is how a baseline and a ceiling drift apart. 618 lines kept with their tests, 971 deleted along with `scripts/sim/evaluate_suite.py`. Gate: `tests/test_net_flow.py` |
| **⭐⭐⭐ where to develop** | |
| `design/rename_identity.py` | ⭐⭐⭐ **EVERY RENAME STAGE PROVEN BIT-IDENTICAL — the gate the bulk rename lands on.** A rename is a numeric NO-OP and this makes the claim falsifiable: `--freeze` captures a reference ONCE, `--check` compares after every stage. ⭐⭐ **It can exist only because `total_threads=1` makes the scan bit-reproducible** — measured, the shipped all-cores default differs run-to-run on SIX float banks by ~3.5e-14, so the harness pins the thread count, the BGZF threads, `OMP_NUM_THREADS` and the EM seed. ⛔⛔ **It compares CONTENT, never NAMES, because the names are the thing under test**: ① the sorted `(dtype, shape, sha256)` multiset of every array, names discarded, and ② the transcript table, which carries no region/boundary vocabulary at all. Neither alone suffices — ① is blind to two arrays SWAPPING names and ② is not — and the hole is pinned as a hole in `--self-test` (8/8) so it is not mistaken for coverage. ⛔ No name map: that is the compatibility hack the owner refused. ⚠ The reference is FROZEN, never rolling — `--freeze` refuses to overwrite — because the renames COMPOUND and a rolling baseline lets a stage-2 defect become the accepted truth for stages 3-9; this is the one place `TRAPS: re-record-the-baseline` is deliberately inverted. ⭐ Falsified on real data: a ONE-ULP nudge to one element of a 13,482-element array fires it, and a restore reads clean |
| `design/rename_census.py` | ⭐⭐⭐ **THE BULK RENAME'S LANDSCAPE, RE-DERIVED — every name the vocabulary ruling touches, by KIND.** REGION (was `node`) · BOUNDARY (a single position BETWEEN TWO REGIONS — was `edge`/`cut`/`line`/`seam`) · SPLICE JUNCTION or `sj`, never bare `junction`. Measures ~10,900 occurrences: 5,192 Python identifiers over 845 distinct names, 978 in C++, and ⭐ **4,733 in PROSE — 43 % of the job**, so a code-only rename leaves half the vocabulary stale. ⛔ **It REPORTS and never renames**, because the danger is not the volume but the names carrying TWO senses, where a correct-looking global replace silently corrupts the other one — four found so far: `cut` (35,229 positions vs 35,041 interior boundaries — three populations, and mapping it wholesale would recreate `TRAPS: two-masks-one-name`), `line` (120 boundaries vs **60 ordinary lines of text**), `donor` (a splice donor vs the toy harness's SOURCE CONDITION), and `node` (vs `ast` nodes, in this very file). ⛔⛔ **THE HISTORICAL TOKENS ABOVE ARE BACKTICKED BECAUSE THIS ROW WAS ITSELF EATEN BY THE RENAME IT DOCUMENTS** — it read *"REGION (was region)"* and *"60 ordinary boundaries of text"* until 2026-08-17, restored from `391aa5eb^`. That is `TRAPS: a-rename-resolves-an-ambiguity-silently` committed against the index entry of the instrument built to prevent it. ⭐ `--sense <token>` dumps every site of one token with context, for a per-site ruling. Plan: `docs/dev/rename_plan.md` |
| `design/module_census.py` | ⭐⭐⭐ **WHERE DOES A CHANGE GO?** — the calibration package's real shape, re-derived from the AST: the layering with every upward import, the graph with each module's importers inside and outside, ⭐ **docstrings that name a sibling with no import boundary** (14 measured, 6 genuinely stale), and dead public surface. ⛔ It REPORTS, it does not judge — an entry point looks dead, and a gated reference implementation looks duplicated |
| **⭐⭐⭐ the backbone** | |
| `design/arm_identity.py` | ⭐⭐⭐ **IS THIS ARM BYTE-IDENTICAL TO THAT ONE? — the A5 gate.** `arm_score.py` AGGREGATES, so a difference that cancels between two fields is invisible there; this compares **every scored field of every row** and exits nonzero on any difference. ⛔ Both of A5's recorded lies are gated: the row-key sets must be EQUAL (an arm with ZERO rows once scored "32/32 IDENTICAL") and both files' mtimes are printed (a stale baseline is what A5 warns about — re-record it, B8). ⭐ Falsified by perturbation: it resolves a **1-ULP** nudge and refuses an empty file |
| `design/backbone_parity.py` | ⭐⭐⭐ **WHAT DOES ONE MESSAGE OPERATOR DO, PER SLOT?** — two policies, one real 70,176-slot chain, one process, every output array compared element by element, plus every diagnostic `_capture` key and ⭐ the five backbone assertions with their ELIGIBLE sets (A14). Strictly stronger per condition than the panel and strictly weaker across conditions, so run it FIRST. `--arm-a head --arm-b no_<switch>`. ⛔ Its verdict on the restructure: **421,056 output elements and 18,245,830 diagnostic elements, zero differences** |
| **⭐⭐ the toy harness** | |
| `design/toy_panel.py` | ⭐⭐ **one toy spec × EVERY cached condition × an RNA-density ladder, scored per object** — what you run once a structure is the *target* rather than a probe. Prior-free pass-0 by default; the RNA density is a multiple of each donor's OWN gDNA density. ⛔ Reports which object carries the error, whether the messages helped or hurt it, and which objects are CONFIDENTLY wrong. ⚠ 13 s per condition; shard with `--conditions` |
| `design/verify_toy_substrate.py` | ⭐⭐⭐ **IS THE INPUT CORRECT? — no solver runs.** Every accumulator bank re-derived from the simulator's per-fragment TRUTH by an independent implementation, plus the splice combinatorics and the length marginal. ⛔ Run this on any new toy spec BEFORE reading a solver number off it |
| `design/verify_capture.py` | ⭐⭐ **hybrid capture, probes ON vs OFF on identical geometry** — the gDNA landscape, the length selection and the sj depletion, each gated on the direction the knobs predict |
| `design/toy_trace_error.py` | ⭐⭐⭐ **the full error trace for ONE scenario** — every region and boundary ranked by error MASS, the four init sources, the relay hop by hop with the licence state, then a ψ channel ABLATION (fidelity-checked to 1e-6) and ⭐ a RELAY-level arm that withholds the composition licence in both twins at once. ⛔ This is what "dissect the error" means |
| `design/zero_controls.py` | ⭐⭐⭐ **THE TWO ZERO CONTROLS — and the owner requires them on every experiment.** ZERO RNA (silent genes, truth `f_g = 1.000` everywhere) and ZERO gDNA (the `g00` donor, truth `0.000`). The truth is a CONSTANT, so every deviation is a false positive with nothing to cancel it. ⭐ The zero-RNA arm is the **biologically dominant** case: most annotated transcripts are OFF. ⛔ Prints the three-rung ladder (strand → self-solve → final), what the MESSAGES delivered as an implied `f_g`, and flags any EMPTY object, because a zero arm can be degenerate and then it tests nothing (`TRAPS.md` A14) |
| `design/reframe_walk.py` | ⭐⭐⭐ **THE REFRAME, WALKED — every count, both flank totals, and every hop in BOTH directions, on one two-exon transcript.** Not a summary statistic: per hop it prints `r` as used, `r` as the predecessor's single total would have given it, what TRUTH says the same ratio is, and the true gDNA-density ratio beside it — so the include/exclude decision at each `(BOUNDARY, side)` is a number. ⭐ Three gDNA arms (zero / high / very high), capture-OFF × unstranded, HIGH RNA so nothing is sparse. ⛔ The geometry is rebuilt and GATED against the solver's published frames to 1e-12 |
| `design/toy_dissect.py` | ⭐⭐ **one scenario opened to every slot and every channel** — the solver's own ladder (strand-only → self-solve → final) beside `cm_g` / `c_tau` / `cg`, so a dead channel is visible as a zero rather than inferred. ⛔ Step 3 of the debug loop, for a toy |
| `design/certified_rna_audit.py` | ⭐⭐⭐ **IS THE CERTIFIED-RNA CHANNEL WIRED?** — `edge_spliced` is a spliced fragment, so it cannot be gDNA and needs no deconvolution. Audits (a) is the bank populated, (b) does it have a divisor, (c) **is a precision emitted** — any one failing looks identical. Runs a TA × TB abundance grid on `tes_readthrough`. ⛔ Scores the **mass**, never `f_g`: the solver's `f_g` is the UNSPLICED fraction and the oracle's is spliced-inclusive, and confusing them reads a 1.5 % defect as a half-unit error. Its docstring carries the 24/24 verdict |
| `design/certified_q_census.py` | ⭐⭐ **CAN A CERTIFIED COUNT SPEAK ABOUT THE UNSPLICED SPLIT? — the answer is NO, and this is why.** Measures the splice-visibility `q = S/(S+C_R)` directly off the origin-split oracle on every condition (36 when measured, 16 since the 2026-08-13 rebuild); no solver runs, so the whole ladder is seconds. ⛔ Its verdict is a NEGATIVE and the docstring says so first: `q`'s mass-weighted median is 0.19–0.71, so the term §2d drops is the same size as the one it keeps, and the raw-count λ term is **worse than the uninformative reference on 12 of 36 conditions** (worst +0.4578). ⭐ Its two extreme rungs ARE the two zero controls, which is what makes the diagnosis certain |
| `design/vertex_ceiling.py` | ⭐⭐ **the re-solve ceiling on the REAL ladder** — pins oracle truth at a chosen object class in `build_node_init` and re-solves, with a `noop` arm that MUST be byte-identical (`TRAPS: a-cancelling-defect-pair`'s sibling discipline). ⛔ Its own verdict is a NEGATIVE and the docstring says so first: the 24.4 % it measures is the value of MISSING INFORMATION, not headroom. ⭐ Reusable override plumbing for any `build_node_init` prototype. ⛔⛔ **IT WAS DEAD AND THE SUITE WAS GREEN — repaired 2026-08-15.** It wrapped `CAL.region_sweep`, a name that vanished when the sweep was renamed to `solve_chain`, so it raised on its FIRST arm; and its `_rna_arm` monkeypatch took one argument after the arms were made symmetric. `TRAPS: a-green-suite-hid-five-dead-instruments`. ⭐⭐ **`--arm ref=A,B` now drives ψ's TWO reference exponents independently** — they are the Beta(a,b) pseudo-counts (`a·log f_g + b·log(1−f_g)` on the λ grid IS Beta(a,b)), so this is the harness for any composition-reference work; `ref=C` keeps the historical single-knob behaviour |
| `design/toy_ceiling.py` | ⭐⭐ **the RE-SOLVE ceiling** — hand one object class a different own belief and re-solve the whole chain, six arms sharing one simulation. ⛔ This is what `TRAPS.md` B17 demands instead of a substitution, and `--arms base noop` is its own falsification (must be byte-identical). Its docstring carries what it measured |
| `design/psi_channel_ablation.py` | ⭐⭐⭐ **WHICH ψ CHANNEL IS DOING THE WORK?** — one level below `worst_objects.py`: that says which OBJECTS carry the error, this says which CHANNEL put it there. Records every argument of the final ψ combine, then re-solves with one `*_imp` channel nulled at a time, so an attribution is a re-solve of the REAL call. ⛔ A5 — the `as-is` arm must reproduce the run bit-identically, and that means reproducing the WRITE-BACK too (the first version read `max\|Δ\| = 1.0` for exactly that reason). ⭐ Read the per-slot table, not the totals: D4h says all-small-singly + large-jointly means the channels share an upstream quantity |
| `design/arm_score.py` · `design/arm_sweep.py` | ⭐⭐ score two `ladder_arm_ab` arms against each other, and sweep a family, **per stratum**. ⛔ Never pooled — the panel total hides a sign flip between strata, and both print `abs_err_all_final` (the deliverable) beside `abs_err_all` (pass-0) on every row, per `TRAPS.md` A15. `$RIGEL_ARMS` points at the `--out` directory |
| `design/ladder_arm_ab.py` | ⭐⭐ **the same override on the REAL ladder** (36 conditions when these numbers were measured; **16** since the 2026-08-13 rebuild), scored by `solvability_audit.py`. ⛔⛔ **Run this before writing a mechanism into `src/`** — FOUR times now a toy- or single-condition-positive change has been panel-negative (`TRAPS.md` B18). ⭐⭐ **`--jobs 8` runs the whole 36-condition panel in ~2.2 min** (was ~20 min): the two unused `c_input_*` arms are no longer built and the MAIN scan payload is cached beside the oracle cache. ⭐ Carries the eight `zc_*` arms plus ⭐⭐ `msgfree_p0` / `msgfree_all` / `msgscale_<k>` / `onesided_rna`, which between them price the whole message layer. ⛔⛔ **`--messages {off,on}` IS PART OF THE ARM AND IS STAMPED INTO EVERY ROW** — under the shipped `off`, **22 of the 26 arms cannot move a number** and 16 of those score byte-identical to `base` while their fire counter reads healthy (`TRAPS: an-ablation-that-never-ran`, third form). An arm the policy cannot express is now REFUSED up front. ⭐⭐ **`--self-test` is the harness's own falsification**: every arm under both policies, gated INERT-under-one / MOVED-under-the-other — 25/25 in ~6 min. `--arm zc_noop` must still be BYTE-IDENTICAL to `base`. ⚠ Every arm file predating 2026-08-11 lacks the `messages` field, so `arm_identity.py` correctly refuses to compare one against a fresh run |
| `design/length_ceiling.py` | ⭐ **what a PERFECT length model is worth, on the ladder, ONE pmf at a time.** ⚠ This prices the fl PMFs the OPPORTUNITY model uses — NOT the length-likelihood composition channel, which is deferred past 0.8.0 and must not be proposed. ⛔ Pricing the pair together hid a 14× split between them (`TRAPS.md` B21). Reports the pass-0 solvable yardstick beside the mass-weighted one, because they disagree |
| `design/toy_harness.py` | ⭐⭐ **a mini chromosome you define, calibrated in 0.1–5 s** (`docs/TESTING.md` §0b), with every object's answer beside per-object truth. The library-level priors a toy cannot fit are harvested from a real cached condition and INJECTED; the gDNA depth is DERIVED to match that donor, never set by hand. `--list` for the spec ladder. ⭐ Reach for this FIRST when a mechanism needs isolating — it found C1's mechanism candidate in one sweep |
| **the substrate — is the panel and the index sound?** | |
| `design/simulator_gates.py` | ⛔ the simulator's own gates G-S1…G-S6, scored on per-fragment truth. Run before trusting the panel |
| `design/suite_resolves.py` | ⛔ can the suite resolve the axis you are changing? Run before quoting any suite number |
| `design/index_census.py` | re-derive an index's census — never quote a stored table, run this |
| `design/verify_index_rebuild.py` | regions byte-identical, boundaries only in contiguous reach |
| **Stage A — the accumulator** | |
| `design/fl_anchor_gap.py` | the anchor and both length models vs truth. `--drain` measures before and after the second pass |
| `design/gdna_pool_census.py` | ⭐ the four gDNA pools, each against its own opportunity and against truth |
| `design/second_pass_accuracy.py` | the second pass scored PER FRAGMENT against the oracle BAM's read names |
| `design/observable_efficiency.py` | what fraction of the length information a storage choice keeps |
| `design/region_density_derivation.py` | the reciprocal-opportunity theorem, T0–T6, each perturbed |
| **diagnostics** | |
| `design/anchor_opportunity_census.py` | ⭐⭐ **IS A ZERO-COUNT ANCHOR'S DENSITY CLAIM TRUE OF ITS NEIGHBOURHOOD? — no solver runs.** Every slot's TRUE gDNA density re-derived from the origin-split oracle through the SHIPPED `build_node_geometry`, then the empty structurally-pure-gDNA population against what its own chain neighbours measure. ⛔ Its verdict is that the claim is false by **346×** under capture and true off it — which is NECESSARY for the "empty means no probe here" hypothesis and, as the panel arms then showed, not sufficient. ⛔⛔ **THAT VERDICT IS ABOUT 1,312 ANCHORS, NOT ABOUT `struct_lock`** — the docstring AND the mask both claimed "exactly `strand_evidence`'s `struct_lock`" until 2026-08-11 and the solver's mask is **15–23× larger** (30,423 vs 1,312 at `g00`), because it also holds every zero-count REGION. This instrument's population is `g1_locked ∧ REGION` — what the standing xfail wants `struct_lock` rescoped TO. ⭐ Both sizes now print on every row, and the measured-redundant `∧ intergenic` term is ASSERTED rather than deleted |
| `design/composition_evidence_census.py` | how much library mass reaches the solver with NO composition evidence. `--inject-kappa 0.5` is its falsification handle |
| `design/held_flux_census.py` | how often a held candidate has ZERO flux evidence, by cause |
| `design/prior_units_check.py` | the EM prior in fragment units vs the old incidence sum |
| `design/boundary_q_population.py` | ⭐⭐⭐ **IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND? — DRAINED, no solver runs.** `assemble_priors` divides a boundary's crossing INCIDENCE by `q = mass/count` to get a fragment count — but `q` is measured on the POOLED population and applied to gDNA and RNA separately. Feeds a PERFECT `f_g` from the origin-split oracle so the `q` conversion is isolated from calibration's own error, and scores `Σ count_c·q_pooled` against the conserved truth `Σ mass_c`, EDGE-only and on the TOTAL prior the EM actually consumes. ⛔ **Its verdict: the defect is real, is NOT driven by fragment length (the equal-length null shows `q_g` 0.633 vs `q_r` 0.523, a LARGER error than a 102 bp gap), and is bounded at ≤ 0.6 pp of composition — so record the bound, build nothing.** ⭐ The driver is PLACEMENT: gDNA crosses boundaries in long intergenic regions where `q → 1`, RNA in short exon regions. ⭐ Gated on the pooled `q` reproducing the shipped `mass_per_crossing` bit-for-bit and on the partition summing to the whole; ⚠ its Ⓖ4 records that ladder capture-ON is NOT an equal-length null, because capture manufactures a +15.33 bp gap |
| **plumbing** | |
| `design/native_parity_on_real_data.py` | the native-parity gate on real cfRNA at full scale |
| `design/scan_profile.py` | ns/fragment, regressed over several BAMs |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
