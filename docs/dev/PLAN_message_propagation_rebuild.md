# PLAN OF ATTACK — rebuild message propagation as a third policy, simple by construction

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** Rulings go to `DESIGN.md`, lessons to
    `TRAPS.md`, derivations to `EQUATIONS.md`, current numbers to `ROADMAP.md` §0. ⛔ Delete this file when
    the plan is executed.

    Written 2026-08-18 (session 2). The DESIGN INPUT is the owner's synopsis in `docs/dev/message_notes.md`
    plus the owner's answers below; where they and this file differ, they win.

## 0. ⛔ VOCABULARY FIRST — the owner's ruling, 2026-08-18

> *"We MUST fix our vocabulary and terminology. This is THE HIGHEST PRIORITY … I would review our naming
> conventions and clean them up as the first step in the plan."*

RULED: **SPLICE IN** replaces `splice_in`. **SPLICE OUT** replaces `deconvolve`. The land/sea analogy the owner used
to explain the geometry is **for understanding only and must not enter code or docs** — this file uses the
predicate instead: an object is **gDNA-MEASURING** where no mature transcript crosses it (the solver's own
`¬mrna_active`) and **IMPUTED** where one does. §5 is the naming review this ruling asks for.

## 1. THE PRINCIPLE — two currencies with complementary invariances

The synopsis gives a currency per hop type. One reason sits under all of them, and stating it turns a table
to memorise into a rule to apply:

| currency | what it is | INVARIANT to | DESTROYED by |
|---|---|---|---|
| **LEVEL** `rho_g` | gDNA fragments per placement | **POPULATION** changes — TSS, TES, strand flips, splicing. gDNA is genomically continuous and knows nothing about transcripts | **ENRICHMENT** changes — a probe edge between the two objects |
| **COMPOSITION** `f_g` | the gDNA share of the crossing population | **ENRICHMENT** changes — a probe enriches everything overlapping it, both components alike, so the ratio survives | **POPULATION** changes — the denominator is a different set of molecules |

⭐ **They are complementary, so at every hop at least one is intact: use that one.** The licence question
*is* the currency question, asked once per hop type. No licence hierarchy, no three-case rule, no
per-strand patch.

⛔ **This is what the shipped relay has backwards.** `RelayPolicy` transports a COMPOSITION by default (the
reframe `r = rho_tot(dst)/rho_tot(src)`) and bolts population licences on top — but the population is
exactly what changes at the `exon|exon` boundaries carrying the largest destination mass. Every licence bug
of the last two sessions repaired a currency CHOICE instead of changing it.

### 1a. ⭐⭐⭐ THE POPULATION OF A MESSAGE IS DIRECTION-DEPENDENT (owner, 2026-08-18)

> *"Always answer the question: what crosses INTO this region?"*

| direction | what crosses into the destination | spliced fragments |
|---|---|---|
| BOUNDARY → EXON — **SPLICE IN** | unspliced boundary crossings (gDNA + nascent RNA) **+ spliced fragments, which splice IN to this exon** | **INCLUDED** |
| EXON → BOUNDARY — **SPLICE OUT** | unspliced fragments crossing contiguously | **EXCLUDED** — they splice OUT and land elsewhere |

⭐ So SPLICE IN and SPLICE OUT are not two bolted-on operators; they are **one rule — the message's
population is whatever physically enters the destination — evaluated in the two directions.** That is the
whole justification for both, and it is why the pair must exist.

✅ **Measured as stated (2026-08-19, §2c)**: with the SPLICE-IN population — unspliced + contiguous spliced +
the sj flux on the exon's side — the composition into an exon from a splice site is at its noise floor off
capture and the clear winner on it; the prototype's "mis-measured row" was this population defect PLUS the
terminus/splice-site conflation §2b names.

### 1b. ORDERING — forward-backward stays (owner, 2026-08-18)

> *"Some objects can receive two messages from two neighbours. The forward-backward message propagation
> ensures that every object gets BOTH its messages. Once the messages arrive, the objects can solve."*

⛔ My depth-ordered single pass is REFUSED and the reason is structural: an ordered pass delivers from one
side only. Forward-backward is what guarantees both arrive.

### 1c. THREE ARMS, NATIVELY (owner, 2026-08-18)

`{gDNA, RNA+, RNA−}` — AXIOM 0, and non-negotiable. Collapsing the two RNA arms into a total would destroy
the strand channel, which is the primary intron deconvolution on stranded data. ⛔ Asking the question at
all was a failure to apply AXIOM 0; it is recorded here so the next session does not re-ask it.

## 2. WHAT THE MEASUREMENTS SAY — ✅ STAGE 2 RAN 2026-08-19 ON THE REBUILT CACHES

⭐ **This section used to carry pre-fix numbers and said so; it now carries the MEASURED map, and the
numbers themselves live in `ROADMAP.md` §0's "HOP CURRENCY MAP" row (moved, not copied).** The instrument
is `scripts/design/hop_currency.py` (all 32 conditions in ~18 s; `--out` writes every condition × hop-type
row; `--self-test` 36/36, falsified by perturbation). Read its docstring for the full table.

### 2a. Imputed stretches are DEEP, so the chain stays — CONFIRMED

Capture-OFF: 27–29 % of imputed mass is ≥ 9 hops from any measured gDNA (g05/g50; 14 % at g98) and
depth 1 reaches only 28–34 %. Capture-ON inverts the measured/imputed ratio (62 % → 18 % measured at g50)
with 40–43 % at depth 1 and 10–15 % at ≥ 9. Forward-backward over the chain stays (§1b).

### 2b. ⛔⛔ THE HOP TYPES ARE NOT THE SEVEN STRATA — the data named the key

Classified by object class alone, `R exon <- B exon|intron` read COMPOSITION 221 k vs LEVEL 27 k on one
capture-OFF condition — as if the SPLICE-IN population fix had failed. Split by the boundary's structural
flags (`is_terminus` / `is_splice_site`, either strand) **all of it sits on the TERMINUS boundaries** — a
TSS/TES lying inside another transcript's intron (3,867 `[term]` vs 15,480 `[sj]` exon|intron boundaries;
`B exon|exon` splits 7,850 `[term]` vs 4,838 `[sj]`). **The hop-type key is `object class × {sj, term, sj+term}`**
(`TRAPS: an-object-class-does-not-see-a-terminus`), and the static table Stage 3 builds in `prepare()`
keys on it. So "six hop types" becomes ~22 ordered types, of which the policy needs about five RULES.

### 2c. The currency oracle — the principle held, and it named one residual

Every error beside its Monte-Carlo noise floor (`TRAPS: the-floor-must-reproduce-the-selection`), in
fragments, per condition:

| hop | currency | evidence |
|---|---|---|
| exon ← SPLICE-SITE boundary (SPLICE IN) | **COMPOSITION** of what ENTERS (unspliced + spliced + sj flux on the exon's side) | both at floor off capture; under capture COMP 1.2–1.5 % at g05/g98, 3.7–6.5 % at g50, against LEVEL 3 → 30 → 58 % |
| exon ← TERMINUS boundary (gene edge, or a TSS/TES inside a gene) | **LEVEL** off capture; ⛔ **NEITHER on capture** | COMP 40–98 % off capture (the RNA originates there); on capture LEVEL dies by the enrichment (3 / 34–37 / 62–64 %) and COMP by the population (92 / 44 / 1.4 %) |
| boundary ← EXON (SPLICE OUT direction), any class | **LEVEL** | 0.0 % on every arm; COMP 3–70 % (the exon's RNA ends or splices out); under capture the claim is CLIPPED to "all gDNA", which is true of the boundary |
| exon\|intron / gene edge ← its own intron / intergenic | **COMPOSITION** | exact on `nrna_none` where LEVEL is off by 78–98 % of the destination under capture — the synopsis's central claim, to the fragment. ⚠ 5–37 % on `nrna_mid × capture-ON` is the PANEL (`TRAPS: the-panel-enriches-nascent-by-its-own-probes`), not the rule |
| intron / intergenic ← its boundary | either | the region solves itself |

⛔ **The residual census is one thing: a TERMINUS into an EXON UNDER CAPTURE** — 0.57–0.65 M fragments per
hop type at g50 (gene edge, exon|intron[term], and 0.24–0.27 M at exon|exon[term]), every strand, every
nascent level. The enrichment changes AND the population changes, so no currency survives; that is
`ROADMAP.md` §1 rank 3's spike-and-slab — a third MECHANISM, re-established here rather than inherited.
⚠ A g98 COMP "win" at those hops (1.1–1.5 %) is `f_g = 1` restated, the mirror of `nrna = 0` at the
controls; read the row across the three gDNA levels, never one cell.

### 2d. ✅ THE SEQUENCING BLOCKER IS CLEARED (2026-08-19)

`calibration_oracle.py`'s FIELD certification passed on only **24 of 32** conditions — the 8 that failed
were exactly the dense capture-OFF ones, i.e. where unstranded × capture-OFF carries its gDNA mass, so the
currency oracle could not run where it matters most.

⭐ **Cause found and repaired**: a transcript-level chimera predicate silently dropped 4,087 gDNA fragments
per condition — the heaviest multi-boundary crossers. The ruling is `DESIGN.md` §3.1b-ii
(compatibility before chimera), the lesson is
`TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`, and the numbers are `ROADMAP.md` §0.
✅ **The rebuild is DONE** (2026-08-19): 32/32 scan caches, 32/32 oracle caches, all 32
`slot_truth.npz` stamped COMPOSITION + FIELD, zero class-bias flags. ⚠ **Read that stamp honestly**: the
uniformity gate is **VACUOUS on the 16 capture-ON conditions** by design (capture makes the field
deliberately non-uniform), so the real verification is **16/16 capture-OFF, up from 8/16** — which is
still exactly the win, since all eight failures were capture-OFF. **§2's currency oracle is now
unblocked on the capture-OFF strata**, which was the whole point of Stage 1.

## 3. THE PLAN

| stage | what | acceptance |
|---|---|---|
| **0 · VOCABULARY** | §5's review, ruled by the owner, then ONE rename (never piecemeal), gated by `rename_identity.py` byte-identity | every stage proven a numeric NO-OP; suite green; `DESIGN.md` §0 carries the ruling |
| **1 · UNBLOCK** | ✅ **DONE 2026-08-19** — a silent chimera drop (§2d); panel fully rebuilt | ✅ the field gate on the conditions where it RUNS: **8/16 → 16/16 capture-OFF**, zero class-bias flags (the 16 capture-ON rows are vacuous by design) |
| **2 · THE MAP** (no `src/`) | ✅ **DONE 2026-08-19** — `hop_currency.py`, all 32, §2 above | ✅ the currency per hop type is MEASURED (§2c), the hop-type key is `object class × {sj, term, sj+term}` (§2b), and the one place both currencies fail is NAMED: a terminus into an exon under capture |
| **2b · THE MAP PER ARM** (no `src/`) | ⛔ **OWNER-RULED ESSENTIAL 2026-08-19**: extend `hop_currency.py` to `{gDNA, RNA+, RNA−}` — the 2026-08-19 map pooled the RNA strands | three currencies per hop type, re-run on the re-simulated caches (step 4 below) |
| **2c · THE SIMULATOR** (owner's step 2) | ✅ **DONE 2026-08-19** — index-backed transcriptome incl. nascent entities, one RNA multinomial, capture by genomic overlap (`ROADMAP.md` §0) | ✅ gDNA 1162× / nascent 1003× under one probe; four perturbations fire; suite 3,538 passed |
| **2d · THE TOY HARNESS** (owner's step 3) | build/revive `toy_harness.py` on the index-based simulator so a toy is solved in seconds with the FL / capture / strand behaviour harvested from a cached full scenario | `TA+` single exon solves; then `TB+` multi-exon beside it |
| **2e · RE-SIMULATE** (owner's step 4) | all 32 conditions, caches and certification rebuilt (`panel.py`, `calibration_oracle.py`), the per-arm map re-run | 32/32 certified; the map re-measured |
| **3 · THE POLICY** (owner's step 5) | ⭐ **TOY-FIRST, ONE TRANSCRIPT AT A TIME** — a THIRD policy beside `SilentPolicy` and `RelayPolicy` on the same gated backbone, its arithmetic `EQUATIONS.md` §3.5e (three densities; rescale-by-totals at a compatible hop with SPLICE OUT / SPLICE IN; level only across a terminus, `sj+term` a terminus) | each added transcript solved exactly, everything before it still solved, both zero controls; the cached full scenarios as the benchmark whenever wanted |
| **4 · THE PANEL** | all 32, three policies, RAW COUNTS (est / true / misplaced fragments), bar written first | never a bare ratio — the ratio framing hid a 1,978,148-fragment win for eleven days |
| **5 · RETIRE** | delete `RelayPolicy` and everything the panel proves dead | §3's predictions, checked not assumed |

### 3a. ✅ STAGE 2 IN FULL — EXECUTED 2026-08-19; kept as the record of the spec the instrument meets

⛔ **NO `src/` CHANGE. NO SOLVER RUNS.** Everything here reads the certified caches and the annotation.
✅ Built as specified, with two additions the data forced: the hop-type key carries the boundary's
structural class (§2b), and every error is printed beside a noise floor that reproduces the scorer's own
selection (§2c). Both defects named below were fixed as written.

**Build ONE instrument**, `scripts/design/hop_currency.py`, to the repo's standard: a docstring that opens
with what it is FOR and its verdict, and a `--self-test` falsified by perturbation with no I/O. It is
picked up by `tests/test_scripts_index.py` — ⚠ a new `scripts/design/` file is worth **+4** collected tests.
⚠ It was prototyped in a session scratchpad that no longer exists; the spec below is complete, and §2's
tables are the SHAPE its output should have (not numbers to reproduce — those were pre-fix).

It answers three questions, on each of the 32 conditions, per stratum:

**① THE HOP CENSUS.** Classify every ordered adjacent pair `(dst ← src)` by the two slots' classes
(`object_composition.strata`'s seven labels — already partition-asserted). Report hop count and
DESTINATION MASS per type. ⭐ Also report the DEPTH distribution: hops from each mrna-active slot to the
nearest slot where no mature transcript crosses (`¬mrna_active`), mass-weighted. That is what says whether
the chain is needed at all — the pre-fix answer was emphatically yes, with ~28 % of that mass ≥ 9 hops
from any measured gDNA at capture-OFF.

**② THE CURRENCY ORACLE — the decisive one.** For every hop, transport the source's TRUE value BOTH ways
and score both at the destination in the SAME units, `Σ|f_hat − true_f_g|·M` in FRAGMENTS::

    LEVEL:        f_hat(dst) = rho_g_true(src) · E_g(dst) / M(dst)
    COMPOSITION:  f_hat(dst) = f_g_true(src)

The source is PERFECT, so what remains is the RULE's error and the smaller column IS that hop type's
currency. ⛔ **Fix the prototype's two known defects**: (a) in the SPLICE IN direction the source population
must be *unspliced crossings PLUS the spliced fragments that splice in* (`DESIGN.md` §0c.0's direction
rule), not the raw `slot_truth` composition — the prototype used the latter and that row is unusable;
(b) run the `nrna_mid` rows too — on `nrna_none` every gDNA-measuring object is pure gDNA, so `f_g = 1 = 1`
and a COMPOSITION "win" there is `nrna = 0` restated rather than a measurement.

**③ THE RESIDUAL CENSUS.** Where `min(LEVEL, COMPOSITION)` is STILL large, neither currency works and a
third MECHANISM is required. The pre-fix run put that at the hops into an exon under capture (~35–43 %),
which is `ROADMAP.md` §1 rank 3's spike-and-slab — ⛔ re-establish it, never inherit it.

**ACCEPTANCE**: one table, per hop type × stratum × nascent level, naming the currency and its residual,
`--self-test` green, instrument committed. ⛔ **No policy code is written until that table exists** — the
entire point is that the currency is MEASURED rather than argued.

⚠ **Scoring discipline that has already cost one false finding**: score per gDNA-BEARING REFERENCE, never
pooled (`TRAPS: a-ratio-needs-a-population-that-can-supply-its-numerator`), and report FRAGMENTS, never a
bare ratio (`TRAPS: never-pool-the-strata`).

⭐⭐ **Why a third policy and not a rewrite.** `DESIGN.md` §6.1's backbone/policy seam already makes
`SilentPolicy` (5 lines) and `RelayPolicy` (1,200) peers behind one gated interface, with five backbone
assertions a policy cannot violate. `RelayPolicy` is NOT touched, so nothing that works can break and the
A/B is a config value rather than a diff. ⛔ The 2026-08 clean rebuild that came out **+103 %** replaced the
whole solver at once; this replaces one peer.

**The new policy's core is one static table.** In `prepare()`, from the annotation alone and no beliefs:
classify every ordered adjacent pair into one of §2b's hop types, look up its currency and its direction's
population rule (§1a), and ASSERT the classification partitions the hops — the discipline
`object_composition.strata` already applies to slots. Capture it so an instrument can census it and so
`TRAPS: an-ablation-that-never-ran` cannot recur.

**The rungs.** ⭐ **Owner's framing (2026-08-19): start with ONE single-exon transcript `TA+ (1000, 2000)`
and develop the policy around it; add a multi-exon `TB+ (5000, 6000), (10000, 11000)` and make sure `TA+`
still solves; keep adding transcripts and let the policy and the map evolve with each addition.** The
table below is the shape that ladder takes; each rung is a toy spec, seconds to solve, scored per object
against truth, and each must pass BOTH zero controls (zero-RNA ⇒ `f_g = 1` everywhere; zero-gDNA ⇒
`f_g = 0`) before the next begins.

| rung | spec | what it proves |
|---|---|---|
| R1 | `TA_single_exon` (exists) | a gDNA LEVEL reaches a single-stranded exon from its two gene-edge boundaries |
| R2 | `spliced_exons` (exists) | COMPOSITION `R intron` → `B exon\|intron`, then SPLICE IN to the exon |
| R3 | `nascent` (exists) | the intron deconvolution (density model vs the intergenic background, NB) is correct UNSTRANDED |
| R4 | **new** — ≥ 9 hops of `exon ↔ exon\|exon ↔ exon` | LEVEL propagation inward from both ends; §2a's 28 % |
| R5 | **new** — the owner's locus verbatim (TA+ TB+ TC+ TD− TE− TF+) | overlapping opposite strands, alternate TSS/TES, nested isoforms |
| R6 | R1–R5 × capture-ON | where §2c says both currencies fail; expected to need rank 3's mechanism |

### What §1's principle predicts dies at Stage 5 — predictions to check, not a plan to execute

| switch | why it should not survive |
|---|---|
| `terminus_population`, `strand_population`, `gdna_level_scale`, `rna_level_scale` | these four ARE the hop-type table, badly factored — one lookup replaces them |
| `mass_rescale` | see §5's entry: it forces a message's three densities to add up to the destination's observed count. Under either currency the claim already accounts for `M` by construction, so there is nothing to restore |
| `conservation_var` | prices the residual `mass_rescale` discards; dies with it |
| `flank_pair`, `transfer_var`, `splice_in_frame_var` | all are properties of `r`, so they survive only on COMPOSITION hops and may be far smaller scoped there |
| `lam_emission_gate` | "this source has no composition to lend" is a hop-type fact, i.e. the table |
| **SURVIVES** | SPLICE IN (the sj-flux RNA level — the operator `DESIGN.md` §0c says the exon needs), SPLICE OUT, the intron deconvolution, ψ's four channels, every precision |

## 4. WHAT THE OWNER HAS ALREADY RULED (do not re-ask)

* SPLICE IN / SPLICE OUT, and the direction rule that justifies them (§1a).
* Forward-backward stays (§1b).
* Three arms, natively (§1c).
* The land/sea analogy is an explanation, never code or docs (§0).
* The 2.5 % deficit is likely an accumulator bug, worth a dedicated session (§2d adds only that it
  sequences before the design measurement, not that it is more important).

## 5. THE NAMING REVIEW — measured, with a proposal per name

Occurrences across `src/` · `tests/` · `scripts/` · `docs/` · memory, and how many are CODE IDENTIFIERS
(the rest are prose, which `rename_census.py` shows is ~43 % of any rename):

| name | total | code ids | what it actually does | proposal |
|---|---:|---:|---|---|
| `splice_in` | 163 | 30 | a flanking BOUNDARY's measured sj flux joins the RNA claim entering an EXON | ✅ **RULED → SPLICE IN** (`splice_in`) |
| `deconvolve` | 96 | 7 | an EXON's RNA leaving into a BOUNDARY is scaled by the continuing share | ✅ **RULED → SPLICE OUT** (`splice_out`) |
| `pin` / `mass_rescale` | 207 | 19 | rescales a message's three densities by a common factor so they exactly account for the destination's OBSERVED fragment count | ⛔ **NEEDS A RULING.** Proposal: `budget_rescale`, or delete it (§3 predicts it is unnecessary) |
| `head` / `RelayPolicy` | 145 | 22 | the name of the shipped policy. Carries no meaning at all | ⛔ **NEEDS A RULING.** It is being retired; the NEW policy should be named for what it does, e.g. `ChainPolicy` |
| `reframe` / `frame` | 91 / 192 | 3 / 26 | multiply a source's claim by `rho_tot(dst)/rho_tot(src)` — i.e. transport a COMPOSITION | proposal: `composition_transport` / `share_transport`; `r` → `total_ratio` |
| `lend` | 71 | 17 | "may this source supply a composition to this destination?" | proposal: `may_share_composition` |
| `hop_logvar` | 25 | 20 | `Var(log r)` — the transport's own scale-sampling variance | proposal: `transport_logvar` (it already has a spelled-out twin, `transfer_logvar`) |
| `premise` | 68 | 0 | the variance of the assumption a hop rests on | prose only; keep or fold into the above |
| `arm` | 1,358 | 235 | ⚠ **TWO SENSES** — an experiment arm (A/B) and a component arm (gDNA/RNA+/RNA−). `TRAPS: two-masks-one-name` | ⛔ **NEEDS A RULING.** Proposal: keep `arm` for experiments, use `component` for {gDNA, RNA+, RNA−} |
| `relay` | 344 | 25 | the message-passing machinery | probably fine; owner call |
| `hop` | 159 | 0 | one step between adjacent objects | prose only, and clear; keep |
| the land/sea analogy | 26 | 0 | an explanation of the geometry | ✅ **RULED OUT of code and docs.** All 26 removed 2026-08-18 — 17 introduced by this file, 3 elsewhere in mine, 1 pre-existing in `TRAPS.md`; the rest are `high-water mark`, a standard memory-profiling term, correctly left alone |

⛔ **One rename, not piecemeal** (`CLAUDE.md`), gated stage-by-stage on `rename_identity.py`'s byte-identity
— a rename is a numeric NO-OP and that must be falsifiable, since the stages COMPOUND.
