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

⛔ **It also says my §2c oracle row is mis-measured rather than contradicting the synopsis** (§2c caveat 2).

### 1b. ORDERING — forward-backward stays (owner, 2026-08-18)

> *"Some objects can receive two messages from two neighbours. The forward-backward message propagation
> ensures that every object gets BOTH its messages. Once the messages arrive, the objects can solve."*

⛔ My depth-ordered single pass is REFUSED and the reason is structural: an ordered pass delivers from one
side only. Forward-backward is what guarantees both arrive.

### 1c. THREE ARMS, NATIVELY (owner, 2026-08-18)

`{gDNA, RNA+, RNA−}` — AXIOM 0, and non-negotiable. Collapsing the two RNA arms into a total would destroy
the strand channel, which is the primary intron deconvolution on stranded data. ⛔ Asking the question at
all was a failure to apply AXIOM 0; it is recorded here so the next session does not re-ask it.

## 2. WHAT THE MEASUREMENTS SAY (2026-08-18; no solver in any of them)

### 2a. Imputed stretches are DEEP, so the chain stays

Depth = chain hops from an IMPUTED slot to the nearest gDNA-MEASURING slot. `g50 ss0.50 nrna_none`:

| | capture-OFF | capture-ON |
|---|---|---|
| gDNA-measuring mass | 4,839,750 (54.3 %) | 1,578,748 (15.3 %) |
| imputed mass | 4,075,084 (45.7 %) | 8,717,307 (84.7 %) |
| imputed mass at depth 1 | 29.0 % | 41.3 % |
| **imputed mass at depth ≥ 9** | **28.4 %** | 12.4 % |

⛔ A one-hop local imputation reaches 29 % of imputed mass and no more — the relay is required, exactly as
the owner said. ⭐ Capture INVERTS the ratio (54 % → 15 %): far less measured gDNA to propagate from, which
is why capture-ON is the hard regime.

### 2b. SIX hop types carry the whole problem

Ordered hops into an IMPUTED slot, by destination mass, `g50 ss0.50 nrna_none capture_off`:

| dst ← src | hops | dst mass |
|---|---|---|
| `B exon\|exon` ← `R exon` | 25,622 | 3,249,058 |
| `R exon` ← `B exon\|exon` | 10,483 | 3,244,078 |
| `R exon` ← `B exon\|intron` | 10,766 | 1,094,494 |
| `R exon` ← `B gene edge` | 1,597 | 524,608 |

plus the two the synopsis names among gDNA-measuring objects (`B exon|intron` ← `R intron`,
`B gene edge` ← `R intergenic`). **Six.** Against nineteen switches.

### 2c. The currency oracle — the principle tested against certified truth

Transport the source's TRUE value both ways; score both at the destination in the same units
(`Σ|f_hat − true_f_g|·M`, fragments). The source is perfect, so what remains is the RULE's error.
`g50 ss0.50 nrna_none capture_on`:

| dst ← src | dst mass | LEVEL err | COMP err | winner |
|---|---|---|---|---|
| `B exon\|intron` ← `R intron` | 861,057 | 854,338 | **0** | **COMP** ⭐ |
| `B gene edge` ← `R intergenic` | 156,156 | 155,286 | **0** | **COMP** ⭐ |
| `R intron` ← `B exon\|intron` | 107,030 | 203 | **0** | COMP ⭐ |
| `B gene edge` ← `R exon` | 125,065 | **480** | 20,910 | LEVEL ⭐ |
| `B exon\|intron` ← `R exon` | 885,239 | **35,024** | 194,511 | LEVEL |
| `R exon` ← `B exon\|exon` | 5,577,039 | **560,426** | 636,431 | LEVEL |
| `B exon\|exon` ← `R exon` | 3,818,045 | 506,398 | 494,196 | tie |
| `R exon` ← `B exon\|intron` | 3,173,458 | 1,152,716 | 1,353,964 | ⛔ both 36–43 % |
| `R exon` ← `B gene edge` | 1,925,819 | 676,760 | 918,261 | ⛔ both 35–48 % |

⭐ The principle predicts the winner on every row that has one, and the synopsis's central claim —
COMPOSITION, not level, from an intron into its `exon|intron` boundary — is exact to the fragment where the
level is off by 854,338.

⛔⛔ **Two caveats travel with this table.**
1. `nrna_none` makes every gDNA-measuring object pure gDNA on both sides, so `f_g = 1 = 1` and those COMP
   zeros are partly `nrna = 0` restated. **Re-run on `nrna_mid` before treating them as measurements.**
2. ⚠ **The `R exon ← B exon|intron` row is MIS-MEASURED, and §1a says exactly how.** That direction is
   SPLICE IN, so the source population must be *unspliced crossings + the spliced fragments that splice in*.
   The oracle used the boundary's `slot_truth` composition, which is neither. **Re-derive it as the
   splice-in composition before drawing any conclusion from that row.**

⭐ **Not a caveat: at capture-ON both currencies fail into an exon (35–43 %).** That is the capture problem
— a partially-enriched boundary feeding a fully-enriched exon — and no currency choice fixes it. It is
`ROADMAP.md` §1 rank 3 (spike-and-slab): a third MECHANISM, not a third currency.

### 2d. ✅ THE SEQUENCING BLOCKER IS CLEARED (2026-08-19)

`calibration_oracle.py`'s FIELD certification passed on only **24 of 32** conditions — the 8 that failed
were exactly the dense capture-OFF ones, i.e. where unstranded × capture-OFF carries its gDNA mass, so the
currency oracle could not run where it matters most.

⭐ **Cause found and repaired**: a transcript-level chimera predicate silently dropped 4,087 gDNA fragments
per condition — the heaviest multi-boundary crossers. The ruling is `DESIGN.md` §3.1b-ii
(compatibility before chimera), the lesson is
`TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`, and the numbers are `ROADMAP.md` §0.
⚠ **What remains is bookkeeping**: rebuild the scan and oracle caches and re-stamp `slot_truth.npz`, after
which FIELD certification should read 32/32 and §2's currency oracle is unblocked on capture-OFF.

## 3. THE PLAN

| stage | what | acceptance |
|---|---|---|
| **0 · VOCABULARY** | §5's review, ruled by the owner, then ONE rename (never piecemeal), gated by `rename_identity.py` byte-identity | every stage proven a numeric NO-OP; suite green; `DESIGN.md` §0 carries the ruling |
| **1 · UNBLOCK** | ✅ **DONE 2026-08-19** — the deficit was a silent chimera drop (§2d) | FIELD certification 24/32 → 32/32, pending the cache rebuild |
| **2 · THE MAP** (no `src/`) | (a) depth + hop census, all 32, per stratum; (b) currency oracle, all 32, with §2c caveat 2's splice-in fix, per hop type × stratum × nascent level; (c) the residual census — hop types where `min(LEVEL, COMP)` is still large | one measured currency table; every place a third mechanism is REQUIRED is named before one is built |
| **3 · THE POLICY** | a THIRD policy beside `SilentPolicy` and `RelayPolicy` on the same gated backbone | the toy rungs below, each exact, both zero controls, before the next starts |
| **4 · THE PANEL** | all 32, three policies, RAW COUNTS (est / true / misplaced fragments), bar written first | never a bare ratio — the ratio framing hid a 1,978,148-fragment win for eleven days |
| **5 · RETIRE** | delete `RelayPolicy` and everything the panel proves dead | §3's predictions, checked not assumed |

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

**The rungs.** Each is a toy spec, seconds to solve, scored per object against truth; each must pass BOTH
zero controls (zero-RNA ⇒ `f_g = 1` everywhere; zero-gDNA ⇒ `f_g = 0`) before the next begins.

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
