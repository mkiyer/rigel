# DESIGN — what is built, and the rulings behind it

**Everything here is SETTLED.** These are decisions that have been made, measured, and are not to be
re-litigated. Where a decision was reached by refuting an alternative, the refutation is in `TRAPS.md`,
not here. Derivations are in `EQUATIONS.md`.

---

## 0. ⭐⭐ VOCABULARY — one word per concept, owner-set 2026-08-04

⛔ **This table is binding on prose, comments, docstrings and identifiers alike.** Every banned word below
was in use for a concept that already had a name, and the ambiguity cost a reader real time.

| ✅ the term | what it means | ⛔ the banned synonym |
|---|---|---|
| **REGION** | a **contiguous genomic interval** — a region of the partition. Has a length in bp | `region` is tolerated: it is the index's own word for the same thing |
| **BOUNDARY** | a **single genomic position** separating two adjacent REGIONs. Zero bp wide | `line` · `seam` · `crossing` (as a noun for an object) |
| **edge** | ⭐ KEPT, as the sanctioned GRAPH-EDGE sense — a BOUNDARY and/or an SJ (owner ruling in the rename's stage 8, `391aa5eb`), which is what makes `edges_df`, `EDGE_KIND_CONTIGUOUS` and `EDGE_KIND_SJ` correct rather than deferred | ⛔ never as a synonym for BOUNDARY in prose |
| **SPLICE JUNCTION** (`sj`) | ⚠ a splice junction — a *different* object on a *different* axis, directed `src → dst` in **GENOMIC** order (`src < dst`, always). Never a synonym for BOUNDARY | — |
| **the sj's LOW end / HIGH end** — `_lo` / `_hi`, or `boundary_left` / `boundary_right` | the two BOUNDARIES a sj is anchored at, named in GENOMIC order, which is the order the code stores them in | ⛔ `donor` · `acceptor` — see below |
| **slot** | one entry of the chain, which alternates REGION, BOUNDARY, REGION, BOUNDARY … A slot is a REGION **or** a BOUNDARY | — |
| **step** | one adjacency move along the chain: REGION→BOUNDARY or BOUNDARY→REGION. So REGION→BOUNDARY→REGION is **two steps** | `hop` |
| **structurally pure-gDNA object** (**G1 object**) | a slot at which no RNA strand is admissible ACCORDING TO THE ANNOTATION, so its composition is certain *given the annotation*: an intergenic REGION, or an `intergenic\|exon` BOUNDARY. Its gDNA density is directly observed, with nothing to deconvolve. The predicate is `region_geometry.g1_locked`. ⚠ On real data the certainty is annotation-derived rather than physical — intergenic space carries unannotated transcription the tool does not model (§0b) — which is safe for a POOLED level and is not a licence to treat one such slot as ground truth | `anchor` — ⛔ that word had two meanings at once and now has only the one below |
| **the mass pin** | the operator that rescales a message so that `Σ_c ρ_c·E_c = M` at the destination (`messages/relay.py`'s `_rescale_v` and its scalar twin inside the scan kernel). "Pin" because the function is named for it | `the mass anchor` |
| **counts** | discrete integer fragment counts | — |
| **density** = **abundance** | counts per base. The two words mean the same thing | ⛔ not the simulator's molar `abundance=` field, which is a per-transcript weight |
| **crossing fragment** | ⭐ a **fragment** that spans a BOUNDARY. Legitimate and necessary — `crossing_eff_length` is the opportunity for exactly this — and it stays | ⛔ only the *object* sense is banned: objects are REGIONs and BOUNDARIES, never "crossings" |
| **switched off** | an A/B in which one code path is disabled and the run repeated, to establish that it is the cause | `ablated` |
| **splice-out** | ⭐ owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA departs** through the sj there. The fragments that splice away leave the contiguous population; the residual continues | ⛔ `splice_out` |
| **splice-in** | ⭐ owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA arrives** through the sj there — the spliced flux joins the destination's population | ⛔ `splice_in` |

⚠ **THE TWO BOUNDARY ROWS OF THE TABLE ABOVE WERE CORRUPTED BY THE RENAME THIS TABLE ORDERED, AND THAT
IS `TRAPS: the-rename-that-corrupted-a-diagram` MET INSIDE THE RULING ITSELF** (repaired 2026-08-17). The BOUNDARY
row's banned-synonym cell read *"`boundary` · `boundary` · `crossing`"* — it banned the very word the row
endorses, twice — because stage 8b mapped both `line` and `seam` onto `boundary`; and the row beneath it
read *"BOUNDARY is an accepted synonym for BOUNDARY"*, which was true of `EDGE` and is vacuous of
`BOUNDARY`. The banned list is restored as measured; the second row now carries the `edge` ruling the same
commit made. ⛔ Nothing else in this section was re-decided.

⛔⛔ **`donor` / `acceptor` ARE BANNED, AND WHAT BANS THEM IS A MEASUREMENT** (owner ruling, carried here
2026-08-14). They are 5′/3′ names and therefore **strand-dependent**, while every structure in this tree is
stored in GENOMIC order. Measured on the suite index: `src < dst` holds on all **13,482** sj including all
**6,527** minus-strand ones — the code is CORRECT — so an identifier such as `donor_cut` holds the
**ACCEPTOR** for **48.4 %** of sj. ⛔ That is worse than a vague name: **the name lies about correct code**,
and a reader who trusts it writes a sign error. Say `boundary_left` / `boundary_right`, or `_lo` / `_hi`.

⚠ **The index's structural flag bits are the one place the words survive, and they are the hazard itself
rather than an exemption.** `FLAG_DONOR_s` marks the genomic-LOW end of an `s`-strand intron on **both**
strands, so on `−` it sits at the transcript's biological ACCEPTOR. `region_geometry.py` therefore states
this at every consumer and names its arrays `_lo` / `_hi`; `test_splice_flux_reframe` gates it on a
`−`-strand sj specifically. ⚠ `TESTING.md` §0b records the same trap from the toy side, where the bits
decide the sign of a whole scenario.

⭐⭐ **`splice-out` / `splice-in` are DIRECTIONAL, and that is the whole reason for the rename.** "Deconvolution"
and "graft" named two operators; the same BOUNDARY is a splice-out for a message travelling one way and a
splice-in for a message travelling the other, so the pair names one thing seen from two sides.
⭐⭐⭐ **LANDED 2026-08-18 (owner: vocabulary is the highest priority).** `graft` → `splice_in` tree-wide
(verified unambiguous); `deconvolve` → `splice_out` **SENSE-SCOPED**, because `deconvolve` also names the DECONVOLUTION
verb — *"deconvolve the gDNA off, RNA is the residual"* — in `density_deconv`, `calibrate`, `object_composition`
and in the named rule `TRAPS: the-deconvolution-is-as-good-as-the-density-it-is-handed`. ⛔ That second sense is
DELIBERATELY UNCHANGED: a global replace corrupts it, which is `TRAPS: two-masks-one-name` and is exactly
what `rename_census.py --sense` exists to catch. ⚠ The deconvolution sense has NOT been ruled on and still
carries a retired word — an owner call, and the next naming decision.
⭐ The rename was gated as a numeric NO-OP: full suite green, **no `tests/golden/` file moved**, and the
diff exactly 290 insertions / 290 deletions.

⭐⭐⭐ **AND THE DERIVATION SHARPENED THAT INTO A RULING (2026-08-05, `EQUATIONS.md` §3.6c).** The two words
name the SEMANTICS of a step and **not** its arithmetic. An BOUNDARY presents one total to its genomic-LOW
neighbour and another to its genomic-HIGH one, and a hop between adjacent slots always uses the low slot's
HIGH-flank total against the high slot's LOW-flank total — *whichever* of them is the source. So the
message direction decides only whether the step is called a splice-out or a splice-in; it never decides
which number is used. ⛔ **A predicate on the message direction is therefore the wrong shape for anything in
this family, and one on the SIDE is the right one.** ⚠ The corollary that matters when writing the plumbing:
this is still not expressible as one array per pass, because within a single forward pass a BOUNDARY presents
its LOW-flank total as a destination and its HIGH-flank total as the very next hop's source.

⚠ **`docs/TESTING.md` §0b carries the counts/density half of this table** for readers who arrive there
first; this section is the canonical one.

⚠⚠ **The banned words are still widespread in code written before this ruling** — ~1,500 occurrences of
`line` / `seam` / `crossing` / `anchor` across `src/`, `tests/` and `scripts/`, almost all in comments and
docstrings. They are being converged as each area is touched, not in one sweep: a blanket regex over that
much prose mis-fires (`crossing_eff_length` and "crossing fragment" are correct; the fragment-length
"anchor" is a third, unrelated meaning), and a 1,500-line mechanical diff is not reviewable.
⚠ **That ~1,500 is a 2026-08-04 measurement and PREDATES the nine-stage rename** (`391aa5eb` and its
siblings), which converged most of it — and did the damage recorded above. Re-derive with
`scripts/design/rename_census.py`; do not quote the old number as current.

---

## 0b. ⭐⭐⭐ THE 0.8.0 RELEASE SCOPE — owner ruling, 2026-08-14

⛔ **This is the organising frame for the whole phase, and it is settled.** The shipped version is
`0.7.1` (`pyproject.toml`); the target is **0.8.0**. `ROADMAP.md` ranks the work *inside* this scope;
what is in it and what is out of it is decided here.

### Three strata are the optimisation target, and the fourth is DEFERRED

The panel crosses **strandedness × capture**, and the four cells are not equally healthy.

> ⭐ **IN SCOPE FOR OPTIMISATION — three strata:** unstranded × capture-OFF, stranded × capture-OFF,
> stranded × capture-ON.
>
> ⛔ **DEFERRED — unstranded × capture-ON.** It is not a development target until the other three are
> fully optimised.

⛔⛔ **DEFERRED IS NOT DROPPED, and that distinction carries the ruling.** The deferred stratum **stays in
every benchmark and every measurement, and must keep being reported**: a panel that cannot see the cell
the tool is worst at cannot tell a real win from a re-labelling, and it is the one cell where a mechanism
that helps everywhere else can hide a catastrophe. ⭐ If it improves as a side effect of work on the other
three, that is a free win — it is never the *justification* for a change. ⛔ It follows that every score
is read **PER STRATUM and never pooled**: a panel total that averages the deferred cell into the three
in-scope ones states neither quantity.

⭐ **Why this shape.** Measured 2026-08-13/14 on the rebuilt 16-condition ladder (noise floor
0.996–1.013): the error is CONCENTRATED — unstranded × capture-ON carries **64.5 % of transcript error
and 90 % of gene-level error**, and the three in-scope strata carry the rest. ⛔ And that cell is not a
gradient anyone can descend yet: on it the tool emits a near-ZERO gDNA fraction *regardless of truth*
(exon `f_g` 0.040 / 0.0016 / 0.0021 at `g05` / `g50` / `g98` against truths 0.054 / 0.518 / 0.982), so it
looks acceptable at low gDNA **by coincidence**. A stratum whose answer barely responds to the truth is
one an optimiser overfits rather than fixes.

### ⛔⛔ The LENGTH CHANNEL is deferred until after 0.8.0

> ⛔ **The fragment-length / "length likelihood" channel as a CALIBRATION COMPOSITION channel is future
> work. 0.8.0 ships WITHOUT it. Do not propose it, do not list it as a candidate, do not rank it.**

⚠ **It does not exist in `src/`.** It was A/B'd once, on 2026-08-10, and never shipped — so this is a
scope ruling and not a code removal, and any doc presenting it as "the next thing to try" is stale.

⛔ **`length_likelihood` in `src/rigel/second_pass.py` is a DIFFERENT thing and is NOT affected** by this
ruling: it is the per-fragment second-pass assignment factor of §4 below, the `f(L)` in
`score = ρ × f(L) × s`. Same words, other axis.

⚠ The two rulings interlock: the optimisation panel gives both origins the SAME length distribution (the
ladder ruling below), so it is structurally incapable of pricing a length composition channel
(TRAPS: equal-lengths-carry-no-composition). The panel and the scope agree.

### CALIBRATION is the focus, and the metric is calibration against ORACLE CALIBRATION

> ⭐ **The focus of development is CALIBRATION.** The metric is the **calibration result scored against
> ORACLE calibration** — not the end-to-end transcript number.

⭐ The end-to-end number stays a thermometer rather than a target (`SUCCESS.md` says why): it moves for
reasons that have nothing to do with calibration — the EM's own evidence, the fragment-length scorer,
depth — so a change to calibration must be scored where it acts. ⛔ The oracle is the origin-split truth,
and the instruments that already own that comparison are `solvability_audit.py`, `pass0_vs_oracle.py` and
`prior_vs_oracle.py`; a calibration change is not measured until one of them has spoken.

⭐ **Scenarios must be CACHED so calibration can be re-run extremely efficiently.** A scan is minutes and
a re-solve is seconds, and every mechanism this project has found came out of the fast loop —
`scripts/design/build_scan_cache.py` and the `cache` stage of
`scripts/sim/panel.py`, which builds the oracle cache beside the scan one.
⭐ **The pattern to copy from any future cache is the KEY**: the scan manifest plus a content hash of
every source file that PRODUCES what is stored, with the module under test deliberately left OUT so its
edit loop stays one second — ⛔ which is sound only while nothing that module produces is stored
(`TRAPS: a-hash-that-misses-its-artifact`, second form).

### ⭐⭐⭐ WHY THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS

> ⭐⭐ **THE EM ALREADY USES THE FRAGMENT-LENGTH DISTRIBUTION.** A large gDNA-vs-RNA length difference lets
> the EM assign fragments on LENGTH ALONE — **bypassing calibration entirely and MASKING bugs in it**.
> Equal lengths FORCE the calibration phase to be exercised.

⭐ It is structural rather than a preference: `scoring.cpp` carries **two** fragment-length lookup tables,
`fl_log_prob_` for RNA and `gdna_fl_log_prob_` for gDNA, so where the two distributions are far apart a
fragment's own length separates the components one fragment at a time and it barely matters what
calibration supplied. The panel would then read green with calibration broken.

⚠ **This reason is STRONGER than the one previously written down**, which was "the length channel is
neutralised, so residual error is attributable to density and strand". That is true and it says what the
panel MEASURES; it does not say what the panel PREVENTS. Where a doc still carries only the weaker form
it is incomplete rather than wrong.

⭐ **What calibration is left with is exactly STRAND and DENSITY** — plus belief propagation across
objects, which is ON (`RelayPolicy`, §0c.2) — and on an unstranded library it is density
alone, since the strand λ-term is exactly 0 at κ = ½ (§7 below). That is the substrate the three in-scope
strata are optimised on.

⛔ **"Equal" is a CONFIGURATION and never a guarantee**: identical parameters still leave a small realised
gap, because an mRNA fragment must fit inside its transcript and gDNA need not. `TESTING.md` carries the
measured residual and what it was priced at; ⛔ score the length axis against the simulator's own truth
table, never against a nominal parameter.

### ⭐⭐⭐ THE NASCENT SCOPE RULING — owner, 2026-08-22

**Real RNA-seq has nascent RNA at very low levels — sparse and rare.** Rigel models it, simply — one
synthetic transcript spanning every multi-exon transcript — so the tool is ROBUST to it. That is the
whole of nascent RNA's place in the design.

* **The default assumption and the default behavior: nascent RNA is absent unless there is an
  abundance of evidence otherwise.**
* **Experiments designed to capture nascent RNA explicitly are a DIFFERENT experiment and are OUT OF
  SCOPE for Rigel.**
* **A number measured at the panel's nascent fragment share (20.2 % capture-OFF and 2.6 % capture-ON) is a STRESS reading, not an
  expected-case reading.** Robustness may be judged at stress level; design decisions may not be
  driven from it. The over-rotation this ruling corrects: the tool accumulated design constraints from
  simulations carrying nascent at levels real data does not have, and started designing AROUND nascent
  instead of being robust TO it.
* ⚠ **STATE THE UNIT, BECAUSE THE MOLAR AND FRAGMENT NUMBERS DIFFER BY ~24× AND ONLY ONE OF THEM IS
  LARGE.** A nascent entity spans a whole gene (mean 40,667 bp on the ladder's index) while a mature
  transcript is spliced (mean 1,708 bp), so a molecularly SPARSE population still supplies a large
  share of FRAGMENTS. On the rebuilt panel the molar ratio is **0.00895** and the fragment share is
  **20.2 %** — the same population described two ways. ⛔ So "the panel is unrealistic" is precise only
  in fragment terms, and how unrealistic depends on the protocol: poly(A) selection strongly depletes
  unspliced RNA, total-RNA/rRNA-depleted protocols retain it. Say which unit a nascent number is in,
  every time, and expect INTRONIC fragment shares in real total-RNA libraries to be non-trivial even
  where the molar population is sparse.
* **Consequences licensed by this ruling:** a global gDNA background pooled over intergenic + intron
  slots is safe on real data — sparse nascent cannot move a pooled megabase-scale level. ⛔⛔ **BUT NOT ON
  THIS PANEL, AND THE NUMBER IS NOT A CONSTANT**: measured on the rebuilt ladder the intergenic+intron
  pool is inflated **1.18× at `g50` and 4.49× at `g05`** capture-OFF (it scales with the RNA:gDNA ratio),
  while intergenic-only is inflated **exactly 1.0000×** on every condition. So an intron-inclusive pool
  is licensed by the RULING for real data and is CONTAMINATED as a measurement here: the excess IS the
  intron nascent count, to the fragment. ⭐ And with strand-specific data most exons solve directly and
  intronic nascent is solved directly by strand — exon-imputation machinery earns its keep on UNSTRANDED
  data.
* **The symmetric honesty this ruling demands:** real intergenic regions carry unannotated
  transcription the tool does not model at all. The sim models nascent and not that — which is exactly
  why nascent got outsized attention: a modeled contaminant is worried about, an unmodeled one is
  glossed over. Neither background pool is "clean" on real data, and both are safe for the same
  reason — sparse RNA cannot move a pooled level.

⛔ **What this ruling is NOT:** it is not a license to break AXIOM 0 — unspliced RNA at an intron is
still RNA, and the three-population set is unchanged — and it does not delete robustness: the
synthetic nascent entity, the `--nrna` harness arms, and the zero controls all stay. It re-ranks
concerns; it does not remove the model.

### ⭐⭐⭐ WHAT "SPARSE" MEANS AS A MODEL — owner, 2026-08-22, and it is the simulator's `sparse` mode

**Sparsity is a per-gene-span ON/OFF pattern, not a low global level.** Transcription and degradation
set the equilibrium, so nascent RNA is absent from most spans and *present and measurable* in a
minority. The three rulings that define the mode (`sim.whole_genome.apply_sparse_nrna`, gated by
`tests/test_whole_genome_sim_config.py`):

* **The unit of sparsity is the GENE SPAN (the nascent entity), not the transcript.** Isoforms share a
  span, so a per-contributor draw would give a 5-isoform gene `1 − 0.9⁵ = 41 %` chance of nascent at
  `on_fraction = 0.1` — the intron slots calibration reads would be four times less sparse than
  configured.
* **The level is LOG-UNIFORM** over its range: where nascent is present it spans decades, and a linear
  draw on (1, 1000) puts 90 % of its mass in the top decade.
* ⭐⭐ **The level is drawn INDEPENDENTLY of the mature level, so `nascent > mature` is a real case the
  tool must survive.** The retired ratio modes set `nascent = mature × ratio`, which made the two
  perfectly rank-correlated and capped the interesting case. The steady state says independence is the
  better null: mature is synthesis/degradation and nascent is synthesis, so the nascent:mature ratio is
  a STABILITY parameter, not an expression one — an abundant stable transcript shows almost no nascent
  signal, an unstable rare one can show more nascent than mature.

⚠ **The fragment share is EMERGENT, and the length geometry dominates it**: an entity spans a whole gene
(mean 40,667 bp on the ladder's index) against a spliced 1,708 bp, so `E[level]` and `on_fraction`
together set the share through a 23.8× factor. Measured for mature ~ logU(1, 10⁴): level ~ logU(1, 100)
gives **4.2 %** at `on_fraction = 0.10` and **20.2 %** at 0.50 — measured on the ladder's own index and
then confirmed by the rebuilt panel to the fragment — while logU(1, 10³) gives ~24 % and ~61 %. ⛔ Price
the share from the parameters BEFORE simulating; a range chosen by eye lands decades off.

---

## 0c. ⭐⭐⭐ CALIBRATION IS MESSAGE PASSING, AND THE EXON IS WHY — owner, 2026-08-17

⛔⛔⛔ **THIS SECTION EXISTS BECAUSE THE DERIVATION BELOW HAS BEEN RE-DERIVED FROM SCRATCH IN THREE OR
FOUR SEPARATE SESSIONS, EACH TIME OVER MANY TURNS, EACH TIME ARRIVING AT THE SAME ENDPOINT.** Owner,
2026-08-17: *"we keep reinventing the wheel, and I keep having to convince the session agent that we are
just re-deriving message passing."* ⭐ **So it is SETTLED, and it is written to be READ rather than
re-derived.** If a design you are considering ends at *"then the exon's gDNA level has to come from its
neighbours"*, you have arrived here again — stop, and read this instead of deriving it.

### The structural fact

⭐⭐ **gDNA is DIRECTLY MEASURABLE only where no MATURE transcript crosses.** Three places, and the
predicate is the solver's own `mrna_active` — the same one §6b's boundary-axis ruling turns on:

* **intergenic REGIONs**,
* **intron REGIONs**,
* and the **BOUNDARIES against them** — `exon\|intron`, `exon\|intergenic`.

At an **EXON** the unspliced mass is `gDNA + unspliced RNA`, and that is **the quantity being solved
for**. An exon therefore has no observation of its own gDNA level, at any depth, in any library.

> ⛔⛔⛔ **AN EXON'S gDNA LEVEL CANNOT BE MEASURED. IT CAN ONLY BE IMPUTED FROM ITS NEIGHBOURS. AND
> IMPUTATION ACROSS THE CHAIN *IS* MESSAGE PASSING. THERE IS NO FOURTH OPTION.**

⭐ This is a statement about the ANNOTATION and the deposit rule, not about any estimator — so no better
estimator, no extra channel and no cleverer prior removes it. It is the reason the solver is a
belief-propagation sweep over the region chain (§1, stage 2) rather than a per-object fit.

### The object system — what each object contributes

```
  intron REGION  <->  intron|exon BOUNDARY  <->  exon REGION  [ <->  exon|intron BOUNDARY  <->  intron REGION ]
       measures                measures              THE TARGET
```

| object | what it contributes |
|---|---|
| **intron REGION** | deconvolve → **gDNA + unspliced RNA**. A direct gDNA observation, subject to the nascent level (§6b.1 ruling ①) |
| **`intron\|exon` BOUNDARY** | its UNSPLICED mass is **gDNA + unspliced RNA**; its **SPLICED** mass is **certified RNA** — gDNA cannot splice, so that arm needs no deconvolution at all |
| **exon REGION** | **gDNA + RNA — THE TARGET**, and unmeasurable alone |

⭐⭐ **Both flanking objects measure something the exon needs, and each measures a DIFFERENT one**: the
flank's gDNA density bounds the exon from BELOW, and the flank's **sj FLUX** measures the exon's **RNA**,
which bounds it from ABOVE (§0c.3). The exon is the only slot in that picture with nothing of its own.

### ⛔⛔ THE THREE ESCAPE HATCHES, ALL CLOSED BY MEASUREMENT (2026-08-17)

Sessions keep reaching for one of three ways around the solve. All three have been measured and all three
fail, and the cause is the same in each case: **hybrid-capture enrichment is PER EXON, nonlinear and
arbitrary** — it depends on which probes the panel happens to contain, which is not derivable from the
annotation, from the chain, or from any pooled statistic.

⭐ **The anchor ladder makes that concrete.** gDNA densities measured at `g98 ss0.99`:

| rung | capture OFF | capture ON |
|---|---|---|
| intergenic REGION | 0.100521 | 0.00418 |
| intron REGION | 0.100452 | 0.00422 |
| `exon\|intergenic` BOUNDARY | 0.098343 | 0.510 |
| `exon\|intron` BOUNDARY | 0.098391 | 0.478 |
| **span across the four rungs** | **1.0×** | **122×** |

⭐⭐ Off capture the four rungs are **one number**. On capture they are a LADDER **ordered by probe
proximity**, and rungs 3 and 4 are only **PARTIALLY** enriched — probes overlap exons, and a fragment
crossing an `exon\|intron` BOUNDARY only partially overlaps a probe (owner, 2026-08-17). ⛔ Which is why
those rungs **UNDER-READ a true exon by 2.6–3.6×**; §6b.1 records the same number as *"the in-gene anchor
is a DETECTOR, not a calibrated level"*.

| ⛔ the hatch | why it is closed |
|---|---|
| **① a POOLED scalar reference** — one library-wide gDNA level | **3.90× WORSE than `base` on stranded × capture-ON.** A scalar cannot express a per-exon enrichment that spans 122× within one condition |
| **② a LOCAL anchor — "just use the nearest measured rung"** | **1.27× worse (nearest rung); 1.50× worse (`density_model.region_gdna_density`, the shipped local boundary-anchored imputation)**, both on stranded × capture-ON. ⭐ And the SAME local form scores **0.4037–0.4977** — far BETTER than `base` — off capture. ⛔ It is not a bad estimator; it is a good estimator whose premise capture destroys. Partial enrichment is exactly the 2.6–3.6× under-read above |
| **③ `capture_eff_length`** | **REFUSED STRUCTURALLY, not by measurement.** `capture_eff_length.transcript_capture_eff_lengths(calibration: CalibrationResult, …)` **takes a solved system as its first argument**. It cannot be the thing that produces the solve; proposing it is a circular dependency, and the circularity is visible in the signature |

⛔ **So the escape hatches are not a matter of effort.** ① and ② fail because the quantity is genuinely
per-exon; ③ fails because it is downstream of the answer. What is left is the neighbours, and using the
neighbours is message passing.

### 0c.0 ⭐⭐⭐ A HOP CARRIES A **LEVEL** OR A **COMPOSITION**, AND THE TWO INVARIANCES ARE COMPLEMENTARY

⭐⭐ **This is the rule the message layer is being rebuilt around (owner, 2026-08-18/19).** The owner's
synopsis gives a currency per hop type; ONE reason sits under all of them, and stating it turns a table to
memorise into a rule to apply:

| currency | what it is | INVARIANT to | DESTROYED by |
|---|---|---|---|
| **LEVEL** `rho_g` | gDNA fragments per placement | **POPULATION** changes — TSS, TES, strand flips, splicing. gDNA is genomically continuous and knows nothing about transcripts | **ENRICHMENT** changes — a probe edge between the two objects |
| **COMPOSITION** `f_g` | the gDNA share of the crossing population | **ENRICHMENT** changes — a probe enriches everything overlapping it, both components alike, so the ratio survives | **POPULATION** changes — the denominator is a different set of molecules |

⭐ **They are complementary, so at every hop at least one is intact: use that one.** The licence question
*is* the currency question, asked once per hop type. No licence hierarchy, no three-case rule, no
per-strand patch.

⛔ **THE SHIPPED RELAY HAS THIS BACKWARDS, AND THAT IS WHY IT IS BEING REPLACED RATHER THAN REPAIRED.**
`RelayPolicy` transports a COMPOSITION by default (the reframe `r = rho_tot(dst)/rho_tot(src)`) and bolts
population licences on top — but the population is exactly what changes at the `exon|exon` boundaries that
carry the largest destination mass on the panel. Every licence bug of 2026-08-18 was an attempt to REPAIR a
currency CHOICE instead of changing it.

⭐⭐⭐ **AND THE POPULATION OF A MESSAGE IS DIRECTION-DEPENDENT** (owner, 2026-08-18) — *"always answer the
question: what crosses INTO this region?"*

| direction | what crosses into the destination | spliced fragments |
|---|---|---|
| BOUNDARY → EXON, i.e. **SPLICE IN** | unspliced boundary crossings (gDNA + nascent RNA) **plus the spliced fragments, which splice IN to this exon** | **INCLUDED** |
| EXON → BOUNDARY, i.e. **SPLICE OUT** | the unspliced fragments crossing contiguously | **EXCLUDED** — they splice OUT and land elsewhere |

⭐ So SPLICE IN and SPLICE OUT are not two bolted-on operators: they are **one rule — the message's
population is whatever physically enters the destination — evaluated in the two directions.** That is the
whole justification for the pair, and it is why both must exist.

⭐ **MEASURED, 2026-08-19 — the principle held on all 32 conditions, and the data named ONE thing the
synopsis did not** (`scripts/design/hop_currency.py`, Stage 2 of the rebuild; the numbers are
`ROADMAP.md` §0). A terminus — a TSS/TES, whether at a gene edge or lying inside another transcript's
intron or exon — is a POPULATION change and carries a LEVEL; a splice site into an exon carries the
SPLICE-IN COMPOSITION; an intron into its own boundary carries a COMPOSITION (exact to the fragment where
a LEVEL is off by 78–98 % under capture); and a hop OUT of an exon carries a LEVEL on every arm. ⛔ **So a
hop type is `object class × {splice site, terminus, both}` read off the boundary flags — the class alone
conflates a TSS/TES inside a gene with a splice site, and the two have OPPOSITE currencies**
(`TRAPS: an-object-class-does-not-see-a-terminus`). ⛔ And at a TERMINUS INTO AN EXON UNDER CAPTURE
neither currency survives — the enrichment changes AND the population changes — which is the residual
§0c.3's spike-and-slab exists for, now re-established on the rebuilt panel rather than inherited.

⭐⭐ **OWNER RULINGS, 2026-08-19, ON THE MAP AND HOW THE POLICY IS BUILT** (the arithmetic is
`EQUATIONS.md` §3.5e, worked in numbers):
* **the map MUST be per ARM — `{gDNA, RNA+, RNA−}`, three currencies measured, not one gDNA share with
  the strands pooled.** The 2026-08-19 map is the pooled-strand first pass and is to be extended before
  it is used.
* **a message is the three DENSITIES; a compatible hop RESCALES by the ratio of totals** (the sj flux
  counted in the boundary's total, removed at the source's scale, the rescale undone — SPLICE OUT; the
  inverse — SPLICE IN). **A terminus is where the rescale is not licensed** and only the level crosses;
  ⛔ a terminus AT a splice junction is still a terminus (`sj+term` rules with `term`, never with `sj`).
* ⛔ **messages DO flow into gDNA-measuring objects** — an exon sends its densities into its
  `exon|intron` boundary and the boundary does the SPLICE OUT arithmetic; "measuring objects are sources
  only" was proposed and is REFUSED.
* **the policy is built TOY-FIRST, one transcript at a time** (`docs/TESTING.md` §0b): a single-exon
  `TA+`, then a multi-exon `TB+` beside it, then more, the policy and the map evolving with each
  addition and the cached full scenarios run as the benchmark whenever wanted. ⛔ Not built straight
  from the map: the owner's judgement is that the logical complexity is not yet held firmly enough for
  that, and the record of this project says solving one thing while breaking another is the failure
  mode.

⛔ **THREE THINGS ARE SETTLED AND MUST NOT BE RE-ASKED**: forward-backward STAYS (an object may need a
message from BOTH neighbours; an ordered single pass delivers from one side only); the three arms
`{gDNA, RNA+, RNA−}` are carried NATIVELY (AXIOM 0 — collapsing the RNA pair would destroy the strand
channel, which is the primary intron deconvolution on stranded data); and the land/sea analogy the owner
used to explain the geometry is an EXPLANATION and never enters code or docs.

### 0c.0b ⭐⭐⭐ THE TWO STRATEGIES ARE ONE CONTINUUM, AND THE POINT ON IT IS FITTED — owner, 2026-08-20

⚠ **NO IMPLEMENTATION IN `src/` AND NO GATE since the 2026-08-27 tear-down.** The policy this ruling was built for was deleted; the reasoning is kept so it is not re-derived, and `EQUATIONS.md` §3.5f carries the arithmetic. It is NOT a description of shipped behaviour and NOT an instruction to rebuild it.

⛔ **§0c.0's table is a statement about INVARIANCES, not an instruction to pick one of two mechanisms**,
and reading it as the latter produced a policy with a strategy SWITCH whose measurement immediately
refuted the shape: the composition rescale helped under capture (0.61x) and hurt off it (**8.7x**
stranded), both correctly. ⭐ **The owner's reading, and it is the settled one:** *"absolute abundance
transfer and composition transfer are on two ends of a spectrum … there is actually a knob that connects
them. Absolute abundance transfer is possible when total abundance is uniform … composition transfer is
required when total abundance is non-uniform."*

⭐⭐ **THE KNOB IS DERIVED, NOT TUNED, AND IT HAS NO CONSTANT.** ABUNDANCE says the enrichment is 1;
COMPOSITION says it is exactly what the totals report. Those are two point hypotheses about ONE unknown,
the log enrichment. Fusing them by inverse variance — the observed `log r` with its counting variance `v`
against the no-enrichment premise, whose error if wrong is `log r` itself — is a shrinkage estimator:

    w = (log r)² / ((log r)² + v)        the claim crosses as   rho · r^w

so a disagreement indistinguishable from counting noise is shrunk away (the capture-OFF end) and one that
dwarfs it is believed in full (the capture-ON end), **per hop, chosen by the data**. ⭐ The residual
`(1 − w)·log r` is the abundance strategy's own transfer variance `((1−w)·log r)²`, which vanishes at the
composition end — so ONE expression spans the continuum and there is no switch to set.
The derivation is `EQUATIONS.md` §3.5f.

### 0c.0c ⭐⭐⭐ AN IMPUTATION IS A WEAK PREDICTOR BY DEFINITION, AND ITS PRECISION MUST SAY SO

⭐ **Owner's ruling, 2026-08-20**: *"message propagation … is a weak predictor by definition. It's an
imputation, not a measurement. The strand model is a measurement — we're using the counts at that object
to measure it … so message propagation needs to be weak."* ⛔ **The consequence is a design constraint,
not a preference: a regression on strand-specific data is a PRECISION defect until proven otherwise.**

⛔⛔ **AND THE GOAL IS NOT A SWITCH.** Strand specificity is a SPECTRUM — a library may be weakly or
strongly stranded — so the tool must not have cliff behaviour at some specificity threshold, and the
answer is never "turn message propagation off when stranded". The answer is honest precision, so that an
excellent stranded dataset drives pass-0 and is not weakened by imputation, and a poor one is still
carried. ⚠ The ladder simulates specificity as a BINARY (0.50 / 0.99) and gDNA at a handful of levels;
both are conveniences, and a mechanism that works only because the panel is coarse is not a mechanism.

⭐ The measurement that makes this checkable is a SPLIT, and it belongs on every message-layer
experiment: score the error separately over destinations that HAVE their own composition evidence and
those that do not. "The messages help" and "the messages trample a measurement" are different findings
that a pooled number cannot distinguish (`TRAPS: an-imputation-must-cost-something-every-hop`).

### 0c.0d ⭐⭐ WHAT PASS-0 IS ACTUALLY FOR — the circular bootstrap, stated so it is not re-litigated

⭐ **Owner, 2026-08-20.** Calibration is circular from the beginning: the tool does not know whether the
library is capture-enriched, nor how much gDNA it holds, and the reference prior is standing in for a
real prior it does not have yet. ⛔ **So "fix the reference prior" is NOT the goal, and treating it as
one forgets what pass-0 is for.** The goal is to solve or impute ENOUGH exons — not all of them — that
the gDNA landscape can be TRAINED on them; once trained, that per-object prior does the work, and it
subsumes much of what message propagation is doing.

⭐ Two things follow, and both are rulings rather than options:
* **Feed the data's background into ψ as a LIKELIHOOD wherever the data allow it** — intergenic space
  is pure gDNA in the annotation's terms (§0b: on real data it carries unmodelled transcription, which a
  POOLED level absorbs). ⛔ NOT as a reference location: the location concept was refuted and deleted
  (§6b.1, owner 2026-08-24) — the channel is the density λ-factor, whose precision scales with counts.
  ⛔ And do not extend it to objects whose composition the annotation does not determine.
* **Message propagation's irreducible job is the deep chain**: long runs of `exon|exon` boundaries where
  imputation is the only source of information there is. Strand specificity would solve those, and often
  is not available — which is exactly why the next best thing has to work, and has to be weak.

### 0c.1 ⭐⭐⭐ THE MECHANISM IS ALREADY BUILT AND IS NOW ON — do not build it again

⛔⛔ **The hop the derivation above asks for EXISTS IN `src/`, is individually switchable, and is behind
one config flag.** It is the **SPLICE IN** — `messages/relay.py`, switch `splice_in` on `RelayPolicy`:

> **SPLICE IN (BOUNDARY → EXON): the BOUNDARY's measured sj flux is a density AT THE SOURCE**, which joins the
> RNA claim entering the destination EXON.

⭐⭐ **Four properties, and each one is precisely what the exon problem needs:**

| the SPLICE IN | and why that is the property required |
|---|---|
| **only an EXON receives it** | an intron carries no sj flux, so there is nothing to SPLICE IN there — the operator is already scoped to the one slot class that has no observation of its own |
| **it is a MEASUREMENT (a COUNT), not an imputation** | so it carries **its own precision**, rather than inheriting the source's belief. This is the difference between an upper bound with a variance and a guess |
| **its transfer variance `hop_logvar` is identically 0** | on a matched-set SPLICE IN the reframe ratio `r` is common-mode and cancels, so the hop adds no scale-sampling variance of its own |
| ⭐⭐ **it is explicitly NOT tau-gated** | the source's PREDICTION precision is 0 on unstranded data, and a tau gate would drop the SPLICE IN on the floor there. ⛔ So the SPLICE IN **survives exactly the stratum where the strand channel is dead** — which is the stratum the exon problem is hardest on |

⚠ **VOCABULARY.** §0 names this operator **SPLICE IN**, and `src/` now spells it that way
(`splice_in`, `splice_in_frame_logvar`, `splice_in_premise_logvar`). ⚠ `graft` was the pre-ruling identifier
and survives only in commit messages and in the history recorded here — a grep for it finds this note.

⛔ **Where it is switched off:** `CalibrationConfig.message_propagation = False` installs
`messages/silent.py`'s `SilentPolicy`, which sends nothing (§6.1). Turning the relay on is one flag; every
operator inside `RelayPolicy` remains behind its own named switch, so the SPLICE IN can be priced alone.

### 0c.2 ⚠ HISTORY — WHY IT *WAS* OFF, AND WHY THAT TABLE MAY NOT BE QUOTED

⛔⛔ **THE RELAY IS ON. `CalibrationConfig.message_propagation = True` (owner, 2026-08-18) and
`calibrate.py` installs `RelayPolicy()`.** Everything in this subsection is the record of the ~11-day
mute, kept because its NUMBERS are still quoted at people; it is not the current state. ⚠ And the
policy it describes is itself being replaced — §0c.0 is the ruling, and the relay is REBUILT rather
than repaired.

⭐ **The mute is a MEASUREMENT and it stands.** Measured 2026-08-07 (`config.py`'s
`message_propagation` is the home of these digits) — ⚠ **on the 36-condition ladder, RETIRED and rebuilt
at 16 conditions on 2026-08-13, so they stand as recorded and are not reproducible as written**:

| stratum | muting the relay | rows better |
|---|---|---|
| stranded × capture-ON | **−58.3 %** | 16/16 |
| stranded × capture-OFF | **−43.7 %** | 16/16 |
| unstranded × capture-OFF | **−32.1 %** | 14/16 |
| ⛔ unstranded × capture-ON | **+154.8 %** | 0/16 |

⭐⭐ **The three the mute IMPROVES are exactly 0.8.0's three IN-SCOPE strata and the one it hurts is the
DEFERRED one**, so the shipped `off` is aligned with §0b rather than in tension with it. ⚠ The `/16` is
SCORED ROWS, not conditions (8 scored conditions × 2 axes); on the 16-condition ladder the same arithmetic
gives **`n/6`**, and that is a smaller panel rather than a coverage regression.

⛔⛔⛔ **BUT THE +154.8 % IS NOT EVIDENCE ABOUT THE RELAY AS IT WOULD BE FIXED, AND INHERITING IT IS THE
ERROR THIS RULING EXISTS TO PREVENT.** It was measured **with a confirmed bug live**:

> `RelayPolicy`'s **composition licence** implements §4.1's rule — *a composition may be imputed across a
> step iff the source SUPPLIED both components AND the two objects measure the same RNA POPULATION* — for
> transcript **TERMINI ONLY** (`region_geometry.terminus_flank_gain`). ⛔ **`mrna_active` flipping across
> a hop is a population change too, and it is NOT checked.** So a **CORRECT** pure-gDNA claim at an intron
> is relayed into the adjacent EXON and drives it to a confident **wrong vertex**.

⭐⭐ **That is a BUG, not a verdict on message passing.** The licence RULE is right and settled (§4.1);
the predicate set implementing it is incomplete by exactly one predicate. ⛔ It was recorded as a strict
xfail in `tests/calibration/test_structural_reference.py`; that file died with the reference-location
deletion (§6b.1, 2026-08-24), so THIS paragraph is now the record — and the defect is localised **to the λ-message**: the same
intron goes **0.9006 → 0.7661**, the wrong way, and nulling `lam_channel` restores **0.9006** exactly while
`cm_g` / `cm_p` stay 0. ⛔ **That localisation is NOT a licence to close it by tuning the reference**: the
xfail records that softening the prior was MEASURED WORSE and that topping up `τ_λ` is refused outright.
The repair is the missing predicate, not a strength.

⛔⛔ **THEREFORE THE PRICE MUST BE RE-PRICED AND MAY NOT BE INHERITED.** Two things differ between the
recorded number and any future one — **the panel** (36 conditions, retired) and **the bug** (live at the
time) — so the +154.8 % is a measurement of a defective relay on a panel that no longer exists. ⛔ A
session that quotes it as *"message propagation was tried and cost +155 %"* has stated something the
measurement does not support. ⭐ The honest form is: **the relay has never been priced with the licence
extended to `mrna_active`, on the current panel.**

⚠ **§6b.1 records the same defect from the reference's side** (*"AND IT EXPOSED A RELAY DEFECT IT DOES NOT
CAUSE"*) and that paragraph is the measurement's home; this is the ruling. ⚠ Two of the suite's xfails are
the recorded price of the switch rather than defects — they go green the instant the flag flips — so the
xfail list must not be read as uniform.

### 0c.3 ⭐⭐⭐ THE SHAPE OF THE REFERENCE UNDER CAPTURE IS SPIKE-AND-SLAB

⭐ **A ruling on the SHAPE. The maths is `EQUATIONS.md` §9d.1–§9d.4; what is settled here is what the
reference has to look like once capture is admitted.**

⛔⛔ **AND THE TWO DOCS SHARE THESE TABLES ON PURPOSE, SO THE HOMES ARE NAMED RATHER THAN LEFT TO
DIVERGE** (`CLAUDE.md`'s MOVE rule): **`EQUATIONS.md` §9d owns the DERIVATION and the DIGITS** — §9d.1 the
structural fact, §9d.2 the anchor ladder and the two failed imputations, §9d.3 the two bounds off §3b's
conserved identity, §9d.4 the mixture and its three limits — and **this section owns the RULING**: that
the shape is settled, that the mechanism already exists (§0c.1), and that the price of turning it on may
not be inherited (§0c.2). ⭐ **Re-measure into §9d and cite it from here; a number edited in one place and
not the other is the failure this note exists to catch.**

Write `rho_g,i = rho_0 · eps_i`, with `rho_0` the un-enriched density and `eps_i ≥ 1` the enrichment at
object `i`. ⭐⭐ **A probe panel makes the distribution of `eps` across exons a SPIKE at `eps = 1`
(unprobed) plus a SLAB at high `eps` (probed) — that is the PHYSICAL STRUCTURE OF A CAPTURE PANEL, not an
analogy.** The 122× ladder above is the same fact seen as a marginal.

So the reference on `log rho_g` is a two-component mixture, and every part of it is set by something
already measured or already ruled:

| part | what sets it |
|---|---|
| **SPIKE** | `rho_0` and its width, measured **EXACTLY** from the **OFF-TARGET anchors** — intergenic + intron REGIONs, §6b.1's calibrated level |
| **SLAB — lower endpoint** | the **MONOTONE bound from the adjacent BOUNDARY**: enrichment is monotone in probe proximity, so a flank's anchor density is a floor for the exon beside it |
| **SLAB — upper endpoint** | the object's **OWN total density** — it cannot hold more gDNA than the mass it holds |
| **the mixing weight** | the **unprobed fraction**. At pass-0 nothing has been observed about the probe indicator, so it is the reference's OWN symmetric Jeffreys exponents `a = b = _JEFFREYS_REF` with no observation added — mean and median **½**, which is `simplex_logodds._NEUTRAL_LOCATION` exactly (`EQUATIONS.md` §9c.1). ⛔ **NO NEW CONSTANT IS INTRODUCED** (`TRAPS: no-magic-numbers`). ⚠ **This is NOT §6b.1 ruling ②'s 0.75** — that is the same exponents *plus one pseudo-observation OF gDNA*, `(a+1)/(a+b+1)`, and it is a claim about a DIFFERENT quantity (ψ's location where `¬mrna_active`). Same convention, different posterior; do not reconcile the two digits |

⭐⭐ **IT DEGENERATES TO THE SHIPPED FORM WHERE ITS EXTRA ASSUMPTION IS INERT**, which is the property
§6b.1 demands of any reference. Three limits, each checked:

* **capture-OFF** — the measured anchor ratio is **0.98** (§6b.1), so the slab **collapses onto the
  spike** and the mixture degenerates to the single measured density: exactly the capture-OFF reference
  already validated.
* **`g00`** — the anchors measure `rho_0 = 0`, so the spike and the slab's lower endpoint are both 0 and
  `f_g → 0`.
* **`g98` capture-ON** — the slab's upper endpoint is **0.9918** against a truth of **0.9817**.

#### The two bounds, and what they ARE

⭐ **Both endpoints come from NEIGHBOURS, which is §0c's structural fact restated in arithmetic:**

    lower:  f_g >= rho_anchor · E_g / M          the adjacent exon|intron flank; enrichment is MONOTONE in probe proximity
    upper:  f_g <= 1 − (rho_r · E_r − S) / M     rho_r from the adjacent BOUNDARY's sj FLUX; S = the certified spliced count

using the conserved identity `rho_r · E_r = unspliced RNA + S`, which `EQUATIONS.md` §3b already carries.

⛔⛔ **THEY ARE A SUPPORT CONSTRAINT, NOT AN ESTIMATE, and reading a bound as an estimate is
`TRAPS: an-upper-bound-is-not-an-estimate`.** Measured mass-weighted over exon REGIONs against the
origin-split oracle, upper bound only:

| condition | true `f_g` | upper bound | mass in violation |
|---|---|---|---|
| `g00 ss0.50` off | 0.0000 | 0.6039 | **0.0 %** |
| `g50 ss0.50` off | 0.0627 | 0.6150 | 5.6 % |
| `g98 ss0.99` off | 0.7672 | 0.8332 | 19.5 % |
| `g50 ss0.99` ON | 0.5220 | 0.8272 | 9.8 % |
| ⭐ `g98 ss0.99` ON | 0.9817 | **0.9918** | 4.6 % |

⭐ **The bound is LOOSE where RNA is abundant and TIGHT where RNA is scarce.** On the last row it sits
**0.0170** from the truth, against **0.4817** for the neutral ½ — on the stratum where a per-exon gDNA
level is otherwise unobtainable. ⭐⭐ **The two bounds are COMPLEMENTARY**: the lower one is tight at LOW
`f_g` and fails under capture; the upper one is tight at HIGH `f_g` under capture. That is why the slab is
bracketed by BOTH and not by either.

#### ⭐⭐⭐ AND IT ESCAPES THE VERTEX THEOREM — which is a real result, not a convenience

`EQUATIONS.md` §9a proves that a simplex vertex is unreachable without evidence: *every **PROPER** prior
on `[0,1]` has a median strictly **INSIDE** `(0,1)`*. ⭐ **That argument needs the prior to have a
DENSITY** — it requires the CDF continuous and strictly increasing. **A SPIKE-AND-SLAB HAS AN ATOM**, so
its median sits **exactly at the spike** whenever the spike carries at least half the mass.

⛔ **The theorem is not violated; its PREMISE does not hold.** A Beta cannot reach `f_g = 0`; a
spike-and-slab can.

⛔⛔ **THE §9a THEOREM ITSELF IS UNTOUCHED AND MUST BE KEPT — it is correct FOR PRIORS WITH A DENSITY.**
What changes is its **SCOPE**. The recorded reading that the vertex shortfall is *the value of MISSING
INFORMATION rather than headroom* — priced by `scripts/design/vertex_ceiling.py`, whose re-solve arm reads
mean `|error|` 0.0538 → 0.0407, **−24.4 %** — is a property of the prior **FAMILY**, not a limit of the
data. ⭐ **Add the scope; do not delete the theorem, and do not re-rank the shortfall as a defect on the
strength of this alone.**

#### ⭐⭐⭐ RECONCILING THIS WITH §6b's "RNA IS THE RESIDUAL AND IS NEVER PREDICTED" — by SCOPE

⛔ **A session WILL trip on this, so it is settled here.** §6b's ruling stands exactly as written, and it
is a ruling about the **OFF-TARGET** case: gDNA is near-uniform over the genome and measurable before any
solve, RNA spans six orders of magnitude with essentially no genomic autocorrelation, so a **pooled** RNA
density is not a population parameter and pooling sj flux across the genome is inadmissible. ⭐⭐ **At an
EXON UNDER CAPTURE the two roles are EXCHANGED**: capture makes gDNA per-exon, nonlinear and arbitrary —
the 122× ladder — so **gDNA is the unpredictable one there**, while **RNA is MEASURED, locally, by the
adjacent BOUNDARY's sj FLUX**. ⛔ That is not a pooled RNA density and never becomes one: the flux is that
exon's own neighbour's COUNT, and pooling flux across the genome remains inadmissible. **Which component
is the residual is decided by which one THE DATA MEASURES AT THAT OBJECT — and capture moves that line.**

---

## 1. The shape of the tool

Three stages, `pipeline.py`:

1. **Scan** (`scan_and_buffer`) — a C++ htslib single-pass reader. Resolves fragments against the index,
   trains the strand and fragment-length models from unique mappers, buffers fragments into a columnar
   `FragmentBuffer`, and deposits per-object tallies into the C++ **accumulator** → `AccumulatorPayload`.
2. **Second pass** (`_drain_side_buffer`) — drains the fragments pass 1 held, then
   **Calibrate** (`calibration.calibrate`) — deconvolves each object's unspliced mass into gDNA vs RNA by
   a belief-propagation sweep over the region chain, and fits the library-level parameters.
3. **Quantify** (`quant_from_buffer`) — scores fragments, builds loci by connected components, runs a
   per-locus EM with `n_transcripts + 1` components. The calibration prior enters as **two per-locus
   Dirichlet scalars**, never per-transcript.

⚠ **Calibration cost is depth-INDEPENDENT** — every region in the index is solved regardless of read depth,
so a 971 k-fragment targeted BAM pays full genome-scale cost. On one real run: index load ~7 s, BAM scan
~2 s, calibration ~66 s, per-locus EM ~24 s. **The scan is ~2 % of runtime — accumulator work is
essentially free and calibration is the whole budget.**

---

## 2. The index

`INDEX_FORMAT_VERSION 8`, shipped as `regions.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`.

> The genome is a graph. **Regions** are intervals; **boundaries** connect them. A fragment is a **path**.
> Regions count fragments *contained*; boundaries count fragments *crossing* (a 0-bp boundary, no width).

- **Regions** tile each reference, cut at **every exon endpoint** of every non-synthetic transcript, with
  **no merging** (`EQUATIONS.md` §1.7). 1 bp regions are legal and common — nothing may assume length > 1.
- **Two boundary kinds.** *Contiguous* boundaries are the boundary between genomically adjacent regions: bidirectional,
  carrying gDNA + RNA, endpoints **implicit** (boundary `i` sits between region `i` and `i+1`). *Junction* boundaries
  are directed `src → dst` in GENOMIC order (§0 — never donor→acceptor), **pure mature RNA** by
  construction, need explicit `(src, dst, strand)`,
  and carry no unspliced channel and no structural flags.
- A splice jump deposits on its **sj boundary only** — never on the contiguous boundaries it splices over.
- Boundaries always run `src < dst`, so **genomic order is a topological order and there is no graph traversal
  anywhere.** The graph is a DAG but **not** a polytree (every sj boundary closes one undirected loop),
  so sj boundaries must be *factors on their endpoint regions*, never message channels (TRAPS: splicing-makes-the-graph-cyclic).
- **8 structural flag bits per contiguous boundary**: TSS / TES / DONOR / ACCEPTOR × {+,−}, not mutually
  exclusive. Carry the raw bits to the consumer; do not pre-derive predicates in the plumbing.
- Validated by invariants I1–I13, two of which re-derive the answer by a **different algorithm**
  (TRAPS: self-checking-validator).
- `manifest.json` records the sources, their sha256 and the build flags. The build is deterministic.

⛔ **Never quote a stored census.** Region and boundary counts are properties of an *annotation*, not of the
tool: run `scripts/design/index_census.py`.

⚠ **`reach` is covered by no existing hash.** A rebuild moved 38 % of human reaches with both
`partition_hash` and `graph_hash` byte-identical, so any calibration-facing cache needs a third key
(TRAPS: a-hash-that-misses-its-artifact).

---

## 3. The accumulator

`tests/native/_accumulator_reference.py` **is the executable specification.** The C++ is gated on
byte-identity to it; where it and a document disagree, the reference wins.

### 3.1 What every object stores

⭐⭐ **A CHANNEL IS STORED WHERE A NAMED CONSUMER READS IT, AND NOWHERE ELSE** (owner-set, 2026-08-08).
The populations therefore do NOT all carry the same channels, and that asymmetry is the design:

    region_contained    count  inv_opportunity_sum
    region (per path)   start_count[2]  end_count[2]  span_count[2]   ← 2026-08-21, see below
    boundary_unspliced  count  inv_length_sum       mass
    boundary_spliced    count                       mass      certified RNA — nothing deconvolves it
    sj                  count  inv_length_sum       mass[2]   inv_length_sum is LIVE in second_pass

⭐⭐ **THE START/END/SPAN REGION BANKS (owner's taxonomy, 2026-08-21).** Every accepted path books its
FIRST covered base (`region_start_count`), its LAST (`region_end_count`), and every region it STRICTLY
spans (`region_span_count` — every base covered by one segment, neither endpoint inside; a region a
spliced path JUMPS over books nothing), each by genome strand. START and END have opportunity `ℓ` for
EVERY fragment length — the composition-free TOTAL at a REGION — and are wall-blind only at the
template's downstream/upstream end respectively, so a consumer side-selects ("use the side whose wall
does not bind"). SPAN's opportunity is `(w−ℓ−1)₊`, a per-component pmf functional BY DESIGN — it is
consumer-gated (the parked length-solve channel). ⭐ The ledger closes TWICE over:
`ΣS = ΣE = qc.deposited`; per region `contained ≤ min(S, E)`; `span ≡ 0` wherever `ℓ ≥ w_max − 1`.
⚠ `region_span_count` is NOT the removed `region_spanning_*` family: those carried the mis-weighted
`1/L` deposit; this is an integer count whose divisor lives with its (gated, future) consumer.

#### 3.1a-i ⭐⭐⭐ THE FOUR FORMS AN ABUNDANCE CAN TAKE — the classification rule, ruled 2026-08-21

A census of 240 sites found that the tree forms a density four ways, and **only two of them are
defects**. This table is the rule that decides whether a given repair applies to a given site; it was
written because a repair was proposed for a site the rule excludes.

| form | expectation | verdict |
|---|---|---|
| the REGION contained reciprocal bank, `Σ 1/(ℓ−w+1)` | `ρ·P(w ≤ ℓ)` — **TRUNCATED** | model-free but biased; §2's defect. 11.6× at a 98 bp exon |
| `count / E_contained(ℓ, pmf)` | `ρ` — **unbiased** | correct, but the fragment-length pmf enters the DIVISOR, so a distorted pmf distorts the level |
| the BOUNDARY/sj reciprocal banks, `Σ 1/(w−1)` | `ρ` exactly, every library | ⭐ already model-free AND unbiased — leave them alone |
| `mass / E_c` with numerator and divisor COMPONENT-MATCHED | `ρ_c` | ⛔ **NOT an abundance site at all.** A deconvolved gDNA mass over the gDNA opportunity is correct |

⭐⭐ **AND THE SCOPE RULE THAT FOLLOWS (owner, 2026-08-21): a TOTAL abundance is a PRE-solve instrument.**
It is needed in exactly three places — ① the density model (the gDNA background and the intron
deconvolution), ② the total-abundance LANDSCAPE, ③ ENRICHMENT RATIOS, which are ratios of two totals.
⛔ **Everything POST-calibration is out of scope by construction**: once beliefs exist, each component
has its own fragment-length distribution and its own opportunity, so a per-component estimate is the
right instrument and a total is the wrong one. ⚠ `capture_eff_length` and `priors` are the sites that
tested this rule — they consume a solved `CalibrationResult`, so their `mass/E_g` is component-matched
and a total does not belong there, however wrong their number may be for other reasons.

#### 3.1a-ii ⭐⭐⭐ THE WALL RULE AND THE SIDE SELECTION — the consumer half, ruled 2026-08-21

`rigel/calibration/total_abundance.py` turns those banks into a per-slot TOTAL, and four things about
it are settled so they are not re-litigated.

**① A side is EXACT iff its template distance clears `w_max − 1`, and `w_max` is READ.** The exact
start opportunity is `A_start(w | d) = min(ℓ, (d + ℓ − w + 1)₊)`, which equals `ℓ` for every `w` iff
the template continues `w_max − 1` bases past the region's genomic-HIGH bound; the END bank mirrors at
the LOW bound. `w_max` comes from the support end of `deposited_lengths` — never a quantile, because a
quantile lets a conservative-looking choice mark a binding wall exact. At `d = 0` the start rule
degenerates into the CONTAINED rule while still depositing the flat `1/ℓ`, which is why the pair exists
rather than a taper on the divisor.

**② The distance is the COMPONENT MINIMUM over the populations AXIOM 0 admits at that slot.** gDNA is
always admitted and its template is the contig; RNA on strand `s` only where `free_s`, taking the
SPLICED distance where an exon covers the region and the genomic contiguous-boundary reach where none
does. ⛔ The mature arm must be able to bind BELOW a long nascent reach — that is the measured
component-differential wall, and a genomic distance at an exon marks a binding wall exact.

**③ A DOUBLE-WALLED slot is honestly NOT MODEL-FREE and reads NaN.** No deposit rule is model-free
where both sides bind; the mask says so and the total refuses to invent a number, rather than papering
it over with an average of two depressed estimates. Measured coverage (`1 − double-walled`, START-mass
weighted, 16-condition ladder): **94.7 %** at capture-OFF, **84.3 %** under capture.

**④ The population for the MATURE distances EXCLUDES synthetic spans, and the filter is written.** A
synthetic entity is a manufactured nascent template and a nascent molecule extends genomically, so it
belongs to the reach arm. ⚠ On every shipped index the EXON rows happen to contain no synthetic
transcript, so this was once true by accident and cost 57 diverging distances between two
implementations (`TRAPS: state-the-population-rule-do-not-inherit-it-from-a-table`).

#### 3.1a-iii ⭐⭐⭐ WHICH OF THE LANDSCAPE'S OUTPUTS A CONSUMER MAY READ — ruled 2026-08-21, from the grid sweep

`abundance_landscape.AbundanceLandscape` publishes `rho_0`, `span_R`, `w_slot`, the mode list and an
anchor verdict. ⛔ **They are NOT equally trustworthy, and the difference is measured rather than
assumed** (`landscape_head_to_head.py --grid-sweep`, 16 conditions, `_N_GRID` swept over a 16× range):

**① `rho_0` AND THE ANCHOR VERDICT ARE CONSUMABLE.** `rho_0` moves 8–25 % across that whole range and
the anchor-consistency verdict holds **12/12 on every contaminated row at every grid**. Against certified
gDNA truth (the `gdna` partition's pooled intergenic rate through the identical side selection) the
depleted mode reads **0.0056 / 0.0264 nats** on the two capture-OFF strata.

**② `span_R` IS NOT CONSUMABLE AS IT STANDS.** On `g50 ss0.99` capture-OFF it reads 58 → 77 → 95.6 →
94.7 → **1.9** as the grid refines. `split_basins` selects the enriched mode by BASIN MASS, so
over-resolution fragments the bulk into sub-bumps and a nearby sub-bump of a huge bulk outweighs the
distant exon mode. ⭐ A consumer that needs the span must either pin the resolution or select the
enriched mode by something other than mass — and `TRAPS: a-mode-count-is-not-a-well-posed-quantity` is
why the mode COUNT may not be read at all.

**③ THE ESTIMAND IS WHAT MAKES THIS LANDSCAPE BEAT THE NPMLE, NOT THE ESTIMATOR.** Both are the same
family of fit; the NPMLE is fitted on `mass / eff_gdna` — a total over one component's OPPORTUNITY
MODEL — so its axis carries the divisor's per-region spread, measured as an offset IQR of **0.12 nats
off capture and 1.66 under it**, comparable to the mode separation itself and removable by no bandwidth.
The landscape's divisor is a geometry. ⚠ On generic held-out predictive likelihood the two TIE off
capture (−5.36 vs −5.35) and the landscape wins by 0.47 nats under it: the tie says the NPMLE is not a
bad density estimate, and the nats column says it is the wrong quantity.

⛔⛔ **AND THIS MODULE IS NOT `ROADMAP.md` §4.3's REFUSED DROP-IN — the distinction matters because the
tokens look alike.** §4.3 refused `region_start_count / ℓ` substituted for `RegionGeometry.inv_abundance`
INSIDE the currency channel, with no wall rule, no side selection and no model-free flag; the knob was
deleted after pricing. This module is a separate per-slot quantity with all three, wired to NO consumer,
and it does not touch `inv_abundance`. §4.3's bar still stands for any consumer swap: re-opening needs a
policy whose DELIVERABLE improves under a better level channel.

| | | |
|---|---|---|
| `count` | `Σ 1` | statistical power — a count is a count |
| `inv_length_sum` / `inv_opportunity_sum` | `Σ 1/A(w)` (float64) | ⭐ TWO deposit rules under two names, each cancelling its own opportunity **on its own support** (`E[Σ] = ρ·P(A>0)`, `EQUATIONS.md` §2): the BOUNDARY/sj form `1/(w−1)` has `P(w≥2) = 1` on any real library — an exact model-free density — while the REGION form `1/(ℓ−w+1)` reads `ρ·P(w≤ℓ)`, a density SHAPE truncated by a per-component pmf functional (`TRAPS: a-cancellation-is-conditional-on-its-support`). Not called `density` because at a REGION it is not one |
| `mass` | `Σ (slice/L)/n_bounds` (float64) | ⭐ the CONSERVED fragment count — sums to **one per fragment**, where `count` is `+1` on each of `max(K,1)` objects. ⭐⭐ A SJ BOUNDARY is a boundary exactly like a contiguous one, so a spliced fragment shares its one unit across every object it crosses, sj included. `EQUATIONS.md` §3b |

⛔ **Six banks were REMOVED on that rule** (2026-08-08): three `region_spanning_*`, the two
spliced-boundary length moments, and `sj_length_sum`. Structs went `Region` 80 → **24 B**, `Boundary`
80 → **48 B**, `SpliceJunction` 40 → **16 B**. ⚠ Deleting `region_spanning` has one structural consequence
worth knowing: **no spliced fragment touches the region axis at all** now — a spliced fragment can never be
*contained*, because both endpoints of an annotated intron are region bounds.

⛔⛔ **AND TWO MORE WENT ON 2026-08-13 — `region_contained_length_sum` and `boundary_unspliced_length_sum`
— BECAUSE THEIR STATED JUSTIFICATION WAS FALSE, WHICH IS WHAT MADE THE DELETION CLEAN RATHER THAN LOSSY.**
This section used to carry a `length_sum  Σ L` row claiming it *"carries the only information about the
gDNA/RNA split when the two components share a mean length"*. It does not: at `mu_g = mu_r = mu` the third
row is `(mu, mu)`, proportional to `(1, 1)` exactly like the other two, so the 3×2 system is still rank
one and `length_sum` removes nothing. It is an independent tilt **only when the means already differ** —
precisely when the first two rows already suffice (`TRAPS: equal-lengths-carry-no-composition`).
⚠ **The claim survived because nothing consumed it**: the banks reached `PopulationView.length_sum` and
stopped, so no test could disagree — a justification with no consumer is a claim with no gate. ⭐ The
retraction is written where the claim was, in `scan_payload.py`'s module docstring, not merely deleted.

⛔ **THE COUNTS KEEP BOTH GENOME-STRAND COLUMNS; THE INV-LENGTH SUMS AND THE BOUNDARY MASSES KEEP ONE —
AND `sj_mass` IS THE ONE REVERSAL OF THAT, RECORDED BELOW.** Which strand
a read aligned to says nothing about whether the molecule was gDNA or RNA, and every consumer summed the
two columns — so a strand axis there is half the bank wasted, and worse, a claim that the split is
meaningful. The counts keep both because the strand model is a Beta-Binomial over them, per strand.

⛔⛔ **AND `sj_count` KEEPS BOTH FOR A REASON THAT IS NOT THE STRAND MODEL** (owner, 2026-08-08). A
sj is stranded by its genomic splicing MOTIF, so the *fragments'* strand looks redundant, and every
consumer today sums the two. It is retained for **ALIGNER-ARTIFACT DETECTION**: aligners emit
false-positive `N` CIGAR ops from plain genomic DNA, `rigel.splice_blacklist` catches only those
`alignable` has enumerated by coordinate, and the empirical detector is that in a **stranded** library a
real sj inherits the global strand specificity while an artifact deposits onto BOTH strands.
⚠ Unstranded data cannot use it (κ = ½ leaves nothing to deviate from) — a property of the detector, not
a reason to drop the column. ⭐ The discriminating information lives ONLY in the split: a clean sj
and an artifactual one at the same depth carry the same total. Gated by
`test_the_junction_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION`.

⛔⛔⛔ **AND THAT SAME PREMISE REVERSES THE ONE-VALUE MASS RULING — ON THE SJ AXIS ALONE** (owner,
2026-08-12; landed 2026-08-13). `sj_mass` is `mass[2]`, and it is the only mass with a strand column.
⭐⭐ **The reversal is admissible because the PREMISE CHANGED, and the premise is recorded WITH the ruling
so that neither direction is re-litigated.** The original ruling was *"nothing reads a mass per strand"*.
That is now false for sj and only for sj: an **artifactual** splice junction accumulates SYMMETRICALLY on
both strands, exactly as gDNA does, so the strand model the tool already has can detect one — but only if it
is given a per-strand OBSERVABLE, and the count alone is not enough, because a count cannot separate a sj
used by many short fragments from one used by few long ones. ⚠ The second reason is structural: without the
bank, artifact filtering needs TWO passes over the BAM (tally, filter, re-accumulate the mass), which is the
one thing the single-pass architecture exists to avoid.
⛔ **The column is `col` — the SAME genome-strand column the count is deposited at** — so `mass[c]/count[c]`
is a per-strand mean and not a ratio of two different populations. ⚠ Summed over columns it is
byte-*comparable* to the single accumulator it replaced but **not bit-identical**: float addition is not
associative and the per-column deposit order differs (~1e-15 relative).
⛔ **SCOPE IT DELIBERATELY: the reversal is for the SJ axis, because that is where the consumer is.**
Whether the BOUNDARY axis's `unspliced_mass` / `spliced_mass` should also go per-strand is a SEPARATE
question with **no consumer named**, and `TRAPS: one-thing-varied` says do not bundle it — the one-value
ruling still STANDS there.
⭐ **`substrate` folds the strand axis at its own boundary**, so `PopulationView.mass` stays
strand-agnostic and no consumer below it changed at all. The per-strand values are deliberately not
re-exported there: their consumer is artifact filtering, which reads the payload.
⭐ **The memory ledger of the bundle, and it goes DOWN.** `SpliceJunction` 16 → **32 B** for the second
mass, against `Region` 24 → **16 B** and `Boundary` 48 → **40 B** for the two `length_sum` banks deleted
with it — **net ~19 MB lower at genome scale**. `sj_mass` moved `SINGLE_COLUMN_AXES` → `BANK_AXES`, and
since the shape digest is DERIVED from those two tables it moved itself.

⛔ **THE FIXED POINT IS GONE** (owner ruling, 2026-08-11: one numeric convention, and mixing two in one
schema is worse than either). **A COUNT is an integer; a FRACTION is float64.** There is no scale
constant and nothing decodes a bank. Headroom is no longer a question anyone has to answer.
⭐ It is also the more ACCURATE choice, measured against exact rational arithmetic on the
reciprocal-opportunity theorem: the fixed point missed the answer by 7.0e-10 … 2.0e-07 where float64
misses by 5.8e-15 … 2.8e-13 — 1e5–7e5× closer. Memory is unchanged **by this ruling** and still flat —
~85 MB at human scale on the struct sizes of the day, no hash map. ⚠ The 2026-08-13 bundle above took a
further ~19 MB off that; the current sizes are `Region` **16 B**, `Boundary` **40 B**, `SpliceJunction`
**32 B**, and the `static_assert`s beside them are the only place worth reading them from.

⭐⭐ **AND ONE OF THEM IS ALREADY LOAD-BEARING IN A WAY WORTH STATING**: `boundary_unspliced_inv_length_sum` is LIVE in `second_pass` (via `pipeline`). ⚠ It is NOT "the one channel whose opportunity and deposit cancel" — there are THREE reciprocal-opportunity channels (`sj_inv_length_sum` takes the identical `1/(w−1)` deposit at `accumulator.cpp:711`, and the REGION bank cancels within its own support) — but the two BOUNDARY-form ones are the UNCONDITIONALLY robust ones: `E[Σ] = rho·P(w≥2) = rho` for any real library, any length gap, where the REGION form reads `rho·P(w≤ℓ)` (`EQUATIONS.md` §2, §3c; `TRAPS: a-cancellation-is-conditional-on-its-support`).

⛔ **THE TALLY IS NOT BIT-REPRODUCIBLE ACROSS WORKER COUNTS, and the owner has signed that off**
(2026-08-11). Integer addition is associative, so every COUNT bank still reproduces exactly — measured
0.000e+00 spread over 8 shipped multi-threaded runs. Float addition is not, so the FRACTION banks are
re-associated by the per-worker merge and wander by ~1e-15; that reaches `posterior_mean` at **1.503e-15**
and the deliverable at ~1e-11, five orders below `EMConfig.convergence_delta`.
⭐ **Tests validate the float banks within a DERIVED tolerance, bracketed from both sides** — a nudge just
past it must fail and one just inside must pass, so it is a tolerance rather than a threshold that never
binds. ⚠ This is a SECOND source of irreproducibility, independent of `EMConfig.seed`
(TRAPS: the-deliverable-is-not-reproducible-by-default); pinning the seed alone is no longer sufficient
for byte-identical output. `TRAPS: integer-channels-reproduce` carries the numbers.

### 3.1b ⭐⭐⭐ WHO OWNS A FRAGMENT — and nothing is ever re-attributed

**Owner ruling, 2026-08-08.** It was already true of the accumulator, the `CalibrationResult` and the
sweep; it is stated here because ONE consumer had drifted from it for a year.

> **A REGION owns the fragments CONTAINED in it. An BOUNDARY owns the fragments that CROSS it. No object's
> mass is ever moved onto another object.**

⭐ **A locus therefore collects BOTH kinds of object**: its REGIONs by genomic overlap, and **its BOUNDARIES are
the boundaries that TOUCH its regions**. Every region contributes both of its boundaries, so a locus of `k` contiguous
regions carries `k + 1` boundaries — **its two outer ones included**, because a fragment crossing a locus's
boundary overlaps the locus and is therefore one of its EM candidates.

⭐ The outer boundaries are unambiguous, and structurally so: a locus is bounded by intergenic sequence and
intergenic regions carry no transcripts, so no two loci contend for a boundary. ⚠ Contention is
*rare rather than impossible* — 20–34 boundaries per flgap condition, carrying ~0.01 % of the mass — so
`priors.contended_boundaries` REPORTS it and an assertion would have died on real data.

⛔ **What this replaced, and why it must not come back.** `assemble_priors` folded each boundary's mass into
ONE flank region, because `_project_regions_to_loci` divides by `region_size_bp` and so cannot see a 0-bp
object. The fold then needed a second heuristic — an intergenic RE-KEY — to stop a locus's far-left boundary
vanishing into its dropped flank. ⚠ The fold predates the accumulator rewrite: shipped **v0.7.1** stored
the two boundary sides as separate banks and used them directly for RNA (no pooling, no owner) while
sending gDNA through `_pooled_seam_arrays`, which rejoined them and had to pick an owner. Its stated
reason is the CONTRACTION — right there, wrong for the COUNT. One object served two purposes.
`TRAPS: a-fold-grows-a-heuristic`.

⭐ **The prior's target, stated once**: `n_gdna` in `em_solver.cpp:apply_grouped_prior_update` is a soft
count of the gDNA fragments that are candidates in ONE multi-locus, each counted once — a multi-locus
being a connected component of transcripts linked by shared fragments. ⛔ **A first-base count of a
locus's fragments is therefore NOT this quantity**: it drops every fragment that overlaps the locus but
starts outside it.

### 3.1b-ii ⭐⭐⭐ COMPATIBILITY IS CHECKED BEFORE A FRAGMENT IS CALLED A CHIMERA — owner, 2026-08-19

> *"To be called a chimera, the paired reads would align to an incompatible orientation (not facing inward)
> or the implied fragment length would be outside the max fragment length parameter. … gDNA should be
> checked for compatibility before a fragment is called a chimera."*

A fragment is a chimera **only if its mates are genomically INCOMPATIBLE**:

    a different reference   OR   not facing inward   OR   implied fragment length > max_frag_length

⛔ **Transcript-set disjointness alone is NOT evidence of a rearrangement.** gDNA is genomically contiguous
and routinely spans two transcripts that share nothing — that is what the annotation looks like there. Such
a fragment carries no junction anywhere in it, so there is nothing to infer a chimera FROM.

⭐⭐ **THE ORIENTATION TEST NEEDED NO NEW PREDICATE, WHICH IS WHY THE RULE IS THREE LINES.** `build_fragment`
keys blocks by `(ref_id, ref_strand)` with R2's orientation FLIPPED, so both mates of an inward-facing pair
carry the SAME strand. `unique_strands.size() == 1` — the condition the code already called
`CHIMERA_CIS_STRAND_SAME` — *is* "facing inward".

⚠ **The length is the implied FRAGMENT LENGTH — outermost start to outermost end — and never `min_gap`.**
They differ by the blocks' own lengths, and the gap would admit a molecule longer than the library can
contain. Both twins gate that distinction (`test_resolution.py`), because breaking the C++ to use the gap
fired no native gate until one was written for it.

⛔ **What it cost while it was wrong** is `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`:
4,087 gDNA fragments per condition, 0.04 % of fragments carrying 2.4 % of every boundary crossing, invisible
to `region_contained` (every dropped fragment was a crosser) and to the conserved-mass gate.

### 3.1c ⭐⭐⭐ THE PRIOR ARBITRATES THE UNSPLICED POOL, AND ITS STRENGTH HAS NO KNOB

**Owner ruling, 2026-08-09/10.** Two statements, both settled, both measured.

> **A spliced fragment is pure RNA and never receives a gDNA candidate**
> (`em_solver.cpp`: `has_gdna = !is_spliced && isfinite(gdna_ll)`), so it does not compete with gDNA and
> **must not enter the pool-level prior**. An unspliced intergenic fragment has no transcript candidate
> and never becomes an EM unit at all. What the prior arbitrates is exactly the unspliced fragments
> inside transcript bounds.

⭐ **AND THE CROSSING→FRAGMENT CONVERSION IS POPULATION-BLIND, MEASURED AND BOUNDED** (2026-08-10). The
boundary term is converted by `q = mass/count`, measured per boundary on the POOLED population and applied to each
component. gDNA and RNA have different true `q` — not because of fragment length (the equal-length null
carries the LARGEST error) but because gDNA crosses boundaries in long intergenic regions where `q → 1` while RNA
sits among short exons. On the TOTAL prior the region term dilutes it to `Δphi` **+0.00013 … +0.00596**, at
or below calibration's own noise floor, and it is not repairable in production because the driver is
placement rather than length. ⛔ Recorded, not fixed: `TRAPS: a-pooled-conversion-applied-per-component`.
⚠ A DIFFERENT quantity from the `Δphi ≤ 5e-4` below, which is the shipped `f_g` against the unspliced
pool; this one isolates the `q` conversion with a perfect `f_g`.

⛔ Putting spliced RNA into `a_r` would penalise gDNA with fragments it could never have won. ⭐ Measured:
against the population it describes — gDNA units + UNSPLICED RNA units — the shipped claim
`phi = a_g/(a_g+a_r)` is exact to **≤ 5e-4** on all four flgap conditions, drained and undrained.
⚠ **Both flgap numbers on this page — the 20–34 contended boundaries above and this ≤ 5e-4 — were measured
on the `flgap_short`/`flgap_long` pair, which was DELETED on 2026-08-13** — and its configs and both
instruments that read it followed on 2026-08-17. They are kept as recorded and **cannot be re-run by
anything on disk**; nothing on the 16-condition ladder reproduces them. ⚠ An
entry once reported a "+0.07…+0.10 tilt" here; it divided by ALL RNA units
(`TRAPS: score-the-consumers-own-count`). The injection is population-matched too, by algebra:
`n_rna = S_r + U_r`, so `R = n_rna + a_r = S_r + (U_r + a_r)` puts the pseudo-count on the unspliced RNA
and leaves the spliced mass untouched; and `out[i] = raw[i]·(1 + a_r/n_rna)` is a UNIFORM scale, so it
changes no transcript's share of RNA — which is exactly "calibration says nothing about which
transcript".

> **The prior's STRENGTH is exactly one pseudo-fragment per real unspliced fragment, by construction, and
> there is deliberately NO knob** (owner, 2026-08-10: a strength parameter is another tunable and we are
> not adding one).

⭐ It follows from the conservation identity — `a_g + a_r` IS the locus's conserved unspliced count, which
IS the pool — and is measured at **Σa/Σpool = 0.999–1.000**, mass-weighted median `w` 1.000, with 0–1 loci
of ~1,300 above 2. So the MAP posterior is a 50/50 blend of calibration and the EM's own evidence.
⚠ Consequence, recorded not acted on: where calibration is wrong by half a unit, half the answer is too.

### 3.2 One strand convention

Everything is stored by **genome** strand (`CHANNEL_PLUS` / `CHANNEL_MINUS`). *Sense* / *antisense* is the
transcript-relative notion and is **derived, never stored**. Two strands exist and they are independent:
`align_strand` and `sj_strand`. An inferred splice is `sj_implicit`, never "inferred".

### 3.3 Two components only

**gDNA and RNA.** "RNA is RNA" — no mature/nascent split in the accumulator. Owner ruling.

### 3.4 Fragment length — ONE definition

`L` = genomic span minus cut introns. The scanner's rival histogram, `FragmentLengthModels` and the
transcript-space definition are **deleted**, and every histogram `build_fl_models` reads comes from the
payload, so a mixed-frame call is unrepresentable.

⭐ **A gap intron is cut on EVERY fragment**, not only unspliced ones, with the gaps the CIGAR already
explained excluded by **exact `(start, end)` equality** — overlap would let a different nearby intron
answer for one and make `L` too short.

⚠ `FragmentLengthModel` **singular** is the scorer and stays.

### 3.5 The five length pools

Each is pure **by construction**, and purity is what removes the circularity: a model is fitted only from a
population known to be one component, so nothing is estimated from the fragments it will later explain.

| pool | rule | component |
|---|---|---|
| `DNA_INTERGENIC` | contained in exactly one intergenic region | gDNA |
| `DNA_INTRONIC` | contained in exactly one intronic region | gDNA |
| `DNA_INTRON_EXON` | crosses exactly one boundary, flanks {intron, exon} | gDNA |
| `DNA_INTERGENIC_EXON` | crosses exactly one boundary, flanks {intergenic, exon} | gDNA |
| `RNA_SPLICED` | used an annotated sj, splice **observed** | RNA |

⭐ There is deliberately **no pool** for an exonic contained fragment or a multi-line crossing — those are
mixtures, and an impure pool is worse than a missing one.

⭐ **The pool is keyed on DETERMINACY, not provenance**: a fragment enters when exactly **one** hypothesis
survived, however it got there (TRAPS: a-purity-filter-is-a-length-filter).

⭐ **The two exon-crossing pools are gDNA**, because mature RNA never crosses an exon↔intron boundary
(TRAPS: mature-rna-never-crosses-a-boundary). They were filed for two milestones as "not gDNA"; they are.

### 3.6 Each pool is divided by its OWN opportunity

The gDNA model is fitted from **all four** gDNA pools, each divided by its own opportunity and then
combined — `calibration/gdna_opportunity.py`, derivation in `EQUATIONS.md` §4.3. The RNA model is the
sj pool divided by its own — `calibration/junction_opportunity.py`, `EQUATIONS.md` §4.2.

⛔ **The four gDNA pools must never be pooled raw** (TRAPS: opposite-tilts-must-not-pool), and every divisor is a
**probability**, not a count (TRAPS: divide-by-a-probability).

⭐ **The contained pair dominates off capture; the crossing pair dominates under it.** Under hybrid capture
the surviving off-target gDNA sits beside a probe, and a fragment beside a probe *reaches* the exon
boundary — so it stops being contained and becomes crossing. Fitting from the contained pair alone measures
the short half of one population.

⚠ **The second pass's scorer reads the same de-tilted pools calibration reads**, so there is one definition
of each length distribution rather than two.

### 3.7 The deposit weight is 1/opportunity

Not `1/length`. `EQUATIONS.md` §2 has the derivation, including the support factor `P(A > 0)` that
bounds where each form is model-free: the BOUNDARY/sj forms unconditionally (`P(w ≥ 2) = 1`), the REGION
form only within `w ≤ ℓ` (`TRAPS: a-cancellation-is-conditional-on-its-support`), and §2.3 the terminus
placement loss that a reach taper repairs.

---

## 4. The two-pass structure

**The accumulator arbitrates.** A fragment arrives with its hypothesis *set*; exactly one survivor
deposits, two or more are held WHOLE in the **side buffer**. A fragment's unsequenced mate gap may hold no
intron, one, or several, and *which* cannot be observed — the bases are not there.

**The second pass drains it**, between the scan and calibration:

1. **score** from pass-1 evidence alone — `score = ρ × f(L) × s`, `EQUATIONS.md` §10;
2. **one multinomial draw** per fragment;
3. **re-deposit through the SAME `deposit`**, with the chosen hypothesis alone — a set of size one, so
   arbitration is degenerate and the ordinary rules decide.

⭐ **One tally path.** There is no second deposit implementation and no duplicated crossing logic, so
byte-identity with the executable specification holds for free rather than by argument.

⭐ **It runs before calibration, and that is the structural decision.** Every input the score needs comes
from pass 1, so **calibration runs exactly once, on the complete tally.**

⚠ **The models the SCORER uses are pass one's; the models CALIBRATION uses are the drained tally's.** Fit
once, score once, drain once, stop. The confident set is biased short, so feeding the drained anchor back
into the score would prefer the shorter — more-spliced — path, and that loop can run away.

⭐ **The FRAGMENT is stored, never its consequences.** Object ids are large, derived, and would have to be
kept consistent with the partition; the fragment is small and replays exactly.

⛔ **The side buffer is the one bank whose ORDER is observable** — it is a list, not a sum — so the C++
export sorts it on the record's own content before it crosses the ABI (TRAPS: integer-channels-reproduce).

**After the drain nothing is held**: the bank is empty and the deferred counters are 0, so "the counter and
the fragments it counts are the same population" stays absolute. Pass one's numbers live in `DrainQC`.

### 4.1 Settled sub-decisions

| | |
|---|---|
| ⭐⭐ **a composition may be imputed across a step iff the source SUPPLIED both components AND the two objects measure the same RNA POPULATION** | owner, 2026-08-04: *"is the source of the message measuring the same thing that I am measuring?"* — if yes, attribute the density discrepancy to capture enrichment; if no, you cannot tell enrichment from a population difference. Derivation + the genomic form of the predicate: `EQUATIONS.md` §3.5b. ⛔ **Termini only** — DONOR/ACCEPTOR change the population too, but their flux is *measured* and the SPLICE IN and the deconvolution route it |
| ⭐ **the population test is written in GENOMIC terms, never TSS/TES** | the strand flips which flank a terminus implicates; `region_geometry.terminus_flank_gain` is the one home, gated on two mirror-image annotations (TRAPS: substitute-the-definitions-first family) |
| ⭐ **the mass pin carries the same licence, plus the structural case** | it fires iff no BELIEF can reach its budget: the composition was supplied, **or** the destination is a structurally pure-gDNA object whose `f_g = 1` is structure and whose `M/E_g` is therefore an observation. `EQUATIONS.md` §3.5c, TRAPS: the-pin-had-a-fixed-point/TRAPS: no-belief-not-no-numbers |
| gDNA's strand term is **0.5** | double-stranded, no sense direction. A fitted mixture marginal was implemented and refuted (`EQUATIONS.md` §5.3) |
| a flat-zero factor is **skipped**, not multiplied | TRAPS: an-all-zero-factor-is-inert |
| the draw is keyed on **queue position**, never content | `EQUATIONS.md` §10 |
| ambiguous assignment stays **integral** | TRAPS: fractional-mass-is-the-problem |

### 4.2 Still open

| | |
|---|---|
| **D-4** — should a density carry its weight of evidence? | a density of 0.01 from 1000 fragments and from 1 are not the same statement |
| the traffic term's **Poisson likelihood** | `ρ` enters as a hard zero; `P(0|λ,E) = e^(−λE)` is not zero (`EQUATIONS.md` §10). Deferred on a measurement, not on a missing piece — the exposure it needs is now available |

---

## 5. nRNA components

Not per-transcript shadows: unique nascent spans keyed by `(ref, strand, start, end)` are shared across
transcripts and materialized as ordinary transcript rows in `index.t_df`, flagged `is_nrna` /
`is_synthetic`.

⚠ On a **non-synthetic** row `is_nrna` means "single-exon, so mature ≡ nascent" — **not** "manufactured
span". **The real-transcript filter is `~is_synthetic`, alone** (TRAPS: nrna-does-not-mean-synthetic).

---

## 6. Code layout

**Python.** Top level: `cli` `pipeline` `config` `index` `scan` `scoring` `buffer` `scan_payload`
`scan_cache` `locus` `locus_partition` `scored_fragments` `estimator` `strand_model` `frag_length_model`
`second_pass` `splice` `splice_blacklist` `native` `gtf` `transcript` `annotate` `stats` `types`, plus the
`report/` and `sim/` subpackages. ⚠ `_accumulator` is GONE — the row-view façade the native `Tally` used
to come through was deleted (`native.py` says so at its import).

`calibration/`: `calibrate` (orchestrator) · `splice_graph` (the v8 index) · ⭐ **`sweep` (THE BACKBONE)**
and **`messages/` (the POLICY: `silent` `head` `variance`)** · `region_chain` `region_geometry`
`region_init` · `substrate` `region_arrays` `signature` · `effective_length` `capture_eff_length` `fl`
`sj_opportunity` `gdna_opportunity` · `strand_likelihood` `gdna_strand` `strand_balance`
`strand_deconv` `strand_summary` · `density_deconv` `density_model` `landscape`
`abundance_landscape`
`simplex_logodds` `derive` `run_fill` · `priors` `result` `errors`
`diagnostics` `track` · `_layers` (the layering the imports already had).
⚠ **Re-derive this list rather than trusting it** — `scripts/design/module_census.py` reads it off the
AST. It named `node_chain` / `node_geometry` / `node_init`, `junction_opportunity` and a `simplex` — five
modules with no file on disk — until this was repaired on 2026-08-17.

**C++** (`src/rigel/native/`, nanobind, C++17, `-O3`, LTO, OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | BAM parsing, fragment grouping, model training, the accumulator |
| `_em_impl` | `em_solver.cpp` | per-locus EM, connected components (Kahan summation, SIMD `fast_exp.h`) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring (`-ffast-math`, no SIMD) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges |
| `_cgranges_impl` | vendored | interval overlap |


### 6.1 ⭐⭐⭐ THE SOLVER IS A BACKBONE PLUS A POLICY — settled 2026-08-07, gated on BYTE-IDENTITY

The belief-propagation solver was ONE 1,635-line function in which the shape of the solve and every
argument about what a message should say were interleaved. It is now two things, and the split is
**structural rather than stylistic**:

| | | |
|---|---|---|
| `sweep.py` | ⭐ **THE BACKBONE.** Two directional scans, one combine, one ψ solve, one write-back, **five assertions**. | It knows nothing about capture, splices in, reframes, pins or enrichment — ⛔ and `test_sweep_backbone.py` asserts those words appear in none of its IDENTIFIERS, read from the AST rather than grepped, because grepping matches the docstring that says they are absent |
| `messages/silent.py` | ⭐ `SilentPolicy` — sends nothing. **THE DEFAULT**, five lines. | A reader who holds `sweep.py` plus this holds the entire working system |
| `messages/relay.py` | `RelayPolicy` — every operator the evolved solver carried, each behind a NAMED switch (**17** of them) | So the panel prices them ONE AT A TIME rather than as a block |
| `messages/variance.py` | was `enrichment_frame.py` — the policy's variance toolbox | ⚠ `count_logvar` is also imported by `region_init`; it has ONE home and this is it |

⛔⛔ **`SilentPolicy` BEING THE DEFAULT IS NOT A STATEMENT THAT MESSAGE PASSING IS UNNECESSARY — §0c
proves the opposite, and `RelayPolicy` ALREADY CONTAINS THE OPERATOR THE EXON NEEDS** (the SPLICE IN / SPLICE-IN,
§0c.1). Read §0c.2 before quoting the price of turning the flag back on.

⭐⭐⭐ **HOW IT WAS ACCEPTED, and this is the part that generalises: a restructure is gated, a rewrite is
not.** Two `TRAPS: byte-identity-gate` gates, opposite in direction, both passed:

| | must be | result |
|---|---|---|
| `--arm backbone_relay` (every switch ON) | **byte-identical to `base`** | ✅ **1,872 / 1,872 scored fields, 72 / 72 rows** |
| `--arm backbone` (`SilentPolicy`) | **byte-identical to `msgfree_all`** | ✅ **1,872 / 1,872 scored fields, 72 / 72 rows** |

and, per-array on one real 70,176-slot chain (`scripts/design/backbone_parity.py`): **421,056 output
elements and 18,245,830 diagnostic elements, zero differences.** ⭐ The second gate is not only plumbing —
one arm runs the whole relay and discards it while the other never runs it at all, so passing it **PROVES
the relay reaches the answer ONLY through ψ's four imputed channels**.

⛔ **The alternative was tried and is why this route was chosen.** A clean rebuild came out **+103 %** and
took two sessions to decompose into a correct derivation, a UNIT ERROR (a log-odds delivered into an
angle's grid) and a STRUCTURAL disconnection (a one-slot-step channel on a bipartite chain reaches 0
REGIONs). One of the three was a typo-class error no amount of derivation review would have caught. **A
refactor gated on byte-identity has exactly zero of that risk: any difference is a bug, full stop.**

#### The interface, and its one contract

```python
relay = policy.prepare(ctx)                  # one working object per sweep
step, publish = relay.scan(backward=False)   # the per-hop kernel; None ⇒ nothing to relay
msg = relay.deliver(left_state, right_state) # -> PsiMessage, from the NEIGHBOURS only
```

⛔⛔ **TRAPS: a-message-from-the-destinations-belief, and the backbone enforces the enforceable half BY CONSTRUCTION.** `deliver` is handed
:class:`NeighbourState`, whose relayed arrays are **already gathered at the source slot** — so a policy
holding one has values FOR THE SOURCE and no way to ask the same array about the destination. TRAPS: a-message-from-the-destinations-belief's nine
costumes were all a message built from the destination's own relayed belief, and none of them is
expressible through that type. ⭐ A gather is exact, so making the rule structural costs no bits.

`StepContext` splits its fields under three headings — **observations** (either end), **geometry /
structure** (either end), **beliefs** (SOURCE-SIDE ONLY) — and the heading is what turns TRAPS: a-message-from-the-destinations-belief from a
discipline into something a reader can check. ⚠ **One belief field is read at the destination by the
shipped policy and it is a named, measured DEBT rather than an oversight**: `belief_fg` reaches the
reframe's frame pair, so the frame at a hop is a function of the destination's belief, at slots carrying
**57–77 % of library mass**.

#### The five assertions, and why they live in the backbone

⭐⭐ **A future policy can be as wrong as it likes and still cannot commit any of these** — each has shipped
at least once.

| the backbone asserts | it would have caught | on HEAD, `g50 ss0.50 capture_on`, pass-0 |
|---|---|---|
| `deliver` sees only the two NEIGHBOUR states | **TRAPS: a-message-from-the-destinations-belief — nine recurrences in nine costumes** | structural — inexpressible, not checked |
| every message mode lies inside its coordinate's own grid | **TRAPS: off-grid-message-mode** — the tilt bug, 74 % of `g00`'s error | ✅ tilt **0 / 4,795** · ⛔ gDNA **15,240 / 50,984 (29.9 %)** · RNA+ 15.0 % · RNA− 14.0 % · λ 0.46 % |
| every delivered share is in `[0, 1]` | the over-unit certified-RNA claim | ⛔ **3,013 / 15,629 (19.3 %)**, and the SUM ⛔ **31,174 / 50,984 (61.1 %)** |
| `\|T\| ≤ 3` | **AXIOM 0**, made executable | ✅ 0 / 70,176, and **9,912** slots reach 3 so it is not vacuous |
| the write-back touches only `solvable` slots | the basis mismatch that made an TRAPS: byte-identity-gate read `max\|Δ\| = 1.0` | ✅ 0 / 32,817 |

⭐⭐⭐ **AND THE ASSERTIONS PAID FOR THEMSELVES ON THEIR FIRST RUN — six of the ten counts above are new
facts about the shipped message layer.** The two that matter most:

* ⛔⛔ **61.1 % of live message packets assert that the three components TOGETHER account for MORE fragments
  than the slot observed.** That is the identity `Σ_c ρ_c·E_c = M` the mass pin exists to restore, and the
  pin is licensed in only two states, so everywhere else the residual is *delivered* rather than fixed. ⭐ It
  is consistent with a number already in the tree — `messages/variance.py` records the over-claim on 52–71 %
  of regions — but nothing had surfaced it as a **checkable invariant**, so nothing could rank it.
* ⛔ **29.9 % of live gDNA level modes are outside ψ's own log-share grid**, whose domain is
  `[log σ(−L), log σ(+L)] = [−10.000045, −4.54e-5]` and **not** `(−inf, 0]`. The cause is the `_EPS = 1e-9`
  floor: `log(1e-9) = −20.723`, **10.72 nats** below the low end. ⚠ **This is not 29.9 % of the error, and
  the difference from TRAPS: off-grid-message-mode is the honest part**: TRAPS: off-grid-message-mode's tilt pinned at the *wrong* corner, whereas a low-side
  share pins at "as little of this component as the grid can express", which at a slot that genuinely holds
  almost none of it is the *right* answer. What is certain is TRAPS: off-grid-message-mode's mechanism — no interior minimum, so
  precision buys a CORNER rather than a location. Whether the corner is right is **unmeasured**.

⚠ **The low-side check is the one the first draft of this section did not have**, and it was found by an
adversarial read rather than by the panel: the gates were already passing byte-identically, so nothing about
the restructure would ever have surfaced it. ⭐ `share_sum_at_most_one` is a SUPERSET of
`share_in_unit_interval` by construction — report both, never add them.

⛔⛔ **AN ASSERTION THE SHIPPED POLICY VIOLATES IS WAIVED WITH ITS MEASUREMENT, NEVER WIDENED.** The two
that fire are recorded in `sweep._KNOWN_VIOLATIONS` with the number beside each, and `test_sweep_backbone.py`
asserts that every waiver carries a written reason. Widening the predicate to fit the defect is how a gate
becomes vacuous (TRAPS: perturb-every-gate/TRAPS: a-gate-that-reconstructs); a waiver keeps the count VISIBLE and rankable.

⭐ **Each assertion also reports how many slots were ELIGIBLE for it**, because TRAPS: could-the-arm-have-fired: an assertion reporting
zero violations where its predicate can never fire is not evidence of anything.

⚠ **One tolerance in the family is derived and not tuned** (TRAPS: no-magic-numbers): the grid-domain check allows
one **grid SPACING** beyond the outermost grid point, because within one spacing the penalty's minimum is
still that boundary CELL — for the tilt, the legitimate answer "all RNA on one strand". It is load-bearing:
the shipped tilt mode overshoots `π/2` by **2 ULP** at 63 of 4,795 live slots, being a convex mean of two
messages that AGREE on `τ = ±1`, and TRAPS: off-grid-message-mode's real overshoot is **57 spacings** — four orders of magnitude
apart, so no threshold is being bought.

---

## 6b. ⭐⭐⭐ ψ's REFERENCE HAS A PER-OBJECT MEAN — the rulings behind it (2026-08-16)

⭐ **The reference does NOT have to be library-wide, and the machinery was already per-slot.** ψ solves
one REGION or one BOUNDARY at a time and the gDNA arm's fitted term is already `(n_slots, K)`; the
reference was the only scalar left in it. `CompositionPriors.location` carries the per-slot mean and
`simplex_logodds._location_term` writes it; `None` means the term is not written at all, which is the
shipped behaviour. The derivation is `EQUATIONS.md` §9c.

⛔⛔ **RNA IS THE RESIDUAL AND IS NEVER PREDICTED.** gDNA is near-uniform over the genome and measurable
before any solve; RNA spans six orders of magnitude with essentially no genomic autocorrelation, so a
pooled RNA density is not a population parameter and pooling splice-junction flux across the genome is
inadmissible (owner, 2026-08-16). The reference's location is therefore `m_i = ρ_g·E_g,i / M_i` — the
gDNA an object's own density predicts, as a share of what it holds — and whatever is left is RNA.

⚠ **THIS RULING IS SCOPED TO THE OFF-TARGET CASE AND §0c.3 CARRIES THE SCOPE.** At an EXON under hybrid
capture the two roles are exchanged — capture makes gDNA per-exon and arbitrary, while RNA is *measured*
locally by the adjacent BOUNDARY's sj flux — and that is a NEIGHBOUR'S OWN COUNT, never a pooled RNA
density, so nothing above is weakened. Read the two together before designing a reference.

⭐⭐ **THE BOUNDARY AXIS SPLITS ON WHETHER MATURE RNA CAN CROSS, NOT ON WHETHER A SPLICE JUNCTION
ATTACHES** (owner, 2026-08-15). A BOUNDARY owns the fragments that cross it CONTIGUOUSLY, and mature RNA
can do that only where the template is contiguous exon on both sides. So an `exon|intron` boundary is
near-pure gDNA under sparse unspliced nascent, while an `exon|exon` boundary — an alternative splice site
inside a contiguous exonic stretch — is crossed freely by mature RNA and is not an anchor. ⛔ The
predicate is the solver's own `mrna_active`. A pool defined as "a splice junction attaches here" measured
true `f_g` **0.0000 over 955,428 fragments** at the zero-gDNA control, because it lumps `exon|exon` in.

⭐⭐ **THE REFERENCE AND THE gDNA LANDSCAPE PARTITION THE OBJECT UNIVERSE RATHER THAN COMPETING.** Where
the annotation determines the answer — intergenic and intron REGIONs, `exon|intron`, `intron|intron`,
gene-edge and opposite-strand `exon|exon` BOUNDARIES, 47.5 % of slots — the reference carries it before
any solve exists. Exons, same-strand `exon|exon` and AMBIG have no structural claim and are the
population the landscape is fitted to serve. ⛔ **The overlap is bounded ONLY where evidence exists**: the
reference is worth one pseudo-fragment, so it is swamped by any evidence channel — but on an
evidence-free object posterior = prior at any depth, and there the reference and the landscape are the
only two voices.

### 6b.1 ⛔⛔⛔ THE REFERENCE LOCATION — REFUTED AND DELETED (owner, 2026-08-24)

⛔⛔⛔ **ψ HAS NO REFERENCE LOCATION.** `_location_term`, `structural_reference_location` (the
annotation-derived `m = 0.75` at `¬mrna_active`), `measured_reference_location` (the measured
background location at ss-intron REGIONs) and both config flags (`structural_reference`,
`measured_intron_reference`) were **deleted outright** on 2026-08-24. The reference is the symmetric
Jeffreys measure and asserts nothing; the location concept may not come back as a flag, a constant,
or a "weakened" variant.

**The refutation, in one paragraph.** A location is a prior ASSERTION at fixed strength. Where the
strand channel carries information it was worth ~one fragment as documented; where it does not
(κ = ½, or composition near a vertex, where the λ-axis Fisher information `∝ [f_g(1−f_g)]²`
collapses) it was **the entire answer at any depth** — measured 0.7471 at every N from 10 to 10⁶ at
κ = ½, and at zero refits it decided 100 % of the claimed-boundary error. Pass-0's job is to LEARN
the prior; asserting one first is circular, and the measured intron form (fit on the data being
deconvolved) was circular twice over. The concept was refuted by the owner with external review on
2026-08-24; the measurements above are re-derivable from `pass0_claimed_ab.py` and
`calibration_vs_oracle.py` on the cached panel.

**What replaces it.** Background information enters as LIKELIHOOD terms whose precision scales with
counts — the intron-factory λ-factor (`density_lambda_factor`) at ss-intron REGIONs today, and the
planned learned-enrichment extension of that channel to the other claimed populations (measured at
boundaries: the intergenic rate predicts certified boundary gDNA to 0.5 % capture-OFF and is 38–50×
low under capture — the enrichment must be LEARNED, never gated on a capture boolean; owner ruling,
2026-08-24: capture is a spectrum).

**The measured price of the deletion** (`calibration_vs_oracle.py`, all 16, same-session before/after,
never pooled): the large in-scope error masses (stranded × capture-ON `g50`/`g98`) +4–8 %; smaller
rows up to +42 % (`g98` ss.99 OFF region); `g00` capture-OFF controls improve; `g00` capture-ON
+30–45 % (the deleted measured form's zero-control win — the density channel is the recorded path to
winning it back). Accepted: the score was partly measuring the term, not the algorithm.

⛔⛔ **THE PRICE CONCENTRATES IN THE FAN-OUT ARM, and that is §H.8's requirement 5 measured again**
(`pass0_claimed_ab.py`, claimed populations, shipped refits, same-session before/after): under
SILENCE most rows move ≤ 10 % and both `g00` unstranded controls IMPROVE ~2–4×, but under FANOUT the
transported flank beliefs lose their anchor — the deferred stratum 2.8–3.1× worse (`g50`/`g98`
unstranded ON), `g98` unstranded OFF boundaries 2.0×, stranded ON boundaries +24–42 % — and on
stranded × capture-ON the fan-out now LOSES to silence on boundaries. The reference and the messages
were priced jointly ON; deleting one re-prices the other. *(Historical: `FanOutPolicy` itself was
deleted on 2026-08-24 once the RNA anchor restored the slots' own evidence — §6b.2.)*

⚠ **What SURVIVES from the old section**: the two refutation mechanisms (strand asymmetry; the
intron-vs-intergenic density factor, `τ_fac = 161.4` at an intron, alive unstranded) are likelihood
channels and remain shipped; `TRAPS: a-priors-curvature-is-not-the-datas-information` and
`TRAPS: a-refutability-test-needs-the-refuting-channel-in-the-fixture` stand; the relay defect the
old reference EXPOSED (a claim carried across an exon↔intron population change) is a LICENCE bug
recorded as a constraint on the relay rebuild (§0c.2) — its strict-xfail fixture died with the
deleted gate file and the constraint carries the record.

### 6b.2 ⭐⭐⭐ THE RNA-ANCHORED EVIDENCE FACTOR — BUILT AND SHIPPED ON (owner design + ruling, 2026-08-24)

**The mechanism** (`calibration.rna_anchor`, default `CalibrationConfig.rna_anchor = True`; the
derivation, the estimator ledger and the recorded residuals live in the module docstring; gate
`tests/calibration/test_rna_anchor.py`, 19 cases, every perturbation verified firing): the RNA side
of the unspliced count is anchored on quantities hybrid capture cannot mis-scale — certified splice
flux at complete-flank exons, the adjacent intron's excess-over-background nascent rate at eligible
ss-intron boundaries — as a likelihood summed into the intron factory's per-slot factor array. No gDNA rate for an enriched slot appears anywhere; anchor and target share
each exon's own probe footprint, so enrichment cancels pair-locally (measured: the RNA-frame
boundary→exon ratio is capture-invariant, median 1.06–1.24 both states, where the gDNA frame
shifts 2.9–3.4×). ⛔ Capture GATING remains refuted (capture is a spectrum); this factor needs no
gate by construction.

**The estimator, hardened same-day through three measured iterations** (each defect caught by the
release metric and pinned by a new gate): ① the transport's CENTER and SPREAD are fitted jointly
from the two lower quantiles of the pair residuals (gDNA can only inflate a residual upward, so
the lower quantiles are RNA-only; two log-normal quantiles give both parameters); ② the fit
carries a SELF-CONSISTENCY GUARD — it predicts its own negative-residual fraction Φ(−m/s) and
refuses when the observed fraction disagrees beyond the binomial noise band, which is what stops a
gDNA-saturated tail (`g98`) from fitting a spurious center (measured: unguarded, the center gave
back most of the high-gDNA win); ③ on refusal the CENTER goes to zero but the WIDTH survives via
the MAD left-tail fallback, max-combined with the two-flank disagreement (measured: dropping the
width with the center over-tightened a zero control 207k → 402k). ⛔ The deep lesson, three
iterations paid for: ANY in-sample residual estimate of the transport has validity that depends on
the unknown gDNA level itself — the guard makes that dependence explicit and refusable instead of
silent.

**ROUND 2 (same day) — derived under adversarial review, prototyped, priced, shipped.** Three
changes, each priced alone: ① ⛔ **the ROUTE-SUM pooling, a bug fix** — a flank's junctions are
disjoint routes (each molecule crosses exactly one), so the flank rate is `Σ_J flux_J / A_J`;
round 1's ratio-of-sums under-predicted k-route exons ~k× and its route-count asymmetry
MANUFACTURED the heavy-tail "transport dispersion" (the deep-flux floor halved under the fix) —
the single largest lever; ② the count-scale Gaussian replaced at EXONS by the intron factory's own
pattern — the Gamma⊗Poisson **NegBinomial marginal (`size = flux + ½`) averaged over the
multiplicative transport scatter by a median-preserving equal-mass quantile quadrature** (bulk-soft,
tails bounded; Gauss–Hermite nodes were tried and measured bulk-TIGHT — a priced wrong shape);
③ the nascent term MARGINALIZED — both point-estimate forms are positively biased at nascent-free
truth (Jensen: the plug-in clamp by ~0.40σ, the truncated-posterior MEAN by ~0.56σ, measured as a
`g98`-ON regression), so quantile nodes of the truncated-excess posterior enter the quadrature and
a clean intron keeps its atom at exactly zero. The excess converts with the intron's RNA
opportunity (panel-invisible by the equal-length design; the fl-gap side panels falsify it).
⛔ **The BOUNDARY factor keeps the round-1 guarded-Gaussian family**: its prediction is near-zero
wherever capture enriches the boundary crossing against the intron the anchor reads (the cliff),
and a quadrature there asserts that near-zero with counting-only width — priced as a 734 → 27.5k
zero-control leak, now pinned by a capture-cliff flat gate and a builder-dispatch sentinel gate
(the revert once shipped as dead code behind 17 green gates). At capture-ON the boundary factor is
therefore ~flat BY DESIGN and the boundary wins arrive via the exon factor and the messages.

**Priced on 0.8.0's own metric** (whole-library misplaced fragments vs oracle, the graduation
table in the round-2 review): the round-1 clean-library residual falls 2.7–3.9× — `g00`-OFF relay
148.8k → 55.4k (claimed E 11.0k), `g05`-OFF silent 124.7k → 31.9k (claimed E 80.7k → 5.5k), the
`g50`-OFF stress-nascent row restored to the pre-anchor base (claimed E 22.9k → 6.6k vs base
6.4k) — while every high-gDNA win holds (`g98`-ON silent 80.7k → 68.5k region with claimed E
5.7k; relay 182k → 146k; `g98`-OFF claimed B 40.9k → 21.5k) and the zero controls sit on the
anchor-OFF base (`g00`-ON claimed B 718 vs 710). Remaining `g00`-OFF relay gap vs the pre-anchor
base: 3.2× (55k vs 17k); the named candidates are the review's unaddressed items (single-flank
center bias, sj strand-column matching, shared-fragment pair correlation, the capture-ON +12 %
mature offset).

⛔⛔ **OBSOLESCENCE TRACKING (owner directive, 2026-08-24: remove what the anchor obsoletes).**
Candidates, each pending the policy re-contrast ON THE ANCHORED TREE and none deleted yet:

* ✅ **EXECUTED 2026-08-24**: `messages/fanout.py` (the whole policy, its stage-3/4 transfers and
  the receiver-side mismatch deflation), `test_fanout_policy.py`, the `"fanout"` config value, the
  `policy_fanout` ladder arm, `pass0_claimed_ab`'s `--dissect`/`--sweep` survey, and
  `StepContext.claims` (dead surface once no policy read it). The final 16-condition re-contrast
  on the anchored tree measured the fan-out DOMINATED on every row — never uniquely best:
  stranded rows belong to silence+anchor, unstranded-exon and deferred wins to relay, zero
  controls to silence. No shipped behavior moved (`message_policy` was already `"relay"`).
* `simplex_logodds.ONE_SIDED_RNA` and the one-sided certified-RNA bound machinery — the anchor IS
  the two-sided certified-RNA statement at complete flanks, made honest by measured variance.
* The messages' certified-RNA transfer channel (`rna_imp_*` at stage-4 destinations) — same
  information, now delivered as own evidence.

⭐ **HYBRID CAPTURE NEEDS NO DETECTION STEP.** Measure the gDNA density at both the off-target anchors
(intergenic + intron REGIONs) and the in-gene `exon|intron` boundaries; their RATIO is the enrichment —
measured **0.98** without probes and **113–114** with, a 116× separation with no threshold and no flag.
⚠ The in-gene anchor is a DETECTOR and not a calibrated level: it under-reads the on-target gDNA density
by **2.6–3.6×** because it sits at the EDGE of the probe footprint.

### 6b.3 ⭐⭐⭐ THE ANCHOR IS A MESSAGE — the citizenship ruling and the integration (owner, 2026-08-25)

⭐⭐⭐ **The ruling.** The spliced-fragment anchor is message propagation — a one-hop imputation from
the flank boundary into the exon — and it may not exist as an appendage beside the message
framework. The SENDER publishes its spliced-fragment observation unchanged; the RECIPIENT decides
how to use it, and the recipient's arithmetic is the route-sum + NB marginal §6b.2 describes. The
2026-08-24 build had it riding the local λ-factor of EVERY policy — the silent control included —
which was refuted as architecture: "we witnessed building an appendage rather than using the
framework that already existed."

**The decision made under the ruling: NO new policy.** The anchor is a new *stream* — the
CERTIFIED-FLUX stream — not a new composition scheme, so it lands in `RelayPolicy` behind one
switch (`RelaySwitches.certified_flux`), `SilentPolicy` stays silent. ⭐ The survey behind the decision found this is a HOMECOMING: the
framework already carried the observation (`StepContext.inv_sj_lo/hi`, per-fragment
reciprocal-opportunity — route-sum-correct by construction) and the relay's splice-in stream
already delivers a certified-RNA claim to exons — built from `sj_count/eff_sj` as a per-face
density, which is the RATIO-OF-SUMS pooling §6b.2's review refuted. ⛔ So the relay's own
splice-in carries the k-route under-read, and the stream is the UPGRADE of that claim: the same
certified observation, correctly pooled, delivered as a count likelihood instead of a Gaussian.
The Gaussian splice-in claim, `ONE_SIDED_RNA` and the `rna_one_sided` plumbing are its
obsolescence candidates, measured at the weak-message re-price.

**The contract** (gates in `tests/calibration/test_rna_anchor.py`, every perturbation verified
firing): `PsiMessage.lam_rows` — an `(n_slots, K)` λ-factor row array in ψ's general evidence
currency — is the one new channel; the backbone asserts its shape/finiteness and sums it into the
FINAL solve only, so phase-A (`build_region_init`) and the own-evidence precision
(`density_factor_precision`) never see an imputation. `calibrate` prepares the evidence once
(`rna_anchor.prepare_flux_evidence`, route table included) and constructs ONE relay instance with
it; the relay's grid-keyed memo builds the rows (`rna_anchor.flux_rows`) per bracket.
`config.rna_anchor` gates the stream: live iff `message_propagation ∧ message_policy == "relay"`.
⛔ **The behavioural consequences are deliberate**: the silent arm reverts to the pre-anchor
CONTROL (its panel rows regress by construction — that is the control being restored); anchored
slots lose own-evidence status, which changes the mismatch-deflation inputs and the refit
training populations honestly; the shipped candidate is relay-with-stream, accepted per stratum
against the §6b.2 graduated numbers. ⚠ Deferred deliberately: folding the route-sum into
`region_geometry` (fixing the relay splice-in's pooling at the source) and per-strand column
matching — both ride the weak-message re-price.

## 6c. ⭐⭐⭐ ψ's COMPOSITION IS A POINT ON THE SIMPLEX, AND CLOSURE IS STRUCTURAL (2026-08-17)

⭐ **THE COMPOSITION HAS TWO DEGREES OF FREEDOM, NOT THREE.** ψ solves a point on the 2-simplex,
parametrised by `λ` (the gDNA-vs-RNA LEVEL) and `θ` (the RNA-internal TILT, a SHARE with no absolute
scale). The composition is their IMAGE — `simplex_logodds._compose`:

    f_g  = the ½-quantile of the λ posterior        RNA total := 1 − f_g
    f_pos = (1 − f_g)·w        f_neg = (1 − f_g)·(1 − w)

so `f_g + f_pos + f_neg = 1` identically. ⭐ **Measured on real conditions: 100.00 % of published objects
close, on both axes and every annotation class, min/p5/p95 all exactly 1.0000** — against 74.7 % of REGIONs
and 77.2 % of BOUNDARYs before.

⛔⛔ **WHAT IT REPLACED WAS NEVER A DECISION.** `f_g` was the posterior MEDIAN and `f_pos`/`f_neg` were
independent posterior MEANS of the grid quantity `1 − f_g`, so `SUM = 1 + median − mean` — the closure
error was exactly the posterior's SKEW. The median was argued for; the RNA fractions fell out as
expectations of an array nobody revisited.

⛔ **TAKING MEANS EVERYWHERE ALSO CLOSES AND IS REFUSED.** It closes by linearity of expectation and scores
**1.352 / 1.573 / 3.756** on the three in-scope strata and **1.801** on the zero control
(`vertex_ceiling.py --arm psi_mean`, 16 conditions): the median is closer to truth at both simplex
vertices, where 49–83 % of in-scope error lives. ⛔ Nor is this renormalisation at publication — nothing is
rescaled to hide a residual; `f_g` and the RNA total are exact complements BY PARAMETRISATION, and a tilt
is estimated as a share because a share is what it is.

⭐⭐ **THE ½-QUANTILE IS CONTINUOUS AND IS READ ON λ.** Snapping to a lattice point put up to half a grid
step into `f_g`, and deriving the RNA total from it propagated that into the whole composition
(`TRAPS: deriving-one-coordinate-propagates-its-error`). ⛔ And the interpolation must be on `λ`, where the
lattice is uniform: in `f_g`-space a concentrated posterior comes back biased toward ½ by **2.71e-03** at
`n_grid` 60, on λ it returns its own grid point to **2.2e-16**
(`TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`). ⭐ That also restores the property the
estimator exists for: `|median(1−f) − (1 − median(f))|` went **8.45e-02 → 3.3e-15**.

⚠ **ADMISSIBILITY IS ENFORCED INSIDE THE MAP**, not by the caller: an unclamped or strand-blind share
produces a NEGATIVE fraction that still sums to 1 — a composition that passes every closure check and is
nonsense. ⚠ **The two degenerate slots are DIFFERENT and this passage used to fuse them.** A slot with **no
counts** (`u_pos + u_neg == 0`) publishes `(0, 0, 0)` — "no data", not a composition claim — and the
sweep's `solvable` write-back leaves it at `region_init`'s signature-binary init, which is a simplex point.
A slot with **no admissible strand** has no RNA to place: `_compose` returns `(0, 0)` and its composition
is `f_g` alone, which is the honest statement and not a closure failure — and nothing dispatches such a
slot to a ψ solve in the first place.

## 7. Where the error is, structurally

Two measurements that shape every design choice above, both worth keeping because they are about the
*structure* rather than about a particular library:

⚠ **Both decompose the error by OBJECT CLASS.** The decomposition by STRATUM — which cell of
strandedness × capture carries it, and therefore which three strata are the 0.8.0 optimisation target —
is §0b above, and the two are read together: an object class carrying the error is a mechanism, a stratum carrying it
is a scope.

**Regions and boundaries measure different components.** The gDNA/RNA opportunity ratio is **0.25 at a crossing
point** (4× RNA-selective) against **115.7 at a 100 bp region**, 25.5 at 147 bp, 2.5 at 300 bp, 1.19 at
1000 bp. A 200 bp RNA fragment cannot fit in a 147 bp region, so a short region is a good gDNA measurement and
says nothing about RNA. **Carry per-component precision, not one scalar.**

**Most of the error is in objects with no evidence of their own.** ⚠ **Measured over a 32-condition sweep
on a panel that no longer exists** — the ladder was rebuilt to 16 conditions on 2026-08-13 — so read the
SHAPE, which the re-derivation below reproduces on a different panel with a different instrument, and not
the digits as current. Objects carried entirely by neighbours were 54.1 % of mass and **91.9 % of all
error** (6.2× the rate); objects
with own evidence 29.6 % / 8.1 %; structurally-locked pure-gDNA objects 16.4 % / 0.0 %. And 80.5 % of
total error was honest under-determination — only 1.9 % was confidently wrong.

⭐ **Re-derived per object against the origin-split oracle** (`pass0_vs_oracle.py`, 4 contaminated
conditions, both axes): the relay carries **93.1 %** of pass-0's error, at 1.0–3.7× its mass share on
every arm, and structurally-locked objects carry **0.0 %** — the same shape, on a different panel, with a
different instrument. Two things the older sweep could not say:

* ⭐ **On an unstranded library the DENSITY model carries the entire own-evidence budget.** At κ = ½ the
  strand λ-term is exactly 0 (`EQUATIONS.md` §5.2), and the intron factory still reaches **100 % of
  intron-REGION mass — both-stranded as well as single-stranded**. Density, not strand, is what makes an
  unstranded library solvable at all.
* ⛔ **But the factory is wired to `(REGION, INTRON)` and to nothing else**, so on an unstranded library
  the relay-only set is exactly *exon regions plus contiguous boundaries*: 27.6 % + 20.6 % of mass off capture,
  53.3 % + 44.6 % on it. ⚠ **Boundaries are the smaller half** — the earlier framing of this as a
  boundary-axis property was wrong.
* ⭐⭐ **100.0 % of the relay-only mass sits on slots that HAVE a count and a gDNA opportunity.**
  Relay-only never meant "no information"; it means no channel is wired to read the information present.
  ⛔⛔ **AND "THE RELAY CARRIES 93.1 % OF THE ERROR" IS NOT AN INDICTMENT OF MESSAGE PASSING — READ IT WITH
  §0c.** `relay-only` names the objects that have no own observation to be right or wrong FROM; §0c proves
  that an exon is structurally one of them, so this figure locates the error at exactly the population the
  relay is the only mechanism for. ⛔ Reading it as *"the relay is what breaks it, so mute it"* is the
  re-derivation §0c exists to stop.
  ⛔ So the coverage of the density deconvolve is a *frame and region-type restriction*, not a limit of what is
  knowable — and `EQUATIONS.md` §3.2 (density is frame-invariant to ~0.36 %) is the argument that ρ_bg
  should transport to the crossing frame at all. What it cannot survive is non-uniform gDNA placement,
  which is precisely what capture creates at exons (TRAPS: capture-is-1000x-on-exons).
⛔⛔ **AND ALL OF THE ABOVE SCORES HONEST IGNORANCE AS ERROR, WHICH IS THE WRONG QUESTION FOR PASS-0.**
The prior-free solve exists to produce a **substrate the gDNA hyperprior is fitted against**, not to be
accurate everywhere; an object with no own evidence reporting `f_g ≈ ½` at *zero precision* is stating a
true fact about itself. Excluding the undetermined population, pass-0's error on the objects it CAN
solve is **0.0456** where the all-objects figure said 0.3150 — 99.5 % of the difference is the
undetermined class. **The measurement that matters is solvable → right / wrong → confidently wrong**,
because a wrong value with a tight variance outvotes correct neighbours and anchors the prior
(`scripts/design/solvability_audit.py`; the retired `node_error_attribution.py` and
`confident_fp_trace.py` had this before the refactor and it was lost).

⚠ On real cfRNA most confident-gDNA regions have **zero** counts (64–94 % across libraries), and genome-wide
80.5 % of regions carry no fragments at all. A density-space estimator floors at `1/E` and discards most of
the evidence.
