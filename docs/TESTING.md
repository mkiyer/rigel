# TESTING — the substrate, the gates, and what they can judge

**Two rules govern this whole file**, and both came from real failures (TRAPS: prove-the-substrate, TRAPS: can-the-benchmark-resolve-it):

> ⛔ **Prove the substrate before you prove the code.** When a simulated axis is the axis you are judging,
> gate the simulator on it.
>
> ⛔ **Prove the suite can resolve the axis you are changing, before quoting a number from it.**

---

## 0. ⭐⭐ THE RANKING PANEL, TWO SIDE PANELS AND A TOY HARNESS — THEY ANSWER DIFFERENT QUESTIONS

⛔⛔ **`pilot/`, `flgap_short/` and `flgap_long/` WERE DELETED ON 2026-08-13 (owner) AND THE LADDER WAS
REBUILT FROM SCRATCH AT 16 CONDITIONS.** ⚠ **Three panels are on disk** under
`~/Downloads/rigel_runs/suite/`, beside the reference and index that were kept: `ladder/`, the only panel
the tool is RANKED on, and the two fl-gap SIDE panels `flgap_rna_long/` and `flgap_rna_short/`, built
under those new names on 2026-08-21 (5.3 GB and 5.7 GB, measured 2026-08-22).
⛔⛔ **The two side panels carry a DIFFERENT nascent model from the ladder** — read the fl-gap section
below before comparing any number across them. ⚠ **The `pilot/` column below is retained as a RECORD of
the Stage-A question no panel currently answers** — its config file survives, nothing is simulated from
it — because deleting the column would hide that the question went away with it.
⭐⭐ **And §0b's TOY HARNESS is the other substrate: the panel says how much error there is and where, a
toy says WHY.** Localise on the panel, isolate on a toy, re-measure on the panel.

| | ⛔ `pilot/` — 8 conditions, DELETED 2026-08-13 | ⭐ `ladder/` — 16 conditions |
|---|---|---|
| **judges** | **Stage A**, the accumulator | **Stage B**, the calibration solver |
| fragment lengths | REALISTIC and DIFFERENT per origin (RNA 206 ± 98, gDNA 157 ± 125) | ⭐ **IDENTICAL for both origins** (206 ± 98, [50, 500], read 100) |
| why that way | a length model can only be judged where the components' lengths actually differ | ⭐⭐ **because the EM ALREADY USES THE FL DISTRIBUTION** — a large gDNA-vs-RNA length gap lets it split the two origins on LENGTH ALONE, BYPASSING calibration and MASKING bugs in it. Equal lengths FORCE calibration to be exercised. See below |
| gDNA axis | {0, 100 %} rate ⇒ f_gdna ∈ {0, 0.5} | ⭐ **4 rungs, f_gdna 0 / 0.05 / 0.50 / 0.98** |
| depth | 10 M RNA, gDNA **added on top** (so 10–20 M total) | ⭐ **10 M TOTAL, fixed**; the rate decides only the SPLIT |
| nascent RNA | `nrna_none` on every condition | ⭐⭐ **SPARSE on every row (2026-08-22)**: `mode: sparse`, `on_fraction 0.50` of gene SPANS on, level ~ **logU(1, 100)** absolute and INDEPENDENT of the mature level, so `nascent > mature` really occurs — **23.6 % of ON spans carry more nascent than the SUMMED mature abundance of their contributors** (measured 2026-08-22 off the index and this panel's own `truth_abundances_nrna_mid.tsv`). Realised: **1,916 of 3,662 ELIGIBLE spans on** (eligible = named by at least one expressed multi-exon transcript), **molar ratio 0.00895, 20.2 % of RNA fragments** off capture — the same total burden as the retired uniform 20 % model, distributed sparsely, so a behaviour change between the two is attributable to the DISTRIBUTION and not to the amount. ⛔⛔ **0.50 is a DEVELOPMENT STRESS level, not real data** (`DESIGN.md` §0b, THE NASCENT SCOPE RULING); a realistic setting is `on_fraction 0.10` ⇒ 4.2 % |
| config | `scripts/sim/configs/pilot.yaml` | `scripts/sim/configs/gdna_ladder.yaml` |

### ⭐⭐ The fl-GAP SIDE PANEL — two arms, opposite signs (added 2026-08-21)

⛔⛔ **A SIDE PANEL, NEVER A LADDER RUNG.** The ladder equalises fragment lengths deliberately (the row
above says why), so it is a NULL for anything whose mechanism is a gDNA-vs-RNA length difference. This
panel exists to exert exactly that one mechanism, and it does not replace the ladder for ranking
anything.

| | `flgap_rna_long/` | `flgap_rna_short/` |
|---|---|---|
| config | `scripts/sim/configs/flgap_rna_long.yaml` | `scripts/sim/configs/flgap_rna_short.yaml` |
| gDNA fl configured | **75 ± 20** | **250 ± 60** |
| RNA fl configured | **250 ± 60** | **75 ± 20** |
| ⭐ gDNA / RNA **MEASURED** | **78.58 / 247.62** | **249.59 / 78.43** |
| ⭐ **realised GAP** | **+169.04 bp** | **−171.15 bp** |
| conditions | `g50` × ss {0.50, 0.99} × capture {off, on} — 4 | the same 4 |
| nascent RNA | ⛔ the **RETIRED UNIFORM** model — `mode: fragment_share`, `shares: [0.20]` | ⛔ the same retired uniform model |

⛔⛔ **BOTH SIDE PANELS PREDATE THE SPARSE NASCENT REBUILD AND STILL CARRY THE RETIRED UNIFORM MODEL, AND
A READER MUST KNOW THAT BEFORE COMPARING EITHER ONE AGAINST THE LADDER.** Their configs ask for
`nrna.mode: fragment_share` at a flat `shares: [0.20]` — nascent RNA on EVERY expressed multi-exon span,
at one solved molecular level — while the ladder was regenerated on 2026-08-22 with `mode: sparse`, where
half the spans carry none at all (the ladder table's nascent row, above). ⭐ **Each panel's data on disk matches its own config, so
both are internally accurate**: this is a difference in what was simulated, not a stale file, and every
measurement already taken on either arm stands. ⚠ "Retired" here means retired as a PANEL's nascent model,
not removed — `fragment_share` and `additive_ratio` are still live modes in `src/rigel/sim/`, so these two
configs run exactly as written. ⛔ What it costs is the CROSS-PANEL reading: a
ladder-vs-side-panel difference confounds the fl gap with the nascent DISTRIBUTION, and the fl gap is the
one thing these panels exist to vary — so name the panels behind any such number, or re-simulate first.
⚠ **Whether to re-simulate them is an open OWNER decision and has not been made.** It is ~11 GB of output,
and the mechanism the arms probe — `E_r/E_g − 1` on a MIXED pool — is exerted by the fl gap rather than by
the nascent distribution, so the arms remain usable for the question they were built for.

⭐ **WHY BOTH DIRECTIONS.** RNA is **not** reliably longer than gDNA — true for cfRNA, false elsewhere —
so a one-sided panel lets the tool overfit to one library type. And the quantity under test,
`E_r/E_g − 1`, **flips sign between the arms**, so a real repair must move the two in OPPOSITE
directions; one that merely trades one bias for another will not.

⚠ **THE STD MOVES WITH THE MEAN, AND THAT IS FORCED.** The sampler is a rejection draw from
`Normal(mean, std)` truncated to `[frag_min, frag_max]`, so a configured mean of 75 at the ladder's std
of 98 **realises at 137.65 bp** — truncation dominates and the arm is not the arm. ⛔ So the configured
"75 / 250" is a CONFIGURATION and the gap is a MEASUREMENT: read it off the partitions' own deposit
histograms, which `total_abundance_audit.py` prints on every row. ⚠ The RNA side lands SHORT of its
configuration (247.62 against 250, and 78.43 against 75 in the other arm) because an mRNA fragment must
fit inside its transcript while gDNA need not — the same truncation the ladder records for its equal
arms. ⭐ The two arms come out near mirror-symmetric, which is what makes the sign-flip test clean. ⚠ Two different "gaps" are quotable
and they differ — the raw deposit histograms give the ladder **−6.79 bp** while the fitted (de-tilted,
EB-shrunk) fl pmfs give **−12.38 bp**. The histogram describes the LIBRARY; the pmf describes what the
tool BELIEVES. Say which one you mean.

⚠ Reads are capped at the fragment length by the simulator, so a 100 bp `read_length` against a 78 bp
mean is safe (no adapter run-off); 88 % of the short arm is under 100 bp.

⛔ **WHAT MAY AND MAY NOT BE READ OFF IT.** The TRANSCRIPT-level number is not a calibration result here
— the EM reads the fl distribution and a gap hands it the answer. Everything that stops BEFORE the EM is
valid: `total_abundance_audit.py` (banks only) and `calibration_vs_oracle.py` (it runs `calibrate`, never
the EM). ⭐ The one END-TO-END row that IS valid is the LIBRARY gDNA FRACTION, which is a composition
deliverable rather than an EM assignment.

### ⭐⭐⭐ AND THERE IS A THIRD fl PANEL — ON THE TEST CHROMOSOME, ALREADY ON THE CURRENT NASCENT MODEL

⛔⛔ **READ THIS BEFORE RE-SIMULATING THE TWO SUITE ARMS ABOVE.** The fl gap is *also* exerted on the
method-development test chromosome (§0a), by three panels that are fully simulated and cached **30/30 scan
and oracle** and that were documented nowhere until 2026-08-30:

| | `scenarios_fl_rna_long/` | `scenarios_fl_gdna_long/` | `scenarios_fl_equal200/` |
|---|---|---|---|
| config (`scripts/sim/configs/`) | `test_reference_fl_rna_long.yaml` | `test_reference_fl_gdna_long.yaml` | `test_reference_fl_equal200.yaml` |
| gDNA / RNA fl configured | 100 / 250 | 250 / 100 | 200 / 200 — ⭐ **the CONTROL** |
| conditions | the full 30 (`g00…g98` × ss `0.50/0.70/0.99` × capture) | the same 30 | the same 30 |

⭐⭐ **They carry the CURRENT nascent model, and that is structural rather than lucky**: the test reference
is `abundance.mode: file`, so each config's `nrna:` block is **dead by design** and nascent comes from
`test_abundances.tsv` — the same file the main test panel reads. ⛔ **So the suite arms' stale nascent model
does NOT block work on the fl gap**, and treating it as a blocking gate has already queued ~11 GB of
unnecessary re-simulation once. Verify rather than trust the config: the condition names carry `nrna_file`,
and `fl_pool_purity.py` shows nascent in the intronic pool and *mature* fragments in the intergenic pool
(the shadow transcripts, §0a).

⭐ **The equal-200 arm is the control that makes the other two a measurement** — an fl-gap result that also
moves on the equal-length arm is an artefact, not the length gap.

### The gDNA ladder, in detail

`f_gdna = rate/(1+rate)`, and the RNA share thins as gDNA rises because the **total** is held:

| rung | g00 | g05 | g50 | g98 |
|---|---|---|---|---|
| rate | 0.0 | 0.052632 | 1.0 | 49.0 |
| f_gdna | 0 | 0.05 | 0.50 | 0.98 |
| n_rna | 10.0 M | 9.50 M | 5.00 M | 0.20 M |

⚠ **The rebuild dropped `g01`, `g10`, `g25`, `g75` and `g90`** (2026-08-13). Three things fixed which
four rungs survived, and none of them is "spread them evenly": **`g00`** is the owner-required
zero-gDNA control, **`g98`** is the top of the range, and ⛔ **`g05` is there because an INSTRUMENT
demands it** — `suite_resolves.py`'s requirement (f) is `0 < gdna_rate <= 0.10 WITH capture on`,
degenerate value 0 conditions, and a `{g00, g50, g98}` panel scores the degenerate value exactly. ⚠ The
test is on the RATE, not on `f_gdna`: the retired `g10` had rate 0.111 and would not have satisfied it
either. ⭐ Three levels is itself a floor, not a preference —
`TRAPS: a-single-level-panel-cannot-see-a-constant`.

× strand {0.50, 0.99} × capture {off, on} = **16**. Owner ruling: real libraries run from almost zero
gDNA to **> 98 %** of all fragments and Rigel must be robust across all of it; 10 M total per condition
is the ceiling, and the RNA-side abundance accuracy that thins at the top rungs is an accepted
trade-off — it is a property of such libraries, not an artifact.

### ⭐⭐ The four strata — THREE are in scope for 0.8.0 and ONE is DEFERRED

⛔⛔ **OWNER RULING, 2026-08-14, AND IT SUPERSEDES THE EARLIER "the focus condition is `ss_0.50 +
capture_on`"** (owner, retained here as the record of what the focus used to be). `ss_0.50` is the
UNSTRANDED half of the strand axis, and that cell is the deferred one:

| stratum | 0.8.0 |
|---|---|
| unstranded (`ss_0.50`) × capture **off** | ⭐ **IN SCOPE** — a development target |
| stranded (`ss_0.99`) × capture **off** | ⭐ **IN SCOPE** — a development target |
| stranded (`ss_0.99`) × capture **on** | ⭐ **IN SCOPE** — a development target |
| ⛔ unstranded (`ss_0.50`) × capture **on** | ⛔⛔ **DEFERRED** — not a development target until the other three are fully optimised |

⛔ **DEFERRED IS NOT DROPPED, AND THE DIFFERENCE IS THE WHOLE POINT.** The deferred stratum REMAINS in
every benchmark and every measurement and must keep being REPORTED on every table; it simply is not what
work is aimed at. If it improves as a side effect of the other three, that is a free win.

⛔⛔ **AND IT IS WHERE THE ERROR IS, SO A POOLED PANEL NUMBER IS MOSTLY A STRATUM NOBODY IS OPTIMISING.**
Measured 2026-08-13/14 on the rebuilt 16-condition ladder (noise floor 0.996–1.013): unstranded ×
capture-ON carries **64.5 % of transcript error and 90 % of gene-level error**; the three in-scope strata
carry the rest. ⛔ On it the tool emits a near-zero gDNA fraction REGARDLESS of truth — exon `f_g`
**0.040 / 0.0016 / 0.0021** at `g05` / `g50` / `g98` against truths **0.054 / 0.518 / 0.982** — so it
looks acceptable at low gDNA **by coincidence**, and a low-gDNA pass there is not evidence of anything.
⭐ Report per stratum, never pooled (TRAPS: a-single-level-panel-cannot-see-a-constant is the same
disease one axis over).

⚠ The `g00` rung is the false-positive check on every stratum, in scope or not. ⛔ Do not tune on a
control.

### ⭐⭐⭐ WHY THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS

⛔⛔ **THE REASON IS NOT THE ONE THIS FILE USED TO GIVE, AND THE WEAKER VERSION IS WHY THIS HEADING
EXISTS** (owner, 2026-08-14). It used to say the length channel was "neutralised, so residual error is
attributable to density and strand". That is a *consequence*, and stated alone it reads as tidiness.

⭐⭐ **The reason is that the EM ALREADY USES THE FRAGMENT-LENGTH DISTRIBUTION.** Give gDNA and RNA
visibly different lengths and the EM can assign fragments on **LENGTH ALONE** — it separates the two
origins downstream, **BYPASSING calibration entirely**, and a panel built that way scores well while
calibration is broken. ⛔ **A length gap MASKS calibration bugs.** Equal lengths remove that shortcut and
**FORCE the calibration phase to be exercised**: with the two origins length-indistinguishable, the only
channels left are **density** and **strand** — plus belief propagation across objects, which is currently
off (`CLAUDE.md`, message propagation).

⚠ **This is about the EM's use of the FL pmf, which is real and shipped.** The fragment-length channel as
a **CALIBRATION composition channel** is a different thing: it does **not** exist in `src/`, it was A/B'd
once and never shipped, and it is **DEFERRED POST-0.8.0** by owner ruling — do not read "the length
channel is neutralised" as a claim that calibration has one. ⚠ `length_likelihood` in
`src/rigel/second_pass.py` is per-fragment second-pass assignment and is a third, unaffected thing.

⛔ **"EQUAL" IS CONFIGURED, NOT ACHIEVED — AND THE RESIDUAL IS MEASURED, NOT ASSUMED.** Identical
parameters still leave a realised gap of **+3.56 to +4.09 bp off capture and +10.80 to +11.26 bp under
it** (re-measured 2026-08-22 with `simulator_gates.py` over all 12 gDNA-bearing ladder conditions;
TRAPS: configured-lengths-are-not-realised: a mature fragment must fit inside its transcript, gDNA need
not, so transcript truncation pulls RNA down). ⚠ **This line used to read `+4.68 / +3.57`, and the
capture-ON half was stale by ~3×** — the discrepancy is not attributable to the sparse nascent rebuild,
since capture's length selection did not change. ⭐ Nascent fragments run **~4.7 bp longer than mature**
(216.4–216.7 against 211.7–212.3 off capture) because a nascent span is not transcript-truncated, so the
RNA mean is a MIXTURE whose weights the nascent model sets — one more reason to quote the gap per
population rather than pooled. What
that residual is WORTH was measured before it was accepted, one thing varied: replacing both pmfs by
their pooled average moves the per-object error by **−0.0002** off capture and **−0.0054** under it —
under 2.5 % of the error, with a *fitted* gap 3–5× larger than the residual. Owner ruling: carry it.

⭐ **An ORACLE CACHE sits beside the panel** at `ladder/oracle_cache/`, built by `panel.py cache`.
⛔⛔ **ITS CONDITION COUNT IS NOT STABLE AND MUST NOT BE READ AS EVIDENCE IN EITHER DIRECTION — this
paragraph has now asserted `12/16` and `16/16` in turn, and both were true on the day they were written.**
Re-derived 2026-08-22 on the rebuilt panel: every condition holds all **FIVE** partitions — the three
`ORIGINS` (`gdna`/`mrna`/`nrna`) plus the per-strand `rna_pos`/`rna_neg` pair — and a stamped
`slot_truth.npz`, while the undrained `_main` payload is present on **12 of 16**: the four `g00` rows have
none. ⛔ `panel.py status` counts `(gdna, mrna, nrna, _main)`, so it reads **oracle cache 12/16** and
prints **✘ on a complete panel**. Two facts pull against each other and a reader needs both:

* ⭐ `pass0_vs_oracle.py`'s own sweep **HOLDS OUT every zero-gDNA condition** ("N zero-gDNA row(s) held
  out as false-positive checks") and builds no cache for one — so a panel cached only by that route
  legitimately reads **12/16**, and `panel.py status` then prints a **✘ on a complete panel**, because it
  counts conditions and the hold-out is invisible to it. ⛔ That ✘ is a false alarm, not a half-built
  panel, and it has cost a session. ⚠ It is also why the four `g00` rows on the rebuilt panel carry five
  partitions and no `_main`: their partitions came from the per-condition prewarm §2 documents, which
  fills a cache without scoring the row.
* ⛔ But `pass0_vs_oracle.measure_condition` **WRITES** `<oracle_cache>/<condition>/_main` whenever an arm
  runs with `--oracle-cache`, and other instruments do sweep the `g00` rows — so any such run fills the
  missing four and the count goes to 16/16. **A `_main` beside a `g00` row therefore proves nothing about
  what has been measured there** (`TRAPS: shard-an-arm-sweep-by-condition`).

⭐ Either way nothing is lost: `quant_accuracy.py`'s `base`, `base_reseed`, `warm_uniform` and every
`alloc_*` arm do not load an oracle at all, so `g00` scores normally; only the `oracle_*` arms need one,
and those are the arms a zero-gDNA condition has no truth-split to feed. It holds
the per-origin split payloads, so `--oracle-cache` turns a 4-minute-per-condition oracle build into
seconds. ⭐⭐ **It stays valid across every CALIBRATION change** — the oracle depends only on the
accumulator and the index — so one cache serves an entire solver-debugging campaign. It is keyed by the
scan cache's own key, so a stale one is refused rather than silently used. ⛔ Rebuild it after any
**accumulator** change (delete the directory; it repopulates).

### ⭐⭐ EVERY SCENARIO MUST BE CACHED — a requirement of the 0.8.0 loop, not an optimisation

⛔ **Owner ruling, 2026-08-14.** The focus of development is **CALIBRATION**, so the loop is *change
calibration → re-run calibration → score it against oracle calibration*, and **nothing in that loop may
re-scan a BAM**. A scan is minutes per condition and an oracle build was 4 minutes; a calibration re-run
off a cache is seconds. Two caches carry it and ⭐ **`panel.py cache` builds BOTH**:

| cache | holds | invalidated by |
|---|---|---|
| **scan** (`build_scan_cache.py`, `panel.py cache`) | the accumulator payload — scan once, calibrate many times | ⛔ any **accumulator** change (`payload_schema_digest`), the **index** (`graph_hash`, `reach_digest`) |
| **oracle** (`ladder/oracle_cache/`) | the per-origin split truth every scoring instrument reads | ⛔ the **accumulator** or the **index**, and nothing else |

⭐⭐ **Neither is invalidated by a CALIBRATION change** — that is precisely what makes a one-second
calibration iteration possible, and it is why one cache serves an entire solver-debugging campaign.
⛔ A scenario without both caches is not usable for development: `panel.py status` names which are
missing, and an instrument fed a stale one REFUSES it rather than silently using it.
⚠ The toy harness's donor bundle is the deliberate exception and for the opposite reason (§0b): it is a
function of the calibration code that fit it, so caching it would serve a stale answer.

---

## 0a. ⭐⭐⭐ THE METHOD-DEVELOPMENT TEST CHROMOSOME — where the message policy is developed

⭐⭐⭐ **A 1 Mb chromosome the owner designs one structure at a time.** Owner: *"over the course of
this method development, we will add transcripts to a 'test chromosome' … collectively, the test
reference chromosome, test transcript GTF, and test abundances comprise a critical method development
benchmark."* It was cleared to zero on 2026-08-27 (the 42 transcripts it had accumulated are in
`git log -- scripts/sim/test_reference/test_chr.gtf`) and rebuilt on 2026-08-28 as:

⭐⭐ **THE ANCHORED TWIN BLOCK** — ONE shape everywhere: exon 1 kb · intron 7 kb · exon 1 kb ·
intron 7 kb · exon 1 kb (span 17,000 bp), 20 kb intergenic between genes, all `+` strand. FIVE gene
types varying **exactly one bit each**, repeated in FIVE abundance blocks = **25
transcripts**, plus **8 SHADOW transcripts** on the blank contig (unannotated transcription the
simulator draws from and the index never sees). ⭐ The five types ARE the message layer's controls,
which is why this substrate is where the policy is developed:

| type | probed | nascent | what it controls |
|---|---|---|---|
| `clean` | no | no | intron + flanks **exactly pure gDNA** — the clean message SOURCE; one certified sj per intron |
| `nasc` | no | yes | the **contaminated** source — the intron is not pure gDNA |
| `cap` | yes | no | the **enrichment cliff** between source and destination |
| `capnasc` | yes | yes | both at once |
| `silent` | yes | `mrna = 0` | every object pure gDNA — **any claimed RNA is a false positive** |

Abundances are molar ladders in half-decade steps with **mature UP the blocks and nascent DOWN**
(10/30/100/300/1000 against 100/30/10/3/1), i.e. independent levels, exactly as the big ladder draws
them. ⛔ `test_abundances.tsv` takes no comments — both its readers key on the first line — so that
design is recorded in `test_chr.gtf`'s header, which is the file to read before editing either.

⭐⭐ **WHY THIS SUBSTRATE EXISTS: the whole loop is SECONDS.** One calibrate on the cached test
chromosome measured **0.04 s** (74 slots), against minutes on the ladder — so a policy change can be
scored on every condition between edits, and the ladder is kept for the shipping judgement.

| piece | |
|---|---|
| `scripts/sim/test_reference/test_chr.gtf` | ⭐ **HAND-EDITED** — one `exon` line per exon. This is where a transcript is added |
| `scripts/sim/test_reference/test_abundances.tsv` | ⭐ **HAND-EDITED** — `transcript_id / mrna_abundance / nrna_abundance`, relative molecular sampling weights |
| `scripts/sim/test_reference/test_probes.bed` | ⭐ **HAND-EDITED** — BED12 capture probes, or empty for no probes |
| `scripts/sim/test_reference/test_shadow.gtf` | ⭐ **HAND-EDITED** — ⭐⭐ **SHADOW transcripts on the BLANK chromosome `test_blank`** (owner design, 2026-08-29): unannotated transcription the SIMULATOR draws from (`shadow_gtf:` in the panel config) and the INDEX never sees. Their abundances live in `test_abundances.tsv` like any transcript (nascent always 0 — a shadow IS the unannotated RNA); a shadow id the annotation declares is REFUSED by the simulator. What it controls: anything that trusts the annotation's notion of "pure gDNA" — the gDNA strand-overdispersion fit above all — must be unmoved by RNA the tool cannot know about |
| the FASTA (⚠ `test_chr.fa` now carries BOTH contigs), the index, the reads and the caches | ⛔ **DERIVED** — never hand-edited; regenerated by the commands below |

⛔ **WHY THE FASTA IS DERIVED.** A spliced transcript needs a GT..AG at every intron or the aligner and
the simulator disagree with the annotation. Generating the chromosome from a fixed seed (1,000,000 bp,
`test_chr`) and injecting a motif at every intron the GTF declares makes that impossible to get wrong; a
versioned FASTA would need a hand edit on every addition and the first missed one is a silent failure.

⛔ **THE BUILDER REFUSES A REFERENCE IT CANNOT SIMULATE** and reports every problem at once: an exon off
the end of the chromosome, overlapping exons within a transcript, an intron under 4 bp (no room for a
motif), a duplicate id, a transcript with no abundance row, an abundance row for a transcript that is
not in the GTF, the wrong reference name. `--self-test` 16/16, each check perturbed.

⚠ **κ IS FITTED FROM SPLICED READS**, so a transcriptome of only single-exon transcripts cannot
calibrate — `calibrate` raises `CalibrationStrandError`. At least one multi-exon transcript carrying
real depth is needed before any benchmark number means anything. (This blocks the STRAND-MODEL FIT, not
calibration itself: given κ, single exons calibrate fine.)

⭐⭐ **NASCENT RNA NEEDS NO DECLARATION** — `rigel index` creates a single-exon nascent ENTITY spanning
each multi-exon transcript, which is exactly the owner's *"just add a single-exon transcript spanning the
multi-exon transcripts"*. Give it abundance through its CONTRIBUTOR's `nrna_abundance` column; the
simulator pools it onto the entity. ⚠ An ANNOTATED single-exon transcript is already its own nascent
equivalent (`is_nrna = True`, no synthetic made).

### ⭐⭐⭐ THE SWEEP — 30 CONDITIONS (owner, 2026-08-27)

`configs/test_reference.yaml` sweeps the chromosome over three axes, and two more requirements are
carried by the SUBSTRATE the owner authors rather than by the config:

| axis | values | where it lives |
|---|---|---|
| **gDNA fraction** | `g00` `g05` `g25` `g50` `g98` — 0 %, 5 %, 25 %, 50 %, 98 % | the config (`rate = f/(1−f)`) |
| **strand specificity** | `0.50` unstranded · `0.70` · `0.99` strand-specific | the config |
| **capture** | off · on | the config |
| **which transcripts are captured** | ⭐ about **HALF** probed, half with NO probe | `test_probes.bed` |
| **nascent RNA** | ⭐ SPARSE — up to about **HALF** the transcripts carry it, the rest exactly none | the `nrna_abundance` column |

5 × 3 × 2 = **30 conditions.** ⭐ `g00` is the zero control (truth is exactly 0, so every gDNA fragment
reported is a false positive with nothing to cancel it) and `ss 0.70` is the transition rung — neither
of the two regimes the bars are written for, and therefore the rung that says whether a policy degrades
gracefully between them.

⛔ **BOTH SUBSTRATE REQUIREMENTS EXIST SO NO CONDITION IS UNIFORM.** A library where every transcript is
probed, or where every intron carries nascent RNA, tests the policy on a world it will never meet; the
half-and-half split puts captured and uncaptured structures, and nascent-free and nascent-bearing
introns, inside every single condition. ⭐ Vary the nascent LEVEL as well as its presence — below, near,
and above mature — and hand-author the pattern rather than drawing it: on a chromosome this small a
Bernoulli draw is mostly sampling noise, and the point of a debug substrate is that you know WHICH
structure carries nascent RNA before you read a number off it.

⚠ **THE `nrna:` BLOCK IN `configs/test_reference.yaml` IS DEAD BY DESIGN.** `abundance.mode: file` makes
the TSV the one source of both weights, so the ladder's `mode: sparse` draw never runs here. Editing
that block changes nothing.

### ⭐⭐⭐ THE COMMANDS — from a hand-edited GTF to a scored benchmark

Everything derived is rebuilt from the three hand-edited files. Run in order; each stage is resumable
and `panel.py status` names the next one.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

# 1. FASTA + index, derived from the GTF (refuses a reference it cannot simulate)
python scripts/sim/build_test_reference.py

# 2. simulate every condition, then build BOTH caches (scan + oracle, with certified slot truth)
python scripts/sim/panel.py status   --config scripts/sim/configs/test_reference.yaml
python scripts/sim/panel.py simulate --config scripts/sim/configs/test_reference.yaml
python scripts/sim/panel.py cache    --config scripts/sim/configs/test_reference.yaml

# 3. score every policy on every condition — seconds
python scripts/design/policy_benchmark.py --panel test
```

⛔ **REBUILD AFTER EVERY EDIT TO THE THREE SOURCE FILES.** The FASTA, index, reads, caches and
`slot_truth.npz` all describe the annotation they were built from; a benchmark run against stale caches
answers a different question and says nothing about it (`TRAPS: a-green-suite-hid-five-dead-instruments`
is the same failure in another dress). ⚠ The 42-transcript derived artifacts were moved to
`~/Downloads/rigel_runs/test_reference_STALE_42tx_2026-08-27/` on the same day rather than left in place
for exactly this reason.

### ⭐⭐ HOW TO READ THE BENCHMARK — two halves, two different bars

⛔ **THE BAR IS NOT "BEAT SILENT".** `SilentPolicy` is the measured floor, and on strand-specific data
it is very hard to improve on: a sighted exon's own strand solve is excellent, so a message can mostly
only disturb it. Message propagation exists for UNSTRANDED data, where the strand channel is dead and
the local answer is a default rather than a measurement.

| rows | the bar |
|---|---|
| **unstranded** (`ss 0.50`) | ⭐ the policy must **WIN** — this is what the message layer is for |
| **stranded** (`ss 0.99`) | ⭐ the policy must do **as little harm as possible** — near 1.00× of silence is a pass |

⛔ **NEVER POOL THE TWO HALVES**, and never pool conditions within them: a total hides a sign flip
between strata, and the halves are being judged against different bars.

⚠ **A TOY AND THE PANEL CAN DISAGREE IN RANK** (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`): a
small substrate's objects mostly lack own evidence, so a message is nearly free there and a policy can
look better than it is. Develop on the test chromosome; confirm on the ladder before believing it.

## 0b. ⭐⭐ THE TOY HARNESS — a mini chromosome you define, calibrated in under a second

`scripts/design/toy_harness.py`. ⭐ **This is the third substrate and the one to reach for when a
mechanism needs isolating below the level of a whole chromosome.**

⚠ **IT PIGGYBACKS ON A REAL CACHED CONDITION AND §0a DOES NOT.** A toy is too small to fit the
library-level quantities calibration needs, so the harness HARVESTS them from a panel condition
(`DonorGlobals`, ~30 s per session) — which means a toy result is conditional on a donor, and a donor
must exist. The test chromosome of §0a has no such dependency: it is a real simulated library of its
own, self-contained, and its whole calibrate is 0.04 s. ⭐ **Prefer §0a for developing a policy**; use
a toy when you need to isolate one mechanism on geometry you fully control. The panel above tells you *how much* error there is and *where*; a toy
tells you *why*, because it is small enough to read every object and cheap enough to sweep one variable
seven times.

⛔ **What it replaces.** Every defect found before it was found by sifting a 36-condition, 10-million
fragment panel for badly-behaved objects and reasoning backwards from a region id. That is slow, the
example is never quite the one you wanted, and a fix can only be measured in aggregate where two errors
cancel. The harness root-caused TRAPS: a-purity-filter-is-a-length-filter in **five objects and 0.1 s** after weeks of panel-level work.

### The one idea that makes it work

A toy cannot fit the library-level quantities calibration needs — one transcript has no population to
estimate a strand balance, an enrichment landscape or an intergenic background from. So **a real cached
condition acts as DONOR**: it is calibrated once, and its fitted bundle is injected
(`InjectedCalibrationPriors`, whose docstring already specified this use). The toy supplies only the
controlled per-region **geometry**, which is the thing under study.

| the donor supplies | so the toy never invents |
|---|---|
| κ, both strand overdispersions, the Fisher noise-floor sample sizes | the strand deadband behaves exactly as it does on real data |
| the intron background, and the pre-solve TOTAL-density landscape | a handful of regions cannot fit these. ⚠ The row used to name the *enrichment NPMLE* and *ρ_bg*: both were converge-and-deleted 2026-08-21, and the injectable field is now `abundance_landscape` |
| both fragment-length pmfs | passed as `calibrate` kwargs, not part of the priors bundle |
| capture on/off + its numeric knobs | reproduced in the toy's own **simulation**, with probes written from the spec |
| frag mean/sd/min/max, read length, strand specificity | read from the donor's post-capture truth |
| ⭐ **gDNA density per base** | see below — this one is not a knob |

⭐⭐ **THE gDNA LEVEL IS DERIVED, NOT CHOSEN, AND THAT IS NOT NEGOTIABLE.** The injected enrichment
landscape is an **absolute** log-density model, so a toy at the wrong depth is a *different* library, not
a small one. The harness measures the donor's gDNA counts-per-base on the donor's own
structurally-pure-gDNA regions (`Σcount / ΣE`, `EQUATIONS.md` §7.1) and simulates the toy to match it.
`ToySpec` therefore has **no gDNA field**, and a gate asserts it never grows one.

### ⭐ Terminology — one word per concept (owner, 2026-08-04)

| | |
|---|---|
| **counts** | discrete **integer** fragment counts. What the accumulator stores, and the solver's Poisson `n` |
| **density** = **abundance** | **counts per base.** The two words mean the same thing and are used interchangeably |

⛔ Not the simulator's molar `abundance=` field, which is a per-transcript weight, not a density.
⚠ And a region's stored counts are **contained** counts, so `counts = density × effective_length`, never
`density × bp`. The sweep reports the density it *asked* for and the counts each object *received*, side
by side, and never converts one into the other.

### Running it

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite

python scripts/design/toy_harness.py --list            # the spec ladder, simplest first

# one spec against one donor — prints EVERY object beside per-object truth   (~0.1-5 s)
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g50_ss_0.50_nrna_mid_capture_off

# ⭐ sweep the transcript's RNA density; the gDNA background stays PINNED, so one variable moves
python scripts/design/toy_harness.py --spec TA_single_exon \
    --donor gdna_g50_ss_0.50_nrna_mid_capture_off --sweep-density

python scripts/design/toy_harness.py --spec all --donor <cond>     # the whole ladder
```

⚠ Any of the 16 ladder conditions can be the donor, and **which one you pick is an experimental
variable** — it sets capture on/off, the strand regime and the gDNA level. ⛔ **It therefore selects a
STRATUM**: `ss_0.50 × capture_on` is the DEFERRED one (§0), so a mechanism isolated on that donor alone is
not aimed at a 0.8.0 target. Harvesting costs one scan +
one calibrate (~30 s) and is deliberately **not cached**: the bundle is a function of the calibration
code that fit it, so a stored copy would go stale on exactly the changes the harness exists to test.
⭐ Harvest once per session and run many toys against it (`harvest()` then `run_toy()` in a loop).

### ⭐⭐⭐ `spliced_exons` — the CURRENT TARGET RUNG, and the one to understand first

**Owner's spec: ONE gene, ONE transcript, `TA+ (1,000, 2,000) (9,000, 10,000)` on 12 kb.** Deliberately
`nested_exons`' **twin** — same chromosome, same gene boundaries — with an INTRON and an SJ where the
nesting was, so the two rungs differ by exactly one structure.

```
REGION intergenic [0, 1000)         BOUNDARY @1,000   intergenic|exon, pure gDNA (TSS+)   G1 object
REGION exon  [1000, 2000)   TA e1   BOUNDARY @2,000   intron|exon, the DONOR+ side        ⭐ NOT a G1 object
REGION intron [2000, 9000)  TA i1   BOUNDARY @9,000   intron|exon, the ACCEPTOR+ side     ⭐ NOT a G1 object
REGION exon  [9000, 10000)  TA e2   BOUNDARY @10,000  intergenic|exon, pure gDNA (TES+)   G1 object
REGION intergenic [10000, 12000)
SJ BOUNDARY 2,000 → 9,000 (+), pure mature RNA, ⚠ NOT a chain slot
```

⭐ **What it adds, and why it is the hard rung: the two exon↔intron BOUNDARIES.** Mature RNA cannot cross an
exon↔intron boundary contiguously (TRAPS: mature-rna-never-crosses-a-boundary), so their truth is pure gDNA — but the solver's own
continuity gate says a strand IS admissible there (RNA that has not spliced *there* could cross), so they
are **not** G1 objects (`DESIGN.md` §0) and the
solver has to *derive* what the structure already implies. ⚠ Every condition of the ladder this rung was
written against was `nrna_none`, so their truth was exactly 1.000 and the panel could not distinguish "no
RNA crosses" from "no *mature* RNA crosses". ⭐ The ladder has carried nascent on every row since
2026-08-19 and does distinguish the two. ⛔⛔ **But under the SPARSE model it distinguishes them only where
a span is switched ON, and that is a property to exploit rather than a limitation**: at
`g50 ss0.99 capture_off` the mass-weighted true intron `f_g` is **0.7019** — unmoved from the retired
uniform panel's 0.7045 — while **39.4 % of live intron slots now read true `f_g` EXACTLY 1.0000**, which
the uniform model gave at ~0 % of them. So the same panel now contains both faces of this rung at once:
introns with nascent RNA in them, and introns whose truth is the exact pure-gDNA 1.000. ⚠ `on_fraction
0.50` is a DEVELOPMENT STRESS level and not real data (`DESIGN.md` §0b, THE NASCENT SCOPE RULING), so read
a number from this face as a robustness check.

⭐ And unlike `nested_exons` there **is** own evidence inside the gene: the 7,000 bp intron REGION is where
the intron factory lives (τ ≈ 0.18–0.38, |Δf_g| 0.05–0.08), so the gDNA level need not travel from the
gene ends. The two exons have essentially **no** own evidence on an unstranded library (τ ≈ 5e−6) and are
carried entirely by messages.

⭐⭐ **The two FACES of an `intron|exon` BOUNDARY** — the derivation this rung exists to land is
`EQUATIONS.md` §3.6. Read it before touching the solver.

### ⭐⭐⭐ `splice_both_strands` — the rung the SPLICE-FLUX REFRAME must be derived against

**Owner's spec, 2026-08-05.** Four transcripts, both strands, overlapping exons AND overlapping introns,
two sj pointing opposite ways, on the same 12 kb chromosome:

```
TA+ (2,000, 3,000) (9,000, 10,000)     2 exons, + strand, intron 3,000–8,999
TB+ (2,000, 10,000)                    1 exon,  + strand, spans TA's intron
TC− (1,000, 11,000)                    1 exon,  − strand, spans everything
TD− (1,000, 2,500) (8,500, 11,000)     2 exons, − strand, intron 2,500–8,499
```

⭐⭐ **What it breaks, and it breaks it immediately.** Every earlier rung let a BOUNDARY ask "is my
neighbour an exon?" and get a yes or a no. Here the question has no answer — the 4-bit signature the
index stores is `{intron₊, intron₋, exon₊, exon₋}` and three regions carry both kinds at once:

| region | signature | what it is |
|---|---|---|
| [1,000, 2,000) | `0001` | exon₋ |
| [2,000, 2,500) | `0011` | exon₊ + exon₋ |
| [2,500, 3,000) | `0111` | ⭐ **intron₋** + exon₊ + exon₋ |
| [3,000, 8,500) | `1111` | ⭐⭐ **intron₊ + intron₋ + exon₊ + exon₋ — all four** |
| [8,500, 9,000) | `1011` | ⭐ **intron₊** + exon₊ + exon₋ |
| [9,000, 10,000) | `0011` | exon₊ + exon₋ |
| [10,000, 11,000) | `0001` | exon₋ |

⛔⛔ **And `coarse_type_array` reports EVERY one of them as `exon`** — exon wins over intron, and the
strand is collapsed. So on this rung the string `intron|exon` never appears, and any rule phrased on the
coarse region type is silent exactly where it is needed. ⭐ The information is not missing; it is discarded
one layer above the place that needs it.

⭐ **The sj axis and the boundary flags carry the rest of it, already plumbed:**

| | src → dst | genomic | strand |
|---|---|---|---|
| sj 0 | region 2 → region 5 | 2,500 → 8,500 | − |
| sj 1 | region 3 → region 6 | 3,000 → 9,000 | + |

| boundary | position | flag bits |
|---|---|---|
| #0 | 1,000 | `FLAG_TES_NEG` |
| #1 | 2,000 | `FLAG_TSS_POS` |
| #2 | 2,500 | ⭐ `FLAG_DONOR_NEG` |
| #3 | 3,000 | ⭐ `FLAG_DONOR_POS` |
| #4 | 8,500 | ⭐ `FLAG_ACCEPTOR_NEG` |
| #5 | 9,000 | ⭐ `FLAG_ACCEPTOR_POS` |
| #6 | 10,000 | `FLAG_TES_POS` |
| #7 | 11,000 | `FLAG_TSS_NEG` |

⚠⚠ **CHECK THE `_NEG` CONVENTION BEFORE RELYING ON IT.** Boundary 2,500 is flagged `DONOR_NEG` and boundary
8,500 `ACCEPTOR_NEG` — but on a **−** transcript the molecule runs right-to-left, so the biological
splice **donor** of TD−'s intron is at 8,500 and the **acceptor** at 2,500. Whether these bits mean
"genomic-low end of a − intron" or "the transcript's actual donor" decides the sign of the whole
derivation. ⛔ `EQUATIONS.md` §3.5b already ruled that this family of predicates must be written in
GENOMIC terms and never in transcript terms, and it is exactly where a sign gets flipped silently.

⭐ **The substrate verifier COVERS this rung, and this rung is where its own falsification lives.**
`verify_toy_substrate.py` takes any number of transcripts on either strand. Re-verified 2026-08-13 on the
rebuilt panel: every gate passes and **all six `--perturb` arms fire** — ⛔ **on `splice_both_strands`,
and the SPEC is part of that claim.** On the one-transcript, `+`-only `spliced_exons` rung only `drop_sj`
fires and the other five are SILENT, because there is no `−` fragment to mirror, no second geometry to
confuse, and the structural-set gate is vacuous there. Reading those silences as a pass is
`TRAPS: could-the-arm-have-fired`, and it was read that way once (TRAPS: self-checking-validator:
re-derive by a different algorithm before trusting).

### ⭐⭐ `toy_panel.py` — one spec × EVERY cached condition × an RNA-density ladder

`toy_harness.py` answers "what does this structure do on ONE library at ONE RNA level". Once a structure is
the **target** rather than a probe, you want the whole space:

```bash
# all 16 conditions x 7 RNA rungs, PRIOR-FREE pass-0, per-object.   ~13 s per condition
python scripts/design/toy_panel.py --spec spliced_exons --out rows.jsonl

# ⛔ capture-ON MUST be lengthened or it measures an empty chromosome
python scripts/design/toy_panel.py --spec spliced_exons --genome-length 120000 \
    --conditions $(ls $SUITE/ladder | grep capture_on)

python scripts/design/toy_panel.py --report rows.jsonl        # re-aggregate, no re-measurement
```

* **Prior-free pass-0 by default** (`--refit-iters 0`) — that is the substrate the gDNA hyperprior is later
  fitted against, so an object confidently wrong here anchors the prior. `--refit-iters 3` is the shipped
  solve; ⛔ the two answer different questions, do not quote one for the other.
* **The RNA density is a MULTIPLE of each donor's own gDNA density**, never absolute — the donors' gDNA
  rates span ~100×, so at rung `m` the exon's true `f_g` is roughly `1/(1+m)` on every row alike.
* Four tables: the gene's mwae per stratum × rung; **per object** the error, its share of the error MASS
  and whether the messages helped or hurt (`loc` vs `pred`); the sweep SHAPE (does the answer track the
  data?); and ⛔ **who is CONFIDENTLY wrong**, which is the population that corrupts a prior.
* Shard it with `--conditions`; harvesting is 30 s per donor and is deliberately not cached.

⚠ **Its per-object CEILING is a SUBSTITUTION, not a re-solve** — it replaces one object's answer with the
truth and re-scores. That is honest for a *sink* and **understates a message SOURCE**, whose value is what
it carries to its neighbours. For a source, run a real arm.

### Writing a new spec

Add a `ToySpec` to `SPECS` in the harness. The ladder is ordered simplest-first and **each rung adds
exactly ONE structure to the one before it** — that is the point: when a row goes wrong, the thing that
changed is the thing to look at.

```python
"my_case": ToySpec(
    name="my_case",
    what_it_probes="one sentence: which mechanism this isolates, and why this structure isolates it",
    genome_length=5_000,
    genes=[{"gene_id": "TA", "strand": "+",
            "transcripts": [{"t_id": "TA", "exons": [(1_000, 3_000)], "abundance": 100.0}]}],
    n_rna_fragments=1_000,      # the RNA knob; gDNA is pinned by the donor
    nrna_abundance=0.0,         # nascent RNA; the ONLY way an intron carries RNA
    captured=None,              # transcript ids to probe when the donor is capture-ON; None = all
    seed=7,
),
```

⚠ **There is no gene-free rung**: `TranscriptIndex` requires at least one transcript, so an unannotated
chromosome cannot be indexed. Use a **silent** gene (`abundance=0.0`) instead — which is a better first
rung anyway, since every object is then structurally pure gDNA and any deviation from `f_g = 1` is a
false positive with nothing to cancel against it.

### ⚠⚠ CAPTURE-ON needs care, and the harness will tell you when it is starved

Under capture the donor's **off-target** density is ~24× lower than its capture-OFF twin
(0.00106 vs 0.02566 counts/bp on `g25` — ⚠ a donor of the **RETIRED** 36-condition ladder, deleted
2026-08-13 and not re-derived on the 16-condition rebuild; see the STAMP below), and capture actively
*depletes* intergenic space. Three
consequences, all measured, and the sweep prints a ⛔ STARVED banner naming which one has bitten:

| starved object | does `--genome-length` help? | why |
|---|---|---|
| intergenic **region** | ⚠ barely | counts scale with bp at fixed density, but capture depletes off-target so hard that 57 kb yields **2 counts**. You would need megabases — which is exactly why the real panel is 93 Mb of chr21+chr22, and exactly why the donor **injects** `background` / `intron_background` / ρ_bg. ⭐ Under capture the toy's own intergenic regions are decoration; the background comes from the donor by design |
| an **boundary**, capture **OFF** | ⛔ **no, not at all** | a 0-bp boundary's counts are `density × mean_FL`, **independent of the chromosome length**. The only lever off capture is library depth |
| ⭐⭐ an **boundary**, capture **ON** | ⭐⭐ **YES — LENGTHEN IT** | and this is the OPPOSITE of the capture-OFF row, which is why the two are listed separately. The gDNA budget is `rate × genome_length` while the probe footprint is **fixed**, and the sampler's on-probe share is `binding·overlap / (off_target·L + binding·overlap)` — so a longer chromosome hands capture a bigger budget to concentrate onto the same probes, and a boundary's count grows with `L` until that ratio saturates |
| … and what else lifts a boundary? | ⭐ **probe geometry, then binding strength** | probes must **tile PER EXON** so each one abuts the exon↔intron boundaries and none straddles a sj (below). Raising `binding_per_base` also works but **un-matches the toy from the donor's chemistry**; lengthening keeps every harvested global intact, so it is the clean lever |

⛔ **STAMP: the three tables that follow were measured 2026-08-04 on the RETIRED 36-condition ladder, and
the donors `g01` / `g10` / `g25` / `g75` / `g90` NO LONGER EXIST** (§0 — the 2026-08-13 rebuild kept
`g00` / `g05` / `g50` / `g98`, and the panel is now 16 conditions, **8** of them capture-ON). They are
retained because the *mechanism* they establish is not panel-specific — capture moves the gDNA signal off
the intergenic and intronic REGIONs and onto the BOUNDARIES abutting an exon, and lengthening the
chromosome is the clean lever. ⚠ Re-measure the counts before quoting one; on the surviving rungs `g50`
and `g98` bracket the solvable end and `g05` / `g00` the starved end.

⭐⭐ **Measured, `spliced_exons` × `g75 ss0.50 capture_on`, 12 kb → 120 kb, same donor and same
chemistry** — this is the table that says capture-ON depletion is *the signal moving*, not starvation:

| object | 12 kb | 120 kb | gDNA density at 120 kb |
|---|---|---|---|
| `intron\|exon` BOUNDARY @2,000 | 2 | **20** | 0.0788 |
| `intron\|exon` BOUNDARY @9,000 | 5 | **36** | 0.1418 |
| `intergenic\|exon` BOUNDARY @1,000 | 1 | **41** | 0.1615 |
| exon REGION interiors | 16 / 12 | 121 / 118 | 0.162 / 0.158 |
| intron REGION (7,000 bp) | 1 | **1** | 0.00015 |
| intergenic REGION | 0 | 13 / 110 kb | 0.00012 |

⭐ So under capture the gDNA signal **leaves the intergenic and intronic REGIONs and arrives at the BOUNDARIES
abutting the exon**, which end up within 2× of the exon interior and ~1000× above the intron interior.
⛔ **That makes the capture-ON rung the one the `intron|exon` BOUNDARY actually matters in** — it is then a
well-counted object at the exon's own capture stratum whose component set excludes mature RNA.

##### ⛔⛔ …BUT THAT TABLE IS ONE CONDITION (`g75`), AND 20–40 COUNTS IS **NOT** TRUE ACROSS THE LADDER

Re-measured 2026-08-04 on all 18 capture-ON conditions × 7 RNA rungs at 120 kb. The `intron|exon` BOUNDARY's
count range, min–max over the rungs:

| donor | `@2,000` | `@9,000` | | donor | `@2,000` | `@9,000` |
|---|---|---|---|---|---|---|
| `g00` | 0–0 | 0–0 | | `g50` | 13–32 | 14–25 |
| `g01` | 0–1 | 0–2 | | `g75` | 20–35 | 24–36 |
| `g05` | 0–1 | 1–3 | | `g90` | 25–43 | 26–47 |
| `g10` | 2–6 | 1–8 | | `g98` | 29–44 | 29–50 |
| `g25` | 5–13 | 6–15 | | | | |

⭐ **`g50` and up are a solve; `g25` and below are not, and `g00` never can be — a zero-gDNA library has no
gDNA fragments to place.** ⛔ So `--genome-length 120000` earns the claim on **8 of 18** capture-ON
conditions, and any capture-ON aggregate over the whole ladder is averaging 8 solves with 10 empty
objects. Say which rows carried it.

⭐ **Lengthening still works at the low-gDNA end — it has not saturated by 1.08 Mb.** Same donor, same
chemistry, m=1, one variable moved (the RNA count is a multiple of the donor's own gDNA DENSITY, so it
does not depend on `L`):

| | 120 kb | 360 kb | 1,080 kb |
|---|---|---|---|
| `g25` `intron\|exon` @2,000 / @9,000 | 5 / 9 | **28 / 18** | **50 / 46** |
| `g05` `intron\|exon` @2,000 / @9,000 | 1 / 1 | 4 / 5 | 7 / 12 |
| ⛔ the intron REGION (7,000 bp), `g25` | **0** | **1** | **0** |

⭐⭐ **And the row that matters most is the last one: the intron REGION is dead under capture at EVERY
chromosome length.** Its bp is fixed and its density is off-probe, so lengthening cannot reach it — while
the BOUNDARY beside it grows without bound. ⛔ **The well-counted side therefore INVERTS with capture**
(TRAPS: capture-inverts-the-counted-side): off capture the intron holds ~315 counts against the BOUNDARY's 12–13, on capture it holds 1
against the BOUNDARY's 20–40.

⚠ **Two length slips in the harness itself, found on the way and worth fixing** (`_donor_sim_params`):
`frag_mean` is read from `truth_summary.json`'s **`all`** row — the gDNA+RNA MIXTURE — and handed to the
toy's RNA draw, which is +2.25 bp on `g50`; and that value is a **post-truncation realised** mean fed back
in as the **pre-truncation generating** mean, which the `[50, 500]` truncation then re-inflates by a
further +5.6 bp. Together they make the toy a ~7.8 bp longer-fragment library than its donor. ⭐ They are
small next to the solver-side gap they sit beside (`EQUATIONS.md` §3.6b, −8.4 bp) but they are the same
units and they add.

#### ⛔ PROBES TILE PER EXON, and the reason is a real suppression

Probes are written in **transcript** space, so a probe spanning an internal sj offset has a genomic
footprint in **two blocks** — and `sim/capture/sampler._split_scale` then multiplies every gDNA fragment
overlapping it by `gdna_split_penalty`. Tiling across the whole transcript put such a split probe over
every internal sj and thereby suppressed exactly the population that spans an `intron|exon` BOUNDARY.
`_toy_probes` now tiles **within each exon**, so every probe is unsplit and ends **on** the boundary.
⚠ The split-probe case is real in a real panel and deserves its own rung; it should not be an accident of
how a tiling loop was written.

### ⚠ What a toy CANNOT judge

* **Magnitudes do not transfer between donors.** The same `nrna` contrast is a factor of **31** against
  the `g25` ladder donor (⚠ RETIRED with the 36-condition ladder on 2026-08-13; the *contrast* is the
  lesson, and neither factor has been re-derived on a surviving donor) and **1.25** against a six-gene
  synthetic donor: direction preserved, size not.
  A small donor cannot determine the enrichment landscape or the intron background, so its injected
  globals are a different regime. ⛔ Gate on direction and ordering; quote magnitudes with their donor.
* **It cannot rank defects.** Five objects cannot tell you what fraction of a real library's error a
  mechanism owns. Localise on the panel, isolate on a toy, then re-measure on the panel.
* **It shares TRAPS: toys-rank-hotspots-backwards' warning** — a toy ranks performance hotspots backwards. It is a
  *correctness* instrument, never a profiling one.

Gates: `tests/calibration/test_toy_harness.py` — 7, each carrying its own perturbation, and the donor is
a scenario the gates **build** so none of them silently skips when the panel is absent.

---

## 1. The simulated panel (shared backbone)

**A real human backbone, not a generated mini-genome:** `chr21` + `chr22` + the **92 ERCC** spike-in
references, carved out of the same GRCh38 + GENCODE v46 sources the production index is built from.

⚠ **Neither the ERCC references nor chr21 is filler.** A single-reference synthetic index once hid a
reference-id-space mismatch that silently dropped 476,719 of 476,732 fragments inside `deposit()` while
every golden test passed (TRAPS: one-reference-hides-refid-bugs). 92 spike-ins make the reference-id space non-trivial **for
RNA** — but they are RNA-only, and gDNA takes a *different* branch through the scanner. **Two genomic
chromosomes is what makes the space non-trivial on the gDNA path too.**

### The grid

The owner's 8-condition minimum, a clean 2³: **gDNA {none, 100 %} × strand {0.50, 0.99} × capture
{off, on}**, nascent held at `none` so every `gdna_none` cell is a pure false-positive test.
Config: `scripts/sim/configs/pilot.yaml`. ⛔ **This grid is the DELETED pilot panel** (removed from disk
2026-08-13; the config file survives, but nothing is simulated from it). What is on disk is §0's
16-condition ladder.

⚠ **`10,000,000` is the RNA depth**, and gDNA is added **on top** at `rate × n_rna` — so a `gdna100`
condition is 20 M fragments and a `gdna_none` condition is 10 M. The depth is set so that fragments per
annotated base — and therefore per **region** — matches the chr22-only panel it replaced; calibration
accuracy depends on per-region counts, so holding that fixed is what keeps a re-baseline attributable to a
code change.

⚠ **Every fragment-length parameter is measured, not chosen** — the count-weighted mean and sd of a real
cfRNA library's own pools: **RNA 206.1 ± 98.3**, **gDNA 156.5 ± 124.6**.

⭐⭐ **Those are PRE-capture parameters and the truth files are POST-capture.** Hybrid capture selects for
length, so the simulator draws the length marginal proportional to the capture-weighted opportunity,
`f_post(w) ∝ f_pre(w) · total_eff(w)`. ⛔ **Score against `truth_fragment_lengths.tsv`, never against
`frag_mean`** — the configured parameters describe a library that was never sequenced (TRAPS: capture-selects-for-length).

⛔ **AND THAT IS TWO DIFFERENT pmfs, WHICH BELONGS TO THE DELETED PILOT — the Stage-A panel.** The ladder
on disk hands **both** origins the RNA pool's parameters, because a gDNA-vs-RNA length gap lets the EM
split the origins on **length alone**, bypass calibration and mask its bugs — §0's *why the ladder gives
gDNA and RNA equal fragment lengths*.

### ⭐⭐⭐ THE SIMULATOR'S TRANSCRIPTOME IS A RIGEL INDEX — including the NASCENT ENTITIES (owner, 2026-08-19)

⛔ **`rigel sim` does not read the GTF into transcripts any more. It builds (or is given) a rigel index and
simulates the index's transcript list**, so what is simulated is exactly what `rigel quant` will read:
annotated transcripts PLUS the **synthetic nascent-RNA entities** `index.create_nrna_transcripts` makes —
one single-exon transcript over each multi-exon span, TSS/TES clustered within `NRNA_MERGE_TOLERANCE`
(20 bp), an annotated single-exon transcript adopted where one already covers the span, gene-neutral.
Config key: `index:` (absent ⇒ built into `<outdir>/rigel_index`).

⭐ **Three consequences, and each removes a divergence between the simulator and the tool:**

| | |
|---|---|
| **nascent RNA is a TRANSCRIPT, not a parallel space** | its molecules are `entity.nrna_abundance`, sampled on its own single-exon template, and its reads keep the `nrna_` origin tag so the oracle's origin split is unchanged. ⭐ **How that number is SET is the nascent MODE and nothing else about the machinery depends on it**: under the panel's `mode: sparse` it is drawn per ENTITY — off entirely with probability `1 − on_fraction`, else log-uniform over `abundance_ranges`, ABSOLUTE and independent of the mature level; under `additive_ratio` and `fragment_share` it is `Σ over contributors of abundance × nrna_ratio`, which is why nascent could never exceed mature under either of them |
| **ONE multinomial over every RNA row** | mature rows and entity rows together, `prob ∝ abundance × capture-aware effective length`. ⛔ The mature/nascent FRAGMENT split is no longer imposed as two pools — it follows from molecules and lengths, and is read off the realised origin counts. ⚠ So a nascent-on condition and its nascent-off twin no longer share a bit-identical mature stream: turning nascent on re-allocates the RNA budget, as it physically must |
| **capture binds by GENOMIC OVERLAP** | every probe → genomic blocks → gDNA, and projected onto EVERY transcript whose exons it touches, any gene, any isoform, **either strand** (the library is ds-cDNA at capture). Pieces separated by an intron take `gdna_split_penalty`, as gDNA does. ⭐ Measured after the change: under one probe, gDNA and nascent are enriched **1162× and 1003×** — the same rate, which is the physics (`TRAPS: the-panel-enriches-nascent-by-its-own-probes` records the defect this replaced) |

⭐⭐⭐ **TWO CONSEQUENCES, BOTH MEASURED AND BOTH RULED CORRECT BY THE OWNER (2026-08-19).**

1. ⭐ **CAPTURE DEPLETES NASCENT RNA ~8x, AND THAT IS THE RIGHT BEHAVIOUR** — nascent is **20.27 %** of
   RNA fragments off capture and **2.60 %** under it, re-measured 2026-08-22 at `g50 ss0.99` off that
   condition's own `truth_summary.json` origin counts (1,013,538 and 129,985 of 5 M RNA fragments).
   ⛔ **Neither figure is a knob and the off-capture one is no longer "exactly the knob"**: under
   `mode: sparse` the fragment share is EMERGENT from the per-span draw, and the 20.2 % it lands on was
   PRICED in advance with `expected_rna_weights` rather than requested.
   Probes tile EXONS: a mature molecule is entirely probe-covered, a pre-mRNA is mostly intron, and
   there are ~100x more mature molecules. ⭐ Owner: *"capture should consolidate the evidence of nascent
   RNA to the intron|exon boundaries where fragments can partially overlap probes within exons. The
   introns should be depleted."* ⚠ Investigators who want to profile nascent RNA under capture design
   probes INTO INTRONS; this panel deliberately does not. ⛔⛔ **AND THE REFUSED ALTERNATIVE IS RECORDED
   SO IT IS NOT REBUILT: adding nascent RNA AFTER the capture simulation — i.e. stating the share
   post-capture — is WRONG** (owner, 2026-08-19). Molecules exist first; capture acts on them.
2. ⭐ **THE POST-CAPTURE FL DISTRIBUTION AND ABUNDANCES ARE THE GROUND TRUTH, NOT THE PRE-CAPTURE ONES**
   (owner, 2026-08-19): *"hybrid capture rescales everything, and the post-capture abundances and FL
   distributions are our ground truth baseline."* The truth files already are post-capture
   (`truth_kind: post_capture_empirical`). ⛔ What this ruling KILLED is `simulator_gates.py`'s old
   G-S5, which pooled mature and nascent into one `mu_rna` and asserted the gDNA gap NARROWS under
   capture: that verdict tracked the nascent MIXTURE, and it passed 12/12 only because the old
   simulator imposed 20 % nascent in every condition. ⭐ G-S5 now asserts the LAW — capture lengthens
   **each** RNA population — and cannot be moved by a mixture. ⭐⭐ **That is exactly why it survived the
   sparse rebuild untouched**: re-run 2026-08-22 on the rebuilt ladder it reads **6/6 gates pass** and
   G-S5 itself **12/12 (condition × pool)**; at `g50 ss0.99` mature goes 211.68 → 228.98, nascent
   216.49 → 238.53 and gDNA 216.67 → 240.40. A gate stated on the mixture would have had to be re-derived
   the moment the mixture's weights changed.

⚠ **`transcript_filter` is refused** when the transcriptome comes from an index: the index IS the transcript
set, so filter the GTF before building it.

### ⚠ What the simulator still does NOT do

| | consequence |
|---|---|
| counts are **Poisson by construction** — a multinomial at fixed abundance, measured ω < 5e-5 | nothing dispersion-dependent validates here. Real sj overdispersion is ≤ 0.02–0.03 |
| the PANEL is all **R1-antisense** | ⭐ The ENGINE can emit either since 2026-08-11 (`ReadSimConfig.r1_sense`, gated both ways in `test_strand_sense_convention.py`), but `orchestrator.run_condition_grid` does not expose it — so no ladder condition is R1-sense. ⚠ The gap moved from the simulator to the panel; it did not close. Real cfRNA is dUTP, so this is not urgent |
| the tool's **gDNA reach** assumption is untested | `region_geometry` says gDNA's template is the chromosome, so `taper_g = 1`. True for 50 Mb, false for a 273 bp contig. Latent: it goes live only when a short reference has two regions, and gDNA is no longer simulated on the spike-ins at all |
| ⭐ each **population is written as one contiguous BLOCK** of read names, not interleaved | measured 2026-08-08: a 10 M-fragment condition has **15** origin transitions in BAM order. So any per-fragment truth JOIN that is checked by "does an impossible label appear?" is nearly blind — a one-fragment slip in the `frag_id` key mislabels ~15 fragments in the whole file and need not produce a single impossible one. ⛔ Gate such a join on a COUNT IDENTITY against the scanner's own `stats.total` / `stats.n_read_names` (`_oracle.check_walk_alignment`) and keep the impossible-label check as a secondary that catches a gross slip. `tests/calibration/test_prior_vs_oracle.py` pins both halves, including that the small slip is invisible |

---

## 2. Building it — ⭐⭐⭐ `panel.py`, ONE COMMAND PER STAGE

⛔⛔ **THIS SECTION USED TO STOP HALFWAY, AND THE MISSING HALF WAS THE POINT.** It documented five manual
steps ending at "cache the scans". **Running the tool and scoring it against truth appeared in no recipe
anywhere** and had to be reassembled from `CLAUDE.md`'s 56-row instrument table. Worse, the ORACLE cache
— the origin-split truth every scoring instrument reads — had no step at all: it was a *side effect* of
`pass0_vs_oracle.py --oracle-cache`, so a reader who followed these steps ended up with a panel that
every scorer refused. `scripts/sim/panel.py` is the whole loop (2026-08-11).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
CFG=scripts/sim/configs/gdna_ladder.yaml       # the panel the tool is RANKED on

python scripts/sim/panel.py status   --config $CFG              # ⭐ ALWAYS FIRST
python scripts/sim/panel.py build    --config $CFG              # index + capture probes
python scripts/sim/panel.py simulate --config $CFG --jobs 8     # ~21 min, ~16 GB
python scripts/sim/panel.py cache    --config $CFG --jobs 8     # ⛔ BOTH caches — MANDATORY, see §0
python scripts/sim/panel.py score    --config $CFG --jobs 8 --arms base oracle noop
python scripts/sim/panel.py report   --config $CFG --arms base oracle
```

⭐ **`status` is the command to run when you do not know where you are.** Every stage is expensive and
resumable, and it prints what exists, what is missing, and which stage to run next. It reads real state
— pointed at the pilot and the ladder on the same day it correctly reported the *inverse* cache pattern
(pilot: scan cache 8/8, oracle 0/8; ladder: scan 0/36, oracle 36/36). ⚠ Both panels in that record are
pre-rebuild: the pilot is deleted and the ladder is now 16 conditions.

### ⛔⛔ TWO WORKFLOW GOTCHAS, BOTH MEASURED 2026-08-22 WHILE REGENERATING THE LADDER

Neither is a defect in the recipe above; both are places where the recipe alone leaves you with a panel
that looks rebuilt and is not.

* ⛔⛔ **`--force` DOES NOT REACH THE `simulate` STAGE, SO A REGENERATION SILENTLY DOES NOTHING.**
  `cmd_cache` passes `--force` through to `build_scan_cache.py` and `cmd_build` honours it for the index
  and the probes, but `cmd_simulate` shells out to `simulate_reads.py` with no such flag — and the
  simulator skips a condition whose oracle BAM already exists (`skip_existing`, keyed on the BAM). ⭐ So
  changing a config's `nrna:` block and re-running `panel.py simulate --force` reports success and
  reproduces the OLD reads. **To regenerate, remove the condition outputs first** (delete the per-condition
  directories, or the whole panel directory) and then simulate. ⚠ The existence key is deliberately the
  oracle BAM rather than the FASTQs; that ruling is in `orchestrator.run_condition_grid` and is not the
  thing to change here.
* ⛔⛔ **`panel.py cache` DOES NOT BUILD THE FOUR ZERO-gDNA ORACLE CACHES.** The stage's oracle half is
  `pass0_vs_oracle.py`, which HOLDS OUT every exactly-zero-gDNA row as a false-positive check — so on the
  16-condition ladder it caches 12 and prints `4 zero-gDNA row(s) held out`. The four `g00` rows then have
  no partitions at all, and `calibration_oracle.py` REFUSES to certify a condition whose five partitions
  are absent. ⭐ Fill them one at a time with the same script's per-condition prewarm, which is a pure cache
  fill and scores nothing (⚠ `--_prewarm` is the worker half of `--jobs` and is hidden from `--help`; it is
  the honest way in because it reuses the SHIPPED loader, so a stale cache is still refused):

  ```bash
  for C in $(ls $SUITE/ladder | grep '^gdna_g00_'); do
      python scripts/design/pass0_vs_oracle.py --suite $SUITE/ladder --index $SUITE/rigel_index \
          --oracle-cache $SUITE/ladder/oracle_cache --_prewarm $C
  done
  ```

  ⚠ Do **not** also pass `--conditions <that row>`: the hold-out runs before the prewarm, so a run whose
  only named condition is a zero-gDNA one exits with "no contaminated conditions found" and builds nothing.
  ⛔ The prewarm writes the five partitions and no `_main`, which is why a complete ladder reads
  **oracle cache 12/16** in `status` (§0).

⚠ **`build` cannot carve the reference** — that needs the source genome/GTF the panel was region bound from, which
the panel config does not name. It prints the exact `build_suite_reference.py` command and stops:

```bash
python scripts/sim/build_suite_reference.py \
    --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf \
    --refs chr21 chr22 --ercc -o $SUITE/reference
```

⛔ **THE GATES — run BOTH before quoting anything.** They are not yet stages of `panel.py`:

```bash
python scripts/design/simulator_gates.py --suite $SUITE/ladder --reference $SUITE/reference
python scripts/design/suite_resolves.py $SUITE/rigel_index --suite $SUITE/ladder
```

⚠ **`status` counts a cached oracle condition complete only when all FOUR of the parts it knows about are
present** — `gdna`, `mrna`, `nrna` and the undrained `_main` payload; counting directories would call a
half-written condition done and fail deep inside an instrument later. Gated by
`tests/test_panel_workflow.py`. ⛔ **The parts a SCORING instrument needs are FIVE**, because
`calibration_oracle.py` also requires the per-strand `rna_pos`/`rna_neg` pair and REFUSES rather than
dropping to three arms — and `_main` is the one part a zero-gDNA row never gets (§0). So `status`'s count
and an instrument's requirement are two different questions; ask the instrument.

⚠ **Every panel config states `gdna.genomic_refs: [chr21, chr22]` explicitly.** The engine does **not** infer
which references carry genomic DNA, and a config that asks for gDNA without stating it is rejected — "has
an annotation" is not "is genomic" (TRAPS: annotated-is-not-genomic).

⚠ **The scan cache is refused, not silently accepted, when it does not describe the index it is loaded
against** — keyed on `graph_hash`, `reach_digest` and `payload_schema_digest`. Any accumulator change
invalidates every cache by design, so re-run step 5.

---

## 3. The simulator's own gates — G-S1…G-S6

`scripts/design/simulator_gates.py`, scored on the panel's per-fragment truth (the oracle BAM's read
names). ⭐ **Every one is directional or an absolute count. Not one carries a threshold** — a pass mark on
"how much longer is a captured fragment" would be inventing the capture efficiency curve.

| | gate | form |
|---|---|---|
| **G-S1** | gDNA fragments on an RNA-only reference | absolute count, must be **0** |
| **G-S2** | genomic references carrying gDNA | **≥ 2**, each non-zero, on every gDNA condition |
| **G-S3** | gDNA mean length, capture off → on | **strictly greater** under capture |
| **G-S4** | on-target vs off-target gDNA mean length | on-target **strictly longer**. ⚠ A regression guard, not a falsification — it passed *with* the capture defect present, because the conditional was right and only the marginal was discarded (TRAPS: a-gate-that-already-passed) |
| **G-S5** | mean length of **EACH** RNA population (mature, nascent) separately, capture off → on | **strictly greater** under capture, per condition × pool |
| **G-S6** | gDNA fragments longer than their own reference | **0** |

⛔⛔ **G-S5's FORM MATTERS AND THIS TABLE USED TO STATE THE RETIRED ONE.** It read `|μ_g − μ_r|` narrowing
under capture — a claim about the mature/nascent MIXTURE, which passed only because the old simulator put
20 % nascent in every condition and would have had to be re-derived the moment the nascent model changed
(§1's second consequence). The shipped gate asserts the LAW instead, on each population separately, so it
is invariant to how much nascent there is. ⭐ Re-run 2026-08-22 on the rebuilt sparse-nascent ladder:
**6/6 gates pass**, G-S5 itself **12/12 (condition × pool)**.

⛔ **"On-target" means OVERLAPS A PROBE, not "its start lands in an exon."** The start-territory version is
geometry-confounded and stays inverted under any correct capture model — TRAPS: on-target-by-start-is-geometry. The script prints
the start-territory table underneath as the diagnostic it is.

---

## 4. Can the suite resolve the axis? — `suite_resolves.py`

⭐ **There are no tuned thresholds in it, deliberately.** Every requirement is scored against its
**degenerate value** — the number a structurally blind suite scores — and passes iff it lands strictly on
the non-degenerate side. Those are boundaries, not choices: an unresolvable partition is *exactly* 1.000×,
a suite with no length variation has variance *exactly* 0, a Poisson simulator has ω *exactly* 0, and a
capture arm that is length-neutral narrows the length gap by *exactly* 0.00.

The requirements: (a) a capture **density step**, (b) fragment-length **variance**, (c) non-Poisson
**counts**, (d) termini strictly **inside an exon**, (e) ample **single-stranded** regions, (f) a low-gDNA ×
strong-capture **corner**, (g) **partition resolution**, (h) a **narrowed length-gap** regime.

⚠ **One is known-failing and it is named work**, not a suite defect (the second, `(f)`, was closed by the
2026-08-13 rebuild):
* **(c)** needs the overdispersion mechanism built in *and* replicate conditions.
* ✅ **(f)** PASSES since 2026-08-13. It needs one gDNA RATE in `0 < rate ≤ 0.10`, and the rebuilt ladder
  carries `g05` (rate 0.052632) precisely for it — 2 conditions, against a degenerate value of 0. ⚠ It was
  known-failing on the pilot, and a `{g00, g50, g98}` rebuild would have re-broken it. It opens exactly the
  regime where TRAPS: capture-is-1000x-on-exons says the hardest failure mode lives.

⭐ **The gate's teeth are proven on three degenerate inputs**, each failing for its own reason: a reference
shaped like the deleted 10 Mb suite (121 regions == 121 merged regions), a "starved toy" with both-strand
regions and no single-stranded ones, and truth files written to sd 0 / capture off / no replicates.

---

## 5. How results are evaluated

⛔⛔ **THE 0.8.0 DEVELOPMENT METRIC IS THE CALIBRATION RESULT SCORED AGAINST ORACLE CALIBRATION, NOT THE
END-TO-END TRANSCRIPT NUMBER** (owner, 2026-08-14). The focus of development is calibration; the
transcript table is a thermometer, and a change is judged on what it does to calibration against truth.

⭐ **There are THREE evaluation questions and they need different instruments.**

| question | instrument |
|---|---|
| ⭐⭐ how wrong is **CALIBRATION**, against **ORACLE calibration**? — **the 0.8.0 metric** | `prior_vs_oracle.py` (calibration's endpoint, the `LocusPriors` the EM actually reads) · `solvability_audit.py` (pass-0) · `pass0_vs_oracle.py` for the T/C/P decomposition |
| how WRONG is the end-to-end answer? — per transcript, per pool, against per-fragment truth | ⭐⭐ `panel.py score` / `report`, i.e. `quant_accuracy.py`. **The only end-to-end scorer**, and a thermometer rather than the target |
| WHERE did the misassigned fragments go? | `rigel.sim.net_flow` (`analyze_net_flow`, `FlowData`) |

⛔⛔ **AND A CEILING ONLY PRICES WHAT ITS ARM CAN REACH — measured 2026-08-13/14, this invalidated a
whole family of them.** `effective_lengths_em` is built inside `_setup_geometry_and_estimator`, called at
`pipeline.py:816` — **before** `assemble_priors` at `pipeline.py:839` — and every existing measurement arm
patches `assemble_priors`, so the
effective-length shrinkage has **never been priced by any ceiling**, and every ceiling number to date was
measured with a wrong ruler installed. ⭐ Say which call your arm patches, and check it sits downstream of
everything you mean to price.

⛔ **`rigel.sim.analysis` and `scripts/sim/evaluate_suite.py` were DELETED on 2026-08-11** (owner). That
was a 1,589-line second scorer running the tool and rendering its own accuracy tables beside
`quant_accuracy.py`'s, against its own definition of truth — and two scorers is how a baseline and a
ceiling drift apart (`TRAPS: score-the-consumers-own-count`). The 618-line flow decomposition, which has
no duplicate, moved to `net_flow.py` with its tests; the other 971 lines went.

**Hard per-fragment label recovery is the wrong target.** An unspliced RNA fragment and a gDNA fragment
from the same locus can be sequence-identical and genuinely unrecoverable. Build the per-locus
`flow[true][assigned]` matrix and reduce to `net(a→b) = flow[a][b] − flow[b][a]`: symmetric,
unrecoverable misassignment cancels, and **only systematic bias survives**.

⚠ **Report absolute per-transcript error alongside the net**, because net cancels.
⚠ And hard-label metrics are nearly blind to a calibration-prior change — TRAPS: hard-labels-miss-soft-change.

⛔ **Missing:** the **soft** 3-pool surplus. The surviving code computes only the hard-label version, and
the soft one is the metric that actually sees a prior change.

---

## 6. The test suite

```bash
python -m pytest tests/ -q                     # ⛔ never bare `pytest` — the repo root must be on sys.path
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
```

⭐ **THE STANDING BASELINE IS GREEN: 0 failures, 3,497 passing, 0 skipped, 11 xfail** (re-derived here
2026-08-17 with `python -m pytest tests/ -q`; it read 3,404 / 10 xfail on 2026-08-14, 3,397 on
2026-08-13 and 3,235 / 7 xfail on 2026-08-11). ⛔ **`CLAUDE.md` is the HOME of this number and carries
every delta ACCOUNTED rather than adjusted; this line is a convenience copy and had gone stale by
+92 / +1 xfail before it was re-derived** — when the two disagree, re-derive rather than pick one.
⚠ It moved **3,496 → 3,497 inside that same session**, and the `+1` is ACCOUNTED rather than adjusted:
one new **dev doc**, which `test_no_jargon_labels` parametrises over and the docs-boundary gate does not
(`pytest --collect-only -q | grep <stem>` prints exactly the one case). A doc-only edit does not move the
count; adding or retiring a FILE does.
**Any failure is a regression** — which is a stronger
and cheaper rule than counting expected ones. ⚠ Re-derive it rather than adjusting it; several tests are
parametrised over doc, source and script files, so adding or retiring one moves the count by a few.

⚠ **The predecessor of this paragraph said "21 `test_golden_output` failures plus one", and both halves
have since resolved**: the goldens were regenerated, and
`tests/scenarios_aligned/test_paralogs.py::test_gdna_sweep[gdna_100]` — a real EM unidentifiability, not
flakiness (TRAPS: identical-paralogs-are-bimodal) — now passes. ⛔ If that row starts failing again, do
not fix it by moving a seed.

⛔⛔ **AND THE OLD WARNING ABOUT THE GOLDENS WAS RIGHT ABOUT THE MECHANISM, WHICH IS NOW MEASURED.** It
said the goldens "run under the default sampling mode, so a flaky expectation baked in now is permanent".
The default is `EMConfig.seed = None` with `assignment_mode = "sample"`, and two runs of the identical
pipeline on the identical BAM genuinely return different transcript counts — 4 transcripts differing by
up to 43 fragments on one small toy. `TRAPS: the-deliverable-is-not-reproducible-by-default`.
⭐ So: **regenerate the goldens twice and diff**, and pin `EMConfig.seed` in any instrument that compares
two end-to-end runs (`quant_accuracy.py` does, and prints a `base_reseed` noise floor beside the effect).

---

## 7. Development discipline for test substrates

**Develop on controlled toys, validate on real data.** A big suite has confounds that hide mechanisms; a
toy ranks hotspots backwards (TRAPS: toys-rank-hotspots-backwards). Both, in that order.

### ⛔⛔ PROFILING TARGETS A HIGH-DEPTH REAL RNA-seq LIBRARY, **NOT** cfRNA — owner, 2026-08-17

⭐⭐ **This REVERSES a standing instruction**, so read it as a correction rather than an addition. The old
rule — *"profile on real cfRNA, never a small synthetic suite"* — got the second half right and the first
half wrong, and it may still be written that way at older sites. ⚠ Where two disagree, **this is the ruling**;
`ISSUES: performance-memory-bounded-solve` carries it forward.

⛔ **Why cfRNA is the wrong substrate for a COMPUTE measurement, even though it is a fine one for an
ACCURACY measurement.** The cfRNA libraries on disk are **sparse and small**: most confident-gDNA regions
carry zero counts (64–94 % across libraries, `DESIGN.md`), so they under-fill exactly the structures whose
cost the performance work is about — the grid solve, the per-slot arrays, the fragment buffer. A library
that never fills them cannot price them, and a profile taken there ranks hotspots by a distribution the
production target does not have. ⭐ This is `TRAPS: toys-rank-hotspots-backwards` applied one rung up: the
lesson was never "cfRNA", it was **"profile on the substrate whose SHAPE you are optimising for"**, and a
small sparse library is a toy in the only sense that matters here.

⭐ So the profiling substrate is, in order of preference: a **deep, real, high-complexity RNA-seq BAM**;
never a panel condition (`TRAPS: toys-rank-hotspots-backwards` outright); never cfRNA. ⚠ The two rules
this does **not** touch, because they are about design inputs rather than profiling substrates:
`TRAPS: real-data-is-a-test-input` (real data is a TEST input, never a DESIGN input — still binding, and
it is *why* the profiling target is a domain call rather than a tuning input) and the accuracy panel,
which stays the 16-condition ladder.

⭐ The instruments are `scripts/profiling/profiler.py` (whole pipeline; wall clock **and** per-phase peak
RSS across `scan` / `calibrate` / `quant`) and `scripts/profiling/scan_profile.py` (the scan alone, swept
across thread and chunk budgets). `scripts/README.md` is their index and both are gated by
`tests/test_scripts_index.py`. ⛔ Set `OMP_NUM_THREADS` deliberately: the shipped default is all cores and
a threaded number answers a different question from a single-thread one.

**Any both-strand stress test needs ample single-stranded regions** — the population prior trains on them,
and a "starved toy" is one of the three degenerate inputs `suite_resolves.py` is proven against.

**How to A/B honestly:** in-process, opposite extremes, never on a saturated condition, one thing varied,
and both arms sharing their random input (TRAPS: perturb-every-gate) — a byte-identical hard-label result is **no
evidence**.
