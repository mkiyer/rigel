# SESSION HANDOFF — ⭐⭐⭐ CERTIFIED RNA IS PER-(EDGE, SIDE) EVIDENCE, AND ψ READS NONE OF IT

> ## ⛔⛔ §0 READ THIS FIRST — §2f WAS SETTLED ON 2026-08-05 AND THE ANSWER IS A **NEGATIVE**
>
> **§2d's raw-count λ term is refuted and MUST NOT be built.** Everything below §2d that assumes it —
> §3's plan, §4's gates T1–T7, the "five sites" — is superseded by **§8**, which is the only current
> section of this file besides §1, §2a–§2c, §5 and §7.
>
> The one-line reason: the term `S·log(1 − f_g)` is honest only if the term §2d drops is small, i.e.
> only if the splice-visibility `q → 1`. Measured against origin-split truth on all 36 ladder
> conditions, `q`'s mass-weighted median is **0.19–0.71** and the term is **worse than the
> uninformative reference on 12 of 36 conditions**, worst **+0.4578** mwae.
> `scripts/design/certified_q_census.py` · `tests/calibration/test_certified_rna_licence.py` · §8.
>
> ⭐ `CERTIFIED_RNA_REPORT.md` beside this file is the same finding written for an **external reviewer** —
> minimal jargon, the derivation from scratch, and §7's next steps in full. Act from §8; hand over the
> report. Both are ephemeral and both are deleted together.

    Written 2026-08-05. ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** The six permanent docs are
      `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **DELETE THIS FILE when its task lands**, promoting §2 into `EQUATIONS.md`, §7's lessons into
      `TRAPS.md`, the measured numbers into `ROADMAP.md`, and the scope ruling into `DESIGN.md`.

---

## §1 THE TASK, AND THE STATE OF THE TREE

**A spliced fragment cannot be gDNA.** So every spliced count is *certified RNA* — origin-labelled by the
deposit rule, needing no length model, no strand model and no deconvolution. It is the only observation in
the tool that is already solved.

⛔ **ψ reads none of it.** `simplex_logodds`'s own docstring: *"There is NO spliced term: `mass_spliced` is
consumed only by the returned `rna_mass`, never by ψ."* The counts reach the final mass sum and one
special-cased slot of `bp_solver`'s MEASUREMENT stream, and they never become an object's **own** evidence.
Measured consequence on the toy: an EDGE holding **11,026 certified-RNA fragments** reports `f_g = 0.5098`
— the uninformative reference — and over-calls gDNA **17×** (TRUE 11 fragments, PRED 186).

**The task: make certified RNA into per-(EDGE, side) composition evidence, derived, and land it.**

### State of the tree — branch `fragment-length-gold-standard`, on `87a5468e`

Nothing is committed. `ruff check src/ tests/ scripts/` clean.
**Suite: 22 failed / 2240 passed / 3 xfailed** — the 22 are the standing set (21 `test_golden_output` + the
paralog row). A 23rd failure or any other name is a regression.

| file | state |
|---|---|
| `src/rigel/second_pass.py` | ✅ **`lift_choices` LANDED** (§5) |
| `src/rigel/pipeline.py` · `tests/calibration/_oracle.py` · `toy_harness.py` | ✅ **§5's WIRING LANDED 2026-08-05** — `_drain_side_buffer(_lift=)` publishes `(undrained, choices, node_types, junctions)`; `OracleTruth.from_bam(drain_with=)` replays them per partition; `run_toy` drains the whole and prints the ambiguity count. On `tes_readthrough`: 36 held, 36 re-deposited, `edge_spliced` **224 → 269 (+20 %)**, ambiguous **0**. ⭐ Both perturbations fire — undrained partitions trip sum-to-full, and the DRAINED whole in `drain_with` trips `lift_choices`' key pool |
| `scripts/design/certified_q_census.py` | ⛔ **NEW — the instrument that refuted §2d.** Verdict in its docstring; §8 |
| `tests/calibration/test_certified_rna_licence.py` | ⛔ **NEW — 28 gates, 6 perturbations, all firing.** The licence: what a certified count may and may not claim |
| `tests/test_lift_choices.py` | ✅ 8 gates, 6 perturbations, all firing |
| `scripts/design/toy_harness.py` | ✅ new rung **`tes_readthrough`** (§6) |
| `scripts/design/certified_rna_audit.py` | ✅ the (a)/(b)/(c) audit on a TA × TB grid; verdict in its docstring |
| `scripts/design/vertex_ceiling.py` | ✅ the re-solve ceiling harness — reusable override plumbing at `build_node_init` |
| `tests/calibration/test_vertex_reference.py` | ✅ 7 gates on ψ's reference. ⭐ **G6 asserts ψ is blind to the spliced channel — it MUST FAIL when this lands, and must be re-pointed, not widened** |

⚠ The oracle cache at `~/Downloads/rigel_runs/suite/ladder/oracle_cache` is **VALID — do not delete.**

---

## §2 ⭐⭐⭐ THE DERIVATION — the full and correct design

### 2a. The three banks, and what each certifies

At one line the accumulator stores three populations with **different component sets**
(`EQUATIONS.md` §3.6):

| bank | the molecules | components |
|---|---|---|
| `unspliced_count` | crossed the line **contiguously**, spliced nowhere in the fragment | ⛔ gDNA + RNA — **this is the deconvolution problem** |
| `spliced_count` (`edge_spliced`) | crossed the line **contiguously**, having spliced *elsewhere* in the fragment | ✅ certified RNA |
| `junction_count` (`sj_count`) | **jumped** from this line | ✅ certified RNA |

### 2b. ⭐⭐⭐ THE RULE — certified RNA lands on the flank where the molecule's BODY lies

An intron spans `[lo, hi)`. A mature molecule that uses it has its body in the exon **outside** the intron
at each end — below `lo`, above `hi`. So:

        junction_count_lo   this EDGE is an intron's genomic-LOW  end  ⇒  the RNA is on its LOW  flank
        junction_count_hi   this EDGE is an intron's genomic-HIGH end  ⇒  the RNA is on its HIGH flank
        spliced_count       crossed contiguously                      ⇒  the RNA is on BOTH flanks

⭐⭐ **This is STRAND-INDEPENDENT in genomic terms, and that is exactly why the arrays are named `_lo`/`_hi`
rather than donor/acceptor.** `EQUATIONS.md` §3.6c: *"the index's `FLAG_DONOR_s` bit marks the genomic-LOW
end of an `s`-strand intron on BOTH strands, so on `−` it sits at the transcript's biological ACCEPTOR."*
The biological names (donor/acceptor, splice-in/splice-out) flip with strand; the **flank does not**.

⭐ **Worked example** (owner, 2026-08-05). `TA+ (1,000, 2,000) (9,000, 10,000)`, a spliced fragment with
blocks `(1,900, 2,000) (9,000, 9,100)`. At the line at 9,000 — the intron's genomic-HIGH end, biologically a
`+`-strand acceptor — the molecule arrives and **continues upward** into node `[9,000, 10,000)`. It is
present on the HIGH flank and absent on the LOW (intron) flank. So the HIGH flank carries certified RNA and
must be solved as such; the LOW flank must not, and its own population is genuinely a mixture (an intron can
hold unspliced RNA, which is RNA and must be deconvolved).

### 2c. ⛔ THE ERROR THIS CORRECTS, so it is not re-made

An earlier version of this file argued that `junction_count` *"cannot constrain this EDGE"*, citing
`intron|exon @2,000` on `alt_splice`: junction flux **14,492** with the EDGE's true `f_g` = **1.0000**.

⭐ **Both facts are true and they are about different things.** That `f_g` is ONE number for the whole line —
the composition of its *unspliced crossing* population, which is pure gDNA because mature RNA cannot cross
an exon↔intron seam contiguously (`TRAPS.md` F9). The LOW flank of the same EDGE *does* carry certified RNA,
which is precisely why `rho_lo` = 47.72 while `rho_hi` = 0.062 there. ⛔ **The mistake was conflating "the
EDGE's composition" with "the (EDGE, side) composition."** The reframe already knows an EDGE presents two
totals, one per flank; the composition evidence must be two-sided in the same way.

⭐⭐ **And the consequence is that the channel is far larger than first thought.** It is not 5,991 terminus
EDGEs — it is **every junction in the annotation**, on its exon-facing flank, plus every contiguous crossing
by a molecule that spliced elsewhere.

### 2d. ⭐⭐⭐ THE λ TERM — one coefficient, and it is the term the reference stands in for

At a slot the unspliced total is `M` with `f_g = σ(λ)`, so the RNA density the bank implies is
`ρ_R(λ) = (1 − f_g)·M/E_r`. A certified count `S` is Poisson on that RNA, so its λ-dependence is

        log L(λ)  =  S · log(1 − σ(λ))   +   (a term scaled by q)

where `q` is the chance a crossing RNA fragment shows a *visible* splice. ⭐ **The first term carries no `q`
at all.** Keeping only it makes the claim **one-sided** — it can push `f_g` down, never up — which is the
honest reading of a floor and a shape this project has already licensed (`EQUATIONS.md` §7.2 uses *"a flat
one-sided prior on the RNA excess, truncated at the observed total"*).

Now compare ψ's own RNA arm, which is `_JEFFREYS_REF · log(1 − f_g)` = `½ · log(1 − f_g)`:

        no evidence:   ψ_RNA  =  (½)     · log(1 − f_g)      ← the Jeffreys reference
        S certified:   ψ_RNA  =  (½ + S) · log(1 − f_g)      ← the SAME term, now measured

⛔ **This is not a new channel bolted on. It is filling in the term the reference was standing in for** — the
reference IS the `S = 0` limit. One coefficient goes from a constant to a per-slot array.

### 2e. ⭐⭐ THE θ TERM — certified RNA is STRANDED, and that is a second axis

Certified RNA carries a strand, and ψ has exactly the axis for it: λ is the two-group split (gDNA vs
RNA-total) and **θ = arcsin(τ) is the tilt between RNA+ and RNA−**. So one certified count makes **two**
claims and they belong on two axes:

| claim | axis | term |
|---|---|---|
| "there is at least this much RNA here" | λ | §2d — coefficient `½ + S` on `log(1 − f_g)` |
| "and it is on THIS strand" | θ | a Gaussian/Beta term on `θ`, strength from `S_pos` vs `S_neg` |

⚠ **The two keyings differ and must not be crossed:** `spliced_count` is keyed by **genome** strand,
`junction_count` by **transcript** strand (`node_geometry`). One conversion, in one place, gated.
⚠ **Build λ first and measure it alone** (`TRAPS.md` B21 — pricing two halves together hid a 14× split
between them). θ is a second, separately-priced step.

### 2f. ⛔ THE TWO OPEN DERIVATION QUESTIONS — settle these before writing the term

1. **Does `S` enter as a raw count, or as `S` scaled by an opportunity ratio?** §2d's coefficient is `S`
   because the term is a Poisson log-likelihood in the same frame as the unspliced bank —
   `node_geometry` gives `spliced_count` the divisor `eff_rna` on the ground that *"a contiguous crossing is
   a contiguous crossing whatever the molecule did elsewhere"*, which IS the same-frame claim. ⛔ **That
   argument is stated for `spliced_count` and NOT for `junction_count`, whose divisor is `eff_junction`.**
   So the junction half needs its own frame statement, and getting it wrong is `TRAPS.md` C0's trap (two
   divisors from one pmf can disagree with opposite sign). ⭐ Gate whichever answer against brute-force
   enumeration — the project's standing rule for any divisor.
2. **Where does the flank-side evidence attach?** `NodeInit` is per-SLOT, one number per object, while the
   evidence is per-(slot, side). Two candidates: (i) attach it to the flank TOTAL the reframe already
   computes, alongside `rho_lo`/`rho_hi`; (ii) emit it as a directed message on the hop that leaves that
   flank. ⭐ (i) keeps one home for "what this slot presents to that neighbour" and is likely simpler;
   (ii) is closer to how the graft works today. Decide by reading `node_total_density` — do not guess.

---

## §3 THE IMPLEMENTATION PLAN

### 3a. Order of work — ⛔ not optional

| | step | why this order |
|---|---|---|
| 1 | **Wire the drain into the toy** (§5) | the toy currently understates the very channel being built, by a measured 13.5 % |
| 2 | **Settle §2f's two questions**, on paper, then as a brute-force gate | a wrong divisor is silent |
| 3 | **Write the gates** (§4), verified failing | `TRAPS.md` A2 |
| 4 | **Implement the λ term** (§3b) | one thing at a time (B21) |
| 5 | **Measure on the toy**: `tes_readthrough` + `certified_rna_audit.py` + ⛔ `zero_controls.py` **both arms** | owner's standing requirement |
| 6 | ⛔⛔ **The LADDER arm, BEFORE `src/`** | `TRAPS.md` B18 — three toy-positive changes have been panel-negative |
| 7 | land in `src/`, re-point G6, then price the **θ** term as a separate step | |

### 3b. The sites, and the observable for each (`TRAPS.md` B14 — a count of gates says nothing)

| | site | change | the observable that moves |
|---|---|---|---|
| 1 | `node_geometry` | **none** — `spliced_count`, `junction_count_lo/hi`, `eff_rna`, `eff_junction_lo/hi` all exist | — |
| 2 | new `certified_rna_lambda_factor(...)` | build the `(m, K)` factor `S·log(1 − σ(λ))`, per flank | the array is nonzero EXACTLY where a certified count is |
| 3 | `node_init.build_node_init` | **one line**: `tau_lam += density_factor_precision(cert, lam_grid)`, ungated, exactly like the existing `tau_len` | `tau_lam` rises where `S > 0`; `prec_pos`/`prec_neg` go positive |
| 4 | `simplex_logodds._rna_arm` | take a per-slot coefficient, mirroring `_gdna_arm`'s `global_logprior` | ψ's `log(1−f_g)` coefficient; ⭐ **G6 flips to failing** |
| 5 | `calibrate.py` | build + thread the factor, mirroring `_build_length_loglik` | the kwarg arrives at `node_sweep` |

### 3c. ⭐⭐ WHAT IT DELETES — the simplification, which is the point

* **ψ's two arms become symmetric.** `_gdna_arm` already accepts a per-slot prior added to `C·log f_g`;
  `_rna_arm` accepts nothing and is `(1,K)`. After this, both are `(C + evidence_c)·log f_c` — same shape,
  same signature. ⭐ The asymmetry a whole investigation chased (§7) is closed by the same change.
* **Zero new precision arithmetic.** `density_factor_precision(factor, lam_grid)` is already generic — its
  docstring says the maths is generic and it is *"named for its first caller"* — and **a FLAT row returns
  exactly 0**. So the channel **self-gates**: `S = 0` ⇒ flat ⇒ no contribution. **No branch, no licence, no
  threshold**, which matters because a threshold here has been refused three times (`TRAPS.md` B11, D4f,
  D4g).
* **`bp_solver`'s MEASUREMENT stream** currently special-cases *"the anchor gDNA count + the spliced RNA
  count"*. Once certified RNA is own-evidence at its object, that path is a **deletion candidate** —
  converge and delete (`TRAPS.md` G3). ⚠ Measure before removing.
* It reaches NODEs for free: the banks are EDGE-only and the relay already carries EDGE evidence to flanking
  NODEs. No new plumbing.

---

## §4 THE GATES

| | gate |
|---|---|
| T1 | `S = 0` everywhere ⇒ output **byte-identical** to HEAD. ⭐ the noop falsification |
| T2 | `density_factor_precision` on the certified factor equals **`S`** to floating point, gated against a brute-force Fisher computation |
| T3 | ⭐ **the FLANK rule**: a junction's genomic-LOW-end EDGE gets the evidence on its LOW flank and **nothing** on its HIGH flank; the mirror for a HIGH-end EDGE; `edge_spliced` gets **both**. Assert per flank, not per object |
| T4 | ⭐ **strand-independence in genomic terms**: the same geometry on a `−` transcript puts the evidence on the same flank. `nested_exons_neg` / `splice_both_strands` are the rungs |
| T5 | the term is **ONE-SIDED** — it can only lower `f_g`, never raise it |
| T6 | on `tes_readthrough`, @9,050's `prec_pos`/`prec_neg` > 0 **and** its unspliced gDNA over-call falls from 17× toward 1× |
| T7 | `_rna_arm` at coefficient `½ + 0` is byte-identical to today's `_JEFFREYS_REF · log1m` — the refactor is a no-op without evidence |

**Perturbations, each naming the gate that must fire:** force `S ≡ 0` → T1 green, T6 fails · delete the
`density_factor_precision` line → T6 fires, G6 stays green · put `S` on the `log f_g` arm → T5 fires ·
**swap `_lo`/`_hi`** → T3 fires · **use donor/acceptor instead of genomic lo/hi** → T4 fires on the `−`
rung · use only `S_pos` → the strand-symmetry gate fires.
⛔ A perturbation that fires nothing is a hole in the gate set. ⚠ And check it had opportunities at all
(`TRAPS.md` A14) — two of the four `verify_toy_substrate` perturbations are `−`-strand-only and are inert on
a `+`-only spec.

---

## §5 ✅ LANDED: `rigel.second_pass.lift_choices` — and the ~30 min left to wire it

**The problem.** Splitting a BAM by origin and re-scanning reconstructs the pass-one payload exactly, which
is what makes an origin-split oracle a valid truth source. The second pass breaks it: its multinomial is
scored against the *whole* payload's densities, so `Σ partitions ≠ whole` (`TRAPS.md` B9). Every calibration
number in this campaign is therefore an **undrained** tally, which is why those instruments say so.

**The repair.** B9 constrains **where** the choice is made, not whether. Score and draw ONCE on the whole,
then replay each fragment's chosen hypothesis inside whichever partition holds it:

```
scores  = score_held_fragments(whole, ...)
choices = choose_hypotheses(scores, whole, seed=...)          # <- ONCE, on the whole
drained = drain(whole, choices, ...)                          # the tally the solver reads
lifted, ambiguous = lift_choices(whole, partitions, choices)  # the truth, per origin
for p, ch in zip(partitions, lifted): drain(p, ch, ...)
```

The key is `_DEFERRED_RECORD_FIELDS`, imported not restated, resting on `DeferredFragments`' own guarantee
that *"two records that tie on that key are identical records."*

⭐ **It takes ALL partitions in one call, and that is load-bearing.** The per-key queue's state must be
shared across partitions or two of them take the same entry and leave another unused. ⛔ The first version
took a single partition with a per-call queue and had exactly that bug — perturbation **P2** found it, the
gate set could not.

### ✅ THE WIRING LANDED 2026-08-05 — and one A14 fact it turned up

`_drain_side_buffer` gained an out-parameter `_lift` publishing `(undrained, choices, node_types,
junctions)` — the four things `lift_choices` needs, all of which exist only inside that function, so a
caller cannot re-derive them and drift (`TRAPS.md` A11). `OracleTruth.from_bam(drain_with=)` consumes it,
lifts ALL partitions in one call and drains each. ⭐ `from_parts`' existing sum-to-full then IS the drain's
end-to-end identity gate, for free. `run_toy` drains the whole, hands the drained payload to `calibrate`,
and prints held / deposited / chose-spliced / **ambiguity count**.

**Measured on `tes_readthrough`:** 36 held, 36 re-deposited, all 36 chose a spliced path,
`edge_spliced` **224 → 269 (+20 %)**, ambiguous **0** (distinct spans ⇒ the lift is exact).
**Both perturbations fire:** drained-whole with UNDRAINED partitions trips sum-to-full
(`node_spanning_count`, max|diff| = 5); the DRAINED whole passed as `drain_with[0]` trips `lift_choices`
with "36 choices for 0 held records" — the loud failure the contract note predicted.

⭐⭐ **AND THE A14 CHECK CHANGED WHAT THE FIX MEANS.** On `spliced_exons`, `TA_single_exon` and
`alt_splice` the side buffer holds **0** fragments — deferral needs two intron paths through one mate
gap, and only `tes_readthrough` has two junctions off one shared donor. So the drain is **inert on every
earlier rung**: the zero controls and every prior toy number are unchanged *by construction*, not by
luck, and the undrained campaign was undrained only on the one rung that could hold anything.

### ⛔ THE ORIGINAL NOTE, kept because the contract subtlety is the load-bearing part

`OracleTruth.from_bam` (`tests/calibration/_oracle.py`:111) does `_split_bam` → scan each → `from_parts`,
and **`from_parts` already validates sum-to-full** — so wiring the drain turns the existing validator into
the drain's own end-to-end identity gate for free.

Add an optional `drain_with=` and drain each partition with its lifted choices before `from_parts`.
⛔ **The subtlety: the drained payload's deferred bank is EMPTY by design** (*"after it nothing is held"*), so
it cannot supply the key pool. `drain_with` must carry the **UNDRAINED** whole alongside the choices —
`(undrained_whole, choices, node_types, junctions)` — while `full_payload` receives the **DRAINED** whole.
Backwards yields an empty pool and every partition raises, which is at least loud.

Then `toy_harness.run_toy` (~line 447): drain the whole, hand the drained payload to `calibrate`, and both
payloads plus the choices to `OracleTruth.from_bam`. ⚠ **Print the ambiguity count** — it bounds the truth
error, and dropping it is the defect gate L6 exists to expose.
⚠ `pass0_vs_oracle.measure_condition` wants the same treatment, but that is a separate change and the
panel's oracle **cache** must be invalidated with it.

---

## §6 THE MEASURED FACTS THIS RESTS ON

**The rung — `tes_readthrough`**, in `toy_harness.SPECS`. `TA+ (1,050, 2,000) (9,000, 9,100)` and
`TB+ (1,000, 2,000) (9,050, 11,000)`; two junctions from one shared donor at 2,000.

| structure | why |
|---|---|
| EDGE @9,100 | TA's TES, **no junction on it** — yet TB's spliced fragments cross it contiguously. ⭐ The population the new TSS/TES EDGEs create |
| EDGE @9,050 | TB's junction acceptor **and** a plain contiguity line for TA, whose exon 2 spans 9,000–9,100 unbroken. ⭐ One line, both banks at once |
| NODE [9,000, 9,050) | TA exon AND TB intron, same strand, 50 bp — below one fragment length, so no resolvable density (`TRAPS.md` D8) |

⛔ **No earlier toy can make this fragment.** On every other rung the exons ARE the nodes, so a spliced
molecule never crosses an interior line contiguously and `edge_spliced` is **structurally zero everywhere** —
measured 0 on `alt_splice`, including at exons holding 68,000 RNA fragments. **So this channel could never
have been found on the substrate the project reaches for first.**

**Substrate: verified.** Every accumulator bank byte-identical to the truth-derived deposit, junction axis
exactly the spec's introns. ⚠ 0.65 % of fragments are **deferred, correctly** — TA's and TB's junctions
share the donor at 2,000 and differ only in acceptor, so an ambiguous unsequenced mate gap implies two
different `L` values. That is the owner's own ruling in `tests/native/_accumulator_reference.py`, and it
understates the certified channel at @9,100 by 13.5 % until §5's wiring lands.

**The audit — 24 of 24 grid cells** (`certified_rna_audit.py`, TA × TB ∈ {0, 30, 300, 3000}):

* ✅ the bank is populated and has a divisor (`eff_rna` = 202.8), is held OUT of the deconvolution, and is
  added back as RNA: `PRED rna` tracks `TRUE rna` to <2 %;
* ⛔ **no precision is emitted** — `cm_p = cm_n = 0` in every cell — so it informs nothing about the
  unspliced split at its own line, and the gDNA is over-called there:

| cell | object | unspl n | spliced | TRUE gDNA | PRED gDNA |
|---|---|---|---|---|---|
| TA 3000 / TB 30 | @9,050 | 365 | 11,026 | **11** | **186** — 17× |
| TA 300 / TB 30 | @9,050 | 288 | 8,739 | 9 | 26 — 2.9× |
| TA 3000 / TB 30 | @9,100 | 116 | 364 | 12 | 12 — exact |

⭐ **@9,050 is the diagnostic and @9,100 its control.** @9,100's flank includes the 1,900 bp exon
[9,100, 11,000), which the relay can speak from — exact. @9,050's flanks are two **50 bp** nodes, so the
relay has nothing to offer and the answer never leaves ½, whatever the object's own certified count. ⚠ And
it moves with **TB**, not with its own count: at TB = 0 it holds 11,649 certified fragments and reads
0.0725; at TB = 30 it holds 11,026 and reads 0.5098. More own evidence, worse answer.

⛔ **Score the MASS, not `f_g`.** The solver's `f_g` is the gDNA fraction of the **unspliced** population;
the oracle's per-object fraction is spliced-INCLUSIVE. Comparing them directly reads a 1.5 %-of-mass defect
as a half-unit error — which is how the first run of the audit mis-reported it.

**Panel scale** (`g50 ss0.50`, substrate census only): 5,991 of 35,038 live EDGEs carry `edge_spliced > 0`,
totalling **85.3 %** of EDGE unspliced mass off capture (24.2 % under capture). ⭐ And per §2b the real
population is larger again, because every junction's exon-facing flank qualifies.

---

## §7 ⛔ WHAT IS DELIBERATELY NOT NEXT — one line each, so it is not rebuilt

| | verdict |
|---|---|
| **the simplex-vertex problem** | ⛔ **NOT A BUILD — a measured negative.** The evidence-free objects are HONEST: `\|Δ\|/sd(f_g)` has median **z = 0.5–0.6** on both vertices and every geometry, so their error is inside their own 1σ. And no prior-free solve can do better: every proper prior on [0,1] has a median strictly inside (0,1), and an object with zero composition information has posterior = prior. The 24.4 %-of-deliverable ceiling `vertex_ceiling.py` measures is **the value of missing information, not headroom.** Its 7 gates stay green and cost nothing |
| **log space vs linear space** | ⛔ settled: **it does not matter.** A posterior median is transform-invariant — re-solving on a uniform-`f_g` grid converges to the log-odds answer at every message precision. And log-odds is *numerically far better*: Beta(½,½) is singular at the vertices in linear space, so linear needs ~10⁷ grid points to match what λ does with 256 |
| **fitting an RNA population prior (`logP_r`)** | ⛔ circular at pass-0, whose entire purpose is to produce the substrate a prior is fitted ON |
| **`length_likelihood`** | ⛔ **stays OFF**, and adding unequal fragment lengths to the ladder stays undone. The ladder's equal lengths and κ = ½ *deliberately* strip every per-object channel so message propagation is the only thing under test. ⭐ The certified channel is different in kind — origin-labelled, no length model, no strand model, no deconvolution — so it adds evidence without conflating error sources |
| **the confidently-wrong population** | ⚠ **2.8 %** of pass-0 node error, and **0.0 %** in the stratum carrying 65 % of it. Real but not first: 92.3 % of the error is on objects with no evidence channel at all |
| **the capture level residual** | ⚠ **the panel's #1 error and a separate build.** Under capture an exon's own total density equals its true gDNA density to **0.2 %** (923× vs 921× background), while the level the relay delivers is **270×** — because it is measured at the flanking EDGE, whose fragments straddle the probe boundary (281×/292×). Toggling capture off on the same condition drops Σ\|err\| **14.4×**. `EQUATIONS.md` §3.5c names this residual; the dissection reranks it from a footnote to first. ⭐ The certified channel is plausibly its discriminator too — `edge_spliced ≈ 0` beside a 921×-enriched exon is the evidence that the enrichment is gDNA — but that is a second measurement, after this one |

### Lessons earned here, for `TRAPS.md` when this lands

1. ⭐⭐ **A channel that is structurally absent from every toy will never be found by the instrument the
   project reaches for FIRST.** `edge_spliced` is identically zero on every rung where the exons are the
   nodes. A7/A8's shape, one level up: not "can the benchmark resolve the axis" but "can the substrate
   *contain* the population".
2. ⭐⭐ **"This bank cannot constrain this object" is a claim about the OBJECT; the evidence may still
   constrain one of its FLANKS.** §2c — the same EDGE was simultaneously pure-gDNA in its unspliced
   crossings and RNA-bearing on its low flank, and reading the first as settling the second cost the whole
   junction population.
3. ⭐ **A fixture that cannot express the thing being gated makes the gate vacuous.** Three of six
   `lift_choices` perturbations initially fired nothing: the record key was untested because every fixture
   record already differed in earlier fields, and the queue and over-consumption guard were never reached.

---

## §8 ⭐⭐⭐ THE SETTLEMENT OF §2f (2026-08-05) — WHAT IS SOUND, WHAT IS REFUTED, WHAT IS BLOCKED

### 8a. Q1 — "raw count, or scaled by an opportunity ratio?" — **RAW COUNT, and the question was the
wrong one**

Write the two banks at one line as ONE population of contiguous crossings split by whether a splice was
VISIBLE, with `q` that visibility probability:

        E[S | f_g, q]  =  (q / (1 − q)) · (1 − f_g) · M

        log P(S)  =  S·log(1 − f_g)   −   (q/(1−q))·(1 − f_g)·M   +   const
                     \___ retained ___/    \_______ dropped _______/

⭐ **The retained coefficient is the raw count and NO opportunity ratio can reach it.** Every frame factor
— `eff_rna` vs `eff_junction` included — multiplies `E[S]`, so in log space it is an additive constant.
⛔ **`TRAPS.md` C0's trap is therefore structurally ABSENT here rather than avoided**: neither divisor
appears in the retained term at all, so there is nothing for them to disagree about, with either sign.
**The junction half's missing "frame statement" is answered: it does not need one.** Gate `C1`, brute-forced
against `scipy.stats.poisson`, with the `c·S` alternative firing it.

### 8b. ⛔⛔ AND THAT IS WHY THE TERM IS WRONG — the real precondition is `q`, not the frame

With `q` unknown, `q/(1−q)` spans `(0, ∞)`, so for **every** `f_g < 1` some `q` makes `E[S]` exactly `S`.
⭐⭐ **The profile likelihood in `f_g` is therefore EXACTLY FLAT on [0,1): the certified count carries
precisely ONE BIT about the unspliced split — "`f_g ≠ 1`".** Everything the `S` coefficient claims beyond
that bit is an assumption about `q`, not evidence. Gate `C2` (brute-force over 26 decades of `q/(1−q)`,
stated as convergence so a loosened tolerance cannot hide a tilt); `C2b` for the `S = 0` mirror.

⛔ **And `S = 0` is flat on the CLOSED interval, which inverts G6's own premise.** `test_vertex_reference`
G6 argues a silent gene's zero certified count is "the strongest possible evidence for the vertex". It is
not: with `q` free, `S = 0` is perfectly explained at `f_g = 1` too (take `q → 0`). The information lives
entirely in `S > 0` and points **away** from the vertex. So §7's vertex population gets nothing here either.

### 8c. THE MEASUREMENT THAT DECIDED IT — `scripts/design/certified_q_census.py`, all 36 conditions

`q` is an identity on the origin-split oracle: `q = S/(S + C_R)`. No solver runs.

| | |
|---|---|
| `q` mass-weighted median | **0.19 – 0.71**; **60–98 % of the mass below 0.9** in every condition |
| the term vs the uninformative reference | ⛔ **WORSE on 12 of 36 conditions** |
| worst | **+0.4578** mwae, `g90 ss0.50 capture_on`, where **98.3 %** of the certified-bearing mass sits on EDGEs whose truth is `f_g` = 0.84 and the term answers **0.04** |
| best | −0.4923, `g00 ss0.50 capture_on` |
| `capture_off` / `capture_on` | mean Δ −0.2937 (worse on 4/18) / −0.0289 (worse on **8/18**) |

⭐⭐ **The two zero controls are the diagnosis.** At `g00` the truth is all-RNA and the term is the best
result on the panel; at `g98` the truth is all-gDNA and it is the worst. **The sign of the effect is the
sign of the truth.** A channel whose benefit tracks the ANSWER rather than the evidence has not added
information — it has added a prior toward RNA (`TRAPS.md` A12's shape).

⭐ And `q` is RNA GEOMETRY, not composition: per-EDGE `q(g00)` vs `q(g50)` is Spearman **ρ = 0.9257** over
5,241 shared EDGEs. A 50-point swing in the gDNA fraction barely moves it.

⛔⛔ **This is `TRAPS.md` B18 honoured in its strongest form: the panel arm ran BEFORE a single line went
into `src/`.** It cost one script and no solver time, because the oracle cache makes it free.

### 8d. Q2 — "where does the flank-side evidence attach?" — **it already does, for the junction half**

Read as instructed, `node_total_density` already carries the junction half per-(EDGE, side) as
`rho_lo`/`rho_hi` (landed at `87a5468e`, `EQUATIONS.md` §3.6c). §2f's option (i) is **built**. The only
certified-RNA gap left on the flank axis is that `spliced_count` is in NEITHER flank total, which that
function's own docstring records as deliberate — a LEVEL change, and §3.6c's verdict is that the reframe
"is priced jointly with [the capture level residual] or not at all". It is not this build.

⚠ And §2b's rule is unchanged and still right: `spliced_count` lies on BOTH flanks, so it has no side to
choose — the λ term it would feed was always per-SLOT, never per-flank. §4's T3/T4 as written gate a
mechanism that the settlement removes; a `_lo`/`_hi` swap is already gated by `test_splice_flux_reframe`.

### 8e. ⛔ AND THE PRECISION THE PLAN REACHED FOR IS OUT OF CONTRACT TOO

§3b site 3 prescribes `tau_lam += density_factor_precision(cert, lam_grid)`, exactly as `tau_len` is
wired. That function reads `1/Var_λ` under the NORMALIZED factor — the Laplace precision, correct for a
PEAKED factor (gated at exactly 1.0 and 25.0, unchanged from `L = 6` to `L = 20`). The certified factor
is **monotone**: no peak, no scale, so its normalized variance is the WINDOW's. Measured, at `S = 1e4`:
**756.6 at `L = 6` → 0.106 at `L = 20`, a 7,100× swing.** ⭐ `simplex_logodds`' own stated acceptance test
is **L-invariance**; a τ that scales with `L` is disqualified by it. The honest form is analytic and
window-free — `I_cert = −∂²/∂λ²[S·log σ(−λ)] = S·f_g·(1−f_g)`, the same shape `strand_evidence` returns —
and gate `C5` pins both halves.

### 8f. ⭐ WHAT SURVIVES, AND WHAT THE CHANNEL IS BLOCKED ON

**Survives** (gate `C4`): ψ's reference plus the term is **exactly Beta(½, ½+S)** on `f_g` — the Jacobian
`df_g/dλ = σ(λ)σ(−λ)` supplies the missing −1 on each exponent. So §2d's headline is literally true: the
count IS the Bayes update of the reference and `S = 0` returns today's Beta(½,½) identically. Whatever is
built here needs no new normalisation and no `S = 0` branch. ⛔ It is the COEFFICIENT that is unlicensed,
not the algebra.

**Blocked on:** `1 − q` is the share of the crossing opportunity that fits inside the unbroken **EXONIC**
reach either side of the line. ⛔ That array does not exist. `splice_graph.build_contiguous_edge_reach_arrays`
is deliberately GENOMIC — its own docstring: *"a contiguous line is also crossed by nascent RNA, which is
genomic; taking the exonic reach here would declare an intronic nascent fragment impossible"* — and the
exonic version is per-TRANSCRIPT and abundance-weighted, which pass-0 does not know.
⭐ **This is the owner's domain call**, and the shape of it is: build a per-(EDGE, side) exonic-reach
opportunity (with its own brute-force enumeration gate, C0's standing rule), or leave the channel to the
mass accounting it already has. With `q` in hand the term becomes two-sided with a finite Fisher
information and an interior maximum at `1 − f_g = S·E_reach/((E_any − E_reach)·M)` — gate `C2`'s
perturbation shows exactly that, and it recovers the true split to 1e-2.

⭐ **The prize is worth the build if it can be made to work**: `edge_spliced` covers **85.3 %** of EDGE
unspliced mass off capture (§6), and where the term is licensed it is worth **−0.49 mwae** (`g00`).

⛔ **And the FIRST measurement of that build is a cheap one, so do it before writing any of it.** `q_geom`
= `1 − E_reach/E_any`, with `reach_lo`/`reach_hi` walked out from each contiguous EDGE to the nearest
junction endpoint or terminus and `crossing_eff_length` on the RNA pmf. Regress the measured `q` — which
`certified_q_census.py` already produces per EDGE — on `q_geom`. ⭐ If a purely annotation-derived
`q_geom` does not explain most of `q`'s variance, the abundance weighting is load-bearing, the quantity is
not available at pass-0, and the channel is closed rather than merely blocked. That is one script and no
solver time.
