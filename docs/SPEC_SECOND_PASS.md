# SPEC — the second pass: draining the side buffer

    Status: SPECIFICATION, 2026-08-02. Nothing here is built yet.
    Entry point: `docs/PLAN_TWO_PASS.md` — this file is its **S3**, in full.
    Prior art it completes: `SPEC_GAP_PATHS.md` (pass 1, landed as S1), `PLAN_TWO_PASS.md` §5 (the
    sketch this replaces), `JUNCTION_OPPORTUNITY.md` (S4, which runs AFTER this).
    Owner rulings, 2026-08-01/02, recorded in §8.

---

## §0 ⭐ THE RULE, IN ONE PAGE

Pass 1 holds every fragment whose gap has more than one surviving explanation — **2–3.5 % of a library**,
and systematically the **long** ones, because a longer gap admits more hypotheses. They are held WHOLE in
the side buffer, so nothing is lost, but nothing is tallied either. This pass resolves them.

> For each held fragment, score every hypothesis it arrived with, draw **one** by multinomial sample, and
> re-offer the fragment through the **same** `Accumulator::deposit` with that hypothesis alone. One
> hypothesis wins the whole fragment; no fractional deposit; integers stay integers.

```
    score(h)  =  rho(h)  x  f(L_h)  x  strand(h)
                  |            |            |
                  |            |            +-- library strand specificity; an OBSERVED junction pins it
                  |            +-- RNA pmf for a spliced path, the unconditional anchor for the genomic one
                  +-- the pass-1 accumulator's own density at the objects the path uses
```

⭐ **Every one of those three comes from pass 1 alone.** No transcript abundance, no calibration output,
no second calibration, no iteration. §1 and §2 are why that is possible, and they are the two decisions
this whole design turns on.

---

## §1 ⭐ THE PRIOR IS THE ACCUMULATOR ITSELF

`PLAN_TWO_PASS.md` §5.2 said `prior(h)` "needs transcript abundance for a spliced path and the gDNA
density for the empty path", and D-C proposed the EM's own abundances with a uniform first pass. ⛔ **Both
are worse than what is already on disk**, and the owner's 2026-08-02 ruling replaces them:

> *"If we have a fragment that maps to different intron combinations, we would extract the nodes/edges
> from the accumulator for each candidate mapping. We can use the first pass accumulator state
> density/count information … So we may not need shortcuts here."*

A hypothesis is a **path**, and a path is a set of accumulator objects. Pass 1 already counted how many
molecules used each of those objects. That is the prior, and it is local, evidence-based, and complete
before the second pass starts.

### 1.1 Why the densities are directly comparable — and why this needs no model

The competing hypotheses disagree about **one thing**: did the molecule cross the gap *contiguously*, or
*jump* it via a junction?

* ∅ says it crossed the gap's contiguous edges. Evidence: `edge_unspliced_inv_length_sum[e]`.
* `h` says it used junction `J`. Evidence: `sj_inv_length_sum[J]`.

⭐ **Both objects are 0-bp lines, and both are deposited with the SAME quantum** — `inv_length_quantum(L−1)`,
junction edges included (`accumulator.cpp`, the `sj_ids` loop). And at an edge

```
    E[inv_length_sum]  =  rho * E_f[placements * (1/placements)]  =  rho     exactly, for ANY f
```

so both are estimates of the same physical quantity, a **start density in fragments per base**, on the
same scale, with the opportunity cancelled identically and **no fragment-length model anywhere**. This is
the reciprocal-opportunity design (`NODE_DENSITY_DERIVATION.md`) paying for itself: the one place the
model-free channel is exact is an edge, and the objects that distinguish these hypotheses are all edges.

⚠ **`count` is NOT interchangeable with `inv_length_sum` here.** A count is `rho * E[placements]`, so
comparing counts across objects compares different opportunity scales. Use `inv_length_sum`. The counts
are still worth carrying as the *weight of evidence* behind each density (§8, D-4).

### 1.2 ⭐ Why this does NOT couple to S4 (junction opportunity)

An obvious objection: `JUNCTION_OPPORTUNITY.md` says the junction pool is opportunity-biased, so surely a
junction's number needs correcting before it can be compared to an edge's?

**No, and the distinction matters.** S4's subject is the **fragment-length distribution** measured on
junction-crossing fragments — `A_j(w)` varies with `w`, so `f` is length-biased. `rho` is not affected:
the `1/(L−1)` deposit weight cancels the opportunity for every `w` identically. So this pass reads `rho`
from junctions safely and leaves `f` on junctions entirely alone. ⛔ **S4 still runs after this pass**,
because S4 needs the drained tally, not the other way round.

### 1.3 What was rejected, and why

| rejected | why |
|---|---|
| **uniform transcript abundance** | breaks the circularity by discarding the evidence. Owner: *"it's actually very easy to use a form of the transcript abundances!"* |
| **the EM's own abundances** (D-C) | ⛔ circular: calibration feeds the EM, and the drained tally feeds calibration |
| **`CalibrationResult`'s deconvolved densities** | ⭐ tempting — `mass_{gdna,rna}_{node,edge} / *_eff_len` is exactly a per-component density, opportunity-normalised, already emitted. ⛔ But it is a calibration **output**, and this pass's product is calibration's **input**. Using it forces a second calibration, which is the loop §7.1 forbids. §2 is the way out |

---

## §2 ⭐ WHERE THE DRAIN SITS — the structural decision

```
    scan  ──►  payload₁ (tally + side buffer)  ──►  SELECT  ──►  DRAIN  ──►  payload₂  ──►  calibrate  ──►  EM
                          │                            │
                          └── FL models, strand model ─┘        ⭐ calibrate runs ONCE, on the whole tally
```

**The drain runs between the scan and calibration.** That is what makes §1 possible, and it rests on one
observation:

> ⭐ **The second pass picks a PATH, not a COMPONENT.** The accumulator stores by structure — contained,
> crossing, spliced, junction — never by gDNA-vs-RNA. Deconvolution is calibration's job and happens
> afterwards. So the drain never needs to know whether a molecule is gDNA or RNA; it needs to know which
> objects it touched, and that is a structural question answerable from pass 1.

Everything the selection needs is available at that point:

| input | source | needs calibration? |
|---|---|---|
| the held fragments and their hypotheses | `payload.deferred` (S1) | no |
| `rho` at every object | `payload.{edge_unspliced,sj}_inv_length_sum` | no |
| `f_rna`, `f_gdna` | `build_fl_models(payload)` — payload only, since C2 | no |
| the unconditional anchor | `payload.deposited_lengths` | no |
| library strand specificity | the scan's `StrandModels` | no |
| the partition (to map an intron → junction id) | `build_junction_edge_arrays(index)` | no |

⚠ **The junction CSR is NOT in the payload** — it carries `ref_sj_offsets` and the banks, but not
`sj_offsets`/`sj_acceptor_cut`, so an intron cannot be resolved to a junction id from the payload alone.
The selection takes it from the index, exactly as `calibrate` already does. That is a dependency to wire,
not a gap to close.

---

## §3 THE SCORE

For held fragment `x` with hypotheses `h ∈ H(x)`:

```
    score(h)  =  rho_hat(h)  *  f_h(L_h)  *  s_h
```

and the three factors are §3.1–§3.3. The scores are normalised over `H(x)` and nothing else; they are a
**relative** statement about one fragment, never a rate.

### 3.1 `rho_hat(h)` — the density term

Let `G` be the fragment's unsequenced gaps (pass 1's own gap list, minus the gaps the CIGAR explained).

* **`h = ∅`** — the molecule crosses `G` contiguously. Its evidence is the *unspliced* crossing density on
  the contiguous edges under `G`.
* **`h` spliced** — the molecule uses junction edges `J_h`. Its evidence is those junctions' density.
  ⚠ A path may also cross contiguous edges *between* its introns (the owner's §1 example: a path can
  cross an exon no read touched). Those are part of the path and count.

⛔ **OPEN — D-1, the aggregation.** A path touches several objects and the molecule must be present at
**all** of them, so the path's rate is not the sum. The candidates are `min` (the bottleneck: a molecule
crossing the whole gap is present at every edge, so the scarcest object bounds it) and the geometric mean.
⚠ `min` is the one with an argument behind it; it is also the one that will be dominated by a single
zero. **Measure both before choosing** — this is a tunable if picked by taste, and a derivation if picked
by measurement.

### 3.2 `f_h(L_h)` — the length term

Owner ruling, 2026-08-02:

> *"Any of the fragments with introns are automatically RNA, so the likelihood will be checked against the
> RNA fragment length distribution. Any of the fragments that are unspliced with the full genomic span
> could be DNA or RNA. And that's the tricky part."*

| hypothesis | pmf | why |
|---|---|---|
| spliced | `fl_models.rna_pmf` | a spliced molecule is certified RNA — gDNA cannot splice |
| `∅` | ⭐ `payload.deposited_lengths`, the **unconditional anchor** | the component is unknown, so it must be *marginalised* — and the anchor is exactly that marginal, measured: every deposited fragment binned at its own `L`, no purity condition, i.e. the library's actual gDNA/RNA blend |

⭐ This is a proper likelihood rather than a fudge. For a spliced path the component is known and we
condition on it; for `∅` it is unknown and we integrate over the library's own composition. And the anchor
is a payload field, so §2's no-calibration property survives.

⚠ **The alternative worth measuring later**: a *locally* composition-weighted mixture instead of the
global anchor. More power, but it reintroduces the calibration coupling. ⛔ Not in v1.

⚠ **The anchor is short-biased and this pass cannot fix that.** The held fragments are the long ones and
they are excluded from the anchor by construction. That is a known, bounded bias — recorded here so it is
not rediscovered as a surprise — and it is one of the reasons §7.1 forbids iterating.

### 3.3 `s_h` — the strand term

Owner ruling, 2026-08-02, and ⭐ **audited 2026-08-02: the first half is already the implemented
behaviour** (`LEDGER.md`, the strand-enumeration audit):

> *"Splice junctions are stranded and asymmetric. So if we detect one splice junction in a fragment, we
> know this fragment strand … Fragments without a splice junction are unspliced and could be either
> strand, but the library strand specificity could constrain this."*

| the fragment | how the strand enters |
|---|---|
| has an **observed** junction | ⭐ the strand is **pinned**, and pass 1 already enforces it: measured, a fragment whose motif is `+` is offered only the `+` transcript's gap intron, and `−` only the `−` one. `s_h` is then constant across `H(x)` and **cancels** |
| **unspliced** | ⛔ both strands are live — measured: 3 hypotheses (`+` intron, `−` intron, `∅`). `s_h` is the discriminator, and without it these are separated by length alone |

For an unspliced fragment, `s_h` compares the fragment's `align_strand` against the strand `h` implies
(`hypothesis.sj_strand`, already carried in the deferred record):

* `h` spliced ⇒ certified RNA ⇒ `s_h = P(align_strand agrees | RNA)` from the library's sense fraction.
* `h = ∅` ⇒ the molecule may be gDNA, which is **symmetric** ⇒ `s_h` is strand-uninformative.

✅ **UNBLOCKED — there was no strand sign bug** (`LEDGER.md` P0, 2026-08-02). `rna_sense_frac` is the
Beta posterior mean of `p_r1_sense = P(align_strand == the junction's strand)`, correctly named and
correctly valued; the apparent mirror was a collision with the simulator's *direction-agnostic* protocol
fidelity. ⭐ So the second pass can consume it directly — it is exactly the
`P(align_strand agrees | RNA)` this term needs.

⚠ **But mind the direction.** On an R1-antisense library `rna_sense_frac ≈ 0.01`, so an RNA hypothesis
whose implied strand *disagrees* with `align_strand` is the LIKELY one. ⛔ A scorer written as
"agreement ⇒ multiply by `rna_sense_frac`" would be exactly backwards on every real cfRNA library. Use
the value as the probability of agreement it is, and let it be small.

---

## §4 ⭐ THE CENSUS ALREADY SAYS WHICH TERM DOES THE WORK

S1's umbrella census partitions the held set by *which question is open* — which turns out to be exactly
*which factor discriminates*. It was not designed for this; it falls out.

| class | what differs | the discriminating factor |
|---|---|---|
| `gap_deferred_rna_or_gdna` | `∅` vs one spliced path | ⭐ `rho`: contiguous-crossing density vs the junction's, plus the FL term across two different pmfs |
| `gap_deferred_which_introns` | ≥2 spliced paths, no `∅` | certified RNA either way, so the component cancels: **junction `rho` ratio + `f_rna`** |
| `gap_deferred_both` | both at once | all three |

⚠ Use this to **read the diagnostics**, not to branch the code. One scorer, evaluated over `H(x)`; the
classes tell you where a regression came from, and give the per-class before/after the gates in §9 need.

---

## §5 SELECTION

### 5.1 The draw

One multinomial draw per fragment over the normalised scores. **One hypothesis wins the whole fragment.**
No fractional deposit — integers are preserved and no transcript is fractionated (owner ruling, and it is
what keeps the accumulator's integer-only contract intact).

### 5.2 ⭐ Seeding — settled, and it needs no new machinery

⭐ **S2.1 answered this** (`LEDGER.md`). The rule is **order by identity, then one stream** — never hash
the fragment's content, because a content key ties exactly on the duplicates it would harm, so 100
identical fragments would draw identically and a 60/40 posterior would collapse to 100/0.

S1 already gives the deferred bank a **canonical order** (sorted on the record's own content, bit-identical
at any worker count), so `index in the drained queue` is well defined. Seed each draw from
`(global_seed, queue index)`.

⚠ `global_seed` lives on the config alongside `EMConfig.seed`. ⛔ The gate is the one S2.1 established for
the whole pipeline and this pass must not weaken it: **byte-identical output at 1/2/4/8 workers**, and two
runs at one seed identical. `PLAN_TWO_PASS.md` §5.6 conceded that calibration "stops being bit-identical"
after this pass — ⛔ **that concession is withdrawn**; reproducibility is achievable here and is gated.

---

## §6 THE DRAIN

### 6.1 One tally path

⭐ The drain re-enters `Accumulator::deposit` with the chosen hypothesis **alone** — a set of size one, so
arbitration is degenerate and it deposits (or is rejected `TOO_LONG` / `EMPTY` by the ordinary rules). There
is no second deposit implementation, no duplicated crossing logic, and byte-identity with
`tests/native/_accumulator_reference.py` is preserved for free. The **specification gains a selection
step, not a second tally** — the reference gets the same `drain()`, so parity stays cheap.

### 6.2 The bookkeeping, and the conservation identity it must satisfy

The drain **consumes** the bank. After it:

```
    deferred bank                 empty
    qc.drained_deposited  +  qc.drained_dropped_*   ==   qc.deferred_undetermined_gap   (pass 1's count)
    qc.deposited(final)   ==   qc.deposited(pass 1)  +  qc.drained_deposited
```

⚠ **New counters, not reused ones.** Folding the drained fragments into `deposited` alone would make the
pass-1 and post-drain numbers indistinguishable, and every before/after measurement in §9 needs them apart.
⭐ And `deferred_undetermined_gap` is **kept** as pass 1's count rather than zeroed — it is the denominator
the drain's own conservation is checked against.

### 6.3 Where the code lives

**Selection in Python, drain in C++.** The scores need the FL models and the strand model, which are
Python-side; the drain must re-enter `deposit()`.

⛔ **OPEN — D-2, the plumbing.** The accumulators are owned by `BamScanner` and `build_result` gathers them
once, inside `scan()`. Two shapes:

| | shape | trade |
|---|---|---|
| **(a)** | keep the scanner alive; add `scanner.drain(choices)` re-entering the live `AccumulatorSet`, then re-gather | one tally, one gather, no new export. ⚠ `scan_and_buffer` must return the scanner or a handle |
| **(b)** | build a **second** `AccumulatorSet` from the index, drain into it, `merge_from` the first | ⭐ the merge is already exact and gated; the drain's **delta is observable**, which is a far better diagnostic than a changed total. ⚠ needs an `AccumulatorSet` binding with a payload export — a second gather unless `build_result` is factored out |

⭐ Recommend **(a)** for v1, with the delta recovered by differencing the two payloads, and (b) noted as
the refactor if the delta ever needs to be first-class. Decide before P3.

---

## §7 WHAT THIS DELIBERATELY DOES NOT DO

### 7.1 ⛔ It does not iterate — owner ruling, and it is a safety property

The confident set is biased **short**. A short-biased FL prior prefers the shorter `L`, which is the
**more-spliced** path; those fragments then re-enter the fit shorter still. **That loop can run away.**

⭐ **Broken by construction: fit once, score once, drain once, stop.** Then *measure* whether a second
iteration would move anything (P5). ⛔ Designing the loop before knowing it is needed is how a tuned
constant gets in.

### 7.2 It does not calibrate twice

§2. Calibration runs once, on the complete tally.

### 7.3 It does not touch junction opportunity

§1.2. S4 runs after, on the drained tally.

### 7.4 It does not resolve unannotated junctions

Deferred by owner ruling (`SPEC_GAP_PATHS.md` §10). They are a real source of long fragments on real data,
are not enumerable from the annotation, and ⚠ **the simulator cannot measure them at all** — so any
residual tail on real cfRNA that survives this pass is the first place to look, and the pilot will be
silent about it.

---

## §8 OPEN DECISIONS

| | | |
|---|---|---|
| **D-1** | how to aggregate `rho` over the several objects a path uses — `min` (bottleneck) or geometric mean | §3.1. ⛔ **Measure both.** Picked by taste it is a tunable; picked by measurement it is a derivation |
| **D-2** | the drain's plumbing: live scanner (a) or second set + merge (b) | §6.3. ⭐ Recommend (a) for v1 |
| **D-3** | ⛔ **a junction with ZERO flux** — impossible, or merely unobserved? | A hard zero makes the hypothesis unselectable and can empty the score vector. ⛔ A pseudocount is a magic number and will not be invented. ⚠ **Measure first**: the held fragments are excluded from the tally they are scored against, so a junction used *only* by deferred fragments reads zero — how often that happens is the fact that decides this |
| **D-4** | should the density carry its **weight of evidence** (the integer `count`) as well as its level? | A density of 0.01 from 1000 fragments and from 1 are not the same statement. The Beta-Binomial machinery calibration already uses is the honest route; ⚠ it is also scope. Measure whether it matters |
| ~~D-5~~ | ✅ **FIXED (P1).** A merged path whose supporters disagree about the strand is marked **AMBIGUOUS** | `LEDGER.md` P1. Taking the first supporter's strand made the answer depend on GTF line order. ⚠ Unreachable on human data; fixed because the alternative is order-dependent |

---

## §9 ⭐ THE PHASED EXECUTION PATH

⛔ **`P0`–`P6` HERE MEAN THIS FILE'S PHASES.** `SOLVER_OBSERVABLES_PLAN.md` has an unrelated P0/P1/P2
(the solver-observables work, where P2 is the fragment-length likelihood and is gated OFF). They share
nothing but the letter. When citing one, name the file.

⚠ A **gate** here means: a test written *before* the code, verified failing first, whose target was fixed
in advance and cannot be tuned; then the code is deliberately broken to confirm the gate fires
(**perturbation**). A phase is not done until its gate passes and its perturbations fail.

| | phase | ⭐ the gate |
|---|---|---|
| ~~**P0**~~ | ✅ **CLOSED 2026-08-02 — there was no bug.** Two quantities were both called strand specificity; `StrandModel.strand_specificity` already recovers the simulated knob (1.00→1.0000, 0.75→0.7701, 0.50→0.5020). No code changed; the deliverable is the gate that stops it being filed a third time | a zero-gDNA library still reads `f_gdna ≈ 0.003`, not 0.499 — the S5.f-addendum measurement, re-run. ⭐ **And the exported scalar now means what it says**, checked against a library simulated at a known sense fraction. ⛔ Blocks §3.3 |
| ~~**P1**~~ | ✅ **CLOSED 2026-08-02.** D-5 fixed (a merged path whose supporters disagree is AMBIGUOUS); both audited behaviours gated through the scan; 3 perturbations, of which Z3 does not fire for a recorded reason | the audit's two cases as tests, through the **scan** path — ⚠ not the resolver, where `t_strand_arr_` is empty and the hypothesis strand silently reads NONE |
| **P2** | ⚠ **PARTIAL — built, UNGATED.** `src/rigel/second_pass.py` exists and runs; `Accumulator.length_under` exposes `L` so the scorer does not reimplement it. ⛔ **No test covers the module.** Remaining, in order: **measure D-3 on the pilot**, then a discriminating fixture + the gate, then measure D-1 and delete the parameter | ⛔ **the discrimination is real**: on a hand-built locus where the truth is known, the correct hypothesis takes the larger share. Plus D-1 and D-3 **measured**, not chosen |
| **P3** | the drain: re-enter `deposit()`, the new counters, the bank consumed | ⭐ **conservation** — `drained_deposited + drained_dropped_* == deferred_undetermined_gap`, exactly, and the bank is empty. Byte-identity C++↔reference. **Byte-identical at 1/2/4/8 workers and across two runs at one seed** (§5.2) |
| **P4** | wire into the pipeline, ahead of calibration | ⭐ **THE TAIL** — no fragment above the library's true longest molecule, **713 bp** on the pilot, read from `truth_fragment_lengths.tsv`. ⛔ **And the gDNA control must not move one digit** — gDNA has no introns to miss, so any movement there means the drain reached fragments it cannot explain |
| **P5** | ⭐ **measure** whether a second iteration would move anything | §7.1. A measurement, not a feature. If it would move things, that measurement tells you what damping and what stopping rule are needed — and only then is a loop designed |
| **P6** | re-measure everything the hold-out invalidated | ⛔ D3's residual above the ceiling; **every** number in `JUNCTION_OPPORTUNITY.md`; `PLAN_TWO_PASS.md` §2.4's anchor-vs-pool table. All are statistics of the deposited set and all currently describe a tally with 2–3.5 % held out |

Then, and only then: **S4** (junction opportunity, on the drained tally) and **S5** (regenerate the
goldens — once, at the very end, twice, and diff).

⛔ **P0 → P1 → P2 is a hard order.** P2's strand term is wrong until P0 lands, and P2's fixtures are
untrustworthy until P1's scan-path gate exists.

---

## §10 THE MEASUREMENTS THIS SPEC IS BUILT ON

Everything asserted above that is a number, and where it came from — so the next session can re-run rather
than trust:

| | |
|---|---|
| the held population is **2–3.5 %** of a library, and systematically long | `PLAN_TWO_PASS.md` §2.3, 8 pilot conditions |
| `inv_length_sum` is an exact density **at an edge**, for any `f` | `NODE_DENSITY_DERIVATION.md`; the reference's own docstring |
| junction edges use the **same** quantum as contiguous edges | `accumulator.cpp`, the `sj_ids` deposit loop |
| an observed junction already pins the candidate strand | audited 2026-08-02, two-strand fixture, `+` → tP only, `−` → tM only |
| an unspliced fragment offers both strands | same fixture, no observed junction → 3 hypotheses |
| **0 of 404,168** human junction coordinates are annotated on both strands | `CLAUDE.md`; the index warns when it happens |
| `rna_sense_frac` is mirrored but consistently, so it cancels **today** | `LEDGER.md` S5.f-addendum; forcing the nominal truth is 166× worse |
| the anchor is the library's own gDNA/RNA length mixture | C1; `deposited_lengths` is unconditional given deposit |
