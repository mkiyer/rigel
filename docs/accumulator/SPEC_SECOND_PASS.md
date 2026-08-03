# SPEC — the second pass: draining the side buffer

    Status: LIVE SPEC, 2026-08-02. ⭐ **P0–P4 ARE ALL CLOSED** — the pass is built, gated, wired into
    the pipeline and measured (`docs/WIP.md` P2.1, P2.2, P3, P4). ⭐⭐ **The anchor's bias against truth
    is now +0.00 % on the zero-gDNA falsification condition**, from −1.61 %. D-1/D-2/D-3/D-5/D-6 decided.
    ⭐⭐ **P0–P5 ARE ALL CLOSED; P6 is partial.** ⛔ **THE RESUME POINT IS NOT IN THIS FILE** — it is
    C3 / S4, junction opportunity, `TODO.md` rank 0. `docs/WIP.md` **B4** is why: calibration fits from the
    fragment-length POOLS, this pass made the *anchor* exact, and the RNA pool it feeds is
    opportunity-biased — so composition got 23.9 % worse. ⚠ §8's **D-4** is the one decision left here,
    and it is the same question as `TODO.md` rank 8 (a zero observation is a likelihood, not an
    impossibility) — both want the EXPOSURE that C3 supplies.
    Entry point: `docs/accumulator/PLAN_TWO_PASS.md` — this file is its **S3**, in full.
    Prior art it completes: `docs/accumulator/SPEC_GAP_PATHS.md` (pass 1, landed as S1), `docs/accumulator/PLAN_TWO_PASS.md` §5 (the
    sketch this replaces), `docs/accumulator/JUNCTION_OPPORTUNITY.md` (S4, which runs AFTER this).
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

`docs/accumulator/PLAN_TWO_PASS.md` §5.2 said `prior(h)` "needs transcript abundance for a spliced path and the gDNA
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
the reciprocal-opportunity design (`docs/accumulator/NODE_DENSITY_DERIVATION.md`) paying for itself: the one place the
model-free channel is exact is an edge, and the objects that distinguish these hypotheses are all edges.

⚠ **`count` is NOT interchangeable with `inv_length_sum` here.** A count is `rho * E[placements]`, so
comparing counts across objects compares different opportunity scales. Use `inv_length_sum`. The counts
are still worth carrying as the *weight of evidence* behind each density (§8, D-4).

### 1.2 ⭐ Why this does NOT couple to S4 (junction opportunity)

An obvious objection: `docs/accumulator/JUNCTION_OPPORTUNITY.md` says the junction pool is opportunity-biased, so surely a
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

⭐ **WHICH contiguous edges ∅'s evidence is — settled by the deposit rule, D-6.** A line is crossed iff it
lies **strictly inside a contiguous segment**, so over a fragment `[s, e)` the ∅ path crosses every line
in `(s, e)` and a path splicing `[a, b)` crosses those in `(s, a)` and `(b, e)`. The difference — and
therefore ∅'s evidence — is the lines at cuts **`a <= c <= b`, endpoints INCLUDED**. Both endpoints of an
annotated intron are cuts by construction, so those two lines always exist and always discriminate.
⛔ The scorer asked for `a < c < b` until 2026-08-02, which drops exactly those two and is **empty when
the intron spans one node**; `docs/WIP.md` P2.1 measured the damage and P2.2 fixed and gated it.

✅ **THE AGGREGATION IS `min`, THE BOTTLENECK — D-1 is closed and the parameter is deleted** (P2.2). A
path touches several objects and the molecule was present at **all** of them, so the path's rate is not
the sum; the scarcest object bounds it. Measured against the geometric mean over the 8 pilot conditions:
a different winner on **0.47–0.59 %** of held records, with the zero mask **bit-identical**. No accuracy
argument separates them at that scale, so the one with a derivation wins.

### 3.1a ⭐ HOW THE THREE FACTORS COMBINE — and the one rule that is not just "multiply"

`score(h) = rho(h) * f(L_h) * s(h)`, normalised over the candidate set and nothing else. ⚠ The product is
**not** a rate and not a calibrated likelihood — `rho` carries units of fragments per base while `f` and `s`
are probabilities — and it does not need to be, because the normalisation makes any common scale cancel.

⛔ **AN ALL-ZERO FACTOR IS UNINFORMATIVE, NOT DECISIVE** (P4.1, 2026-08-03). That same normalisation is
why: a factor taking the *same* value for every candidate cancels and cannot affect the answer. Zero is the
one value where the arithmetic loses that — it destroys the product instead of cancelling, the
normalisation cannot run, and the record collapses to a coin toss that throws away the factors that did
know. Measured: of 10 such records, the length term alone would have been right **8 of 8** times it could
decide. So an all-zero factor is dropped for that fragment.

⛔ **The PARTIAL-zero case is untouched** — a zero for *some* candidates says "no evidence for this path"
and stays decisive, which is D-3 exactly as ruled. ⭐ And when the product is zero everywhere because two
factors **contradict**, the coin toss is restricted to candidates whose `L` the length term has not ruled
out: `f = 0` is already the statement `max_fragment_length` makes, read off the measured distribution.
`rigel.second_pass.combine_factors` is the one place all of this lives.

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
behaviour** (`docs/WIP.md`, the strand-enumeration audit):

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

✅ **UNBLOCKED — there was no strand sign bug** (`docs/WIP.md` P0, 2026-08-02). `rna_sense_frac` is the
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

⭐ **S2.1 answered this** (`docs/WIP.md`). The rule is **order by identity, then one stream** — never hash
the fragment's content, because a content key ties exactly on the duplicates it would harm, so 100
identical fragments would draw identically and a 60/40 posterior would collapse to 100/0.

S1 already gives the deferred bank a **canonical order** (sorted on the record's own content, bit-identical
at any worker count), so `index in the drained queue` is well defined. Seed each draw from
`(global_seed, queue index)`.

✅ **DECIDED (P3): one stream in queue order**, which is S2.1's rule verbatim. The bank's canonical sort
makes "the i-th record" well defined, so `np.random.default_rng(seed)` consumed in that order is
reproducible by construction — and it vectorises, which is what makes the correspondence between stream
position and queue index a property of the code rather than of a loop body. ⛔ The per-index alternative
was rejected as machinery for a capability (drain a subset, re-draw one record) nothing needs.

⚠ `global_seed` is still a parameter of `drain`, not yet a config field — that is **P4**, where the pass is
wired in. ⛔ The gate is the one S2.1 established for
the whole pipeline and this pass must not weaken it: **byte-identical output at 1/2/4/8 workers**, and two
runs at one seed identical. `docs/accumulator/PLAN_TWO_PASS.md` §5.6 conceded that calibration "stops being bit-identical"
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

Deferred by owner ruling (`docs/accumulator/SPEC_GAP_PATHS.md` §10). They are a real source of long fragments on real data,
are not enumerable from the annotation, and ⚠ **the simulator cannot measure them at all** — so any
residual tail on real cfRNA that survives this pass is the first place to look, and the pilot will be
silent about it.

---

## §8 OPEN DECISIONS

| | | |
|---|---|---|
| ~~**D-1**~~ | ✅ **CLOSED (P2.2) and the parameter DELETED.** `min`, the bottleneck | §3.1. Measured on the 8 pilot conditions: the two pick a **different winner on 0.47–0.59 %** of held records and leave the **zero mask bit-identical**, so D-1 never interacted with D-3. No accuracy argument separates them at that scale, so the derivation decides. ⚠ Shares move > 0.10 on 1.9–3.4 % of records — the place to look if P4's tail gate misses |
| ~~**D-2**~~ | ✅ **DECIDED 2026-08-02 — NEITHER (a) nor (b) as written.** The drain is a pure function over the payload | §6.3, and it went against that section's own recommendation. ⛔ (a) forces the drain inside the scan, which bakes one draw into every cached payload and turns each P5/P6 re-measurement into a 61 s rescan of 8 conditions. ⭐ And (b)'s objection — "needs an `AccumulatorSet` binding with a payload export" — is false: every tally channel is already a `def_prop_ro` on the bound per-reference `Accumulator`. So the shape needs **no new C++**, drains a **cached** payload, and keeps the delta observable |
| ~~**D-3**~~ | ✅ **CLOSED 2026-08-02 — owner ruling: NO fallback.** A hard zero stays hard. Measured (P2.1): | 8 pilot conditions. The empty score vector is **≤ 1 %** of records; the zeros instead **decide** 56–80 % of them. D-3's population is the spliced arm: **17.9–22.8 %** of spliced hypotheses read zero, **100 % of them `annotated_empty`** (zero `unannotated`), and against truth **69–89 % are correctly zero** — nothing expressed the junction. ⛔ **The predicted self-exclusion is not the mechanism**: the heaviest empty slot has **1,509 held claimants and zero truth molecules**, so a pseudocount would hand those a route the library never used. ⚠ The ∅ arm's raw rate is **not data** — see D-6 |
| **D-4** | should the density carry its **weight of evidence** (the integer `count`) as well as its level? | A density of 0.01 from 1000 fragments and from 1 are not the same statement. The Beta-Binomial machinery calibration already uses is the honest route; ⚠ it is also scope. Measure whether it matters |
| ~~D-5~~ | ✅ **FIXED (P1).** A merged path whose supporters disagree about the strand is marked **AMBIGUOUS** | `docs/WIP.md` P1. Taking the first supporter's strand made the answer depend on GTF line order. ⚠ Unreachable on human data; fixed because the alternative is order-dependent |
| ~~**D-6**~~ | ✅ **FIXED AND GATED (P2.2).** `_lines_inside` → `_distinguishing_lines`; arm 3 of the P2 gate was verified failing first, and Y1/Y2/Y2b all fire | §3.1. `second_pass._lines_inside` asks for lines **strictly between** an intron's endpoints; the deposit rule says the lines distinguishing ∅ from a spliced path are those at cuts **`a <= c <= b`**, endpoints included — and both endpoints are always cuts, so the two guaranteed discriminators are exactly the two dropped, and the set is **empty when the intron spans one node**. Brute-forced against the reference. Measured: on gdna100-capture-off the shipped rule reads **0.4344** zero and the derived rule reads ⭐ **0.0000**. ⚠ Not a free-standing fix — on zero-gDNA the derived rate is still 0.72, which is **correct** (no gDNA ⇒ nothing crosses an intron contiguously), so it must be judged with D-3, not before it |

---

## §9 ⭐ THE PHASED EXECUTION PATH

⛔ **`P0`–`P6` HERE MEAN THIS FILE'S PHASES.** `docs/calibration/SOLVER_OBSERVABLES_PLAN.md` has an unrelated P0/P1/P2
(the solver-observables work, where P2 is the fragment-length likelihood and is gated OFF). They share
nothing but the letter. When citing one, name the file.

⚠ A **gate** here means: a test written *before* the code, verified failing first, whose target was fixed
in advance and cannot be tuned; then the code is deliberately broken to confirm the gate fires
(**perturbation**). A phase is not done until its gate passes and its perturbations fail.

| | phase | ⭐ the gate |
|---|---|---|
| ~~**P0**~~ | ✅ **CLOSED 2026-08-02 — there was no bug.** Two quantities were both called strand specificity; `StrandModel.strand_specificity` already recovers the simulated knob (1.00→1.0000, 0.75→0.7701, 0.50→0.5020). No code changed; the deliverable is the gate that stops it being filed a third time | a zero-gDNA library still reads `f_gdna ≈ 0.003`, not 0.499 — the S5.f-addendum measurement, re-run. ⭐ **And the exported scalar now means what it says**, checked against a library simulated at a known sense fraction. ⛔ Blocks §3.3 |
| ~~**P1**~~ | ✅ **CLOSED 2026-08-02.** D-5 fixed (a merged path whose supporters disagree is AMBIGUOUS); both audited behaviours gated through the scan; 3 perturbations, of which Z3 does not fire for a recorded reason | the audit's two cases as tests, through the **scan** path — ⚠ not the resolver, where `t_strand_arr_` is empty and the hypothesis strand silently reads NONE |
| ~~**P2**~~ | ✅ **CLOSED 2026-08-02 (P2.1 + P2.2).** All three steps done: **D-3 measured** (no fallback — owner ruling), the **gate built** (`tests/test_second_pass_scoring.py`, 7 mirror-paired loci, 10 tests, **no threshold anywhere**), and **D-1 measured and its parameter DELETED**. ⭐ Along the way the measurement found **D-6** — ∅'s evidence set was the wrong set of contiguous edges — and it is fixed and gated. ⚠ Two perturbations do not fire and are recorded: §3.2's pmf choice and ∅'s contested-union scoping are both still ungated | ⭐ **the discrimination is real** — ✅ and it is, on `rho` (arms 1/2), the **length** term (arm 5) and the **strand** term (arm 7), each isolated by a mirror rather than a constant. **4 of the 7 loci exist only because a perturbation passed.** D-1 and D-3 **measured**, not chosen |
| ~~**P3**~~ | ✅ **CLOSED 2026-08-02.** `second_pass.drain()` is a **pure function, payload in / payload out, with ZERO new C++** — D-2 decided against §6.3's recommendation because (a) breaks the scan cache, and because (b)'s stated export cost does not exist. 16 gates, 18 of 19 perturbations fire. ⭐ Two findings: the census would have depended on the RNG, and §6.2's bookkeeping conflicted with a payload door check — the door won. Plus two latent cache holes the new field exposed | ✅ conservation exact, bank empty, **1/2/4/8 workers and two runs at one seed byte-identical**, and ⭐ §6.1's claim checked directly: a drained fragment gives byte-identically the tally that offering only that hypothesis would have given |
| ~~**P4**~~ | ✅ **CLOSED 2026-08-02.** `_drain_side_buffer` runs between the scan and `build_fl_models`, so calibration sees the complete tally and runs once. Seed on `PipelineConfig.second_pass_seed`, ⭐ deliberately **not** `em.seed`. 5 gates, all 6 perturbations fire | ⭐⭐ **THE ANCHOR IS EXACT**: mean −1.61 % → **+0.00 %**, sd −1.48 % → **+0.02 %** against truth on zero-gDNA. ⭐ **The gDNA control holds EXACTLY** — all four pure pools move by **0**, and it is derived, not a tolerance: a drained fragment is either a multi-line crossing (no pool) or an annotated splice (`RNA_SPLICED`). ✅ **THE TAIL HOLDS EXACTLY as of P4.1** — **0** fragments above the true ceiling on all 8 conditions. It read 0–3 until the annihilation bug was fixed: an all-zero factor was destroying the other two, and the length term had already scored every impossible answer at zero |
| ~~**P5**~~ | ✅ **ANSWERED 2026-08-03, by measurement rather than by building the loop** | §7.1. ⭐ **NO ITERATION.** Three measurements settle it: (1) the anchor is **exact** after one drain (+0.00 % / +0.02 %), so there is nothing on it for a second pass to improve; (2) the held fragments' assigned lengths are right to **+0.12 bp** on average, so the assignment is not what moves the pool; (3) the pool's residual +6 % is a **selection** effect — the annotated-junction-using population genuinely is +3.8 % longer than the library (P4.2's subpopulation table) — and **iteration cannot fix a selection effect**. ⛔ Only S4's opportunity correction can, which is why `TODO.md` rank 0 is C3 and not a loop |
| **P6** | ⚠ **PARTIAL** — re-measure everything the hold-out invalidated | ✅ **D3's residual: 0 fragments above the true ceiling on all 8 conditions** (P4/P4.1). ✅ `docs/accumulator/PLAN_TWO_PASS.md` §2.4's table: superseded by `docs/WIP.md` B4. ✅ The FL subpopulations measured against truth (P4.2). ⛔ **STILL OPEN: every number in `docs/accumulator/JUNCTION_OPPORTUNITY.md` §3**, which is C3's own first task — its §1 derivation is untouched |

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
| the held population is **2–3.5 %** of a library, and systematically long | `docs/accumulator/PLAN_TWO_PASS.md` §2.3, 8 pilot conditions |
| `inv_length_sum` is an exact density **at an edge**, for any `f` | `docs/accumulator/NODE_DENSITY_DERIVATION.md`; the reference's own docstring |
| junction edges use the **same** quantum as contiguous edges | `accumulator.cpp`, the `sj_ids` deposit loop |
| an observed junction already pins the candidate strand | audited 2026-08-02, two-strand fixture, `+` → tP only, `−` → tM only |
| an unspliced fragment offers both strands | same fixture, no observed junction → 3 hypotheses |
| **0 of 404,168** human junction coordinates are annotated on both strands | `CLAUDE.md`; the index warns when it happens |
| `rna_sense_frac` is mirrored but consistently, so it cancels **today** | `docs/WIP.md` S5.f-addendum; forcing the nominal truth is 166× worse |
| the anchor is the library's own gDNA/RNA length mixture | C1; `deposited_lengths` is unconditional given deposit |
