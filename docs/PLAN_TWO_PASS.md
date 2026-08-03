# THE PLAN — one two-pass structure, serving both open problems

    Status: ⚠ **HISTORY, as of 2026-08-03.** Its S1, S2 and S3 have ALL LANDED — the side buffer, its
      persistence, and the second pass that drains it. ⭐ Read it for §0–§1, which are why the gap-path and
      junction-opportunity problems are ONE problem and still the clearest statement of that.
      ⛔ **§5 is superseded by `SPEC_SECOND_PASS.md`** and **§2.4's anchor-vs-pool table by `LEDGER.md` B4**,
      which measured the drained numbers §2.4 said could not be read yet. S4 is now `TODO.md` rank 0.
      Originally: CONVERGENCE PLAN, 2026-08-01. This supersedes the separate tracks.
    Replaces as the entry point: `SPEC_GAP_PATHS.md` (the gap-path track) and
    `JUNCTION_OPPORTUNITY.md` §7–§8 (the opportunity track). Both remain valid as
    DERIVATIONS; this file is the order of work and the gates.
    ⭐ Written to be handed to ONE session with clean context. Two sessions were working
    these as separate problems; §1 is why they are not.

---

## §0 ⭐ THE ONE-PAGE VERSION

Two problems were being worked in parallel:

| | problem | was being called |
|---|---|---|
| **A** | a fragment's unsequenced mate gap may hold no intron, one, or several — and *which* is a likelihood question | the intron / gap-path problem |
| **B** | the RNA fragment-length distribution is measured from junction-crossing fragments, which is a length-biased sample | junction opportunity (C3) |

⭐ **They are the same structure.** Both need a quantity that does not exist until a first pass is over:
A needs a fragment-length distribution to score paths; B needs transcript abundance to weight the
opportunity. Neither can be answered in one pass, and **both are answered by the same two-pass shape**.

> **Pass 1** deposits every fragment whose answer is determined, and holds the rest **whole** in a side
> buffer. **Pass 2** fits on the confident set, scores the held-back fragments, picks one answer each by
> **multinomial sample**, and drains them through the *same* deposit path.

⛔ **Do not build two two-pass machines.** Build one; B rides on it.

### The owner's rulings this plan encodes

* *"we should be searching for gap introns within every fragment"* — landed (C2.6).
* *"a fragment compatible with multiple transcripts is essentially the same thing as deciding which of
  several inferred intron combinations is most likely. Every transcript has exactly one exon-intron
  combination."* — ⭐ this is why A and B share a latent variable, §1.
* *"a simple multinomial sampler will be sufficient. Nondeterministic but correct and doesn't randomly
  fractionate transcripts."* — §5.
* *"first pass and then draining the side buffer, assigning likelihoods, and doing multinomial sampling
  in a quick second pass"* — ⛔ **no iteration in v1**, §5.4.
* *"our algorithm and our code is fine … this is a small fraction of the total problem … all calibration
  does is provide a prior for the EM."* — the measured share is **2–3.5 %**, §2.3. Build proportionately.

---

## §1 WHY A AND B ARE ONE PROBLEM

A transcript fixes exactly one exon/intron structure. So a fragment's candidate transcript set **is** its
set of possible splice configurations — *"which transcript?"* and *"which intron combination?"* are one
latent variable seen from two sides.

⚠ The map is **many-to-one**: several isoforms can imply the *same* intron set across one fragment's span
while differing elsewhere. What a fragment really has is a **partition of its candidate set by implied
configuration**, and the accumulator already acts on exactly that distinction:

* candidates disagree on the transcript but **agree on the intron set** → `L` determined → **deposits**;
* candidates **imply different intron sets** → `L` undetermined → **deferred**.

⭐ So the architecture is already aligned. What is missing is the **resolution step** for the
more-than-one-block case — and that step needs a fragment-length model, which is B's subject. Hence one
plan.

---

## §2 WHERE THE WORK ACTUALLY IS — measured 2026-08-01, not assumed

### 2.1 Landed and committed

| | |
|---|---|
| `40006126` | **C2** — one definition of fragment length; the scanner's rival histogram deleted; `build_fl_models(payload)` takes the payload and nothing else |
| `e2f41cd0` | **the junction-opportunity derivation** (`JUNCTION_OPPORTUNITY.md`) and the gap-intron spec |

### 2.2 ✅ S1 AND S2 HAVE LANDED — the migration is finished and the side buffer exists

**The suite is 21 failed / 1868 passed, and all 21 are `test_golden_output`** — the clean baseline.
`LEDGER.md`'s **S1 + S2** entry is the record. What was half-done at `750cc8ee` is done:

| was | now |
|---|---|
| the reference had moved to the hypothesis interface; the C++ had not | ⭐ byte-identical over **every** `Tally` field, the deferred CSR and the umbrella census included |
| the binding did not expose the deferred bank | `Accumulator.deferred` / `.gap_resolution`, and `build_result` gathers both through the same exporters |
| the side buffer existed in C++ and nowhere else | it crosses the ABI, lives on `AccumulatorPayload`, and survives the scan cache byte-identically (**S2**) |
| old call sites passed `introns=` | migrated, along with the two design scripts |

⭐ **Two things S1 found that were not on anyone's list**, both in the LEDGER entry:

1. **`GapResolution.RESOLVED_UNSPLICED` could not be entered by any fragment** — the genomic path is
   always the longest, so it can never be the sole survivor. Deleted; the *ordering* is pinned instead.
2. ⛔ **The reference stamp was untestable.** Two perturbations — a constant `ref` in the accumulator, and
   a constant in `AccumulatorSet` — each passed **all 1860 tests**, because every fixture was
   single-contig or deferred only on reference 0. Two fixtures moved to close it.

⚠ **And the owner's known one-line defect is fixed and gated**: `payload_schema_digest` recurses into the
nested banks now. It had no gate until a perturbation said so.

⭐ **The resume point is S3, the second pass** (§5). ✅ **D-D is answered** (§7): the second pass needs no
new seed machinery — S2.1 established the rule as *order by identity, then one stream*, and S1 already
gives the deferred queue its canonical order.

⚠ `C2.6` (cut a gap intron on every fragment, not only unspliced ones) **did land and is measured** —
`LEDGER.md`. It is the reason `dropped_too_long` is now 0 and the anchor's sd error fell from +27 % to
+2 %. The in-flight work is the *generalisation* of it from "one intron per gap" to "a hypothesis set".

### 2.3 The size of the problem, measured on the 8 pilot conditions

| library | deferred | share of offered | `dropped_too_long` |
|---|---|---|---|
| zero-gDNA, capture off | 172,242 | **3.44 %** | ⭐ **0** |
| zero-gDNA, capture on | 104,745 | 2.09 % | **0** |
| gdna100, capture off | 172,771 | 1.73 % | **0** |
| gdna100, capture on | 108,359 | 1.08 % | **0** |

⭐ **Nothing is discarded any more.** ⚠ The gdna100 rate is lower only because gDNA fragments have no gaps
to be ambiguous about — they dilute the denominator; it is not evidence that gDNA libraries are easier.

### 2.4 ⛔ B's target CANNOT be read right now, and this is a trap

Measured this session against the clean ruler, zero-gDNA, vs the simulator's truth:

| | mean | sd |
|---|---|---|
| anchor `deposited_lengths` | −1.6 % | −1.5 % |
| RNA_SPLICED pool, raw | +2.4 % | −6.4 % |
| pool corrected with the **TRUE** θ | **−3.2 %** | −3.3 % |

Read naively: the opportunity correction overshoots. ⛔ **Do not bank that.** The 3.44 % deferred
fragments are held out of **both** the anchor and the pool, and they are systematically the **long** ones
— a longer gap admits more hypotheses. Everything above is measuring the hold-out.

⚠ **And it retires an earlier result.** Before C2.6 the correction appeared to take the pool from +8.0 %
to +0.3 % — near-perfect. That was **two errors cancelling**: uncut introns inflated the pool, and the
opportunity correction happened to undo the inflation. **B's real target is only knowable after the side
buffer drains.** That is why B is sequenced after A and not before.

---

## §3 THE ARCHITECTURE

```
   PASS 1  ─ every fragment arrives with its HYPOTHESIS SET, not one path and a flag
             (the empty set is the genomic hypothesis)
             │
             ├─ exactly one hypothesis survives  ──▶  DEPOSIT it
             └─ two or more survive              ──▶  SIDE BUFFER, whole
                                                       │
   fit the fragment-length model on the DEPOSITED set only ─┐
                                                       │    │
   PASS 2  ─ score each held hypothesis:  P(L | FL) × prior(path)
             pick ONE by multinomial sample
             drain through the SAME deposit()  ◀───────┘
```

⭐ **One tally path.** The drain re-enters `Accumulator::deposit`; there is no second deposit
implementation, so byte-identity to `tests/native/_accumulator_reference.py` is preserved for free. The
specification gains a **selection** step, not a second tally.

⭐ **The accumulator decides, the adapter only supplies.** Today the scanner collapses the answer into a
bool before handing it over, which is why the accumulator "cannot decide". Give it the set.

---

## §4 THE STEPS, AND WHAT EACH ONE HAS TO PROVE

⚠ A **gate** here means: a test written *before* the code, verified failing first, whose target was fixed
in advance and cannot be tuned. Then the code is deliberately broken to confirm the gate fires
(**perturbation**). A step is not done until its gate passes and its perturbations fail.

| | step | ⭐ the gate, in one line |
|---|---|---|
| ~~**S1**~~ | ✅ **LANDED** — C++ to the reference, deferred bank across the binding, old call sites | ⭐ **conservation** — `deposited + deferred + dropped_* == offered`, exactly. Plus byte-identity C++↔reference and the same answer at 1/2/4/8 workers. All met; 9 perturbations, two of which found real holes |
| ~~**S2**~~ | ✅ **LANDED** — the side buffer through `scan_payload` and the scan cache | the queue survives a cache round-trip **byte-identically** and comes back TYPED; its merge is deterministic at any worker count. ⚠ `payload_schema_digest` now recurses into the nested banks, so every existing cache is invalidated by design |
| **S3** | ⭐ **the second pass** — ⭐ **SPEC'D IN FULL: `docs/SPEC_SECOND_PASS.md`**, which supersedes §5 below and carries the phased path P0–P6 | ⭐ **the tail** — no fragment longer than the library's true longest molecule. On the pilot that is **713 bp**, read from `truth_fragment_lengths.tsv`. ⛔ And the **gDNA control must not move one digit** |
| **S4** | **B** — the junction-opportunity correction, on the drained tally (§6) | on a zero-gDNA library the corrected RNA distribution must agree with the anchor, because there they describe **one population** |
| **S5** | regenerate `tests/golden/` | ⛔ regenerate **once**, at the very end, **twice, and diff** |

⛔ **S1 must NOT be judged by a calibration A/B.** Between S1 and S3 the tally is deliberately *thinner* —
the ambiguous mass is held, not yet deposited. Its gate is **conservation**, not accuracy.

---

## §5 THE SECOND PASS — the design

### 5.1 Fit

Fit the fragment-length model on **directly deposited fragments only** — the set with exactly one
surviving hypothesis. That is what the pools already are; no new estimator.

### 5.2 Score

For each held fragment, for each surviving hypothesis `h`:

```
    score(h)  =  P(L_h | FL)  ×  prior(h)
```

`L_h` is that hypothesis's implied length. `prior(h)` needs transcript abundance for a spliced path and
the gDNA density for the empty (genomic) path. ⭐ Both are quantities **B also needs** — one shared
source, per §0.

### 5.3 ⭐ Assign — multinomial, and the seeding is the part that will bite

One draw per fragment from the normalised scores. **One hypothesis wins the whole fragment.** No
fractional deposit: integers are preserved and no transcript is fractionated. This is also what makes it
compatible with the accumulator's standing integer-only ruling.

⛔ **Seed per fragment by its INDEX IN THE MERGED QUEUE — not from a global or per-worker stream.**

The gate is bit-identity at 1/2/4/8 workers. A stream RNG makes the draw depend on scheduling and that
gate fails immediately. The merged deferred queue is already required to be deterministically ordered, so
seed each draw from `(global_seed, index_in_sorted_queue)`. The draw then depends only on the queue's
*content*.

⚠ **Do not content-address the seed from the fragment's own fields.** Identical fragments would then draw
identically — 100 copies of one fragment would go 100/0 instead of 60/40. The queue index is what breaks
that tie correctly.

### 5.4 ⛔ Do NOT iterate — owner ruling, and it is a safety property

The confident set is biased **short**, because a longer gap admits more hypotheses. A short-biased
fragment-length prior prefers the shorter `L`, which is the **more-spliced** path; those fragments then
re-enter the fit shorter still. **That loop can run away.**

⭐ **Break it by construction: fit once, score once, drain once, stop.** Then *measure* whether a second
iteration would move anything. If it would, that measurement tells you what damping and what stopping
rule are needed. ⛔ Designing the loop before knowing it is needed is how a tuned constant gets in.

### 5.5 Where the code lives

**Selection in Python, drain in C++.** The scores need the fragment-length model and the abundances,
which live Python-side after calibration; the drain must re-enter `deposit()` so there is one tally path.
The reference gains the same selection step, so byte-identity stays cheap.

### 5.6 ⚠ Calibration stops being bit-identical — say it out loud

`CLAUDE.md` currently states that calibration is bit-identical run to run and only the EM samples. After
S3 that is **false**. With a fixed seed it is *reproducible*, which is enough — but every A/B note that
assumes determinism must be updated, and the honest check is the one already prescribed for the EM: **run
two seeds and report the spread.**

---

## §6 B — JUNCTION OPPORTUNITY, ON THE DRAINED TALLY

The formula is **derived and proven** (`JUNCTION_OPPORTUNITY.md`; 48,648 exhaustive configurations
against an independent oracle). For a transcript with exon lengths `e_i` and total `L`:

```
    A_j(w) = (L − w + 1)₊ − Σ_i (e_i − w + 1)₊
```

— because a fragment crosses **no** junction exactly when it fits inside a single exon, and the exons are
disjoint, so the complement partitions with no inclusion–exclusion.

The library aggregate is `Σ_t θ_t · A_j(w, t)`, and **θ is transcript abundance** — the estimand. Measured
sensitivity: across an adversarial sweep the corrected mean spans 214–301 bp, an **86.9 bp ambiguity
against a 16.6 bp tilt**. So the correction is *not* unconditionally safe and must use a θ that is
actually estimated, not assumed. ⭐ That θ is the same one §5.2's prior needs.

⚠ **Re-measure B's target first** (§2.4). Every number quoted for B predates the drain.

---

## §7 OPEN DECISIONS

| | | |
|---|---|---|
| **D-A** | Does `SPLICE_IMPLICIT` survive as a reported class once such fragments are deferred? | ⭐ **Recommend keep, purely descriptive** — it is the scanner's census of what it saw, and C2 put QC where it is generated |
| **D-B** | May a *non-annotated* unspliced path be a hypothesis for a certified-RNA fragment? | ⭐ **Recommend no** — that is the unannotated-junction problem, deferred by owner ruling |
| **D-C** | Which θ estimator feeds §5.2 and §6? | ⭐ **Recommend the EM's own**, with a uniform first pass. It is the only option whose error is bounded by something measurable |
| **D-D** | Where does the `global_seed` live, and is it exposed on the CLI? | ✅ **ANSWERED, and it needs no new machinery** (`LEDGER.md` S2.1). The rule is **order by identity, then one stream** — never hash the content, which ties on exactly the duplicates it would harm. S1 already gives the deferred queue a canonical order, so `(global_seed, index in that queue)` is well defined. ⚠ The old worry that `EMConfig.seed` was a broken pattern to copy is **retired**: the seed was always fine; the ORDER it was consumed in was not, and that is now fixed for the whole pipeline |

### ⛔ What S1 settled, so it is no longer open

| | |
|---|---|
| **D1** (was `TODO.md` rank 0a) | **deleted.** The RNA pool is keyed on determinacy, not provenance — a purity filter on a length pool is a length filter |
| **C1** (`SPEC_GAP_PATHS.md` §3) | **implemented and gated.** A transcript whose intron overlaps an aligned block's interior is not a candidate; `test_h2_a_CONTRADICTED_transcript_is_not_a_candidate` |
| **C3** (`SPEC_GAP_PATHS.md` §8) | **measured, two-sided, on a real scan.** The `max_fragment_length` prefilter does change classification, and both arms are now fixtures |
| **D3 / C2.7** (was rank 0b) | **solved** by grouping candidates on the whole-fragment path. ⛔ Its residual cannot be re-measured until the drain |

⚠ **Deferred by owner ruling:** unannotated splice junctions. They are a real source of long fragments in
real data, are not enumerable from the annotation, and **the simulator cannot measure them at all**. Any
residual tail on real cfRNA that survives S3 is the first place to look, and the pilot will be silent.

---

## §8 THE TWO SESSIONS, AND THE CONVERGENCE

| session | did | state |
|---|---|---|
| **the first** | C2 (landed, committed); the junction-opportunity derivation and proof; found the uncut-intron bug and specced it; measured the deferred population; designed the second pass (§5) | ⭐ retired |
| **the parallel one** | C2.6; `SPEC_GAP_PATHS.md`; the hypothesis machinery in the reference and partly in C++ | ⭐ retired |
| **the convergence one** | ✅ **S1 and S2** — the migration finished, the side buffer landed and persisted, 9 perturbations, two real holes closed | |

✅ **The tree is green again at the baseline: 21 failed / 1863 passed, all 21 goldens.**

### Reading order for the next session

1. **this file** — the order of work and the gates; §2.2 is the state, §5 is what to build
2. `LEDGER.md`, the **S1 + S2** entry — what landed, what it found, and what it deliberately did not measure
3. `SPEC_GAP_PATHS.md` §0–§4 — the enumeration rule and the interface, in detail
4. `JUNCTION_OPPORTUNITY.md` §1, §4 — the proven formula, and the evidence that found the bug

⛔ **Start with D-D (§7), then S3.** And re-measure B's target *after* the drain — every
junction-opportunity number in the docs predates it.
