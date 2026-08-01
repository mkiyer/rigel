# SPEC — the mate gap is a PATH problem, and the accumulator arbitrates it

    Status: SPECIFICATION, 2026-08-01. Partly implemented — see `PLAN_TWO_PASS.md` §2.2 for the
            exact resume point (reference done, C++ partial, call sites unmigrated).
    ⭐ ENTRY POINT IS NOW `docs/PLAN_TWO_PASS.md` — it carries the order of work, the gates, and the
      second-pass design. This file remains the detailed enumeration rule and interface (§0–§4).
    Supersedes: `SPEC_GAP_INTRONS.md` (C2.6, landed) wherever they differ — that spec cut ONE intron
                per gap and decided ambiguity in the adapter. Both are corrected here.
    Owner rulings, 2026-08-01, recorded verbatim in §8.
    Prior art it completes: `ACCUMULATOR_DESIGN.md` §9 (the deferred queue) and §9.1 (distinct intron SETS).

---

## §0 ⭐ THE RULE, IN ONE PAGE

A fragment's unsequenced mate gap may hold **no intron, one intron, or several** — and *which* is not
observable in the first pass. It is a **likelihood** question, and the likelihood needs a fragment-length
distribution that does not exist until the first pass is over.

> **A fragment arrives at the accumulator with its HYPOTHESIS SET, not with one path and a flag.**
> Each hypothesis is a set of introns cut from the fragment's genomic extent. **The empty set is the
> genomic hypothesis.** The accumulator filters hypotheses it can rule out, then:
>
> * **exactly one survives → DEPOSIT it.** The path is determined; nothing is in doubt.
> * **two or more survive → SIDE BUFFER.** There is no answer to deposit, and picking one is picking a
>   fragment length at random.

⭐ **That single rule generates the whole of §8's decision tree.** Nothing below is a special case — §2
shows each branch falling out of it, including the `max_fragment_length` rule, which turns out not to be
a separate rule at all.

⛔ **The accumulator decides. The adapter only supplies.** Today `bam_scanner.cpp` computes a bool
`path_ambiguous` and the accumulator obeys it — and `accumulator.h:210` says why: *"the accumulator
cannot decide it — only the caller has the candidate list."* It cannot decide it **because the caller
collapses the answer before handing it over**. Give it the set and the decision lives exactly where the
outcome is already reported (`DepositOutcome.DEFERRED`). ⭐ Owner ruling: the rule for what gets
deferred **will change as this work proceeds**, so it must live in one place, behind one interface.

---

## §1 WHY THE GAP IS A PATH PROBLEM — the owner's example, kept exactly

```
    TA  exons (1000,2000)      (3000,3050)      (4000,5000)      introns (2000,3000) (3050,4000)
    TB  exons (1000,2000)                       (4000,5000)      intron  (2000,4000)

    fragment    blocks [1800,1950)                    [4050,4200)
                       ==========|~~~~ unsequenced ~~~~|=========
                                        gap [1950,4050)
```

| hypothesis | introns cut | L |
|---|---|---|
| **TA** | (2000,3000) + (3050,4000) = 1950 bp | 2400 − 1950 = **450** |
| **TB** | (2000,4000) = 2000 bp | 2400 − 2000 = **400** |

Both are compatible with every sequenced base. ⭐ **TA's path crosses an exon — (3000,3050) — that no
read ever touched.** So a gap hypothesis is not "an intron"; it is a **path through the annotation**, and
a path may contain several introns and several exons.

⛔ **This is what `SPEC_GAP_INTRONS.md` got wrong twice.** It searched each gap for the *first* matching
intron of the *first* matching transcript, so TA reads as "(2000,3000)" — a path **no molecule has** — and
the unanimity test then compares that against TB's "(2000,4000)" per gap rather than comparing SETS.
`ACCUMULATOR_DESIGN.md` §9.1 already required the set test; the implementation approximated it, and the
approximation is exact only while a transcript implies **≤ 1 intron per gap**. C2.6's D3 residual is that
approximation failing, measured at **98.5 %** of the remaining fragment-length tail.

---

## §2 ⭐ ONE RULE, AND THE DECISION TREE FALLS OUT OF IT

The hypothesis set for a fragment is:

```
    H = { intron-set(t, gaps)  :  t ∈ compatible transcripts }        ∪  { ∅ if the molecule may be genomic }
```

and `∅` — cut nothing — is **the genomic hypothesis**: an unspliced molecule whose gap is real DNA.

⭐ **The nascent shadow is not a candidate and does not need to be.** Owner ruling: *"The nascent shadow
is unspliced and gDNA is unspliced… Nascent RNA and gDNA always compete within the boundaries of a
transcript span. If there is an intron gap between two reads, we are by definition within a transcript
span."* So the shadow row **is** `∅`, already present. ⚠ This retro-justifies the `~is_synthetic` filter
with a better reason than `ACCUMULATOR_DESIGN.md` §9.1's — that filter was defended as "not conflating
path with component"; the truth is simpler, the shadow is *redundant*.

⭐ **And it explains the old measurement.** §9.1 records that counting the shadow *"deferred 100 % of
implicit fragments on all three real cfRNA libraries and deposited none"*, and treated that as proof the
test was unsatisfiable. Under this spec **100 % deferral is the correct answer** for that population —
they genuinely have two explanations. It only looked like a failure because there was nowhere to put them.

### When is `∅` available?

| the fragment | is `∅` a hypothesis? | why |
|---|---|---|
| ≥ 1 **annotated** CIGAR-N junction in the reads | ⭐ **only if some candidate transcript is unspliced across the gap** (a retained-intron isoform) | gDNA cannot be spliced, so the molecule is **certified RNA** and the genomic explanation is dead. `FragmentPool.RNA_SPLICED`'s purity argument, reused |
| no annotated junction | ✅ **always** | the molecule may be gDNA or nascent, and then the gap is real template |

### The `max_fragment_length` filter, and ⭐ why the span rule is not a separate rule

Short-read chemistry does not sequence molecules much beyond 1 kb; `--max-fragment-length` (default
**1000**) already states that, and the accumulator already rejects `L` above it as `TOO_LONG`. Apply it to
the hypotheses:

> **Drop every hypothesis whose `L` exceeds `max_fragment_length` — unless that empties the set, in which
> case keep the survivors and let the accumulator reject them as `TOO_LONG` as it does today.**

⭐ **The owner's rule "if the genomic span > 1 kb, assume it is RNA" is this filter applied to `∅`.**
`∅`'s `L` *is* the genomic span, so a span above the limit deletes the genomic hypothesis automatically.
There is no second rule and no new constant.

⚠ **It is not purely a cost gate**, and this must be said plainly because the owner called it one. It
changes *classification*: hypotheses at `L = 400` and `L = 1200` are "determined" with the filter and
"ambiguous" without it. That is defensible — it is the second pass's likelihood applied early, at a point
where the FL distribution has no mass — but its **false-positive rate is a measurement, not an
assumption**, and §9's G1 gate takes it.

### The tree, re-derived

| fragment | hypotheses | outcome |
|---|---|---|
| no gap, contiguous span | `{∅}` | ✅ **deposit** — unspliced |
| gap, no candidate intron inside it | `{∅}` | ✅ **deposit** — confirmed unspliced |
| certified RNA, all candidates agree, no retained-intron isoform | `{S}` | ✅ **deposit**, cut `S`, RNA pool |
| certified RNA, candidates disagree | `{S₁, S₂, …}` | ⛔ **deferred queue** |
| not certified, gap intron, span ≤ 1 kb | `{∅, S₁, …}` | ⛔ **deferred queue, always** |
| not certified, gap intron, span > 1 kb | `{S₁, …}` (`∅` filtered) | ✅ deposit if one survives, else deferred queue |

---

## §3 THE ENUMERATION

```
    enumerate_gap_paths(exons, observed_introns, candidate_t, scratch) -> GapPaths
```

1. **Gaps** = holes between merged aligned blocks, ⛔ **minus every gap that exactly matches an observed
   CIGAR-N intron**. This is C2.6 §2.1 and it survives unchanged: exact `(start, end)` equality, never
   overlap, or a *different* annotated intron within the ±K anchor tolerance answers for one the CIGAR
   already stated and `L` comes out too **short** (measured: 500 bp against a 502 bp molecule).
2. **Per candidate transcript** `t`, collect **every** intron of `t` lying inside **every** gap (±K). That
   ordered list is `t`'s path. ⭐ There is no choice within a transcript — a transcript determines its own
   path — so the enumeration is over transcripts, and this is what makes the set finite and small.
3. **Group by path.** Distinct paths are distinct hypotheses; the transcripts sharing one are its
   supporters (and carry the abundance prior the second pass needs).
4. **Add `∅`** per §2, if not already produced by a retained-intron candidate.
5. **Filter** by `max_fragment_length`, unless that empties the set.

⚠ **The path is the whole fragment's, never one gap's.** §9.1's warning is live: a per-gap union can take
gap A's intron from T1 and gap B's from T2 and emit a path **no single molecule has**. Group by the
full-fragment path and that is structurally impossible rather than merely avoided.

### ⚠ CONCERN C1 — `t_inds` is overlap-compatible, not structure-compatible

The owner's rule is *"any transcript with a compatible combination of introns + exons is a candidate."*
`cr.t_inds` is not that set. It comes from `merge_sets`, which **falls back to `MC_UNION` when the
intersection is empty**, and the nRNA shadow is injected into every block's exon set. So a transcript the
reads *contradict* — one with an intron inside an aligned block — can be in `t_inds` and would contribute
a spurious hypothesis.

⭐ **There is an existing quantity that is exactly right when the candidate is compatible.**
`compute_frag_lengths_aligned` already computes, per candidate, `|tx_pos(gend) − tx_pos(gstart)|` — every
intron of `t` between the endpoints, cut. For a structurally compatible `t` every such intron is either
observed, or in a gap, or inside an aligned block (which is the contradiction), so **B ≡ the per-hypothesis
`L`** exactly when there is no contradiction. ⚠ That makes it a *test*, not a shortcut: a candidate whose
definition-B length disagrees with its enumerated path length **is incompatible and must be dropped**. §9
G1 gates on it.

---

## §4 THE INTERFACE — what the adapter supplies

```
    FragmentPaths {
        ref, start, end          // the extent on ONE reference, mate gap included
        align_strand             // the column; the accumulator still rejects an undefined one
        observed_introns         // CIGAR-N: cut under EVERY hypothesis
        hypotheses[]             // { implied_introns, sj_strand, supporting_t_inds[] }
    }
```

⭐ **`observed_introns` sits outside the hypothesis list because it is not in doubt.** Every hypothesis
cuts it. That keeps each hypothesis to what actually varies, and it means a fragment with no gap carries
exactly one empty hypothesis — the degenerate case is the general case, not a branch.

```
    Accumulator::deposit(FragmentPaths, scratch) -> DepositOutcome
```

`DepositOutcome` is unchanged; `DEFERRED` stops being an instruction and becomes a **verdict**.

---

## §5 ⭐ WHAT THIS DELETES — the cleanliness argument, stated as a diff

The owner asked for code cleaner than what is there. It is a net **deletion**:

| gone | why |
|---|---|
| `FragmentPath.sj_implicit` | ⭐ **and D1 with it.** It means "part of this `L` was inferred". A *directly deposited* fragment now has exactly one possible path, so its `L` is not in doubt at all — **determinacy replaces provenance**. ⭐ C2.6 measured that trade: the pool reads **+0.67 % / +2.40 %** against truth on determinacy and **−9.58 % / −22.46 %** on provenance. A deletion that improves a number |
| `FragmentPath.path_ambiguous` | an input that was already an outcome |
| `RawResolveResult.implicit_ambiguous` | the verdict is the accumulator's now |
| the per-gap `seen_none` / `seen_intron` unanimity scan | ⭐ ~40 lines of subtle bookkeeping, and three comment blocks defending it, replaced by "group candidates by path" — which is what §9.1 asked for in the first place |
| the `is_real` / `~is_synthetic` filter *as a special case* | the shadow is `∅`, which the enumeration produces anyway |
| `SPEC_GAP_INTRONS.md`'s D1, D2, D3 | D1 deleted, D2 moot (`sj_strand` is per hypothesis), D3 **solved** rather than deferred |

⚠ **`splice_type` does not move, again.** It classifies what was **observed** and feeds scoring, the
deferred queue, strand training and the report's census. ⭐ **OPEN DECISION D-A (§8):** `SPLICE_IMPLICIT` now
describes a fragment that is *deferred*, not deposited — is the class still meaningful, or does it become
purely descriptive?

---

## §6 THE SIDE BUFFER

`ACCUMULATOR_DESIGN.md` §9 specifies it as `(candidate object ids, channel, L)` and it **does not exist**.

⭐ **Store the fragment, not its consequences.** Object ids are large, derived, and would have to be kept
consistent with the partition; the fragment is small and replays exactly. One record:

```
    ref, start, end, align_strand, observed_introns[], hypotheses[]
```

⭐⭐ **The drain re-enters `Accumulator::deposit` with the chosen hypothesis.** There is no second deposit
path, no duplicated crossing/containment logic, and byte-identity with
`tests/native/_accumulator_reference.py` is preserved for free — the specification gains a *selection*
step, not a second tally.

⚠ Per-worker deferred queues merge like the tally does, and the merge must be **deterministic at any worker
count** — `test_accumulator_worker_determinism.py` extends to cover it. ⛔ `payload_schema_digest`
**moves**; all 8 pilot caches rebuild (~46 s).

---

## §7 THE SECOND PASS

1. **Fit** the fragment-length distribution on **directly deposited** fragments only — unique mappers with
   exactly one surviving hypothesis. That is the confident set, and it is what the pools already are.
2. **Score** each deferred hypothesis: `P(L | FL) × prior(path)`, where the prior needs **transcript
   abundance** for a spliced path and the **gDNA density** for `∅`.
3. **Assign** — §9's multinomial across matching paths.
4. **Drain** through `deposit()`.

### ⛔ CONCERN C2 — the bootstrap is length-biased, and it is positive feedback

The confident set is **biased short**, because a longer gap admits more hypotheses. C2.6 measured that
bias directly on the nearest available proxy: excluding the ambiguous long fragments moves the pool from
**+0.67 % / +2.40 %** to **−9.58 % / −22.46 %** against truth.

A short-biased FL prior prefers the **shorter** `L`, which is the **more-spliced** path; those fragments
then re-enter the fit even shorter. ⛔ **That loop can run away and needs a stopping rule and a measured
convergence check** — the same requirement `JUNCTION_OPPORTUNITY.md` §7(c) has for θ.

⭐ **And that is the argument for building it once.** C3 needs post-capture transcript abundance; this
needs transcript abundance; both want the EM's own. **One two-pass structure, shared.**

---

## §8 ⛔ DECISIONS AND CONCERNS

| | | |
|---|---|---|
| **D-A** | Does `SPLICE_IMPLICIT` survive as a class once such fragments are deferred rather than deposited? | ⭐ **Recommend: keep, purely descriptive.** It is the scanner's census of what it saw, and C2's ruling put QC where it is generated |
| **D-B** | `∅` for a certified-RNA fragment requires a retained-intron candidate. Should a *non-annotated* unspliced path also be allowed? | ⭐ **Recommend: no.** The molecule is RNA; an unannotated retained intron is the unannotated-junction problem, deferred by the owner |
| **C1** | `t_inds` is overlap-compatible, not structure-compatible — §3 | needs the compatibility predicate; **G1 gates on it** |
| **C2** | the FL bootstrap is short-biased and self-reinforcing — §7 | needs a convergence rule and a measured check |
| **C3** | the `max_fragment_length` prefilter changes classification, not just cost — §2 | measure its FP rate; **G1 gates on it** |
| **C4** | between G1 and G2 the tally is **thinner** — the ambiguous mass is retained in the deferred queue but not yet deposited | ⛔ **G1 must NOT be judged by a calibration A/B.** Its gate is conservation, not accuracy |

---

## §9 ROADMAP

| | step | gate |
|---|---|---|
| **G1** | enumeration + the `FragmentPaths` interface + accumulator arbitration + deferred-queue **storage**, landed **together** | ⭐ **conservation: `deposited + deferred + dropped_* == offered`, exactly.** Plus: the owner's §1 example reproduced to the base pair (450 and 400, both deferred); a compatibility-predicate test (C1); the FP rate of the `max_fragment_length` prefilter (C3); byte-identity to the reference; determinism at 1/2/4/8 workers |
| **G2** | the second pass: fit on the confident set, score, assign, **drain through `deposit()`** | the anchor's mass above the library's true ceiling — **713 bp** on the pilot, read from `truth_fragment_lengths.tsv` — reaches **0**, and `dropped_too_long` collapses. ⛔ C2.6 left this at **0.00137** and **38,309**; M3 showed the mechanism is exactly what G1 enumerates |
| **G3** | the shared abundance prior + the convergence check | ⭐ the gDNA control must stay bit-identical throughout (it has no introns to miss), and the FL fit must converge with a **measured** stopping rule, not a chosen one |

⚠ **G1 and G2 both move the goldens.** ⛔ Still regenerate **once**, at the end, **twice, and diff**.

---

## §10 DEFERRED, BY OWNER RULING

**Unannotated splice junctions** are a real source of long fragments in real data and absent from the
simulator. They are not enumerable from the annotation, so they are out of scope here — ⚠ and that means
the pilot **cannot** measure this population at all. Any residual tail on real cfRNA that survives G2 is
the first place to look, and the pilot will be silent about it.
