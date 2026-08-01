# The junction pool's opportunity function — the derivation, and what measuring it found

    Status: DERIVATION COMPLETE, 2026-08-01. ⛔ NO CODE WRITTEN.
    ⭐ ENTRY POINT IS NOW `docs/PLAN_TWO_PASS.md`. This work is step S4 there, and it rides on the
      SAME two-pass structure as the gap-path work — §0 of that file is why.
    ⚠ §3's numbers PREDATE the side buffer and must be re-measured after the drain — PLAN_TWO_PASS §2.4.
    ⭐ §8 D1 IS ANSWERED: the uncut intron jumped the queue and LANDED as C2.6 (`LEDGER.md`).
    ⛔ **§3's MEASUREMENTS ARE NOW STALE AND MUST BE RE-RUN.** Every number in §3 was taken against the
      contaminated anchor (sd +27.0 %, now +1.98 %) and against an `RNA_SPLICED` pool that has since
      moved from +8.00 % / +67.35 % to −9.58 % / −22.46 % against truth. ⚠ In particular §3.2's θ
      control — "corrected, TRUE θ → +0.3 % mean / +59.8 % sd" — is scored on both of those, so
      **C3's real target is not knowable from this file any more.** §4's diagnosis stands; §1's
      derivation and its 48,648-configuration proof are untouched.
    Asked for: `ACCUMULATOR_DESIGN.md` §8.1(b) gives the opportunity for two of the three pools and
    leaves the junction pool as "the transcript-level count" — and §"what the accumulator will not
    decide" defers it as *"a fragment-length-model question, downstream of the tally."* So C3's
    central correction had no formula. This supplies it.

---

## §0 ⭐ THE ANSWER, AND THE SURPRISE

**The formula exists, it is exact, and it is proven.** For a transcript with exon lengths
`e_1 … e_K` (transcript space) and total length `L = Σ e_i`, the number of start positions at which a
length-`w` fragment crosses **at least one** junction is

```
    A_j(w) = (L − w + 1)₊  −  Σ_i (e_i − w + 1)₊
```

⭐ **48,648 exhaustive configurations agree with an independent enumerating oracle.**

**But measuring it overturned the premise C3 was scoped on.** The audit records the zero-gDNA residual
as **+7.7 % mean / +32.0 % sd** and calls the whole of it "the junction-opportunity tilt, and it is
C3's target". Correcting the pool with the **true** abundances — which the simulated pilot supplies —
against the **true** fragment-length distribution:

| vs simulated TRUTH | mean | sd |
|---|---|---|
| RNA_SPLICED pool, raw | +8.0 % | +67.3 % |
| ⭐ **corrected, TRUE θ** | ⭐ **+0.3 %** | ⛔ **+59.8 %** |
| the anchor `deposited_lengths` | +0.4 % | ⛔ **+27.0 %** |

> ⭐ **C3 fixes the mean almost exactly and closes only ~19 % of the sd gap.** The mean error is
> junction opportunity and the formula removes it. **The sd error is not, and no opportunity function
> can remove it** — it is present in the *anchor* too, and it is a **defect upstream of all of this**.

⛔ **The sd excess is fragments the library never contained.** The simulator's realized support is
**50–713 bp**. The accumulator reports **1.35 % of mass ≥ 600 bp**, **0.97 % ≥ 700 bp**, and rejects a
further **6.1 % as longer than 1000 bp**. Truncating at 600 bp removes **95.5 %** of the sd excess.

⭐ **And the cause is located, to one line.** §4.

---

## §1 THE DERIVATION

### 1.1 What "opportunity" means here

`NODE_DENSITY_DERIVATION.md` §4.1: a population's opportunity `A(w)` is **the number of integer start
positions producing that event**, and because starts are Poisson,

```
    E[count(w)]  ∝  f(w) · A(w)          ⟹          f(w) = count(w) / A(w)
```

That is the whole content of §8.1(b)'s "divide each pooled histogram by its own opportunity". Two of
three pools already have theirs: `(ℓ − w + 1)₊` contained, `placements(w)` crossing.

### 1.2 The junction pool

RNA lives on transcripts, so the start positions are **transcript-space**. A fragment
`[s, s+w)` on transcript `t` enters `RNA_SPLICED` iff it crosses ≥ 1 annotated junction.

⭐ **Work with the complement, and it decomposes with no inclusion–exclusion.** A window crosses **no**
junction **iff it lies wholly inside a single exon** — and the exons are disjoint, so the placements
partition by *which* exon contains the window:

```
    placements crossing nothing  =  Σ_i (e_i − w + 1)₊
    total placements             =  (L − w + 1)₊
    ⟹   A_j(w)  =  (L − w + 1)₊ − Σ_i (e_i − w + 1)₊
```

⚠ Attempting this by inclusion–exclusion **over junctions** is where §8.1(b)'s "not `w−1`" warning
leads — the union of "crosses junction `j`" events is messy and gets worse with every exon. The
complement is a partition, so it is exact in one line.

### 1.3 It is proven, not asserted

`scratchpad/junction_opportunity_proof.py`. The oracle **enumerates every start position** and tests
`s < c < s+w` per junction — it shares no code with the formula (`CARRY_FORWARD.md` §3 trap 1).

| coverage | |
|---|---|
| exhaustive, 1–4 exons × lengths 1–7, every `w` up to `L+2` | **48,552** configurations |
| realistic scales (1 bp exons, 2 kb exons, 20 × 50 bp, `w` at every boundary) | 96 more |

And the properties the derivation claims, **checked rather than stated**:

| | property | why it matters |
|---|---|---|
| **P1** | a single-exon transcript has `A_j ≡ 0` | an unspliced transcript can never populate this pool — the formula knows it |
| **P2** | `w = 1 ⇒ A_j = 0` | a 1 bp fragment crosses no 0-bp line |
| **P3** | `A_j` **rises** with `w` up to the longest exon | ⭐ **this IS the tilt** — the measured +0.70 log-ratio-vs-length correlation |
| **P4** | `w >` longest exon ⇒ `A_j = (L−w+1)`, i.e. every placement crosses | the tilt **saturates**; it is not unbounded |
| **P5** | two transcripts of **equal length** but different exon structure have different `A_j(w)` | ⛔ **the aggregate depends on WHICH transcripts are expressed** |

---

## §2 ⛔ WHY THIS POOL IS THE HARD ONE — P5, STATED PROPERLY

The library-level opportunity is a **θ-weighted sum**:

```
    A_pool(w)  =  Σ_t θ_t · A_j(w, t)              θ_t = molar abundance (copies) of transcript t
```

⚠ **θ is the quantity the tool exists to estimate.** And by P5 it does not cancel: a transcript of many
short exons saturates early (steep tilt), one of few long exons saturates late (shallow tilt), so the
*shape in w* of the aggregate depends on the expression profile.

⭐ **This is exactly why the other two pools are easy and this one was deferred.**

| pool | weights | known? |
|---|---|---|
| gDNA contained | gDNA density along the genome — **modelled as uniform** | ✅ annotation only |
| crossing | same | ✅ annotation only |
| **junction (RNA)** | **transcript abundance** | ⛔ **the estimand** |

`ACCUMULATOR_DESIGN.md`'s deferral — *"a fragment-length-model question, downstream of the tally"* — is
therefore correct and not an oversight. The tally cannot know θ.

---

## §3 THE MEASUREMENTS

`scratchpad/theta_sensitivity.py`, `scratchpad/theta_truth_control.py`, on the chr22 + ERCC index
(5,440 non-synthetic transcripts, 640 single-exon and so contributing nothing to this pool) and the
8 cached pilot conditions.

### 3.1 How much does the answer depend on θ? — a lot, if θ is unconstrained

`gdna_none ss0.99 capture_off`, corrected pool vs the anchor:

| θ regime | corrected mean | corrected sd |
|---|---|---|
| uniform (annotation only) | +4.0 % | +36.8 % |
| length-weighted | +0.4 % | +27.1 % |
| log-normal expression, sd 2 | +4.1 % | +37.5 % |
| log-normal expression, sd 4 | +19.1 % | +64.0 % |
| ⛔ **adversarial: short-exon quartile only** | **+38.0 %** | **+106.7 %** |
| ⛔ **adversarial: long-exon quartile only** | **−1.8 %** | **+20.7 %** |

⛔ **Across the full sweep the mean spans 214–301 bp — an 86.9 bp ambiguity against a 16.6 bp tilt.**
**5.2× larger than the thing being corrected.** A correction applied with a badly wrong θ is worse than
no correction at all.

⚠ The adversarial rows are not realistic libraries; they are there to show the correction is **not
unconditionally safe**, which a sweep over only plausible θ would have hidden.

### 3.2 ⭐ THE CONTROL: correct with the TRUE θ, score against the TRUE distribution

Both are on disk for the pilot (`truth_abundances.tsv` → `mrna_abundance`, the molar abundance;
`truth_fragment_lengths.tsv` → the realized mRNA length distribution).

⚠ θ is the **molar** abundance, not the observed fragment count: `A_j` already counts start positions,
and observed fragments ∼ copies × effective length. Using the count applies the length weighting twice.

| `gdna_none ss0.99 capture_off`, vs TRUTH (mean 217.13, sd 87.41) | mean | sd |
|---|---|---|
| anchor `deposited_lengths` | +0.4 % | +27.0 % |
| RNA_SPLICED pool, raw | +8.0 % | +67.3 % |
| ⭐ **corrected, TRUE θ** | **+0.3 %** | **+59.8 %** |
| corrected, uniform θ | +4.4 % | +73.6 % |
| corrected, length-weighted θ | +0.8 % | +61.4 % |

⭐ **The derivation is right.** With θ known, the mean error collapses **+8.0 % → +0.3 %**. §8.1(b)'s
claim that this feeds `μ_g − μ_r` and that a 10 % error costs 0.010–0.026 of composition makes the
8 % mean error worth **≈ 0.008–0.021** — real, and C3 removes it.

⛔ **And the sd is untouched: +67.3 % → +59.8 %.** Against the anchor that is 31.8 % → 25.9 %, i.e.
**C3 done perfectly closes ~19 % of the sd gap the audit assigned to it.**

### 3.3 The capture question, answered

C2.1 flagged that the residual is capture-dependent — **+7.7 % / +32.0 %** off, **+11.0 % / +43.4 %**
on — and asked whether §8.1(b) needs a capture-aware term. **It does not.**

| vs TRUTH | capture_off | capture_on |
|---|---|---|
| pool raw, mean | +8.0 % | +11.5 % |
| corrected TRUE θ, mean | +0.3 % | **+1.6 %** |
| transcripts with non-zero abundance | 2,592 | **1,845** |

⭐ Capture changes **which transcripts are expressed**, and that is a change in θ, not a new physical
effect. Correcting with the *post-capture* θ removes nearly all of it. ⚠ The requirement this creates is
that θ must be the **post-capture** abundance — the thing the library actually contains, not the input
molarity.

---

## §4 ⛔ THE FINDING: THE SD GAP IS AN UNCUT INTRON, AND THE GATE IS ONE LINE

### 4.1 The anchor is measuring lengths the library does not contain

`scratchpad/sd_excess_diagnostic.py`. The simulator's realized mRNA support is **50–713 bp**
(5,000,000 fragments). The accumulator's anchor:

| ≥ bp | truth | anchor | ratio |
|---|---|---|---|
| 500 | 0.0014 | 0.0188 | 13× |
| 600 | 0.000029 | 0.0135 | **463×** |
| 700 | **0** | 0.0097 | ⛔ **∞** |
| 800 | **0** | 0.0060 | ⛔ **∞** |

…and `qc.dropped_too_long = 280,558`, **6.1 % of deposits**, rejected for exceeding 1000 bp on a
library whose longest molecule is 713 bp.

⭐ **Truncating the anchor at 600 bp removes 95.5 % of its sd excess.** The +27 % is a tail, not a shape
error.

### 4.2 ⭐ THE CONTROL THAT NAMES THE MECHANISM

`gdna100 ss0.99 capture_off`. gDNA fragments have **no introns to miss**; RNA fragments do.

| series | ≥600 | ≥700 | ≥800 |
|---|---|---|---|
| **DNA_INTERGENIC pool** (pure gDNA) | 0.00023 | 0.00001 | 0.00000 |
| ⭐ **TRUTH gdna** | **0.00024** | **0.00001** | **0.00000** |
| **RNA_SPLICED pool** (pure RNA) | **0.04182** | **0.03002** | **0.01853** |

> ⭐ **`L` is exact to five decimal places when the fragment has no introns, and badly wrong when it
> has them.** That excludes every generic explanation — a span bug, chimeras, mate-gap handling — all of
> which would hit gDNA identically. ⚠ And gDNA truth runs to 785 bp, so the RNA pool's 1.85 % above
> 800 bp is not gDNA leaking in either.

⚠ **C0's proof still stands.** `L` is correct *given its inputs* — 333,684 configurations say so. The
inputs are wrong: the **intron list is incomplete**.

### 4.3 The line

`src/rigel/native/resolve_context.h:1442`:

```cpp
if (cr.splice_type == SPLICE_UNSPLICED &&          // ⛔ HERE
    cr.chimera_type == CHIMERA_NONE &&
    collect_implicit_splice_introns(...)) {
    cr.splice_type = SPLICE_IMPLICIT;
}
```

⛔ **Implicit-splice detection runs ONLY on fragments the resolver already called UNSPLICED.** A
fragment with an observed CIGAR-N splice *and* a second intron sitting in its unsequenced mate gap
never has the gap intron detected, so that intron is never cut and `L` over-counts by its full length.

**It predicts precisely the measured pattern:**

| population | detection runs? | tail ≥600 |
|---|---|---|
| gDNA, no introns | n/a | 0.0002 — exact |
| anchor, mostly unspliced | ✅ yes | 0.0066 |
| **RNA_SPLICED, always has an observed splice** | ⛔ **never** | **0.0418** |

The pool with an observed splice — the one where the gate *guarantees* detection is skipped — is
**6× worse than the anchor**. Such fragments must span ≥ 2 introns, so they are long by construction:
exactly the tail.

⚠ **Stated as the leading hypothesis with its evidence, not as proven.** The falsification is direct
and cheap: count fragments that have an observed splice **and** a candidate transcript's intron wholly
inside a mate gap, and check that their `L` distribution is the tail. That is the first step of any fix,
and it is not written here because this document writes no code.

### 4.4 ⚠ A correction to the record: C1 read the support ceiling backwards

C1 recorded the anchor's support moving from the scanner's **713** to the accumulator's **1000** as
*"the support ceiling mismatch closed entirely"* — a headline win.

⛔ **713 was right.** It is the library's true maximum, to the base pair. The scanner's definition **B**
took its length from the *transcript*, so an intron is removed whether or not it was sequenced — that
definition **cannot** produce an uncut intron. Matching at 1000 is matching `max_frag_length`, the
clamp, not the library.

⚠ **This does not undo C2.** Definition B was deleted for reasons that all still hold: it was summed
into one array with a genomic-span definition, it silently discarded 4.6 % of its own candidates, and
its population was never stated. But on this one question it was right and `L` is wrong, and the
ledger's C1 entry should not be read as evidence that `L`'s support is correct.

---

## §5 WHAT THIS MEANS FOR C3

C3 as scoped in `FRAGMENT_LENGTH_AUDIT.md` §4 targets *"the residual +7.7 % / +32 %"*. Measured, that
residual is **two unrelated defects**:

| | | fixed by |
|---|---|---|
| the **mean** error, +8.0 % vs truth | junction-opportunity tilt | ⭐ **C3**, essentially exactly, given θ |
| the **sd** error, +27 % in the anchor and +67 % in the pool | ⛔ **uncut introns** — §4 | **not C3.** A new step, upstream of it |

⛔ **The sd defect is upstream of the anchor**, so it is in every consumer C2 just wired up: the scorer,
the calibration divisors, and — by D7 — **every transcript's effective length in the EM**. It is not
confined to the length channel.

⭐ **And it inverts the priority.** C3's prize is ~0.008–0.021 of composition. The uncut-intron tail
misstates 1.4 % of fragments by hundreds of bp and throws away another 6.1 % as too long, on a library
where nothing is longer than 713 bp.

---

## §6 G3, RESTATED SO IT CAN BE MET

§5's gate reads: *"On a zero-gDNA condition, `unconditional_pmf` vs de-tilted `rna_fl_pmf`: mean and sd
agree within sampling error, and the support ceilings match."*

⛔ **As written it cannot pass after C3, and it has no threshold.** "Within sampling error" is
undefined, and inventing a tolerance would be `CARRY_FORWARD.md` §3 trap 12 — a tuned constant hiding
an upstream defect. It splits:

| | gate | target |
|---|---|---|
| **G3a** | the corrected pool's **mean** against the anchor's | ⭐ a **measured** target, not a tolerance: with true θ the residual is **+0.3 %**, so C3 must land within the θ-estimation error of that. Score the θ estimator against the pilot's known θ and quote the induced mean error |
| **G3b** | the **sd**, and the support ceiling | ⛔ **not C3's gate.** It belongs to the uncut-intron work, and its target is the true library maximum — **713 bp**, which the simulator states exactly |

⚠ Neither introduces a constant. G3a's target is measured per-condition from data already on disk;
G3b's is read from the truth file.

---

## §7 THE θ ESTIMATOR — the remaining piece of C3, and its options

Given §3.1, C3 needs a θ good enough that the induced mean error is well under the 8 % it removes.

| | option | cost | ⚠ |
|---|---|---|---|
| **a** | **uniform** over non-synthetic transcripts | free, annotation only | +4.4 % residual — halves the error only |
| **b** | **length-weighted** | free, annotation only | +0.8 % measured here, ⛔ **but that is one condition and it has no mechanism behind it** — it may be coincidence of this simulator |
| **c** ⭐ | **the EM's own θ, one extra pass** — calibrate with an uncorrected pmf, run the EM, recompute `A_pool` from its abundances, re-calibrate | one extra calibration (~2.6–6.4 s/condition cached) | the honest answer; needs a convergence check |
| **d** | **observed junction counts** (`sj_count`) as weights | free, already in the payload | RNA-pure by construction, but junction-level not transcript-level, and it is the same data the pool is drawn from — circularity needs proving harmless |

⭐ **Recommendation (c)**, with (a) as the first-pass seed. It is the only option whose error is
*bounded by something measurable* rather than by a hope about the annotation.

⚠ **(b) must not be adopted on the strength of the +0.8 %.** One condition, one simulator, no
mechanism — that is the shape of a coincidence, and `CARRY_FORWARD.md` §3 trap 19's lesson applies.

---

## §8 ⛔ THE TWO OWNER DECISIONS

### D1 — ⭐ Does the uncut intron jump the queue ahead of C3?

The evidence says it should. It is a **larger defect**, it sits **upstream** of C3, it contaminates the
**anchor** that C3 is scored against, and by D7 it reaches every transcript's effective length in the
EM. Fixing C3 first means de-tilting a pool against a reference that is itself wrong by +27 % in the sd.

| | order | consequence |
|---|---|---|
| **1** ⭐ | **uncut introns first**, then C3 | C3 is then scored against a clean anchor and G3b becomes meetable. ⚠ Touches the resolver and reopens the S3 byte-identity gate; the goldens move a third time |
| **2** | **C3 first**, as planned | keeps the recorded plan; ⛔ C3's own gate G3 cannot pass, and the mean correction is scored against a contaminated reference |

### D2 — Which θ estimator (§7)?

Recommend **(c)**, the EM's own θ with a uniform first pass. It is the only option with a measurable
error bound. ⚠ It makes calibration two-pass, which is a structural change worth your explicit sign-off.

---

## §9 WHAT WAS NOT DONE

⛔ **No code.** No formula is implemented, no gate is written, nothing is wired.
The three measurement scripts live in the session scratchpad and are referenced above; they are
diagnostics, not deliverables, and none of them is in `scripts/`.
