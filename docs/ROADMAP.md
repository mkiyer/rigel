# ROADMAP — what to do next, and why, in priority order

    Status: 2026-08-03. ⭐ **Fragment length is DONE.** It was the blocker on this whole area for months.
    Everything below is ranked by what it buys, and every claim here is a MEASURED number with the entry
    in `docs/WIP.md` that produced it.

**Sibling documents.** `docs/TODO.md` is the exhaustive ranked list including the small self-contained
items; this file is the *narrative* — why the order is what it is. `docs/WIP.md` is what landed, with its
gates. `docs/SESSION_HANDOFF.md` §3 is the trap list, and it is the most useful thing in the repo to read
before designing anything.

---

## §0 ⭐ WHERE THE TOOL ACTUALLY IS — four numbers

| | measured against | |
|---|---|---|
| fragment length, the **anchor** | the simulator's own record, 8 pilot conditions | ⭐ **+0.00 % mean / +0.02 % sd** |
| the second pass, **per fragment** | 171,534 held fragments vs per-fragment truth | ⭐ **90.5 % exact**, mean error **+0.12 bp** |
| fragments longer than any molecule in the library | all 8 conditions | ⭐ **0** |
| ⛔ **the deliverable** — library gDNA fraction | 4 contaminated conditions vs truth | ⛔ mean \|error\| **0.0472** |

⚠ **The first three are upstream plumbing. The fourth is the product.** The tool exists to deconvolve a
library into gDNA and RNA; everything else is in service of that. Read the roadmap in that light.

---

## §1 ⛔ THE CENTRAL FINDING — and it is why the order below is what it is

`WIP.md` **B4**: draining the side buffer made the composition estimate **worse**, mean \|error\| on the
contaminated conditions **0.0381 → 0.0472 (+23.9 %)**. Decomposed:

* ⭐ **~55 % is the drain being RIGHT.** Truth says **99.8 %** of held fragments are RNA, so depositing them
  genuinely lowers the gDNA *fraction* — and the estimate was already **too low** on 3 of 4 conditions, so
  adding real RNA moves it further out. **The drain exposes a pre-existing gDNA under-call rather than
  creating one.**
* ⚠ **~45 % is the fragment-length POOLS moving.** `calibrate` fits from the two **pure pools**, not from
  the anchor — and the anchor is what the second pass made exact. `RNA_SPLICED` went the *other* way,
  **+2.4 % → +6.2 %** against truth.

⭐⭐ **THE MECHANISM, IN ONE SENTENCE.** The `RNA_SPLICED` pool is selected on *"used an annotated
junction"*, and that population is genuinely **+3.8 %** longer than the library — because seeing a splice
is roughly length-independent while having an unsequenced mate gap is a pure length threshold. The drain
feeds long junction-using fragments into that pool, so it makes a known, uncorrected bias **bigger**.

⛔ **So the next step is not more accuracy upstream. It is the correction for that pool.**

---

## §2 RANK 1 — C3: the junction pool's opportunity correction

`docs/accumulator/JUNCTION_OPPORTUNITY.md`. ⭐ **The formula exists and is proven**; no code is written.

For a transcript with exon lengths `e_1 … e_K` in transcript space and total `L = Σ e_i`, the number of
start positions at which a length-`w` fragment crosses **at least one** junction is

```
    A_j(w) = (L − w + 1)₊  −  Σ_i (e_i − w + 1)₊
```

⭐ Verified against an independent enumerating oracle over **48,648 exhaustive configurations**.

**What it buys.** It is the only correction for the bias §1 identifies, and it is on the path from "the
length model is accurate" to "the composition estimate is accurate".

**What it also unblocks.** `docs/calibration/SOLVER_OBSERVABLES_PLAN.md` **P2** — the fragment-length
likelihood in the per-node solve — has been **built and gated OFF since 2026-07-31, blocked on the FL
pools**, with its mechanism already proven (blind mass **100 % → 0 %**). C2 removed half that block; C3
removes the rest. That is the payoff the whole fragment-length track was scoped for.

⛔ **Before scoping it, re-measure.** Every number in that file's §3 is **stale** — taken against the
contaminated anchor and the pre-drain pool. ⭐ `WIP.md` **P4.2** supplies the replacement targets: the
observed-splice population is **+3.8 %** longer than the library, and after the drain the pool reads
**+6.2 … +8.1 %**.

⚠ **And note what C3 cannot fix.** Measured with the true abundances, the correction fixes the pool's
**mean** almost exactly and closes only ~19 % of the **sd** gap. The sd error is a different problem and
naming it as C3's target would be a mis-scope.

---

## §3 RANK 2 — turn on the per-node fragment-length likelihood

Immediately after C3, because it is already written and gated off. Mechanism proven: it takes the share of
library mass reaching the solver with **no composition evidence of its own** from **100 % → 0 %** on the
conditions where the strand channel is blind.

⚠ Its own A/B harness exists (`scripts/design/length_likelihood_ab.py`) and its docstring already names the
sharpest target: `gdna100 ss0.50 capture_on`.

---

## §4 RANK 3 — the worst row in the suite: `gdna100 ss0.50 capture_on`

Reads `f_gdna` **0.3704 against a truth of 0.5** — a −0.13 error, by far the worst, **unexplained since
S5.f and unmoved by everything since.** Unstranded *and* capture-enriched: the strand channel carries no
information at all, and capture destroys the intron signal **75×**. ⭐ That combination is the design's
hardest failure mode, and it now has a name, a number and a reproducible condition.

⚠ A related datum worth keeping in view: on an **unstranded zero-gDNA** library the tool reads
`f_gdna = 0.0854` against a truth of **0**. With no strand information composition comes from length alone
and over-calls gDNA by 8.5 pp. The stranded zero-gDNA rows read 0.0018–0.0032.

---

## §5 RANK 4 — regenerate the goldens, once, at the end

21 `test_golden_output` failures are the **standing baseline**, and they have moved **six times**: P1 (the
EM prior's units), C2 (the FL models), C2.6 (`L` itself), S1 (the hold-out), P4 (the drain wired in), P4.2
(the combination rule). C3 will move them again.

⛔ **Regenerate twice and diff.** The goldens run under the default `sample` assignment mode, so a flaky
expectation baked in now is permanent. ⚠ And regenerating is **not** validating: it records whatever the
tree does. Adjudicating the 16–52 % `gdna_em_count` drop from the index change is a separate comparison
against the suite's truth.

---

## §6 WHAT IS DELIBERATELY NOT NEXT, each with its measurement

| | why not now |
|---|---|
| **iterating the second pass** | ✅ **Answered by measurement, not by building the loop.** The anchor is already exact, the assignment is right to **+0.12 bp**, and the pool's residual is a **selection** effect — iteration cannot fix one |
| **a Poisson likelihood for the traffic term** | ⭐ Correct in principle: zero observations is `P(0 \| λ, E) = e^(−λE)`, not zero. And `rho` is *partially* zero on **62 %** of held records, so a sampling zero hard-eliminates a candidate on most of them. ⛔ But it is right **99.9 %** of the time here — **15** misassignments in 171,534 — because the hard zero is the **large-exposure limit** of that likelihood and the pilot has 5 M confident fragments. ⚠ It needs each object's **exposure**, which *is* junction opportunity. **Do C3 first and this becomes cheap** |
| **D-4, the density's weight of evidence** | The same question from the other side, deferred on the same 15-fragment measurement, wanting the same exposure |
| **a fitted strand mixture for the genomic candidate** | ⛔ **Implemented, measured and REFUTED.** The orientation discrimination is `(1−p)/p` and any *constant* for ∅ cancels out of it, so an orientation-dependent marginal from a global genic average destroys **78 %** of the signal. gDNA is biologically 50/50 and that is what the code uses |
| **the two "splash" pools** | ⛔ `DESIGN.md` §8 calls them pure gDNA and **neither is** — one is statistically indistinguishable from the pure RNA pool. Nothing consumes them, so nothing is mis-fitted today, but the design's claim is wrong as stated and should be marked unverified until settled |
| **`POOL_EB_PRIOR_ESS = 1000`** | It now **dominates** the gDNA length model — the pure pool is 4,467 fragments, so the prior pulls the fitted mean 79.4 → 100.3, a **+26 %** shift toward RNA. ⚠ It is a magic number and changing it is the owner's call |

---

## §7 ⛔ THE RULES THAT MADE THE NUMBERS ABOVE TRUSTWORTHY

Not process decoration — each of these caught something real in the last session.

1. **A falsification test written FIRST, verified failing, and then the code broken to watch the gate
   fire.** ⭐ In the second-pass work, **four of seven fixture loci exist only because a perturbation
   passed** — the gate was blind and said nothing until the mirror case was added.
2. **No magic numbers.** Stop and discuss. Two candidate constants were proposed in the last session and
   both were replaced by something derived: the combination rule and the ∅ strand term.
3. **Score against truth, not against the previous run.** ⭐ The simulator writes per-fragment ground truth
   into the oracle BAM's read names; that is what found the annihilation bug, and it is what turned "the
   deliverable improved" into "the deliverable got 23.9 % worse".
4. ⛔ **Real data is a TEST input, never a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
   spectrum. Sweep the plausible space, report the worst case, bring the domain call to the owner.
5. **The code does not reference these documents.** Docs evolve and rot; source must stand alone. A
   docstring may cite a *test* or the executable specification, because those cannot rot silently.
