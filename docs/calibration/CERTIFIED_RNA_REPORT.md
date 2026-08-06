# Can a spliced read tell us how much DNA contamination is nearby? — a negative result

    Written 2026-08-05 on branch `fragment-length-gold-standard` (base `87a5468e`).
    ⚠ **EPHEMERAL WORKING DOC**, like `SESSION_HANDOFF.md` beside it. The six permanent docs are
      `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⚠ **Audience: an external reviewer.** `SESSION_HANDOFF.md` §8 carries the same finding in the
      project's own vocabulary and is the version to act from. This document explains it from scratch.
    ⛔ Delete both when the question is closed.

---

## 1. What the tool is doing, in plain terms

Rigel measures how much of each RNA transcript is present in a sequencing library. Its complication is
that a real library is a **mixture of two kinds of molecule**:

* **RNA** — what we want to measure;
* **genomic DNA (gDNA)** — contamination that was never RNA at all, but which sequences into reads that
  land in exactly the same places.

Nothing in a single read says which it is. So before Rigel can quantify transcripts it has to run a
**calibration** step that estimates, at every position in the genome, what fraction of the local reads
came from DNA rather than RNA. That fraction is called `f_g` ("fraction gDNA"). Getting it wrong in
either direction corrupts every downstream number.

### The unit of measurement

The genome is chopped into **regions** (an exon, an intron, an intergenic stretch) with **lines** between
them. For each line, the tool counts fragments in three separate buckets:

| bucket | what it holds | can it be DNA? |
|---|---|---|
| **unspliced** | fragments that ran straight across the line in one piece | **yes — this is the mixture we must untangle** |
| **spliced** | fragments that ran straight across the line *but were spliced somewhere else along their length* | **no** |
| **junction** | fragments that *jumped* from this line to a distant one | **no** |

The second and third buckets are special. **Splicing is something only RNA does.** DNA has no splicing
machinery, so a fragment that shows a splice is RNA with certainty — no statistics required. The project
calls these **certified RNA**.

`f_g` is defined only for the **first** bucket. It is the fraction of the *unspliced* count that is DNA.
That definition matters a great deal below, and confusing it with "the DNA fraction of everything at this
line" has already caused one mis-report in this codebase.

---

## 2. The opportunity, and why it looked compelling

The core solver (`ψ`, a small per-position Bayesian model) **never looks at the certified buckets at all**.
Its own source comment says so. The certified counts are used for bookkeeping — added back into the final
RNA total — but they never inform the estimate of `f_g` at the position where they were observed.

That looks like free information left on the table, and the size is not marginal: on the 36-condition
benchmark panel, **5,991 of 35,038 active lines carry certified spliced reads, together covering 85.3 %
of the unspliced read mass** (off hybrid capture).

And there is a concrete failure that the missing channel appears to explain. On a small test scenario
(`tes_readthrough`), one line carries **365 unspliced reads and 11,026 certified-RNA reads**. Of those 365,
the truth is that **11 are DNA**. The tool predicts **186** — a 17-fold over-call — because with no
evidence of its own it falls back to "I don't know", which is `f_g ≈ 0.5`.

Eleven thousand reads at that position are *known* to be RNA, and the tool answers "coin flip". So the
proposed fix was: let the certified count `S` speak.

---

## 3. The proposed fix, and the algebra behind it

The solver already contains a placeholder for exactly this. Its prior has a term
`½ · log(1 − f_g)` — a conventional weak default meaning "no information". The proposal was to replace the
½ with `½ + S`:

```
no evidence:      ψ_RNA  =  (½)     · log(1 − f_g)      ← today
S certified:      ψ_RNA  =  (½ + S) · log(1 − f_g)      ← the proposal
```

This is not a bolt-on. It is precisely the Bayesian update you would write if the `S` certified reads
were `S` observations of "this position is RNA". **I verified that this identity is exact** (gate `C4`):
the resulting distribution for `f_g` is exactly `Beta(½, ½+S)`, and at `S = 0` it returns today's
behaviour bit for bit. There is no new normalisation to get wrong and no special case for `S = 0`.

### Where the coefficient comes from

Write down the expected number of certified reads. Let `M` be the unspliced count at this line, so
`(1 − f_g)·M` of them are RNA. Let **`q`** be the probability that an RNA fragment crossing this line
*also shows a visible splice* somewhere along its length — i.e. that it lands in the certified bucket
rather than the unspliced one. Then:

```
E[S]  =  (q / (1 − q)) · (1 − f_g) · M
```

Taking logs of the Poisson likelihood:

```
log P(S)  =  S·log(1 − f_g)   −   (q/(1−q))·(1 − f_g)·M   +   constant
             \___ retained ___/    \_______ dropped _______/
```

The proposal keeps the first term and drops the second, on the argument that the first is the honest
"one-sided" reading — evidence that RNA is *present*, which can push `f_g` down but never up.

### One thing this derivation gets right, and it is worth stating

There was a live worry (recorded as trap **C0**) that the two certified buckets use *different* divisors —
spliced reads are normalised one way, junction reads another — and that a term built from both could be
inconsistent.

**That worry does not apply here, and the reason is structural.** Any such normalisation factor
*multiplies* `E[S]`. Multiplying inside a logarithm adds a constant, and a constant cannot change the
coefficient of `log(1 − f_g)`. So the retained term's coefficient is the **raw count `S`**, and neither
divisor appears in it at all. There is nothing for them to disagree about, in either direction.

I gated this by brute force against SciPy's own Poisson implementation across a grid of normalisation
factors and read counts (gate `C1`), with the wrong alternative (`c·S`) confirmed to break it.

---

## 4. What I actually did, in order

| | step | outcome |
|---|---|---|
| 1 | **Wire the second pass into the test harness** | ✅ done and gated |
| 2 | **Settle the open derivation questions before writing code** | ⛔ **the derivation is refuted** |
| 3 | **Write the gates** | ✅ 28 gates, 6 perturbations, all firing |
| 4 | Implement the term | ⛔ **not done — step 2 forbids it** |
| 5 | Measure on the toy scenario | superseded by step 6 |
| 6 | Measure on the full benchmark panel | ✅ **ran first, and is what produced the negative** |

### Step 1 — a data-quality fix that had to come first

A small fraction of reads (~1 %) can't be assigned to a unique genomic path from the sequence alone,
because the unsequenced gap between paired reads admits more than one explanation. A "second pass" resolves
these by drawing one explanation at random with the right probabilities.

The problem: the tool's **ground-truth oracle** works by splitting the reads by their known true origin
(DNA / mature RNA / nascent RNA), re-counting each group separately, and checking the three add back up to
the whole. That check is what makes the oracle trustworthy. But the second pass breaks it — its
probabilities are computed against the *whole* library, so running it separately on three subsets is a
different operation and the three no longer sum.

The repair (already written in a previous session) is to make the random draw **once** on the whole
library and then *replay each fragment's already-chosen answer* inside whichever subset holds it. I wired
that through:

* `pipeline._drain_side_buffer` now publishes what the replay needs;
* the oracle accepts it and drains each subset with the replayed choices;
* the existing "three must sum to the whole" check becomes the end-to-end proof that this is correct, for
  free.

**Result on the target scenario:** 36 fragments resolved, certified-RNA count **224 → 269 (+20 %)**, and
**zero** cases where the replay was ambiguous. Two deliberate sabotages both fail loudly, as designed.

⭐ **An unexpected finding here.** I checked whether this fix could have changed any *earlier* result — the
discipline being that "nothing changed" is only meaningful if something *could* have. On every previous
test scenario the count of unresolved fragments is **zero**: you need two competing splice paths through
one gap, and only the newest scenario has that. So every earlier number in this campaign was already
correct, and this fix matters for exactly one scenario — the one it was built for.

---

## 5. Step 2 — the negative result

### The question I had to answer

Dropping the second term is safe **only if that term is small**, and it is small only if `q → 1` — i.e.
only if nearly every RNA fragment crossing a line reveals a splice.

`q` is not something the tool knows. But it is directly *measurable* from the ground-truth oracle, because
the oracle knows which unspliced reads were really RNA:

```
q  =  S / (S + C_R)          C_R = the TRUE number of RNA reads in the unspliced bucket
```

So I measured it. All 36 benchmark conditions. **No solver runs at all** — this reads cached truth and
does arithmetic, so the whole panel takes seconds.

### Finding 1 — `q` is nowhere near 1

Read-weighted median `q` per condition: **0.19 to 0.71**. And **60–98 % of the read mass sits at
`q < 0.9`** in *every* condition.

So the dropped term is not a correction. It is the same size as the term being kept. The "one-sided floor"
reading is not conservative — it is simply an unmodelled half.

### Finding 2 — the term is worse than doing nothing on a third of the panel

I scored what the term alone would answer (the median of `Beta(½, ½+S)`) against per-position truth,
weighted by read mass, over exactly the positions the term would fire on. "Reference" is the
uninformative `f_g = 0.5` such a position gets today.

| condition | median `q` | term's error | reference's error | change |
|---|---|---|---|---|
| `g00 ss0.50 capture_off` | 0.200 | 0.0107 | 0.5000 | **−0.4893** |
| `g50 ss0.50 capture_off` | 0.234 | 0.0695 | 0.4540 | −0.3845 |
| `g50 ss0.50 capture_on` | 0.356 | 0.5245 | 0.3285 | **+0.1960** |
| `g90 ss0.50 capture_on` | 0.600 | 0.8654 | 0.4076 | **+0.4578** |
| `g98 ss0.50 capture_on` | 0.694 | 0.9230 | 0.4779 | +0.4451 |

(`g00 … g98` is the DNA-contamination level, 0 % to 98 %. `capture_on/off` is whether hybrid-capture
enrichment was used. Lower error is better.)

**Worse than the uninformative default on 12 of 36 conditions.** Worst case `+0.4578`, where **98.3 % of
the affected read mass** sits on positions whose true `f_g` is **0.84** and where the term answers **0.04**
— confidently, and almost exactly backwards.

### Finding 3 — the tell, and this is the part that generalises

Look at the extremes. At `g00` there is **no DNA at all**, so the truth is "everything is RNA" — and the
term is the best result on the panel (−0.49). At `g98` the library is **almost entirely DNA** — and the
term is the worst (+0.45).

**The sign of the benefit is the sign of the answer.** A channel that helps precisely when the truth
happens to be what it always says has not added information; it has added a *prior* — a standing bias
toward "RNA". It looks like a large accuracy gain on half the panel and it is not one.

This is the same shape as a trap the project has already recorded (`A12`), and the two extreme conditions
of the panel function as the two "zero controls" the owner requires on every experiment: one where the
truth is constant-RNA, one where it is constant-DNA. Having both is what makes this a diagnosis rather
than a suspicion.

### Why it fails — the mechanism, stated exactly

With `q` unknown, `q/(1−q)` can take **any positive value**. So for *any* `f_g` less than 1, there is some
`q` that makes the model predict exactly the `S` we observed. The observation is therefore compatible with
every possible answer.

Formally: the profile likelihood in `f_g` — the best the model can do at each candidate `f_g`, optimising
over the unknown `q` — is **exactly flat** on the whole interval `[0, 1)`. I verified this by brute force
over 26 orders of magnitude of `q`, stated as a convergence test so that a loose tolerance can't hide a
real tilt (gate `C2`).

> **A certified read count carries exactly one bit of information about the DNA fraction:
> "it isn't 100 % DNA."**
> Everything the `S` coefficient claims beyond that one bit is an assumption about `q`, not evidence.

That is why the failure is concentrated under hybrid capture and at high DNA. Capture concentrates DNA
onto exons, so an exon-interior position genuinely *is* DNA-rich **and** carries plenty of certified RNA at
the same time. The term reads the RNA, has no way to see the DNA, and answers "RNA".

### A corollary that undercuts the original motivation

The investigation began partly because a *zero* certified count at a silent gene looked like the strongest
possible evidence that a position is pure DNA. **It isn't.** With `q` free, `S = 0` is equally well
explained at every `f_g` **including 1** — just take `q → 0`. The profile is flat on the *closed* interval
(gate `C2b`). All the information in this channel is in `S > 0`, and it points *away* from "pure DNA",
never toward it.

### A separate defect found on the way

The implementation plan also specified reusing an existing helper, `density_factor_precision`, to compute
how *confident* the new term should make the solver. That helper measures the width of a probability
distribution on a numerical grid, which is correct for a distribution with a peak. The certified term has
no peak — it is monotone — so what the helper actually measures is **the width of the grid**.

Measured at `S = 10,000`: the reported "confidence" is **756.6** on a narrow grid and **0.106** on a wide
one — a **7,100× swing** from changing a numerical setting that is supposed to be irrelevant. On a
genuinely peaked input the same helper returns exactly the right answer (1.0 and 25.0), unchanged across
the same range, which confirms the helper is fine and the *input* is out of contract.

The solver module's own stated acceptance test is invariance to that grid setting, so this is disqualified
by the project's existing rule. The correct form is available in closed form and involves no grid at all
(gate `C5`).

---

## 6. What survives

Not everything here is negative, and the surviving parts are what a future attempt should build on.

1. **The coefficient is the raw count, and no normalisation factor can reach it.** The `C0` worry is
   structurally inapplicable. (gate `C1`)
2. **The Bayesian identity is exact.** Prior + term = `Beta(½, ½+S)`; `S = 0` reproduces today's behaviour
   bit for bit. No new normalisation, no special case. (gate `C4`)
3. **`q` is a property of RNA geometry, not of contamination.** Per-position `q` measured at 0 % DNA
   versus 50 % DNA agrees at Spearman **ρ = 0.9257** across 5,241 shared positions. A 50-point swing in
   contamination barely moves it.

   ⭐ **This is the most important constructive result in the report.** It means `q` is in principle a
   *computable* quantity — a function of where the splice sites are and how long the fragments are — not
   an unknowable biological nuisance. The channel is **blocked, not impossible**.
4. **With `q` supplied, the term works.** The same likelihood, with `q` pinned, has a sharp interior
   maximum that recovers the true DNA fraction to within 0.01. (gate `C2`'s perturbation)
5. **The prize is large.** Certified reads cover 85.3 % of the relevant read mass off capture, and where
   the term *is* legitimate it is worth −0.49 in error — the largest single move available.

---

## 7. What to do next

### 7.1 The immediate decision: one cheap measurement, then a go/no-go

`1 − q` is "the chance a fragment crossing this line fits inside the unbroken stretch of exon on either
side of it". If it fits, it never had to splice; if it doesn't, it did. So:

```
q_geom  =  1  −  E_reach / E_any
```

where `E_any` is the number of ways a fragment can cross the line at all, and `E_reach` the number of ways
it can cross while staying inside the unbroken exon. Both come from one fragment-length distribution via
the existing `crossing_eff_length` function; `E_reach` needs the distance from the line to the nearest
splice site (or transcript end) in each direction.

**The measurement:** compute `q_geom` for every position, and regress the *measured* `q` on it.
`certified_q_census.py` already emits the measured `q` per position, so this is one script and no solver
time.

* **If `q_geom` explains most of the variance** → the channel is buildable. Proceed to 7.2.
* **If it does not** → the missing ingredient is transcript-level abundance, which calibration does not
  know at this stage. The channel is **closed**, not merely blocked, and this document becomes the reason
  it is never attempted again.

⚠ There is a specific reason to expect difficulty, and it should be checked rather than assumed. The
existing "reach" array in the codebase (`build_contiguous_edge_reach_arrays`) is **genomic**, not exonic —
deliberately, because these lines are also crossed by *nascent* (unfinished, unspliced) RNA, which runs
into introns freely and does not splice. The exonic version needed here is per-transcript, and where
several transcripts overlap, the right answer is an abundance-weighted blend of their different reaches.
Abundance is what the downstream step computes, so there is a potential circularity to resolve.

### 7.2 If it is buildable — the shape of the work

1. **Build the exonic-reach array**, per line and per side, as new geometry. Under the project's standing
   rule, any new divisor must be verified against brute-force enumeration — literally counting the valid
   fragment placements for small cases and checking the formula reproduces them.
2. **Restore the dropped term** rather than approximating it away. The estimator becomes two-sided with a
   finite confidence and an interior maximum at

   ```
   1 − f_g  =  S · E_reach / ( (E_any − E_reach) · M )
   ```

   which is a genuine point estimate of the RNA content of the unspliced bucket, from the certified count
   and geometry alone.
3. **Use the closed-form confidence**, `S · f_g · (1 − f_g)`, not the grid helper — see §5's last
   subsection. This matches the form the existing strand-evidence channel already uses.
4. **Re-run `certified_q_census.py` as the acceptance test.** It scores the proposed answer against truth
   on all 36 conditions in seconds, with no solver in the loop. The requirement is that the two extreme
   conditions (`g00`, `g98`) move in the *same* direction — that is the specific failure this report
   documents, and it is cheap to check continuously rather than at the end.
5. **Only then** run the full solver arm, and only then write into `src/`.

### 7.3 A separate, smaller item that is now unblocked

Independently of `q`, the certified count is currently absent from the "level" each position reports to its
neighbours (`node_total_density` — its own comment records this as a deliberate deferral). That is a
different quantity from the one this report is about: it needs no `q`, because the certified reads share
their normalisation with the unspliced ones by construction. It is worth pricing on its own, but the
project has already ruled that this must be measured jointly with a known hybrid-capture defect rather
than alone.

### 7.4 What must NOT be done

* ⛔ **Do not implement the `(½ + S)` coefficient.** It is measured at +0.4578 error on the worst
  condition. The gates in `test_certified_rna_licence.py` will fail if someone tries.
* ⛔ **Do not re-point gate `G6`** in `test_vertex_reference.py`. It asserts the solver is blind to this
  channel, and the solver *is* still blind — that assertion is currently true and should stay green until
  a legitimate fix lands. I added a note in its docstring so the next attempt starts from this constraint
  rather than rediscovering it.
* ⛔ **Do not add a threshold** ("only apply the term where the DNA level is low"). It would work on the
  panel and it is exactly the kind of tuned constant this project has refused three times. The problem is
  a missing quantity, not a missing cutoff.

---

## 8. Artefacts, and how to reproduce

| file | what it is |
|---|---|
| `scripts/design/certified_q_census.py` | **the instrument.** Measures `q` and scores the proposed term against truth on all 36 conditions. No solver. Verdict in its docstring |
| `tests/calibration/test_certified_rna_licence.py` | **28 gates** pinning what a certified count may and may not claim. Scored against SciPy's own Poisson and Beta rather than re-derived arithmetic |
| `docs/calibration/SESSION_HANDOFF.md` §8 | the same finding in project vocabulary — the version to act from |
| `docs/TRAPS.md` C0b | the transferable lesson |
| `src/rigel/pipeline.py`, `tests/calibration/_oracle.py`, `scripts/design/toy_harness.py` | step 1's wiring |

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1

python scripts/design/certified_q_census.py          # the measurement, whole panel, seconds
python -m pytest tests/calibration/test_certified_rna_licence.py -q     # 28 gates
python -m pytest tests/ -q                            # 22 failed / 2268 passed / 3 xfailed
```

The 22 failures are the project's standing baseline (21 stored-output comparisons plus one known
multi-mapping case), unchanged by this work. `ruff check src/ tests/ scripts/` is clean.

### On confidence in the gates

Writing a test that passes proves nothing on its own. Every gate here was also **deliberately broken** and
watched to fail:

| sabotage | gate that fired |
|---|---|
| coefficient scaled by the normalisation factor | `C1` |
| the unknown `q` restricted to a narrow band | `C2` |
| the zero-count claim asked about a non-zero count | `C2b` |
| the coefficient zeroed | `C3` |
| the count added to the DNA side instead of the RNA side | `C4` |
| the grid helper handed a peaked input instead of a monotone one | `C5` |

All six fired, each on exactly its intended gate and nothing else. No sabotage passed silently — which,
on this project's own history, is the failure mode worth checking for.

---

## 9. Summary for a reviewer

A spliced sequencing read cannot come from DNA, so it is proof that RNA is present. The tool was ignoring
that proof, and a natural fix was to feed the count of such reads into its estimate of local DNA
contamination.

**That fix does not work, and the reason is not an implementation detail.** The count of *visibly* spliced
reads depends on two things that are multiplied together: how much RNA is there, and how likely a fragment
was to reveal its splice. The second factor is pure geometry — it depends on how far the nearest splice
site is compared with the fragment length — and it varies from near 0 to near 1 across the genome. Without
knowing it, the observation is mathematically compatible with **any** contamination level; it constrains
only the extreme case of 100 % DNA.

Assuming the factor away amounts to assuming the answer. Measured against ground truth on 36 benchmark
conditions, the resulting estimate is better than doing nothing where the truth happens to be RNA-rich and
substantially worse where it is DNA-rich — worse overall on a third of the panel.

The cost of finding this out was one script and no solver time, because the benchmark's ground truth is
cached. **It was found before any of it was written into the production code**, which is the outcome the
project's own process rule was designed to produce.

The channel is not dead. The missing geometric factor was shown to be a property of the annotation rather
than of the contamination (correlation 0.93 across a 50-point contamination swing), so it is computable in
principle. §7.1 gives the one cheap measurement that decides whether it is computable in practice.
