# PLAN — the model-free TOTAL abundance at a REGION, and the bank that already carries it

    ⚠ **A DEV DOC. Nothing outside `docs/dev/` may cite it, and it is NOT the state.** The state is
    `ROADMAP.md` §0/§1; rulings go to `DESIGN.md`, derivations to `EQUATIONS.md`, lessons to
    `TRAPS.md`, numbers to `ROADMAP.md`. §9 below is the MOVE list — execute it and delete the
    corresponding section here in the same edit.

    Written 2026-08-20 (session 3) at the owner's request: *"my hope is that we can find or re-derive
    the model free abundance for regions and for boundaries."* No `src/` change was made. Nothing was
    committed. Every measurement below is a read, a brute-force enumeration, or an arm over a cached
    payload.

    ⛔ **This is a DIFFERENT THREAD from `session_2026-08-20_dissection_loop.md`.** That one is the
    owner's iterative panel-dissection protocol and it remains the primary loop. This one is an
    ACCUMULATOR / DEPOSIT-RULE thread that came out of a question about what the accumulator deposits
    before calibration runs. They meet in exactly one place — §6, where the defect found here is a
    named mechanism for a result that loop measured.

---

## 0. THE WHOLE THING IN SIX LINES

1. The owner remembered a previously-derived *model-free abundance* for REGIONS, denoted with a
   capital `A`. **It was found, it is correct, and it already shipped** (§1).
2. It does **not** give a model-free TOTAL, because its opportunity is zero for `w > ell` — a
   fragment longer than the REGION can never be contained. `E[Sum 1/A] = rho * P(w <= ell)`, and
   `P(w <= ell)` is a functional of the fragment-length pmf, **per component** (§2).
3. At a 98 bp exon on the ladder's own fragment lengths that is a **11.6x under-read**; below 50 bp it
   is exactly **zero at any depth** (§2).
4. The truncation-free answer is a different RELATION: **the fragment's first covered base lies in the
   REGION**. Its opportunity is `ell` for every `w` — no fragment length appears in the weight at all
   (§3).
5. **That bank is already deposited** as `region_start_count` and is read by nothing downstream (§4).
6. It has exactly one real defect — the TEMPLATE WALL — and the exact repair is the mirror bank
   `region_end_count`, whose value is already computed and discarded (§5).

⭐ **FIRST COMMAND OF THE SESSION**, before any of this:

    python scripts/design/preflight.py            # ~2 min, or --fast
    python -m pytest tests/ -q                    # ANY failure is a regression

⛔ Do not start on a ✘. A missing DERIVED artifact is a command you have not run, not damage.
⚠ Suite baseline measured **3,614 collected** when this note was written (this file makes it 3,615 —
one `test_no_jargon_labels` case; a `docs/dev/` file gets `+1`, never the `+2` a `tests/calibration/`
file gets). ⛔ Re-derive, never adjust (`TRAPS: re-record-the-baseline`).

---

## 1. WHAT WAS FOUND — the capital-`A` derivation, and it already shipped

**Provenance.** Written `8c673324` (2026-07-30) as `docs/NODE_DENSITY_DERIVATION.md` §4.1–4.2; moved to
`docs/accumulator/` at `5d838580`; that doc tree was deleted at `c6b2ea89`. The same derivation survived
in a dev doc, `docs/dev/counts_densities_paradigm.md` §8, until `a5ed752b`. The adoption plan is
`be9b7efb` (*"the node deposit becomes 1/A"*), preceded by `e9bc8bf3` and `83644a7d`. **It landed in the
accumulator at `69a85be2` (2026-08-10).**

⚠ Historical vocabulary in all of those: `node` = REGION, `edge` / `line` / `seam` = BOUNDARY. Do not
carry those words forward (`DESIGN.md` §0).

**The derivation, as written.** Capital `A` is the OPPORTUNITY — the number of admissible start
positions for a length-`w` fragment at that object:

    A(w) = (ell - w + 1)_+     at a REGION of length ell     relation: CONTAINED
    A(w) = (w - 1)_+           at a BOUNDARY                 relation: CROSSES the point

    deposit  h(w) = 1 / A(w)

    E[ Sum h ] = Sum_w rho * f(w) * A(w) * (1/A(w)) = rho      "for ANY length pmf"

**Where it lives today.**

| object | deposit | site |
|---|---|---|
| BOUNDARY | `inv_boundary = length >= 2 ? 1.0/(length-1) : 0.0` | `accumulator.cpp:622` |
| BOUNDARY, unspliced crossing | `boundary.unspliced_inv_length_sum += inv_boundary` | `accumulator.cpp:643` |
| sj | `sj.inv_length_sum += inv_boundary` | `accumulator.cpp:711` |
| REGION, contained only | `region.contained_inv_opportunity_sum += 1.0/(region_len - length + 1)` | `accumulator.cpp:733` |

Executable specification: `tests/native/_accumulator_reference.py:927`, `:990`, `:1013`.
⚠ `69a85be2` also renamed the bank `..._inv_length_sum` -> `..._inv_opportunity_sum`, which is what
invalidated every cache. That rename was luck rather than design — see §8's cache-key note.

⭐ **The one correct statement of the derivation's LIMIT anywhere in the tree** is
`scripts/design/region_density_derivation.py:20`:

    T2  E[sum h] = rho * P(A > 0). Model-free up to the population's own support truncation.

**Nothing cites that line**, and every other home dropped the factor (§9).

---

## 2. THE DEFECT — the cancellation is conditional on its own support

`A_contained(w) = 0` for `w > ell`. Such a fragment is not merely down-weighted, it **deposits nothing
at all**, so the cancellation holds only on `{A > 0}`:

    E[ Sum 1/A ]  =  rho * Sum_{w : A(w) > 0} f(w)  =  rho * P(w <= ell)          <- not rho

`accumulator.cpp:728` states the conditioning as a reassurance — *"`A >= 1` is structural: the fragment
IS contained, so `w <= ell`"* — which is exactly the sentence that hid it.

**Measured on the ladder's own fragment lengths** (truncated normal mean 206, sd 98, on `[50, 500]`,
`gdna_ladder.yaml:131-134` / `:156-159`; realised deposited mean 217):

| `ell` | 50 | 98 | 120 | 151 | 250 | 400 | 500 |
|---|---|---|---|---|---|---|---|
| `P(w<=ell)` | 0.0012 | **0.0861** | 0.1445 | 0.2479 | 0.6571 | 0.9764 | 1.0000 |
| under-read | 823x | **11.6x** | 6.9x | 4.0x | 1.5x | 1.02x | 1.00x |

Median exon `ell` on the panel index is **98 bp**. Of 35,135 REGIONs, 24,018 are exons, and **7,123 of
those have `ell < 50 = frag_min`, so their contained bank is exactly ZERO by construction at any
depth.**

**Brute-forced on an exactly uniform field** — one fragment at every admissible start of every length,
which is `rho = 1` per length by construction, the same falsification style as
`tests/native/test_fragment_length_proof.py`:

    ell =   30   w in [40,80]  ("RNA-like")     contained 0.0000     start 1.000000
    ell =   30   w in [150,300] ("gDNA-like")   contained 0.0000     start 1.000000
    ell =   98   w in [40,80]                   contained 1.0000     start 1.000000
    ell =   98   w in [150,300]                 contained 0.0000     start 1.000000
    ell = 1000   either                         contained 1.0000     start 1.000000

⛔⛔ **AND THIS IS THE gDNA-vs-RNA PROBLEM RESTATED, NOT A SEPARATE ONE.** `P(w <= ell)` is per
COMPONENT. On a library with a real length gap (gDNA mean ~330, RNA ~120, lognormal) the two components'
truncation factors differ by **orders of magnitude at exon scale** — the ratio runs from ~10x at
`ell = 200` to three digits at `ell = 98`, and closes to 1.0 only above ~500 bp. Two independent
parameterisations of the gap were tried and both give that shape; the exact multiplier is a function of
the assumed pmf and must be re-derived from a real library rather than quoted from here.

So the circularity `messages/currency.py:194-195` was written to eliminate — *"NEVER
`mass / effective_length`. That divisor is a function of the composition being solved for."* — is still
present. It moved out of the DIVISOR and into the SUPPORT.

⭐⭐ **THE EQUAL-LENGTH LADDER SEES THE TRUNCATION IN FULL, AND THAT IS THE OPPOSITE OF WHAT IS WRITTEN
DOWN.** The bias is `w` against `ell`, not `w_gDNA` against `w_RNA`. Confirmed through the shipped banks
on the real cached payloads: `Sum(contained_inv) / Sum(region_start_count/ell)`, per class per `ell`
band, divided by the condition's own `P(w<=ell)` from its own `deposited_lengths`, reads **0.93–1.03 on
eleven of twelve live capture-OFF bins across all three REGION classes** at `g98 ss0.99`. The law is a
measurement on shipped code, not a derivation.

⚠ **Name the outlier rather than leaving it in the residual** (`TRAPS: read-the-whole-failure-list`):
the one bin outside the band is `exon` at `ell` in `[50,150)`, which reads **0.853** — the band where
`P(w<=ell)` is smallest and the quotient is therefore worst-conditioned. ⛔ And the twelfth bin
(`intergenic` at `ell` in `[50,150)`) was NOT reported at all; re-derive it, because an unreported bin
is not a passing bin.

⛔ Therefore `messages/currency.py:200-201` — *"An equal-fragment-length panel cannot see the
difference, which is why this had to be reasoned about rather than measured"* — **is false for this
half, and it is the sentence that will stop the next session from measuring it.** Delete it (§9).
⚠ It remains TRUE for the per-component half of §2: with `E_g == E_r` the ladder cannot separate
`P_g` from `P_r`. Two different claims under one sentence.

---

## 3. THE DERIVATION ASKED FOR — the truncation-free total

### 3.1 The family theorem, stated so the proviso is visible

REGION `R` = the half-open interval `[0, ell)`. A fragment of length `w` occupies `[s, s+w)`. Let
`lambda(s,w)` be the intensity of fragment starts. For a relation with admissible start set `S(w)`, put
`A(w) = |S(w)|`. Let `g(s,w) >= 0` be any placement functional computable at deposit time and
`C(w) = Sum_s g(s,w)`; deposit `h = g/C`. Then **by linearity alone — no Poisson assumption and no
independence assumption**:

    E[ Sum_frag h ]  =  Sum_s Sum_w lambda(s,w) * g(s,w) / C(w)

and if `lambda(s,w) = rho * f(w)` for every `s` in the support of `g(.,w)` and every `w` in the support
of `f`:

    E[ Sum h ]  =  rho * Sum_w f(w) * C(w)/C(w)  =  rho     for any f,  PROVIDED C(w) > 0 wherever f(w) > 0.

⭐ **The truncation is exactly the failure of that last proviso.** So the question is precisely: *is
there a relation whose `A(w)` is strictly positive across the whole support?*

### 3.2 The opportunity table

Closed forms verified against an independent integer-set oracle by exhaustive enumeration over
`ell` in `[1,12]` x `w` in `[1,15]` and again over `ell` in `[1,60)` x `w` in `[1,80)`, zero mismatches.
⚠ `ell = 0` is EXCLUDED: the interval degenerates and the formula/relation split there is §3.5's subject,
not a mismatch the oracle can see.
`o` denotes the overlap in bases.

| relation | `A(w)` | support | `E[Sum 1/A]` |
|---|---|---|---|
| CONTAINED in R | `(ell - w + 1)_+` | `w <= ell` | `rho*P(w<=ell)` — **SHIPPED at a REGION** |
| CONTAINS every base of R | `(w - ell + 1)_+` | `w >= ell` | `rho*P(w>=ell)` |
| SPANS R strictly | `(w - ell - 1)_+` | `w >= ell+2` | `rho*P(w>=ell+2)` |
| CROSSES a designated interior point | `(w - 1)_+` | `w >= 2` | `rho*P(w>=2)` — **SHIPPED at a BOUNDARY** |
| COVERS a designated base | `w` | all `w` | `rho` |
| OVERLAPS at least one base of R | `ell + w - 1` | all `w` | `rho` |
| MIDPOINT in R | `ell` | all `w` | `rho` |
| ⭐ **STARTS in R** (first covered base) | **`ell`** | **all `w`** | **`rho`** |
| ENDS in R (last covered base) | `ell` | all `w` | `rho` |

Atom-averaged forms, for completeness: COVER is `h = o/(ell*w)`, CROSS is `h = (o-1)_+/((ell-1)(w-1))`.
Composite identity, zero exceptions: `(ell + w - 1) = (ell-w+1)_+ + 2*(w-1)_+ - (w-ell-1)_+`.

### 3.3 Which one, and why it is unique

**FIVE** relations are truncation-free across the whole support — COVERS, OVERLAPS, MIDPOINT, STARTS-IN,
ENDS-IN — six if CROSSES is counted, whose `w >= 2` is vacuous at `frag_min = 50`. They separate on the
**UNIFORMITY WINDOW** `W`, the exact set of bases over which `lambda` must be constant for §3.1's proviso
to hold. Widths over the ladder's support `w` in `[50, 500]` — ⚠ two of them depend on `w_min` and not
`w_max`, and that matters: CONTAINED's 49 is `ell - w_min + 1`, so at `w_min = 1` it would read 98 and
the uniqueness premise below would evaporate; ENDS-IN's 548 is `ell - w_min + w_max`. MIDPOINT
(`ell + w_max//2 - w_min//2`) is omitted from the table because it is strictly dominated by STARTS-IN on
both axes:

| relation | window at `ell = 98` | at `ell = 500` |
|---|---|---|
| CONTAINED (shipped) | 49 | 451 |
| ⭐ **STARTS in R** | **98** | **500** |
| CROSSES a point | 499 | 499 |
| COVERS a base | 500 | 500 |
| ENDS in R | 548 | 950 |
| OVERLAPS R | 597 | 999 |

**UNIQUENESS — and state it as the MAXIMUM-OPPORTUNITY claim, which is what is true.** Confine the
assumption to the REGION itself: require `S(w) ⊆ [0, ell)`. Then `A(w) <= ell`, with equality for every
`w` only when `S(w) = [0, ell)` — i.e. *the fragment's first covered base lies in R*.

> ⭐⭐⭐ **STARTS-IN is the unique MAXIMUM-OPPORTUNITY relation confined to R, and among UNBIASED
> deposits supported on R it attains the minimum variance `rho/ell`** (Cauchy–Schwarz:
> `Sum_s h(s)^2 >= 1/|S(w)| >= 1/ell`, with equality iff `h = 1/ell` on the whole of R).

⛔ **Do NOT write "the unique truncation-free relation with `S(w) ⊆ R`" — that is false.** "The first
covered base lies in the first half of R" has `S(w) = [0, floor(ell/2))` for every `w`: a subset of R,
constant in `w`, truncation-free. There is an infinite such family; what singles STARTS-IN out is that it
is the largest of them. ⚠ And the variance ranking needs the word UNBIASED (`h = 0` has variance zero)
and a second-moment model; the unbiasedness itself needs linearity alone (§3.1).

### 3.4 The rule, and the proof

    A_start(w) = ell   for every w.        h(w) = 1/A(w) = 1/ell   is a CONSTANT.

    E[ Sum h ]  =  Sum_{s in R} Sum_w lambda(s,w) * (1/ell)
                =  (1/ell) * Sum_{s in R} lambda_bar(s)          lambda_bar(s) = Sum_w lambda(s,w)
                =  the REGION's MEAN START DENSITY over its own ell bases

    and under lambda(s,w) = rho*f(w):   =  (1/ell) * ell * rho  =  rho      EXACTLY.

Two properties make this stronger than "one more candidate":

1. ⭐ **No `w` appears in the weight.** There is no fragment length in the estimator, so there is
   nothing for a pmf to change. Model-freeness is STRUCTURAL rather than a cancellation that can fail
   on a support edge — which is the whole failure mode of §2.
2. ⭐ **`f(w)` never had to factor out of `lambda`.** The middle line holds for an *arbitrary*
   `lambda(s,w)`, so the estimator is exact for the REGION's mean observed start density under any
   field, hybrid capture included. Every other candidate needs `lambda` uniform over a window strictly
   wider than the REGION.

⭐ **And the deposit degenerates to an integer count** — the `/ell` happens at read time. That means
bit-identical merges across worker threads for free (`TRAPS: integer-channels-reproduce`), no float
bank, no new memory in the per-REGION struct.

⭐ It is also the estimand the consumer already assumes: `contained_eff_length`'s own docstring
(`effective_length.py:83`) calls it *"the count of starts placing the whole molecule inside"*, so
Rigel's `rho` at a REGION is already defined as fragments per admissible START position. A start-rate
REGION face is commensurate with a start-rate BOUNDARY face; a containment-weighted interior average is
not.

### 3.5 The BOUNDARY case — and `ell -> 0` is NOT a limit relation

The shipped BOUNDARY deposit is its own relation: *crosses a designated inter-base point*,
`A(w) = (w-1)_+`, `E = rho * P(w >= 2)`, which is `1.000000` on the ladder. Its uniformity window is the
`w_max - 1` bases immediately **LEFT** of the point: a 0-bp object has no interior, so its estimator is
one-sided by construction.

    ell -> 0 :   CONTAINED  -> A = (0 - w + 1)_+ = 0  for all w >= 2      <- NOT w-1
                 OVERLAPS   -> the FORMULA gives w-1; the RELATION gives 0
                 CROSSES    -> A = w-1 at EVERY ell; no limit is needed
                 STARTS-IN  -> A = 0

⛔ **`accumulator.cpp:726` and `_accumulator_reference.py:1008` both attach the `ell -> 0` limit claim
to the CONTAINED rule, whose limit is 0.** The claim is false as attached. The BOUNDARY deposit needs no
limit argument; it is a different relation at every `ell`. (The only relation whose `ell -> 0` limit is
both `w-1` and the boundary relation is SPANS-STRICTLY.)

⭐ **The useful exact result about BOUNDARIES**, which is what kills the "just use the flanks" idea in
§7: every admissible start of a crossing fragment is strictly LEFT of the point, so the RIGHT-adjacent
REGION contributes exactly ZERO. A BOUNDARY's `inv_length_sum` estimates the start density of its
LEFT-adjacent REGION, blended with whatever lies BEYOND THAT ON THE LEFT:

    E[ boundary deposit, point ell bases from the left REGION's far edge ]
        =  rho_adj * (1 - c(ell))  +  rho_beyond_left * c(ell)
        c(ell) = E_f[ (1 - ell/(w-1))_+ ]                <- A FUNCTIONAL OF THE PMF

⚠ `c(ell)` is the fraction of the start window that falls off the left end of the adjacent REGION.
§7.2's `boundary-right = 0.5752` is `1.0*(1 - 0.4720) + 0.1*0.4720` — the `0.1` is the LEFT flank, not
the right one.

    ladder fl:  c(50)=0.7155  c(98)=0.4720  c(151)=0.2681  c(250)=0.0624  c(400)=0.0016  c(>=499)=0.0000

Reproduced by two independent implementations and by exact enumeration.

---

## 4. ⭐⭐⭐ IT IS ALREADY DEPOSITED, AND IT IS READ BY NOTHING

```cpp
// src/rigel/native/calibration/accumulator.cpp:602-606
const std::int64_t first_base = scratch.segments.front().first;
const std::int64_t last_base  = scratch.segments.back().second - 1;   // :603 — computed, then discarded

const std::int64_t first_region = region_of_pos(first_base);
region_start_count_[static_cast<std::size_t>(first_region)] += 1u;
```
⭐ `:603` is shown deliberately: §5's mirror bank is `region_of_pos(last_base)`, and that value is
already here.

⚠ Note the comment already there: *"The path's own first and last COVERED base, not the fragment's
extent"* — so it is the right quantity for a spliced path too.

| stage | site |
|---|---|
| deposit | `accumulator.cpp:606`; merge `:865`; init `:266` |
| executable specification | `_accumulator_reference.py:900` |
| export | `bam_scanner.cpp:2224`, prop at `:2910` |
| typed payload | `scan_payload.py:509`, `:312`, `:761-764`, `:814` |
| into calibration | `substrate.py:149`, `:200` |

⛔ **A grep over `src/` returns those sites and the invariant `sum(region_start_count) == qc.deposited`
(`accumulator.h:475`, `scan_payload.py:119`) — plus `scan_payload.py:614`, which reads only its SHAPE
to get `n_regions` — and NOTHING ELSE. It never reaches `region_geometry.py`.** Every one of the 16 cached ladder payloads and all **16** test-chromosome scenarios
already contain it (`test_reference.yaml:56-57` is 4 rates x 2 strand x 2 capture; the `8` in earlier
notes is the pre-2026-08-20 grid).

The read side is one line, and `region_len` is already in scope at `region_geometry.py:257`:

```python
# region_geometry.py:296 today:
inv_abundance[is_region] = np.asarray(_r_inv, np.float64)[obj[is_region]]
# the swap:
inv_abundance[is_region] = (substrate.region_start_count[obj[is_region]]
                            / region_len[obj[is_region]])
```

---

## 5. THE ONE REAL DEFECT — the TEMPLATE WALL, and its exact repair

`A_start(w) = ell` holds only while every start slot in R admits a length-`w` molecule. Within distance
`d` of the template's genomic-HIGH end — the chromosome wall for gDNA, the TES or TSS whichever is
genomically higher for RNA, the reference end for an ERCC spike-in — the true opportunity becomes

    A_start(w | d)  =  min( ell , (d + ell - w + 1)_+ )        <- w-DEPENDENT again

and ⛔ **at `d = 0` that is exactly `(ell - w + 1)_+`, the CONTAINED opportunity.** So at a flush wall
the start rule DEGENERATES INTO THE RULE IT REPLACES while still depositing the flat `1/ell`, and reads
worse than it. Brute-forced, `ell = 98`, ladder fl, truth 1.0 everywhere:

| `d` | 0 | 25 | 50 | 100 | 200 | 400 | 600 |
|---|---|---|---|---|---|---|---|
| start | **0.0196** | 0.0501 | 0.1000 | 0.2507 | 0.6482 | 0.9926 | 1.0000 |
| end | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| contained | 0.0861 | 0.0861 | 0.0861 | 0.0861 | 0.0861 | 0.0861 | 0.0861 |

At a flush 98 bp exon the start rule is **4.4x worse** than the bank it replaces, and its error factor
is composition-dependent — gDNA has no transcript wall and reads 1.000, RNA does not. Equal fragment
lengths do not rescue it: what differs is the WALL, not `f`.

⭐⭐ **THE REPAIR IS EXACT AND IT IS THE MIRROR RELATION.** `region_end_count` at
`region_of_pos(last_base)` — and `last_base` is **already computed** at `accumulator.cpp:603` and used
only for the containment predicate at `:719`. Its opportunity is also `ell` and it is bounded by the LOW
wall instead. The two fail at opposite walls and each is exact at the other.

> **RULE: use the side whose wall does not bind. The estimator is exact iff `d_high >= w_max - 1`
> (start form) or `d_low >= w_max - 1` (end form).**

This is the REGION analogue of the `reach_lo`/`reach_hi` trapezoid Rigel already applies at BOUNDARIES
(`EQUATIONS.md` §1.5, `effective_length.crossing_eff_length`). ⛔ When BOTH walls bind — `d_low` and `d_high` both `< w_max - 1 = 499`, which is GUARANTEED for any
REGION inside a template shorter than 499 bp: **seven** of the 94 references (273, 274, 274, 274, 493,
494, 494), four of them `<= 274 bp` — **no deposit rule is model-free there** and the bound must be
RECORDED rather than papered over. ⚠ `template < 2*w_max - 2 = 998` is the NECESSARY condition, not this
one, and it admits 48 of 94.

⭐ **Exposure, RE-DERIVED 2026-08-20 in SPLICED-TRANSCRIPT coordinates (RANK 3 executed), and the
genomic 27.30 % is retired** — it used distance to the NEAREST transcript end, which is the wrong
collapse. With the named conservative collapse (`d_high(r) = MAX` downstream spliced template over the
templates covering `r`, so the wall binds only if it binds for ALL of them), weighted by
`region_start_count`:

| population of templates | FLUSH `d_high = 0` | BINDING `d_high < 499` | condition |
|---|---|---|---|
| ALL (nascent entities included) | **3.69 %** of exonic starts | **7.38 %** | `g00 ss0.99 off` |
| ALL | 7.21 % | 12.17 % | `g50 ss0.99 on` |
| MATURE templates only | **9.12 %** (7.47 % of library) | **22.62 %** (18.53 %) | `g00 ss0.99 off` |

⛔ **The wall is COMPONENT-DIFFERENTIAL, which is the finding**: gDNA's template is the chromosome and
never binds; the MATURE population is wall-bound at up to ~23 % of exonic starts; the nascent entities'
genomic reach rescues the pooled number to ~7 %. A per-component differential bias is the same
composition circularity as §2's truncation, relocated to TERMINAL exons — so RANK 4 is a CORRECTNESS
requirement for the mature population there, not an optimisation, and it is bounded (the two rules fail
on complementary populations). The exposure curve is smooth in `d_high` (no threshold artifact); the
enumeration is ~80 lines off `regions_df` + `get_exon_intervals`, re-derivable from this paragraph.

⚠ Also: the shipped REGION divisor `contained_eff_length(region_len, pmf)` (`region_geometry.py:268`)
takes **no reach argument at all**, which is correct for the CONTAINED relation and would be a hole for
the START relation. The mirror bank is the right repair; a taper on the divisor is not.

---

## 6. WHAT THIS EXPLAINS — a named mechanism for the CurrencyPolicy result

`region_geometry.py:288-299` (a block COMMENT, not a docstring) builds ONE flat `inv_abundance` array in which REGION slots carry the
**truncated** bank and BOUNDARY slots the **untruncated** one, under the docstring *"No divisor is
applied here and none may be: the bank IS the density"* — which is false at a REGION.

The chain strictly alternates REGION / BOUNDARY, so **every hop in `enrichment_ratio`
(`messages/currency.py:191`) is a REGION-to-BOUNDARY ratio, and every one carries a factor
`1/P(w <= ell_region)` with no enrichment in it at all** — 11.6x at a 98 bp exon, 823x at 50 bp,
undefined below `frag_min`. ⚠ **But the factor multiplies only the `inv` TERM of each face, not the
face total**: `currency.py:204-206` builds each face as `inv + inv_sj_lo.sum(axis=1)` (and `_hi`), and
the sj flux is a BOUNDARY-form `1/(w-1)` bank that enters UNDIVIDED — so on an expressed exon it dilutes
the truncated term and the in-situ distortion at any single hop is far smaller than 11.6x. That is why
the measured median below is 0.7424 rather than `1/11.6 = 0.086`. That feeds `rescale_weight`
(`currency.py:222`),
`w = (log r)^2 / ((log r)^2 + v)`; at exon depths `v` is small, so `w -> 1` and the policy believes the
artefact in full.

Measured per hop at `g98 ss0.99 capture_off`, where the true **gDNA** field is uniform so the
gDNA-level ratio is 1.0 everywhere: exon-destination median `r` = **0.7424** under the shipped bank
against **1.0490** under the start bank, with **28,028 of 47,852 live exon hops (58.6 %) having the
REGION bank at exactly 0** — where `currency.py:219` returns its `1.0` "no evidence of enrichment"
default.

⛔ **Two caveats on that reading, and the next session must resolve both before quoting it.** ① The TOTAL
field is NOT uniform at `g98`: that rung is 98 % gDNA and 2 % RNA (`gdna_ladder.yaml:149-150`) and every
rung carries nascent RNA at a 20 % fragment share, so on expressed exons — the population this section is
about — the total ratio departs from 1.0. Score against `slot_truth.npz`, never against an assumed 1.0.
② If 58.6 % of hops return the hard-coded `1.0`, the median over ALL live exon hops is exactly 1.0, so
the 0.7424 must be restated with its population named (hops with a NON-ZERO REGION bank, or all live
hops).

⭐ That is a mechanism for the recorded finding that `CurrencyPolicy` is the WORST arm on every in-scope
stratum while winning only where the local solve is blind (`f295a313`, `13f14e97`): it hops at exons,
and exons are where the ratio is wrong. ⛔ It is a mechanism, **not** a prediction that fixing it makes
the policy win — see §8 RANK 2's ceiling note.

---

## 7. HONEST LIMITS — each as a bound or a mechanism, not a worry

**7.1 Model-free buys a LEVEL, never a SPLIT.** A model-free channel gives all three populations
`{gDNA, RNA+, RNA-}` the same coefficient, so its row is PROPORTIONAL TO `(1,1)` and it supplies exactly
one equation: the TOTAL. Measured discrimination efficiency of `Sum 1/A` alone: 0.113 at
`ell = 151`, 0.000 at `ell = 1000`. The owner asked for a LEVEL, so this is a licence rather than an
obstacle.

⛔ But the recovered doc's phrasing — *"its determinant against any other row is zero"* — is **false as
written**. With rows `[count ; Sum 1/A]` the matrix is `[[E_g, E_r],[1,1]]` and `det = E_g - E_r`, zero
only against another row proportional to `(1,1)`. `accumulator.h:32` carries the correct qualified form.

⚠⚠ **THEREFORE A SCOPE HAZARD EXISTS AND NEEDS AN OWNER RULING.** Pairing a truncation-free total
`T = rho_g + rho_r` with the conserved unspliced mass `M = rho_g*E_g + rho_r*E_r` gives
`det = E_g - E_r != 0` — **that pair IDENTIFIES the split**, and it is a fragment-length composition
channel by construction. Using the total as a LEVEL (what `enrichment_ratio` does) is clear of the
post-0.8.0 length-channel retirement; using it to improve `f_g` is not.
⚠ **And the hazard is UNMEASURABLE ON THIS PANEL BY CONSTRUCTION**: the ladder sets `E_g == E_r`, so
`det = 0` there. Pricing it needs §10.3's side panel — which is the same substrate the per-component
half of §2 needs.

**7.2 The tiling problem is what refuses the OVERLAP family.** Rigel's REGIONs tile the genome, so a
fragment hanging off the end samples a DIFFERENT `rho`. Exact expectations, a 98 bp REGION at
`rho = 1.0` between two 1 kb flanks at `rho = 0.1`:

    start 1.0000 | contained 0.0861 | overlap 0.4047 | cover 0.3473 | cross ~0.347
    midpoint 0.2174 | end-in 0.1176 | spans 0.0912 | boundary-left 0.1000 | boundary-right 0.5752

⚠ MIDPOINT uses `mid = s + w//2`; CROSS's third decimal is sensitive to the interior-point convention
and is quoted to three figures for that reason. ⛔ An earlier draft of this row read `midpoint 0.3509`,
which is CROSS evaluated at the single centre point 49 — MIDPOINT has `A = ell`, so its window slides
wholly into the left flank for every `w >= 2*ell` and it cannot sit near 0.35.

⭐ In the start-rate currency there IS no tiling: a fragment belongs to whichever REGION its first
covered base lies in, and the map is a partition — `sum(region_start_count) == qc.deposited` exactly
(`accumulator.h:475`).

⛔ And OVERLAP is worse than it looks. Its window extends into the flanking REGION, and the ANNOTATION
does not admit RNA of every strand there — `statics.free_pos` / `free_neg` are false for an intron flank
that no transcript occupies contiguously, while gDNA is genomically continuous and always admitted. So
the true overlap opportunity **differs between the populations sharing the slot**, and a channel whose
OPPORTUNITY depends on which population supplied the fragment is not model-free at all. ⚠ Stated as an
OPPORTUNITY on purpose (AXIOM 0): the question is what each population's chance of being there IS, not
which species a fragment belongs to.

**7.3 The flanking BOUNDARIES do not give this for free.** By §3.5 the right-flanking BOUNDARY returns
the REGION's own `rho` exactly iff `ell >= w_max - 1` — the same population on which containment already
works — and the left one carries the NEIGHBOUR's `rho`. On short REGIONs both are blended by `c(ell)`, a
pmf functional, so inverting the blend re-imports the exact length model the channel exists to
eliminate. **The two channels do not complement; they fail on the same population.**
⚠ The opportunities DO compose exactly
(`A_overlap = A_contained + A_leftBoundary + A_rightBoundary - A_spans`, zero exceptions) — but the
DEPOSITS do not, for the one-sidedness reason above. Do not mistake the identity for a recipe.

**7.4 Under capture no deposit rule delivers the REGION's latent enrichment.** `lambda_bar(s)` genuinely
varies inside the REGION, so "the abundance at this REGION" is not a scalar. The start rule delivers the
mean OBSERVED start density exactly — verified by exact enumeration at 113x enrichment in three probe
geometries (probe = the REGION; probe abutting; probe covering only the first 200 of 1000 bases):
`1.000000` in every cell, where overlap / cover / cross read 0.26–0.40 and the BOUNDARY bank reads
0.0089. Everything beyond that is a PLACEMENT model, which `EQUATIONS.md` §4.4 already rules out of
scope. That ruling stands untouched.

⛔ **One agent argued the reverse — that under capture CONTAINED is the placement-correct estimator and
START is the leaky one. The data refuses it.** CONTAINED's window `[0, ell-w]` structurally excludes the
REGION's last `w-1` start slots, which is exactly where capture concentrates mass at an intron abutting
an exon probe. On the real cache at `g98 ss0.99 capture_ON`, in bands where `P(w <= ell) = 1.0000`
exactly so the truncation is provably absent: long introns read `contained/start = 0.0803`, long
intergenic `0.1458`, long exons `1.0470` — against `1.0000 / 1.0005 / 1.0092` for the same bins at
capture-OFF. ⚠ Recorded because it is the one place the six independent passes disagreed.

**7.5 Two smaller bounds.** ① The deferred population is lost identically by both banks
(`accumulator.cpp:542` returns before the start count at `:606`) — 1.37 % at `g50` capture-OFF, 1.84 %
capture-ON, and every cached manifest reads `drain: null`, so every `scripts/design/` arm is missing them
equally. ② `length` is the path's covered-base count, not its genomic span
(`accumulator.cpp:438-458`), so for a fragment with an UNANNOTATED intron inside one REGION the
contained divisor uses the spliced `L` where the placement space is governed by the span. Zero on the
panel (`unannotated_introns: 0`), live on real data. `region_start_count` is immune — a fragment has one
first covered base whether it splices or not.

---

## 8. ⭐⭐⭐ WHAT THE NEXT SESSION DOES — ranked, with exact steps

⛔ **Read `TRAPS: panel-before-src` first.** Nothing in §3–§5 has been priced against the 0.8.0 metric.
This section is ordered so that everything free happens before anything irreversible.

### RANK 1 — land the doc and gate corrections. Cost: zero, and do it first.

§9 is the list. It costs nothing, it is unambiguously right, and until it lands the next session reads
`EQUATIONS.md` §3.5g and rebuilds this error. ⛔ In particular delete `currency.py:200-201`; it is the
sentence that prevents the measurement in RANK 2.

⛔⛔ **RANK 1 LANDS ONLY THE CORRECTIONS TO TEXT DESCRIBING SHIPPED CODE** — §9.1 rows 1–5, §9.3, §9.4,
§9.5, §9.6. **§3's derivation STAYS in `docs/dev/` until RANK 2 or RANK 4 lands a consumer.**
`docs/dev/README.md` sends a derivation to `EQUATIONS.md` when *"the code depends on"* it, and no code
depends on §3: `region_start_count` reaches `substrate.py:200` and stops. Writing it in now would make a
permanent doc describe code that does not exist — and §9.1's last row says so.

⚠ Two of the §9 items are TEST defects, not prose, and one of them is a dead gate — fix that one with a
perturbation, per `TRAPS: could-the-arm-have-fired`.

### RANK 2 — ✅ EXECUTED 2026-08-20, AND THE VERDICT IS A MEASURED NEGATIVE

The A/B is landed as `CalibrationConfig.region_abundance_bank` ("contained" shipped / "start"),
threaded `calibrate → build_region_geometry`, gated by an ABSOLUTE-truth fill gate written first,
verified failing, then fired by three perturbations (banks swapped; `ell+1` divisor; divisor
forgotten). `relay_pool_ab.py --region-bank` and `exon_solvability.py --region-bank` carry it, every
row stamped. Inertness verified on all 16 (`off` rows byte-identical across banks; base reproduces
session 2's table 48/48 by column name). ⚠ Latent instrument defect found: `relay_pool_ab.py`'s
report crashes (`StopIteration`) when `--arms` omits `off`, and the `--out` write sits AFTER the
report, so a single-arm run writes nothing.

**The verdict — currency arm, gdna |err| in FRAGMENTS, `start ÷ contained`, per stratum:**

| stratum | contained | start | ratio |
|---|---|---|---|
| unstranded × OFF ⭐ | 583,028 | 582,931 | 1.000× |
| stranded × OFF ⭐ | 635,702 | **536,176** | **0.843×** |
| stranded × ON ⭐ | 848,458 | 1,016,135 | 1.198× |
| DEFERRED unstr × ON | 7,890,982 | 14,537,395 | **1.842×** |
| the four zero controls | 71,826 | 156,857 | **2.184×** (each of the four regresses, 2.0–3.3×) |

⛔ **Fixing the level channel's real defect makes the losing policy WORSE nearly everywhere it was
best** — the zero controls and the capture rows. The per-condition shape names the suspect mechanism
(NOT yet dissected — it is a hypothesis with one measured leg): under `contained`, **58.6 % of live
exon hops had a REGION bank of exactly 0** and `enrichment_ratio` returned its 1.0 default — an
ACCIDENTAL MUTE — while under `start` the bank is alive at every slot, so the rescale machinery now
fires on hops it previously skipped, and the currency arm's zero-control standing partly RESTED on
the truncated bank's silence. The two capture-OFF g98 rows, where the truncation is largest and the
wall smallest, are the only clear wins (0.76–0.80×). ⭐ The channel defect (§2) and its repair remain
real and landed; what this prices is the POLICY's use of it — consistent with session 2's finding
that the currency gap is not one mechanism, and further evidence for §0c.0c that the policy's
delivered precision, not its inputs, is the binding problem.

⭐ **And the pass-0/deliverable SPLIT is measured on the worst in-scope condition**
(`exon_solvability --refits 0`, currency arm, `g98 ss0.99 ON`): at pass-0 the start bank IMPROVES the
exon solve — misplaced 149,736 → 138,136, and the NO-own-evidence bucket's ≤0.05 mass coverage goes
45.8 % → 66.1 % — while the FINAL whole-chain error on the same condition regresses 404,197 → 475,476.
The corrected level channel helps exactly where §6 predicted, and the propagation/refit machinery
converts the better input into a worse deliverable
(`TRAPS: the-intermediate-is-not-the-deliverable`'s shape, measured inside one arm).

#### the original RANK 2 spec, kept as the method record

1. Make the A/B a CONFIG VALUE rather than a diff (the precedent is `message_policy`; `DESIGN.md` §0c).
   ⛔ **It is not a two-line change, and the obvious path does not exist**: `build_region_geometry`
   (`region_geometry.py:198-206`) takes no config argument, and there is no
   `src/rigel/calibration/config.py` — `CalibrationConfig` is `src/rigel/config.py:247` and
   `message_policy: str = "relay"` is `src/rigel/config.py:450`. The shape that works:
   * add `region_abundance_bank: str = "contained"` to `CalibrationConfig`, beside `message_policy`;
   * add a KEYWORD-ONLY `region_abundance_bank="contained"` to `build_region_geometry` so the 15
     non-`calibrate` call sites are untouched (there are 16 in total — `calibrate.py:337`,
     `tests/calibration/test_region_geometry.py` x10, `tests/calibration/_synthetic.py:284`,
     `hop_currency.py:486`, `reframe_walk.py:111`, `calibration_oracle.py:189`,
     `object_composition.py:394`, `anchor_opportunity_census.py:132`);
   * thread `config.region_abundance_bank` from `calibrate.py:337`;
   * branch at `region_geometry.py:295-296`.
2. ⛔ **Write the falsification gate FIRST and watch it fail**, then break the fixed code three ways and
   watch each gate fire. **No test in the suite gates how `inv_abundance` is FILLED** — every occurrence
   is constructed synthetically (`test_currency_policy.py:81-113`, `test_sweep_backbone.py:36`,
   `test_density_model.py:102`), so this would otherwise land in an ungated line.
   ⭐⭐ **The gate must assert the ABSOLUTE truth, not flatness across two length sets.** Perturbing 14
   off-by-ones: an accuracy gate caught 14/14; a flatness gate caught only 10/14, missing every pure
   SCALE error (`1/(ell+1)`, predicate `s < ell+1`) — perfectly flat and perfectly wrong.
3. Score with `relay_pool_ab.py --arms off on currency` and
   `exon_solvability.py --arms relay currency --refits 0` — ⛔ **these are the only two instruments in
   the tree that can install `CurrencyPolicy`.** `solvability_audit.py` CANNOT: `:660` hard-codes
   `config = CalibrationConfig()` and its argparse has no policy flag, so scoring the swap there needs
   an `--arms` flag built first. `grep -rn message_policy scripts/` returns exactly
   `relay_pool_ab.py:116,173` and `exon_solvability.py:212-213`.
4. ⛔ **Per stratum, never pooled; the DEFERRED unstranded x capture-ON row reported as its own row;
   BOTH zero controls — which on the panel means the four `g00` rows of that same
   `relay_pool_ab.py` run.** ⚠ `zero_controls.py` is the TOY control: it reaches policy only through
   `toy_harness.with_messages`, which is `dataclasses.replace(config, message_propagation=...)`
   (`toy_harness.py:175-178`) — RelayPolicy vs SilentPolicy — and it runs toy specs rather than panel
   conditions. It cannot express this arm.

⚠⚠ **STATE THE CEILING UP FRONT.** `inv_abundance` has exactly ONE consumer in the tree —
`currency.py:204` inside `enrichment_ratio`, called once at `currency.py:445`. `RelayPolicy` and
`SilentPolicy` never read it, and `config.py:450` defaults to `relay`. **The swap cannot move the
shipped tool at all.** It can only move the `currency` arm, which `ROADMAP.md` §0 records as the worst
arm on every in-scope stratum. Land it as a DEFECT FIX INSIDE A LOSING POLICY, not as a 0.8.0 metric
move — and if it is presented as the latter, `TRAPS: an-ablation-that-never-ran` applies.

⛔ **It also changes the channel's POPULATION, and that must be argued rather than smuggled in.**
`region_start_count` includes SPLICED starts; the REGION slot today is unspliced-contained. Spliced
fragments are 0.61 / 15.5 / 20.3 / 31.4 % of the library at `g98 off / g50 off / g50 on / g00 off`. This
moves the REGION face TOWARD the BOUNDARY face, which is already spliced-inclusive via the hand-added
`inv_sj_lo`/`inv_sj_hi` (`currency.py:205-206`) — defensible, and a population change, not a bias fix.

### RANK 3 — ✅ EXECUTED 2026-08-20; §5 carries the table. What remains below is the method record.

§5's 27.3 % is genomic and is a lower bound for internal exons. Redo it against each REGION's
downstream SPLICED template length, weighted by `region_start_count`, on `g00` (the RNA-only control) and
on a contaminated row.

⛔⛔ **IT IS NOT A DROP-IN SUBSTITUTION, BECAUSE `region_start_count` HAS NO TRANSCRIPT AXIS.**
`substrate.py:149` is `int64[n_regions]` — one number per REGION — while `d_high` is per
`(REGION, transcript)`: an internal exon sits on many transcripts with different genomic-high ends and
different downstream spliced lengths. **Name the collapse rule before measuring: take the MAXIMUM
downstream spliced template over the transcripts containing the REGION, so the wall binds only if it
binds for ALL of them.** That is the conservative direction and it is the one that can refuse RANK 4.
⭐ The machinery exists: `transcript_truth.py` already computes interior cumulative-exon boundaries in
spliced-transcript coordinates, and `exon_solvability.py` already carries a region-to-transcript map.

⛔ **This number decides whether RANK 4 is a correctness requirement or an optimisation.** If it is
small, the start bank alone may be enough; if it is ~27 %, the mirror bank is mandatory and RANK 2's
result is not interpretable without it.

### RANK 4 — `region_end_count`, the mirror bank. One integer bank, one full re-scan. Only after RANK 3.

⭐ Treat this as the CONDITION OF CORRECTNESS for §3, not an enhancement.

* One `uint32[n_regions]` at `region_of_pos(last_base)`. ⚠ `last_base` is already computed at
  `accumulator.cpp:603`, but `region_of_pos(last_base)` is evaluated only inside the unspliced guard at
  `:719` — the mirror bank needs that lookup on SPLICED paths too.
* ⛔ **Trace the plumbing against `region_start_count`, NOT against
  `region_contained_inv_opportunity_sum`.** The latter is float64 and sits in `SINGLE_COLUMN_AXES`
  (`scan_payload.py:290-300`, *"every row here is a fraction, so every row is float64"*), so it is
  auto-plumbed by that table and an integer bank cannot follow it. The integer path is hand-plumbed:
  `ADDITIVE_AXES` (`scan_payload.py:312`), `:509`, `:614`, `:761-764`, `:814`, then `substrate.py:149`,
  `:200`. ⚠ `region_start_count` has 102 occurrences across 26 files and the float bank 33 across 15;
  the MANDATORY set is far smaller than either. Re-derive the count before quoting one.
* ⭐ INTEGER, so merges are bit-identical and it gets no `FLOAT_BANK_DEPOSITS` entry
  (`rescan_panels.py:78-85`) — it stays at byte-identity rather than a derived budget.
* CPU ~2–10 ns/fragment against a per-fragment scan budget in the high hundreds of ns. ⚠ **Not measured
  this session** — one reviewer timed 586 ns/fragment on `gdna_g50_ss_0.99_nrna_mid_capture_off`; re-time
  with `scan_profile.py` before quoting. Either way this is under 2 %.
* Memory +4 B/REGION: 141 KB on the suite index, ~4.5 MB/worker at human scale.
* ⛔ Land it WITH the executable specification (`_accumulator_reference.py`) in the same change; the C++
  is gated on byte-identity to it.
* RE-SCAN, **ESTIMATED — no cache was touched this session, so re-time with `scan_profile.py` before
  quoting any of these**: ~10 s/condition for the scan cache (16 conditions ~2.7 min), ~3.3 min/condition
  for the six oracle partitions (~53 min), ~2 min for the 16 test-chromosome scenarios — **about one
  hour.** The ~21 min simulate stage is NOT needed; the BAMs do not change.
* ⛔ **`panel.py cache --force --jobs 8` parallelises the SCAN caches ONLY, and `--force` does not reach
  the oracle layer.** `panel.py:196` passes `--force` to `build_scan_cache.py` alone;
  `pass0_vs_oracle.py` has no `--force` and picks its prewarm workers by file EXISTENCE (`:1044-1045`),
  so a bank ADDITION leaves every file in place, `todo` comes out empty, no worker spawns and the oracle
  stage runs SERIAL. Delete `<suite>/oracle_cache/*` by hand first, or drive that layer with
  `rescan_panels.py --jobs 8`.
* ⭐ **`rescan_panels.py` IS the right tool for a bank ADDITION and needs no change.** `:185` returns the
  symmetric difference, so a pure addition lands in `only_new` with `differing` empty, and two self-test
  cases already gate exactly that (`:399`, `:417`). ⚠ Its `PAYLOADS` is `("_main", *ORIGINS)` (`:61`) —
  four payloads — so `rna_pos` / `rna_neg` need one `panel.py cache` pass afterwards.
* ⛔ Re-run `preflight.py` and the instruments after the re-scan, not just the suite
  (`TRAPS: a-green-suite-hid-five-dead-instruments`).

### ⛔ REFUSED, with reasons — do not rebuild these

| candidate | why not |
|---|---|
| OVERLAP `1/(ell+w-1)` | Window straddles both flanking density steps, and the flank's OPPORTUNITY differs between the populations sharing the slot (`statics.free_pos`/`free_neg` gate RNA there; gDNA is always admitted), so it is not model-free at all. Localises at 0.30 of a 100x step at `ell = 98`. Needs a new float bank plus a companion count bank. |
| COVER / CROSS | Strictly between START and OVERLAP on both axes; beat neither. CROSS yields `0/0 = NaN` on the **508 REGIONs with `ell == 1`**, and `currency.py:218`'s `ok = (nbr >= 0) & (t_src > 0.0) & (t_dst > 0.0)` evaluates `NaN > 0.0` as False, **laundering the NaN into the 1.0 "no evidence" default** rather than crashing. |
| SPANS-STRICTLY, or CONTAINED+SPANS as a pair | Closes the truncation only to `1 - f(ell+1)` (0.9961–0.9981 measured); costs a full re-scan; placement-non-local (0.043 at 113x); reads 0.0912 in the tiling test. ⚠ The `region_spanning_*` banks were deleted at `5591cc01` on measured evidence — with the weight `1/L`, never `1/(w-ell-1)`, so the two halves of the identity have never been in the tree together. |
| `contained_inv / P_hat(w<=ell)` | Available exactly where it is not needed. `0/P_hat = 0` resurrects NONE of the dead exons; `P_hat` is the pooled mixture's CDF where the correct divisor is per component, re-importing the composition dependence; `P_hat = 0` for the 7,123 REGIONs with `ell < frag_min`. |
| BORROW — use the flanking BOUNDARY totals as the REGION's level | ⚠ Scores best for a DEGENERATE reason: the chain alternates, so `r` becomes boundary/boundary `== 1` by construction. It does not fix the level, it **switches the enrichment channel off** (`log r -> 0` so `w -> 0`, pure ABUNDANCE transport). Given `SilentPolicy` currently beats `CurrencyPolicy` on every in-scope stratum that may be the best OUTCOME — but land it under that name, as a channel mute, not as a fix. |

### The cache key — NOTHING needs renaming, and the received rule is stale

⛔ **`CLAUDE.md`'s *"THE SCHEMA DIGEST DOES NOT SEE A DEPOSIT-RULE CHANGE"* is no longer true.**
`scan_cache.py:220-224` hashes `deposit_digest()` into `payload_schema_digest()`, and `deposit_digest()`
(`scan_cache.py:227-298`) runs the NATIVE accumulator over a fixed 7-fragment fixture and hashes every
bank's raw bytes. Gated by `tests/test_scan_cache.py:302` and `:283`.

⭐ **Restate the rule as: the key covers `accumulator.cpp`'s DEPOSIT and not `resolve.cpp`'s FRAGMENT
CONSTRUCTION.** The 2026-08-19 "cached / skip on all 32" incident was the chimera repair (`386ad09e`),
which changed which fragments are OFFERED and no deposit rule — so `deposit_digest` was right not to
move and the caches were wrong to be reused. `--force` remains the tool for that class.
⚠ The `69a85be2` rename-that-moved-the-key was luck; it has since been replaced by a real key.

### 8.9 Scope — an OWNER CALL is needed before RANK 4

The 0.8.0 metric is the calibration result against ORACLE CALIBRATION, per stratum, and the
fragment-length COMPOSITION channel is retired until after 0.8.0. Reading this as:

* a **LEVEL** fix to an existing model-free channel that is measurably wrong — in scope, and RANK 2 is
  free;
* a **new bank plus a full re-scan** — that is an accumulator change in a calibration-focused phase, and
  §7.1's split hazard sits next to it.

⛔ Bring the owner the RANK 3 number and the §7.1 hazard together, and let them rule. Do not land RANK 4
on this note's authority — nothing in `docs/dev/` is authoritative.

---

## 9. THE MOVE LIST — EXECUTED 2026-08-20 (RANK 1), remaining rows only

⛔ The corrections to text describing SHIPPED code were MOVED to their permanent homes and their rows
deleted here per the MOVE rule: §9.1 rows 1–5 → `EQUATIONS.md` §1.2/§2/§3.5g; §9.2 → `TRAPS.md`
`a-cancellation-is-conditional-on-its-support` (150/150 re-derived); §9.3 → `DESIGN.md` §3.1-table /
§3.1-note / §3.7, `ROADMAP.md` §0's enrichment-ratio row, `CLAUDE.md`'s cache-key rule (restated:
covers the DEPOSIT, not fragment construction); §9.4 → the source and spec comments themselves
(including the seven `EQUATIONS.md` citations replaced with test/spec references); §9.5 → both gates in
`tests/native/test_fragment_length_proof.py` repaired/extended, each fired by TWO deposit-rule
perturbations (a `1/ell` rule and a perfectly-flat `0.5×` scale error) and re-greened; §9.6 → the
substrate rename (`PopulationView.inv_opportunity_sum` beside `inv_length_sum`, one rule per name).

### 9.1 — the one row deliberately NOT moved

| site | defect | correction |
|---|---|---|
| new | the whole of §3 above, plus `c(ell)` from §3.5 | ⛔ **NOT in RANK 1.** The truncation-free derivation belongs in `EQUATIONS.md` only once code depends on it — `docs/dev/README.md`'s criterion. Until RANK 2 or RANK 4 lands, it stays in this dev doc. ⚠ `scripts/design/region_density_derivation.py:20` is the only correct statement of the LIMIT on disk (T2), now cited from `EQUATIONS.md` §2 |

### 9.7 To `docs/TESTING.md` — the panel, the harness and the two gates (NOT in RANK 1)

⚠ `docs/dev/README.md` sends *"a panel, a harness, a gate"* there, and four things in this note are
currently unhomed: §10.3's two-condition fl-gap SIDE PANEL (⛔ never a ladder rung), §10.2's enumeration
harness, RANK 2 step 2's `inv_abundance` fill gate, and §9.5's repaired
`test_fragment_length_proof.py` gate (now landed; home its TESTING.md row when RANK 2's gate lands
beside it). Land each with its owner section.

## 10. WHAT WOULD FALSIFY THIS

1. **The truncation law.** `Sum(contained_inv) / Sum(region_start_count/ell)` per REGION class per `ell`
   band, against the condition's own `deposited_lengths` CDF. If the quotient is not ~1.0 off capture,
   the law in §2 is wrong. It read 0.93–1.03 on 11 of 12 live bins; the outlier is `exon`
   `ell` in `[50,150)` at 0.853 and the twelfth bin was never reported (§2). ⭐ **The equal-length ladder
   is the right substrate for this** — no new panel is needed.
2. **The start rule's unbiasedness.** The enumeration harness of §2: place one fragment at every
   admissible start of every length and require the estimate to equal the true `rho` ABSOLUTELY, for at
   least one length set entirely above `ell` and one entirely below. ⛔ Assert the absolute value, not
   flatness.
3. **The per-component half of §2 — the ladder is blind BY DESIGN and the substrate must be built.**
   With `E_g == E_r` there is nothing to see, and `flgap_short`/`flgap_long` were deleted 2026-08-13.
   ⭐ It is a CONFIG EDIT, not new code: `gdna_ladder.yaml` already carries independent
   `simulation.frag_*` (`:131-134`) and `gdna.frag_*` (`:156-159`) blocks, deliberately identical. **One
   two-condition SIDE PANEL with a real gap prices the whole thing.**
   ⛔⛔ It must be a SEPARATE side panel, never a ladder rung — the ladder's equal lengths are a forcing
   function (`TESTING.md` §0) and diluting them re-opens the bug class they exist to expose
   (`TRAPS: equal-lengths-carry-no-composition`).
4. **The wall exposure** (§8 RANK 3), in transcript coordinates. If it is near zero, RANK 4 is optional.

---

## 11. WHAT WAS NOT DONE

* No `src/` change, no test change, no commit, no cache touched. `git status` is unchanged apart from
  this file.
* The RANK 2 A/B was NOT run. §6's per-hop numbers are a read of the shipped banks against truth, not a
  re-solve, so they say the ratio is wrong — not what fixing it is worth.
* §5's wall exposure is genomic and is a lower bound; §8 RANK 3 exists to replace it.
* The per-component truncation ratios in §2 are computed from ASSUMED lognormal pmfs, not from a real
  library. Treat the SHAPE as established and the multiplier as illustrative.
* ⛔ No ceiling was measured. If one is wanted, read the ceiling caution in `CLAUDE.md` first — every
  measurement arm patches `assemble_priors`, and the effective-length shrinkage sits outside it.
* ⛔ **The evidence trail must be recoverable from the REPOSITORY, not from a session scratchpad.** Both
  primary sources resolve by git ref:

      git show 8c673324:docs/NODE_DENSITY_DERIVATION.md
      git show a5ed752b^:docs/dev/counts_densities_paradigm.md

  Every enumeration harness behind §3.2, §5, §7.2 and §9.5 is ~30 lines and is re-derivable from the
  formulas as written; ⚠ any claim here that CANNOT be re-derived that way should be treated as
  unverified until it is.
