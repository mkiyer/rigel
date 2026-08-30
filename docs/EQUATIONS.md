# EQUATIONS — the derivations the code depends on

**Every formula here is implemented somewhere and gated by a test.** This file exists so the derivation
does not have to be reconstructed from the code, and so two modules cannot come to disagree about one
quantity (TRAPS: two-docstrings-one-quantity).

⚠ **ONE BLOCK IS AN EXCEPTION AND IT IS BANNER-HEADED AS ONE.** §3d–§3e derive the fragment-length
**composition** channel, which is **DEFERRED POST-0.8.0** (owner, 2026-08-14) and is not in the tree —
nor is the §3.1 2×2 it rests on solved anywhere in calibration. They are kept because a derivation is a
record and it will be wanted after the release. ⛔ Nothing else here is optional in that way, and the
deferral does **not** reach the length **model**: the opportunities, effective lengths and length moments
of §1, §3.6b and §4 are live and shipping, and so is the second pass's per-fragment length term (§10).

⚠ **A SECOND EXCEPTION, AND IT IS THE OPPOSITE KIND — §9d is a derivation AHEAD of the tree rather than
behind it.** The capture reference's spike-and-slab is derived, its three limits are checked and every one
of its parameters is measurable off the shipped payload today — but **nothing in `src/` builds the
mixture**, and what ships is §9c/§9c.1's Beta with the structural per-object location. It is banner-headed
as such, and the numbers inside it measure the BOUNDS and the gDNA FIELD, never a shipped estimator.

⛔ **No measurements here.** Where a number appears it is a property of the *formula* (a limit, a
degenerate case, a worked value that makes the shape concrete), never a property of a library.

⚠ **THAT RULE IS THE INTENT AND THE FILE DOES NOT FULLY KEEP IT.** Several derivations below are anchored
by the measurement that verified them, and those stay — a derivation with no verification is worth less,
and deleting one to satisfy the rule loses the finding rather than the untidiness. ⛔ What is required of
them instead is a **PANEL**: a number describing a library must name what it was measured on.

⛔⛔ **STAMP, APPLYING TO THE WHOLE FILE: THE 36-CONDITION LADDER AND THE `flgap_short` / `flgap_long`
PAIR WERE RETIRED ON 2026-08-13.** The ladder was rebuilt from scratch at **16** conditions
(`g00/g05/g50/g98` × ss `0.50/0.99` × capture off/on) and both `flgap` panels were deleted from disk —
only their configs survive, and they are a PAIR that works only as one. So every number below attributed
to *"the 36-condition ladder"*, to a `g01` / `g10` / `g25` / `g75` / `g90` donor, or to a `flgap`
condition was measured on a panel that **no longer exists**, and **none of them has been re-derived on
the panel now on disk**. ⭐ They are kept because a measurement is a record. **Re-derive rather than
trust — a number that has moved is a result, not a documentation bug.**

---

## 1. The deposit rule

**Notation.** A fragment is an interval `[s, s+w)` of *molecule* length `w`. A reference is region bound into
**regions** at every exon endpoint of every non-synthetic transcript; the 0-bp boundaries between adjacent
regions are **boundaries** (contiguous boundaries). `ell` is a region length.

**1.1 Fragment length `L` — ONE definition.** `L` = genomic span **minus region bound introns**. A paired-end mate
gap counts toward `L`; an intron does not. Whatever counts toward `L` must also count as coverage for
crossing, or the estimator is biased.

**1.2 Crossing count — the core identity.** A fragment `[s, s+w)` crosses an inter-base point `p` iff
`s < p < s+w`, i.e. `s ∈ [p−w+1, p−1]` — exactly `w − 1` start positions (the shipped rule,
`accumulator.cpp`'s `1/(length−1)` deposit is its reciprocal). So

    E[count at p]  =  ρ · Σ_w f(w)·(w−1)  =  ρ · (mean_FL − 1)

exactly, for any fragment-length distribution, and **independent of both flanking region sizes**. This is
why the partitioning problem dissolves: a count at a 0-bp boundary depends on nothing but the density and the
length distribution. ⚠ §1.4's formula and §3c already carry the `w − 1` form; this section was the odd
one out (corrected 2026-08-20 — the `w`-position version was never true of shipped code).

**1.3 Contained opportunity.** The starts at which a length-`w` fragment fits wholly inside a region of
length `ell`:

    (ell − w + 1)₊

and over a set of regions, `A(w) = Σ_n (ell_n − w + 1)₊`. Its fragment-length expectation is
`E_f[(L−w+1)₊] = (L+1)·F(L) − S(L)` with `F` the FL CDF and `S(L) = Σ_{w≤L} w·f(w)`; beyond the support it
is `L + 1 − mean_FL`.
⛔ **The `+1` is the discrete count of start positions, not a fudge** — drop it and the divisor is exactly
0 when a region is one fragment long.

**1.4 Crossing-exactly-one-line opportunity.** For a boundary with flanking region lengths `a` (left) and `b`
(right):

    A(w)  =  (w−1)₊  −  (w−1−a)₊  −  (w−1−b)₊  +  (w−1−a−b)₊

⭐ **The two nearest boundaries are the only ones that need excluding**, so this is exact rather than a
truncation: a fragment is an interval containing the boundary, so if it reaches any boundary beyond `p−a` it must
also cross `p−a`. ⚠ Reference ends need no special case — the partition region bounds at `0` and `L_ref`, so the
outermost region's length *is* the distance to the wall.

**1.5 General crossing divisor — ONE formula for both boundary kinds.** With `R_lo`, `R_hi` the molecule's own
remaining sequence either side:

    E_J  =  E_f[ min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)₊ ]

Mean fragment length is its **large-reach limit**, not a separate case. On RNA N(200,50): 199.0 at R=550
(so mid-transcript sj are exact), 160.1 at 200, 87.8 at 147, 19.6 at 100, 50.0 at R=50 — a **4×
error** if mean length is used blindly at a first exon.

**1.6 Reach.** At position `p` inside exon `e`: `reach_lo = exonic_bases_before(e) + (p − e.start)`,
`reach_hi = total_exonic − reach_lo`; maximised over transcripts independently per side AND per strand.
gDNA is unbounded on a chromosome; mature RNA stops at the polyA site. **A reach of 0 is meaningful, not a
sentinel.**

**1.7 Partition rule.** `region bounds = unique({exon.start, exon.end over all non-synthetic transcripts} ∪
{0, ref_length})`, region `i` = `[region bounds[i], region bounds[i+1])`, **no merging**. Introns and termini are already exon
endpoints, so only the merge step was removed. Boundaries always run `src < dst`, so **genomic order is a
topological order** and there is no graph traversal anywhere.

---

## 2. Reciprocal opportunity, and where it is model-free

⭐ **The general rule (shipped 2026-08-10, `69a85be2`): deposit `1/A(w)` where `A(w)` is that
population's own OPPORTUNITY** — the number of admissible start positions for a length-`w` fragment at
that object. Then, by linearity alone,

    E[Σ 1/A]  =  Σ_c ρ_c · P_c(A > 0)

⛔⛔ **THE CANCELLATION IS CONDITIONAL ON ITS OWN SUPPORT** (`TRAPS:
a-cancellation-is-conditional-on-its-support`): a fragment with `A(w) = 0` does not deposit a small
number, it deposits **nothing**, so the factor `P(A > 0)` survives — and it is a functional of exactly
the fragment-length pmf the channel claims independence from, per component. The executable statement is
`scripts/design/region_density_derivation.py`'s T2, perturbed.

**2.1 At a BOUNDARY it is model-free for every real library.** The crossing opportunity is
`A(w) = (w−1)₊` (§1.2), the shipped deposit `1/(L−1)`:

    E[Σ 1/(w−1)]  =  Σ_c ρ_c · P_c(w ≥ 2)  =  Σ_c ρ_c      whenever frag_min ≥ 2

so the support factor is exactly 1 for any library whose fragments are at least 2 bp — unconditionally
robust in practice, across a mixture with different component lengths. The same holds for the sj bank,
which takes the identical deposit.

**2.2 At a REGION it is model-free only WITHIN its support — a density SHAPE, not a TOTAL.** The
contained opportunity is `A(w) = (ell − w + 1)₊`, the shipped deposit `1/(ell − w + 1)`
(`region_contained_inv_opportunity_sum`), and

    E[Σ 1/A]  =  ρ · P(w ≤ ell)          ← NOT ρ

`P(w ≤ ell)` is a pmf functional and **differs per component**, so at REGION scale the gDNA-vs-RNA
circularity this channel exists to remove moves out of the divisor and into the support. Measured on the
ladder's own lengths: an **11.6× under-read at a 98 bp exon** (the panel's median exon), exactly **zero
at any depth** below `ell < frag_min = 50` (7,123 of 24,018 exon REGIONs). ⛔ **Do not read the REGION
bank as a model-free TOTAL**; within a fixed `ell` band it is a valid density SHAPE, and the truncation
law `Σ(contained_inv) / Σ(start_count/ell) ≈ P(w ≤ ell)` is confirmed on the shipped banks (0.93–1.03 on
11 of 12 live capture-OFF bins, `g98 ss0.99`).

**2.3 And it fails at a terminus, exactly.** `E[Σ1/L] = ρ · E_f[placements(w)/w]`, which equals `ρ` only
where placements ∝ `w`. At a point 50 bases from an end, placements = 50 for **every** `w > 51`,
independent of `w`. So `Σ1/L` fixes the length bias and a reach taper fixes the placement loss; **neither
substitutes for the other.**

**2.3b ⭐⭐⭐ THE START/END RELATION IS MODEL-FREE AT A REGION, AND THE WALL IS THE ONLY EXCEPTION.**
A fragment's FIRST covered base falls in region `r` iff its start lies in `r`, so the opportunity is
`ℓ` — the same number for every fragment length, which is exactly what the contained relation is not:

    E[S_r]  =  ρ · ℓ            for every w        ⟸ no pmf functional, no support factor
    E[C_r^inv]  =  ρ · P(w ≤ ℓ)                    ⟸ §2.2, truncated per component

so `S_r/ℓ` is a TOTAL where the contained bank is a SHAPE. The one exception is the template WALL:
within `d` of the template's genomic-HIGH end, only starts that still admit a length-`w` molecule
count, so

    A_start(w | d)  =  min( ℓ , (d + ℓ − w + 1)₊ )              ← w-dependent again
    A_start(w | d)  =  ℓ    for every w    ⟺    d ≥ w_max − 1   ← the exactness condition

and at `d = 0` this is `(ℓ − w + 1)₊`, the CONTAINED opportunity: **at a flush wall the start relation
degenerates into the relation it replaces while still depositing the flat `1/ℓ`, and reads worse than
it.** `E_r` is the mirror, exact iff `d_low ≥ w_max − 1`. ⭐ The two therefore fail at OPPOSITE ends,
which is what makes the pair closed: **use the side whose wall does not bind; average where both are
exact** (two counts of one rate at one opportunity, so the pooled rate IS the precision-weighted
combination). Where both bind, no deposit rule is model-free and the honest output is a refusal.
⭐ **The wall is COMPONENT-DIFFERENTIAL, and that is the substantive content**: gDNA's template is the
chromosome and never binds, a nascent molecule's is its genomic span, a mature molecule's is its
SPLICED length — so the distance must be taken at the component MINIMUM over the populations §0's
`T(slot)` admits, and a genomic distance at an exon marks a binding mature wall exact.

⭐⭐ **AND THE PAIR GIVES A FIELD-FREE TEST OF ALL OF THIS, which no comparison against the contained
bank can.** `S_r/ℓ` and `E_r/ℓ` share their opportunity and read the same field at the same slot, so
their ratio has an expectation of exactly 1 wherever both sides are exact — with no uniformity
assumption and no fragment-length distribution entering. Where one side's wall binds, that side is
depressed and the ratio moves in a KNOWN direction. Measured on the 16-condition ladder, every stratum
and both zero controls: both-exact **1.0018–1.0057**, start-exact/end-bound **0.712–0.879** (`< 1`),
end-exact/start-bound **1.242–1.428** (`> 1`) — capture-ON included, where every field-dependent
⚠ **Re-measured on the SPARSE-nascent ladder (2026-08-22): 1.0000–1.0074 / 0.6974–0.9955 / 1.0216–1.4878, 100 gates 0 failures — the DIRECTION holds on every row, and the bound-side margin has nearly collapsed at `g98` capture-OFF (0.9952).**
(the numbers above are the PRE-REBUILD panel's)
comparison is vacuous. ⛔ A comparison against the contained bank is NOT field-free and must not be
read as one: the two weight a region's positions differently, so a probe-shaped field separates them
legitimately (`TRAPS: two-estimators-of-one-rate-weight-the-field-differently`).

**2.4 Superadditivity.** Contained effective lengths are **superadditive**: over 118,195 splitting regions
`Σ E(children)/E(whole)` = 0.7652, and 0.0917 for 305 bp → 145+160. So densities, effective lengths and
variances **cannot be pooled across two partitions**; only truth pools additively
(`f_g = ΣG/(ΣG+ΣR)`).

---

## 3. The two-component deconvolution

**3.1 The 2×2 at one object.**

    N       =  ρ_g·E_g[w]  +  ρ_r·E_r[w]
    Σ 1/L   =  ρ_g         +  ρ_r

Identified **iff the two MEAN lengths differ**; the second row being literally `[1,1]` is what makes the
2×2 well conditioned. `N / Σ(1/L)` is the abundance-weighted mean fragment length.
⛔ **So the identifying quantity is the GAP `μ_g − μ_r`, not either mean.** At equal means the channel
carries exactly zero information at any depth (TRAPS: equal-lengths-carry-no-composition).

⛔⛔ **AND THE CHANNEL THIS 2×2 DESCRIBES IS DEFERRED POST-0.8.0** (owner, 2026-08-14). Fragment length as
a CALIBRATION composition channel is out of scope for the release, so the 2×2 above stands as a
derivation and **no calibration module solves it**: outside `substrate.py` nothing in `calibration/` reads
an `inv_length_sum` bank at all (⚠ `sj_inv_length_sum` *is* live — in the SECOND PASS, §10, which is a
different question). ⛔ Do not propose it, do not rank it — §3d–§3e carry the banner, the verdict and the
mechanism.

⭐⭐ **AND THIS IS WHY THE STAGE-B LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS — the reason is
STRONGER than the one that used to be written down.** It is not "the length channel is neutralised, so the
residual is attributable to density and strand". It is that **the EM ALREADY USES THE FL DISTRIBUTION**: a
large gDNA-vs-RNA length gap lets the EM assign fragments on LENGTH ALONE, **bypassing calibration
entirely and MASKING its bugs**. Equal lengths FORCE the calibration phase to be exercised, and what
calibration then has to work with is strand and density (plus propagation across objects, currently off).
⚠ So equal lengths are a property of the INSTRUMENT and never a claim about real libraries; §3b records
what that choice costs, and why a gapped panel is where that cost is measured instead.

**3.2 Density is the frame-invariant currency; a fraction is not.** `ρ_c = C_c/E_c` agrees across the
contained / crossing / spliced frames to ~0.36 % and does not degrade across 1×/30×/300× capture. The
log-odds shifts between frames by exactly `log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`, and **capture
cancels identically** — a ratio transports across a capture cliff, an absolute density does not.

**3.3 Conservation with unequal effective lengths is `Σ_c ρ_c·E_c = M`**, not `Σ_c ρ_c = M/E`. Its
sensitivity is bounded: a purely compositional error moves `M/Σρ_c E_c` by only ×1.04 on a contained region
and ×1.50 at a crossing, so a large violation is accumulated drift, never one hop.

**3.4 Why an integer count must be stored.** `Var(log ρ_c) = 1/(f_c·n) ≡ Var(log f_c) + 1/n`, exactly.
Mass sums fractional per-fragment shares, so `1/mass` is not a counting variance.

**3.5 ⭐⭐ The reframe is a COMPOSITION IMPUTATION, and a gDNA LEVEL crosses unscaled.** Substituting
`ρ_c(src) = φ_c(src)·ρ_tot(src)` into the reframe `r = ρ_tot(dst)/ρ_tot(src)` gives what it actually
delivers:

    ρ_c^msg(dst)  =  ρ_c(src)·r  =  φ_c(src) · ρ_tot(dst)

— the source's density **share** applied to the destination's observed total. There is **no level transport
in it at all**: it is exact iff `φ_c(src) = φ_c(dst)`, and wrong by exactly `φ_c(src)/φ_c(dst)`.
⛔ When the source carries only gDNA that factor is `1/φ_g(dst)` and the delivered level collapses to
`ρ_tot(dst)`, the destination's **own** total, independently of the source's measurement — TRAPS: a-message-from-the-destinations-belief/TRAPS: a-total-density-ratio.

⭐ So a message makes two claims and they need two scales:

| claim | currency | scale | licence |
|---|---|---|---|
| COMPOSITION (`λ`, `τ`) | a share | `r` | ⭐ **§3.5b, two conjuncts** |
| a component's LEVEL | an absolute rate | its own landscape's ratio | — |

**3.5b ⭐⭐⭐ THE LICENCE — "is the source measuring the same thing I am?"** (owner, 2026-08-04.) A
composition may be imputed across a step iff **both** hold:

| | conjunct | why |
|---|---|---|
| **SUPPLY** | the source **SUPPLIED both components** of the pair — the λ-emission gate's own predicate, a statement about **precision**, not density | a source carrying one component has no share to lend; λ is not "large" for it, it is UNDEFINED |
| **POPULATION** | the two objects measure the **same RNA population** | otherwise enrichment and a population difference are indistinguishable, and the imputation reads one as the other |

⭐ The population conjunct is set algebra, not a model. An BOUNDARY counts what spans it **contiguously**, so

    T(BOUNDARY)  =  T(NODE_left) ∩ T(NODE_right)

and `T(BOUNDARY) = T(right)` fails iff a transcript's body **begins** at the BOUNDARY, `T(BOUNDARY) = T(left)` iff one
**ends** there. So a transcript **terminus** is exactly what makes one flank's population larger, the test
is an **equality** per `(BOUNDARY, side)` pair (both directions of mismatch corrupt `φ_g`), and it is per-step
rather than per-object. `region_geometry.terminus_flank_gain`.

⛔⛔ **WRITE IT IN GENOMIC TERMS, NEVER IN TSS/TES.** TSS/TES is transcript-relative and the strand flips
it:

    the RIGHT flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic LOW  end is there  ⟺  TSS₊ or TES₋
    the LEFT  flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic HIGH end is there  ⟺  TES₊ or TSS₋

⛔ **TERMINI ONLY.** A DONOR/ACCEPTOR BOUNDARY also changes the population, but there the flux is **measured**
(`junction_count`) and the SPLICE IN and the SPLICE OUT exist to route it. A terminus has no flux to measure: a
transcript simply begins. That is the derived boundary between the two treatments.

**3.5c ⭐⭐ THE MASS PIN CARRIES THE SAME LICENCE, PLUS ONE STRUCTURAL CASE.** `Σ_c ρ_c·E_c = M` is an
identity *under the imputation premise*, restored by `k = M/S` with `S` filling every unsupplied component
from the destination's **own** density. TRAPS: a-message-from-the-destinations-belief permits a message to use the destination's constants and
observations, never its beliefs — so the pin is licensed in exactly the two states where no belief reaches
`S`: **(i)** the message supplied the composition (§3.5b, nothing is filled in), or **(ii)** the
destination is a **structurally pure-gDNA object**, where there is no unsupplied component and `f_g = 1` is
structure, so `S = ρ_g·E_g` and the pin hands the object its own **measured** `M/E_g`.
⭐ Case (ii) is how the capture landscape reaches an exon at all — a G1 BOUNDARY has `prec_g = 0` and cannot
originate a level through the fuse, so the pin is its only channel. Unlicensed, the pin was TRAPS: a-message-from-the-destinations-belief at full
strength: `k = 1/(φ_msg + R_own)`, fixed point `(1−R_own)·ρ_tot`, and `R_own` is **exactly ½** at a slot
with no composition evidence — so it drove the delivered gDNA *fraction* to ½ regardless of the truth.

and for **gDNA the level crosses UNSCALED**, because gDNA is uniform along the genome *before* capture.
⭐⭐ **This is one rule for both capture arms with no branch on capture**, because the capture landscape is
carried by the structurally-pure-gDNA population's **own measurements**, not by a scale factor: the relay's
mass pin re-anchors the running level at every object, and at a pure-gDNA object that total *is* its gDNA
density at its own capture stratum. A per-capture-class landscape ratio built to do this explicitly measured
**byte-identical off capture** and 1.2 % of one class on one capture-ON condition, and was deleted.

⚠ **One residual is separate and is not this**, so a reader does not attribute it here: the level an exon
receives is the one measured at its flanking BOUNDARY, and a fragment spanning a BOUNDARY at a gene end lies only
PARTLY under the capture probe while one contained in the exon lies wholly under it — so that level is a
**lower bound** under capture, and closing the gap needs a probe-geometry model the tool does not have.
⭐ The second residual that used to be listed here — the mass pin filling unsupplied components from the
destination's own `f_g ≈ ½` — is **fixed**: §3.5c.

**3.5d ⭐⭐⭐ WHY `r` FROM TOTALS IS NOT A COMPONENT'S OWN RATIO — the share-weighted-average identity.**
(Owner's question, 2026-08-05; measured by `reframe_walk.py` §3b.)

§3.5 says the reframe delivers `φ_c(src)·ρ_tot(dst)` and is exact iff `φ_c(src) = φ_c(dst)`. This is the
same statement read from the other end, and it is the one that says **what to do**:

    r_tot  =  [ρ_g(dst) + ρ_R(dst)] / [ρ_g(src) + ρ_R(src)]  =  φ_g(src)·r_g  +  φ_R(src)·r_R

⭐⭐ **The ratio of TOTALS is the share-weighted average of the components' OWN density ratios, weighted by
each component's share of the SOURCE's mass.** Verified to floating point on every hop where both ratios
exist. Two immediate consequences:

| `φ_g(src)` | what `r_tot` is | is scaling the gDNA arm by it valid? |
|---|---|---|
| **1** (intron, intergenic REGION, a gDNA-only BOUNDARY flank) | `r_tot ≡ r_g`, identically | ⭐ **yes** — nothing foreign enters |
| **≈ 0.0006** (an expressed exon) | `0.9994·r_R + 0.0006·r_g` | ⛔ **no** — it is the RNA ratio wearing a total's name |

Measured, `spliced_exons` × `g75 ss0.50 capture_off`, 200 k RNA: at `exon → intron|exon @2,000`,
`r_g = 1.3545`, `r_R = 1.1026`, `r_tot = 1.1027`. On the eight hops whose both ends are pure gDNA,
`r_tot` and `r_g` agree to 4 dp.

⛔⛔ **AND THE TWO COMPONENTS FAIL IN DIFFERENT WAYS, WHICH IS WHY DEPTH CANNOT SETTLE THIS.** With capture
off every one of these ratios should be exactly **1.0** — the true enrichment is 1, gDNA is uniform along
the genome, mature RNA is uniform along its transcript. Standard errors from 1.0, same run:

* the **gDNA** arm at a BOUNDARY: **0.6–1.8 se**. That is NOISE, and depth shrinks it. An BOUNDARY is 0 bp, so its
  opportunity is ~one mean fragment length whatever the chromosome does (`E_g = 215.7`, `E_J = 202.7`) —
  measured, the 7,000 bp intron REGION holds **518** gDNA counts while the two BOUNDARIES hold **19** and **24**.
  ⚠ Only library DEPTH helps; lengthening the reference does nothing off capture (`TESTING.md` §0b).
* the **RNA** arm across a sj: **13.5–13.7 se**. That is a BIAS. Depth makes it MORE significant, not
  smaller. It is §3.6b's frame gap: `E_J = E[w]−1` rises with the fitted mean fragment length while
  `E_r = e−E[w]+1` falls, so a length-model error is amplified with opposite sign at 0.62 %/bp.
  ⛔ TRAPS: a-variance-cannot-fix-a-bias — a variance cannot fix a biased mode.

⭐ **So the existing `σ²_transfer = Var(log r)` is the right medicine for the gDNA arm and the wrong
medicine for this**: `r_tot` is not a *noisy* estimate of `r_g`, it is a *precise* estimate of a different
quantity. Damping it treats a targeting error as an uncertainty.

⛔⛔ **THE EXTREME CASE, where the identity cannot be written at all.** If the source has ZERO of a
component the destination has plenty of, `r_R` is undefined, `φ_R(src) = 0`, and the product is `0·∞`:
`r_tot` then contains a term with **no counterpart at the source**. Measured at
`intergenic|exon @1,000 → exon`: `r_g = 0.7014` while `r_tot = 1187.87` — **1,694× apart**. ⭐ This is
exactly what §3.5b's SUPPLY conjunct catches, and it does: `lend` is False there and the gDNA level crosses
unscaled. The licence is right; what §3.5d adds is that between that extreme and `φ_g = 1` there is a
**continuum**, and the licence is a boolean across it.

⭐⭐⭐ **WHAT THE IDENTITY LICENSES, since the error is now a closed form.** The error of `r_tot` as an
estimator of `r_g` is

    r_tot − r_g  =  φ_R(src) · (r_R − r_g)

and every term on the right is measurable from the solver's own state. So the choice is no longer
boolean-or-nothing: `r_g` is directly available and unbiased-but-noisy, `r_tot` is precise-but-biased with a
**computable** bias, and two estimators with opposite failure modes fuse by inverse variance — the same
shape §3.6 already uses for its two closures and the-residual-level for the level. ⛔ No threshold, which matters because
a threshold here has been refused three times (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace). ⚠ Unbuilt and unpriced: the
ceiling to measure first is what a PERFECT per-component `r_g` at every hop is worth, re-solved.

**3.5e ⭐⭐⭐ THE TWO OPERATORS AND THE TERMINUS, WORKED IN NUMBERS — owner's derivation, 2026-08-19.**
(The ruling is `DESIGN.md` §0c.0. ⚠ The policy this was built against was DELETED 2026-08-27 — the derivation is kept so it is not re-derived, not a description of shipped behaviour, and the
toy rungs are checked against it. A message is the three DENSITIES `{gDNA, RNA+, RNA−}`, always.)

⭐ **(i) SPLICE OUT — exon REGION → `exon|intron` BOUNDARY.** Single-stranded, + strand. The exon's
contained belief is `{gDNA = 10, RNA+ = 90, RNA− = 0}`; the boundary has a measured sj flux `SJ+ = 17`
and an unspliced belief `{gDNA = 2, RNA+ = 1, RNA− = 0}`.

    1. the exon sends {10, 90, 0} — its densities, as they are
    2. the boundary has a + strand sj, so it knows it must SPLICE OUT
    3. the boundary's total INCLUDING the sj is 2 + (1 + 17) + 0 = 20; the exon's is 100; so the
       boundary RESCALES itself by 100/20 = 5:  {2, 18, 0} -> {10, 90, 0}, and SJ+ 17 -> 85
    4. now it SUBTRACTS the sj density: RNA+ 90 − 85 = 5 unspliced, still at the exon's scale: {10, 5, 0}
    5. with the sj removed it UNDOES the rescale: {10/5, 5/5, 0} = {2, 1, 0}
    6. that is the unspliced message the boundary compares with its own {2, 1, 0}

The rescale is by the **ratio of totals with the sj flux counted in the boundary's total** (§3.6c's
per-flank total), the sj is removed at the exon's scale, and the rescale is undone. ⭐ The exon and the
boundary are *composition-compatible* — the same transcript population, RNA+ partly leaving by the sj —
and the rescale factor absorbs whatever common-mode difference in opportunity or enrichment the two
totals carry. ⛔ Nothing about it branches on capture.

⭐ **(ii) SPLICE IN — intron REGION → `intron|exon` BOUNDARY → exon REGION.** The inverse.

    1. the boundary sends its densities INCLUDING the spliced-in flux, because those fragments continue
       contiguously into the adjacent exon:  {gDNA, unspliced RNA+ + SJ, RNA−}
    2. the exon sees a message of total `rho_tot_src` and has its own total `rho_tot_dst`; the rescale
       factor is `rho_tot_dst / rho_tot_src`
    3. the exon rescales the message (they are compatible) and compares it with its own belief
       {gDNA, RNA+, RNA−}

⛔⛔ **(iii) THE TERMINUS — where the arithmetic above is NOT licensed.** `TA+` exons (1000, 5000) and
`TB+` exons (2000, 5000); the boundary at 2000 is `TB`'s TSS and, by flank, `exon|exon`. An unspliced
fragment crossing 2000 is compatible with `TA+` and with gDNA, **not with `TB+`**. So the boundary at 2000
does not measure `TB+`, and the transcript population in the REGION (2000, 5000) is **not** the population
at the boundary or in (1000, 2000). That is what the terminus bit says: **the population changed, the
rescale arithmetic of (i)/(ii) cannot account for a transcript that ORIGINATES after the boundary, and
only the LEVEL can be propagated** from the boundary at 2000 into (2000, 5000). ⛔ **A terminus that sits
AT a splice junction is still a terminus** (owner, 2026-08-19): the population still changes across it
and the rescale is still incorrect — `sj+term` is ruled with `term`, never with `sj`.

⚠ **Capture, in the same example.** The probes could be anywhere along (1000, 5000): one at (1100, 1220)
enriches `TA+` and not `TB+`; one at (3000, 3120) enriches both — and the gDNA NEAR it, while the gDNA
crossing position 2000 is not enriched by it. So the LEVEL that crosses a terminus under capture is the
boundary's level, and the region's may differ from it by a per-exon enrichment nobody has measured.
⭐ Measured on the rebuilt panel (`hop_currency.py`, 2026-08-19, strands pooled): the LEVEL from a
gene-edge or `exon|intron` terminus under-reads an exon's gDNA by 34–37 % of the exon's mass at `g50`
capture-ON, from an `exon|exon` terminus by 7–9 %, and by ≤ 1 % off capture. §9d.3's two bounds (the
flank's level as a floor, the total minus the sj-measured RNA as a ceiling) are what bracket that level.
⛔ Nothing here is built; the toy rungs meet this case at capture-ON, with ground truth.

⭐⭐ **(iv) THE THREE ARMS ARE NOT OPTIONAL** (owner, 2026-08-19): every one of the steps above is
written in `{gDNA, RNA+, RNA−}`; a map or an operator that pools the two RNA arms has not measured the
thing the policy carries.

**3.5f ⭐⭐⭐ THE TRANSPORT KNOB — the two strategies as ONE shrinkage estimator, with no constant.**
(The ruling is `DESIGN.md` §0c.0b. Derived and gated 2026-08-20.)

A hop's two strategies are two point hypotheses about ONE unknown, the log enrichment `eta` between the
two objects:

    ABUNDANCE     eta = 0                 transport the abundances as they are
    COMPOSITION   eta = log r             transport them rescaled by the observed ratio of totals

Neither is a mechanism to select; both are estimates of `eta`, and they fuse. The observed `log r` carries
the counting variance `v` of the two abundance observations; the no-enrichment premise, if it is wrong, is
wrong by `log r` itself, so its own squared error is `(log r)²`. Inverse-variance fusion of the two gives

    w      =  (log r)² / ( (log r)² + v )
    eta^   =  w · log r                       Var(eta^)  =  w · v
    claim  =  rho · r^w

⭐ **Every quantity is measured**: `r` from the model-free reciprocal-opportunity banks (§3.5g) and `v` from
`count_logvar` of the two slots' counts. **No constant is introduced**, and the two ends are recovered
exactly — `w → 0` where the disagreement is indistinguishable from counting noise, `w → 1` where it dwarfs
it. ⭐ The unabsorbed residual `(1 − w)·log r` is the ABUNDANCE strategy's transfer variance
`((1−w)·log r)²`, which is the full squared disagreement at `w = 0` and vanishes at `w = 1`: one expression
covers the continuum, which is what makes this a derivation rather than a switch.

⛔⛔ **AND THE FACTOR MAY NOT BE DERIVED FROM THE SOURCE'S OWN CLAIM.** The mass-identity form
`k = M_dst / Σ_c ρ_c,src·E_c,dst` is exact under the premise and reproduces §3.5e's worked numbers in both
directions — and is unbounded, because the source's claim is its denominator: measured, a source holding
`ρ_g = 3.9e-4` with no RNA returned `k = 235,800` (`TRAPS: a-rescale-that-reads-the-source-belief-is-unbounded`).
The measured ratio agrees with `k` wherever the source's belief is right and cannot amplify a claim the
source does not have.

**3.5g ⭐⭐⭐ A TOTAL ABUNDANCE MUST NOT BE `mass / effective_length` — the accumulator already deposits
the composition-free quantity.** (Owner's question, 2026-08-20: *"how are we trusting our enrichment
ratios now?"*)

An effective length is a function of the FRAGMENT-LENGTH distribution, and gDNA and RNA have different
ones — so `mass / E` is a function of the composition being solved for and any enrichment ratio built from
it is circular. The owner's own arithmetic: 100 counts in a 500 bp region reads `100/(500−100) = 0.25` as
pure gDNA and `100/(500−200) = 0.33` as pure RNA.

⭐⭐ **The accumulator deposits the RECIPROCAL OPPORTUNITY per fragment** — `1/(w−1)` crossing a BOUNDARY,
`1/(ℓ−w+1)` contained in a REGION — whose expectation cancels the opportunity **on the deposit's own
support** (§2, corrected 2026-08-20):

    E[ Σ 1/(w−1) ]      =  rho · P(w ≥ 2)   =  rho   at a BOUNDARY (frag_min ≥ 2, any pmf, any composition)
    E[ Σ 1/(ℓ−w+1) ]    =  rho · P(w ≤ ℓ)            in a REGION of length ℓ — NOT rho

(`tests/native/_accumulator_reference.py`; the region form is why `1/L` was wrong there and `1/(L−1)`
right at a boundary). ⛔⛔ **So the BOUNDARY form is a model-free TOTAL and the REGION form is only a
model-free density SHAPE**: `P(w ≤ ℓ)` is a pmf functional that differs per component — 11.6× at a 98 bp
exon on the ladder's own lengths, exactly zero below `ℓ < frag_min` — which is the same gDNA-vs-RNA
circularity this section removes from the divisor, surviving in the support
(`TRAPS: a-cancellation-is-conditional-on-its-support`). ⭐ A FACE's total adds the sj flux whose bodies
lie on that side, in the same units, from the sj bank's own reciprocal-opportunity column — without it the
comparison is between an object that can hold mature RNA and one that cannot
(`TRAPS: a-face-total-is-not-a-total-without-its-flux`).

⚠ **What is still model-dependent, stated so it is not assumed away**: the per-component divisors `E_g`
and `E_r` remain fragment-length models, so any statement about a COMPONENT's density still uses one. What
§3.5g buys is that the TOTAL at a **BOUNDARY** — and any enrichment ratio built between boundaries — does
not. ⛔ At a REGION the shipped bank's `P(w ≤ ℓ)` factor means every REGION↔BOUNDARY ratio still carries a
pmf functional with no enrichment in it (the chain alternates REGION/BOUNDARY, so this is
every hop).

**3.5h ⭐⭐⭐ THE PREMISE VARIANCE — why an imputation must cost something on every hop.**
(The ruling is `DESIGN.md` §0c.0c.)

Every variance that scales with counts vanishes between two deeply-counted slots, so a message layer built
only from counting terms transports an imputation for FREE and delivers it at full strength beside a
measurement. The premise of a hop — *"my neighbour's values apply here"* — is not a counting statement and
does not shrink with depth. It is estimable from the chain itself by METHOD OF MOMENTS, which is what keeps
it constant-free:

    Var_obs( log r )  =  premise  +  E[ v_r ]        ⇒     premise  =  max( 0,  Var(log r) − mean(v_r) )

⛔ Floored at 0 rather than at something small: a substrate whose ratios vary no more than Poisson predicts
has exhibited no heterogeneity, and a fit may not manufacture doubt the data do not show. Charged once per
hop, under either strategy, it makes a deep imputation arrive weaker than a shallow one and a measurement
outweigh both. ⚠ ONE library-level scalar is the first form; a per-hop-type fit is the obvious refinement
and needs a substrate with enough hops per type to estimate on.

**3.6 ⭐⭐⭐ THE TWO FACES OF AN `intron|exon` BOUNDARY — component-set matching, derived and verified.**
(Owner, 2026-08-04. Verified against oracle truth on `toy_harness --spec spliced_exons`.)

At **one boundary** the accumulator stores three populations, and their COMPONENT SETS DIFFER — which is the
whole content of this section:

| bank | what it counted | components |
|---|---|---|
| `unspliced_count` = `U` | crossed the boundary **contiguously**, spliced nowhere | gDNA + **nascent** |
| `junction_count` = `J` | never crossed it — it **JUMPED** from here | **mature**, certified |
| `spliced_count` = `S` | crossed contiguously, spliced *elsewhere* | mature, certified |

⛔ Mature RNA cannot cross an exon↔intron boundary contiguously (TRAPS: mature-rna-never-crosses-a-boundary), so it is **absent from `U`**
and present in `J`. Three densities, one of which needs no deconvolution at all:

    rho_g   = C_g / E_g       unknown split of U          U = C_g + C_nas
    rho_nas = C_nas / E_r     unknown split of U
    rho_mat = J / E_J         ⭐ MEASURED — certified RNA, nothing to deconvolve

Now match component sets against each flank:

    T(INTRON) = {gDNA, nascent}            = T(BOUNDARY, U only)        ⇒ shares transferable
    T(EXON)   = {gDNA, nascent, mature}    = T(BOUNDARY, U + J)         ⇒ shares transferable

so, with **the same `rho_g` in both numerators** and only the denominator changing:

    (I)  INTRON face:  phi_g(BOUNDARY) = rho_g / (rho_g + rho_nas)             ==  phi_g(INTRON)
    (II) EXON   face:  phi_g(BOUNDARY) = rho_g / (rho_g + rho_nas + rho_mat)   ==  phi_g(EXON)

⭐⭐ **One BOUNDARY, one gDNA density, TWO composition statements — differing only in whether the sj
term is in the total.** This is what makes the BOUNDARY a *solvable* object rather than a relay: (I) plus the
boundary's own mass identity `U = rho_g·E_g + rho_nas·E_r` is two equations in two unknowns, and the intron's
composition is prior-free (the intron factory). Then (II) delivers the exon's composition with **every
term measured** — no imputation, no prior, no strand, no length model.

**Verified against oracle truth**, `g50 ss0.50 capture_off`, residual in `phi_g`:

| arm | INTRON face | EXON face |
|---|---|---|
| m=100, `nrna_none` | +0.000000 | −0.0011 / −0.0031 |
| m=100, **nrna = 60** (the non-trivial control) | +0.0266 / −0.0198 | +0.0039 / −0.0012 |

⚠ The `nrna_none` arm tests (I) only trivially (both sides are 1.000); the nascent arm is the real one, and
there the boundary's nascent density reads 0.572 / 0.631 against the intron's 0.578 — the two boundaries straddling
it on 15 and 9 gDNA counts against the intron's 349, i.e. Poisson, not bias.

⛔⛔ **THE ESTIMATOR IS NOT THE IDENTITY, AND THIS IS WHERE IT FAILS.** The sj sees only
`E_J/(E_J + Σ E_r)` ≈ **11 %** of a two-exon transcript's mature fragments, so at low RNA `rho_mat` is a
handful of counts: at m=1 with `J = 3`, face (II) returns 0.715 against a truth of 0.458. And a BOUNDARY's
gDNA count is 8–36 where the intron's is ~349. So:

* ⭐ **face (I)'s job is to transport a WELL-COUNTED gDNA level from the intron to the boundary**, not to take
  the boundary's 8 counts as authoritative — measured at m=1, referencing the intron gives the exon
  `f_g = 0.499` and referencing the boundary gives 0.350, against a truth of 0.458 and a shipped answer of
  0.625;
* ⭐ **face (II) closes the composition two ways** — via `rho_mat` (tight at high RNA, over-stated by
  §3.6b, useless at low) or via the exon's own mass identity closed with `rho_g` (strong at low RNA, a
  small difference of large numbers at high). They are strong in opposite regimes, so they FUSE by
  inverse variance rather than by precedence.

⛔⛔ **AND THE WELL-COUNTED SIDE IS NOT A FIXED SIDE.** The 349-vs-8 ranking above is capture-OFF. Under
capture the intron REGION holds **1** count — median 0, max 3, unchanged from 120 kb to 1.08 Mb, because its
bp is fixed and its density is off-probe — while the `intron|exon` BOUNDARY holds **20–40** at the exon's own
capture stratum. So the direction of transport inverts with capture, and only the SHARE survives that:
transferring the intron's DENSITY instead is measured at **+0.207** mass-weighted `|Δf_g|` on capture-ON ×
unstranded. TRAPS: capture-inverts-the-counted-side.

⛔⛔ **WHAT THE CEILING SAID, AND IT IS THE REASON NOTHING WAS BUILT** (2026-08-04, `toy_ceiling.py` then
`ladder_arm_ab.py`). Handing both `intron|exon` BOUNDARIES the ORACLE truth and **re-solving the whole chain**
is worth **−0.000** on capture-OFF × unstranded and −0.033 on capture-ON × unstranded; the achievable
form (the intron's share) captures 82 % of that. Carried to the 36-condition ladder it is **negative**:
solvable mwae 0.0413 → 0.0426, confidently-wrong 20,173 → 22,336. ⭐ The identity above is not what
failed — improving that BOUNDARY is simply not worth anything where the tool is wrong. §9a is the general form of that reading.

**3.6c ⭐⭐⭐ THE SPLICE-FLUX REFRAME — AN BOUNDARY HAS TWO TOTALS, ONE PER FLANK.** (Owner's framing
2026-08-05; derived and gated the same day. `test_splice_flux_reframe`, `region_total_density`.)

§3.6 gives the two faces of an `intron|exon` BOUNDARY as a property of the OBJECT. Made per-STEP it becomes a
statement about which of the two flanks a hop is talking to, and then it applies at every BOUNDARY with
sj flux rather than only where the coarse types happen to read `intron|exon`.

The reframe is a composition imputation, so numerator and denominator must be totals over the **same
component set** — the intersection of what the two slots can carry. At a BOUNDARY the accumulator holds two
disjoint populations: what crosses **contiguously** (`U`, a gDNA/RNA mixture) and the **sj flux**
(`J`, certified mature). A molecule counted in `J` spliced *at this position*, so its body lies in the
exon on exactly **one** side of it and it never enters the other flank:

    at a sj's genomic-LOW  end   its exon is on the LOW  side
    at a sj's genomic-HIGH end   its exon is on the HIGH side

⭐⭐ Hence one total per flank, and the split is **per sj**, summed over the sj attached that
way:

    rho_lo  =  rho_U  +  Σ_{j : low end here}   J_j / E_J,j          — used against the LOW  neighbour
    rho_hi  =  rho_U  +  Σ_{j : high end here}  J_j / E_J,j          — used against the HIGH neighbour

and at a REGION both sums are empty — a REGION stores only CONTAINED fragments and a contained fragment used
no sj — so `rho_lo = rho_hi = rho_U` there and every junction-free chain is unchanged.

⭐⭐⭐ **THE PAIRING RULE, AND DIRECTION DOES NOT ENTER IT.** A hop joins adjacent slots `(k, k+1)`.
Whichever is the source, the pair is the same pair, so

    r  =  rho_lo[k+1] / rho_hi[k]

always. Travelling low→high (mature departing — a **splice-out**, `DESIGN.md` §0) versus high→low (mature
arriving — a **splice-in**) changes only which is numerator; it never changes which total each slot
presents. ⛔ **But it is NOT one array per direction**: within ONE forward pass a BOUNDARY at a sj's low
end is the DESTINATION of the hop from its low flank (flux INCLUDED) and the SOURCE of the next hop into
its high flank (EXCLUDED). Two arrays indexed by ROLE is what expresses that; one per pass is not.

⛔⛔ **WHAT THE PREDECESSOR DID AND WHAT IT COST.** One junction-inclusive total per slot, on every hop, in
both directions and in both twins. Measured against origin-split truth on `g50 ss0.50 nrna_none
capture_off`, the intron-facing side of the two `intron|exon` BOUNDARIES is inflated by exactly `J/E_J`:

| step | shipped `r` | flank-pair `r` | TRUTH ratio | shipped/true | pair/true |
|---|---|---|---|---|---|
| intron → BOUNDARY @9,000 | 1.4595 | **1.1789** | 1.1421 | ⛔ 1.28× | ✅ 1.03× |
| intron → BOUNDARY @2,000 | 1.0061 | **0.7255** | 0.7029 | ⛔ 1.43× | ✅ 1.03× |

⚠ It inflates the two hops of a two-hop pair in opposite directions and therefore **cancels in a
compounded ratio** (2.159 used vs 2.153 true), which is why no endpoint, conservation or aggregate check
saw it — TRAPS: recompute-from-the-oracle. ⛔ And the opposite error is not available either: using `rho_U` on both sides was
measured WORSE, because at a BOUNDARY→EXON step the exon genuinely contains the spliced population.

⛔⛔ **WRITE THE PREDICATE IN GENOMIC TERMS — §3.5b's ruling, and here is where it bites.**
`splice_graph`'s `FLAG_DONOR_s` marks the genomic-LOW end of an `s`-strand intron on **both** strands
(`don_bit ← intron_start`, inside a loop over strand that changes nothing else), so on `−` it sits at the
transcript's biological **acceptor**. The names are a misnomer on `−` and the data is uniform, which is
exactly the arrangement in which a genomically-phrased predicate is safe and a biologically-phrased one
flips sign silently. The derivation above never asks the biological question, the fields are named
`_lo`/`_hi`, and keying the split on the sj's strand instead fires **exactly one** of twelve gates.

⭐ **And a rule keyed on the coarse region type provably cannot do this.** On `splice_both_strands` three
regions are simultaneously intron and exon and `coarse_type_array` reports **every** gene-body region as
`exon`, so `intron|exon` never appears. Nor is "does the neighbour admit mature RNA on strand `s`?"
enough: at BOUNDARY @3,000 the flanking region carries mature⁺ (from a single-exon transcript spanning the
intron) but not the mature⁺ population *that sj's* flux belongs to.

⛔⛔ **BUILT AND MEASURED, AND IT DOES NOT PAY ON ITS OWN (2026-08-05).** The derivation above is not what
failed — it reproduces the truth ratio to 3 % where the shipped total is 28–43 % off. What the panel says is
that the **BOUNDARIES improve and the REGIONs get worse**: on a 6-condition shard spanning all four strata, the
boundary axis moves mwae 0.16046 → 0.16000 with confidently-wrong −7.7 % and the shipped solve better on 6 of 6,
while the region axis moves 0.12232 → **0.12297** with confidently-wrong **+36.9 %** and Σ|err| +20,884
fragments. `solv%` is byte-identical in both arms, so this is not a moved denominator.
⭐ **The reason is TRAPS: a-cancelling-defect-pair and it was predictable from TRAPS: recompute-from-the-oracle**: an evidence-free exon is fed through
`intron → BOUNDARY → exon`, the two hops' errors cancel, and the second hop carries a *different* defect (§3.5's
composition ratio applied to a level). Correcting one of a cancelling pair is worse than correcting neither.
⛔ **So this is priced jointly with that defect or not at all** — `EQUATIONS.md` §3.6c.

**3.6b ⭐⭐ THE SJ AND CONTAINED FRAMES ARE A LEVER ON THE FRAGMENT-LENGTH MEAN, AT 0.62 %/bp.**

⚠ **This is the length model as a DIVISOR — live, shipping, and NOT the deferred composition channel of
§3d–§3e.** A fitted-mean error here moves an *opportunity*, and it does so whether or not anything ever
reads length as evidence about composition.

`E_J` and the exon's `E_r` are built from one pmf and are exactly consistent — measured 202.8 and 797.2 on
a 1,000 bp exon, summing to 1,000.0. But they differentiate with **opposite sign**:

    E_J   = E[w] − 1          →   dE_J/dE[w]  = +1
    E_r   = e − E[w] + 1      →   dE_r/dE[w]  = −1
    ⇒  d log(rho_mat / rho_R) / dE[w]  =  1/E_J + 1/E_r  =  0.0062 per bp

so a length model that is wrong by `Δ` reports the sj estimator as `1 + 0.0062·Δ` times the exon's
own. Measured: the fitted RNA pmf's mean is **203.80 against a true 212.20 on all 36 ladder conditions**
(−8.4 bp off capture, −3.5 bp under it), and −8.4 bp predicts **1.1125** against a measured 1.103 / 1.111.
✅ **AND THERE IS NO SECOND TERM.** A "finite-transcript placement" contribution of 1.024 was recorded
here and is **WITHDRAWN**: it assumed the simulator draws a length and then places it uniformly among the
`L − w + 1` legal starts. It does not — `wgs_engine._post_capture_length_allocation` reweights the length
MARGINAL by that same opportunity (`f_post(w) = f_pre(w)·total_eff(w)/Z`), so the realised crossing count
at length `w` is `f_pre(w)·(L−w+1)·(w−1)/(L−w+1) = f_pre(w)·(w−1)` and the factor cancels identically.
Computed exactly from the pmf, `k = 1.000000` at every transcript length. ⛔ The whole gap is the
length-model mismatch, with no geometric component. TRAPS: two-divisors-opposite-sign.

✅ **This WITHDRAWS the sign correction this section used to carry.** `calibration/messages/variance.py`'s
`splice_in_premise_logvar` (P1d) asserts
`rho_R(exon) ≥ rho_nas(B) + rho_mat(B)`, a LOWER bound, and that is **right**: the measured ratio of
**1.103** (no nascent) and **1.049** (with) is `1 + (1−s)(k−1)` with `k` the frame gap above and `s` the
nascent share of the exon's RNA — the nascent arm is measured in the exon's own frame and dilutes it,
which is precisely what identifies the inflation as `rho_mat`'s rather than the bound's. ⚠ What does
survive is that `splice_in_premise_logvar` is fitted on fluxes carrying the same inflation, so it inherits it.

---


### 3.7 ⭐⭐⭐ THE FRAME-FREE COMPOSITION COORDINATE — and what it does and does not buy

    ⭐ Promoted verbatim from a working doc (`variance_model_notes.md` §4–§5) when it was
    deleted (2026-08-07). ⛔ **The DERIVATION below is sound and is worth keeping; the REBUILD that was
    based on it measured +85–103 % on the panel and is deleted.** `DESIGN.md` §6.1 records why those two
    facts are compatible: the rebuild's regression decomposed into a correct derivation, a UNIT ERROR and a
    STRUCTURAL disconnection, and only the first of the three is written here.
    ⚠ So read this as *"what a composition coordinate would have to look like"*, not as a design that is
    waiting to be re-attempted.


⛔ **The first draft of this document claimed "composition transfer is the identity on λ". THAT IS FALSE**,
and the error matters because it is the difference between a belief-free exact operation and a wrong one.

**What actually transfers is the DENSITY, not the share.** "gDNA is uniform" and "the same transcripts run
through both slots" are statements about *densities*. Shares are derived, and they move whenever the two
slots' opportunity ratio differs — which between a REGION and a BOUNDARY it always does. Verified:

| `E_g,E_r` at source → destination | `φ_g` | `λ` | `η` |
|---|---|---|---|
| (100, 200) → (250, 100) | 0.3333 → 0.7143 | −0.6931 → +0.9163 | **0 → 0** |
| (3000, 1500) → (254, 254) | 0.6667 → 0.5000 | +0.6931 → +0.0000 | **0 → 0** |
| (50, 400) → (900, 300) | 0.1111 → 0.7500 | −2.0794 → +1.0986 | **0 → 0** |

Since `λ = log(φ_g/φ_R) = log(ρ_g/ρ_R) + log(E_g/E_r)`, define

$$\boxed{\;\eta \;\equiv\; \lambda - \log\frac{E_g}{E_r} \;=\; \log\frac{\rho_g}{\rho_R}\;}$$

**`η` is the frame-free composition coordinate, and it transfers as the IDENTITY across any boundary with
`T_s = T_d`.** The tilt `θ` — built from `(ρ_p − ρ_n)/ρ_R`, in which the common `E_r` cancels — likewise
transfers as the identity.

**Four consequences, and together they are the whole simplification:**

1. **The conversion is a known geometric constant.** `λ_d = λ_s + [log(E_g/E_r)_d − log(E_g/E_r)_s]`. No
   belief, no count, no unknown — read off the index.
2. **The mass identity holds by construction.** Rebuild at the destination as
   `ρ_c(d) = φ_c(d)·M_d/E_c(d)` with `φ` from `λ_d`; then `Σ_c ρ_c E_c = M_d Σφ_c = M_d`, exactly, for any
   `η` (verified to 1e-10). ⭐⭐ **So there is nothing for the mass pin to restore.**
3. **No ratio of totals is ever formed** — so no reframe `r`, no `framed` guard, no undefined ratio, no
   `rho > _EPS` guards.
3b. ⭐⭐ **AND THIS — NOT §2 — IS WHAT DISSOLVES TRAPS: a-ratio-cannot-carry-zero.** TRAPS: a-ratio-cannot-carry-zero says a MULTIPLICATIVE transport
   cannot carry zero, so "there is none here" is unrepresentable by construction. Under §4/§5 the
   conversion is an **additive constant in log space** and the shared density crosses as the **identity**:
   `rho_g = 0` needs no `rho > _EPS` guard, is never a denominator, and arrives as `f_g = 0` intact. ⛔ The
   first draft credited this to §2's Jeffreys location, which is the opposite of true — §2's positive floor
   is what *destroys* the zero, and it cost +96,299 % on the control to find out (§6).
4. ⭐ **It explains what the reframe was for.** `r = ρ_tot(d)/ρ_tot(s)` is the code compensating, with
   *beliefs*, for working in `λ` where it should have worked in `η`. The correct compensation is exact and
   belief-free.

⚠ **The premise `η` inherits, stated plainly:** `η` is invariant only where gDNA density and RNA density are
each unchanged — i.e. within one capture stratum, with the same transcripts. `T_s = T_d` is a proxy for the
second condition; the first is the same uniform-gDNA assumption the whole tool rests on. This is not a new
weakness, but §4 does not remove it.

⚠ The **spliced** bank is untouched: a spliced fragment is certified RNA, entering as a *measurement* of the
RNA side. The identity above is over the **unspliced** bank, which is what needs deconvolving.

---

### 3.8 THE MISMATCHED BOUNDARY — under the null there is NO UNKNOWN

`T_s ⊊ T_d`: the source (say an `intergenic|exon` boundary, `T_s = {g}`) knows only `ρ_g`.

⭐ **A density is already frame-free.** So under the null — gDNA uniform, no capture enrichment —
`ρ_g` crosses **unchanged**, and the destination's own observations convert it:

$$f_g(d) \;=\; \frac{\rho_g(s)\cdot E_g(d)}{M_d}, \qquad \text{RNA mass} = M_d - \rho_g(s) E_g(d)$$

**Fully determined. No interval, no free parameter, no scale to integrate out.** The newly-active
population is not "unconstrained by the message" — the null *implies* it as the residual of the
destination's own count. This is the owner's formulation and it is a strict improvement on the first draft's
interval framing.

⭐⭐ **And this is exactly what HEAD already computes** — `mo_g = log(c_g·E_g/M) = log f_g`.

#### 3.8.1 The general case — TSS/TES and strand changes, in one rule

Let `C = T_s \cap T_d` (shared) and `N = T_d \setminus T_s` (newly active). Under the null:

1. **Within `C`, the density RATIOS are frame-free and transfer as identities** (§4) — so the relative
   composition *among the shared populations* crosses exactly.
2. **Each shared density crosses unchanged**, so the destination-frame mass they account for is determined:
   $$A \;=\; \sum_{c \in C} \rho_c(s)\,E_c(d)$$
3. **`N` takes the residual** `M_d − A`, and how it splits depends only on `|N|`:

| `\|N\|` | example | the split |
|---|---|---|
| 0 | `T_s = T_d` | §4 — nothing new; the message is complete |
| 1 | a `+`-strand TSS: `{g} → {g, R⁺}`; or `{g,R⁺} → {g,R⁺,R⁻}` | **determined** — one unknown, one equation |
| 2 | `{g} → {g, R⁺, R⁻}` (an intergenic boundary into a both-strand exon) | the RNA **total** is determined; the **tilt** `θ` is not, and falls to ψ's reference |

⭐⭐⭐ **So the only thing a component mismatch ever leaves undetermined is the TILT, and only when both
strands are newly active.** The gDNA-vs-RNA split — the quantity the tool exists to estimate — is
determined in every case. That is a much smaller residual unknown than the first draft supposed, and it is
`θ`, a nuisance parameter, not `λ`.

⚠ `A > M_d` is possible from sampling noise (the shared populations "account for" more than the destination
observed). Then the null is inconsistent with the destination's own count and the honest reading is `f_N = 0`
with the excess absorbed by the shared claim's own uncertainty. ⛔ Do **not** introduce a shift here; §6.

⛔⛔ **So the open problem is not "what crosses". It is: HOW MUCH SHOULD THE DESTINATION DOUBT THE NULL?**
Capture makes `ρ_g(d) = γ·ρ_g(s)` with `γ ≥ 1` possible; nothing at pass-0 is evidence that it *happened*.
§6 shows that every rule that introduces doubt is refused by the zero control, and §8 says where the doubt
must come from instead.

---

---

## 3b. ⭐⭐⭐ THE CONSERVED MASS, AND WHY ONE SHARE FOR TWO COMPONENTS IS A BIAS

Derived 2026-08-08 when `boundary_unspliced_mass` landed. The accumulator deposits `+1` on every boundary a
fragment crosses — `max(K, 1)` of them — so a sum over objects is an object-INCIDENCE count while the EM
adds a FRAGMENT count. The mass bank closes that gap: it sums to ONE per fragment.

**The deposit.** A fragment of length `w` is region bound by the crossed boundaries into slices; a slice of length
`s` bounded by `n_cross ∈ {1,2}` boundaries deposits `s/(w·n_cross)` at each. Every slice of a
single-segment path has `n_cross ≥ 1`, so `Σ = Σ s/w = 1` exactly.

**The opportunity.** Put the boundary at 0 with flanking regions of length `a` (left) and `b` (right). A
crossing fragment with `u ∈ [1, w−1]` bases to the left contributes `g_a(u)/w`, where `g_a(u) = u` for
`u ≤ a` and `a/2` beyond — the second branch because the slice then has a boundary on both sides. Summing:

    Σ_{u=1}^{w−1} g_a(u) = (w−1)w/2      if w−1 ≤ a
                         = a·w/2          otherwise      → both equal  w·min(w−1,a)/2

so, adding the mirror term and dividing by `w`,

    A_mass(w; a, b) = [ min(w−1, a) + min(w−1, b) ] / 2

exact regardless of how many further boundaries the fragment crosses, and `E[mass] = ρ_c · E_{f_c}[A_mass]`.
⚠ **It is a CENSORED functional**, so it is sensitive to the pmf's *shape*, not only its mean.

⭐⭐ **THE POOLING THEOREM — the result that matters.** Define a component's share at a boundary as the mean
conserved mass one of its crossings carries, `share_c = E_c[A_mass] / E_c[w−1] ∈ (0,1]`. The accumulator
can only measure the MIXTURE's, `share_pooled = M/C = φ·share_g + (1−φ)·share_r`. Rescaling both
components by it gives

    a_g + a_r  =  M_g + M_r                    ← the TOTAL is conserved, exactly
    (â_g/â_r) / (a_g/a_r)  =  share_r / share_g   ← and is INDEPENDENT of the true mixing ratio

⛔⛔ **So the error is purely COMPOSITIONAL, and no conservation check can ever see it** — the locus total
is right to the fragment while the split is wrong (`TRAPS: conservation-misses-mis-attribution`). It is
identically zero iff `share_g = share_r`, which holds whenever `f_g = f_r` — which is why an
equal-length panel is *structurally* blind to it (`TRAPS: an-equal-length-panel-defeats-the-lift`).

⚠ **That blindness is a PRICE THE LADDER PAYS ON PURPOSE, not an argument for giving it a length gap**
(§3.1): a gap lets the EM assign on length alone and bypass the phase the ladder exists to exercise. The
gapped `flgap` pair is where this bias **was** measured; the ladder is where calibration is measured.
⛔ **Both `flgap` panels were DELETED on 2026-08-13, so nothing on disk can measure this bias today, and
on 2026-08-17 their configs and the two instruments that read them went too** (`prior_yardstick.py`,
`flgap_study_cache.py` — dark since the panels went, and neither could be repointed at the ladder because
the whole premise is a DRAINED oracle, which the ladder refuses). ⛔ **Re-measuring therefore means
DESIGNING a length-gap panel, not restoring one**, and the paragraph above is the argument it has to
answer first. ⚠ This is a statement about what is reproducible, not a licence to give the ladder a gap.

⚠⚠ **THE SECOND BOUNDARY IS EXACT PER BOUNDARY AND NOT AT A LOCUS, AND THE GAP IS LARGE ENOUGH TO MIS-SIZE A
CORRECTION** (measured 2026-08-08, `tests/calibration/test_prior_units.py`). `share_pooled` cancels from
the ratio only where *every* fragment of both components passes through it. Two things break that at the
scale a prior is assembled on, and they break it by very different amounts:

* **Contained mass never passes through the share at all.** A contained fragment deposits on exactly one
  region, so that term is already a fragment count and is carried through untouched — it *dilutes* the
  tilt, and because `share_pooled` itself moves with `φ`, the dilution is mixture-dependent. Swinging
  the mixture 900× at a fixed length gap moved the locus-level distortion **0.578 → 0.767**, where the
  gDNA component was 53 % contained.
* **Summing over boundaries with different flanks** re-introduces a weak dependence, because the cancellation
  is per boundary and the two components weight the boundaries differently: **0.5030 → 0.5073** over the same
  900× swing.

⭐ So `share_r/share_g` is the right *mechanism* and the wrong *magnitude* for any correction applied at
locus granularity. A repair belongs where the identity is exact — **per boundary, before the contained term
is added** — and a per-locus factor derived from this ratio would be wrong by the contained fraction.

⭐⭐ **MEASURED 2026-08-08, DRAINED, ON ALL FOUR flgap CONDITIONS: THE POOLED SHARE IS 82–99 % OF THE
ASSEMBLER'S WHOLE RESIDUAL, AND THERE IS NOTHING ELSE LEFT TO EXPLAIN.** With the truth masses in, the
assembler misses the EM's own candidate count by `rel` 8.8e-4 … 0.0111; give each component its own true
per-line share as well and that falls to **2.8e-5 … 2.0e-3** — 67 to 9,653 fragments out of 2.4–4.9 M.
⛔ **This paragraph is the ONLY home of that verdict.** It was measured by `prior_yardstick.py`, which was
deleted on 2026-08-17 along with the `flgap` panels its measurement required — so nothing on disk
reproduces these numbers, and re-deriving them means designing a drainable length-gap panel first.
⭐ Also measured there, and recorded where it belongs: the assembler's composition claim is exact against
the UNSPLICED pool (`DESIGN.md`, `Δphi ≤ 5e-4`); scoring it against ALL RNA units instead reads a phantom
`+0.07…+0.10` tilt (`TRAPS: score-the-consumers-own-count`); and the prior's STRENGTH is one
pseudo-fragment per real unspliced fragment by construction, with no knob (`DESIGN.md`).

⛔ **The predecessor of this paragraph said "25–28 % … the other ~72 % survives perfect per-component
shares", and the 72 % did not exist.** It was scored against a per-locus count of the fragments whose
FIRST BASE lands in the locus, and the EM counts every fragment that is a *candidate*
(`TRAPS: score-the-consumers-own-count`). ⚠ Read the ratios rather than the absolutes: `flgap_long`
carries a real spliced leak in its drained gDNA partition (8,641 records against `flgap_short`'s 1), so
its 0.0020 is an upper bound and `flgap_short`'s 4.3e-5 is the admissible figure.

⚠ **And `S` is a diagnostic ceiling, not a design.** It consumes the true per-component share, which
production cannot have. What the measurement establishes is that the assembler's *arithmetic* is right
and that a per-line share correction is the whole of the remaining work — not that the work is done.

## 3c. ⭐⭐ THE MODEL-FREE LOCAL MEAN LENGTH — and its one precondition

At a contiguous boundary the opportunity is `w−1` and the deposit `1/(w−1)`, so they cancel:

    E[count] = ρ·E[w−1],   E[inv_length_sum] = ρ    ⇒    count / inv_length_sum + 1 = E[w]

⭐ **No pmf, no model, per boundary.** Given the two global component means it inverts directly to a local
gDNA fraction: `φ = (E[w] − μ_r) / (μ_g − μ_r)`. ✅ Verified numerically (400 k fragments, μ 326/208,
φ = 0.40): reads 256.6 against a true mixture mean of 255.1, giving `φ̂ = 0.413`.

⛔⛔ **THE PRECONDITION, AND IT IS EASY TO VIOLATE: the identity holds only under NATURAL PLACEMENT.**
A longer fragment is more likely to cross a given boundary, and the `1/(w−1)` deposit is exactly what
cancels that length bias. Evaluated on a length-SELECTED subpopulation the cancellation fails and the
ratio becomes a harmonic-type mean instead: forcing every fragment to cross a boundary in a check of this
identity returned **185.9 — below both component means** — which is impossible for a mixture mean and is
the signature of the violation.

⭐ **And note what the denominator is**: `μ_g − μ_r`, the same 2×2 determinant §3 relies on. The
fragment-length gap that biases every divisor is also the *only* θ-independent composition evidence an
AMBIG slot can get. ⛔ **So the tool may not be made gap-robust by shrinking the estimated gap** — that
destroys the identifying quantity. Robustness must come from per-component divisors.

⛔⛔ **AND THE EXIT THAT LAST PARAGRAPH IMPLIES — hand an AMBIG slot its composition from a local mean
length — IS DEFERRED POST-0.8.0 AND MUST NOT BE PROPOSED FOR THE RELEASE** (owner, 2026-08-14). It is the
same channel §3d–§3e derive, and it has been measured: A/B'd on the flgap pair on 2026-08-10 it is better
on 7 of 8 conditions and *appears* to resolve the blindness (0.0324 → 0.5222 against a truth of 0.507) —
and on the `g00` ZERO-gDNA control it reports **54–57 % gDNA in a library containing none**, because it
returns a near-CONSTANT regardless of the truth (~0.5 on that arm, 0.59–0.72 on the gap sweep of §3d–§3e)
while the flgap panel is all `g50` (`TRAPS: a-single-level-panel-cannot-see-a-constant`).
⭐ **The identity above is untouched by that and is still live** — it is a statement about `E[w]`, which
is a measurement, not about `f_g`, which is the inversion that was refused.

## 3d–3e. ⛔⛔⛔ THE LENGTH CHANNEL — DEFERRED POST-0.8.0, AND THE DERIVATIONS ARE KEPT

⛔⛔ **OWNER RULING, 2026-08-14. Fragment length as a CALIBRATION composition channel is OUT OF SCOPE for
0.8.0, and 0.8.0 ships without it.** Do not propose it, do not list it as a candidate, do not rank it
against anything, and do not open it as "the next thing to try" — it becomes admissible again only after
the release, and only against the objections recorded below.

⭐ **What is deferred, precisely, so the ruling is not over-read:**

| | |
|---|---|
| ⛔ **DEFERRED** | `length_likelihood` as CALIBRATION's composition channel — a length row over the `λ` grid, the precision derived in §3d, the moment repair in §3e, and any per-slot `f_g` inverted from a local mean length (§3c) |
| ⭐ LIVE | the length **model** as an OPPORTUNITY and a DIVISOR — `E_g`, `E_r`, `E_J` (§1), the four gDNA pools and the sj pool (§4), the frame gap §3.6b prices, and the opportunity-tilted moments in `calibration/effective_length.py`, where they were always a geometry rather than a composition claim |
| ⭐ LIVE | `length_likelihood` in `src/rigel/second_pass.py` — a **DIFFERENT thing** with the same word in its name: the per-fragment assignment term of §10. The deferral does not touch it |

⚠ **This is a SCOPE ruling about what may be proposed, not a code removal waiting to happen.** The module
is already out of the tree: `calibration/length_likelihood.py`, its config flag, its wiring and the seven
instruments that priced it went in `b7ed7a0b` (2026-08-10), byte-identical on every production path
because the flag defaulted False.

⭐⭐ **THE TWO SECTIONS BELOW ARE KEPT DELIBERATELY, AND DELETING THEM AGAIN IS THE WRONG MOVE.** They are
a derivation and its measurements; a derivation is a record, and this one will be wanted after 0.8.0, when
the objections should be the starting point rather than a surprise. ⚠ They were written on 2026-08-10
while the module still existed, so read every present tense in them as of that date, and every script they
name as recoverable from `f470a570` rather than runnable now. Module names are the current ones
(`region_init`, not the `node_init` of the day).

⛔⛔ **AND READ THE VERDICT BEFORE THE DERIVATIONS, because it is not a pricing question and §3d does not
answer it.** The channel's answer was **not a function of the pmf gap at all**: at a gap of ~1e-9 bp it
reported 0.72 / 0.59 / 0.72 on libraries whose truth was 0.00 / 0.00 / 0.57, and CLOSING the gap made
every one of them worse. The mechanism is `TRAPS: a-linear-likelihood-emits-a-sign` — a Gaussian
log-likelihood is asymptotically LINEAR in the composition, so its argmax is a sign saturated at a grid
endpoint — with `TRAPS: amplitude-fades-influence-does-not` for why repairing the precision (§3d) could
never have been enough on its own. ⭐ So a post-0.8.0 attempt starts at the MODE, not at the precision and
not at the moments.

### 3d. ⭐⭐⭐ THE LENGTH CHANNEL'S OWN PRECISION — smooth shrinkage, not a gate

⚠ **DERIVED AND MEASURED, NOT LANDED** (2026-08-10). Recorded here so it is not re-derived; the reason it
did not land is at the end.

**The problem.** `length_likelihood` returns a log-likelihood row over the `λ` grid, and `region_init`
registers its curvature as source 4 of a slot's composition precision. That curvature is measured
NUMERICALLY — normalise `exp(row)` over the grid, take the variance, invert. ⛔ A flat row normalises to
the UNIFORM distribution on the grid, so the precision has a FLOOR that is not information:

    tau_flat  =  1 / Var(uniform over the lambda grid)  =  0.029016

Measured at `g00` (zero gDNA, the two fitted pmfs **1.2 bp** apart): the channel reported bias **+0.66**
and median `|Δ| = 1.0000` on 5.1 M fragments with `med tau = 0.029` — the grid's own width sold as
evidence, on 43–99 % of slots. ⛔ And it lands where it does most harm: `tau_len` is deliberately
UNGATED, because the channel exists to speak on AMBIG slots where every other channel is silent, so
there a grid-width claim at an arbitrary grid EDGE is the ONLY evidence and it wins.

⛔ **A wider gate is not the repair** (owner, 2026-08-09: this codebase prefers a channel that fades to
one that switches, and three thresholds in this family have already been refused). The channel's own
Fisher information fades on its own.

**The derivation.** The likelihood is a bivariate Gaussian on `(Σu, Σw)` whose MEAN moves with `pi`:

    mean = N·mu(pi),  mu(pi) = pi·mu_g + (1−pi)·mu_r      cov = N·V(pi)

For a Gaussian whose mean depends on a parameter the information from the mean is
`(dm/dpi)' Sigma^-1 (dm/dpi)`, and with `dm/dpi = N·Delta` where `Delta = mu_g − mu_r`:

    I_pi  =  (N·Delta)' (N·V)^-1 (N·Delta)  =  N · Delta' V^-1 Delta
    I_lam =  I_pi · [pi(1−pi)]²                              (push forward, pi = sigma(lam))

⭐ `Delta' V^-1 Delta` is the **squared Mahalanobis distance between the two components' length
signatures** — how distinguishable gDNA and RNA are by length here, in the metric the noise itself sets.

⭐⭐ **Why it is the right shape.** `Delta → 0 ⇒ I → 0` QUADRATICALLY, so the channel fades smoothly as
the two pmfs converge and `Delta = 0` exactly is merely the limiting case — the exact-inequality guard
stops being load-bearing and becomes numerical hygiene. It is linear in `N`. It never touches the grid, so
it has neither floor nor ceiling. And it introduces no constant.

⛔ **It deliberately omits `½ tr[(V^-1 dV/dpi)²]`**, the heteroscedastic term from the covariance also
moving with `pi` — the same `−½ log det` that displaces the peak by **0.32 at N = 1**. That term is
`O(1)` in `N`, so it does NOT vanish at low depth where it is least trustworthy. Omitting it is the
conservative reading.

**Measured** (`length_channel_census.py` table ④, deleted with the channel), REGION slots at `N >= 50`:

| | tabulated `tau` | analytic `I_lam` |
|---|---|---|
| `g00`, pmf gap **1.2 bp** | 0.0292 (43–99 % of slots pinned at the floor) | **0.0021** |
| `flgap_long/OFF`, gap **110 bp** | 0.2032 | **10.68** |

⭐ `tau` separates "no gap" from "real gap" by **7×**; the analytic form by **5,000×**, growing smoothly
with `N` (0.0011 → 10.68).

⛔⛔ **WHY IT DID NOT LAND.** `tau` is bounded at BOTH ends — a floor at `1/Var(grid)` and a ceiling at
the grid spacing — so `I/tau` runs **0.04 → 52** across the depth range. The analytic form does not only
shrink uninformative slots, it **amplifies confident ones 40–50×**, faithfully amplifying whatever the
moments underneath say including their error (§3e). ⭐ It affects ONLY the precision, never the mode
(`tau_lam` feeds `region_init.own_composition_logvar`; `fg_loc` comes from ψ) — ⛔ which is exactly why it
was never the repair: the verdict in the banner above is a MODE failure, and a precision that is right
about a mode that is wrong is worse, not better.

### 3e. ⚠ THE LENGTH CHANNEL'S MOMENTS MUST DESCRIBE THE POPULATION IN ITS BANK

The channel models the length distribution of fragments landing at an object as `g_c(w) ∝ f_c(w)·A(w)`
— the library pmf times the opportunity. ⭐ The opportunity is exact for a UNIFORMLY placed component and
that is verified rather than assumed: gDNA cannot splice, and off capture its predicted moments match the
realised bank to **1.000** (REGIONs) and **0.996** (BOUNDARIES), at every region-length stratum.

⛔ **The RNA arm is off by 2.3 % (REGIONs) and 7.5 % (BOUNDARIES), and it is NOT the pmf.** `pi(w)`'s
de-tilt of the spliced pool is essentially exact — spliced-pool ÷ `pi(w)` reads **211.77** against a true
library **212.20** with per-bin ratios 1.005 … 0.979 — and feeding the TRUE pmf in place of the fitted one
leaves the residual almost unchanged (193.49 → 194.00; 246.57 → 247.19). What is missing is a
per-population SELECTION term: at a BOUNDARY, `A(w) = w − 1` counts every crossing while the bank
`boundary_unspliced` holds only the unspliced ones (one that spliced elsewhere is in `boundary_spliced`).
At a REGION nothing is missing — containment excludes splicing structurally.

⛔⛔ **AND THE OBVIOUS REPAIR IS NOT ADMISSIBLE.** Bounding the crossing opportunity by the containing
exonic block closes 62–67 % of the BOUNDARY deficit with no free parameter, but it assumes every RNA
molecule at a BOUNDARY is confined to an exonic block. RNA that has not spliced at the intron it meets
runs straight through it and violates that; and the block's extent is per-TRANSCRIPT, which is what the
EM solves, so the constraint is circular (owner, 2026-08-10). ⭐ `A(w) = w − 1` is the CONSERVATIVE
bound — exactly right for gDNA always, exactly right for unspliced RNA always, and too permissive only
for RNA that spliced elsewhere — so it errs in one direction for one sub-population. ⛔ The ladder this was measured on held `nrna_none` on every condition — the substrate most favourable to
the tighter bound; a repair validated there is not validated. ⚠ The rebuilt ladder (2026-08-19) carries
nascent on every row, so the substrate now exercises the permissive direction at STRESS level
(`DESIGN.md` §0b, the NASCENT SCOPE RULING).
⚠ `length_channel_census.py` regenerated every number here and `short_exon_fl_probe.py` priced the
tighter bound; both were deleted in `b7ed7a0b` and are recoverable from `f470a570`.

## 4. Opportunity corrections for length pools

A pool is a length-dependent **selection**, so its raw histogram is not the library's length
distribution. Let `A(w)` be the pool's opportunity and `T(w)` the total opportunity for the same
population.

**4.1 The divisor is a probability.**

    pi(w) = A(w) / T(w)          fitted(w) = count(w) / pi(w)

⛔ Never `count(w)/A(w)` — see TRAPS: divide-by-a-probability.

**4.2 The sj pool** (`calibration/sj_opportunity.py`). For a transcript with exon lengths
`e_1..e_K` and total `L = Σ e_i`, the starts at which a length-`w` window crosses **at least one**
sj:

    A_j(w)  =  (L − w + 1)₊  −  Σ_i (e_i − w + 1)₊

⭐ Derived via the **complement** — a window crosses no sj iff it lies wholly inside one exon, and
the exons are disjoint — so there is no inclusion-exclusion. The library quantities are abundance-weighted
sums, `T(w) = Σ_t θ_t (L_t − w + 1)₊` and `A(w) = Σ_t θ_t A_j(w,t)`.
⚠ `θ` is a **molar** abundance (copies), not an observed fragment count — `A_j` already counts start
positions, so a count applies the length weighting twice. ⭐ Production uses a **uniform** `θ` over the
non-synthetic transcripts: the ratio cancels most of the dependence, so uniform lands within 0.2 pp of the
simulator's own molar abundances.

**4.3 The four gDNA pools** (`calibration/gdna_opportunity.py`). Two contained (§1.3), two
crossing-exactly-one-line (§1.4), with `T(w) = Σ_refs (L_ref − w + 1)₊`. The combination:

    f(w)  ∝  [ Σ_p count_p(w) ] · T(w) / [ Σ_p A_p(w) ]

⭐ **This is the opportunity-weighted average of the four de-tilted pools**, `Σ_p A_p f_p / Σ_p A_p`, and
under Poisson counts `Var(count_p) ∝ A_p`, so those weights are exactly **inverse-variance**. There is no
tunable weight. ⛔ It is *not* the same as pooling the four histograms and applying one divisor
(TRAPS: opposite-tilts-must-not-pool).

**4.4 What an opportunity correction cannot do.** Both forms above assume the population is placed
**uniformly** over its template. Under hybrid capture gDNA is not: placement is proportional to a capture
landscape the tool cannot see. That residual is a *placement* model, not a better divisor.

---

## 5. Strand

**5.1 The likelihood — the only intrinsic gDNA/RNA signal.**

    p    =  ½·f_g + κ·(1−f_g)
    var  =  N·p(1−p)  +  (N·f_g)²·¼·od_g  +  (N(1−f_g))²·κ(1−κ)·od_r
    loglik = −½·(sense − N·p)²/var − ½·log(var)

**5.2 What it can and cannot say.** With RNA tilt `d = f₊ − f₋`, `p = ½ + (κ−½)·d` — **the gDNA fraction
cancels identically** (verified to 5.6e-17). Strand measures the tilt; it reaches gDNA only through the
triangle bound `f_g ≤ 1 − |d|`.

    I(f_g)  =  N·(2κ−1)² · [f_g(1−f_g)]² / (4·p(1−p))

**exactly zero at κ = ½** for any count and any overdispersion, and saturating in `N` at
`(½−κ)²/(p(1−p)·od)` — capped by dispersion, not by depth.

**5.2b ⭐⭐ WHAT THE CODE ADDITIONALLY DOES — the DERIVED noise-floor DEADBAND, which §5.2 above does not
describe.** `κ` is **fitted**, so `(2κ−1)²` is a squared point estimate and is strictly positive on a
genuinely unstranded library. `region_init.strand_evidence` therefore does not evaluate §5.2's `I(f_g)`
verbatim — it evaluates

    N_eff  =  N / (1 + (N−1)·od_r)                          the OVERDISPERSED effective count
    σ²_d   =  ¼·(1/N_rna + od_r)  +  ¼·(1/N_gdna + od_g)     the noise floor on κ̂
    disc   =  4·max(0, (κ−½)² − σ²_d)                        ⭐ REPLACES  (2κ−1)² = 4(κ−½)²
    I_strand  =  N_eff · disc · [f_g(1−f_g)]² / (4·p(1−p))

⚠ `p` here is §5.1's `p = ½·f_g + κ·(1−f_g)` — written in the code as `κ + f_g·(½−κ)`, the same
expression — and **not** §5.2's tilt form `½ + (κ−½)·d`, which is a different parametrisation. §5.2 uses
the letter for both and this is the one place that matters.

⭐ **`disc` is `(2κ−1)²` with the sampling variance of `κ̂` subtracted and floored at 0**: a κ̂ within
`√σ²_d` of ½ is not composition signal. No constant is chosen — each half is the binomial `¼/N` plus that
arm's fitted overdispersion — so this is the derived answer to `TRAPS: a-licence-with-no-floor` and
`TRAPS: a-threshold-on-a-fitted-residue`, both of which refused a tuned floor.

⛔⛔ **AND IT HAS ONE CONSEQUENCE §5.2's "exactly zero at κ = ½" DOES NOT COVER: a gDNA-FREE library
switches the channel off at EVERY κ.** `N_gdna = 0` is guarded as `max(N_gdna, _EPS)` with
`_EPS = 1e-9`, so `1/N_gdna` becomes `1e9`, `σ²_d ≈ 2.5e8`, and since `(κ−½)² ≤ ¼` always, `disc = 0`
**even at κ = 0.99**. Measured 2026-08-17 on the 16-condition ladder: `tau_lam` is exactly **0.0 at all
70,176 slots on all four `g00` rows**.
⚠ **Read that as the derivation's own boundary case, not as a defect**: with no gDNA there is no
gDNA/RNA split for strand to speak about. What it *is* is a gap in this section — §5.2 alone predicts a
live channel on a stranded zero-gDNA library, and the code has never had one.

**5.3 gDNA's strand term is ½.** Double-stranded, no sense direction. A fitted mixture marginal was
implemented and **refuted**: the orientation discrimination is `(1−p)/p` and any *constant* for the
unstranded case cancels out of it, so an orientation-dependent marginal from a global genic average
destroys 78 % of the signal.

---

## 6. Overdispersion and second-moment evidence

Symmetric `Beta(a,a)` gives `od = 1/(2a+1)` (a=2 → 0.200, a=14 → 0.0345). Effective count
`n_eff = n/[1 + (n−1)·od]` — at od 0.2 a 1,523-fragment seed is worth **five coin flips**. Pooled moments:

    od = Σ_s[(k_s − n_s μ_s)² − n_s μ_s(1−μ_s)] / Σ_s n_s(n_s−1) μ_s(1−μ_s)

Its exact null information `I = (Σ n(n−1)pq)² / Σ[2n²p²q² + npq − 6np²q²]` collapses to the pair count
`Σn(n−1)/2` **only at mean ½**. ⭐ **Second-moment evidence is counted in PAIRS of fragments inside one
object — a singleton carries exactly zero.**

### 6a. ⭐⭐⭐ THE AWAY-HALF MOMENT — gDNA overdispersion with NO pure seed (2026-08-29)

The gDNA fit needs seeds whose RNA does not read as strand spread, and no structural class can be asserted
pure (pervasive transcription; the intergenic space is whatever the GTF leaves over). Orient each genic seed
so that RNA of its own gene pulls the residual DOWN, `d = (k − n/2)·sign(½ − κ)`; under pure gDNA `d` is
symmetric about 0 and the moment excess `d² − n/4` is EVEN in `d`, so the pooled moment restricted to the
away half has the null expectation of the full one:

    od = Σ_s a_s·(d_s² − n_s/4) / Σ_s a_s·n_s(n_s−1)/4        a_s = 1[d_s > 0] + ½·1[d_s = 0]

Unbiased for `ρ_g` under ANY distribution of RNA content across the seeds; a contaminated seed reaches the
away side only by noise, with small `d`, so it biases DOWN, never up. ⭐ **At `κ = ½` EXACTLY — reachable,
since `κ = (n_same+1)/(n_obs+2)` is exactly ½ whenever `2·n_same = n_obs`, the modal outcome on an unstranded
library — the orientation is degenerate and the FULL two-sided moment is used instead: RNA at the same mean ½
is symmetric, so it contributes only `n_g(n_g−1)/[N(N−1)] ≤ 1` of the excess and the one-sided guarantee
survives. Without that branch every residual collapses to 0 and the fit returns a hard `od = 0` — perfect
Binomiality, the most confident strand likelihood assertable.** The tie weight ½ is exact: at `n = 2`
under `BetaBinom(2, ½, ρ)`, `P(k=1) = ½(1−ρ)` and `P(k=0) = P(k=2) = (1+ρ)/4`, and the away half returns
exactly `ρ` while full tie weight returns `(3ρ−1)/(3−ρ)`. Half the pairs enter, so the information is half §6's pair count — `I = P/2` for the
TOTAL `P`, since `Var(e_s)|₀ = n(n−1)/8` and `E[a_s] = ½` give `Var(num) = P/8`, `E[den] = P/4`,
`Var(od_mom) = 2/P`. ⛔ NOT half the away half's OWN pair count: that halves twice and overstates the
standard error by √2. ⛔ Requires a gene strand to orient by — intergenic and AMBIG objects cannot enter — and
unannotated ANTISENSE RNA pushes toward the away side, the one recorded way to inflate it.

### 6b. ⭐⭐⭐ INFLUENCE WEIGHTING — why a deep seed is not worth its pair count (2026-08-30)

Pooling `od̂_s = (d_s² − n_s/4)/(n_s(n_s−1)/4)` by pair count is minimum-variance ONLY at `ρ = 0`. Given the
seed's latent rate `p` (`u = p − ½`), `E[od̂_s | p] = 4u²` **exactly** (since `E[d² − n/4 | p] = n(n−1)u²`), so
by the law of total variance

    Var(od̂_s | ρ) = V∞(ρ) + E_p[Var(od̂_s | p)] ≈ 2ρ²(1−ρ)/(1+2ρ) + 2/(n(n−1))

with the between-seed term EXACT from the symmetric Beta's moments (`E[u²] = ρ/4`, `E[u⁴] = 3ρ²/(16(1+2ρ))`;
MC-verified to 4 s.f. at ρ = 0.01/0.05/0.2/0.5). ⭐ **`V∞` does not depend on `n`** — a seed's information
about ρ SATURATES with depth — so the inverse-variance (Gauss–Markov) weights are

    w_s = 1/(½ + c_s·V∞(ρ))          c_s = n(n−1)/4

a constant at ρ = 0 (the pair-count estimator is this one with ρ pinned at 0) and equal-per-seed once
`c_s·V∞ ≫ ½`. ⭐ **No constant is introduced**, and the one approximation — the sampling term at `p = ½`
rather than integrated over `p` — cannot bias the fit, because the weights depend only on `n_s` and ρ and
never on a seed's own data, so the ratio has expectation ρ for ANY weight function. ρ enters its own weights,
and the root of `g(ρ) = clip(moment(ρ)) − ρ` is **bracketed by construction** (`g(0) ≥ 0`, `g(ceiling) ≤ 0`),
so bisection terminates with no iteration limit. Measured on real cfRNA, where one seed carried 77.8 % of a
library's pooled numerator: effective seeds 1.9 → 1241, and the fit came off the ceiling.

⭐ **THE MEAN MATTERS, AND THE TWO COMPONENTS DO NOT SHARE IT.** In general

    V∞(ρ, μ) = 3ρ²·[2ρ + μ(1−μ)(1 − 7ρ)] / [μ(1−μ)(1+ρ)(1+2ρ)] − ρ²
    b_s      = c_s·Var(od̂_s | ρ=0) = (2·n·pq + 1 − 6·pq)/(n − 1)      w_s = 1/(b_s + c_s·V∞(ρ, μ))

which reduces algebraically to `2ρ²(1−ρ)/(1+2ρ)` and `b_s = ½` at μ = ½ (verified against Monte Carlo at
ρ ∈ {0.01, 0.05, 0.2} × μ ∈ {0.5, 0.9, 0.99, 0.0023}). ⛔ **At a real library's κ = 0.0023 the same ρ = 0.05
gives V∞ = 0.285 against gDNA's 0.0043** — a seed at an extreme mean carries far less information about ρ than
its pair count suggests, so the two components may NEVER be compared on pair counts.

### 6c. ⭐⭐⭐ THE TWO COMPONENTS RECONCILE AGAINST EACH OTHER (2026-08-30)

A weighted estimator's precision is `Σ 1/V_s = Σ w_s·c_s` — its own weighted denominator — so each component
reports the precision of the estimate it actually made, at its own ρ and its own μ. The weaker then borrows
its DEFICIT from the better-measured one:

    od_w' = (I_w·od_w + (I_s − I_w)·od_s) / I_s          borrow weight (I_s − I_w)/I_s

0 when the two are equally informed (neither moves), 1 when the weak one measured nothing (it takes the
other's value outright). ⛔ NOT `(I_w·od_w + I_s·od_s)/(I_w + I_s)`, which drags an equally well-measured
component to the midpoint while the other stays put; and NOT a symmetric pooling, which would erase a real
difference. With NEITHER measured both take the ceiling: any common value leaves the strand channel
uninformative, since the composition term reads only the DIFFERENCE of the two dispersions.

---

## 7. Background rate, and deconvolving counts without strand

**7.1 A faint background rate is measurable only in aggregate.** `ρ_bg = Σg/ΣE`,
`Var(log ρ_bg) ≈ 1/Σg`. A region of effective length `E` resolves a rate only above ~`1/E` (Fisher
information = `ρ·E`), so no per-region estimator finds it and **true zero is resolved sharpest**.
⚠ **One-sided:** `ρ_bg > 0` proves DNA present; `ρ_bg ≈ 0` does NOT prove absence, because capture
depletes the off-target floor. Never a denominator or a scale.

**7.2 Counts against a gDNA background with NO strand data.** `P(g|C) ∝ P_bg(g)·1[0 ≤ g ≤ C]` with
`g ~ NegBinom(ρ_bg·E_g, α_eff)`, a flat one-sided prior on the RNA excess, truncated at the observed
total. `α = Σμ²/max(Σ(g−μ)² − Σμ, 0⁺)` (∞ ⇒ Poisson), `1/α_eff = 1/α + 1/(Σg + n₀)`. **No tuned
constant.**

---

## 8. Self-limiting damping

Treat a claim and the destination's own estimate as two studies of one quantity:

    b̂² = max(0, G² − v_msg − v_own)        ⇒     p_eff = 1 / max(v_msg, G² − v_own)

with `G` the log gap. Exact safety property: **a claim can outweigh a region's own belief only if it agrees
within `√2·σ_own`**, and where a region has no evidence `v_own = ∞` and the term switches itself off. No
knob.

---

## 9. Priors on a grid

`p(ρ_c) ∝ ρ_c^(c−1)` gives `Beta(c,c)`. `c = ½` (Jeffreys) is the only grid-width-stable choice — see
TRAPS: no-prior-means-haldane for what omitting the term does instead.

### 9a. ⛔⛔ WHY A SIMPLEX VERTEX IS UNREACHABLE WITHOUT EVIDENCE — and why that is not headroom

ψ lands 5–8 % short of `f_g = 1` on unexpressed genes, and that gap was once ranked as the top defect. It
is a **theorem, not a bug**:

* every PROPER prior on `[0,1]` has a median strictly INSIDE `(0,1)`;
* an object with zero composition evidence has posterior = prior;
* ⇒ a vertex is unreachable there **in any coordinate, at any depth**.

⛔⛔ **THE FIRST BULLET IS TRUE FOR PRIORS WITH A DENSITY AND IS FALSE IN GENERAL — §9a.1 states the
premise, which was never written down and has been read across this project as "no prior can reach the
vertex".** The theorem is kept; what changes is that "irreducible" turns out to be a property of the prior
FAMILY. Read §9a.1 before quoting the bullets or a vertex ceiling.

⭐ The empirical companion says the estimator is honest rather than merely stuck: measured per object,
`|f_g − truth| / sd(f_g)` has median **z = 0.5–0.6** on both simplex vertices, so every wrong answer sits
inside its own 1σ with a variance that is if anything conservative.

⛔⛔ **THE COMPANION IS SCOPED TO THE STRANDED CASE, AND AT `κ = ½` IT MAKES NO PREDICTION AT ALL**
(measured 2026-08-17). "Inside its own 1σ" carries an implied RATE: `sd(f_g)` shrinks like `n^(−1/2)`, so
the vertex shortfall should shrink with depth at the same rate. That rate is bought by the strand
channel, and §5.2 says the channel is identically zero at `κ = ½`. Measured, the shortfall's **log-log
slope against depth**:

| `κ` | slope | window it was measured on |
|---|---|---|
| the stranded arm | **−0.5221** | `n ≤ 50` fragments/slot |
| the real unstranded fit, `κ = 0.500369` | **−0.0000 in every decade** | 5.4 decades, 1.0 → 263,621 fragments/slot |

⭐ So `n^(−1/2)` is the **STRANDED** reading, and on unstranded data the shortfall is **depth-independent**:
there is no rate to quote and no depth at which it closes.
⛔ **Do not quote an `n^(−1/2)` shrinkage on an unstranded stratum, and do not read a flat shortfall there
as a regression.**

⚠⚠ **THE TWO ROWS ARE NOT ONE ARM SWEPT OVER `κ`, SO READ THE MEASUREMENT AND NOT AN ATTRIBUTION.** The
stranded row reproduces a clause that was established by DRIVING the solver on a synthetic pure-gDNA
object; the unstranded row is the shipped pipeline with the fitted landscape LIVE. More than one thing
differs between them. ⛔ In particular *"no strand channel ⇒ posterior = prior"* does **not** follow from
`κ = ½` on its own: `region_init.build_region_init` SUMS `τ_λ` over sources and
`density_deconv.density_factor_precision` is a second one — zero in the prior-free pass-0 the conclusion
below is about, and non-zero once a `λ`-factor is passed. **The depth-independence is measured; what
causes it is not settled here.**
⚠ The `κ = 0.500369` above is one fit on one panel; `SUCCESS.md` records **0.500689** from a different
run, and the point is the *order of magnitude of `κ̂ − ½`*, never either digit.
⛔ **Therefore the ceiling a vertex-pinning arm measures is the value of MISSING INFORMATION, not headroom
for a fix** — `scripts/design/vertex_ceiling.py` prices it and its docstring carries the number. ⛔ And
"fit a prior to fix it" is circular: pass-0 must stay prior-free, because its purpose is to produce the
substrate a prior is fitted ON. ⭐ The defect worth hunting is in the **confidently-wrong** population,
which is a different set of objects.
⚠ The one channel that could have supplied vertex evidence is closed independently: a certified-RNA count
of zero is consistent with `f_g = 1` too, gated by `tests/calibration/test_certified_rna_licence.py`.

### 9a.1 ⛔⛔⛔ THE PROOF NEEDS A CONTINUOUS CDF, AND A SPIKE-AND-SLAB HAS AN ATOM (2026-08-17)

⭐⭐ **THE THEOREM ABOVE IS KEPT AND IS CORRECT. WHAT IS ADDED IS ITS PREMISE.** The argument for *every
proper prior on `[0,1]` has a median strictly inside `(0,1)`* runs through the CDF: `F` is **continuous**
and **strictly increasing** on `(0,1)`, so `F(x) = ½` has a solution there and that solution is the
median. ⛔ **Both of those are properties of a DENSITY, not of properness.** For a general proper prior the
median is the quantile

    median  =  inf { x : F(x) ≥ ½ }

and an ATOM at `x₀` carrying mass `π ≥ ½` makes that infimum exactly `x₀`: `F` *jumps* across ½ at a single
point and there is nothing to solve. Three consequences, and the third is the one that matters:

* a **Beta** — §9c's family, atom-free, `F` continuous and strictly increasing — **cannot** put its median
  at `f_g = 0` or `f_g = 1`, at any strength, in any coordinate, ever. The bullets above hold verbatim;
* a **spike-and-slab** (§9d.4) is a proper prior on the same interval and **can**, whenever its spike
  carries at least half the mass;
* **no theorem is violated** — the premise simply does not hold for the second family.

⛔⛔ **THEREFORE "the value of MISSING INFORMATION", AND THE RECORDED READING THAT ABOUT THREE QUARTERS OF
THE VERTEX SHORTFALL IS IRREDUCIBLE, ARE STATEMENTS ABOUT THE PRIOR FAMILY AND NOT ABOUT THE DATA.**
`scripts/design/vertex_ceiling.py` prices the shortfall under the family the tree ships, which is the right
thing to price and is not headroom *within* that family — but a family with an atom moves the ceiling
itself, so **do not quote the shortfall as a property of the evidence.**

⛔ **This is not a licence to widen a bound, to soften the reference, or to fit a prior at pass-0** — all
three are refused above and stay refused. **The atom has to be EARNED**: §9d.4's spike is *located by a
measurement* (the off-target anchors, available before any solve) and is exactly `ρ_0 = 0` at the `g00`
zero control, which is why it reaches the vertex there rather than being placed there. An atom asserted
without a measurement behind it is the `m = σ(L)` mistake of §9c.1 in a new coordinate.

⭐⭐ **AND THE TWO SCOPES NOW ON THIS SECTION READ TOGETHER RATHER THAN COMPETING — they narrow different
things.** The `n^(−1/2)` companion is scoped to the STRANDED case: a statement about EVIDENCE and the rate
at which depth buys it. This one is scoped to priors with a DENSITY: a statement about the PRIOR. They meet
at the object that has neither — on unstranded data the shortfall is depth-independent (log-log slope
**−0.0000** over 5.4 decades, above), so posterior = prior there and **the prior family is the only thing
left that can reach a vertex**. ⛔ Neither scope licenses the other's claim: more depth still buys nothing
at `κ = ½`, and an atom still buys nothing where evidence exists, because it is worth one
pseudo-observation (§9c.1) and is swamped exactly as the Beta location is.

---

## 9b. ⭐⭐ THE EM's RNA PRIOR, AND WHY A SYNTHETIC NASCENT ENTITY GETS NONE

`assemble_priors` hands the EM two per-locus pseudocounts. The gDNA one lands additively on the single
gDNA component; the RNA one is shared among the RNA components in proportion to the evidence each
already carries (`em_solver.cpp:apply_grouped_prior_update`):

    rna_count      = Σ_{i ≠ g} raw[i]                       the whole RNA pool
    annotated_count = Σ_{i ≠ g, i not synthetic} raw[i]      the prior's eligible recipients

    out[g] = raw[g] + gdna_prior
    out[i] = raw[i] · (1 + rna_prior/annotated_count)        i annotated
    out[i] = raw[i]                                          i SYNTHETIC nascent

⭐ **Summed over the RNA components this is exactly `rna_count + rna_prior` either way**, so the
gDNA:RNA split is unchanged and the rule redistributes strictly WITHIN the RNA pool.

⛔⛔ **THAT IS A PER-M-STEP IDENTITY FOR GIVEN `raw` COUNTS — NOT AN END-TO-END INVARIANT.** The EM
iterates: a different `theta` gives a different E-step, hence different `raw`, hence a different
converged split. Measured on the 36-condition ladder, the rule moved gDNA by **+4.21 M fragments**
panel-wide. Three separate attempts to state this as "the library gDNA fraction cannot move" were wrong;
the sentence above is the correct one.

⚠ **The gate must name the denominator the chosen branch divides by.** It briefly tested
`annotated_count && annotated_carried` while the count branch divides by `annotated_count` alone, so a
locus with zero annotated count and nonzero carried alpha kept a live `rna_prior`, multiplied it by
`inv = 0` and silently dropped it — RNA summing to `rna_count` while gDNA kept `gdna_prior`. Reachable
under VBEM, the shipped default.

### Why the synthetic branch exists

A synthetic nascent entity is a shadow span the INDEX manufactured; no annotation asserts it exists. The
null hypothesis is therefore that it is **absent until the data proves otherwise**, which is a Dirichlet
`alpha = 0` on those components and needs no constant. Two consequences follow from the arithmetic:

* **A zero-count component cannot be revived by the prior**, because `out[i]` is proportional to
  `raw[i]`. So the coverage-weighted warm start is the ONLY spark a synthetic entity ever gets, and it
  must stay strictly positive — with `alpha = 0`, `theta = 0` is a fixed point.
* **A zombie decays geometrically** at `kappa/(1 + rna_prior/annotated_count)` per iteration, with
  `kappa = w_N/w_T < 1` for free: a shadow span is LONGER than the transcript it shadows, so its
  effective length already disfavours it.

⭐ **The survival criterion is derived, not chosen.** With `m` fragments whose only RNA candidate is the
entity, it grows iff `m·w_N > Total·theta_g·w_g` — *a nascent entity survives exactly when its
intron-exclusive evidence exceeds what gDNA alone would explain there.*

⛔ **The price, and it is strand-gated.** ~27 % of the freed mass goes to gDNA rather than annotated RNA
(86/14 at the weakest in-scope strandedness). On a nascent-ENRICHED toy it inverts to 83 % gDNA at
SS=0.65 and **zero at SS=1.0** — gDNA is strand-symmetric and RNA is not, so strand evidence alone
arbitrates the contested intronic fragment when it is strong enough.

⚠ `alpha = 0` is a blunt instrument — it cannot distinguish a nascent entity with real intronic support
from one with none — and under the NASCENT SCOPE RULING (`DESIGN.md` §0b) it is the RIGHT default rather
than a placeholder: nascent is absent until the data proves otherwise. ⛔ Building per-entity nascent
ACCURACY is designing around nascent and is not a 0.8.0 candidate; the per-transcript RNA prior lane
remains valuable for annotated transcripts, which is a different claim.

### 9b.1 ⭐⭐ THE PRIOR IS ALREADY AN ADDITIVE PER-COMPONENT PSEUDOCOUNT — the weights are the design

⚠ **This subsection was numbered `9c` until 2026-08-17, which collided with the ψ-reference `§9c` below —
two different sections carried one number, so a citation could not say which it meant.** It is renumbered
`9b.1` because it is a property of `§9b`'s rule and nothing else changed; `§9c` now means the Beta
reference alone, which is the sense `DESIGN.md` §6b and `ROADMAP.md` cite.

The rule above is habitually described as "multiplicative, hence neutral on the within-RNA split". That
is a description of one CHOICE OF WEIGHTS, not of a different kind of update. Writing `A` for
`annotated_count` and `P` for `rna_prior`:

    out[i] = raw[i]·(1 + P/A)  ==  raw[i] + P·raw[i]/A  ==  raw[i] + a_i,    Σ a_i = P

so the shipped prior is `a_i = P · w_i / Σ_eligible w` at **`w_i = raw[i]`** — an allocation in
proportion to the EM's own current belief. ⭐ **A prior that echoes the posterior carries no information,
which is exactly why it is neutral.** Any per-transcript scheme differs from it in the weights alone.

⛔⛔ **THE ONE PLACE THE WEIGHTS ARE NOT INTERCHANGEABLE IS `raw[i] = 0`, AND IT IS THE CONSEQUENTIAL
ONE.** With `w_i = raw[i]`, `out[i] = 0` is an ABSORBING STATE: no prior magnitude whatsoever revives a
component with no warm-start evidence, since `alpha` floors to `EM_LOG_EPSILON` and `digamma` of that is
`−1e300`. Any strictly positive weight has no such state. That absorbing state is the structural guard
against zombie revival named above, so **a weight vector that lifts every entity off zero removes it**.

⚠ **And the threshold that replaces it is not the one to design against.** With per-row max subtraction,
`EXP_CUTOFF = −708` and `digamma(a) ≈ −1/a`, a component's responsibility underflows to exactly zero
below **0.00142** alpha units. That point is real and is never binding: measured on the shipped native
EM, a component actually activates at **0.4709** alpha units when its candidates are equally likely and
**0.1625** at a 2-nat likelihood advantage — 331× and 115× higher — because the binding constraint is the
VBEM fixed point `alpha = Σ_u resp(alpha)`, not the exponential cutoff. ⛔ Separately, `assign_posteriors`
uses `log(theta + eps)` with `std::exp` and **no** cutoff, so a component whose E-step responsibility is
identically zero still collects posterior mass. Design against ~0.16–0.47, not against 0.0014.


## 9c. ⭐⭐⭐ ψ's COMPOSITION REFERENCE IS A Beta, AND ITS MEAN IS A THIRD TERM THE CODE DID NOT CARRY

`a·log f_g + b·log(1−f_g)` on the λ grid **is exactly Beta(a, b) in `f_g`** — verified numerically to six
decimals — so `a` and `b` are PSEUDO-COUNTS with a STRENGTH `a + b` and a MEAN `a/(a+b)`. The shipped
`_JEFFREYS_REF = ½` fixes BOTH: strength 1, which is the ignorance statement and is correct; and mean ½,
which **asserts the library is half gDNA**.

**The derivation.** Give the two components independent Gamma rate priors `ρ_c ~ Gamma(α_c, β_c)` with
Poisson counts `n_c ~ Poisson(ρ_c E_c)`. Then `X = ρ_g E_g ~ Gamma(a, s_g)` and `Y = ρ_r E_r ~ Gamma(b, s_r)`
with `s_c = β_c/E_c`. Writing `f = X/(X+Y)`, `T = X+Y` and integrating `T` out:

    p(f, T) ∝ f^(a−1)(1−f)^(b−1) · T^(a+b−1) · e^(−T(s_g f + s_r(1−f)))
    ⇒ p(f)  ∝ f^(a−1)(1−f)^(b−1) / ( s_g·f + s_r·(1−f) )^(a+b)

and on the solver's grid, where `|df/dλ| = f(1−f)`,

    log p(λ) = a·log f_g + b·log(1−f_g) − (a+b)·log( f_g + r·(1−f_g) ) ,   r = s_r/s_g

⛔⛔ **The shipped reference is exactly this with `s_g = s_r`** — the code has silently committed to the
two components having MATCHED SCALES, which is the "half gDNA" assertion appearing as a MISSING TERM
rather than as a chosen number. With `m` the object's prior expected composition
`ρ̄_g E_g / (ρ̄_g E_g + ρ̄_r E_r)` and `β_c = α_c/ρ̄_c`, the ratio is `r = (b/a)·m/(1−m)`, and at the
shipped `a = b = ½` the third term collapses to

    − log[ (1−m)·f_g + m·(1−f_g) ]

**Four properties, all structural.** (i) `m = ½` makes the bracket the constant ½, so the term drops — a
strict generalisation agreeing with the shipped constant exactly where its assumption is true.
(ii) ⭐⭐ **The tails stay `e^(−|λ|/2)` for every `m`** — in the `+λ` tail the bracket → `m(1−f_g)` and the
`log(1−f_g)` cancels one Jeffreys half — so **only the LOCATION moves and `L`-invariance is untouched**.
⛔ Moving `a`/`b` instead cannot do this: they set the tails and the location together, and `b = 0.03`
leaves 57 % of the mass outside `L = 10`. (iii) Proper for every `m ∈ (0,1)`. (iv) ⭐ Substituting
`u = f/(f + r(1−f)) ~ Beta(a,b)`, which is symmetric at `a = b`, gives **`median(f_g) = m` in closed
form** — the re-derivation of the relay's hard-coded `R_own = ½`.

### 9c.1 ⭐⭐⭐ THE SHIPPED LOCATION AND ITS STRENGTH — `L − log 2` NATS, AND NO CONSTANT IS CHOSEN

⛔ **THE CLAMP IS NOT A CHOSEN EPSILON, AND IT IS NOT A PER-OBJECT FLOOR EITHER.** Three candidates were
priced on the 16-condition ladder. A flat `eps = 1/Σeff_g` reads `m = 1 − 1e-8`, at which **99.1 % of the
reference's mass falls outside `L = 10`** and the answer becomes a function of the grid
(`TRAPS: a-clamp-at-the-closed-end-escapes-the-window`). The derived per-object floor
`m_i = E[g]_i/(E[g]_i + 1)` — one pseudo-fragment of RNA *here* — keeps 0.909–0.990 of the mass inside the
window but was measured **WORSE on every stratum** (0.609 / 1.045 / 0.580 / 1.000, against 0.381 / 0.659 /
0.363 / 0.800 for the hard clamp): on a structurally pure-gDNA object the truth IS `f_g = 1`, and a soft
floor pulls it back off the vertex it should sit on.

⭐⭐⭐ **THE STRENGTH IS A LOG-ODDS, AND THAT IS WHAT MAKES IT CHOOSABLE WITHOUT A CONSTANT.** The term is
`−log[(1−m)f_g + m(1−f_g)]`; at `f_g → 1` it is `−log(1−m)` and at `f_g → 0` it is `−log m`, so its full
range — how far it can move ψ — is

    strength  =  log( m / (1−m) )  =  logit(m)          ⇒   the claim's ODDS are e^strength

**The location written on the λ scale IS its strength in nats.** So the question "how strong?" is exactly
"what odds does the annotation entitle us to?", and Bayes answers it on the reference's own exponents: one
pseudo-observation of gDNA takes `Beta(a, b)` to `Beta(a+1, b)`, whose mean is

    m  =  (a + 1) / (a + b + 1)  =  0.75   exactly at   a = b = ½        (strength log 3 = 1.0986 nats)

⭐ A composition mean from pseudo-counts — same units as `m`, no new number, and it tracks `_JEFFREYS_REF`
automatically. ⚠ **Not `m = σ(a+b)`**, which an earlier draft used: `a + b` is a pseudo-COUNT and a strength
is a NAT, so equating them is an analogy rather than a derivation. Numerically they barely differ
(0.731 vs 0.75) and the measured optimum is BROAD — 0.69 → 1.50 nats within 6 % — which is what makes the
derived value safe rather than lucky.

⛔⛔ **WHAT THIS REPLACES, AND WHY THE LATTICE WAS THE WRONG PLACE TO GET IT.** *"A prior may not assert
more than the lattice can represent"* is a valid CAP — it is applied — but it was being used to CHOOSE, as
`m = σ(L)`. On the grid that saturates at

    strength  =  log[ (m² + (1−m)²) / (2m(1−m)) ]  =  L − log 2 + O(e^(−L)) ,   m = σ(L)

i.e. **9.31 nats ≈ 10,000:1** at `L = 10` — while the term's own docstring claimed "influence bounded by
one fragment". ⚠ It also breaks above `L ≈ 20.72`, where `1 − σ(L)` falls under the location clamp and the
strength sticks at `log(1/2ε) = 20.723` forever (measured 20.709 at `L = 25` against a claimed 24.307).

**The overturn depth, in FRAGMENTS.** One sense read on a `κ`-stranded library is worth `log(2κ)` nats, so

    n  =  strength / log(2κ)      fragments

which at `κ = 0.99` is **1.46** for the shipped `log 3` and was **13.6** for `σ(L)`. ⚠ The budget assumes
each fragment delivers `log(2κ)` nats ON THE λ AXIS, and near the vertex it does not — `I_strand ∝
[f_g(1−f_g)]²` collapses exactly where this prior points, so at `σ(L)` the residual reached **0.0058 at
5,000 fragments**, 370× the budget. ⭐ At `log 3` that residual never exceeds **0.006**, below one `K = 60`
grid step (0.085): a prior that cannot move the answer by a representable amount cannot outvote anything.
⛔⛔ At `κ = ½` the strand channel carries no composition information, `log(2κ) = 0`, and this budget
DIVERGES — there the refutation must come from the DENSITY mechanism (`fit_intron_background` +
`density_lambda_factor`), which measures `τ_fac = 161.4` at an intron and is what actually overturns the
claim on unstranded data.

⚠ **The strength is one pseudo-fragment on the EXPONENTS, so the location is swamped where evidence exists
and DECIDES where none does.** Measured, swinging `m` 0.02 → 0.98 on a single-strand object at `κ = 0.9`
moves `f_g` by 0.7700 at 6 fragments, 0.0090 at 120 and **0.0000** at 15,000 — but on an AMBIG object, or
once `κ = ½` kills the strand channel, by **0.95 at any depth**, because posterior = prior there.

### 9c.2 ⛔⛔⛔ WHY THE LOCATION MAY NOT ENTER `τ_λ` — the Jacobian, in one line

`I_strand ∝ [f_g(1−f_g)]² / (4p(1−p))` (§5.2). So when the location moves a slot's mode from
`f_g = 0.98576` to `0.99975`, the strand term falls by

    [0.98576·0.01424]² / [0.99975·0.00025]²  =  3,154×

against a MEASURED `tau_lam` fall of 3,227× — **~98 % of it is the Jacobian.** The likelihood is genuinely
that flat on λ at the point the prior chose; nothing was lost and there is nothing to restore.

⛔ And restoring it would not be a small increment. `own_composition_logvar` gives
`Var(log f_g) = (1−f_g)²/τ_λ`, which at the vertex is `~8e-08`, so `own_precision = 1/(Var + count_logvar)`
saturates at the COUNT ceiling: **τ = 0.029 and τ = 1e6 both return 850.44 against a ceiling of 850.50.**
Only `τ > 0` does any work, so feeding a prior in is a BOOLEAN gate flip that releases the whole count
precision — and, carrying no count, it does so at empty slots too (`n = 0` ⇒ `prec_g` 0 → `1/ψ'(½)` =
0.2026). `TRAPS: a-priors-curvature-is-not-the-datas-information`.

## 9d. ⭐⭐⭐ THE REFERENCE UNDER HYBRID CAPTURE — `ρ_g = ρ_0·ε`, AND WHY A PROBE PANEL MAKES `p(ε)` SPIKE-AND-SLAB

⛔ **STATUS, SO NOBODY GREPS FOR IT: THIS IS A DERIVATION AHEAD OF THE TREE.** Nothing in `src/` builds the
mixture below; what ships is §9c's Beta with §9c.1's strength and the structural per-object location
(`DESIGN.md` §6b.1). Everything here is derived, its three limits are checked, and every parameter it
needs is measurable off the shipped payload before any solve runs.

§9c gives the reference a per-object MEAN, `m_i = ρ_g,i·E_g,i / M_i`, and §9c.1 its strength.
⭐⭐ **EVERY TERM IN `m_i` IS EXACTLY KNOWN PER OBJECT EXCEPT `ρ_g,i`**, so the whole of the remaining work
is one question — *what is the gDNA density at THIS object?* This section answers it and changes nothing
else: the mean's form, the strength, the closed-end clamp (§9c.1) and the Jacobian ruling (§9c.2) all carry
over unchanged.

### 9d.1 ⛔⛔ AN EXON'S `ρ_g` IS NOT MEASURABLE — SO IT IS IMPUTED FROM NEIGHBOURS, AND THAT IS MESSAGE PASSING

gDNA is genomically continuous (AXIOM 0), so its density is **directly measurable exactly where no mature
transcript crosses**: intergenic REGIONs, intron REGIONs, and the BOUNDARIES against them (`exon|intron`,
`exon|intergenic`) — the `¬mrna_active` population, `DESIGN.md` §6b. At an EXON the unspliced mass is
`gDNA + unspliced RNA`, which **is the quantity being solved for**, so

    an exon's ρ_g cannot be measured at the exon — at any depth, in any coordinate

and the only remaining source is the objects either side of it along the chain:

    intron REGION  ↔  intron|exon BOUNDARY  ↔  exon REGION  ↔  exon|intron BOUNDARY  ↔  intron REGION

    intron REGION         deconvolve            → gDNA + unspliced RNA
    intron|exon BOUNDARY  unspliced             → gDNA + unspliced RNA ;  SPLICED → certified RNA
    exon REGION           the TARGET            → gDNA + RNA, unmeasurable alone

⛔⛔⛔ **IMPUTING A QUANTITY AT ONE OBJECT FROM ITS NEIGHBOURS ALONG A CHAIN IS MESSAGE PASSING, AND THERE
IS NO SECOND MECHANISM.** Any scheme that gives an exon a gDNA level from the objects around it *is* a
message, whatever it is called; re-deriving it under a new name is this project's most repeated wasted
turn. ⭐ **The mechanism is already BUILT**: `messages/relay.py`'s SPLICE IN is exactly the
BOUNDARY → EXON hop — only an EXON receives it; what it carries is a MEASUREMENT (a COUNT) and not an
imputation; it carries its own precision; its transfer variance is 0; and it is deliberately not `τ`-gated, so it
survives unstranded data where the strand channel is dead (§5.2). It sits behind
`CalibrationConfig.message_propagation = False`.
⚠ **The recorded price of switching it on is not evidence about the fixed relay**: it was measured with a
named, confirmed defect live — the composition licence checks transcript TERMINI (`terminus_flank_gain`)
and **not** `mrna_active` flipping, which is precisely the predicate saying the RNA population differs
across an exon↔intron hop, so a correct pure-gDNA claim is transported into the adjacent exon and drives it
to a confident wrong vertex (formerly a strict xfail in the reference-location gate file, which died
with the 2026-08-24 deletion — `DESIGN.md` §6b.1/§0c.2 carry the record) — and on the retired
36-condition ladder. **Re-price;
never inherit.**
⭐ **What this section supplies is the DISTRIBUTION the imputation is over, not a new transport.**

### 9d.2 ⛔⛔ WHY BOTH OBVIOUS IMPUTATIONS FAIL UNDER CAPTURE — the field is neither smooth nor scalar

Off capture the gDNA field is near-uniform, so *any* imputation works. Under capture it is a per-EXON,
nonlinear, **arbitrary** function of which probes the panel happens to contain. Measured on the
16-condition ladder at `g98 ss0.99` (2026-08-17), the four anchor classes:

| | intergenic REGION | intron REGION | `exon\|intergenic` | `exon\|intron` | span |
|---|---|---|---|---|---|
| capture OFF | 0.100521 | 0.100452 | 0.098343 | 0.098391 | **1.0×** |
| capture ON | 0.00418 | 0.00422 | 0.510 | 0.478 | **122×** |

All three imputations available before a solve therefore fail (16-condition ladder,
stranded × capture-ON, ratio to a `base` re-recorded in the same session, 2026-08-17):

| imputation | stranded × capture-ON | capture-OFF |
|---|---|---|
| a POOLED scalar reference | **3.90× WORSE** | — |
| the NEAREST anchor rung | 1.27× worse | — |
| LOCAL — `density_model.region_gdna_density` | 1.50× worse | **0.4037–0.4977, much BETTER than base** |

⭐ **The last row is the whole difficulty in one line: the same local form is the best thing available off
capture and a regression on it.** A single functional form cannot serve both, which is what forces a
MIXTURE rather than a better point estimate.

⚠ **And the anchor ladder is not a clean ladder.** Rungs 2 and 3 are PARTIALLY enriched — probes overlap
exons, and a fragment crossing an `exon|intron` BOUNDARY only partially overlaps a probe — so the in-gene
anchor under-reads a true exon's gDNA density by **2.6–3.6×**. It is a DETECTOR of enrichment, not a
calibrated level of it (`DESIGN.md` §6b.1).

⛔ **AND THE CAPTURE MODEL ALREADY IN THE TREE CANNOT BE BORROWED TO CLOSE THIS.**
`transcript_capture_eff_lengths(calibration: CalibrationResult, …)` consumes a **solved**
`CalibrationResult`, so using it to supply `ρ_g,i` to the reference that precedes the solve is circular.

### 9d.3 ⭐⭐⭐ THE TWO BOUNDS ON AN EXON — both from neighbours, off §3b's conserved identity

Write `M_i` for the object's conserved UNSPLICED mass (§3b's bank, which sums to ONE per fragment),
`E_g,i` / `E_r,i` for its two opportunities, `ρ_anchor,i` for the gDNA density measured at the adjacent
`exon|intron` flank, and `S` for the certified SPLICED count at the adjacent BOUNDARY — a spliced fragment
cannot be gDNA (AXIOM 0), so it is certified RNA and needs no deconvolution. §3b's conserved bank gives
each component's expected mass as `ρ_c·E_c`, and the RNA component's mass splits over exactly the two banks
it can land in:

    ρ_r·E_r  =  unspliced_RNA  +  S            ← §3b's conserved identity; SPLICED is a SEPARATE bank

⛔ **So `S` SUBTRACTS inside the RNA term rather than bounding `f_g` on its own.** `f_g ≤ 1 − S/M` is a
different statement and it is FALSE — a first draft wrote it and the origin-split truth violated it.

⚠⚠ **`E_r,i` IN THAT IDENTITY IS THE RNA's TOTAL CROSSING OPPORTUNITY — BOTH BANKS — AND IT IS NOT THE
`E_r` INSIDE §9c's `m_i`, WHICH IS THE UNSPLICED ONE.** The two differ by exactly `S`, and `f_g` is the
UNSPLICED fraction everywhere in this file; carrying one symbol for both is
`TRAPS: two-masks-one-name` committed in a derivation. ⭐ The algebra below is self-consistent in the
total-opportunity sense and lands back on the unspliced `f_g` exactly — that is what the `− S` does.
⚠ **And `ρ_r` is ONE ADJACENT BOUNDARY's own count, never a pooled RNA density**: `DESIGN.md` §6b's ruling
that RNA is the residual and is never predicted stands unweakened, and §6b carries the scope for the
capture case.

**LOWER — from the flank, because enrichment is MONOTONE in probe proximity.** The adjacent `exon|intron`
BOUNDARY sits at the EDGE of the probe footprint and the exon interior sits inside it, so `ρ_anchor,i` is a
lower bound on the exon's own `ρ_g,i` (this is the same 2.6–3.6× under-read of §9d.2, used in the one
direction it is sound). The gDNA mass is at least `ρ_anchor,i·E_g,i`, hence

    f_g  ≥  ρ_anchor,i · E_g,i / M_i

**UPPER — from the adjacent BOUNDARY's SPLICE-JUNCTION FLUX.** The unspliced pool is
`M_i = gDNA + unspliced_RNA`; substituting the conserved identity for `unspliced_RNA` gives

    f_g  =  1 − (ρ_r·E_r,i − S) / M_i          exact, if ρ_r were the exon's own true RNA density

The junction flux measures the RNA that actually splices through ONE adjacent junction, which is a SUBSET
of the transcripts crossing the exon, so it UNDER-states `ρ_r` and the equality relaxes to

    f_g  ≤  1 − (ρ_r·E_r,i − S) / M_i

Measured, mass-weighted over exon REGIONs against the origin-split oracle (16-condition ladder,
2026-08-17), the UPPER bound only:

| condition | true `f_g` | upper bound | mass in violation |
|---|---|---|---|
| `g00 ss0.50` off | 0.0000 | 0.6039 | **0.0 %** |
| `g50 ss0.50` off | 0.0627 | 0.6150 | 5.6 % |
| `g98 ss0.99` off | 0.7672 | 0.8332 | 19.5 % |
| `g50 ss0.99` ON | 0.5220 | 0.8272 | 9.8 % |
| `g98 ss0.99` ON | 0.9817 | **0.9918** | 4.6 % |

⭐⭐ **The bound is LOOSE where RNA is abundant and TIGHT where RNA is scarce**, which is the right way
round: it is worth most exactly where the composition question is hardest. At `g98 ss0.99` capture-ON the
endpoint ALONE carries an error of **0.0170** against the neutral ½'s **0.4817**.

⛔⛔ **IT IS A SUPPORT CONSTRAINT, NOT AN ESTIMATE** (`TRAPS: an-upper-bound-is-not-an-estimate`, measured
while refusing the soft-min per-transcript weight). It says where the answer CAN be; reading an endpoint as
the answer is the refused move, and the nonzero violation mass above is the honest record that both
premises — monotone probe proximity, and the junction being a subset — are not universal. ⭐ **The two are
COMPLEMENTARY**: the lower bound is tight at low `f_g` and fails under capture, the upper is tight at high
`f_g` under capture. A bound with a measured violation rate is a support with a tail, which is exactly what
the slab below is.

### 9d.4 ⭐⭐⭐ THE MIXTURE — an ATOM plus a SLAB, and NO NEW CONSTANT

Decompose the object's gDNA density into an un-enriched level and an enrichment,

    ρ_g,i  =  ρ_0 · ε_i ,        ε_i ≥ 1

A probe panel is a FINITE LIST of targets, so `ε` is not a smooth field with a scale: an exon is either in
the panel or it is not. **Unprobed ⇒ `ε = 1` exactly** — the same molecules, no pull-down — which is a
POINT MASS, not a narrow density. **Probed ⇒ `ε` is a large factor** set by probe count, overlap and
efficiency, about which the object itself says almost nothing — a broad SLAB. ⭐ **That is the physical
structure of a capture panel, not an analogy**, and it is what the 122× span and the partially-enriched
middle rungs of §9d.2 look like when read as a distribution rather than as a field. Hence, on `log ρ_g`:

    p(log ρ_g,i)  =  π · N( log ρ_0 , σ_0² )  +  (1 − π) · Unif[ log ρ_anchor,i , log ρ_max,i ]

    ρ_max,i  =  ( M_i − (ρ_r·E_r,i − S) ) / E_g,i    ← §9d.3's UPPER bound, in density coordinates
             →  M_i / E_g,i                          ← where no adjacent junction supplies ρ_r:
                                                       f_g = 1, the structural cap

⛔ **THE SLAB'S UPPER ENDPOINT IS §9d.3's BOUND, AND THE STRUCTURAL CAP IS ITS DEGENERATE CASE — the two
are not the same number and this is where it is reconciled.** *"An object cannot hold more gDNA than the
mass it holds"* gives `ρ_g ≤ M_i/E_g,i`, i.e. `f_g ≤ 1`, which is the endpoint whenever no adjacent
junction supplies `ρ_r`; where one does, §9d.3's flux term subtracts the RNA it proves is there and the
endpoint is strictly tighter. ⚠ **`DESIGN.md` §0c.3's table names this endpoint *"the object's OWN total
density"* while its own limit ③ quotes `0.9918` rather than `1.0000` — those are the two cases above, and
the ruling should say which it means.** Limit ③ below is the flux case.

| term | where it comes from | why it is not a new constant |
|---|---|---|
| `ρ_0` | the OFF-TARGET anchors — intergenic + intron REGIONs | MEASURED exactly, before any solve; the `¬mrna_active` population of §9d.1 |
| `σ_0²` | the same anchors' own spread | an empirical second moment of a measured population |
| slab LOWER | §9d.3's monotone flank bound `ρ_anchor,i` | derived from the geometry, per object |
| slab UPPER | §9d.3's upper bound `ρ_max,i` — the object's own total density LESS the RNA the adjacent junction proves is there, degenerating to that total density where there is no flux to subtract | structural at both ends: **an object cannot hold more gDNA than the mass it holds**, and the flux term is a NEIGHBOUR's own count (§9d.3) |
| the slab's SHAPE | uniform on `log ρ` | the scale-invariant (Jeffreys) prior for a positive RATE restricted to a known interval — the same ignorance statement §9 already makes with `c = ½`, written in the coordinate a rate has |
| `π`, the UNPROBED fraction | at pass-0, nothing has been observed about the probe indicator, and the reference's own Jeffreys exponents (`a = b = _JEFFREYS_REF`) give a symmetric Beta with mean and median **½** | the SAME convention §9c.1 uses to fix the strength — and `π` is MEASURABLE rather than assumed: the ratio of the off-target anchor to the in-gene one IS the enrichment (**0.98** without probes, **113–114** with, `DESIGN.md` §6b.1), so where probes exist it is read off the data, and where a panel's target list is supplied it is COUNTED |

**THE THREE LIMITS, each checked.**

**① capture-OFF ⇒ the shipped form, exactly.** The measured anchor ratio is **0.98**: `ε = 1` within
measurement, so there is no enriched population for the slab to describe, `π → 1`, and the mixture IS the
spike — a single measured density, which is the capture-OFF reference §9c and `DESIGN.md` §6b.1 already
validate. ⭐ This is the *"a strict generalisation agreeing with the shipped constant exactly where its
assumption is true"* property §9c's (i) demands of any reference term, met **by construction
rather than by tuning**.

**② `g00`, the ZERO-gDNA control ⇒ the vertex.** The anchors measure `ρ_0 = 0`, and the flank bound is 0
as well, so the spike and the slab's lower endpoint COINCIDE at 0. ⭐ `m` is a MEDIAN (§9c's (iv)) and a
median commutes with the monotone map `ρ_g ↦ ρ_g·E_g,i/M_i`, so §9a.1's atom argument carries from
`log ρ_g` to `f_g` unchanged: at the pass-0 `π = ½`, and at every `π > ½`, the median IS the atom and
`m_i → 0`. ⛔ **It is NOT "whatever `π` is", and that overstatement was in a draft of this section**: the
slab still spans up to `ρ_max,i` — 0.6039 in `f_g` at this very condition (§9d.3's table) — so a MINORITY
spike does not reach the vertex. That is §9a.1's `π ≥ ½` condition restated, not an extra assumption.
⛔ This is the limit that makes §9a.1 concrete: the reachable vertex is reached **because a measurement put
the atom there**, not because a family was chosen for its reach.
⚠ **AND IT IS THE ONE PLACE THE COORDINATE BREAKS**: `log ρ` is undefined at `ρ = 0`, so at the zero
control the slab must be taken in `ρ` rather than `log ρ`, or its endpoint clamped exactly as §9c.1 clamps
the location (`TRAPS: a-clamp-at-the-closed-end-escapes-the-window`). The LIMIT is well defined — both
components sit at 0 — the log coordinate is not.

**③ `g98` capture-ON ⇒ the slab's upper endpoint is informative on its own.** It reads **0.9918** against a
truth of **0.9817** (§9d.3's table), i.e. the support alone already excludes the neutral ½.

⛔ **WHAT THIS SECTION DOES NOT CHANGE, so it is not re-opened as a side effect.** §9c.1's strength is
untouched: a location is worth ONE pseudo-observation (`log 3` nats, overturned by 1.46 fragments at
`κ = 0.99`), atom or no atom, and it is swamped wherever evidence exists. §9c.2 still forbids the location
from entering `τ_λ` (`TRAPS: a-priors-curvature-is-not-the-datas-information`). And §9d.1's transport
question is separate from this section's distribution question — deriving one does not settle the other.

## 10. The second pass's score

⚠ **`f(L)` here is the SECOND PASS's per-fragment length term (`second_pass.py`'s `length_likelihood`
array), and the 0.8.0 deferral of the length COMPOSITION channel (§3d–§3e) does not touch it.** The two
share a word and nothing else: this one ranks one fragment's already-enumerated candidates against each
other, and never claims a composition.

    score  =  ρ × f(L) × s

normalised within one fragment's candidate set, factors applied in order of the evidence behind them and
**skipped when flat-zero among the survivors** (TRAPS: an-all-zero-factor-is-inert). One multinomial draw per fragment, in the
side buffer's canonical order, from a single RNG stream.
⛔ Never key the draw on the fragment's **content**: a content hash ties on exactly the duplicates it
would harm, so 100 identical fragments would draw identically and a 60/40 posterior would collapse to
100/0.

⚠ **Known approximation:** `ρ` enters as a hard multiplicative zero, but zero observations is
`P(0 | λ, E) = e^(−λE)`, not zero. The hard zero is the **large-exposure limit** of the correct
likelihood, so it is right where the library is deep and wrong where it is shallow.
