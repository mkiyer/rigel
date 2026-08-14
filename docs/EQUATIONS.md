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

⛔ **No measurements here.** Where a number appears it is a property of the *formula* (a limit, a
degenerate case, a worked value that makes the shape concrete), never a property of a library.

---

## 1. The deposit rule

**Notation.** A fragment is an interval `[s, s+w)` of *molecule* length `w`. A reference is region bound into
**regions** at every exon endpoint of every non-synthetic transcript; the 0-bp boundaries between adjacent
regions are **boundaries** (contiguous boundaries). `ell` is a region length.

**1.1 Fragment length `L` — ONE definition.** `L` = genomic span **minus region bound introns**. A paired-end mate
gap counts toward `L`; an intron does not. Whatever counts toward `L` must also count as coverage for
crossing, or the estimator is biased.

**1.2 Crossing count — the core identity.** A fragment `[s, s+w)` spans a point `p` iff
`s ∈ [p−w, p−1]` — exactly `w` start positions. So

    E[count at p]  =  ρ · Σ_w f(w)·w  =  ρ · mean_FL

exactly, for any fragment-length distribution, and **independent of both flanking region sizes**. This is
why the partitioning problem dissolves: a count at a 0-bp boundary depends on nothing but the density and the
length distribution.

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

## 2. Reciprocal length, and where it is model-free

**2.1 At a BOUNDARY it is exactly model-free.**

    E[Σ 1/L]  =  Σ_c Σ_w ρ_c·w·f_c(w)·(1/w)  =  Σ_c ρ_c

the opportunity factor `w` cancels identically, across a mixture with different component lengths.

**2.2 At a REGION it is not.** There the opportunity is `(ell − w + 1)₊`, which `1/w` does not cancel;
`Σ1/L` is then only a better-conditioned second moment. ⛔ **Do not claim region-level model-freeness.**

**2.3 And it fails at a terminus, exactly.** `E[Σ1/L] = ρ · E_f[placements(w)/w]`, which equals `ρ` only
where placements ∝ `w`. At a point 50 bases from an end, placements = 50 for **every** `w > 51`,
independent of `w`. So `Σ1/L` fixes the length bias and a reach taper fixes the placement loss; **neither
substitutes for the other.**

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
rather than per-object. `node_geometry.terminus_flank_gain`.

⛔⛔ **WRITE IT IN GENOMIC TERMS, NEVER IN TSS/TES.** TSS/TES is transcript-relative and the strand flips
it:

    the RIGHT flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic LOW  end is there  ⟺  TSS₊ or TES₋
    the LEFT  flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic HIGH end is there  ⟺  TES₊ or TSS₋

⛔ **TERMINI ONLY.** A DONOR/ACCEPTOR BOUNDARY also changes the population, but there the flux is **measured**
(`junction_count`) and the graft and the peel exist to route it. A terminus has no flux to measure: a
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
2026-08-05; derived and gated the same day. `test_splice_flux_reframe`, `node_total_density`.)

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

✅ **This WITHDRAWS the sign correction this section used to carry.** `bp_solver`'s P1d asserts
`rho_R(exon) ≥ rho_nas(B) + rho_mat(B)`, a LOWER bound, and that is **right**: the measured ratio of
**1.103** (no nascent) and **1.049** (with) is `1 + (1−s)(k−1)` with `k` the frame gap above and `s` the
nascent share of the exon's RNA — the nascent arm is measured in the exon's own frame and dilutes it,
which is precisely what identifies the inflation as `rho_mat`'s rather than the bound's. ⚠ What does
survive is that `graft_premise_logvar` is fitted on fluxes carrying the same inflation, so it inherits it.

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
gapped `flgap` pair is where this bias is measured; the ladder is where calibration is measured.

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
`scripts/design/prior_yardstick.py`.

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
for RNA that spliced elsewhere — so it errs in one direction for one sub-population. ⛔ Every ladder
condition is `nrna_none`, i.e. the substrate most favourable to the tighter bound and the least
representative; a repair validated there is not validated.
⚠ `length_channel_census.py` regenerated every number here and `short_exon_fl_probe.py` priced the
tighter bound; both were deleted in `b7ed7a0b` and are recoverable from `f470a570`.

## 4. Opportunity corrections for length pools

A pool is a length-dependent **selection**, so its raw histogram is not the library's length
distribution. Let `A(w)` be the pool's opportunity and `T(w)` the total opportunity for the same
population.

**4.1 The divisor is a probability.**

    pi(w) = A(w) / T(w)          fitted(w) = count(w) / pi(w)

⛔ Never `count(w)/A(w)` — see TRAPS: divide-by-a-probability.

**4.2 The sj pool** (`calibration/junction_opportunity.py`). For a transcript with exon lengths
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

⭐ The empirical companion says the estimator is honest rather than merely stuck: measured per object,
`|f_g − truth| / sd(f_g)` has median **z = 0.5–0.6** on both simplex vertices, so every wrong answer sits
inside its own 1σ with a variance that is if anything conservative.
⛔ **Therefore the ceiling a vertex-pinning arm measures is the value of MISSING INFORMATION, not headroom
for a fix** — `scripts/design/vertex_ceiling.py` prices it and its docstring carries the number. ⛔ And
"fit a prior to fix it" is circular: pass-0 must stay prior-free, because its purpose is to produce the
substrate a prior is fitted ON. ⭐ The defect worth hunting is in the **confidently-wrong** population,
which is a different set of objects.
⚠ The one channel that could have supplied vertex evidence is closed independently: a certified-RNA count
of zero is consistent with `f_g = 1` too, gated by `tests/calibration/test_certified_rna_licence.py`.

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

⚠ `alpha = 0` is the correct null but a blunt instrument: it cannot distinguish a nascent entity with
real intronic support from one with none. The successor is a per-transcript allocation of the RNA prior.

### 9c. ⭐⭐ THE PRIOR IS ALREADY AN ADDITIVE PER-COMPONENT PSEUDOCOUNT — the weights are the design

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
