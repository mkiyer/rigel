# EQUATIONS — the derivations the code depends on

**Every formula here is implemented somewhere and gated by a test.** This file exists so the derivation
does not have to be reconstructed from the code, and so two modules cannot come to disagree about one
quantity (`TRAPS.md` E13).

⛔ **No measurements here.** Where a number appears it is a property of the *formula* (a limit, a
degenerate case, a worked value that makes the shape concrete), never a property of a library.

---

## 1. The deposit rule

**Notation.** A fragment is an interval `[s, s+w)` of *molecule* length `w`. A reference is cut into
**nodes** at every exon endpoint of every non-synthetic transcript; the 0-bp boundaries between adjacent
nodes are **lines** (contiguous edges). `ell` is a node length.

**1.1 Fragment length `L` — ONE definition.** `L` = genomic span **minus cut introns**. A paired-end mate
gap counts toward `L`; an intron does not. Whatever counts toward `L` must also count as coverage for
crossing, or the estimator is biased.

**1.2 Crossing count — the core identity.** A fragment `[s, s+w)` spans a point `p` iff
`s ∈ [p−w, p−1]` — exactly `w` start positions. So

    E[count at p]  =  ρ · Σ_w f(w)·w  =  ρ · mean_FL

exactly, for any fragment-length distribution, and **independent of both flanking node sizes**. This is
why the partitioning problem dissolves: a count at a 0-bp line depends on nothing but the density and the
length distribution.

**1.3 Contained opportunity.** The starts at which a length-`w` fragment fits wholly inside a node of
length `ell`:

    (ell − w + 1)₊

and over a set of nodes, `A(w) = Σ_n (ell_n − w + 1)₊`. Its fragment-length expectation is
`E_f[(L−w+1)₊] = (L+1)·F(L) − S(L)` with `F` the FL CDF and `S(L) = Σ_{w≤L} w·f(w)`; beyond the support it
is `L + 1 − mean_FL`.
⛔ **The `+1` is the discrete count of start positions, not a fudge** — drop it and the divisor is exactly
0 when a node is one fragment long.

**1.4 Crossing-exactly-one-line opportunity.** For a line with flanking node lengths `a` (left) and `b`
(right):

    A(w)  =  (w−1)₊  −  (w−1−a)₊  −  (w−1−b)₊  +  (w−1−a−b)₊

⭐ **The two nearest lines are the only ones that need excluding**, so this is exact rather than a
truncation: a fragment is an interval containing the line, so if it reaches any line beyond `p−a` it must
also cross `p−a`. ⚠ Reference ends need no special case — the partition cuts at `0` and `L_ref`, so the
outermost node's length *is* the distance to the wall.

**1.5 General crossing divisor — ONE formula for both edge kinds.** With `R_lo`, `R_hi` the molecule's own
remaining sequence either side:

    E_J  =  E_f[ min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)₊ ]

Mean fragment length is its **large-reach limit**, not a separate case. On RNA N(200,50): 199.0 at R=550
(so mid-transcript junctions are exact), 160.1 at 200, 87.8 at 147, 19.6 at 100, 50.0 at R=50 — a **4×
error** if mean length is used blindly at a first exon.

**1.6 Reach.** At position `p` inside exon `e`: `reach_lo = exonic_bases_before(e) + (p − e.start)`,
`reach_hi = total_exonic − reach_lo`; maximised over transcripts independently per side AND per strand.
gDNA is unbounded on a chromosome; mature RNA stops at the polyA site. **A reach of 0 is meaningful, not a
sentinel.**

**1.7 Partition rule.** `cuts = unique({exon.start, exon.end over all non-synthetic transcripts} ∪
{0, ref_length})`, node `i` = `[cuts[i], cuts[i+1])`, **no merging**. Introns and termini are already exon
endpoints, so only the merge step was removed. Edges always run `src < dst`, so **genomic order is a
topological order** and there is no graph traversal anywhere.

---

## 2. Reciprocal length, and where it is model-free

**2.1 At an EDGE it is exactly model-free.**

    E[Σ 1/L]  =  Σ_c Σ_w ρ_c·w·f_c(w)·(1/w)  =  Σ_c ρ_c

the opportunity factor `w` cancels identically, across a mixture with different component lengths.

**2.2 At a NODE it is not.** There the opportunity is `(ell − w + 1)₊`, which `1/w` does not cancel;
`Σ1/L` is then only a better-conditioned second moment. ⛔ **Do not claim node-level model-freeness.**

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
carries exactly zero information at any depth (`TRAPS.md` F3).

**3.2 Density is the frame-invariant currency; a fraction is not.** `ρ_c = C_c/E_c` agrees across the
contained / crossing / spliced frames to ~0.36 % and does not degrade across 1×/30×/300× capture. The
log-odds shifts between frames by exactly `log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`, and **capture
cancels identically** — a ratio transports across a capture cliff, an absolute density does not.

**3.3 Conservation with unequal effective lengths is `Σ_c ρ_c·E_c = M`**, not `Σ_c ρ_c = M/E`. Its
sensitivity is bounded: a purely compositional error moves `M/Σρ_c E_c` by only ×1.04 on a contained node
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
`ρ_tot(dst)`, the destination's **own** total, independently of the source's measurement — `TRAPS.md` D4/D4b.

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

⭐ The population conjunct is set algebra, not a model. An EDGE counts what spans it **contiguously**, so

    T(EDGE)  =  T(NODE_left) ∩ T(NODE_right)

and `T(EDGE) = T(right)` fails iff a transcript's body **begins** at the EDGE, `T(EDGE) = T(left)` iff one
**ends** there. So a transcript **terminus** is exactly what makes one flank's population larger, the test
is an **equality** per `(EDGE, side)` pair (both directions of mismatch corrupt `φ_g`), and it is per-step
rather than per-object. `node_geometry.terminus_flank_gain`.

⛔⛔ **WRITE IT IN GENOMIC TERMS, NEVER IN TSS/TES.** TSS/TES is transcript-relative and the strand flips
it:

    the RIGHT flank gains RNA at an EDGE  ⟺  a transcript's genomic LOW  end is there  ⟺  TSS₊ or TES₋
    the LEFT  flank gains RNA at an EDGE  ⟺  a transcript's genomic HIGH end is there  ⟺  TES₊ or TSS₋

⛔ **TERMINI ONLY.** A DONOR/ACCEPTOR EDGE also changes the population, but there the flux is **measured**
(`junction_count`) and the graft and the peel exist to route it. A terminus has no flux to measure: a
transcript simply begins. That is the derived boundary between the two treatments.

**3.5c ⭐⭐ THE MASS PIN CARRIES THE SAME LICENCE, PLUS ONE STRUCTURAL CASE.** `Σ_c ρ_c·E_c = M` is an
identity *under the imputation premise*, restored by `k = M/S` with `S` filling every unsupplied component
from the destination's **own** density. D4 permits a message to use the destination's constants and
observations, never its beliefs — so the pin is licensed in exactly the two states where no belief reaches
`S`: **(i)** the message supplied the composition (§3.5b, nothing is filled in), or **(ii)** the
destination is a **structurally pure-gDNA object**, where there is no unsupplied component and `f_g = 1` is
structure, so `S = ρ_g·E_g` and the pin hands the object its own **measured** `M/E_g`.
⭐ Case (ii) is how the capture landscape reaches an exon at all — a G1 EDGE has `prec_g = 0` and cannot
originate a level through the fuse, so the pin is its only channel. Unlicensed, the pin was D4 at full
strength: `k = 1/(φ_msg + R_own)`, fixed point `(1−R_own)·ρ_tot`, and `R_own` is **exactly ½** at a slot
with no composition evidence — so it drove the delivered gDNA *fraction* to ½ regardless of the truth.

and for **gDNA the level crosses UNSCALED**, because gDNA is uniform along the genome *before* capture.
⭐⭐ **This is one rule for both capture arms with no branch on capture**, because the capture landscape is
carried by the structurally-pure-gDNA population's **own measurements**, not by a scale factor: the relay's
mass pin re-anchors the running level at every object, and at a pure-gDNA object that total *is* its gDNA
density at its own capture stratum. A per-capture-class landscape ratio built to do this explicitly measured
**byte-identical off capture** and 1.2 % of one class on one capture-ON condition, and was deleted.

⚠ **One residual is separate and is not this**, so a reader does not attribute it here: the level an exon
receives is the one measured at its flanking EDGE, and a fragment spanning an EDGE at a gene end lies only
PARTLY under the capture probe while one contained in the exon lies wholly under it — so that level is a
**lower bound** under capture, and closing the gap needs a probe-geometry model the tool does not have.
⭐ The second residual that used to be listed here — the mass pin filling unsupplied components from the
destination's own `f_g ≈ ½` — is **fixed**: §3.5c.

---

## 4. Opportunity corrections for length pools

A pool is a length-dependent **selection**, so its raw histogram is not the library's length
distribution. Let `A(w)` be the pool's opportunity and `T(w)` the total opportunity for the same
population.

**4.1 The divisor is a probability.**

    pi(w) = A(w) / T(w)          fitted(w) = count(w) / pi(w)

⛔ Never `count(w)/A(w)` — see `TRAPS.md` C3.

**4.2 The junction pool** (`calibration/junction_opportunity.py`). For a transcript with exon lengths
`e_1..e_K` and total `L = Σ e_i`, the starts at which a length-`w` window crosses **at least one**
junction:

    A_j(w)  =  (L − w + 1)₊  −  Σ_i (e_i − w + 1)₊

⭐ Derived via the **complement** — a window crosses no junction iff it lies wholly inside one exon, and
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
(`TRAPS.md` C4).

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

with `G` the log gap. Exact safety property: **a claim can outweigh a node's own belief only if it agrees
within `√2·σ_own`**, and where a node has no evidence `v_own = ∞` and the term switches itself off. No
knob.

---

## 9. Priors on a grid

`p(ρ_c) ∝ ρ_c^(c−1)` gives `Beta(c,c)`. `c = ½` (Jeffreys) is the only grid-width-stable choice — see
`TRAPS.md` D5 for what omitting the term does instead.

---

## 10. The second pass's score

    score  =  ρ × f(L) × s

normalised within one fragment's candidate set, factors applied in order of the evidence behind them and
**skipped when flat-zero among the survivors** (`TRAPS.md` D7). One multinomial draw per fragment, in the
side buffer's canonical order, from a single RNG stream.
⛔ Never key the draw on the fragment's **content**: a content hash ties on exactly the duplicates it
would harm, so 100 identical fragments would draw identically and a 60/40 posterior would collapse to
100/0.

⚠ **Known approximation:** `ρ` enters as a hard multiplicative zero, but zero observations is
`P(0 | λ, E) = e^(−λE)`, not zero. The hard zero is the **large-exposure limit** of the correct
likelihood, so it is right where the library is deep and wrong where it is shallow.
