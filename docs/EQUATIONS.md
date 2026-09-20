# EQUATIONS — the derivations the code depends on

**Purpose.** Every derivation here is implemented by a named module and function in `src/rigel/`, or is
the maths behind a ruling in `DESIGN.md` that the code obeys, and each section opens by naming what
implements it. The file exists so a derivation does not have to be reconstructed from the code, and so
two modules cannot come to disagree about one quantity (`TRAPS: two-docstrings-one-quantity`).

**What does not belong here.** Panel measurements — a derivation is maths and its assumptions; a number
that appears is a property of the formula (a limit, a degenerate case, a worked value), and the
instrument that verified it is named in one sentence. Rulings and their measured prices live in
`DESIGN.md`, open problems and refusals in `ISSUES.md`, lessons in `TRAPS.md`, the panels in
`TESTING.md`. Derivations for mechanisms not in the tree do not live here either: the fragment-length
composition channel is deferred past 0.8.0 (`DESIGN.md` §0b); the length model as an opportunity and a
divisor (§1, §3.6b, §4) is a different, live thing. Section numbers are anchors cited from the code and
the other docs, so a deleted section leaves a gap.

---

## 1. The deposit rule

Implemented by `native/calibration/accumulator.cpp` (the deposits), `calibration/region_geometry.py`
(reach, the partition) and `calibration/effective_length.py` (the opportunities). The executable
specification is `tests/native/_accumulator_reference.py`, which wins over this text.

**Notation.** A fragment is an interval `[s, s+w)` of *molecule* length `w`. A reference is cut into
**regions** at every exon endpoint of every non-synthetic transcript; the 0-bp cuts between adjacent
regions are **boundaries**. `ell` is a region length.

**1.1 Fragment length `L` — one definition.** `L` = genomic span minus the introns it spans. A paired-end
mate gap counts toward `L`; an intron does not. Whatever counts toward `L` must also count as coverage for
crossing, or the estimator is biased.

**1.2 Crossing count — the core identity.** A fragment `[s, s+w)` crosses an inter-base point `p` iff
`s < p < s+w`, i.e. `s ∈ [p−w+1, p−1]` — exactly `w − 1` start positions (the shipped `1/(length−1)`
deposit is its reciprocal). So

    E[count at p]  =  ρ · Σ_w f(w)·(w−1)  =  ρ · (mean_FL − 1)

exactly, for any fragment-length distribution, and independent of both flanking region sizes. This is
why the partitioning problem dissolves: a count at a 0-bp boundary depends on nothing but the density and
the length distribution.

**1.3 Contained opportunity.** The starts at which a length-`w` fragment fits wholly inside a region of
length `ell`:

    (ell − w + 1)₊

and over a set of regions, `A(w) = Σ_n (ell_n − w + 1)₊`. Its fragment-length expectation is
`E_f[(L−w+1)₊] = (L+1)·F(L) − S(L)` with `F` the FL CDF and `S(L) = Σ_{w≤L} w·f(w)`; beyond the support it
is `L + 1 − mean_FL`. ⛔ The `+1` is the discrete count of start positions, not a fudge — drop it and the
divisor is exactly 0 when a region is one fragment long.

**1.4 Crossing-exactly-one-line opportunity.** For a boundary with flanking region lengths `a` (left) and
`b` (right):

    A(w)  =  (w−1)₊  −  (w−1−a)₊  −  (w−1−b)₊  +  (w−1−a−b)₊

The two nearest boundaries are the only ones that need excluding, so this is exact rather than a
truncation: a fragment is an interval containing the boundary, so if it reaches any boundary beyond `p−a`
it must also cross `p−a`. Reference ends need no special case — the partition cuts at `0` and `L_ref`,
so the outermost region's length *is* the distance to the wall.

**1.5 General crossing divisor — one formula for both boundary kinds.** With `R_lo`, `R_hi` the molecule's
own remaining sequence either side:

    E_J  =  E_f[ min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)₊ ]

Mean fragment length is its large-reach limit, not a separate case. On RNA N(200,50): 199.0 at R=550
(so mid-transcript sj are exact), 160.1 at 200, 87.8 at 147, 19.6 at 100, 50.0 at R=50 — a 4× error if
mean length is used blindly at a first exon.

**1.6 Reach.** At position `p` inside exon `e`: `reach_lo = exonic_bases_before(e) + (p − e.start)`,
`reach_hi = total_exonic − reach_lo`; maximised over transcripts independently per side and per strand.
gDNA is unbounded on a chromosome; mature RNA stops at the polyA site. A reach of 0 is meaningful, not a
sentinel.

**1.7 Partition rule.** `region bounds = unique({exon.start, exon.end over all non-synthetic transcripts}
∪ {0, ref_length})`, region `i` = `[region bounds[i], region bounds[i+1])`, no merging. Introns and termini
are already exon endpoints. Boundaries always run `src < dst`, so genomic order is a topological order and
there is no graph traversal anywhere.

---

## 2. Reciprocal opportunity, and where it is model-free

Implemented by `accumulator.cpp` (the `1/A(w)` deposits), `calibration/substrate.py` (the banks) and
`calibration/total_abundance.py` (the start/end pair of §2.3b). The executable statement of the support
factor is gated in `tests/native/test_conserved_mass.py`.

**The general rule: deposit `1/A(w)` where `A(w)` is that population's own opportunity** — the number of
admissible start positions for a length-`w` fragment at that object. Then, by linearity alone,

    E[Σ 1/A]  =  Σ_c ρ_c · P_c(A > 0)

⛔ The cancellation is conditional on its own support
(`TRAPS: a-cancellation-is-conditional-on-its-support`): a fragment with `A(w) = 0` deposits nothing, so
the factor `P(A > 0)` survives — a functional of exactly the pmf the channel claims independence from.

**2.1 At a boundary it is model-free for every real library.** The crossing opportunity is
`A(w) = (w−1)₊` (§1.2), the shipped deposit `1/(L−1)`:

    E[Σ 1/(w−1)]  =  Σ_c ρ_c · P_c(w ≥ 2)  =  Σ_c ρ_c      whenever frag_min ≥ 2

so the support factor is exactly 1 for any library whose fragments are at least 2 bp — unconditionally
robust across a mixture with different component lengths. The sj bank takes the identical deposit.

**2.2 At a region it is model-free only within its support — a density shape, not a total.** The
contained opportunity is `A(w) = (ell − w + 1)₊`, the shipped deposit `1/(ell − w + 1)`
(`region_contained_inv_opportunity_sum`), and

    E[Σ 1/A]  =  ρ · P(w ≤ ell)          ← NOT ρ

`P(w ≤ ell)` is a pmf functional and differs per component, so at region scale the gDNA-vs-RNA circularity
this channel exists to remove moves out of the divisor and into the support (exactly zero below
`ell < frag_min`). ⛔ Do not read the region bank as a model-free total; within a fixed `ell` band it is a
valid density shape.

**2.3 And it fails at a terminus, exactly.** `E[Σ1/L] = ρ · E_f[placements(w)/w]`, which equals `ρ` only
where placements ∝ `w`. At a point 50 bases from an end, placements = 50 for every `w > 51`, independent
of `w`. So `Σ1/L` fixes the length bias and a reach taper fixes the placement loss; neither substitutes
for the other.

**2.3b The start/end relation is model-free at a region, and the wall is the only exception.** A
fragment's first covered base falls in region `r` iff its start lies in `r`, so the opportunity is `ℓ` —
the same number for every fragment length, which is exactly what the contained relation is not:

    E[S_r]  =  ρ · ℓ            for every w        ⟸ no pmf functional, no support factor
    E[C_r^inv]  =  ρ · P(w ≤ ℓ)                    ⟸ §2.2, truncated per component

so `S_r/ℓ` is a total where the contained bank is a shape. The one exception is the template wall: within
`d` of the template's genomic-high end, only starts that still admit a length-`w` molecule count, so

    A_start(w | d)  =  min( ℓ , (d + ℓ − w + 1)₊ )              ← w-dependent again
    A_start(w | d)  =  ℓ    for every w    ⟺    d ≥ w_max − 1   ← the exactness condition

and at `d = 0` this is `(ℓ − w + 1)₊`, the contained opportunity, still deposited as the flat `1/ℓ`.
`E_r` is the mirror, exact iff `d_low ≥ w_max − 1`. The two fail at opposite ends, which closes the pair:
use the side whose wall does not bind, average where both are exact (two counts of one rate at one
opportunity, so the pooled rate is the precision-weighted combination), refuse where both bind. The wall
is component-differential — gDNA's template is the chromosome and never binds, a nascent molecule's is
its genomic span, a mature molecule's its spliced length — so the distance is taken at the component
minimum over the populations `T(slot)` admits. The pair also gives a field-free test: `S_r/ℓ` and `E_r/ℓ`
share their opportunity and read the same field, so their ratio has expectation exactly 1 wherever both
are exact, and a binding wall moves it in a known direction. ⛔ A
comparison against the contained bank is not field-free: the two weight a region's positions differently
(`TRAPS: two-estimators-of-one-rate-weight-the-field-differently`).

**2.4 Superadditivity.** Contained effective lengths are superadditive — `Σ E(children) < E(whole)` for
any split — so densities, effective lengths and variances cannot be pooled across two partitions; only
truth pools additively (`f_g = ΣG/(ΣG+ΣR)`).

---

## 3. The two-component deconvolution

**3.1 The 2×2 at one object.** Stored by `calibration/substrate.py` (the `inv_length_sum` banks); no
calibration module solves it, because fragment length as a composition channel is deferred past 0.8.0
(`DESIGN.md` §0b). Kept because `SUCCESS.md`'s information census is defined on it.

    N       =  ρ_g·E_g[w]  +  ρ_r·E_r[w]
    Σ 1/L   =  ρ_g         +  ρ_r

Identified iff the two mean lengths differ (the second row being literally `[1,1]` is what conditions
it); `N / Σ(1/L)` is the abundance-weighted mean fragment length. The identifying quantity is the gap
`μ_g − μ_r`: at equal means the channel carries zero information at any depth
(`TRAPS: equal-lengths-carry-no-composition`), which is why the ladder gives gDNA and RNA equal lengths —
a gap would let the EM split the origins on length alone and mask calibration (`DESIGN.md` §0b).

**3.2 Density is the frame-invariant currency; a fraction is not.** `ρ_c = C_c/E_c` agrees across the
contained / crossing / spliced frames. The log-odds shifts between frames by exactly
`log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`, and capture cancels identically — a ratio transports
across a capture cliff, an absolute density does not.

**3.3 Conservation with unequal effective lengths is `Σ_c ρ_c·E_c = M`**, not `Σ_c ρ_c = M/E`
(`density_deconv.py`, `messages/transfer.py`). Its sensitivity is bounded: a purely compositional error
moves `M/Σρ_c E_c` by only ×1.04 on a contained region and ×1.50 at a crossing, so a large violation is
accumulated drift, never one hop.

**3.4 Why an integer count must be stored** (``count_logvar`` in `native/transfer_rows.h`, the one home of the
counting term). `Var(log ρ_c) = 1/(f_c·n) ≡ Var(log f_c) + 1/n`, exactly. Mass sums fractional
per-fragment shares, so `1/mass` is not a counting variance. The shipped counting term is
`trigamma(n + ½)` — the Jeffreys posterior's exact `Var(log ρ)` at every count including zero, which is
`1/n` from `n ≈ 10` and `π²/2` at `n = 0`: a zero count is a measurement, not an absence.

**3.5 A composition crosses only where the population is shared; a gDNA level crosses unscaled.**
`native/transfer_kernel.h`: the face maps (`splice_faces`, `terminus_rules`) carry a composition under
§3.5b's licence; `gdna_lane` carries the gDNA level everywhere else. Why the level needs its
own lane is one substitution. Rescaling a source's density by the ratio of totals
`r = ρ_tot(dst)/ρ_tot(src)` and writing `ρ_c(src) = φ_c(src)·ρ_tot(src)` gives

    ρ_c^msg(dst)  =  ρ_c(src)·r  =  φ_c(src) · ρ_tot(dst)

— the source's density *share* applied to the destination's observed total. There is no level transport
in it: it is exact iff `φ_c(src) = φ_c(dst)` and wrong by exactly `φ_c(src)/φ_c(dst)`. ⛔ When the source
carries only gDNA that factor is `1/φ_g(dst)` and the delivered level collapses to `ρ_tot(dst)`, the
destination's own total, independently of the source's measurement
(`TRAPS: a-message-from-the-destinations-belief`). So a message makes two claims on two scales: a
composition (a share, scaled by `r`, licensed by §3.5b) and a level (an absolute rate, carried unscaled
as a profile over `log ρ`, because gDNA is uniform along the genome before capture). Under capture a
boundary's level is a lower bound on the exon inside it — a fragment spanning a gene-end boundary lies
only partly under the probe — which is why every gDNA-lane hop is lower-only.

**3.5b The licence — "is the source measuring the same thing I am?"** (ruling 2026-08-04;
``outside_flank`` and ``boundary_shares_strand`` in `native/transfer_kernel.h` are the predicates). A composition may
be imputed across a step iff both hold:

* **SUPPLY** — the source supplied both components of the pair (a statement about precision): a source
  carrying one component has no share to lend; its λ is undefined, not "large".
* **POPULATION** — the two objects measure the same RNA population: otherwise enrichment and a
  population difference are indistinguishable, and the imputation reads one as the other.

The population conjunct is set algebra, not a model. A boundary counts what spans it contiguously, so

    T(BOUNDARY)  =  T(NODE_left) ∩ T(NODE_right)

and `T(BOUNDARY) = T(right)` fails iff a transcript's body begins at the boundary, `T(BOUNDARY) = T(left)`
iff one ends there. A transcript terminus is exactly what makes one flank's population larger; the test
is an equality per `(BOUNDARY, side)` pair, per step rather than per object. ⛔ Write it in genomic
terms, never in TSS/TES, which the strand flips:

    the RIGHT flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic LOW  end is there  ⟺  TSS₊ or TES₋
    the LEFT  flank gains RNA at a BOUNDARY  ⟺  a transcript's genomic HIGH end is there  ⟺  TES₊ or TSS₋

⛔ Termini only. A donor/acceptor boundary also changes the population, but there the flux is measured
(`junction_count`) and the SPLICE IN and SPLICE OUT maps route it. A terminus has no flux to measure: a
transcript simply begins. That is the derived line between the two treatments.

**3.5e The two operators and the terminus, in `{gDNA, RNA+, RNA−}`.** The ruling is `DESIGN.md` §0c.0
(2026-08-19); `transfer_kernel.h`'s `splice_faces` and `terminus_rules` implement it. A message is the three
densities `{gDNA, RNA+, RNA−}`, always: an operator that pools the two RNA components has not measured
what the policy carries.

* **SPLICE OUT — exon region → `exon|intron` boundary.** The boundary knows the `+` sj leaves here. It
  rescales the exon's densities by the ratio of totals *with the sj flux counted in its own total*
  (§3.6c), subtracts the sj density at the exon's scale, and undoes the rescale. Worked: exon
  `{10, 90, 0}`, boundary `{2, 1, 0}` with `SJ+ = 17` → rescale ×5 → `{10, 90, 0}`, `SJ+ = 85` → subtract
  → `{10, 5, 0}` → unscale → `{2, 1, 0}`, the unspliced message compared with the boundary's own. The
  rescale absorbs any common-mode opportunity or enrichment; nothing branches on capture.
* **SPLICE IN — `intron|exon` boundary → exon region.** The inverse: the boundary sends its densities
  *including* the spliced-in flux, which continues contiguously into the exon; the exon rescales by
  `ρ_tot(dst)/ρ_tot(src)` and compares.
* **THE TERMINUS — where that arithmetic is not licensed.** Two `+` transcripts starting at 1000 and
  2000 and sharing an exon to 5000: an unspliced fragment crossing 2000 is compatible with the first and
  with gDNA, not with the second, so the population in (2000, 5000) is not the population at the
  boundary. The rescale cannot account for a transcript that originates after the boundary; only the
  level can be propagated inward (THE LEVEL RULE, `terminus_rules`). A terminus at a splice junction is
  still a terminus: `sj+term` is ruled with `term`, never with `sj`.

**3.5g A total abundance must not be `mass / effective_length` — the accumulator already deposits the
composition-free quantity** (`calibration/total_abundance.py`; ruling 2026-08-20). An effective length is
a function of the fragment-length distribution, and gDNA and RNA have different ones, so `mass / E` is a
function of the composition being solved for and any enrichment ratio built from it is circular (100
counts in a 500 bp region read 0.25 as pure gDNA and 0.33 as pure RNA at means 100 and 200). The
reciprocal-opportunity deposit (§2) cancels the opportunity on its own support:

    E[ Σ 1/(w−1) ]      =  rho · P(w ≥ 2)   =  rho   at a BOUNDARY (frag_min ≥ 2, any pmf, any composition)
    E[ Σ 1/(ℓ−w+1) ]    =  rho · P(w ≤ ℓ)            in a REGION of length ℓ — NOT rho

So the boundary form is a model-free total and the region form only a density shape (§2.2), and every
REGION↔BOUNDARY ratio still carries a pmf functional — the chain alternates the two, so that is every
hop. A face's total adds the sj flux whose bodies lie on that side, in the same units, from the sj bank's
own reciprocal-opportunity column (`TRAPS: a-face-total-is-not-a-total-without-its-flux`). The
per-component divisors `E_g`, `E_r` stay length models; what this buys is that a boundary's total, and
any enrichment ratio between boundaries, uses none.

**3.5h The premise variance — why an imputation must cost something on every hop.** The ruling is
`DESIGN.md` §0c.0c; ``hop_price`` in `native/transfer_rows.h` and the pair terms in the builders'
`terminus_rules` and `alternative_splice_site` implement it. Every variance that scales with counts
vanishes between two deeply-counted slots, so a layer built only from counting terms delivers an
imputation at full strength beside a measurement. The premise of a hop — "my neighbour's values apply
here" — is not a counting statement and does not shrink with depth; it is estimable from the pair itself
by method of moments, which keeps it constant-free:

    Var_obs( log r )  =  premise  +  counting        ⇒     premise  =  max( 0,  (log r)² − counting )

with `r` the ratio of the two nodes' densities and `counting = 1/n_s + 1/n_x`. ⛔ Floored at 0 rather
than at something small: a pair whose ratio varies no more than Poisson predicts has exhibited no
heterogeneity, and a fit may not manufacture doubt the data do not show. Charged per hop on the pair's
own witness and never pooled across pairs, it makes a deep imputation arrive weaker than a shallow one
and a measurement outweigh both.

**3.6 The two faces of an `intron|exon` boundary — component-set matching.** Ruling 2026-08-04;
`transfer_kernel.h`'s `splice_faces`. At one boundary the accumulator stores three populations, and their
component sets differ:

| bank | what it counted | components |
|---|---|---|
| `unspliced_count` = `U` | crossed the boundary contiguously, spliced nowhere | gDNA + unspliced RNA |
| `junction_count` = `J` | never crossed it — it jumped from here | spliced RNA, certified |
| `spliced_count` = `S` | crossed contiguously, spliced elsewhere | spliced RNA, certified |

Spliced RNA cannot cross an exon↔intron boundary contiguously
(`TRAPS: mature-rna-never-crosses-a-boundary`), so it is absent from `U` and present in `J`. Three
densities, one of which needs no deconvolution:

    rho_g   = C_g / E_g       unknown split of U          U = C_g + C_u
    rho_u   = C_u / E_r       unknown split of U          (RNA unspliced at this position)
    rho_j   = J / E_J         measured — certified RNA

Matching component sets against each flank, `T(INTRON) = T(BOUNDARY, U only)` and
`T(EXON) = T(BOUNDARY, U + J)`, so with the same `rho_g` in both numerators:

    (I)  INTRON face:  phi_g(BOUNDARY) = rho_g / (rho_g + rho_u)            ==  phi_g(INTRON)
    (II) EXON   face:  phi_g(BOUNDARY) = rho_g / (rho_g + rho_u + rho_j)    ==  phi_g(EXON)

One boundary, one gDNA density, two composition statements differing only in whether the sj term is in
the total. (I) plus the boundary's own mass identity `U = rho_g·E_g + rho_u·E_r` is two equations in two
unknowns, and the intron's composition is prior-free (the intron factory); (II) then delivers the exon's
composition with every term measured.

The estimator is not the identity. The sj sees only `E_J/(E_J + Σ E_r)` of a transcript's spliced
fragments (about 11 % for two exons), so at low RNA `rho_j` is a handful of counts; face (I)'s job is to
transport a well-counted gDNA level from the intron to the boundary, and face (II) closes the composition
two ways — via `rho_j` (tight at high RNA, over-stated by §3.6b's frame gap) or via the exon's own mass
identity closed with `rho_g` (strong at low RNA) — which fuse by inverse variance, never by precedence.
⛔ The well-counted side is not a fixed side: under capture the intron is off-probe and nearly empty
while the boundary holds tens of counts, so the direction of transport inverts and only the share
survives it (`TRAPS: capture-inverts-the-counted-side`). The transfer policy carries the certified flux
as the splice-in map's cap rather than solving face (II) outright (`DESIGN.md` §6b.13).

**3.6c The splice-flux reframe — a boundary has two totals, one per flank.** Ruling 2026-08-05;
`ChainView.sj_count_lo` / `sj_count_hi` (`messages/__init__.py`) carry the split and
`transfer_kernel.h`'s `splice_faces` reads it. §3.6 made per step: which flank is a hop talking to? Numerator and
denominator of a composition imputation must be totals over the same component set, and a molecule
counted in `J` spliced at this position, so its body lies in the exon on exactly one side:

    at a sj's genomic-LOW  end   its exon is on the LOW  side
    at a sj's genomic-HIGH end   its exon is on the HIGH side

Hence one total per flank, split per sj:

    rho_lo  =  rho_U  +  Σ_{j : low end here}   J_j / E_J,j          — used against the LOW  neighbour
    rho_hi  =  rho_U  +  Σ_{j : high end here}  J_j / E_J,j          — used against the HIGH neighbour

and at a region both sums are empty (a contained fragment used no sj), so every junction-free chain is
unchanged. Direction does not enter the pairing: a hop joins adjacent slots `(k, k+1)`, so
`r = rho_lo[k+1] / rho_hi[k]` always, and the two totals are indexed by role (a boundary at a sj's low end
is the destination of the hop from its low flank, flux included, and the source of the next hop into its
high flank, excluded). One junction-inclusive total per slot inflates the intron-facing side by exactly
`J/E_J`, in opposite directions on the two hops of a pair, so the error cancels in a compounded ratio and
no aggregate check sees it (`TRAPS: recompute-from-the-oracle`).

⛔ Write the predicate in genomic terms (§3.5b). `splice_graph`'s `FLAG_DONOR_s` marks the genomic-low end
of an `s`-strand intron on both strands, so on `−` it sits at the transcript's biological acceptor; the
names are a misnomer on `−`, the data is uniform, and a predicate keyed on the sj's strand flips sign
silently. A rule keyed on the coarse region type cannot do this either: where regions are simultaneously
intron and exon, `coarse_type_array` reports every gene-body region as `exon`.

**3.6b The sj and contained frames are a lever on the fragment-length mean, at 0.62 %/bp.**
`calibration/effective_length.py` (`crossing_eff_length`, `contained_eff_length`) — the length model as a
divisor, live, and not the deferred composition channel. `E_J` and the exon's `E_r` are built from one
pmf and are exactly consistent (they sum to the exon length), but they differentiate with opposite sign:

    E_J   = E[w] − 1          →   dE_J/dE[w]  = +1
    E_r   = e − E[w] + 1      →   dE_r/dE[w]  = −1
    ⇒  d log(rho_j / rho_R) / dE[w]  =  1/E_J + 1/E_r  =  0.0062 per bp   (E[w] ≈ 200, e = 1,000)

so a length model wrong by `Δ` reports the sj estimator as `1 + 0.0062·Δ` times the exon's own
(`TRAPS: two-divisors-opposite-sign`). There is no second, geometric term: the simulator reweights the
length marginal by the same opportunity it places with (`wgs_engine._post_capture_length_allocation`), so
the realised crossing count at length `w` is `f_pre(w)·(w−1)` and the placement factor cancels. The
inflation is `rho_j`'s: `rho_R(exon) ≥ rho_u(B) + rho_j(B)` is a correct lower bound, diluted by
`1 + (1−s)(k−1)` with `s` the unspliced share of the exon's RNA.

---

## 3b. The conserved mass, and why one share for two components is a bias

`accumulator.cpp`'s mass bank (`boundary_unspliced_mass`), gated by `tests/native/test_conserved_mass.py`;
`tests/calibration/test_prior_units.py` gates the locus-scale consequence. The accumulator deposits `+1`
on every boundary a fragment crosses, so a sum over objects is an object-incidence count while the EM
adds a fragment count. The mass bank closes that gap: it sums to one per fragment.

**The deposit.** A fragment of length `w` is cut by the crossed boundaries into slices; a slice of length
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
It is a censored functional, so it is sensitive to the pmf's shape, not only its mean.

**The pooling theorem.** Define a component's share at a boundary as the mean conserved mass one of its
crossings carries, `share_c = E_c[A_mass] / E_c[w−1] ∈ (0,1]`. The accumulator can only measure the
mixture's, `share_pooled = M/C = φ·share_g + (1−φ)·share_r`. Rescaling both components by it gives

    a_g + a_r  =  M_g + M_r                    ← the TOTAL is conserved, exactly
    (â_g/â_r) / (a_g/a_r)  =  share_r / share_g   ← and is INDEPENDENT of the true mixing ratio

So the error is purely compositional and no conservation check can see it
(`TRAPS: conservation-misses-mis-attribution`). It is identically zero iff `share_g = share_r`, i.e.
whenever the two components' length pmfs agree, which is why an equal-length panel is structurally blind
to it (`TRAPS: an-equal-length-panel-defeats-the-lift`) — a price the ladder pays on purpose (§3.1); the
fl-gap side panel (`TESTING.md`) exercises it. The cancellation is exact per boundary, not at a locus:
contained mass never passes through the share (a contained fragment deposits on exactly one region,
already a fragment count) and summing over boundaries with different flanks re-introduces a weak
dependence. So `share_r/share_g` is the right mechanism and the wrong magnitude for a locus-level
correction; a repair belongs per boundary, before the contained term is added.

---

## 4. Opportunity corrections for length pools

`calibration/sj_opportunity.py` (§4.2) and `calibration/gdna_opportunity.py` (§4.3). A pool is a
length-dependent selection, so its raw histogram is not the library's length distribution. Let `A(w)` be
the pool's opportunity and `T(w)` the total opportunity for the same population.

**4.1 The divisor is a probability.**

    pi(w) = A(w) / T(w)          fitted(w) = count(w) / pi(w)

⛔ Never `count(w)/A(w)` — `TRAPS: divide-by-a-probability`.

**4.2 The sj pool** (`sj_opportunity.sj_opportunity`, `crossing_probability`, `detilt_pool`). For a
transcript with exon lengths `e_1..e_K` and total `L = Σ e_i`, the starts at which a length-`w` window
crosses at least one sj:

    A_j(w)  =  (L − w + 1)₊  −  Σ_i (e_i − w + 1)₊

Derived via the complement — a window crosses no sj iff it lies wholly inside one exon, and the exons
are disjoint — so there is no inclusion-exclusion. The library quantities are abundance-weighted sums,
`T(w) = Σ_t θ_t (L_t − w + 1)₊` and `A(w) = Σ_t θ_t A_j(w,t)`. `θ` is a molar abundance (copies), not an
observed fragment count — `A_j` already counts start positions, so a count applies the length weighting
twice. Production uses a uniform `θ` over the non-synthetic transcripts: the ratio cancels most of the
dependence.

**4.3 The four gDNA pools** (`gdna_opportunity.contained_opportunity`, `crossing_opportunity`,
`total_opportunity`). Two contained (§1.3), two crossing-exactly-one-line (§1.4), with
`T(w) = Σ_refs (L_ref − w + 1)₊`. The combination:

    f(w)  ∝  [ Σ_p count_p(w) ] · T(w) / [ Σ_p A_p(w) ]

This is the opportunity-weighted average of the four de-tilted pools, `Σ_p A_p f_p / Σ_p A_p`, and under
Poisson counts `Var(count_p) ∝ A_p`, so those weights are exactly inverse-variance. There is no tunable
weight. ⛔ It is not the same as pooling the four histograms and applying one divisor
(`TRAPS: opposite-tilts-must-not-pool`).

**4.4 What an opportunity correction cannot do.** Both forms above assume the population is placed
uniformly over its template. Under hybrid capture gDNA is not: placement is proportional to a capture
landscape the tool cannot see. That residual is a placement model, not a better divisor.

---

## 5. Strand

**5.1 The likelihood — the only intrinsic gDNA/RNA signal** (the kernel's `strand_term`, `native/transfer_rows.h`;
its two-component form is the gates' readable reference, `tests/calibration/_psi_reference.strand_loglik`).

    p    =  ½·f_g + κ·(1−f_g)
    var  =  N·p(1−p)  +  (N·f_g)²·¼·od_g  +  (N(1−f_g))²·κ(1−κ)·od_r
    loglik = −½·(sense − N·p)²/var − ½·log(var)

**5.2 What it can and cannot say.** With RNA tilt `d = f₊ − f₋`, `p = ½ + (κ−½)·d` — the gDNA fraction
cancels identically. Strand measures the tilt; it reaches gDNA only through the triangle bound
`f_g ≤ 1 − |d|`.

    I(f_g)  =  N·(2κ−1)² · [f_g(1−f_g)]² / (4·p(1−p))

exactly zero at κ = ½ for any count and any overdispersion, and saturating in `N` at
`(½−κ)²/(p(1−p)·od)` — capped by dispersion, not by depth.

**5.2b The strand channel is live iff the protocol preserves strand — a decision, not a band**
(`region_init.strand_discriminability`, landed 2026-09-14). `κ` is fitted, so `(2κ̂−1)²` is a squared
point estimate and strictly positive on a genuinely unstranded library. The code therefore evaluates

    N_eff  =  N / (1 + (N−1)·od_r)                          the OVERDISPERSED effective count
    disc   =  4·(κ̂−½)²  if BF₁₀ > 1  else  0                 REPLACES  (2κ−1)² = 4(κ−½)²
    I_strand  =  N_eff · disc · [f_g(1−f_g)]² / (4·p(1−p))

`p` here is §5.1's `p = ½·f_g + κ·(1−f_g)` (in the code `κ + f_g·(½−κ)`), not §5.2's tilt form. `BF₁₀`
is the Bayes factor, on the spliced 2×2 the strand fit read (`n_same` sense reads of `N`,
`strand_balance.fit_strand_balance`), of a free `κ` under the fit's own Beta(1, 1) prior (H1: a
strand-preserving protocol) against `κ = ½` EXACTLY (H0: an unstranded protocol — read 1's strand is
independent of the transcript's, so every sj reads ½ and the pooled split is Binomial(N, ½) with no free
parameter). Both marginals are closed form:

    ln BF₁₀  =  N·ln 2 + ln B(a, b),      a = n_same + 1 = κ̂·(N+2),   b = n_opp + 1 = (1−κ̂)·(N+2)
             ≈  ½·[z² − ln(2N/π)],          z = (κ̂ − ½) / √(¼/N)                          (large N)

so the channel is live iff `BF₁₀ > 1` at equal prior odds — no constant. The asymptotic form is the
free parameter's Occam penalty: the excursion a sampling fluctuation must clear to be read as a
protocol grows as `√ln N` (2.0σ at N = 100, 3.0σ at 10⁴, 3.8σ at 3·10⁶), while a stranded library
clears it by orders of magnitude (the ladder: −2.0 … −7.3 nats on the eight unstranded rows,
+3.9·10⁴ … +2.5·10⁶ on the eight stranded ones). With no spliced observation the two marginals are
equal and the channel is dead; `calibrate` raises before that on a real library.

What this replaced (2026-08 → 2026-09-14): `disc = 4·max(0, (κ̂−½)² − σ²_d)` with
`σ²_d = ¼(1/N_rna + od_r) + ¼(1/N_gdna + od_g)`, an unbiased estimate of `(κ−½)²` floored at zero.
Under H0 `(κ̂−½)²/Var(κ̂)` is χ²₁, so the floored estimate is positive with probability 0.32: a coin
toss, not a deadband — the ladder's `g98 ss.50 OFF` (z = 1.22) shipped with a live channel, and
`g00 ss.50 OFF` (z = 1.05) with a dead one only because of the gDNA term, without which it read 21,484
false gDNA fragments. That term had no derivation: gDNA's strand mean is ½ by symmetry (§5.3) and is
not estimated, `od_g` is per-slot noise that lives in §5.1's variance and not in `Var(κ̂)`, and its
`1/N_gdna` switched every gDNA-free library's channel off at any κ (`N_gdna = 0`, the modal real case
and all four ladder `g00` rows). The `od_r` in `Var(κ̂)` was wrong in the same way: the pooled mean
averages the per-sj spread over the number of sj, not once.

**5.3 gDNA's strand term is ½** — double-stranded, no sense direction. A fitted mixture marginal was
implemented and refuted: any constant for the unstranded case cancels out of the orientation
discrimination `(1−p)/p`, so a global genic average destroys most of the signal.

---

## 6. Overdispersion and second-moment evidence

`gdna_strand.overdispersion_for_beta` and `_null_information`. Symmetric `Beta(a,a)` gives
`od = 1/(2a+1)` (a=2 → 0.200, a=14 → 0.0345). Effective count `n_eff = n/[1 + (n−1)·od]` — at od 0.2 a
1,523-fragment seed is worth five coin flips. Pooled moments:

    od = Σ_s[(k_s − n_s μ_s)² − n_s μ_s(1−μ_s)] / Σ_s n_s(n_s−1) μ_s(1−μ_s)

Its exact null information `I = (Σ n(n−1)pq)² / Σ[2n²p²q² + npq − 6np²q²]` collapses to the pair count
`Σn(n−1)/2` only at mean ½. Second-moment evidence is counted in pairs of fragments inside one object — a
singleton carries exactly zero.

### 6a. The away-half moment — gDNA overdispersion with no pure seed

`gdna_strand.away_half_moment` and `fit_gdna_strand_overdispersion` (2026-08-29). The gDNA fit needs
seeds whose RNA does not read as strand spread, and no structural class can be asserted pure (pervasive
transcription; the intergenic space is whatever the GTF leaves over). Orient each genic seed so that RNA
of its own gene pulls the residual down, `d = (k − n/2)·sign(½ − κ)`; under pure gDNA `d` is symmetric
about 0 and the moment excess `d² − n/4` is even in `d`, so the pooled moment restricted to the away half
has the null expectation of the full one:

    od = Σ_s a_s·(d_s² − n_s/4) / Σ_s a_s·n_s(n_s−1)/4        a_s = 1[d_s > 0] + ½·1[d_s = 0]

Unbiased for `ρ_g` under any distribution of RNA content across the seeds; a contaminated seed reaches the
away side only by noise, with small `d`, so it biases down, never up. At `κ = ½` exactly — reachable,
since `κ = (n_same+1)/(n_obs+2)` is exactly ½ whenever `2·n_same = n_obs`, the modal outcome on an
unstranded library — the orientation is degenerate and the full two-sided moment is used instead: RNA at
the same mean ½ is symmetric, so it contributes only `n_g(n_g−1)/[N(N−1)] ≤ 1` of the excess and the
one-sided guarantee survives. Without that branch every residual collapses to 0 and the fit returns a
hard `od = 0`, the most confident strand likelihood assertable. The tie weight ½ is exact: at `n = 2`
under `BetaBinom(2, ½, ρ)`, `P(k=1) = ½(1−ρ)` and `P(k=0) = P(k=2) = (1+ρ)/4`, and the away half returns
exactly `ρ` while full tie weight returns `(3ρ−1)/(3−ρ)`. Half the pairs enter, so the information is
half §6's pair count — `I = P/2` for the total `P`, since `Var(e_s)|₀ = n(n−1)/8` and `E[a_s] = ½` give
`Var(num) = P/8`, `E[den] = P/4`, `Var(od_mom) = 2/P`. ⛔ Not half the away half's own pair count: that
halves twice and overstates the standard error by √2. ⛔ Requires a gene strand to orient by — intergenic
and AMBIG objects cannot enter — and unannotated antisense RNA pushes toward the away side, the one
recorded way to inflate it.

### 6b. Influence weighting — why a deep seed is not worth its pair count

`gdna_strand.between_seed_variance` and `influence_weights`; the root is found by bisection in
`fit_gdna_strand_overdispersion` (2026-08-30). Pooling `od̂_s = (d_s² − n_s/4)/(n_s(n_s−1)/4)` by pair
count is minimum-variance only at `ρ = 0`. Given the seed's latent rate `p` (`u = p − ½`),
`E[od̂_s | p] = 4u²` exactly (since `E[d² − n/4 | p] = n(n−1)u²`), so by the law of total variance

    Var(od̂_s | ρ) = V∞(ρ) + E_p[Var(od̂_s | p)] ≈ 2ρ²(1−ρ)/(1+2ρ) + 2/(n(n−1))

with the between-seed term exact from the symmetric Beta's moments (`E[u²] = ρ/4`,
`E[u⁴] = 3ρ²/(16(1+2ρ))`; Monte-Carlo-verified). `V∞` does not depend on `n` — a seed's information about
ρ saturates with depth — so the inverse-variance (Gauss–Markov) weights are

    w_s = 1/(½ + c_s·V∞(ρ))          c_s = n(n−1)/4

a constant at ρ = 0 (the pair-count estimator is this one with ρ pinned at 0) and equal-per-seed once
`c_s·V∞ ≫ ½`. No constant is introduced, and the one approximation — the sampling term at `p = ½` rather
than integrated over `p` — cannot bias the fit, because the weights depend only on `n_s` and ρ and never
on a seed's own data, so the ratio has expectation ρ for any weight function. ρ enters its own weights,
and the root of `g(ρ) = clip(moment(ρ)) − ρ` is bracketed by construction (`g(0) ≥ 0`,
`g(ceiling) ≤ 0`), so bisection terminates with no iteration limit.

The mean matters, and the two components do not share it. In general

    V∞(ρ, μ) = 3ρ²·[2ρ + μ(1−μ)(1 − 7ρ)] / [μ(1−μ)(1+ρ)(1+2ρ)] − ρ²
    b_s      = c_s·Var(od̂_s | ρ=0) = (2·n·pq + 1 − 6·pq)/(n − 1)      w_s = 1/(b_s + c_s·V∞(ρ, μ))

which reduces algebraically to `2ρ²(1−ρ)/(1+2ρ)` and `b_s = ½` at μ = ½. ⛔ At a real library's
κ = 0.0023 the same ρ = 0.05 gives V∞ = 0.285 against gDNA's 0.0043 — a seed at an extreme mean carries
far less information about ρ than its pair count suggests, so the two components may never be compared
on pair counts.

### 6c. The two components reconcile against each other

`gdna_strand.reconcile_overdispersions` (2026-08-30). A weighted estimator's precision is
`Σ 1/V_s = Σ w_s·c_s` — its own weighted denominator — so each component reports the precision of the
estimate it actually made, at its own ρ and its own μ. The weaker then borrows its deficit from the
better-measured one:

    od_w' = (I_w·od_w + (I_s − I_w)·od_s) / I_s          borrow weight (I_s − I_w)/I_s

0 when the two are equally informed (neither moves), 1 when the weak one measured nothing (it takes the
other's value outright). ⛔ Not `(I_w·od_w + I_s·od_s)/(I_w + I_s)`, which drags an equally well-measured
component to the midpoint while the other stays put; and not a symmetric pooling, which would erase a
real difference. With neither measured both take the ceiling: any common value leaves the strand channel
uninformative, since the composition term reads only the difference of the two dispersions.

---

## 7. Background rate, and deconvolving counts without strand

**7.1 A faint background rate is measurable only in aggregate** (`gdna_density.pooled_log_rate`,
`one_sided_rate`). `ρ_bg = Σg/ΣE`, `Var(log ρ_bg) ≈ 1/Σg`. A region of effective length `E` resolves a
rate only above ~`1/E` (Fisher information = `ρ·E`), so no per-region estimator finds it and true zero is
resolved sharpest. One-sided: `ρ_bg > 0` proves DNA present; `ρ_bg ≈ 0` does not prove absence, because
capture depletes the off-target floor. Never a denominator or a scale.

**7.2 Counts against a gDNA background with no strand data** (`density_deconv.fit_gdna_background`,
`_log_negbinom`). `P(g|C) ∝ P_bg(g)·1[0 ≤ g ≤ C]` with `g ~ NegBinom(ρ_bg·E_g, α_eff)`, a flat one-sided
prior on the RNA excess, truncated at the observed total. `α = Σμ²/max(Σ(g−μ)² − Σμ, 0⁺)` (∞ ⇒ Poisson),
`1/α_eff = 1/α + 1/(Σg + n₀)`. No tuned constant.

---

## 9. Priors on a grid

`simplex_logodds._JEFFREYS_REF`. `p(ρ_c) ∝ ρ_c^(c−1)` gives `Beta(c,c)`. `c = ½` (Jeffreys) is the only
grid-width-stable choice — `TRAPS: no-prior-means-haldane` for what omitting the term does instead.

### 9a. Why a simplex vertex is unreachable without evidence — and why that is not headroom

ψ lands short of `f_g = 1` on unexpressed genes, and that
gap is a theorem, not a bug:

* every proper prior with a density on `[0,1]` has a median strictly inside `(0,1)`;
* an object with zero composition evidence has posterior = prior;
* ⇒ a vertex is unreachable there in any coordinate, at any depth.

The estimator is honest rather than merely stuck: per object, `|f_g − truth| / sd(f_g)` sits inside its
own 1σ on both vertices. That companion is scoped to the stranded case — `sd(f_g)` shrinks like
`n^(−1/2)` only where the strand channel buys it, and §5.2 says the channel is identically zero at
`κ = ½`, so on unstranded data the shortfall is depth-independent. ⛔ Do not quote an `n^(−1/2)` shrinkage
on an unstranded stratum, or read a flat shortfall there as a regression. The ceiling a vertex-pinning
arm measures is therefore the value of missing information, not headroom, and "fit a prior to fix it" is
circular: pass-0 must stay prior-free to produce the substrate a prior is fitted on. A certified-RNA
count of zero is consistent with `f_g = 1` too (`tests/calibration/test_certified_rna_licence.py`).

### 9a.1 The proof needs a continuous CDF, and an atom escapes it

The first bullet runs through the CDF: `F` continuous and strictly increasing on `(0,1)`, so `F(x) = ½`
has a solution there — properties of a density, not of properness. For a general proper prior the median
is `inf { x : F(x) ≥ ½ }`, and an atom at `x₀` carrying mass `π ≥ ½` makes that infimum exactly `x₀`. So
a Beta (§9c's family, atom-free) cannot put its median at a vertex at any strength in any coordinate,
while a spike-and-slab can whenever its spike carries at least half the mass. "The value of missing
information" is thus a statement about the prior family the tree ships, not about the data. ⛔ Not a
licence to widen a bound, soften the reference, or fit a prior at pass-0: an atom would have to be earned
by a measurement that places it, and none ships.

---

## 9b. The EM's RNA prior goes to every RNA component, in proportion to its evidence

`calibration/priors.assemble_priors` hands the EM two per-locus pseudocounts;
`native/em_solver.cpp:apply_grouped_prior_update` applies them. The gDNA one lands additively on the
single gDNA component; the RNA one is shared among the RNA components in proportion to the evidence
each already carries:

    rna_count = Σ_{i ≠ g} raw[i]                     the whole RNA pool, and the whole set of recipients

    out[g] = raw[g] + gdna_prior
    out[i] = raw[i] · (1 + rna_prior/rna_count)      every i ≠ g

Summed over the RNA components this is exactly `rna_count + rna_prior`, so the gDNA:RNA split is
unchanged and the rule redistributes strictly within the RNA pool. ⛔ That is a per-M-step identity
for given `raw` counts, not an end-to-end invariant: the EM iterates, a different `theta` gives a
different E-step, hence different `raw`, hence a different converged split.

**Why every RNA component, with none singled out** (owner, 2026-09-19). RNA is RNA (Axiom 0). Whether
the annotation happens to assert a given RNA component — a synthetic nascent entity is a shadow span
the index manufactured — is a fact about the annotation, not about this locus's composition, and the
pseudocount is a statement about composition. One consequence makes the rule sayable in a line and is
the whole reason to want it: because the weights echo the EM's own current belief, the prior enters
every RNA component as the SAME factor `(1 + rna_prior/rna_count)`, so it moves the gDNA:RNA split —
which is what it is for — and nothing else.

**What the withheld share did, measured.** The rule this replaced held synthetic entities out of the
denominator and left them unscaled, on the null that a manufactured span is absent until the data
proves otherwise. That made the factor un-common over the pool, so the prior ALONE redistributed RNA
between entities the data cannot tell apart. On a locus of two components explaining all 200 fragments
equally well, one of them synthetic, a prior of 500 drove the synthetic component from 100.0 fragments
to 2.79e-298 and handed all 200 to the other — the withheld factor compounds once per M-step
(`tests/test_estimator.py`). On the antisense-intronic scenarios the displaced mass landed on the one
annotated transcript that could also explain it: the leak onto an unexpressed antisense `t2` fell from
70 fragments to 0 at SS 0.65, and the two nested-antisense rungs from 124 → 14 and 24 → 2
(`tests/scenarios/test_antisense_intronic.py`; `ISSUES: nested-antisense-leak-under-the-sane-ruler`,
closed by this rule).

**A zombie decays against the TRANSCRIPT it shadows, and GROWS against gDNA.** The withheld factor was
also an anti-zombie force, and dropping it is the restoration's whole price: the geometric decay rate of
a shadow entity holding nothing of its own falls from `kappa/(1 + rna_prior/annotated_count)` per
iteration — the withheld rule's own denominator, over the components it did admit — to `kappa = w_N/w_T`.
That is still strictly below 1 for free, since a shadow span is longer than the transcript it shadows, so
the decay against the MATURE isoform survives — slower, and now a property of the likelihood rather than
of the prior. ⛔⛔ **THAT GUARANTEE DOES NOT COVER THE OTHER COMPETITOR.** The survival criterion is unchanged and
still derived: with `m` fragments whose only RNA candidate is the entity, it grows iff
`m·w_N > Total·theta_g·w_g`. Written in densities that is `m/L_n > Total·theta_g/L_g` — the entity's
footprint against the gDNA component's average — and in a locus whose gDNA is uniform the footprint holds
its share of it (`m ∝ L_n/L_g` of the total), so the factor is the density ratio and not `L_g/L_n`: a
whole MultiLocus's one gDNA opportunity does not by itself destabilise a shadow. What makes `theta_n = 0`
unstable under capture is the gDNA component priced at `ρ q̄` of its density (§11), and it is measured
with one E-step from the TRUE counts (2026-09-20): off capture the truth is a fixed point (the shadows
return −1.6 %); on capture it drifts +24 % per step, +47 % without the priors, and +9 % with the gDNA
opportunity at its oracle value. The toy below is a locus with ALL of its gDNA under one shadow — the
contest the note of 2026-09-19 solved — and is kept as the record of what the solver does there, not of
the ladder: off capture `L_g/L_fam` reaches 5 with no leak, on capture families over-claim 3.5× where
the ratio is 1, and single-twin families leak 4.3× against 3.7× for families of ten
(`ISSUES: nascent-siphons-gdna-under-capture`). On a pool of `N` fragments that only gDNA and the
shadow can explain, split evenly by genome strand, the fixed point in `r = b/a` is

    r · L_n/L_g = [ ss·r/(ss·r + ½) + (1−ss)·r/((1−ss)·r + ½) ]
                / [    ½/(ss·r + ½) +        ½/((1−ss)·r + ½) ]

whose only root at `L_n ≥ L_g` is `r = 0`. Solved and confirmed against the shipped solver
(`tests/test_estimator.py`), at `ss = 0.99` the shadow's share of `N` runs 0 % at `L_g/L_n = 1`, 35.2 %
at 2, 52.7 % at 6.2 and 82.0 % at 20 — 44.4 % at 6.2 with a `gdna_prior` of `N/2`. ⛔ A length knob is not
a repair here either — it kills the live entities with the dead ones. Under the NASCENT SCOPE RULING
(`DESIGN.md` §0b) nascent RNA is modelled for robustness, which is an argument for treating it like
any other RNA here and not for a null that suppresses it.

### 9b.1 The prior is an additive per-component pseudocount — the weights are the design

`em_solver.cpp`. The rule above is habitually described as "multiplicative, hence neutral on the
within-RNA split". That is a description of one choice of weights, not of a different kind of update.
Writing `R` for `rna_count` and `P` for `rna_prior`:

    out[i] = raw[i]·(1 + P/R)  ==  raw[i] + P·raw[i]/R  ==  raw[i] + a_i,    Σ a_i = P

so the prior is `a_i = P · w_i / Σ w` at `w_i = raw[i]` — an allocation in proportion to the EM's own
current belief. A prior that echoes the posterior carries no information, which is exactly why it is
neutral on the split it does not exist to set. `AggregatePrior::component_rna_prior_weight` generalises
`w`; the lane is built end to end and `pipeline.py` passes none
(`ISSUES: per-transcript-prior-lane`), so `w_i = raw[i]` is what ships.

**The one place the weights are not interchangeable is `raw[i] = 0`, and it is why an EQUAL SHARE was
refused** (owner, 2026-09-19). With `w_i = raw[i]`, `out[i] = 0` is an absorbing state — no prior
magnitude revives a component with no warm-start evidence, since `alpha` floors to `EM_LOG_EPSILON`
and `digamma` of that is `−1e300`. A weight vector that lifts every component off zero removes it, and
the threshold to clear is the VBEM fixed point `alpha = Σ_u resp(alpha)`, at which a component
actually activates (~0.16–0.47 alpha units on the shipped EM), not the exponential cutoff (0.0014).
`P` is a conserved FRAGMENT COUNT — calibration's unspliced RNA mass on the locus
(`priors.assemble_priors`), tens to thousands on an expressed locus — so a flat `P/n` would clear that
threshold by one to two orders of magnitude at every component, reviving any shadow entity outright
and flattening the within-RNA split toward uniform, an assertion nothing measured licenses. Admitting
every RNA component at `w_i = raw[i]` restores fairness and keeps the absorbing state, because the
state is a property of the WEIGHTS and not of the eligibility test that was removed.

### 9b.2 The prior cancels exactly under MAP; under VBEM the digamma residual is what is left

`em_solver.cpp`, gated in `tests/test_estimator.py`. §9b's "it moves the gDNA:RNA split and nothing
else" is exact in the MAP M-step, where `theta_i ∝ alpha_i`: the common factor `c = 1 + P/R` divides
out of every within-RNA ratio and the converged answer is bit-for-bit what it was at `P = 0`. Under
VBEM — the shipped mode — `theta_i ∝ exp(psi(alpha_i))`, which is not scale-equivariant, so the
cancellation is asymptotic rather than exact. With `psi(x) = log x − 1/(2x) + O(x^-2)`,

    psi(c·alpha_i) − psi(alpha_i) = log c + (1 − 1/c)/(2·alpha_i) + O(alpha_i^-2)

and `log c` is common, so it normalises away. What survives is a per-component residual
`(1 − 1/c)/(2·alpha_i)`, largest at the smallest alpha and bounded by `1/(2·alpha_min)` since `c ≥ 1`.
So the allocation's uniformity is intact and only the M-step's own nonlinearity moves a share, by an
amount that vanishes as the locus deepens. ⭐ It is still two orders of magnitude from the thing a
gate must separate: on a three-component locus the residual moves a share by 1.4e-3 against its bound
of 1.5e-2, while the eligibility rule §9b replaced moved a component by 100 % of its mass. ⛔ Do not
read a small VBEM drift here as a defect in the allocation, and do not widen the MAP gate to
accommodate it — they are different statements about different M-steps.

## 9c. ψ's composition reference is a Beta, and its mean would be a third term

`calibration/simplex_logodds.py` (`_gdna_arm`, `_rna_arm`, `_JEFFREYS_REF`). `a·log f_g + b·log(1−f_g)`
on the λ grid is exactly Beta(a, b) in `f_g`, so `a` and `b` are pseudo-counts with a strength `a + b`
and a mean `a/(a+b)`. The shipped `_JEFFREYS_REF = ½` fixes both: strength 1, the ignorance statement;
and mean ½, which asserts the library is half gDNA.

**The derivation.** Give the two components independent Gamma rate priors `ρ_c ~ Gamma(α_c, β_c)` with
Poisson counts `n_c ~ Poisson(ρ_c E_c)`. Then `X = ρ_g E_g ~ Gamma(a, s_g)` and `Y = ρ_r E_r ~ Gamma(b, s_r)`
with `s_c = β_c/E_c`. Writing `f = X/(X+Y)`, `T = X+Y` and integrating `T` out:

    p(f, T) ∝ f^(a−1)(1−f)^(b−1) · T^(a+b−1) · e^(−T(s_g f + s_r(1−f)))
    ⇒ p(f)  ∝ f^(a−1)(1−f)^(b−1) / ( s_g·f + s_r·(1−f) )^(a+b)

and on the solver's grid, where `|df/dλ| = f(1−f)`,

    log p(λ) = a·log f_g + b·log(1−f_g) − (a+b)·log( f_g + r·(1−f_g) ) ,   r = s_r/s_g

The shipped reference is exactly this with `s_g = s_r` — matched scales, the "half gDNA" assertion
appearing as a missing term rather than as a chosen number. With `m` an object's prior expected
composition `ρ̄_g E_g / (ρ̄_g E_g + ρ̄_r E_r)` and `β_c = α_c/ρ̄_c`, the ratio is `r = (b/a)·m/(1−m)`, and at
`a = b = ½` the third term collapses to

    − log[ (1−m)·f_g + m·(1−f_g) ]

Four properties, all structural: (i) `m = ½` makes the bracket the constant ½, so the term drops — a
strict generalisation agreeing with the shipped constant exactly where its assumption is true; (ii) the
tails stay `e^(−|λ|/2)` for every `m` — in the `+λ` tail the bracket → `m(1−f_g)` and the `log(1−f_g)`
cancels one Jeffreys half — so only the location moves and `L`-invariance is untouched, which moving
`a`/`b` instead cannot do (they set the tails and the location together, and `b = 0.03` leaves 57 % of
the mass outside `L = 10`); (iii) proper for every `m ∈ (0,1)`; (iv) substituting
`u = f/(f + r(1−f)) ~ Beta(a,b)`, symmetric at `a = b`, gives `median(f_g) = m` in closed form.

⛔ **ψ ships without this term.** A reference location built on it was refuted and deleted on
2026-08-24 (`DESIGN.md` §6b.1): a location is a prior assertion at fixed strength, and where the strand
channel is dead it was the entire answer at any depth. Background information enters as likelihood terms
whose precision scales with counts (the intron factory's rows, built in the solve's kernel from `density_deconv`'s background; the landscape prior). The
derivation stays because it is what any future location would be judged against.

### 9c.1 The strength of a reference mean is a log-odds, and one pseudo-observation sets it

Record of the refuted location (`DESIGN.md` §6b.1), kept because `DESIGN.md` and `ISSUES.md` cite its
arithmetic. The term `−log[(1−m)f_g + m(1−f_g)]` is `−log(1−m)` at `f_g → 1` and `−log m` at `f_g → 0`, so
its full range — how far it can move ψ — is

    strength  =  log( m / (1−m) )  =  logit(m)          ⇒   the claim's ODDS are e^strength

The location written on the λ scale is its strength in nats. Bayes sets it on the reference's own
exponents: one pseudo-observation of gDNA takes `Beta(a, b)` to `Beta(a+1, b)`, whose mean is

    m  =  (a + 1) / (a + b + 1)  =  0.75   at   a = b = ½        (strength log 3 = 1.0986 nats)

— not `m = σ(a+b)`, which equates a pseudo-count with a nat. The lattice may cap a claim but may not
choose it: `m = σ(L)` saturates at `strength = L − log 2 + O(e^(−L))`, i.e. 9.31 nats at `L = 10`, and
above `L ≈ 20.72` sticks where `1 − σ(L)` falls under the location clamp
(`TRAPS: a-clamp-at-the-closed-end-escapes-the-window`). The overturn depth in fragments is
`n = strength / log(2κ)` — 1.46 fragments for `log 3` at `κ = 0.99` — and it diverges at `κ = ½`, where
posterior = prior and a location of any strength is the whole answer: the mechanism of the refutation.

### 9c.2 Why a reference's curvature may not enter `τ_λ` — the Jacobian

Ruling 2026-08-16, `region_init.py` (ψ's reference contributes nothing to `tau_lam`). `I_strand ∝
[f_g(1−f_g)]² / (4p(1−p))` (§5.2), so moving a slot's mode from `f_g = 0.98576` to `0.99975` lowers the
strand term by `[0.98576·0.01424]² / [0.99975·0.00025]² = 3,154×` — the likelihood is genuinely that flat
on λ near the vertex, and nothing is lost there. `tau_lam` is the data's Fisher information; a prior's
curvature is not the data's information (`TRAPS: a-priors-curvature-is-not-the-datas-information`), and
feeding one in acts as a boolean gate that releases the whole count precision, at empty slots too.

## 9d. sweep_replay: the tolerance budget — why the read-out cannot amplify a rounding error

`scripts/profiling/sweep_replay.py --tolerance` reports, per output array, the slots moved, the largest
absolute and relative move, and a budget beside them. The budget is derived, not chosen, and it has two
halves.

**The read-out is perfectly conditioned.** ψ's answer at a slot is a posterior mean over the grid,
`f = Σ_k w_k σ_k` with `w_k ∝ exp(ψ_k − max ψ)`. Perturb cell k's log-density by `δ_k`; to first order
`w_k` moves by the factor `(1 + δ_k)` and

    δf = Σ_k w_k (σ_k − f) δ_k,   so   |δf| ≤ max_k |δ_k| · Σ_k w_k |σ_k − f| ≤ ½ · max_k |δ_k|,

because σ lies in [0, 1] and a variable's mean absolute deviation about its mean is at most half its
range. K does not appear: a max-normalised log-profile through exp and a normalised sum amplifies
nothing, at any grid size. The log-variance read-out `Var(log f) = Σ w_k g_k² − (Σ w_k g_k)²` has the same
form with `(g_k − ḡ)² − Var` in place of `σ_k − f`, bounded by the squared range of `g = log σ(λ)` on the
window, `L̃² = log²(1 + e^L)`, so `|δ Var| ≤ L̃² · max_k |δ_k|`. The fraction is actually read as the
posterior MEDIAN (§9a): a weight error of relative size `δ` shifts the CDF by at most `δ`, so the median
moves by at most `δ` over the posterior density at the median — the mean's bound wherever the posterior
is unimodal and its density there is not small. A balanced bimodal posterior has no such bound: its
median jumps between the modes under any perturbation, and a slot the report flags beyond the budget is
its way of pointing at one.

**What an implementation perturbs.** Any implementation that keeps the solver's expressions and differs
in libm and summation order rounds each of the `T` terms of a cell within one rounding unit `ε` of that
term's magnitude, so `|δ_k| ≤ ε · Σ_terms |term_k| ≤ ε · T · A` with `A` the largest term magnitude a
cell can carry. Two scales bound it. The strand Gaussian `−½ (u − n p)² / var` scales as
`n / (2κ(1−κ)) = c_κ · n`, since `|u − n p| ≤ n` and its variance is at least `n κ(1−κ)` on the window;
with strand overdispersion `od > 0` the variance grows as `n²·od`, which caps the term at `c_κ / od`
however deep the slot, so on deep slots the count bound is loose by that factor — a bound, still. The
fitted arms are kernel log-densities whose bandwidth is the grid step, so on a `K`-cell grid they are
bounded by the kernel's range, `K²/2`. The rows are profiles of the same likelihoods, blurred, and carry
neither scale beyond them. With `N` the largest slot count on the chain,

    A = max(c_κ · N, K²/2),      B_f = ½ · ε · T · A,      B_var = L̃² · ε · T · A,

with `T = 6` (the strand term, two arms, the factory rows, the message rows, the cube row) and `ε` the
solve's rounding unit — float64's `2⁻⁵³`, the whole of ψ's precision. The
bound is loose by construction (a full cancellation at a high-weight cell is assumed), which is why the
report always shows the actual move beside it. Achievable rounding lands orders inside it: float32
rounding of every term and of the strand mean at disagreeing slots moves a fraction by 10⁻⁴–10⁻⁵ of the
budget (the gate records the scale), because neighbouring cells' rounding errors are incoherent and the
read-out averages them; the chunk-exact reordering of 2026-09-11 moved ≤ 3.1e-15 per slot — the
read-out's own rounding, `K · ε`, with no term re-rounded. Gate: `tests/test_sweep_replay_tolerance.py`
— the budget covers float32 rounding of the terms and of the strand mean at every strength.

## 9e. ψ's θ quadrature — the nodes follow the strand term's peak (`native/psi_kernel.h`'s `tilt_window`)

**The integrand.** At an AMBIG slot ψ is read out on its θ-marginal, `M(λ) = ∫ exp ψ(λ, sin θ) dθ` over
`θ ∈ [−π/2, π/2]` (§9c: the arcsine measure on the tilt `τ` is cancelled by the θ coordinate, so no measure
weight is written). Only the strand term and a delivered cube row depend on θ. The strand term
(`_mixture_strand_loglik`) is a Gaussian quasi-likelihood in the aligned-strand rate `p` with the variance
frozen at the reference composition, and `p = ½·f_g + κ·f₊ + (1−κ)·f₋ = ½ + a(λ)·τ` with
`a(λ) = (1 − f_g)(κ − ½)`, so at fixed λ it is an **exact Gaussian in τ**:

    L(λ, τ) = −(τ − τ̂(λ))² / 2σ_τ(λ)² + const(λ),    τ̂ = d/a,   σ_τ = σ_p/|a|,   d = u₊/n − ½,   σ_p = √V/n.

Its width in θ is `σ_θ = σ_τ / cos θ̂ ≈ 1/(2√n·|κ−½|·(1−f_g)·cos θ̂)`: 0.015 rad at 5k fragments, 0.005 at
50k, 0.0015 at 500k.

**Why a fixed lattice fails, and how.** With `K_t` uniform nodes at step `h = π/(K_t−1)` (0.053 at 60), wherever
`σ_θ < h` the lattice sum at one λ is `≈ exp(−δ(λ)²/2σ_τ²)` for `δ` the distance from `τ̂(λ)` to its nearest node
— a factor anywhere between 1 and `e^{−(h/2)²/2σ_τ²}` (`e^{−100}` at 50k) chosen by where the drifting centre
`τ̂(λ) = d/a(λ)` happens to fall. Across λ that is a COMB: the λ posterior of a deep interior-tilt slot is a set
of spikes at arbitrary λ, and the read-out is a coin toss on the node placement (the ladder's recorded K_t 30
failure was one 25k-fragment slot with a 0.006 rad peak; 60 nodes happened to land on it). At a strand-pure
slot the centre is `τ̂ = 1/(1−f_g)`, on the boundary at `f_g = 0` and beyond it after; there the peak is a quartic
in `δ = π/2 − θ` of width `∝ n^{−¼}` at the boundary and a one-sided quadratic of width `σ_τ/√(τ̂−1)` beyond —
the lattice's error there is λ-dependent too, but smaller, and was not the recorded mechanism.

**The rule.** Per `(slot, λ)`, integrate only where the term has mass. With `τ_m = clip(τ̂, −1, 1)` the term's
maximum ON the domain,

    ρ = √((τ_m − τ̂)² + 2σ_τ²T),    [τ_lo, τ_hi] = [τ̂ − ρ, τ̂ + ρ] ∩ [−1, 1],

is the window where it lies within `T` nats of that maximum — one closed form for an interior peak, a peak on
the boundary and a peak beyond it. Place `K_t` uniform nodes in θ across `[arcsin τ_lo, arcsin τ_hi]` and sum
them with the trapezoid weights: `M(λ) ≈ h·[½g₀ + g₁ + … + ½g_{K_t−1}]`, written into ψ as `log h(slot, λ)`
plus `log ½` at the two end nodes. In θ the integrand `g(θ) = exp ψ(λ, sin θ)` is smooth and EVEN about `±π/2`
(`sin(π − θ) = sin θ`), so it is the restriction of a smooth periodic function and the trapezoid rule is
spectrally accurate in every regime: at an interior window end `g` is `e^{−T}` of its peak and the end
weight is immaterial; at a domain end the `½` is exactly the periodic trapezoid rule's weight on the
reflected window, on the quartic and the one-sided quadratic peak alike. The `log h` term is not optional:
the window scales with `σ_τ(λ) ∝ 1/(1−f_g)`, and that λ-dependence is the marginal's own — the fixed lattice
sampled it wrongly, a rule that dropped `log h` would drop it altogether. A slot with no strand information
(`κ = ½`, `n = 0`, or `ρ ≥ 1` on either side) gets the whole domain, i.e. the uniform lattice with trapezoid
weights — the rule degrades to the lattice exactly where the lattice was right.

**The two constants are derived.** `T = −log ε₆₄ ≈ 36 nats`: the term's mass outside the window is
`erfc(√T) ≈ e^{−T}/√(πT) < ε₆₄` of the peak's integral — the truncation is below double precision, so no
wider window can change a bit. `K_t`: the interior window is `2√(2T)·σ_θ` wide, and the trapezoid rule's error
on a Gaussian of width `σ_θ` at spacing `h` is `2·e^{−2π²(σ_θ/h)²}`, below `e^{−T}` once `h ≤ σ_θ·π√2/√T`, so
`K_t − 1 ≥ 2√(2T)·√T/(π√2) = 2T/π ≈ 23` — **`K_t = 24`** (`_TILT_NODES`). At a domain end the window is only
`T^{¼} ≈ 2.4` widths per side, which 24 nodes over-resolve. Measured against adaptive quadrature over
`n ∈ [500, 500k]`, `f_g ∈ {0, 0.3}`, `τ ∈ {0, 0.5, 0.9, 0.99, 1}`: the λ-shape error of `log M` is ≤ 2·10⁻⁶ nats
at 24 nodes (12 nodes reach 0.07 and are refused); the fixed lattice at 60 reads 0.01–0.5 nats at n = 500,
8–12 at 50k and 90–130 at 500k, and 2,400 lattice nodes still 0.1–0.4 at 500k. The gate is
`test_vertex_reference.test_the_theta_marginal_matches_adaptive_quadrature_at_every_depth`.

**A delivered row is evaluated at the nodes.** The RNA level lanes deliver a row's ingredients
(a row of `simplex_logodds.CubeRows`: the held profile per strand over `u = log(ρ/ρ_ref)`, the slot's total `n` and RNA
opportunity `a_r`, each lane's `ρ_ref`) and ψ evaluates them at its own nodes: at each cell the strand's share
`f_s = (1 − f_g)(1 ± τ)/2` implies the density `f_s·n/a_r`, and the held profile is read at `log(ρ_s/ρ_ref)`
(the kernel's row map, `profile_of_level` with the tilt inside; the readable form is the gates' `_psi_reference.row_at`). No θ lattice exists for a row to be built
on and nothing is interpolated. Before this, rows built on a 60-node lattice equalled rows on 240 on the
shared-exon stress and rows on 24 did not — the lattice's own resolution error, now gone with it.

**What the exact marginal makes visible (an open issue, not the quadrature's).** At an interior tilt the exact
marginal is `∝ σ_τ(λ) = σ_p/((1−f_g)|κ−½|)` up to where the band fills the domain (`1 − f_g ≈ σ_p/|κ−½| ≈
1/√n`): the volume of tilts consistent with the data grows as `f_g → 1`, an Occam factor of order `√n` toward
gDNA at a balanced slot — Bayes-correct under an `f_g`-independent arcsine measure on τ, and contrary to §5.2's
"strand reaches gDNA only through the triangle bound". On a real chain the RNA level lanes and the landscape
prior pin such slots; a slot with no junction and no prior reads `f_g ≈ 1 − 1/√n` from its own solve. See
`ISSUES: strand-marginal-volume-factor`.

## 9f. The tilt atom — the AMBIG tilt's hypothesis space is {pure +, pure −, mixed} (`native/psi_kernel.h`'s `slot_cube`)

**What §9e's exact marginal cannot say.** At an AMBIG slot whose RNA is all on one strand the truth sits AT
the strand cap: with `τ = 1` the split `p̂ = ½ + (1 − f_g)(κ − ½)` identifies `f_g` exactly as at a
single-strand slot. But ψ integrated the tilt over a continuum, and every `f_g` below the cap fits the same
`p̂` with a slightly impure tilt `τ̂(λ) = d/a(λ) < 1`; the θ-marginal is spread over `[0, cap]`, weighted by
the strand term's width `σ_τ(λ) ∝ 1/(1 − f_g)` (§9e's volume factor), and its median lands below the cap.
Prior-free, one slot, at the pure tilt: a truth of 0.50 read 0.31–0.37, 0.85 read 0.69–0.80, 0.97 read
0.85–0.95 (n = 30 … 30k). On the ladder's stranded × capture-ON rows this was the largest AMBIG-class error
in scope (−26k net on `g50 ss.99 ON`), because under capture the prior that repairs it off capture is broad.

**The hypothesis space.** Axiom 0's opportunity geometry admits both RNA strands at an AMBIG slot, but PRESENCE
per strand is discrete: the slot's RNA is all on `+`, all on `−`, or on both. The tilt's reference measure is
therefore a mixture of three hypotheses at equal reference weight,

    P(τ)  =  ⅓·δ(τ − 1)  +  ⅓·δ(τ + 1)  +  ⅓·(arcsine measure on (−1, 1)),

i.e. in `θ = arcsin τ`: two atoms at `±π/2` and the uniform density `dθ/π` between them (§9c's arcsine
measure, which the θ coordinate carries with no weight). The per-λ marginal is

    M(λ)  =  ⅓·[ e^{L(λ, +1)}  +  e^{L(λ, −1)}  +  (1/π)·∫ e^{L(λ, sin θ)} dθ ]

and the ⅓ cancels from the λ posterior. In ψ this is two more θ columns per slot, `τ = ±1` exactly with
log-weight 0, beside the continuum's `K_t = 24` windowed nodes whose trapezoid log-weights carry `−log π`
(the window's share of the domain: with no strand information the window is the whole domain, its weights
sum to π, and the three hypotheses' masses are equal at every λ — the gate
`test_vertex_reference.test_the_three_tilt_hypotheses_carry_equal_reference_weight`). The pure columns are
single-strand solves inside the cube — `f₊ = 1 − f_g, f₋ = 0` and the mirror — the strand term evaluated at
the vertex, no tilt parameter and no width. Where the data are pure the pure hypothesis explains them with
no parameter and wins the Occam contest against the continuum's `∝ σ_θ·e^{L_max}`; where they are not, the
pure column sits `e^{−(1 − τ̂)²/2σ_τ²}` below the peak and vanishes. The read-out is unchanged: `f_g` the
posterior median over the θ-marginal, `w₊` the RNA-mass-weighted share over every column, atoms included.
Every θ-independent term (the arms, the λ-factor rows) is common to the three hypotheses and cancels.

**The witness.** An atom's cost is at an INTERIOR tilt: at `f_g = cap` the pure hypothesis also explains
`p̂` with no parameter, so a both-strand slot is pulled toward the cap (prior-free, a truth of 0.30 at
`τ = 0.5` reads 0.46 → 0.63 with the plain atom). The RNA level lanes carry each strand's PRESENCE as a
delivered lower bound (§6b.13), and a held level on strand `s` is a certified witness that `s` carries RNA,
so it rules the hypothesis "all the RNA is on the other strand" out: in ψ the pure `−s` column is `−∞`
wherever the slot's delivered row holds a profile on `s`. Nothing is pooled and no constant enters: the witness
is the delivery itself (a level on `s` says nothing against "pure `s`"; with nothing delivered both atoms
stand). The structural witness — the per-strand exon bits, a strand whose RNA here could only be nascent
cannot be the pure carrier — adds nothing measurable on top (the θ note §13: within 2 % by presence truth)
and is not written. Where a both-strand slot holds a level on only one strand the residual cost remains; it
is a limit of the information and is accepted as such (`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`).

## 10. The second pass's score

`src/rigel/second_pass.py` (`combine_factors`, `choose_hypotheses`). `f(L)` here is the second pass's
per-fragment length term (its `length_likelihood` array); the deferral of the length composition channel
does not touch it — this one ranks one fragment's already-enumerated candidates against each other and
never claims a composition.

    score  =  ρ × f(L) × s

normalised within one fragment's candidate set, factors applied in order of the evidence behind them and
skipped when flat-zero among the survivors (`TRAPS: an-all-zero-factor-is-inert`). One multinomial draw
per fragment, in the side buffer's canonical order, from a single RNG stream. ⛔ Never key the draw on the
fragment's content: a content hash ties on exactly the duplicates it would harm, so 100 identical
fragments would draw identically and a 60/40 posterior would collapse to 100/0.

Known approximation: `ρ` enters as a hard multiplicative zero, but zero observations is
`P(0 | λ, E) = e^(−λE)`, not zero. The hard zero is the large-exposure limit of the correct likelihood, so
it is right where the library is deep and wrong where it is shallow.

## 11. The ruler — a transcript's bases at their pieces' capture efficiencies (`capture_eff_length`, `capture_efficiency`, `priors.assemble_priors`; the reference, `abundance_landscape.located_enriched_mode`)

Under hybrid capture the EM divides a transcript by its effective length under capture: its
FL-marginal length times the mean capture efficiency of its bases,

    eff_em_t = fl_t · factor_t,    factor_t = Σ_p ℓ_p^τ · c̃_p / Σ_p ℓ_p^τ,

over the pieces `p` (regions of the partition) the transcript's exons overlap, with `c̃_p` the piece's
capture efficiency — its gDNA density against the fully captured level `ρ_ref`, clipped at 1 — and
`ℓ_p^τ` the taper-weighted count of the transcript's bases in the piece. `factor_t ∈ (0, 1]`, and it is
exactly 1 when no piece is depleted relative to the reference. The locus gDNA component takes the same
sum over the locus's pieces, introns included, since gDNA's template is the genome.

**The length is a sum over bases (role one).** A transcript's effective length under capture is the
sum over its start positions of the efficiency of the fragment starting there, and a fragment's
efficiency is the mean per-base efficiency over the bases it covers, so

    eff_t = Σ_x c̃(x) · τ(x),      τ(x) = E_f[ (starts whose fragment covers base x) / w ],

over the transcript's own bases `x` in transcript coordinates. A fragment of length `w` starting at
`s` covers `x` iff `s ≤ x < s + w` with `0 ≤ s ≤ L − w`, which is `min(x + 1, w, L − x, L − w + 1)⁺`
starts, each spreading its unit over `w` bases; so `Σ_x τ(x) = Σ_w f(w)(L − w + 1)⁺ = fl_t` exactly —
the per-base frame partitions the start count, it does not re-derive it — and away from both ends
`τ = 1`, only the bases within a fragment of an end being tapered (`effective_length.BaseTaper`; for
`L ≥ 2 w_max − 1` the four-way min is `min(d, w)` with `d = min(x + 1, L − x)` and one cumulative
table serves every template). With `c̃` constant on a piece the sum is the factor above. No boundary
object, no junction object and no contained support enters the length: a junction-spanning start is
counted through the bases it covers, a piece shorter than a fragment carries its bases' weight, and
`span = fl` holds for every structure by construction — the property the junction objects were built
to secure, which they secured by imputing forty junctions from pieces that had no support on the
ladder's dense annotation.

**The locus gDNA length counts what the count counts.** The EM's gDNA component for a locus carries
a COUNT, the calibration's gDNA mass on the locus's regions and boundaries (one crossing converted by
`q`, the conserved mass per crossing), and the length it divides that count by is those same objects'
starts at their own efficiencies, each start counted ONCE:

    eff_g = Σ_r S_r · c̃_r + Σ_e q_e · S_e · c̃_e,

the contained support of every region at its efficiency and the crossing support of every boundary at
the boundary's own efficiency `c̃_e = E[min(ρ_e/ρ_ref, 1) | k_e, S_e]`, converted by that boundary's `q_e`.
The conversion is the deposit rule's own: a crossing fragment deposits a count of +1 at EVERY boundary it
crosses and a mass summing to 1 across them, so a boundary's crossing support `E_f[w − 1]` (§2; gDNA's
reach is unbounded) counts INCIDENCES, and a start whose fragment spans a piece shorter than itself sits
in the support of both of that piece's boundaries. Under a uniform field of ρ fragments per start every
object reads its own density (`k_e = ρ S_e`, `m_r = ρ S_r`), the count is `ρ (Σ S_r + Σ q_e S_e) = ρ N_starts`,
and the UNconverted length `Σ S_r + Σ S_e` is `N_incidences`, so a component read against it has density
`ρ · N_starts / N_incidences = ρ q̄` — the field's only where no fragment crosses two boundaries. Enumerated
(`tests/calibration/test_priors.py`: three 40-bp pieces inside a long reference, four boundaries counting the
locus's outer two, fragments of 60), 179 distinct starts overlap the locus against 236 incidences,
`q = 49.5/59` at the outer boundaries and `40/59` at the inner ones, and `Σ q_e S_e = 179` exactly. At
`q = 1` — flanks longer than every fragment — the two forms coincide, so the factor-one identity holds as
before and a capture-OFF prior moves only at loci whose fragments span short pieces.

**What the unconverted support cost (shipped 2026-09-16 → 2026-09-20).** It had been refused by a
measurement on the test chromosome's capture-OFF transcript number, before nascent RNA's share of the RNA
prior was restored — 1–2 % worse, a contaminated toy's injection gate insensitive — which ranked a
calibration input on the thermometer and read the unmasking of a cancelling error as harm. Under capture
the introns contribute no opportunity and the probed exons are shorter than a fragment, so the crossing
support is 64 % of the locus length over the ladder's `g50 ss.99 ON` loci and 84 % in the loci that leak
(15 % off capture), and the boundary term reads `1/q̄` times the once-counted crossing gDNA it holds
(correlation 0.988 over 722 loci). A gDNA component priced at `ρ q̄` hands its sense fragments at the
probed exons to the RNA hypotheses, and the synthetic nascent entities, pinned by nothing, take them
(`ISSUES: nascent-siphons-gdna-under-capture`). Converting the support takes that row's siphon from
+541,216 to +32,908 fragments and its gDNA pool from −626,550 to −56,319, with `calibration_vs_oracle.py`,
`zero_controls.py` and `policy_benchmark.py` identical on every metric.

One other form stays refused by measurement: the transcript's per-base form over the locus's bases drops
the boundary objects whose masses the count keeps, and where the calibration's crossing masses sit above
their geometry the gDNA component then reads denser than its objects, over-claims the exonic unspliced
fragments and every probed gene under-calls — the test chromosome's `g50 ss.99 ON` row through the
thermometer, gene-level Σ|Δ| 25,633 → 38,174 against 23,967 with the object form, and the gDNA pool +5.5 %
against −0.1 %.

**The evidence is every unspliced gDNA object over the piece (role two).** gDNA is one template at a
uniform rate before capture, so its density after capture on a piece is that piece's efficiency up to
the one unit `ρ_ref`; and the accumulator counts a fragment at every boundary it crosses, so a piece
too short to contain a fragment is witnessed at its edges. Two kinds of count: the piece's own
CONTAINED count `k_p` on its support `S_p = E_f[(ℓ_p − w + 1)⁺]`, `E[k_p] = ρ_p S_p`; and the CROSSING
count `k_e` at every boundary within a fragment's reach, whose fragments' bases lie in the pieces they
cover in proportions the geometry fixes,

    E[k_e] = Σ_q A_eq · ρ_q,      A_eq = Σ_w f(w) Σ_{a=1}^{w−1} (bases of the placement in q) / w,

with `A_eq = G(c_j) − G(c_{j−1})` for the piece at cumulative distance `(c_{j−1}, c_j]` on a side,
`G(c) = Σ_w f(w) Σ_{a=1}^{w−1} min(a, c)/w`, and `Σ_q A_eq = E_f[w − 1]` exactly, the crossing
opportunity partitioned over the bases it counts (`effective_length.crossing_base_shares`, gated
against enumeration). A crossing at an interior boundary of a long exon has its bases in two pieces of
equal efficiency and reads that efficiency directly; a crossing at the edge of a 40 bp exon has a share
`≈ 35 / 216` of its base-starts in the exon and the rest in the intron, whose efficiency is pinned by
its own long contained count, so the exon's efficiency is identified from the edge counts against the
intron's known level.

**The efficiency is a posterior mean, not a plug-in.** Under the fitted gDNA landscape `P(log ρ)` — the
same `DensityLandscape` ψ's composition arm reads on the refits, the population's own statement of
where gDNA densities sit — and the Poisson counting rule as the likelihood,

    c̃_p = E[ min(ρ_p / ρ_ref, 1) | k_p, S_p, {z_ep, A_ep}_e ],

on the landscape's grid, with the crossings APPORTIONED: a crossing count is a Poisson sum over the
pieces within reach, and the E-step of a Poisson sum attributes it,
`z_eq = k_e · A_eq ρ̄_q / Σ_q' A_eq' ρ̄_q'` with `ρ̄` the pieces' own-count posterior means, each piece
then reading its share as a count on its own exposure `A_eq` (the terms of one piece pool exactly into
one Poisson term). The neighbours are held at what their own counts say — one pass. At high depth the
posterior is the plug-in `min(k/S/ρ_ref, 1)`; at low depth it is the population's mixture weighted by
the piece's own likelihood; a piece with no evidence reads the population's clipped mean; no constant
enters, and no floor: the multimapper floor `w = C/(C+1)` this replaces was a +3.4 to +3.7 nat bias on
every unprobed transcript at every gDNA level, protecting against a plug-in's exact zero that a
posterior mean never produces. `k_p` and `k_e` are deconvolved masses, so the Poisson enters through
the gamma function as a continuation, and `k_p` already carries the landscape through the solve, so the
prior enters twice in a small way — measured against the solve's own posterior (the belief's `f_g` and
`Var(log f_g)` as a log-normal, closed form), which loses badly where a slot is not solved
(+1.25 to +4.96 nat on the unprobed class), and against the clip taken outside the expectation,
`min(E[ρ]/ρ_ref, 1)`, which is indistinguishable on every row measured. Iterating the apportionment
with the updated means is EM on the joint and converges (13 passes on the test chromosome, 59 on the
ladder, at the grid step), changing nothing the truth instrument can see — identical on every test
chromosome row, +0.44 against +0.46 nat on the ladder's `g05 ss.99 ON` unprobed class — so the one pass
ships; a joint update of neighbours instead (each conditioning on the other's mean with the whole count)
never settles, 64–66 pieces of the test chromosome and ~2,000 of the ladder flipping by 8 nat every pass.

**`ρ_ref` is a population quantity.** It is the enriched mode of the population density `P(log ρ_g)`,
and the tool fits exactly that density: `landscape.DensityLandscape`, ψ's composition arm on the refits,
trained on the located compositions and the zero-count anchors (DESIGN §7.1). Its census
(`abundance_landscape._census`) partitions the grid into basins at the minima between interior maxima;
the depleted basin is the largest by rendered mass (for gDNA the unprobed objects outnumber the probed
ones — 0.70–1.00 of the mass on every row of both panels), and the enriched candidate is the basin
above it holding the most located kernels, `None` when nothing lies above.

**The location floor on the mode.** A basin's MEMBERS are the kernels with a location: a count of at
least one fragment, published by the fit as `DensityLandscape.located` beside each kernel's `centre`
(`log(max(count, 1)/E)`, in nats). A zero-count anchor or a sub-fragment kernel is centred at its
resolution wall `1/E`, which is where the kernel could not see and not where a density is, so it is no
member; a human index trains ~250–325 k anchors whose walls span every decade, and on a sparse library a
basin above the bulk can be packed with them around ten measured kernels
(`ISSUES: the-ruler-reference-on-sparse-real-libraries`). The candidate is a mode only if its members
resolve it at the located population's own resolution, `k = √n_located` (`landscape.knn_widths`' k):
each member's width is half the distance to its k-th nearest MEMBER, and

    n_members > k   and   median(width_k)² ≤ _LOCATED_VAR = 1 nat²,

the floor DESIGN §7.1 rule 4 applies to a slot, in the same variable — not a constant chosen but the
identity's value at the one-fragment wall (`Var(log c) = 1/c`), read at the population's own resolution. A
basin with k members or fewer has no k-th neighbour inside itself — the cluster smaller than √n that reaches
outside itself — and is no mode however narrow the rendered density's cut made it; the within-basin spread
is NOT the statement, because a basin cut by the grid's edge is narrow whatever its kernels (a 1-fragment
exon piece on 0.008 bp of support rendered a 0.30-nat basin at the top of the ladder's `g98 ss.50 OFF`
grid). Reading the members at their own √n_members instead would call twenty-one kernels strewn across
three nats a mode. The result publishes the count behind the reference
(`CalibrationResult.gdna_reference_members`). Measured: every enriched mode on both panels has 257–3,606
members at widths at the grid step, at a peak stable to 0.02 decades across an 8× range of the render
resolution; the two sparse real libraries hold 10–16 located kernels in any basin above the bulk against
k = 33 and 123, and read `None`.

**No enriched mode ⇒ no contraction, exactly.** `CalibrationResult.gdna_reference_density` is `None`,
every efficiency is exactly 1 on both axes (`gdna_capture_efficiency_region`, `_boundary`, which the
result refuses otherwise), the ruler returns `fl` verbatim and `assemble_priors` reads the locus's
uncontracted span. This is the capture-OFF field (unimodal, Poisson noise around
one level) and the gDNA-free field (the anchors' wall, with any false-positive basin above it a lone
kernel). The plug-in `min(ρ_n/ρ_ref, 1)` on a noisy uniform field is biased below 1 (Jensen plus the
clip), which is why a per-object reference read from the field itself contracted the oracle's own
counts by 8 % at capture-OFF; a modal decision has no per-object noise to clip.

**What the witness cannot see.** A probe is captured in genomic coordinates on gDNA and in transcript
coordinates on cDNA, and the two differ within a fragment of every junction a probe spans and at every
exon shorter than a fragment: on the ladder's transcript-designed panel gDNA at a split probe is captured
at a fifth of the cDNA's weight, and on the test chromosome's tiny-exon block a probe centred on a 40 bp
exon binds a gDNA fragment over 125 bp while the simulator's non-stacking rule binds a spliced fragment
over one exon's 40. The efficiency reads the gDNA and the transcript's factor inherits the difference
(`ISSUES: ruler-witness-geometry-on-transcript-panels`), declared and not repaired.

## 12. The flux price's witness — the column count on the protocol's share of the opportunity (`transfer_kernel.h`'s `rna_lane`)

The certified flux at one of an exon's junctions is that strand's RNA level at the exon (§6b.13's source):
the spliced count `c_J` at the junction's route rate `r_J = Σ flux / A_route`, and every hop pays the
pair's price (`hop_price`) — both witnesses' counting plus the disagreement of the two densities beyond
what counting explains, `max(0, log(r)² − (1/n_s + 1/n_x))`. The junction's side is in WHOLE-STRAND units:
`sj_count` and `route_rate` are keyed by the junction's own transcript strand, so `c_J` is every spliced
fragment of that strand's routes whatever column its reads landed on. The exon's side must be priced in
the same units. Its witness is the unspliced count on the column strand `s`'s RNA reads on, and under a
protocol whose read rate is `κ_read = max(κ, 1 − κ)` that column holds

    E[c_s] = κ_read·R_s + (1 − κ_read)·R_o + g/2,

so read on the exon's whole opportunity `a_r` the ratio to the junction's rate is `κ_read` when the pair
agrees at a single-strand gDNA-free exon, and the price carried `log(κ_read)²` of disagreement that is not
there — 0.48 nats² at κ = ½, on every flux level of every unstranded library. THE CORRECTION IS TO THE
OPPORTUNITY, NOT THE COUNT: a column is an opportunity of `κ_read·a_r` for the strand's RNA to be counted
on it, so the witness is `(c_s, κ_read·a_r)` — the column's own count, at its own counting precision, on
the protocol's share of the exon's opportunity. At κ = ½ its density is the exon's total, `2c_s/a_r`, a
bound on `R_s` (the two columns are exchangeable there); at κ → 1 it is the column's; at a both-stranded
exon on a strand-preserving protocol it is the strand's own share `R_s + g/2` and not the other strand's.
Nothing chosen: `κ` is the fitted strand model's.

Measured: the same chain at κ = 0.99, 0.7 and ½ delivers the flux level at counting alone whenever the
pair agrees in whole-strand units (the gate in `test_transfer_rna_lanes.py`); the golden
`strand_ss65_multi_iso`'s exon inside the other isoform's intron, with no gDNA in the library, reads
0.008 gDNA against 0.235 under the uncorrected witness.

Refused with numbers (2026-09-14): the exon's TOTAL unspliced count as the witness (whole units, two-sided)
is right on unstranded data but at a both-stranded exon on a strand-preserving protocol it mixes the other
strand in, widens the strand's level, and two AMBIG exons on the ladder's stranded capture-ON zero row read
0.60 and 0.18 gDNA in a gDNA-free library (233 → 435 false fragments on that row); the total as a one-sided
BOUND (only a junction rate above the exon's density contradicts a lower bound) loses stranded capture-ON
8–10 % on the test chromosome and 1.5 % on the ladder (the exon-holds-more widening protects the probed
exons); and the strand's own count from the column split, `R̂_s = (c_s − c_o)/(2κ − 1)` at its
Poisson-equivalent precision `(c_s − c_o)²/(c_s + c_o)`, loosens the golden's ceiling to 0.323 at κ = 0.65
(the asymmetry's precision is low and the price widens the level away) — the failure
`ISSUES: flux-witness-in-strand-units` recorded.
