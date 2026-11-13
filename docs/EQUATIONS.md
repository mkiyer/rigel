# EQUATIONS — the derivations the code depends on

Each section names the module and function that computes its result; where this text and the code disagree, the
code wins, and a gate in `tests/` holds each formula against brute force. Measurements, rulings, open problems,
lessons and history live elsewhere: a number here is a property of the formula. Section numbers are anchors cited
from tests and the other docs, so a deleted section leaves a gap.

---

## 1. The deposit rule

`native/calibration/accumulator.cpp` (the deposits), `calibration/splice_graph.py` (reach, the partition),
`calibration/effective_length.py` (the opportunities). The executable specification,
`tests/native/_accumulator_reference.py`, wins over this text.

**Notation.** A fragment is an interval `[s, s+w)` of molecule length `w`. A reference is cut into **regions** at every
exon endpoint of every non-synthetic transcript; the 0-bp cuts between adjacent regions are **boundaries**. `ℓ` is a
region length.

**1.1 Fragment length** `L` = genomic span minus the introns it spans. A mate gap counts toward `L`; an intron does
not. Whatever counts toward `L` must also count as coverage for crossing.

**1.2 Crossing count.** `[s, s+w)` crosses the inter-base point `p` iff `s ∈ [p−w+1, p−1]`: exactly `w − 1` starts. So
`E[count at p] = ρ · (mean_FL − 1)` for any length law, independent of both flanking regions.

**1.3 Contained opportunity** (`contained_eff_length`). A length-`w` fragment fits in a region at `(ℓ − w + 1)₊` starts:

    E_f[(ℓ − w + 1)₊] = (ℓ+1)·F(ℓ) − S(ℓ),     F the CDF,  S(ℓ) = Σ_{w≤ℓ} w·f(w)

and `ℓ + 1 − mean_FL` beyond the support. The `+1` counts start positions: without it a region one fragment long has
divisor 0.

**1.4 Crossing exactly one boundary** (`gdna_opportunity.crossing_opportunity`). With flanks `a` (left), `b` (right),

    A(w) = (w−1)₊ − (w−1−a)₊ − (w−1−b)₊ + (w−1−a−b)₊

exact: a fragment containing the boundary that reaches any boundary beyond `p−a` also crosses `p−a`. The partition
cuts at `0` and `L_ref`, so an outermost region's length is the distance to the wall.

**1.5 General crossing divisor** (`crossing_eff_length`). With `R_lo`, `R_hi` the molecule's own remaining template
either side, a crossing placement with `x` bases left needs `1 ≤ x ≤ w−1`, `x ≤ R_lo`, `w − x ≤ R_hi`:

    E_J = E_f[ min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)₊ ]

With the reaches ordered `a ≤ b`, the min is `w−1` to `w = a+1`, `a` to `b+1`, `a+b+1−w` to `a+b+1`, then 0: three
sums off §1.3's cumulative sums. At `UNBOUNDED_REACH` both sides it is `mean − 1` exactly.

**1.6 Reach.** At `p` inside exon `e`: `reach_lo = exonic_bases_before(e) + (p − e.start)`,
`reach_hi = total_exonic − reach_lo`, maximised over transcripts per side and per strand. gDNA's is unbounded; RNA's
stops at the transcript's ends. A reach of 0 is meaningful.

**1.7 Partition.** `region bounds = unique({exon.start, exon.end of non-synthetic transcripts} ∪ {0, ref_length})`,
no merging. Boundaries run `src < dst`, so genomic order is a topological order.

**1.8 The offer** (`AccumulatorSet::offer`). A fragment arrives with one PLACEMENT per resolved alignment, each with
hypotheses `H_p`. One placement deposits unchanged. Otherwise the strand check and the clip exclude a placement alone
(none left: `dropped_strand_undefined` if every exclusion was the strand, else `dropped_empty`); identical placements
are one; one left deposits with all its hypotheses; several keep the union of pairs with `L_{p,h} ≤ max_length` (all,
if none survives); one survivor deposits alone, two or more are held whole as `deferred_multiple_placements`. Ledger:
`deposited + deferred_undetermined_gap + deferred_multiple_placements + dropped_* == offered`, in fragments.

**1.9 The observable genome** (`observable_contained_eff_length`, `observable_crossing_eff_length`,
`observable_cut_shares`, `gdna_opportunity.observable_gdna_opportunity`). A start counts only where its fragment
touches no base of the index's unobservable mask (`rigel.observable`): a region's contained support is §1.3 summed
over its observable pieces, and gDNA's template at a boundary ends at the nearest masked base each side — its reach
in §1.5 and §11. With nothing masked each form is its plain one, to the bit.

---

## 2. Reciprocal opportunity, and where it is model-free

`accumulator.cpp`'s `1/(w−1)` banks (`Boundary::unspliced_inv_length_sum`, `SpliceJunction::inv_length_sum`); read by
the second pass's `ρ(h)` (§10). Deposit `1/A(w)`, `A(w)` the population's own opportunity at the object; by linearity

    E[Σ 1/A] = Σ_c ρ_c · P_c(A > 0)

A fragment with `A(w) = 0` deposits nothing, so `P(A > 0)` survives, a functional of the pmf the channel claims
independence from (`TRAPS: a-cancellation-is-conditional-on-its-support`).

**2.1 At a boundary it is model-free.** `A(w) = (w−1)₊`, so `E[Σ 1/(w−1)] = Σ_c ρ_c` whenever every fragment is at least
2 bp, across any mixture of length laws. The sj bank takes the same deposit.

**2.2 At a region it is a shape, not a total.** With `A(w) = (ℓ − w + 1)₊`, `E[Σ 1/A] = ρ · P(w ≤ ℓ)`, per component,
and 0 below `ℓ < frag_min`. No bank deposits it.

**2.4 Superadditivity.** `Σ E(children) < E(whole)` for contained effective lengths under any split, so densities,
effective lengths and variances do not pool across two partitions; only counts pool additively.

---

## 3. The two-component deconvolution

The per-slot solve and the transfer messages: `calibration/region_geometry.py`, `native/transfer_kernel.h`,
`native/transfer_rows.h`, `calibration/messages/transfer.py`.

**3.2 Density is the frame-invariant currency; a fraction is not.** `ρ_c = C_c/E_c` agrees across the contained,
crossing and spliced frames, while the log-odds shifts by `log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`. Capture
cancels from that ratio only when it is one common multiplier on both components at each object.

**Observed strand evidence** (`Chain::strand_profile`). At expectations `a = ρE_g` and `r` the strand columns have
means `a/2 + q r` and `a/2 + (1−q) r`; conditional on their total the `+` column is `Binomial(M, p)`:

    f = a/(a+r),     p = f/2 + (1−f) q,     log L(f) = u log p + v log(1−p) + const

`q = κ` on a `+`-strand slot, `1 − κ` on a `−`; constant at `q = ½`; no row without a strand model. It is the
likelihood of an expected fraction: an inferred fraction times `M` is no observed Poisson DNA count
(`ISSUES: the-reader-reads-false-gdna-as-capture`). Gate: `tests/calibration/test_observed_strand_claims.py`.

**Row arithmetic** (`transfer_rows.h`). A level row is a log-likelihood over `u = log(ρ/ρ_ref)` on a fixed lattice.

- *Coordinates.* A new reference `ρ_ref·e^d` translates every knot by `−d`; resampling at the old numbers would read
  other densities. The RNA lanes' `ρ_ref` (`TransferPolicy.library`) is the unspliced RNA density over single-strand
  exons; where 0 it is the ratio of summed positive count/exposure pairs — `(n, E_r)` per RNA-admitting slot,
  `(C, C/R)` per certified flux — a unit, not an estimate (`tests/calibration/test_message_coordinates.py`).
- *The blur* (`blur_row_at_step`): a discrete Gaussian of variance `v` at spacing `d`, radius `max(⌈4√v/d⌉, 1)`,
  `B_j = log Σ_t w_t exp(ℓ_{j+t})`; below the normal float range, the window-local log-sum-exp
  `m_j + log Σ_t exp(ℓ_{j+t} + log w_t − m_j)` (`m_j` the window's maximum), so tails keep relative likelihood.
- *A known source's footprint* (`flux_level`). A flux `C` at route rate `R` has exposure `A = C/R` and
  `ℓ(u) = C log(ρ_ref A e^u) − ρ_ref A e^u`; it is evaluated on `h` extra cells each side at the original spacing,
  blurred and cut back, so no output cell reads an edge-padded value (`tests/calibration/test_flux_footprint.py`).
- *The edge row* (`edge_level_row`). A pure-DNA boundary with count `N_b` and opportunity `E_b` gives
  `Q_b(ρ) = N_b log(min(ρE_b, N_b)/N_b) − min(ρE_b, N_b) + N_b` (Poisson below its mode, flat above), read at
  `ρ = σ(λ)N_e/E_e` on the recipient. It assumes the recipient's DNA density is at least the boundary's
  (`ISSUES: the-edge-density-floor-under-capture`) and is a projection through the observed total, not an intensity.

**3.2b Component opportunities in a face map** (`face_map_lambda`). A boundary with unspliced count `U`, count fraction
`f`, opportunities `E_g,b`, `E_r,b`, and certified RNA rate `s` on flank `x`:

    g = U f/E_g,b,     u = U(1−f)/E_r,b,     λ_x(f) = log(g E_g,x) − log((u + s) E_r,x)

The rate is the route sum `Σ_J J/A_J`, plus `S/E_r,b` for a contiguous spliced crossing. The forward face pushes the
boundary's profile through this monotone map; the reverse (`splice_out_row`) reads the exon's profile at `λ_x(f)`,
marginalising `s` over its counting width. Profiles, so no Jacobian. With no RNA added it is §3.2's shift. That the
two supplied laws share a capture frame is not established here. Gates: `test_splice_out_opportunities.py`,
`test_transfer_faces.py`; the route rate's width is `ISSUES: message-layer-open-cases`.

**3.3 Conservation with unequal effective lengths is `Σ_c ρ_c E_c = M`**, not `Σ_c ρ_c = M/E` (the message densities
`ρ = f·M/E`).

**3.4 The counting term** (`count_logvar`). `Var(log ρ_c) = 1/(f_c n) = Var(log f_c) + 1/n`, so the count must be an
integer count, not a fractional mass. Shipped: `trigamma(n + ½)`, the Jeffreys posterior's exact `Var(log ρ)` —
`≈ 1/n` from `n ≈ 10`, `π²/2` at `n = 0`.

**3.5 A composition crosses only where the population is shared; a gDNA level crosses unscaled** (`splice_faces`,
`terminus_rules`, `gdna_lane`). Rescaling a source by `r = ρ_tot(dst)/ρ_tot(src)` gives
`ρ_c^msg(dst) = φ_c(src) · ρ_tot(dst)`, exact iff `φ_c(src) = φ_c(dst)`: a gDNA-only source delivers the destination's
own total (`TRAPS: a-message-from-the-destinations-belief`). So a message makes two claims — a composition scaled by
`r` under §3.5b's licence, and a level carried unscaled over `log ρ`, since gDNA is uniform before capture. Under
capture a boundary's level is a lower bound on the exon inside it, so every gDNA-lane hop is lower-only.

**3.5b The licence** (`face_is_licensed`, `outside_flank`, `boundary_shares_strand`). A composition crosses a step iff
the source supplied both components and both objects measure the same RNA population. A boundary counts what spans it,
so `T(BOUNDARY) = T(left) ∩ T(right)`: `= T(right)` unless a transcript begins there, `= T(left)` unless one ends there.
In genomic terms, never TSS/TES:

    RIGHT flank gains RNA  ⟺  a transcript's genomic LOW  end is there  ⟺  TSS₊ or TES₋
    LEFT  flank gains RNA  ⟺  a transcript's genomic HIGH end is there  ⟺  TES₊ or TSS₋

Termini only: at a donor/acceptor the flux is measured (`sj_count`) and the splice maps route it.

**3.5e The operators, in `{gDNA, RNA+, RNA−}`.** SPLICE OUT (exon → `exon|intron` boundary): rescale the exon by the
ratio of totals with the sj flux in the boundary's total (§3.6c), subtract the sj density at the exon's scale,
unscale — exon `{10, 90, 0}`, boundary `{2, 1, 0}`, `SJ+ = 17`: ×5, `SJ+ = 85`, `{10, 5, 0}`, `{2, 1, 0}`. SPLICE IN
(`intron|exon` boundary → exon): the inverse, the boundary's densities including the spliced-in flux. THE TERMINUS:
where a transcript begins inside the step the rescale is unlicensed and only the level propagates (the level rule);
a terminus at a junction is still a terminus.

**3.5g A total is not `mass / effective_length`.** `E` depends on a length law the two components do not share, so
`mass/E` is a function of the composition being solved for. Face arithmetic runs on counts and each component's own
opportunity (§3.2b), and a face's total adds its flank's flux (§3.6c;
`TRAPS: a-face-total-is-not-a-total-without-its-flux`).

**3.5h The premise variance** (`hop_price`; the level rule's pair term; `alternative_splice_site` on two strand-mode
log-odds). "My neighbour's values apply here" does not shrink with depth, and it is estimable from the pair:

    Var_obs(log r) = premise + counting   ⇒   premise = max(0, (log r)² − (1/n_s + 1/n_x))

`r` the ratio of the two densities; the hop's variance is `count_logvar(n_s) + count_logvar(n_x) + premise`. Floored
at 0 (a pair within Poisson shows no heterogeneity), charged per hop on the pair's own witness, never pooled.

**3.6 The two faces of an `intron|exon` boundary** (`splice_faces`). The boundary holds `U` (crossed contiguously,
spliced nowhere: gDNA + unspliced RNA), `J` (jumped from here: certified RNA) and `S` (crossed contiguously, spliced
elsewhere: certified RNA). Spliced RNA cannot cross an exon↔intron boundary contiguously
(`TRAPS: mature-rna-never-crosses-a-boundary`), so with `U = C_g + C_u`:

    ρ_g = C_g/E_g,   ρ_u = C_u/E_r,   ρ_j = J/E_J
    (I)  intron face:  φ_g(BOUNDARY) = ρ_g / (ρ_g + ρ_u)          == φ_g(INTRON)
    (II) exon face:    φ_g(BOUNDARY) = ρ_g / (ρ_g + ρ_u + ρ_j)    == φ_g(EXON)

Given the intron's composition, (I) and `U = ρ_g E_g + ρ_u E_r` fix the boundary's split and (II) maps it to the
exon; the identities do not measure the intron. The transfer carries the flux as the splice-in map's cap. Under
capture the intron is off-probe and nearly empty while the boundary is counted, so transport inverts and only the
share survives.

**3.6c One total per flank** (`ChainView.sj_count_lo` / `_hi`, `route_rate_lo` / `_hi`). A molecule in `J` has its
body in the exon on one side — LOW at a sj's genomic-low end, HIGH at its high end:

    ρ_lo = ρ_U + Σ_{j : low end here} J_j/E_J,j         ρ_hi = ρ_U + Σ_{j : high end here} J_j/E_J,j

and a hop `(k, k+1)` reads `r = ρ_lo[k+1]/ρ_hi[k]`. One junction-inclusive total per slot errs by `J/E_J` in opposite
directions on a pair's two hops, which cancels in a compounded ratio. `FLAG_DONOR_s` marks the genomic-low end of an
`s`-strand intron on both strands (the biological acceptor on `−`), so a predicate keyed on sj strand flips silently.

---

## 3b. The conserved mass, and why one share for two components is a bias

`accumulator.cpp`'s `Boundary::unspliced_mass`; `CalibrationResult.boundary_mass_per_crossing`, read by
`priors.assemble_priors` (§11). Gates: `tests/native/test_conserved_mass.py`, `tests/calibration/test_prior_units.py`.

**The deposit.** The boundaries a fragment crosses cut it into slices; a slice of length `s` bounded by
`n_cross ∈ {1, 2}` crossed boundaries deposits `s/(w·n_cross)` at each, so a fragment deposits exactly 1 in total.

**The opportunity.** With flanks `a`, `b`, a crossing with `u ∈ [1, w−1]` bases left gives `g_a(u)/w` from its left
slice, `g_a(u) = u` for `u ≤ a`, `a/2` beyond; `Σ_{u=1}^{w−1} g_a(u) = w·min(w−1, a)/2` in both cases, so

    A_mass(w; a, b) = [min(w−1, a) + min(w−1, b)] / 2,     E[mass] = ρ_c · E_{f_c}[A_mass]

**The pooling theorem.** A component's share is `share_c = E_c[A_mass]/E_c[w−1]`; the accumulator measures only the
mixture's, `q = M/C = φ share_g + (1−φ) share_r`. Rescaling both crossing counts by `q` gives

    â_g + â_r = M_g + M_r                          the total, conserved exactly
    (â_g/â_r) / (a_g/a_r) = share_r / share_g      independent of the true mixing ratio

a purely compositional error (`TRAPS: conservation-misses-mis-attribution`), zero iff the two length laws agree. The
EM's prior reads the gDNA half, `â_g/a_g = q/share_g` per boundary (`ISSUES: the-pooled-q-in-the-gdna-count`).

---

## 4. Opportunity corrections for length pools

`calibration/sj_opportunity.py`, `calibration/gdna_opportunity.py`, `calibration/fl.py`. A pool is a length-dependent
selection; `A(w)` is its opportunity and `T(w)` the population's total.

**4.1 The divisor is a probability:** `π(w) = A(w)/T(w)`, `fitted(w) = count(w)/π(w)` (`detilt_pool`, which keeps the
pool's total as its evidence weight). Never `count(w)/A(w)` (`TRAPS: divide-by-a-probability`).

**4.2 The sj pool** (`sj_opportunity`, `crossing_probability_from_index`). A transcript with exons `e_1..e_K`,
`L = Σ e_i`, has `A_j(w) = (L − w + 1)₊ − Σ_i (e_i − w + 1)₊` sj-crossing starts (a window crosses no sj iff it lies
in one exon). The library's `T(w) = Σ_t θ_t (L_t − w + 1)₊`, `A(w) = Σ_t θ_t A_j(w, t)`, `θ` molar (copies, not
fragments: `A_j` already counts starts); production takes `θ` uniform over the non-synthetic transcripts.

**4.3 The four gDNA pools** (`GdnaOpportunity.combined_probability`). Two contained (§1.3), two crossing one boundary
(§1.4), `T(w) = Σ_refs (L_ref − w + 1)₊`:

    f(w) ∝ [Σ_p count_p(w)] · T(w) / [Σ_p A_p(w)]

the opportunity-weighted mean `Σ_p A_p f_p / Σ_p A_p` of the de-tilted pools, exactly inverse-variance under Poisson
counts (`Var count_p ∝ A_p`); not one divisor on the pooled histogram (`TRAPS: opposite-tilts-must-not-pool`). It is
gDNA's law where §4.4's contrast declines.

**4.4 The two-pool contrast** (`fl.deconvolved_gdna_counts`). The intergenic and intronic contained pools share one
geometry, so de-tilted by their own opportunities they are two shapes at purities `a_p = min(ρE_p/n_p, 1)` (`ρ` §7):

    f_0 = a_0 g + (1−a_0) r,   f_1 = a_1 g + (1−a_1) r   ⇒   î = [(1−a_1) f_0 − (1−a_0) f_1] / (a_0 − a_1)

`î` is clipped at zero and renormalised, then `g = m + λ(î − m)` with `m = w_0 f_0 + w_1 f_1` (`w_p` the count shares)
and `λ = sep²/(sep² + SE²)`, `sep = a_0 − a_1`, `SE² = a_0²/n_0 + a_1²/n_1`; `g = m` at `sep = 0`. Under capture
`a_0` clips at 1 and `g = f_0` exactly. With no density, an empty pool or nothing after the clip it declines and §4.3
stands in. Both forms assume uniform placement over the template; gDNA under capture is not, a placement residual no
divisor fixes. The law is then shrunk to the all-fragment law, `(counts + ν·global)/(n + ν)`,
`ν = CONSTANTS.fragment_length.pool_prior_ess` (`_smooth_eb`).

**4.5 RNA's law from its own counts** (`fl.own_count_law`). The sj pool's raw count is

    c(l) ~ Poisson(N π(l) law(l)),     law(l) = m(l) e^{d(l)} / Σ_u m(u) e^{d(u)}

with `π` §4.2's, `m` the all-fragment law (`payload.deposited_lengths`), `d` a random walk over the lengths `m` holds,
`d(l+1) − d(l) ~ N(0, q)`. At given `q` the MAP minimises the convex

    f(d) = Σ_{l: π>0, m>0} [μ(l) − c(l) log μ(l)] + Σ_l (d(l+1) − d(l))²/2q,     μ(l) = N π(l) m(l) e^{d(l)}

(Hessian `H = diag(μ) + K/q`, `K` the path's first-difference precision; Newton with step halving), and `q` maximises
`log p(c | q) ≈ −f(d̂) + ½[log n − (n − 1) log q] − ½ log |H|`, `n` the path length. A length with `π > 0` and no count
is evidence; `π = 0` is not; an empty pool gives `law = m`. The variance is the normalised law's by the delta method,
`J = diag(law) − law lawᵀ`, `H = UᵀU`:

    Var law(l) ≈ [J H⁻¹ Jᵀ]_ll = law(l)² |U⁻ᵀ (e_l − law)|²

`d`'s level cancels (`J 1 = 0`); a one-length law reads `law(l)² [H⁻¹]_ll`. Conditional on `q`. Gate:
`tests/calibration/test_fl.py`.

---

## 5. Strand

**5.1 The strand term** (`transfer_rows.h`'s `strand_term`, `strand_variance`; readable form
`tests/calibration/_psi_reference.strand_loglik`). At one cell of the `(f_g, f₊, f₋)` simplex,

    p = ½ f_g + κ f₊ + (1−κ) f₋
    var = n p_ref(1−p_ref)
    log L = −½(u₊ − n p)²/var − ½ log var

`var` is frozen at the slot's reference composition: the count sets precision, never composition. A Gaussian
approximation; the message claims use the exact conditional binomial (§3.2). The variance is binomial (§6).

**5.2 What it can say.** With tilt `d = f₊ − f₋`, `p = ½ + (κ−½) d`: gDNA cancels, reached only through
`f_g ≤ 1 − |d|`. At a single-strand slot `dp/df_g = ½ − κ`, so the Fisher information on `λ = logit f_g` is

    I_λ = n (2κ−1)² [f_g(1−f_g)]² / (4 p(1−p))           exactly 0 at κ = ½

**5.2b The channel is live iff the protocol preserves strand** (`strand_balance.strand_discriminability`; the kernel's
`strand_evidence` is `τ_λ`, read in §9c):

    disc = 4(κ̂−½)² if BF₁₀ > 1 else 0,   τ_λ = N · disc · [f_g(1−f_g)]²/(4p(1−p))

`BF₁₀` compares, on the spliced 2×2 (`n_same` of `N`, `strand_balance.fit_strand_balance`), κ free under Beta(1, 1)
against `κ = ½` exactly:

    ln BF₁₀ = N ln 2 + ln B(a, b),   a = κ̂(N+2) = n_same + 1,   b = (1−κ̂)(N+2)
            ≈ ½[z² − ln(2N/π)],      z = (κ̂ − ½)/√(¼/N)

Live iff `BF₁₀ > 1`, no constant: the excursion a fluctuation must clear grows as `√ln N`. With no spliced observation
the marginals are equal and `κ = ½`. gDNA's strand mean is ½ by symmetry and enters nowhere.

**5.3 gDNA's strand term is ½** — double-stranded, no sense direction.

**5.4 κ, the genuine junctions' rate** (`junction_fit.genuine_sense_fraction`). Junction `j` has `k_j` of `n_j` reads
on the minority side, oriented so the genuine rate is below ½:

    ℓ(κ, w) = Σ_j log[ w_g Bin(k_j; n_j, κ) + w_a Bin(k_j; n_j, ½) + w_r Bin(k_j; n_j, 1−κ) ]
    κ = (m + 1)/(R + 2),   m = Σ_j [r_gj k_j + r_rj (n_j − k_j)],   R = Σ_j (r_gj + r_rj) n_j

at the maximum, `r_gj`, `r_rj` the genuine and reversed responsibilities: the Beta(1, 1) mean over the reads credited
to RNA. For fixed κ the objective is concave in `w`, so KKT finds the maximum at a vertex, an edge or by a log-barrier
ascent; `max_w ℓ` is scanned in `log(κ/(½ − κ))` and refined by Brent's method (`CONSTANTS.junction_fit`). All
junctions genuine gives `(n_same + 1)/(n_obs + 2)`. Consumers see `n_same' = κ(n_obs + 2) − 1` clipped to
`[0, n_obs]`, and the raw `n_same` where the raw split reads no protocol (`StrandModel.genuine_n_same`).

---

## 6. Overdispersion is 0, and strand counts cannot measure it

The strand term has no overdispersion term. With intra-class correlation `od`, `Var(k) = n μ(1−μ)[1 + (n−1) od]`.
**6a.** At gDNA's mean ½, exactly, for `n ≤ 3`: `BB(n, ½, ρ) = (1 − ρ) Bin(n, ½) + (ρ/2) δ₀ + (ρ/2) δₙ`, and `δ₀`, `δₙ`
are what pure single-strand RNA reads as κ → 0. So gDNA with a fraction ρ of objects pure RNA has the likelihood of
pure gDNA at `od = ρ`: composition heterogeneity reads as overdispersion unless composition is known from elsewhere.

---

## 7. Background rate

`gdna_density.one_sided_rate`. Transcription only adds fragments, so every density `n_i/E_i` over-reads gDNA and the
clean objects are on the low side. Under pure gDNA `n_i ~ Poisson(ρE_i)`, and De Moivre's identity is exact:

    E[(λ − N)₊] = ½ E|N − λ| = λ^{⌊λ⌋+1} e^{−λ}/⌊λ⌋!  =: D(λ),        F(ρ) = Σ_i (ρE_i − n_i)₊ − Σ_i D(ρE_i) = 0

solved by bisection: `F(0) = 0`, `F(0+) < 0` once any object is counted, `F` unbounded above, one positive root.
Contamination only lowers `(λ_i − n_i)₊`, so the bias is one-signed, high. The pooled `Σn/ΣE` is reported beside it.
A region resolves a rate only above `~1/E` (Fisher information `ρE`); `ρ ≈ 0` does not prove absence.

---

## 9. Priors on a grid

The kernel's `JEFFREYS_REF`. `½ log f_g + ½ log(1−f_g)` on the λ grid is exactly
Beta(½, ½) in `f_g` — strength 1, mean ½ — and the only grid-width-stable choice (`TRAPS: no-prior-means-haldane`).
A fitted landscape's gDNA arm is added on top of both halves (`ISSUES: the-gdna-prior-enters-psi-twice`).

**9a. A vertex is unreachable without evidence.** A prior with a density has its median strictly inside `(0, 1)`, and
with no composition evidence posterior = prior, so the read-out (the median) cannot reach a vertex at any depth. On an
unstranded library §5.2's information is 0, so that shortfall is depth-independent.

---

## 9b. The EM's RNA prior goes to every RNA component, in proportion to its evidence

`native/em_prior.h`'s `apply_grouped_prior_update`; the pseudocounts are §9b.3's.

    rna_count = Σ_{i≠g} raw[i]
    out[g] = raw[g] + gdna_prior,     out[i] = raw[i] · (1 + rna_prior/rna_count)   for every i ≠ g

The RNA pool sums to `rna_count + rna_prior` and every RNA component, annotated or synthetic, takes the same factor:
the prior sets the gDNA:RNA split and nothing within RNA (per M-step, not end to end). Where `rna_count = 0` the prior
is shared by the carried state (VBEM's α), and dropped if that is 0 too.

**A shadow span against gDNA.** A synthetic span holding nothing decays against the transcript it shadows, being
longer. Against gDNA, with `m` fragments whose only RNA candidate is the span, it grows iff `m w_N > Total θ_g w_g`. On
`N` fragments only gDNA and the span explain, split evenly by genome strand, the fixed point in `r = b/a` (span over
gDNA) at strand specificity `ss` is

    r · L_n/L_g = [ ss·r/(ss·r + ½) + (1−ss)·r/((1−ss)·r + ½) ] / [ ½/(ss·r + ½) + ½/((1−ss)·r + ½) ]

with only root `r = 0` at `L_n ≥ L_g`. At `ss = 0.99` the span's share is 0 at `L_g/L_n = 1`, 35.2 % at 2, 52.7 % at
6.2, 82.0 % at 20 — the closed form `tests/test_estimator.py` pins the solver against.

### 9b.1 The prior is an additive per-component pseudocount — the weights are the design

`out[i] = raw[i](1 + P/R) = raw[i] + a_i`, `a_i = P·w_i/Σw` at `w_i = raw[i]`: an allocation echoing the current
posterior, hence neutral within RNA. `AggregatePrior::component_rna_prior_weight` generalises `w`; `pipeline.py` passes
none (`ISSUES: per-transcript-prior-lane`). At `w_i = raw[i]`, `out[i] = 0` is absorbing: no prior revives a component
with no warm-start evidence, while a flat `P/n` (`P` a fragment count) would lift every component past VBEM's
activation point, reviving every shadow span.

### 9b.2 The prior cancels exactly under MAP; under VBEM the digamma residual is left

Under MAP `θ_i ∝ α_i` and `c = 1 + P/R` divides out of every within-RNA ratio. Under VBEM `θ_i ∝ exp ψ(α_i)`,

    ψ(c α_i) − ψ(α_i) = log c + (1 − 1/c)/(2 α_i) + O(α_i^{−2})

and the residual after `log c` normalises away is bounded by `1/(2α_min)`, vanishing as the locus deepens.

### 9b.3 The pseudocounts' odds are calibration's gDNA share of the locus

`pipeline.em_pseudocounts`, formed in `run_locus_em_partitioned`; gate `tests/test_em_pseudocounts.py`. The
deterministic fragments (one candidate, an annotated junction) skip the E-step but enter every M-step, so θ_g is
gDNA's share of all `N` fragments — as it must be, since a transcript's length counts its spliced starts. The grouped
update is exact EM for `log L(θ) + P_g log θ_g + P_R log(1 − θ_g)`, which leaves the EM's fixed point in place iff
`P_g : P_R = G : (N − G)`. With `G_c` calibration's gDNA count of the locus (§11), `U` the units with a gDNA candidate,
and `N_c` the units neither multimapper nor splice artifact, plus the deterministic fragments:

    P_g = U · min(G_c/N_c, 1),     P_R = U − P_g,     no prior where N_c = 0

neutral iff `G_c/N_c = G/N`; calibration's RNA count never enters. `G_c` counts multimappers at their second-pass draw
and `N_c` does not (`ISSUES: the-em-prior-share-excludes-drawn-multimappers`). The strength `U` is not derived.

---

## 9c. ψ's reference and its curvature

The reference is §9's Beta(½, ½); ψ ships no reference location. **9c.2.** `τ_λ` (§5.2b; `solve_kernel.cpp`'s
`strand_evidence`) reads the strand counts alone: the data's Fisher information, which a prior's curvature is not
(`TRAPS: a-priors-curvature-is-not-the-datas-information`). Near a vertex the likelihood is genuinely flat in λ —
moving the mode from `f_g = 0.98576` to `0.99975` lowers `[f_g(1−f_g)]²` by 3,154× — and a reference's curvature fed
in would release the whole count precision at empty slots.

---

## 9e. ψ's θ quadrature — the nodes follow the strand term's peak

`native/psi_kernel.h`'s `tilt_window`; `simplex_logodds.TILT_NODES`. At an AMBIG slot ψ is read on the θ-marginal
`M(λ) = ∫ exp ψ(λ, sin θ) dθ`, `θ ∈ [−π/2, π/2]`; the tilt's arcsine measure `(1−τ²)^{−½}` is cancelled by
`|dτ/dθ| = cos θ`, so θ carries no weight. With `a(λ) = (1 − f_g)(κ − ½)`, `p = ½ + a τ`, so at fixed λ the strand term
is an exact Gaussian in τ:

    L(λ, τ) = −(τ − τ̂)²/2σ_τ² + const(λ),     τ̂ = d/a,   σ_τ = σ_p/|a|,   d = u₊/n − ½,   σ_p = √var/n

of θ-width `σ_θ ≈ 1/(2√n |κ−½| (1−f_g) cos θ̂)`; a fixed lattice coarser than `σ_θ` reads it as a comb.

**The rule.** With `τ_m = clip(τ̂, −1, 1)` the maximum on the domain, the window within `T` nats of it is

    ρ = √((τ_m − τ̂)² + 2σ_τ² T),     [τ_lo, τ_hi] = [τ̂ − ρ, τ̂ + ρ] ∩ [−1, 1]

for a peak inside, on and beyond the boundary alike. `K_t` uniform θ nodes over `[arcsin τ_lo, arcsin τ_hi]` take
trapezoid weights, written into ψ as `log h(slot, λ)` plus `log ½` at the ends. In θ the integrand is smooth and even
about `±π/2`, so the trapezoid rule is spectrally accurate at an interior window end (`e^{−T}` of the peak) and at a
domain end (the periodic rule's `½` on the reflected window). `log h` carries the marginal's own λ-dependence. No strand
information (`κ = ½`, `n = 0`) gives the whole domain.

**The constants.** `T = −log ε₆₄ ≈ 36` nats (`ε₆₄ = 2⁻⁵²`): the mass outside is `erfc(√T) < ε₆₄`. The interior window
is `2√(2T)σ_θ` wide and the trapezoid error on a Gaussian is `2e^{−2π²(σ_θ/h)²}`, below `e^{−T}` once
`h ≤ σ_θ π√2/√T`, so `K_t − 1 ≥ 2T/π`: **`K_t = 24`**
(`test_vertex_reference.test_the_theta_marginal_matches_adaptive_quadrature_at_every_depth`).

**A delivered row is evaluated at the nodes** (`profile_of_level`; readable form `_psi_reference.row_at`). A
`CubeRows` row holds each strand's profile over `u = log(ρ/ρ_ref)`, the slot's total `n` and RNA opportunity `a_r`; at
each cell strand `s`'s share `f_s = (1 − f_g)(1 ± τ)/2` implies `ρ_s = f_s n/a_r`, read at `log(ρ_s/ρ_ref)`.

**A property of the exact marginal.** At an interior tilt `M(λ) ∝ σ_τ(λ) ∝ 1/(1 − f_g)` until the band fills the
domain: an Occam factor of order `√n` toward gDNA at a balanced slot with no other evidence.

## 9f. The tilt atom — the AMBIG tilt's hypotheses are {pure +, pure −, mixed}

`native/psi_kernel.h`'s `slot_cube`. A slot whose RNA is all on one strand has its truth at the strand cap, but a
continuum of slightly impure tilts fits the same `p̂` at lower `f_g`, and §9e's volume factor pulls the median below
it. Presence per strand is discrete, so the measure is

    P(τ) = ⅓ δ(τ − 1) + ⅓ δ(τ + 1) + ⅓ (arcsine on (−1, 1)),
    M(λ) = ⅓ [ e^{L(λ, +1)} + e^{L(λ, −1)} + (1/π) ∫ e^{L(λ, sin θ)} dθ ]

In ψ: two more columns, `τ = ±1` at log-weight 0, beside the 24 windowed nodes whose log-weights carry `−log π` (with
no strand information the window's weights sum to π and the three hypotheses carry equal mass —
`test_vertex_reference.test_the_three_tilt_hypotheses_carry_equal_reference_weight`). Every θ-independent term cancels
among the hypotheses; `f_g` is still the θ-marginal's median, `w₊` the RNA-mass-weighted share over every column.

**The witness.** A pure column also explains an interior tilt at the cap. A held RNA level delivered on strand `s`
certifies `s` carries RNA, so the pure `−s` column is `−∞` wherever the slot's row holds a profile on `s`; with
nothing delivered both atoms stand (`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`).

---

## 10. The second pass's score

`second_pass.py` (`score_held_fragments`, `combine_factors`, `placement_abundance`, `choose_hypotheses`). It ranks one
held fragment's enumerated candidates; it claims no composition.

    score(h) = ρ(h) × f(L_h) × s(h)

`ρ(h)` is §2's reciprocal deposit, the bottleneck (minimum) of `sj_inv_length_sum` over a spliced hypothesis's
junctions or of `boundary_unspliced_inv_length_sum` over the boundaries the genomic hypothesis crosses where the
others jump — one deposit rule, one scale. `f` is RNA's law for a hypothesis with an intron or annotated junction,
else the all-fragment law; `s` the observed orientation's probability under the hypothesis. Factors are normalised in
the fragment, applied in order, skipped when flat-zero among survivors (`TRAPS: an-all-zero-factor-is-inert`), and
formed in logs with the record's maximum subtracted. One multinomial draw per fragment in the buffer's canonical order
from one RNG stream, clamped onto the last candidate with mass; never keyed on content, which would draw identical
fragments identically and collapse a 60/40 posterior to 100/0.

**Certified RNA.** A placement whose observed intron resolves to an annotated junction (`Accumulator::sj_edge_id`, the
motif not ambiguous) is booked spliced, so every path takes the RNA law and RNA's strand term — by the motif if the
aligner wrote one, else by its junctions' transcripts' strand (½ when both). `ρ`'s hard zero is the large-exposure limit
of `P(0 | λ, E) = e^{−λE}`.

**On a multi-placement record (§1.8)**, `score(p, h) = e(p) × A_p × ρ(h) × f(L_{p,h}) × s(p, h)`: `e(p)` the hit's
sequence at the EM's prices (§16), at the library's own `ε`; `A_p` the pass-one unique traffic where the placement's
observed path deposits — the bottleneck of `sj_inv_length_sum` over the annotated junctions its introns resolve to,
else of `boundary_unspliced_inv_length_sum` over the boundaries it crosses, else
`region_contained_count / contained_eff_length(ℓ, global_pmf)` in its one region, 0 in no object. Each is a start
density per base, and neither factor is applied on a one-placement record. All placements at `A = 0` make the factor
uninformative; a zero beside a positive stays hard.

**gDNA's floor** (`pass_one_bulk`). gDNA is one template at a uniform rate before capture and capture only enriches.
Pass one deposits a start only when its fragment is not held, a fraction `u_p`, so `A_p = k_p/S_p` estimates
`λ_p u_p`; over the intergenic regions the same estimator gives the bulk `ρ̂_b = Σk/ΣS ≈ ρ_G ū`. Every placement gDNA
can be at (none booked spliced) reads `A_p ← max(A_p, ρ̂_b)`, so gDNA-only copies at or below the intergenic
mappability split evenly and an expressed copy keeps its traffic. The exact `k_p/(u_p S_p)` needs mappability.

---

## 11. The conserved frame — every component's shares of its objects at their capture efficiencies

`effective_length.conserved_cut_shares`, `capture_eff_length`, `capture_efficiency.capture_efficiencies`,
`priors.assemble_priors`, `effective_length.observable_cut_shares`. Under capture the EM divides each component — the
locus gDNA, each synthetic span, each transcript — by the sum over the objects its fragments deposit on of its
CONSERVED SHARE of the object times the object's capture efficiency:

    L_T = Σ_pieces S_T(p) · c_p + Σ_cuts M_T(k) · c_k

PIECES are the regions the template covers; between consecutive pieces is a CUT (a boundary or, in a spliced
template's coordinates, a junction). Objects and not bases, because the deposit rule gives a crossing fragment's unit
to the objects it crosses.

**The shares.** `S_T(p) = E_f[(ℓ_p − w + 1)₊]` (§1.3); `M_T(k)` is §3b's deposit summed over placements. With the cut
at 0 between pieces `a` (left) and `b` (right) and reach `near ≥ a`, `far ≥ b`, a placement with `x` bases left is
admissible iff `max(1, w − far) ≤ x ≤ min(w − 1, near)`, `w ≤ near + far`, and gives `[L_a(x) + L_b(w − x)]/w`:

    L_a(x) = x  (x ≤ a: an end slice, bounded by this cut)        L_a(x) = a/2  (x > a: runs through the piece)

**The closed form** (`_side_share`). The left inner sum with no far side is `(w−1)w/2` to `w = a+1`, `a w/2` to
`near+1`, `a(near+1)/2` beyond; the far side removes `(w−far−1)(w−far)/2` to `far+1+a` and `a(w−far)/2` beyond. Over
`f` up to `near + far` each segment is a difference of cumulative sums of `f`, `w f` and `f/w`: O(objects). The right
part is the mirror; a reach excluding its own piece is refused. Unbounded,
`M = ½E_f[min(a, w−1)] + ½E_f[min(b, w−1)]`, §3b's `A_mass` in expectation.

**Conservation.** A placement crossing no cut gives its unit to its piece; one crossing cuts gives its whole unit to
cuts, so `Σ_p S_T(p) + Σ_k M_T(k) = Σ_w f(w)(L_T − w + 1)₊`, the fl-marginal length, exactly. A total never gates the
rule — a `1/K` split conserves the same sum (`TRAPS: conservation-misses-mis-attribution`) — so the gates are per
object against the reference accumulator (`tests/calibration/test_effective_length.py`, `test_capture_eff_length.py`).

**The gDNA component** (`assemble_priors`). Reach `UNBOUNDED_REACH`; a boundary's share `M_e` is
`CalibrationResult.gdna_boundary_conserved_len`. Over the locus's regions and every boundary touching them, at the
count's overlap weights `ω`:

    gdna_eff_len = Σ_r ω_r S_r c_r + Σ_e ω_e M_e c_e,         gdna_count = Σ_r ω_r k_r + Σ_e ω_e k_e q_e

`q` (§3b) turning crossing counts into fragment counts. At every efficiency 1, for a locus owning its regions whole,
`Σ S_r + Σ M_e = L + ½E_f[min(a_L, w−1)] + ½E_f[min(a_R, w−1)]`: the starts whose fragments overlap the locus, each
once. The crossing support `E_f[w−1]` counts a start at every boundary its fragment crosses and would read the density
`ρ q̄`. The share is gDNA's own geometry, never the pooled `q` (`tests/calibration/test_priors.py`).

**Transcripts and synthetic spans** (`transcript_capture_eff_lengths`). Pieces and cuts come
from the annotation and the RNA law; a cut is a junction where its pieces are different exons with a region between.

    eff_em_t = fl_t · (Σ_p S_p c_p + Σ_k M_k c_k) / (Σ_p S_p + Σ_k M_k)

which conservation makes the sum itself wherever `fl_t` is unfloored; a transcript with no share keeps `fl_t`.

**The junction, by conservation of bases** (`_junction_prices`). gDNA never deposits on a junction. With `c_lo`, `c_hi`
its low and high boundaries and `c_intron,lo`, `c_intron,hi` the intron pieces beside them, and capture additive over a
fragment's bases: a gDNA crossing of the low boundary holds the junction fragment's low-exon bases plus `w − x` intron
bases, one of the high boundary `x` intron bases plus its high-exon bases; on average each extra is half an intron
fragment. Capture saturates (a fragment binds its best probe part), so

    c_junction = min( max(c_lo + c_hi − ½(c_intron,lo + c_intron,hi), 0),  max(C, c_exon,lo, c_exon,hi, c_lo, c_hi) )

`C` the junction cap. Where each boundary reads its flanks' mean it is `½(c_exon,lo + c_exon,hi)`. An intron piece
too short for a gDNA fragment reads the boundary on its far side, so a one-region intron gives `½(c_lo + c_hi)`.

**The junction cap, from RNA's own split** (`fit_junction_cap`). The saturation level `S` (below) is gDNA's, and gDNA
never deposits on a junction, so `S` bounds only what gDNA can show. RNA crosses the junction. In an exon region `A`
RNA and gDNA are the same placements under the same probes, so a transcript at rate `ρ_t` leaves `ρ_t S_t(A) c_A`
contained fragments in `A` and `ρ_t M_t(J) p_J` on a junction `J` priced `p_J`: the ratio is free of `ρ_t`. Over the
transcripts using `J`, with `A` covered by them alone and their `M_t(J)` read at the mean `X_J`,

    E[W_J] / E[m_J] = X_J p_J / D_J,        D_J = Σ_{A beside J, its own} S(A) c_A

`W_J` the junction's conserved spliced mass (`sj_conserved_mass`) and `m_J` calibration's unspliced RNA in those
regions (`count_rna_region`). As counts, with `p_J(C)` the rule above at cap `C`:

    W_J | W_J + m_J  ~  Binomial( q_J(C) ),        q_J(C) = X_J p_J(C) / (X_J p_J(C) + D_J)

    C = argmax_{C ≥ S}  Σ_J [ W_J log q_J(C) + m_J log(1 − q_J(C)) ]

over the junctions with a region of their own, a fragment in their footprint and a price the cap can move
(conservation's sum above `max(S, the best object beside)`). A term is constant in `C` below its junction's best
object and above its sum, and concave in `log C` between, so the sum is not concave: it is maximised by search on a
geometric grid from `S` to the largest sum (`CONSTANTS.calibration.junction_cap_grid_points`) and the first maximum
is taken, so a cap the counts cannot tell from `S` is `S`. The cap is one-sided (`C ≥ S`), and with no level above the
bulk (`S = 1`) it is 1. The shares are the ruler's own, carried to gDNA's table as an efficiency is (§14), so the fit
and the lengths it feeds state one model. What it assumes without checking is `ISSUES: the-junction-cap-reads-bookings`.

**The efficiency** (`capture_efficiencies`). Each object reads its gDNA count `k` on its support `S` (a region's
contained, a boundary's crossing) as `ℓ(x) = k log(S eˣ) − S eˣ` over `x = log ρ`, against the fitted landscape as
Efron's two groups — the bulk basin collapsed onto its peak `x_b` with its whole mass, elsewhere the unpadded density
`DensityLandscape.logP_free`:

    log c_o = E[min(x, s(x)) | k_o, S_o] − x_b

`s(x)` the peak of the enriched basin holding `x` (`+∞` outside one): capture saturates per level. A mean of the LOG
density, so a lone level of mass `μ` at `Δ` nats above the bulk lifts an evidence-free object by `e^{μΔ}`, not
`1 + μ(e^Δ − 1)`. `S = 0` reads 1; `k` is a deconvolved mass, the Poisson continued by the gamma function.

**The bulk, the saturation, the unit** (`landscape.census`, `split_basins`). The grid splits into basins at the minima
between maxima; the bulk is the most massive, since gDNA off the probes outnumbers gDNA on them, and every basin above
it is an enriched level. `S = exp(x_sat − x_b)` with `x_sat` the most massive enriched peak
(`CalibrationResult.gdna_capture_saturation`; 1 if none). Efficiencies are relative yields against the bulk
(`gdna_bulk_density`), so a factor common to every component cancels at convergence
(`ISSUES: the-em-stops-before-it-converges`). One basin makes every efficiency exactly 1; no landscape publishes no
bulk, and `transcript_capture_eff_lengths` returns `fl` while `assemble_priors` reads `Σ S_r + Σ M_e`.

**Where the rule is exact.** For gDNA under any capture physics; for another component on the objects where its
fragments are gDNA's — contained pieces and cuts away from junctions, exon edges and template ends. Near an edge or end
its fragments splice or stop where gDNA's run on, and a probe captures in genomic coordinates on gDNA but transcript
coordinates on cDNA (`ISSUES: ruler-witness-geometry-on-transcript-panels`).

---

## 12. The flux price's witness — the column count on the protocol's share of the opportunity

`transfer_kernel.h`'s `rna_lane`, `lane_receive`. A certified flux at an exon's junction — count `c_J` at route rate
`r_J`, in whole-strand units since `sj_count` and `route_rate` are keyed by the junction's strand — is priced by
`hop_price` (§3.5h) against the exon's witness: the unspliced count on the column strand `s`'s RNA reads on. Under
read rate `κ_read = max(κ, 1−κ)`,

    E[c_s] = κ_read R_s + (1 − κ_read) R_o + g/2

so on the exon's whole opportunity `a_r` an agreeing pair would read a false disagreement `log(κ_read)²`. The column is
an opportunity of `κ_read · a_r` for strand `s`'s RNA, so the witness is `(c_s, κ_read · a_r)`: at `κ = ½` the exon's
total `2c_s/a_r`, at `κ → 1` the column's. `κ` is the fitted model's (`κ_read = ½` with none).

**Local strand difference** (`split_live`). The column difference's mean is `(2κ − 1)(R₊ − R₋)`; the DNA cancels. Where
the protocol is live (§5.2b) and the recipient holds a witness, `lane_receive` prices the pair on `n_x = c_s − c_o`
at its Poisson variance `v_x = c_s + c_o`: `count_logvar(n_s) + count_logvar(c_s) + max(0, l² − (v_s/n_s² + v_x/n_x²))`,
`l = log((n_x/a_x)/(n_s/a_s))`, `n_s`, `v_s` the sender's held RNA count and variance.

---

## 13. Whole counts are a transportation problem

`em_assign.h`'s `assign_posteriors` (`assignment_mode = "sample"`); gate `tests/test_estimator.py`'s
`TestWholeCounts`. In one locus fragment `u` has candidates `C_u` with posteriors `p_uc > 0`; `n_c = Σ_u p_uc`,
`Σ_c n_c = m`. An ASSIGNMENT `z` sends each fragment to one candidate, with counts `N_c`.

**13.1 A rounding is always reachable.** In the network source → fragment (exactly 1) → candidate (0 to 1) → component
(`⌊n_c⌋` to `⌈n_c⌉`) → sink the posterior matrix is a feasible flow; every bound is an integer and a network matrix is
totally unimodular, so an integral flow exists: every `N_c ∈ {⌊n_c⌋, ⌈n_c⌉}`, the total exact.

**13.2 One particular rounding need not be.** For targets `q_c` (largest remainder, `Σ q_c = m`) an assignment with
`N = q` exists iff `Σ_{c∈S} q_c ≤ |U(S)|` for every component set `S`, `U(S) = {u : C_u ∩ S ≠ ∅}` (Hall). The
fractional counts satisfy it; rounding several of a small set up can break it (one fragment split ½ : ½ between `B` and
`C` beside four components at ¼ hands `B` and `C` the two largest remainders).

**13.3 The repair is augmenting paths, and it cannot cycle.** From the draw let `e_c = N_c − q_c`. A path
`c_0 → u_1 → c_1 → … → u_r → c_r` from `e_{c_0} > 0` to `e_{c_r} < 0`, each `u_i` on `c_{i−1}` with `c_i` a candidate,
moves every `u_i` one step: `e_{c_0}` falls by one, `e_{c_r}` rises by one, nothing else moves. (1) It terminates: each
path lowers `Σ_c max(e_c, 0)` by one. (2) It succeeds exactly when the targets are reachable: if no under-target
component is reachable from `c_0`, the reachable set `R`'s fragments have all their candidates in `R`, and
`Σ_R N_c > Σ_R q_c` violates 13.2. (3) Each path is one breadth-first search, `O(n)`.

**13.4 The fallback reaches 13.1's rounding.** On unreachable targets each component relaxes to `[⌊n_c⌋, ⌈n_c⌉]`: counts
above their ceiling shed along paths, then counts below their floor draw along backward paths. A path exists: else the
fragments on the reachable set `R` would number `Σ_R N_c > Σ_R ⌈n_c⌉ ≥ |{u : C_u ⊆ R}|`; the fill is the same with
`⌊n_c⌋`. Whole counts end within one of every fractional count, total exact, and ARE the largest-remainder rounding
whenever it is reachable.

**13.6 For identical fragments the draw is an urn.** The draw weights `c` by `p_c · owed_c / to_come_c`. Fragments
sharing a candidate set and posterior have `to_come_c = (left) · p_c`, so the weight is `owed_c/(left)`: an urn
without replacement, each fragment going to `c` with probability `q_c/m`, none stranding. Stranding comes only from
different candidate sets (`test_the_first_fragment_and_the_last_are_drawn_with_the_same_odds`).

---

## 14. The EM's two length tables — one shared law unless the data show a difference capture cannot make

`calibration/shared_length_law.py`; `capture_eff_length.rna_captured_fraction`, `_carried`; the variances in `fl`.
Gate: `tests/calibration/test_shared_length_law.py`.

**14.1 Capture cancels from the length odds.** With `P(x | X) = f_X(l) · w(s, l)/Z_X` and `w` the placement's capture
weight (probes bind sequence: one `w` for both origins), an unspliced fragment's odds are
`θ_R/θ_G · f_R(l)/f_G(l) · Z_G/Z_t`. Capture enters only the normalisers; under `f_R = f_G` equal tables are exact.

**14.2 The two measured laws are read in different frames.** `g` (`FLModels.gdna_pmf`) is uncaptured (§4.4); `r`
(`rna_pmf`) is the spliced census, captured. A law captured under response `m(l)` (non-decreasing, `m(l)/l`
non-increasing: overlap with probes grows with length and saturates) is

    f_m = f(1 + βm)/Z = (1 − κ') f + κ' f m/⟨m⟩_f,     κ' = β⟨m⟩/(1 + β⟨m⟩)

so on equal chemistry `r = g_m ≠ g`: a false gap.

**14.3 κ, RNA's captured fraction.** An object at gDNA yield `c` owes `1/c` of its fragments to the background:
`κ = Σ_exons n_RNA (1 − 1/c)₊ / Σ_exons n_RNA`, `κ < 1`, 0 with no landscape.

**14.4 The capture mixture** (`capture_mixture`). `G_c = (1 − c) g + c g l/μ_g`; `ℓ(c) = Σ_l r(l) log G_c(l)` is concave
with `ℓ'(0) = μ_r/μ_g − 1`, so `c = 0` when RNA's census is no longer than `g`; else `c` maximises `ℓ` on `[0, κ]` by
bisection on the slope.

**14.5 The bound is the worst case over every proportional response** (`capture_bound`). With
`ρ_m = dg_m/dg = 1 − κ + κ m(L)/⟨m⟩`, §14.2's premise gives `m(L)/E m(L) ≤_cx L/E L` (Lorenz order,
Jakobsson–Fellman), preserved by the affine map, so `E_g φ(ρ_m) ≤ E_g φ(1 − κ + κL/μ_g)` for convex `φ`: every
f-divergence is at most the size-biased mixture's. For Hellinger, `H(g, g_m) ≤ B`, `B² = H²(g, (1 − κ)g + κ g l/μ_g)`,
non-decreasing in `κ`. A faster-than-proportional response can exceed it
(`test_a_response_rewarding_length_faster_than_proportionally_can_exceed_the_bound`).

**14.6 The decision, on a metric** (`trust`, `geometric_path`). `H²(a, b) = Σ_l (√a_l − √b_l)²` is a squared L2
distance, so `H` is a metric (the root of the symmetrised KL is not). On equal chemistry `H(g, r) ≤ B`; otherwise the
chemistry distance is at least `D − B`. With `D² = (H²(g, r) − N)₊`,

    w = clip(D/B − 1, 0, 1)   (D = 0: w = 0;  B = 0 < D: w = 1),       T_G = normalize(r^{1−w} G_c^w),   T_R = r

**14.7 The noise term** (`sampling_noise`). The null is two noisy estimates of one law, expanded at `M = (g + r)/2`:
`N = Σ_l (Var ĝ_l + Var r̂_l)/(4 M_l)` — each law's own delta method would divide a near-empty law's noise by nearly
nothing (`test_a_noisy_near_empty_tail_of_one_law_does_not_hide_a_real_gap`). RNA's variance is §4.5's, and an RNA law
no spliced fragment measured shows no difference (`test_an_rna_law_no_spliced_fragment_measured_shows_no_difference`).
gDNA's is its raw law's plug-in `law_l²/raw_l` (`_plugin_variance`) at the empirical-Bayes weight squared,
`(n/(n + ν))²`; through §4.4's contrast it is differentiated as built, with coefficients
`c₀ = (1 − λ) w₀ + [î > 0] λ(1 − a₁)/(s Σî⁺)`, `c₁ = (1 − λ) w₁ − [î > 0] λ(1 − a₀)/(s Σî⁺)`, `s = a₀ − a₁`, and
`Var = c₀² Var f₀ + c₁² Var f₁`; a clipped length carries no inverse variance
(`test_the_contrasts_variance_is_its_own_estimators_through_the_clip`).

**14.8 RNA's normaliser carries gDNA's efficiencies to the EM's tables.** §11's efficiencies are read on `g`, RNA's
normaliser scores with `T_R` and `T_G`, so each object enters at `ĉ_o = c_o · S_o(g)/S_o(T_G)` (`c_o` where
`S_o(T_G) = 0`):

    eff_t = fl_t · Σ_o S_o(T_R) ĉ_o / Σ_o S_o(T_R)

At one table this is `Σ_o c_o S_o(g)`, the normaliser exactly under the null (§14.1); at `T_G = g` the factor is 1.

---

## 15. A multimapper's places

`scoring.cpp`'s `merge_places`, `flush_mm_group`; `capture_eff_length.place_efficiencies`, `place_log_weights`;
`em_locus.h`'s `IntergenicTerm`, `intergenic_log_weight`; `priors.intergenic_gdna_rate`.

**15.1 One fragment, several places.** The molecule came from one component `c`, at one place, with one length:

    P(u | c) = Σ_{p ∈ S_u(c)} w(p) · f_c(ℓ_{c,p}) · s_c(p) · e(p) / Z_c

`w` the place's capture weight (§15.2), `Z_c = L̃_c` (§11), `S_u(c)` the places `c` can be at: a transcript's candidate
places, the locus gDNA's gDNA-explainable places, the background's intergenic places (§15.3). For a transcript or the
locus gDNA, hits whose extents overlap on one reference are alternatives at one place and take the larger term;
distinct places sum; each intergenic hit is its own place. A `−∞` term is skipped, so hit order cannot matter. On a
unique fragment `w` is common to every candidate and is not applied.

**15.2 The place weight is the deposit rule's mass at the reader's efficiencies.** With `m_o(x)` placement `x`'s
deposit mass at object `o` (`Σ_o m_o = 1`) and `c_o` §11's efficiency,

    w(x) = Σ_o m_o(x) c_o,         Σ_x f(ℓ_x) w(x) = Σ_o c_o Σ_x f(ℓ_x) m_o(x) = Σ_o c_o share_o(f)

so `w` sums over a component's placements to its normaliser on the law the efficiencies are read on (a bottleneck,
maximum or per-base mean would not). A junction takes §11's price (RNA only); an implicit place is weighed on its
observed path. On `g`, gDNA's sum is calibration's supports and a transcript's at one table is §14.8's
(`test_a_templates_place_weights_sum_to_its_length_on_the_uncaptured_law`). Where a unit's hypotheses imply one length
the scorer's table factor cancels; where they differ the ratio is the shared law's per-fragment one.

**15.3 The intergenic background: a fixed rate in absolute units.** An intergenic place lies in no locus, so its gDNA
arrives at a rate `λ` no locus parameter carries. The locus emits a Poisson total `M` with composition `θ` (the
grouped prior on θ_g, silent on `M`). With `q_uc` the per-candidate term (place weights included), `h_G` gDNA's
length table and the deterministic fragments `d`,

    b_u = Σ_{p ∈ IG(u)} h_G(ℓ_p) · ½ · e(p)                    (w = 1 at an intergenic place)
    ℓ(M, θ) = Σ_u log[ Σ_c M θ_c q_uc/L̃_c + λ b_u ] + Σ_d log(M θ_{t_d}) − M + P_g log θ_g + P_R log(1 − θ_g)

E-step `r_uc ∝ M θ_c q_uc/L̃_c`, `r_u,IG ∝ λ b_u`; θ by the grouped update (§9b) on `Σ_u r_uc + n_det,c`; and
`∂ℓ/∂M = 0` gives `M = Σ_u Σ_{c ∈ locus} r_uc + N_det = N − Σ_u r_u,IG`, the locus's own data count. In the locus's
frame the background is one more column of weight `log λ − log M + log b_u`, outside the prior and θ's normalisation,
never pruned; `M` sits in SQUAREM's state at its slot (`M/N` under MAP, `M` under VBEM). `λ = gdna_bulk_density`; with
no landscape, the intergenic regions' gDNA count over their contained support
(`test_the_em_reaches_the_reference_fixed_point`).

---

## 16. A hit's sequence

`scoring.sequence_log_penalties`, `scoring.cpp`'s `score_mm_alignment`; the scan's `n_unique_unspliced_mismatches` /
`_bases`. Gate: `tests/test_hit_sequence.py`.

**16.1 Each base, aligned or clipped.** At mismatch rate `ε` an aligned base reads its reference with probability
`1 − ε` and each other base with `ε/3`; a soft-clipped base is placed nowhere and reads `¼`. Against a hit whose every
base aligns and matches,

    log e(p) = NM · log[ε/(3(1 − ε))] + clip · log[1/(4(1 − ε))]

`NM` its edits and `clip` its soft-clipped bases, both mates'.

**16.2 One hit reads no sequence.** It gives every candidate the same `e(p)`, which cancels from responsibilities and
prunes; `e` is read only between a multimapper's places.

**16.3 ε, read on the library.** `ε = Σ NM / Σ aligned bases` over the unique unspliced hits, the binomial MLE; spliced
hits (the aligner's own cut) and multimappers (the question being priced) stay out. No mismatch makes one impossible
(`−∞`); no unique unspliced base gives no price; `ε ≥ ¾` is refused.

**16.4 A spliced reading beside its contiguous one.** `S` is a spliced reading whose junctions transcript `T` holds,
`C` a contiguous reading of the same fragment, which gDNA, the synthetic span and an intron-retaining isoform read:

    P(gDNA | u) / P(T | u) = [θ_G/L̃_G · f(L_C) · ½ · e(C)] / [θ_T/L̃_T · f(L_S) · s_T(o) · e(S)]

- `C` aligns every base `S` aligns: `e(C)/e(S) = [ε/(3(1 − ε))]^{m_s − m_c}`. At an exact tie the odds are density ×
  strand × length, and the split is the calibrated one.
- `C` clips `Δc` bases `S` aligns: `e(C)/e(S)` gains `[1/(4(1 − ε))]^{Δc}`. A +1/−1 aligner clips a tail only where its
  mismatches `m ≥ Δc/2`, so against the reference the clipped bases' ratio is at most `√(ε/(3(1 − ε)))` per base,
  below `¼` for every `ε < 3/19`: `¼` can only favour `C`.
- More than one junction: gDNA reads only a fully contiguous reading; a partial one is RNA's alone (an intron-retaining
  isoform). The aligner reports whole-fragment readings, so no per-junction rule is needed.
