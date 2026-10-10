# RNA-short repair and the next design — independent review

*Sandbox document (`docs/dev/`): provisional, not authoritative, cited by nothing outside the sandbox. 2026-10-06.*

**Read:**

- the working-tree diff against `a7b03103`, including the untracked `test_splice_out_opportunities.py`;
- `RNA_SHORT_IMPLEMENTATION_REVIEW.md`, and `RNA_SHORT_FIX_PLAN.md` §8–10;
- DESIGN §0c.3, §6b.4–6b.8 and §7.1–7.2;
- EQUATIONS §3.2–3.2b;
- the three named ISSUES entries;
- FRAGMENT_LENGTH_POSTMORTEM §2.

**Archive:** `~/Downloads/rigel_runs/prototypes/2026-10-06_rna_short_counts/`. Its receipts and gate logs were read, not re-run.

**Computed here:** two offline numpy passes over the archive's frozen native inputs (no solve; commands in the appendix).

**Not done:** no source change, build, benchmark, commit or push.

Tags: **[V]** verified here from code, receipts or the frozen inputs · **[P]** plausible, decided by the named experiment ·
**[C]** contradicts a statement in the brief.

## 1. Verdicts

### The count-map repair: keep the arithmetic; do not land it as written

The four map families are correct count-frame algebra for given component opportunities [V]. I checked:

- the direction of each map;
- the inverse pairing (TRANSPORT pulls back through the inverse map; SPLICE_OUT evaluates the sender at the forward map);
- the sign of the shared-density shift;
- the route-rate units (per-placement densities, so `u + s` is commensurate);
- the zero-opportunity limits, where each message tends continuously to vacuous.

Three problems stand between that and a landing:

1. **The inputs are in two frames (`two-length-frames`).** RNA opportunities are evaluated on the spliced census law, which capture
   selects. gDNA opportunities are evaluated on the uniform-frame law. On equal-chemistry captured libraries this
   manufactures opportunity ratios that the old identity maps ignored. The receipts show a pool shift whose sign
   follows the capture label on all 12 contaminated ladder rows. The deferred-row loss is the capture-ON half of
   that shift. Whether the two frames cause it is what `law-frame-arm` decides [P].
2. **Map family 3 bundles a second mechanism (`center-is-a-second-mechanism`).** The alternative-splice discrepancy center changes widths at
   equal opportunity, on every stranded library.
3. **The widths were not re-derived (`widths-of-the-old-centers`).** They still belong to the old centers, and the two directions of a face
   are not one factor.

Two regressions in the brief are mis-attributed:

- **g98 ss0.70 ON does not indict the maps.** It is a knife-edge in ψ's prior and readout that the repair trips (`g98-valley-knife-edge`).
- **The test-chromosome capture-OFF transcript losses are noise.** They sit on pools and genes that are identical, so
  they are EM-fork noise (`off-deltas-are-em-fork-noise`).

### The proposed next design: right evidence abstraction, wrong order, and the capture half is not yet an estimator

- **Profiles are the correct interface.** The contract `p_ij ∝ L_i(ρ_j) π_j` is exactly derivable (`profile-derived`).
- **But the first consumer must be ψ's landscape.** It rests on the same point-count-as-Poisson premise the plan has
  just falsified for the ruler, plus training weights that bias its mixing weights (`landscape-has-the-same-premise`). §10 instead fits a second
  population inside the ruler, which would give two priors for one quantity.
- **The prior-once partner is plausibly that population estimator.** It is not yet demonstrated, and it is decidable
  cheaply (`population-vs-oracle`).
- **Profiles do not settle the hidden prior influence (`hidden-prior-influence`) or the unstranded case (`unstranded-needs-messages`).** The plan must decide both.
- **`E_cap` is a specification, not an estimator.** Placement-level capture is not identifiable from object totals
  (`capture-identifiability`). Price the scalar transfer with the arms already designed before building anything (`transfer-pricing`). Any capture-aware
  opportunity inherits `two-length-frames`'s missing RNA uniform-frame law (`missing-rna-uniform-law`).

## 2. Findings, ranked

| # | Severity | Area | Finding | Status |
|---|---|---|---|---|
| `two-length-frames` | High | repair input | RNA opportunity on the census law, gDNA on the uniform law; the pool shift's sign follows the capture label on every contaminated ladder row | V ratios & signs · P cause |
| `g98-valley-knife-edge` | High | ψ prior/readout | g98 ss0.70 ON: identical landscape mass, a 1.3-nat deeper valley at the truth; strand alone recovers the truth | V |
| `center-is-a-second-mechanism` | Med-high | repair scope | the alternative-splice center change is a separate mechanism acting at equal opportunity | V |
| `widths-of-the-old-centers` | Medium | widths | route-sum variance, sensitivity-free forward blur, zero-flux asymmetry, `v_ratio` of the old center | V |
| `what-the-gates-cannot-see` | Medium | tests | the 48 gates cannot tell the route sum from `S/Er_b`; one orientation; deep and mode-only; perturbations bypass the builder | V |
| `zero-opportunity-lane-switch` | Low | topology | at exactly zero opportunity the gDNA level lane takes over the face | V |
| `off-deltas-are-em-fork-noise` | Medium | reading | test-chromosome OFF transcript deltas sit on identical pools and genes | V |
| `shared-crossings-fused-as-independent` | Medium | evidence | adjacent boundaries share most gDNA crossings at short exons, yet the exon fuses them as independent | V geometry · P impact |
| `capture-identifiability` | High | §7 | placement capture is unidentifiable; per-object gDNA averages and a library-level length shape are identifiable | V |
| `missing-rna-uniform-law` | High | §7 | the pre- vs post-capture law: the repair already mixes frames; RNA's uniform-frame law is the missing estimator | V |
| `landscape-has-the-same-premise` | High | §8/§10 | the landscape has the falsified premise plus biased training weights; §10 adds a second population | V code · P size |
| `profile-derived` | High (+) | §8 | exact `L_i` derived: a Gamma–Poisson marginal, not "ψ minus the arm"; every limit falls out | V |
| `hidden-prior-influence` | High | §8 | hidden prior influence: belief-frozen strand variance, the intron factory, training by posterior variance | V |
| `unstranded-needs-messages` | Med-high | §8 | own-evidence-only profiles reproduce the refused strand-only reference when unstranded | V |
| `uninformative-reads-the-mean` | Medium | §10 | an uninformative object reads π's mean, so in a captured library it gets a large weight | V |
| `psi-median-discontinuity` | Medium | ψ | ψ's median readout is discontinuous at bimodal posteriors and already moves calibration counts | V |

### `two-length-frames` — the repaired maps read two length laws in different frames

**Where.**
- [region_geometry.py:208-209](../../src/rigel/calibration/region_geometry.py#L208-L209) builds `eff_gdna` and
  `eff_rna` from `gdna_fl_pmf` and `rna_fl_pmf`.
- [pipeline.py:943,959](../../src/rigel/pipeline.py#L943) passes `FLModels.gdna_pmf` (the uniform frame) and
  `FLModels.rna_pmf`.
- `rna_pmf` is the spliced pool: junction-de-tilted, but capture-selected
  ([fl.py:168](../../src/rigel/calibration/fl.py#L168)). The census lengthening is documented in
  `ISSUES: the-scorer-reads-a-census-length-law`: RNA 212 → 228.7 bp under capture on the ladder.
- Boundaries use one library constant per component, at unbounded reach
  ([region_geometry.py:193](../../src/rigel/calibration/region_geometry.py#L193)).

**Frozen inputs** of the test chromosome, which has equal chemistry: both laws are N(206, 98) in
`test_reference.yaml`. The true ratio in the code's own model is therefore exactly 0 everywhere.

| condition | boundary log(Eg/Er) | regions: median / q95 / max | regions with Er < 50, mean |
|---|---:|---:|---:|
| g50 ss0.50 OFF | +0.007 | −0.002 / −0.000 / −0.049 (min) | −0.032 |
| g98 ss0.70 ON | −0.044 | +0.013 / +0.109 / +1.198 | +0.511 |

**Receipts.** Each line is gDNA est−true, repaired minus current; the nRNA pool moves the other way.

- **Ladder capture-ON** (g05/g50/g98 × ss0.50/ss0.99), gDNA: −9.6k, −1.6k, −86.0k, −5.9k, −108.6k, −5.2k.
  - gDNA falls on 6 of 6 rows.
  - The nRNA pool rises on 6 of 6.
- **Ladder capture-OFF** (same six), gDNA: +1.2k, +0.3k, +4.7k, +1.7k, +2.0k, +1.4k.
  - gDNA rises on 6 of 6.
  - The nRNA pool falls on 6 of 6.
- **The deferred stratum** (ladder 10.42 → 12.40) is the capture-ON half of this.
- **Test chromosome ss0.70 ON** moves real pools: at g50, mRNA −7.2k and nRNA +6.2k.

The sign follows the capture label. On equal chemistry, a correct repair would be a no-op in both halves.

**[C]** EQUATIONS 3.2b and the brief §3 say "All four opportunities are capture-blind geometry, each evaluated on its
component's length law". The geometry is capture-blind; the RNA law it is evaluated on is not. The RNA level lanes
already read `eff_rna` (`rna_lane`); the repair extends that dependence to every composition face.

**Mechanism [P].** Under capture the census RNA law is longer than the uniform gDNA law. That gives:

- `Er_b > Eg_b` at boundaries;
- `Er_x < Eg_x` at short regions.

The reverse splice face then places a boundary's gDNA odds below the exon's, by `log(Eg_x/Er_x) − log(Eg_b/Er_b)`.
The FORWARD face passes that to the intron. The net sign at each object depends on its length relative to the two
laws. The receipts show the net, not the per-object path.

**Falsifier.** On an equal-chemistry captured simulation, every face's map must be the identity to within the laws'
estimation tolerance. Today it fails on g98 ss0.70 ON. Decisive arm: `law-frame-arm`.

### `g98-valley-knife-edge` — g98 ss0.70 ON is a knife-edge in ψ, tripped by the repair

**What the receipts and inputs show.**

- **The landscape alone decides [V].** `cross_prior.py` shows the counts depend only on which landscape is used.
  Slots 2962/3010/2986 read 2,949/2,438/2,989 under the current landscape and 341/319/2,241 under the repaired one.
  This holds with either tree's inputs and under either native build.
- **The three slots** are both-strand exons:
  - `Eg` ≈ 9,781;
  - own `log(Eg/Er)` = 0.001, so the repair does not touch their faces;
  - total density 0.30–0.32/bp, with truth about 99 % gDNA, i.e. 0.32/bp.
- **The two landscapes**, log-density relative to the maximum, current vs repaired:

  | density (/bp) | current | repaired |
  |---|---:|---:|
  | 0.003 | −0.40 | −0.16 |
  | 0.03 | −4.03 | −3.20 |
  | 0.3 | −8.59 | −9.92 |
  | 3.2 | −0.01 | −0.00 |

  Mass below 0.01, in 0.01–0.3, and above 0.3: 0.696 / 0.027 / 0.277 in both, identical to three decimals.
- **The replay on the repaired build:**
  - strand only: 3,106 / 2,822 (truth 3,137 / 2,904);
  - strand + landscape: 30 / 30;
  - shipped: 341 / 319.

**What this means.**

- The truth sits in the population's valley. Partially captured exons at 0.3/bp are real ("capture is a spectrum"), but
  the landscape is a kernel sum whose valley depth is a kernel-tail quantity.
- Local evidence (the tilt's Occam factor) favours the truth by roughly 5–6 nats. That figure is inferred, not
  computed: the median held at a 4.6-nat prior penalty and flipped at 6.7.
- Between truth and the 0.03/bp alternative, the prior's penalty grew from 4.6 to 6.7 nats, and the median flipped.
- **[C]** The brief says "These ambiguous exons have no local evidence that outweighs the population". In fact, local
  evidence favours the truth; a deeper valley in the kernel tails outweighs it.

**Consequences.**

- The repair is the trigger. The causes predate it:
  - a population with no mass at intermediate capture (`landscape-has-the-same-premise`);
  - a median readout of a bimodal posterior (`psi-median-discontinuity`).
- Any change that touches kernel tails can swing this row by about 20 points. It cannot be judged on one A/B (`g98-spread`).
- It must not be excused: after the repair the condition is worse than 0.7.1 (69.5 / 42.4 against 56.0 / 39.6).

### `center-is-a-second-mechanism` — the alternative-splice center change is its own mechanism

The old center, `logit f_b − log((U+S)/U)`, scales the odds. Its own map adds the spliced RNA instead:
`U f/(U(1−f) + S)` at equal opportunity. So the old width disagreed with its own map even with no length gap.

**Counterexample:** `f_b = ½`, `U = S = 100`.
- Old predicted odds: 0.5.
- Map odds: 0.333.
- So `d = 0.405`, and the width is about 0.16 nat² at the truth.

The correction ([transfer_kernel.h:337](../../src/rigel/native/transfer_kernel.h#L337)) is right. But it changes widths
on every stranded library whatever the gap, and DESIGN §6b.8's per-pair width was measured with the inconsistent
center. Part of that protection may have been geometric artefact. A/B it alone (`center-arm`).

### `widths-of-the-old-centers` — the widths belong to the old centers, and a face's two directions are not one factor

**(a) The route-sum variance.**
- `Var(log Σ_J J/A_J) = (Σ J/A_J²)/(Σ J/A_J)² ≥ 1/Σ J` (Cauchy–Schwarz), with equality only for equal `A_J`.
- The reverse, alternative-splice and terminus faces are now centred on route sums, but keep `trigamma(S+½)`
  ([transfer_rows.h:204](../../src/rigel/native/transfer_rows.h#L204)).
- Measured `log(A_eff/Er_b)` over flux-bearing boundary sides:
  - median 0;
  - 5th percentile −0.56 (test g50), −0.60 (test g98 ON), −0.95 (RNA-long);
  - minimum −4.8, from reach-limited junctions.
- **Counterexample:** `J = 200` on `A = 200` and `J = 2` on `A = 0.5`.
  - The rate is 5 and the true `Var(log rate)` is 0.32.
  - `trigamma(202.5)` = 0.005, so the sd is understated 8×.
  - If J2 is noise, the face sends a confident wrong message.

**(b) The forward blur has no sensitivity.**
- `∂λ_x/∂log S = −s/(u+s)` and `∂λ_x/∂log U = s/(u+s)`, so
  `Var λ_x ≈ (s/(u+s))²(Var log U + Var log S)`.
- `transport_row` blurs by the bracket at sensitivity 1
  ([transfer_rows.h:190](../../src/rigel/native/transfer_rows.h#L190)).
- At zero flux the map is an exact shift, yet it is blurred by at least π²/2 ≈ 4.93 nat² (`trigamma(½)`).
- The same face's reverse (`splice_out_row` at `splice_rate = 0`) has zero marginal width. The pair is not one factor
  read both ways.
- This predates the repair. Still, the brief's "they must use the same opportunities and central rate" stops one step
  short: they must be the same factor.

**(c) `v_ratio` belongs to the old center.**
- [transfer_kernel.h:340](../../src/rigel/native/transfer_kernel.h#L340): `v_ratio = Var log((U+S)/U)` is the old
  center's variance.
- The new center's variance weights `v_b` by `((1−f) + f·u/(u+rate))²`, and adds
  `(rate/(u+rate))²·Var log rate`.

**Falsifier.** A shallow-count coverage test:
- draw `U, S, J_k` from Poisson around known densities;
- the delivered profile's 68 % interval must cover the recipient's true log-odds about 68 % of the time;
- check per face kind and per direction.

### `what-the-gates-cannot-see` — what the 48 gates and the perturbations cannot see

- **The route sum is untested against `S/Er_b`.**
  - The fixture sets `sj = J·ar[1]·10⁴` and `route = J·10⁴`, so `A_J = Er_b`
    ([test_splice_out_opportunities.py:68](../../tests/calibration/test_splice_out_opportunities.py#L68)).
  - A builder using `S/Er_b` instead of the route sum passes all 48.
  - The archived rate perturbation fires only 2 of 6 cases (`perturb_rate.log`).
  - Add a fixture with `A_J = 0.2·Er_b`.
- **One orientation only.** The fixture is one chain, with the exon on the left, `DON_POS`/`TSS_POS`, + strand, and
  `route_rate_hi = 0`. A builder that always reads `route_lo` passes. Add the mirror: `ACC_POS`, `TES_POS`,
  − strand, `route_hi`.
- **Deep and mode-only.** Counts are ×10⁴, and `_check_map` reads only the argmax, so the width and the marginal are
  never exercised.
- **Perturbations bypass the builder.** They mutate the prepared tables, never the builder's choices:
  - the route side;
  - the `S/Er_b + route` composition;
  - the `s_out`/`rate` pairing;
  - the new `a_r > 0` guards.

  Add builder mutations.
- **The face tests check wiring, not arithmetic.** In `test_transfer_faces.py` the "independent" expectations call
  `R.splice_out_row`, `R.face_map_lambda` and `R.transport_row`, and restate the builder's rules. The numpy-composed
  intron → exon expectation is the exception.
- **Also missing:**
  - the zero-opportunity topology (`zero-opportunity-lane-switch`);
  - shared-fragment fusion at both sides (`shared-crossings-fused-as-independent`);
  - an equal-chemistry identity gate (`two-length-frames`);
  - the coverage test (`widths-of-the-old-centers`).

### `zero-opportunity-lane-switch` — zero opportunity switches the face to the gDNA level lane

At any nonpositive opportunity, `shared_density` declines to install FORWARD
([transfer_kernel.h:195](../../src/rigel/native/transfer_kernel.h#L195)). The gDNA lane then opens a level face
there, because `F.kind == NONE` ([:377](../../src/rigel/native/transfer_kernel.h#L377)). At `a_r = 1e−300` the face
installs as a vacuous FORWARD and blocks the lane instead. So the message is discontinuous at exactly zero (regions
shorter than a law's minimum). The old code always installed FORWARD.

The impact is low, but pin the intended behaviour with a test.

### `off-deltas-are-em-fork-noise` — the capture-OFF test-chromosome losses are not calibration effects

Take g00 ss0.99 OFF:

- transcripts 7.00 → 7.71;
- genes 0.835 in both;
- mRNA −9,060 in both;
- nRNA +90 in both;
- gDNA 2 → 1.

Every other ss0.99 OFF row also moves pools by under 50 fragments. `ISSUES: benchmark-noise-floors-unmeasured`
records about 8k fragments of transcript Σ|Δ| from a last-bit change, roughly 0.7 points at 1.1 M fragments. Under
CLAUDE.md's rule (read genes and pools beside the transcript table), these deltas count neither against the repair
nor, as "beats 0.7.1", for it.

### `shared-crossings-fused-as-independent` — correlated boundary evidence at short regions

Boundary counts are incidences: every crossed boundary gets the full count (`_accumulator_reference.py`). For a region
of length ℓ, a fraction `(w−1−ℓ)₊/(w−1)` of each boundary's crossings cross both ends:

- 60 % for 250 bp gDNA at ℓ = 100;
- 0 for 78 bp RNA.

`solve_block` adds the two held compositions as "independent witnesses", so the shared gDNA evidence is counted
twice — at the short exons of the over-call class. This predates the repair. For the profile design: two boundaries'
`L` are not independent, so never multiply them.

### `capture-identifiability` — what capture information is actually identifiable

**Identified.**
- Per object, a gDNA count identifies one weighted average of capture:
  `ρ0 · Σ_w P0_g(w) Σ_{a∈A_i(w)} C(a)`, where gDNA is separable at that object.
- RNA counts confound capture with expression.
- A per-object length tilt would need per-object gDNA length moments. The `Σ 1/A` banks are one column mixing both
  origins, so the tilt is separable only at gDNA-pure objects.
- At library or class level, gDNA's census law against its uniform law identifies a length-selection shape.

**Scale is not a problem.**
- The EM's weights are relative within a locus, so the `ρ0`/`c` scale cancels.
- No reference, clip or `None` is needed for the EM.

**Count vs mass is benign if weights are intensive.**
- A per-object density is intensive and may be applied per placement in the conserved frame (EQUATIONS §11).
- Summing incidence counts across objects is not valid. It was already refuted: "pooling a region's count with its
  boundaries' crossings".

**The smallest defensible model:**
- a per-object scalar `c_i` from gDNA evidence;
- transferred to RNA under an explicit "uniform within the object" assumption;
- both components' geometry on pre-capture laws;
- junction placements keep the capped rule;
- a library-level length shape only if `transfer-pricing` shows the transfer error is material (Q7).

### `missing-rna-uniform-law` — the law frame is the missing estimator

`E_cap` needs pre-capture laws with selection applied once. The implemented repair (`two-length-frames`) and the RNA lanes already
apply a selected RNA law to capture-blind geometry.

The missing piece is RNA's uniform-frame law. The alternative is a matched-placement pair: within exons the two
census laws agree, 219.00 vs 218.98 bp (`ISSUES: the-scorer-reads-a-census-length-law`).

Two shortcuts are off the table:
- Building RNA's law from gDNA's times a ratio was refused (`ISSUES: a-length-table-built-from-the-other-origins-law`).
- A capture de-tilt of RNA's own law is different in kind, but carries the same tail risk and needs its own falsifier.

### `landscape-has-the-same-premise` — the landscape has the premise the plan falsified for the ruler

**The premise.** `_poisson_kernels` puts a Poisson kernel at ψ's point count
([landscape.py:138](../../src/rigel/calibration/landscape.py#L138)).

**The training weights** ([:223-250](../../src/rigel/calibration/landscape.py#L223-L250);
[calibrate.py:189-228](../../src/rigel/calibration/calibrate.py#L189-L228)):
- `_reliability` gives `w = ref/(v+ref)`, with `ref = 1/max(c,1) + 0.119`.
  - An exon at 100 fragments with `Var(log f) = 0.5` weighs about 0.2.
  - An intron at 1,000 fragments with `Var(log f) = 0.05` weighs about 0.7.
- Non-exon zero-count anchors weigh 1.
- Both-strand regions and boundaries are excluded.

The docstring says the weight "separates classes … every enriched region is an exon". Under capture this pushes the
mixing weights toward the depleted mode by construction.

That is a candidate mechanism for `ISSUES: the-gdna-prior-enters-psi-twice`'s finding that the extra +½ tilt "was
holding enriched exons up" [P]. If it holds, the missing partner is the population estimator, not a ruler.

§10 instead fits its own population inside the ruler, which would give two priors for one quantity. Recommendation:

- one population, fitted from profiles with unit object mass;
- ψ applies it once;
- the ruler reads ψ's posteriors.

Unit object mass is not a free choice: ψ applies the prior per object. An opportunity vote estimates a per-base
population. Measure first (`population-vs-oracle`).

### `profile-derived` — the profile, derived

ψ's Beta(½,½) reference is exactly independent Jeffreys `Gamma(½)` priors on the two component rates (the
Beta–Gamma construction). So, with `M` an object's unspliced count:

```text
L_i(ρ) = Σ_{k=0..M}  Pois(k; ρ·Eg_i) · Γ(M−k+½)/Γ(M−k+1) · S̄_i(k/M)      [× received rows, if message-inclusive]
p_ij  ∝  L_i(ρ_j) · π_j          π_j = the landscape's mass on atom j of its log-ρ grid
```

Here `S̄_i` is the strand likelihood marginalised over ψ's tilt measure. The continuous limit is
`L_i ∝ exp(strand(λ(ρ)))·(1−f)^(−½)`, with `f = ρEg/M`.

**Consequences:**

- **Dividing the posterior by π does not recover `L_i`.** `L_i` is neither `exp(ψ − arm)` nor
  `exp(ψ − arm − Jeffreys)`: it carries the RNA reference and a Jacobian. "Prior-free" means free of π only; the RNA
  reference stays.
- **The limits are derived, not chosen:**
  - `M = 0` → `e^{−ρEg}`;
  - G1 → `Pois(M; ρEg)`;
  - `Eg → 0` → constant in ρ, so no opportunity vote and no `eff ≥ 1` guard are needed;
  - the vertex `k = M` is finite (√π);
  - `ρEg > M` gives the Poisson tail.
- **A short exon should not seed a false mode.** Take `Eg = 0.16` with one wrong-strand read in a hundred at κ = 0.99.
  Across the plausible ρ range `k` moves by at most a fragment or two, and the strand term by about half a nat. So `L`
  is nearly flat in ρ. A point-count kernel is located there instead, and that is how the original false mode formed.
  `profile-false-mode` tests this.
- **The same algebra confirms Part 1.** ψ with the landscape replacing `½ log f` is exactly this posterior on the λ
  lattice. Shipped ψ carries an extra `f^½ ∝ ρ^½`.
- **Storage:** a 260-vector per object, computed in bounded blocks. Use the continuous form away from the vertex.

### `hidden-prior-influence` — prior influence hidden in today's inputs

Each item breaks the plan's own gate ("`L_i` invariant when only π changes").

1. **The strand variance is frozen at the incoming belief.** `Chain::strand_profile` reads `belief[x]`
   ([transfer_kernel.h:124](../../src/rigel/native/transfer_kernel.h#L124)).
   - `var = n·p_ref(1−p_ref)`, which at κ = 0.99 spans 0.0099n to 0.25n.
   - So a claim's width changes up to 25× as the belief moves.
   - Fix: freeze at the own data instead (Q4).
2. **The intron factory is a prior from other objects.**
   - It is `log NegBinom(f·C; ρ_bg·E, α)`, with `ρ_bg` the intergenic background (`density_deconv.py`).
   - It is sent as the intron's "own claim".
   - It is applied in ψ beside the landscape ([solve_kernel.cpp:529,712](../../src/rigel/native/solve_kernel.cpp#L529)),
     and the landscape also trains on intergenic regions.
   - So at an intron the gDNA density prior enters three times: the Jeffreys half, the landscape and the factory.
   - In unstranded libraries it is the only source of composition rows (Q3).
3. **Training reads posterior variances.** Both admission (`located_var`) and the reliability weights use them.

### `unstranded-needs-messages` — unstranded libraries need message-inclusive profiles

**Own evidence alone fails.** In an unstranded library, own evidence comes only from the factory and from G1 objects,
both off-target. A population trained only on own evidence then reads no enrichment under capture. That is the refused
strand-only reference by another route.

**Including messages costs independence.** Message-inclusive profiles put each datum in several objects' terms.
- **Training π:** point estimates stay consistent as a composite likelihood, provided each `L_i` is a correct
  marginal. So this is an efficiency and weighting question.
- **One object's posterior:** fusing correlated messages makes it over-confident (`shared-crossings-fused-as-independent`).

The plan must state which `L` trains π and which `L` is read out.

### `uninformative-reads-the-mean` — an uninformative object reads the population mean

Take the g98 ss0.70 ON landscape (27.7 % of its mass above 0.3/bp) and an object with no gDNA opportunity, such as a
short exon under the RNA-short gap. Its posterior is π itself:

- read by the arithmetic mean, it gets about 300× background;
- read by `exp E log ρ`, about 7× background.

That weight is then applied to RNA placements on the object. Borrowing from neighbours fixes it, at `shared-crossings-fused-as-independent`/F14's price.
The plan must state the readout functional and its uninformative limit.

### `psi-median-discontinuity` — ψ's readout has the discontinuity the ruler's median was rejected for

`posterior_median` ([psi_kernel.h](../../src/rigel/native/psi_kernel.h)) on a bimodal posterior jumps between modes at
a mass tie. `g98-valley-knife-edge` shows that this already moves calibration's counts by thousands, and so the EM's priors. The spectrum
ruling's continuity requirement should cover ψ's readout, not only the ruler's.

## 3. The brief's claims

| Claim | Status |
|---|---|
| The four map families are correct count-frame algebra; the forward and reverse faces share opportunities and center | **V** |
| "All four opportunities are capture-blind geometry, each on its own law" | **C**: the RNA law is the capture-selected census (`two-length-frames`) |
| "No fragment-length estimator was changed" | **V**, but calibration now *depends* on `rna_pmf` at every composition face |
| The 48 tests do not call the production map functions | **V**; but `A_J = Er_b`, one orientation, deep, mode-only (`what-the-gates-cannot-see`) |
| The four families failed 4 / 4 / 18 / 8 cases before their repairs, and the perturbations fire | **V** for the intron, alternative and terminus logs and every `perturb_*.log`; the reverse-splice record is a census and agrees with the brief |
| The g98 loss is wholly a refitted-landscape effect | **V**, but incomplete: identical mass, deeper valley (`g98-valley-knife-edge`) |
| "These ambiguous exons have no local evidence that outweighs the population" | **C**: strand alone recovers the truth (`g98-valley-knife-edge`) |
| Every transcript stratum beats 0.7.1 | **V**; per condition, g98 ss0.70 ON is now worse than 0.7.1 |
| The stranded OFF test-chromosome losses are real | **C** as calibration effects: pools and genes identical (`off-deltas-are-em-fork-noise`) |
| The ladder deferred-row loss comes from these maps | **P**: the capture-signed shift (`two-length-frames`); `law-frame-arm` decides |
| RNA-short ON 13.33 % refutes a ~20 % capture-length floor | **V** as logic; part of the ON movement may be `two-length-frames`, so read `law-frame-arm` on the gap panels |
| Calibration/oracle per-object totals agree within 7.3e-12; replays bit-identical | not re-checked (receipts only) |
| Profiles are the partner for prior-once | **P**: unproven; `landscape-has-the-same-premise`'s estimator bias is the competing explanation (`population-vs-oracle`) |
| `E_cap` with one learned capture serves both components | **P/C**: not identifiable as stated (`capture-identifiability`); needs the RNA uniform law (`missing-rna-uniform-law`) |

## 4. Revised implementation sequence

- **Retain:** the four families' map arithmetic, the 48 gates (extended per `what-the-gates-cannot-see`), the census, replay and `cross_prior`
  instruments, and the oracle-count ruler diagnostics.
- **Remove:**
  - the ruler-first order of §10;
  - any population fitted inside the ruler;
  - `E_cap` as an implementation target (it stays as a specification);
  - the Jeffreys identity as the sole prior gate (`profile-invariance` joins it).
- **Split:**
  - family 3's maps from its center;
  - prior-once from the population estimator, and both from the ruler;
  - the law frame from the maps.

| Stage | Content | Output before advancing |
|---|---|---|
| A (hours, no `src/`) | `law-frame-arm` law frame · `center-arm` center · `g98-spread` g98 · `transcript-noise-floor` noise floor | which repair losses are real, and their cause |
| B | land the maps per Q1/Q2, with the missing falsifiers (`what-the-gates-cannot-see`) and an equal-chemistry identity gate | the suite re-derived; calibration-vs-oracle per class and stratum |
| C | widths: one pair factor read both ways; carry `Σ J/A²` per side for the route-sum variance | its own A/B on the calibration metric, including the coverage gate |
| D | the population: `population-vs-oracle` → exact `L_i` offline on frozen dumps (`profile-false-mode`; `profile-invariance` after Q4) → a unit-mass landscape from profiles, replacing the kernels and training weights → prior-once on top → the factory per Q3 | each one alone, on calibration vs oracle, with zero controls separate |
| E | ψ readout continuity (`psi-median-discontinuity`), gated with a tie sweep; the fix lands with D's population | g98 stable under `g98-spread`'s perturbations |
| F | ruler = a continuous functional of D's posteriors, with spatial borrowing for no-opportunity objects; delete reference/clip/`None` | owner checkpoint |
| G | capture transfer: `transfer-pricing` first; a reduced model only if material (Q7) | the priced gap per probe layout |

## 5. The cheapest decisive experiments

| | Uncertainty | Experiment | Cost | Rejects the hypothesis if… |
|---|---|---|---|---|
| `law-frame-arm` | `two-length-frames` causes the capture-signed shift | Re-calibrate cached scans with the repaired build. Override only `calibrate()`'s `rna_fl_pmf` argument, not the EM's `FLModels`: `rna_fl_pmf := gdna_fl_pmf` on the ladder and the test chromosome (equal chemistry); the simulator's true laws on the gap panels. Census and calibration vs oracle; then end-to-end only on the ladder's capture-ON rows | seconds per calibration; about 10 min end-to-end, sharded | the ON shift persists with commensurate laws |
| `center-arm` | the center change is harmless | Repaired maps with the old center vs the new, on the test chromosome's stranded rows and ladder stranded ON; `policy_benchmark.py --by-class` at the alt-ss boundary class | minutes | alt-ss boundary error rises with the new center (then `center-is-a-second-mechanism` unmasks a premise bias, the "other half") |
| `g98-spread` | g98 is a knife-edge | (a) Diff the two landscapes' training sets (slots in/out, kernel centres, weights) from frozen dumps. (b) On the current tree, K small unrelated perturbations (drop one random located slot; jitter one kernel width); count flips | seconds | flips do not occur under perturbations of the repair's size (then the repair moved something material: find it) |
| `transcript-noise-floor` | OFF transcript deltas are noise | A last-bit perturbation of calibration output on all 30 test conditions, on the current tree | minutes | the g00 ss0.99 OFF spread is ≪ 0.7 points |
| `population-vs-oracle` | the prior-once partner is the population estimator | On stranded and deferred ON rows, compare enriched-mode mass and exon-class posterior mass: the shipped landscape vs a unit-mass NPMLE on the same kernels vs the oracle population (`slot_truth.npz`) | seconds | the shipped enriched mass already matches the oracle (then exchangeability is the problem → Q5) |
| `profile-false-mode` | exact profiles prevent the false mode | Offline on the archived pre-repair RNA-short stranded OFF dump: `L_i` (`profile-derived`) for single-strand regions, NPMLE | seconds | an enriched atom forms from the short exons |
| `transfer-pricing` | scalar capture transfer is good enough | The B_match vs B_transfer1 arms already specified in `ISSUES: the-scorer-reads-a-census-length-law`, on the test-chromosome gap panels and the flgap ON rows | minutes | the B_match − B_transfer1 gap exceeds the release margin in any in-scope stratum (then a component-specific opportunity is needed) |
| `profile-invariance` | `L_i` is prior-free | Frozen messages; swap the two archived g98 landscapes; extract `L_i`; compare | seconds | `L_i` moves (expected today, via `hidden-prior-influence`(1); a passing `profile-invariance` is the gate for D) |

## 6. Decisions for the owner

1. **Q1 — land the maps before the law frame is fixed?** Only if `law-frame-arm` confirms `two-length-frames`. Options:
   - land now and record the capture-ON shift as release-blocking under the length-frame entry; or
   - hold until RNA's law is in gDNA's frame.

   Landing now keeps RNA-short stranded OFF at 42.3 → 3.53 %. It costs a systematic gDNA → nRNA shift on every
   captured library: stranded ON gDNA −1.6k to −6k, transcripts flat; deferred +2 points.
2. **Q2 — which frame do composition maps read?**
   - Both uniform laws: needs a capture de-tilt of RNA's spliced law, with tail risk (`missing-rna-uniform-law`).
   - A placement-matched pair: needs gDNA's length-resolved mass on RNA's placements; stranded-only, via within-object
     strand contrast.
3. **Q3 — is the intron factory a prior or data?**
   - As a class prior for introns, it should replace π there. That is a class-conditional prior already shipped,
     which needs reconciling with "one population".
   - As borrowed data, intergenic evidence must not also train π, or it enters twice.
4. **Q4 — move the strand variance freeze from the incoming belief to an own-data reference?** This is required for a
   prior-free profile. It revisits the count-zero-information freeze ruling.
5. **Q5 — if `population-vs-oracle` shows exchangeability, not estimator bias, is the problem:** is a class-dependent mixing weight
   admissible? The component shapes stay shared, the weights are learned and continuous, and the split is inert off
   capture by construction. Or is that the refused class split?
6. **Q6 — ψ's count readout: keep the median, or a continuous functional?** The median is robust at vertices and
   discontinuous at bimodality. The mean of the count is the alternative; or hand the EM the mixture.
7. **Q7 — if `transfer-pricing` shows the scalar transfer is material:** is a library-level length-selection *shape* admissible under
   the spectrum ruling? Intensities stay per object, with no on/off.

Settled by derivation, so not a question:
- the population's unit is the object (ψ applies the prior per object);
- the RNA reference measure is part of `L_i`;
- the EM needs no capture reference (`capture-identifiability`).

## Appendix — reproducing the two offline computations

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
cd ~/Downloads/rigel_runs/prototypes/2026-10-06_rna_short_counts
python - <<'EOF'
import pickle, numpy as np
for p in ["test_g50_terminus/scenarios_0.50_native.pkl", "test_g98_terminus/scenarios_0.70_on_native.pkl"]:
    d = pickle.load(open(p, "rb")); ag, ar, b = d["eff_gdna"], d["eff_rna"], d["is_boundary"]
    ok = (ag > 0) & (ar > 0); x = np.log(ag[~b & ok] / ar[~b & ok])
    print(p, np.unique(np.round(np.log(ag[b & ok] / ar[b & ok]), 4)), np.quantile(x, [.5, .95]), x.max())
    for side in ("lo", "hi"):                          # junction opportunity against Er_b
        sj, rr = d[f"sj_count_{side}"].sum(1), d[f"route_rate_{side}"].sum(1); m = b & (sj > 0) & (rr > 0)
        print(side, np.quantile(np.log(sj[m] / rr[m] / ar[m]), [.05, .5]))
D = {n: pickle.load(open(f"test_g98_{n}/scenarios_0.70_on_native.pkl", "rb")) for n in ("current", "terminus")}
for n, d in D.items():                                 # the landscape at the AMBIG exons' competing densities
    lr, lp = d["gdna"]; lp = lp - lp.max()
    print(n, [round(float(np.interp(np.log(r), lr, lp)), 2) for r in (0.003, 0.03, 0.3, 3.2)])
EOF
```

The pool deltas in `two-length-frames`/F7 are read from `current_{ladder,test}.jsonl` against `terminus_{ladder,test}.jsonl`
(`axis == "library"`, `*_est − *_true`).
