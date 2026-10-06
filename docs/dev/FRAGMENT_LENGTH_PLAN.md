# The fragment-length fix: the plan

*Sandbox document (`docs/dev/`): provisional, not authoritative, and cited by nothing outside the sandbox. Fifth
version, 2026-10-02, against `bdfd8709`. Nothing in it is built. The problem statement is
[`FRAGMENT_LENGTH_REVIEW.md`](FRAGMENT_LENGTH_REVIEW.md); the reviews that shaped this plan are in
[`FRAGMENT_LENGTH_REVIEWS.md`](FRAGMENT_LENGTH_REVIEWS.md); earlier versions are kept outside the repo.*

## Where it stands

**Status, 2026-10-03:**
- Stage 0 tier 1, Stage 1, Stage 3 and the bench-only Part 2 prototype have run. The verdicts and the owner's open
  questions are in `~/Downloads/rigel_runs/prototypes/HANDOFF_2026-10-03.md`.
- In short:
  - Part 1 is correct but needs a class-conditional partner.
  - Part 2 is right as the ruler's reference source and not ready as ψ's prior.
  - On Part 2's nuisance, a shared profiled RNA amount collapses at κ 0.7. Conditioning on the region total is
    consistent there, but reads RNA-free one-strand regions low.
- The flgap and junction panels are being re-simulated on the current physics.

**The principle.** Calibration learns from observations with their likelihoods, never from an estimate recycled as
data, and never from a prior counted twice.

**Ready now:**
- **Stage 0:** the falsifier bench.
- **Stage 1:** the prior arithmetic.
- **Stage 3:** the capture × length experiment as ISSUES specifies it.
- **A bench-only prototype of Part 2.**

The fourth and fifth reviewers agree on Stages 0, 1 and 3; the fifth also on the bench-only prototype. The prior
arithmetic is a correction; the training redesign remains a candidate.

**Not yet:**
- **Part 2 in `src/`:** how the capture reference is read, and whether evidence-only training is adequate for
  unstranded capture-OFF, are open (§6).
- **Part 3:** it has an identifiability gap and waits on Stage 3.

## 1. The design

### Part 1: the prior enters ψ once

**The defect.** The gDNA rate prior has one place in ψ. It holds the reference's gDNA half, `½·log f`, when there is no
fitted prior, and the fitted landscape when there is. Today the kernel adds both (`psi_kernel.h:182`), which tilts
gDNA by `+½` per nat of `log ρ`. The derivation is in the reviews (second review, item 1).

**The changes.**
- Split the reference into its gDNA and RNA halves; use the landscape *instead of* the gDNA half.
- Evaluate the landscape at the unclipped, stably computed `log σ(λ)`, not at the `10⁻¹²` clamp.
- Below the landscape's grid, continue from its end value at the reference's slope of `½` per nat.
- Build the grid's domain from the same floored opportunities and masses the arm evaluates, for every slot the arm
  is read at. Today the domain stops at `10⁻⁹` while the arm floors at `10⁻¹²`, so a counted slot in between reads
  above the grid. Equal constants are not enough: a domain filter could still drop a slot the arm processes. Its
  gate enumerates every consumer's actual coordinates and checks that none exceeds the grid's top.
- Make the bracket demand (`required_logodds_window`) cover the completed lower tail, not only the landscape's floor.
- Leave the arm's `T = N` plug-in alone here. It is a separate approximation, and training changes do not settle it.
  A quadrature check on sparse counts and on narrow or bimodal landscapes decides whether it changes.

### Part 2: the prior is fitted from evidence

Each region contributes its own likelihood `L_o(z, Y_o)` over `z = log ρ_g`, where `Y_o` is the region's RNA amount.
The population fit maximises

```
Σ_o log Σ_z q(z) · L_o(z, Y_o)      over q and one Y_o per region
```

by plain EM, with the objective checked at every step.
- **No fragment is labelled.** For each candidate density, `L_o` scores how well it explains the region's counts.
- **Each region has one RNA amount, shared by every candidate density.** It is maximised inside the fit, neither
  integrated nor re-maximised per candidate. The two rejected forms each manufacture gDNA:
  - **Integrating `Y` against its reference** leaks that reference into the evidence. A count-only region then
    "prefers" gDNA near its whole count, by +0.6 nats at `N = 10` and +1.2 at `N = 100`.
  - **Re-maximising `Y` per candidate** lets every candidate fit its own RNA amount. In a pure-RNA null (expected RNA
    10, `κ` 0.99, no gDNA), a false 0.1 % weight at one gDNA fragment then gains 2.1×10⁻⁵ nat per region by exact
    enumeration, and the gain grows with the number of regions.
  - **One shared amount** loses 5.1×10⁻⁵ in the same null. Neither surviving form is declared adequate before the
    bench's non-null cases.
- **Each `Y_o` update is a one-dimensional concave maximisation.** The likelihoods are products of Poissons, and there
  is no quadrature.

**The first prototype takes two region classes, chosen by structure, never by outcome.** Write `X = e^z·E_g`.

| region | `L_o(z, Y)` |
|---|---|
| admits no RNA | `Poisson(n | X)`, with no `Y` |
| admits one RNA strand, stranded library | `Poisson(u_anti | ½X + (1−κ)Y) · Poisson(u_sense | ½X + κY)` |

A zero count gives `e^{−X}` in both rows.

**Added in step two**, with the same shared-nuisance treatment:
- **unstranded regions:** `Poisson(N | X + Y)`;
- **both-stranded regions:** one RNA amount per strand, so the minority strand bounds gDNA, as a 90:10 split should.

**Excluded:** boundaries, because a crossing count books one fragment once per boundary; and messages, the factory row
and any earlier landscape, none of which is the region's own evidence.

**Three rules the fit needs:**
- **Zero support.** Zero counts pooled over gDNA-only regions resolve densities near `1/ΣE`, far below today's grid
  floor of `1/max E`. That is a resolution scale, not a floor: their likelihood is still `e^{−1}` there and is
  largest at zero. So the lower support is extended until predictions and the consumed quantities stop changing.
  The fitted mass at the floor is not the test, since it can stay at 100 % while the implied rate keeps moving. Zero
  states stay closed for 0.8.0, so there is no zero atom.
- **An unidentified fit is reported, never published.** If the likelihoods carry no information, EM returns its
  starting point; the fit then says so and installs no landscape.
- **Training no longer needs a sweep.** Every input is in the scan payload: counts, strand counts, opportunities and
  κ. So `q̂` is fitted once, and one sweep runs with it. Pass-0's training role and the three refits go: with a fixed
  prior every refit repeats the same calculation.

### Part 3: the length laws (waits on Stage 3)

**The idea.** At a shared placement and length, capture is common to both origins. Within one-RNA-strand objects,
`logit q_ob = α_o − log r_b`, with `r_b` the RNA:gDNA length ratio.

**What is open:**
- `r_b` is identified only up to a scale factor. An independently established RNA law fixes that scale on shared
  support, but a capture-biased spliced histogram cannot supply it.
- The capture response without a probe BED.

**When it is prototyped:** use exact integer lengths stored sparsely, with no bins until their error is measured, and
keep conserved mass jointly with length and strand.

## 2. Stage 0: the falsifier bench (2 days)

Synthetic regions with known truth, counts drawn from the exact model, run through each candidate. It has three tiers,
in this order.
- **Tier 1, analytic unit tests:**
  - **the population null, first:** a pure-RNA population must not reward weight at a false gDNA alternative. Run it
    for both surviving nuisance forms, and also on non-null cases;
  - `N` of 0 and 1;
  - `κ = ½`;
  - both RNA strands;
  - zero opportunity;
  - zero RNA;
  - all-flat input, which must report "unidentified";
  - the count-only region, which must be flat below its total;
  - pooled zero counts, whose predictions must converge as the lower support extends.
- **Tier 2, population falsifiers:**
  - **RNA depth at fixed tiny gDNA opportunity** (30, 300 and 3,000 RNA reads; `E_g` of 0.1 to 10 positions). The
    false mode must not form. It fails on the shipped tree.
  - **An enriched minority,** with its size, exposure and RNA depth varied, including short probed exons. It must be
    recovered wherever the observations can distinguish it, and not where they cannot.
  - **`g00`.**
  - **A mis-specified strand error** (1.5× the model's).
- **Tier 3, at the panel stage:** repeated realisations, held-out loci and oracle-prior controls, run on the
  genome-scale cache rather than the bench.

Arms: the shipped training, the guard (a comparison arm only), and Part 2's prototype. What survives becomes permanent
tests in `tests/calibration/`, verified failing first and broken on purpose afterwards.

## 3. Stage 1: Part 1 (2–3 days)

- **Where:** a C++ worktree.
- **The gates:**
  1. A Jeffreys-shaped landscape reproduces the reference-only ψ.
  2. The consumed quantities converge as the omitted tail mass goes to zero. This replaces equality at fixed brackets:
     even the reference alone reads `Var(log f)` 2.804 at `L = 10` against 3.290 at `L ≥ 40`, because 0.9 % of its
     mass lies outside the shipped bracket.
  3. Every arm consumer's actual coordinates lie inside the grid's represented domain.
  4. A slot below `log 10⁻¹²` still sees the tail slope.
- **The A/B:** the test chromosome, then the ladder. Read per stratum: the zero controls and unstranded rows, the
  opportunity bins, and a pinned `quant_accuracy.py` pair.

## 4. Stage 2: Part 2 (prototype about 1 week; panel A/B about 1 week)

1. A Python prototype outside the tree, built on the cached payload and judged on the bench.
2. **Two measurements, kept apart:**
   - one sweep against repeated, fully reset sweeps under the same fixed prior, which should be identical;
   - the new prior's accuracy per stratum, with **unstranded capture-OFF read explicitly** (it loses message-derived
     training) and the deferred stratum reported.
3. The genome-scale cached loop, read with `calibration_vs_oracle.py`, `prior_vs_oracle.py` and the opportunity table.
4. The pinned end-to-end A/B, then the candidate broken on purpose.

## 5. Stage 3: capture × length (2–3 days, in parallel)

- Run `ISSUES: the-scorer-reads-a-census-length-law`'s arms `B_match` and `B_transfer1` against that entry's controls,
  unchanged.
- A partial improvement establishes a contributing mechanism. Failing to close the whole gap does not prove the
  mechanism lies elsewhere.
- Part 3 is designed only after this.

## 6. Open questions

1. **For you:** start Stages 0, 1 and 3, and the bench prototype of Part 2?
2. **For reviewers:** is one shared RNA amount per region the right nuisance treatment? It passes the pure-RNA null,
   where re-maximising per candidate fails it. Tier 1's non-null cases are its next test.
3. **Part 2's capture reference.** How should the capture level be read off `q̂`? Many weak regions can locate a
   population together, so the shipped rule that counts individually located members needs replacing or a fresh
   justification. The prototype measures candidates first.
4. **Unstranded capture-OFF.** Is the evidence-only prior adequate there? If not, a junction's flux enters its exon's
   likelihood, counted once.

## 7. Rejecting a candidate

Reject a candidate that does any of these:
- learns composition from blind regions' reference rather than their data;
- creates locations through numerical boundaries;
- erases an enriched minority the observations can distinguish;
- moves the capture reference away from truth;
- increases in-scope truth error;
- succeeds only by disabling genuine capture.

Capture-ON rows are judged in scope; the deferred stratum is reported. The guard is a comparison arm, not a production
fallback.

## 8. What no design can do

- **A region that cannot hold a gDNA fragment has no gDNA density to measure.**
- **An unstranded both-stranded exon has no composition channel of its own.**
- **Without probe information, capture and abundance can be confounded.**

A coherent model removes artificial information; it cannot create real information. These limits are reported as
uncertainty.

## Symbols

| symbol | meaning |
|---|---|
| `N`, `n` | a region's unspliced count; `n` where the region admits no RNA |
| `u_anti`, `u_sense` | a one-RNA-strand region's counts against and along its RNA strand |
| `E_g` | gDNA's opportunity: contained positions for a region |
| `f_g`, `λ`, `z` | the gDNA share, its logit (the solver's axis), and `log ρ_g` (the landscape's axis) |
| `X`, `Y` | a region's expected gDNA and RNA counts |
| `κ` | the library's RNA sense fraction |
| `q`, `q̂` | the gDNA population prior, and its fit |
| `L_o`, `Y_o` | region `o`'s likelihood over `z`, and its one RNA amount |
| `r_b`, `α_o`, `q_ob` | the RNA:gDNA length ratio at length `b`; object `o`'s log density ratio; its gDNA fraction at `b` |
