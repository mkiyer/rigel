# The gDNA fragment-length model — what it actually costs, and the estimator that survives (2026-08-30)

    ⚠ **A DEV DOC and a SANDBOX.** Nothing here is authoritative and nothing may cite into it. The
    diagnosis it builds on is `DIAGNOSIS_gdna_fragment_length.md` and still stands unchanged.

⭐⭐⭐ **THE HEADLINE, AND IT RE-RANKS THE WORK: this is NOT a calibration defect.** Measured on the
fl-gap arms — the only substrate that can see it — a **perfect** gDNA length pmf moves the calibration
result by less than the benchmark resolves, in an **inconsistent direction between the two arms**. The
exposure that is real is the **EM's per-fragment length likelihood**, which no instrument in the repo
priced until this session promoted `scripts/design/em_fl_ceiling.py`; measured there it is worth **8–32 %
of the library gDNA bias**, the product. ⛔ So rank 1b should not be ranked on the 0.8.0 metric, and the
previous session's sequencing —
"regenerate the arms, price with `length_ceiling.py`, then build Option A" — is answered rather than
pending.

## 0. The gate is CLOSED: the fl-gap arms already exist on the current model

⛔ Rank 1b said the price was "UNRANKABLE until the fl-gap arms are regenerated on the current
sparse-nascent model". That is true of the **ladder-scale** side panels in `suite/`. It is **not** true of
the test chromosome, where three fl panels already exist, fully simulated and oracle-cached 30/30:

    scenarios_fl_rna_long    RNA 250 / gDNA 100   (the +sign arm)
    scenarios_fl_gdna_long   gDNA 250 / RNA 100   (the −sign arm)
    scenarios_fl_equal200    both 200             (the control)

⭐ **They carry the CURRENT nascent model already**, and that is structural rather than lucky: the test
reference is `abundance.mode: file`, so the `nrna:` block in each config is dead and nascent comes from
`test_abundances.tsv` — the same file the main test panel reads. Confirmed by measurement, not by reading
the config: the condition names carry `nrna_file`, pool 1 holds 34,291 nascent fragments, and pool 0 holds
8,423 **mature** fragments, which can only be the shadow transcripts. ⛔ Re-simulating anything was
therefore unnecessary and nothing was re-simulated.

## 1. What a PERFECT gDNA length model is worth to CALIBRATION — measured, both arms

`length_ceiling.py --suite <fl arm>`, `mwae_all` on the shipped column, capture-OFF:

| | region, both pmfs exact | boundary, both pmfs exact |
|---|---|---|
| **rna_long** g05/g50 × ss0.99/ss0.50 | **+0.8 %, −0.3 %, +5.6 %, +2.7 %** | −35.5 %, −5.9 %, −35.4 %, −10.8 % |
| **gdna_long** same four | **−1.2 %, −1.6 %, −3.1 %, −4.7 %** | +15.5 %, −13.7 %, **+32.4 %**, +2.5 % |

⛔ **The sign is inconsistent BETWEEN THE ARMS on both axes.** Making the gDNA pmf exact is a small win on
one arm and a small loss on the other, at regions and at boundaries alike; the deliverable (`mwae_all_final`,
gDNA pmf exact alone) moves +0.10 … +0.30 % on rna_long and −1.95 … +0.90 % on gdna_long. Every magnitude
is at or under the benchmark's own stated resolution (`NEXT_SESSION.md`: do not believe a single row below
~2 %). ⭐ The shape is the one already recorded for the effective-length ruler in `ROADMAP.md` §0 — **both
rulers are wrong and their errors partly cancel**, so correcting one alone is as likely to hurt as help.

**Why calibration is nearly blind, derived then measured** (a session-scoped prototype; the algebra is the two lines above). The
pmf reaches calibration only through opportunity functionals, and the two have completely different
sensitivity:

    REGION   E[(ell−L+1)+]   error enters DIVIDED BY THE OBJECT LENGTH
    BOUNDARY E[L−1]          error passes through UNDIVIDED

Measured at a shipped-vs-true gap of +121.01 bp: the boundary functional is **+106.94 %** wrong, while the
region functional is −65.9 % at ℓ=100 bp, −1.22 % at 10 kb and **−0.06 % at 200 kb**. gDNA is measured at
introns and intergenic regions, which are kb- to Mb-scale, so the error is divided away exactly where the
gDNA evidence lives. ⭐ This is `TRAPS: a-pooled-rate-cannot-see-a-short-object-factor` in a second place.

## 2. Where the exposure IS: the EM — and it took a new instrument to see it

`pipeline.py:816` turns `fl_models.gdna_pmf` into a `FragmentLengthModel` and hands it to `_score_fragments`.
For an unspliced fragment the length term contributes `log g(L) − log r(L)` to the origin call, so a wrong
`g` injects `log g_ship(L) − log g_true(L)` nats **per fragment**, independent of `r`
(a session-scoped prototype; the statistic is the line above, weighted by the gDNA fragments that actually occur):

| panel | condition | median \|error\| | 95th | gDNA mass whose length-term origin preference FLIPS |
|---|---|---|---|---|
| rna_long | g05 ss0.99 off | **1.005 nats** | 1.250 | 7.7 % |
| gdna_long | g05 ss0.99 off | **1.686 nats** | 1.768 | 6.0 % |
| rna_long | g50 ss0.99 off | 0.165 | 0.188 | 0.3 % |
| equal-200 | g05 ss0.99 off | 0.141 | 0.484 | 23.4 % |

⭐ A median of 1.0–1.7 nats is a factor of 2.7–5.4 on the per-fragment origin likelihood ratio. ⚠ **Read the
equal-length row's 23.4 % as noise, not as damage**: when `g` and `r` are the same shape the sign of
`log(g/r)` is arbitrary, so the flip rate is high while the magnitude is negligible. Magnitude is the
column that matters.

⛔ **No instrument measured this** — `length_ceiling.py` and `calibration_truth_ab.py --ceiling` both stop at
`calibrate`, and `quant_accuracy.py`'s arms all inject `LocusPriors` fields, never an fl pmf — so one was
built and PROMOTED: `scripts/design/em_fl_ceiling.py`, which reuses `quant_accuracy.run_condition` rather than
re-implementing a scorer, and carries the harness's own three gates (the injection counts its fires; a
`noop_fl` arm that replaces the pmf **with itself** must come back byte-identical; `base_reseed` is the
noise floor). ⭐ The noop gate FIRED on its first run — renormalising the pmf moved it by 8.7e-19 — which is
the perturbation half of `TRAPS: perturb-every-gate` catching the harness rather than the code.

### ⭐⭐⭐ THE ANSWER: it is worth 8–32 % of the LIBRARY gDNA BIAS, and nothing else

`gdna_frac_est` is `cli.py`'s own `gdna_fraction` — **the product**. Capture-OFF, true gDNA fraction from
`origin_counts` (exactly 0.050000 / 0.500000):

| arm | condition | base | perfect gDNA pmf | truth | bias removed | vs reseed floor |
|---|---|---|---|---|---|---|
| rna_long | g05 ss0.99 | 0.099320 | 0.095205 | 0.05 | **8.3 %** | 22× |
| rna_long | g50 ss0.99 | 0.538155 | 0.529610 | 0.50 | **22.4 %** | 47× |
| gdna_long | g05 ss0.99 | 0.101440 | 0.095065 | 0.05 | **12.4 %** | 29× |
| gdna_long | g50 ss0.99 | 0.535215 | 0.528835 | 0.50 | **18.1 %** | 9× |
| gdna_long | **g50 ss0.50 (UNSTRANDED)** | 0.548280 | 0.533055 | 0.50 | **31.5 %** | 27× |
| **equal-200 control** | g05 / g50 ss0.99 | 0.100200 / 0.551590 | 0.099455 / 0.551460 | — | (−0.0007 / **inside the floor**) | 17× / **0.4×** |

⭐⭐ **All five fl-gap rows move the deliverable TOWARD truth, on BOTH sign arms, far outside the noise
floor** — the consistency the calibration-level ceiling of §1 did not have. The reason it is one-signed is
that the tool **over-calls gDNA at every rung here**, and a correct gDNA length model reduces the over-call
whichever way the length gap points. ⭐ The displaced mass goes to nascent (`nrna_est` +857 / +1714), not to
mature. ⭐ The equal-length control is inert, which is what makes the rest a measurement rather than an
artefact.

⛔ **TRANSCRIPT accuracy is NOT a win and must not be quoted as one.** `count_abs_err` improves on rna_long
(−36 of 13,373; −31 of 7,432) and **regresses** on gdna_long (+48, +107, +146 of ~8–14 k), all outside their
floors. Sign-inconsistent between the arms, exactly like §1. ⚠ And the fl-gap panel's transcript number is
not interpretable anyway (`TRAPS: a-length-gap-bypasses-calibration`) — the LIBRARY row is, because it is a
composition deliverable rather than an EM assignment.

## 3. THE STORAGE QUESTION — answered with arithmetic, and the answer is no

The owner asked whether to store FL histograms in every region and boundary. On the panel index
(chr21+chr22) there are **70,270 regions + 70,082 boundaries = 140,352 objects** at `max_length` 1000:

| | chr21+22 | whole genome (×32) |
|---|---|---|
| per-object FL histogram, float32 | 0.56 GB | **~18 GB** |
| per-object FL histogram, uint16 | 0.28 GB | ~9 GB |
| two per-object length MOMENT banks (ΣL, ΣL²), float64 | 2.2 MB | ~0.1 GB |
| **the current `pool_lengths` bank, 5 × 1001** | **39 KB** | **39 KB — constant in genome size** |

⭐⭐ **THE ASYMMETRY IS THE WHOLE ANSWER: the POOL axis is library-wide and costs nothing to refine; the
OBJECT axis is genome-sized.** Splitting `pool_lengths` by strand takes 39 KB → 78 KB. Storing per-object
histograms takes 39 KB → 18 GB, to serve a composition channel that is **retired until after 0.8.0** and was
refuted on its own terms (`TRAPS: a-linear-likelihood-emits-a-sign` — a MODE failure, so a better pmf does
not rescue it). ⛔ **Do not store per-object FL histograms.** If per-object length information is ever
wanted, the affordable form is moments, and the model-free local mean is already free: `count /
inv_length_sum + 1 = E[w]` (`EQUATIONS.md` §3c), from two banks that already exist.

## 4. THE ESTIMATOR THAT SURVIVES — a two-pool contrast, no template, no schema change

**The derivation.** Both CONTAINED pools share one opportunity geometry, and their contaminants — nascent
RNA in an intron, unannotated transcription in intergenic space — are genomically **contiguous** there
exactly as gDNA is, so the length tilt cancels identically for both components. De-tilted, each pool is a
mixture of the same two shapes at **different weights**:

    f_p(L) = a_p·g(L) + (1 − a_p)·r(L)        p ∈ {0 intergenic, 1 intronic}

Two equations, two unknown functions, and `g` falls out with **no template at all**:

    g = [ (1 − a_1)·f_0  −  (1 − a_0)·f_1 ] / (a_0 − a_1)

⭐⭐ **The conditioning is `1/(a_0 − a_1)` — the SEPARATION of the two purities — not the `1/a` of template
subtraction.** That is precisely what killed the refuted option: at `a = 0.18` a 5 % template error left a
−14 bp residual. Here the separation is measured at **0.42–0.61** off capture, so the amplification is ~2×
rather than ~21×. ⭐ And it is **exactly the identity when the pools agree**: set `f_0 = f_1 = f` and the
numerator is `(a_0 − a_1)·f`, so a pool pair with nothing to say changes nothing.

**Measured against the origin-split oracle** (a session-scoped prototype of the formula above, oracle weights, capture-OFF,
mean L error in bp — POOLED is what ships today):

| arm | condition | POOLED | pool 0 alone | **CONTRAST** | reduction |
|---|---|---|---|---|---|
| rna_long | g05 ss0.99 | +118.96 | +71.37 | **+10.16** | 12× |
| rna_long | g25 ss0.99 | +59.91 | +20.20 | **+2.39** | 25× |
| rna_long | g50 ss0.99 | +27.96 | +7.59 | **+0.89** | 31× |
| gdna_long | g05 ss0.99 | −120.59 | −80.10 | **−4.70** | 26× |
| gdna_long | g25 ss0.99 | −60.84 | −22.33 | **+0.17** | 360× |
| gdna_long | g50 ss0.99 | −27.81 | −8.37 | **−0.01** | 2800× |
| gdna_long | **g50 ss0.50 (UNSTRANDED)** | −27.66 | −8.27 | **−0.26** | 106× |

⭐⭐⭐ **The unstranded row is the reason this beats the strand split.** `Option A` (splitting the length
banks by strand) is undefined at κ = ½ and has no orientation for the intergenic pool at all — and
unstranded × capture-OFF is an IN-SCOPE 0.8.0 stratum. The contrast needs no strand, no schema change and
no new bank.

⭐ **It is self-limiting, and that was measured rather than designed in.** `a_0 − a_1` collapses exactly
where nothing needs correcting: **0.0000 at g00** (the estimator refuses — with no gDNA there is no gDNA
length distribution), 0.0197 at g98, 0.0336 at g50 capture-ON.

## 5. ⭐⭐⭐ THE MIXING WEIGHTS — the away-half again, this time on DENSITY. BUILT AND MEASURED.

The contrast needs `a_0`, `a_1`. By definition `a_p = n_gdna,p / n_p`, and gDNA being genomically UNIFORM,
`n_gdna,p = rho_g · E_p`. So the two weights need exactly ONE scalar, the gDNA density `rho_g`
— and note `rho_g` **cancels from the ratio** `a_0/a_1 = (E_0/n_0)·(n_1/E_1)`, which is index geometry
plus observed counts, so only the LEVEL is unknown.

⭐⭐ **THE ONE-SIDED ARGUMENT: transcription can only ADD fragments to an object, never remove them.** So
every object's observed density `n_i/E_i` is an OVERESTIMATE of `rho_g`, and the clean objects are the
ones on the LOW side. That is exactly `TRAPS: purity-is-a-property-of-the-annotation`'s prescription —
*where a contaminant can only push a statistic ONE WAY, use the half it cannot reach* — transposed from
the strand residual to density, and it asserts no class pure.

    rho_hat(q) = sum_{lowest q by density} n_i  /  sum_{lowest q} E_i          (q = 1 is TODAY's rate)

**Measured** (test chromosome, intergenic + intronic REGIONs, `rho_hat/rho_true`): today's pooled rate is
**8.46× too high** at `g05`; the curve has a clear PLATEAU at 0.890 / 0.960 / 1.023 for q = 0.25/0.50/0.75,
then breaks upward (1.77× at 0.90, 3.34× at 0.95). ⭐ At `g05` the resulting weights are **(1.000, 0.049)
against an oracle (1.000, 0.048)**, where today's rate gives `a_1 = 0.403` — an 8× error.

### ⭐⭐⭐ THE WHOLE CHAIN, WITH NO ORACLE INPUT ANYWHERE

`(1)` per-object density → `(2)` `rho_g` from its low side → `(3)` `a_p = rho_g·E_p/n_p` → `(4)` the
contrast. The oracle is read only to SCORE. Mean gDNA fragment length, capture-OFF, error in bp:

| condition | POOLED (ships today) | **CHAIN** | oracle-weight ceiling |
|---|---|---|---|
| rna_long `g05` | +118.96 | **+8.88** | +10.16 |
| rna_long `g25` | +59.91 | **+0.77** | +2.39 |
| rna_long `g50` | +27.96 | **−0.02** | +0.89 |
| rna_long `g50` **unstranded** | +27.49 | **−0.04** | +0.82 |
| gdna_long `g05` | −120.59 | **+1.17** | −4.70 |
| gdna_long `g25` | −60.84 | **+5.27** | +0.17 |
| gdna_long `g50` | −27.81 | **+2.34** | −0.01 |
| gdna_long `g50` **unstranded** | −27.66 | **+3.02** | −0.26 |

⭐ **89–99.9 % of the error removed, on BOTH sign arms, stranded and unstranded, with no schema change.**
⚠ It sometimes beats the oracle-weight ceiling because `rho_hat` runs a few percent low and that partly
cancels the shared-`r` premise error — **luck, not design**; do not report it as headroom.

**The controls behave.** `g00`: nothing to estimate, and it says so. `g98`: `a_0 − a_1` collapses to 0.008–0.03
and the estimator **REFUSES** — correctly, since the pooled estimate is already within 0.7 bp there.
⭐ **On the equal-length LADDER at realistic scale (70,270 regions) it is inert**: `rho_hat/rho_true`
= **0.999–1.009** and the chain moves the pmf by −0.30 / +0.02 / −0.00 bp. It does not damage a correct
estimate.

### ⛔ THREE THINGS STILL WRONG WITH IT

1. ⛔⛔ **`TRIM_Q = 0.50` IS A MAGIC NUMBER and may not ship** (`TRAPS: no-magic-numbers`). The basin is
   wide — q ∈ [0.20, 0.65] all give +4.7 … +10.3 bp at `g05` — but **q = 0.80 BREAKS it** (+67.7 bp,
   `rho_hat` 1.98× truth), because the safe depth is a function of the CONTAMINATED FRACTION, which is a
   property of the library and not a constant. ⭐ **The derived replacement**: include an object iff its
   count is consistent with `Poisson(rho·E_i)`, iterated to a fixed point in `rho`. That is a self-consistent
   trim with no chosen depth, and it is the same shape as the bracketed bisection `EQUATIONS.md` §6b already
   uses for the strand overdispersion.
2. ⚠ **Small-sample wobble.** On the equal-length TEST CHROMOSOME (152 regions) `g05` costs −6.58 bp. The
   ladder's 70,270 regions cut that to −0.30, so it is object count, not bias — but the object count is a
   property of the annotation and must be reported, not assumed.
3. ⛔ **Capture-ON is untouched and must not be claimed.** `TV(r_0, r_1)` is **0.95** under capture — the
   shared-contaminant premise is dead — and the purities stop being separated. Capture-ON is defects **B**
   and **C**, which are independent and own the already-priced −5.90 %.

⛔ **And it has NOT been A/B'd through the shipped pipeline.** Everything above is scored on histograms
against the oracle. `TRAPS: panel-before-src`.

## 6. The order the next session should work in

1. ⭐⭐⭐ **Derive the trim** (§5 item 1) — the Poisson-consistency fixed point. This is the only thing
   between the chain and a shippable estimator, and it is the one place a constant would otherwise enter.
2. **A/B it end to end with `em_fl_ceiling.py`**, on all three test-chromosome fl panels AND the ladder,
   judged on `gdna_frac_est`. The perfect-pmf arm says the ceiling is 8–32 % of the library gDNA bias; the
   chain recovers 89–99.9 % of the pmf error, so most of that ceiling should be reachable — ⛔ but that
   inference is a PREDICTION until the A/B runs, since the pmf error and the composition error are not the
   same quantity.
3. ⛔ **Watch the equal-length control and both sign arms on every run.** §1 and §2 both show a single arm
   giving a confidently wrong sign.
4. **B and C** are independent, cheap, and own the already-priced −5.90 % — the better buy if step 1 stalls.

## 7. ⚠ Three test-chromosome fl panels exist and are documented NOWHERE

`scenarios_fl_rna_long`, `scenarios_fl_gdna_long` and `scenarios_fl_equal200` are on disk, fully cached,
and appear in no permanent doc — not `TESTING.md` §0/§0a, which owns the panels, and not `CLAUDE.md`'s
table, whose only fl-gap row is the **suite-scale** `flgap_*.yaml` pair. ⭐ That omission is what made the
previous session record "the fl-gap arms must be regenerated first" as a blocking gate: it is true of the
suite panels and false of these, and ~11 GB of re-simulation was queued behind a panel that already
existed. They belong in `TESTING.md` §0.
