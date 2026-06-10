# Robust strand-cleaning: shrink to the global gDNA fraction (g₀=global)

**Status:** proposed fix. 2026-06-09. Builds on `strand_cleaning_robustness_design.md` (the options
analysis, which recommended O2 precision-weighted shrinkage) and the g₀=1 trial result.

## 1. Current behavior

The calibrator separates each node's unspliced mass into gDNA vs RNA *by strand*, via a closed-form
unmix of the observed sense fraction `s = sense/total` between gDNA (sense rate ½, unstranded) and
RNA (sense rate κ = `rna_sense_frac`):

```
ĝ = (s − κ) / (½ − κ)          (method-of-moments gDNA fraction)
```

This is unbiased but its variance is `p(1−p)/[N(½−κ)²]` — it **explodes as κ→½** (no strand
information) or N→0 (sparse). We made it robust with a precision-weighted shrinkage (O2): weight ĝ
by its own Fisher information `τ = 4N(½−κ)²` against a prior precision `τ₀`, shrinking toward a
fallback `g₀`:

```
g = [ 4N(½−κ)(s−κ) + τ₀·g₀ ] / [ 4N(½−κ)² + τ₀ ]
```

This cancels the 1/(½−κ) blow-up and degrades smoothly: when strand is strong (`τ ≫ τ₀`) it →ĝ; when
strand is weak (κ→½ or N→0, `τ→0`) it →**g₀**.

**The problem with the current target g₀ = 1.** The trial set g₀=1 (the "keep raw count" fallback).
The full 20-scenario benchmark showed this is **behaviorally inert** — it reproduces the old cliff,
because "default to g₀=1 at low SS" *is* what the cliff did (return 1 at κ=½). At an unstranded node
this means **"assume the mass is all gDNA."** Consequences, concentrated in the unstranded (ss=0.50)
conditions:

| condition | g₀=1 net leak | interpretation |
|---|---|---|
| gdna_none cap_on ss0.50 (0 % gDNA) | **−4.5 %** | pure-RNA library, 4.5 % of RNA falsely called gDNA |
| gdna400 cap_on ss0.50 (80 % gDNA) | +17.4 % | |
| gdna1000 cap_on ss0.50 (90 % gDNA) | +18.7 % | |

At zero contamination, defaulting unstranded nodes to "all gDNA" manufactures a false-positive
RNA→gDNA siphon. The default is wrong because it ignores *how much gDNA the library actually has*.

## 2. The proposed fix: shrink to the global gDNA fraction

Make the shrinkage target the **library-wide gDNA fraction**, node-type-aware:

- **Count-observable regions** (intron / intergenic) — gDNA-pure by signature ⇒ **g₀ = 1** (their
  unstranded default of "all gDNA" is *correct*; unchanged from today).
- **Non-observable (exon) regions & boundary sides** — a genuine gDNA/RNA mixture whose split is
  unidentifiable at the node when strand is weak ⇒ **g₀ = ĝ_global**, the library's overall gDNA
  fraction.

So an unstranded exon node no longer defaults to "all gDNA" — it defaults to the **actual
contamination rate** of the library.

### Estimating ĝ_global without strand (no circularity)

The global gDNA fraction must be available *before* the per-node clean, and must not depend on the
clean itself. It is — from the **count-observable regions**, which are gDNA-pure *by signature*
(intergenic / intron-only have no exon bit, so they carry no mature RNA). Their unspliced density is
the gDNA density **with no strand cleaning required**:

```
ρ_global = Σ unspliced_count(observable regions) / Σ eff_len(observable regions)
```

The expected gDNA fraction of any node is then `ρ_global · eff_len / count` — i.e. how much of the
node's mass is explained by gDNA at the library density. Aggregated, `ĝ_global` is the library gDNA
fraction. This is **strand-free** (observable regions need no strand), breaking the circularity (we
do *not* need the post-deconvolution global density as the target).

## 3. Theory — why this is the right target

This is **empirical-Bayes shrinkage toward the population mean**. The strand clean is a noisy
per-node measurement of g; its likelihood precision is `τ = 4N(½−κ)²` (the strand Fisher
information). With a prior `g ~ N(g₀, 1/τ₀)`, the posterior mean is the precision-weighted blend
above. The *correct* prior mean for a node about which strand says nothing is the **expectation of g
across the population it is drawn from** — the library gDNA fraction. Not 1 (that asserts every node
is pure gDNA), not ½ (that asserts maximal ignorance even when the library is, say, 2 % gDNA), but
the empirical rate.

The key statistical fact: **at zero strand information, a node's best estimate of its gDNA fraction
is the population gDNA fraction** — exactly what an EB prior provides. g₀=1 and g₀=½ are both *fixed*
priors that ignore the data's own contamination level; ĝ_global is the *empirical* prior that adapts
to it.

## 4. What it does / why it works

- **gdna_none (0 % gDNA):** ρ_global ≈ 0 ⇒ ĝ_global ≈ 0 ⇒ unstranded nodes default to ~0 % gDNA ⇒
  the **false-positive RNA→gDNA siphon disappears** (the −4.5 % → ~0).
- **gdna1000 (90 % gDNA):** ĝ_global ≈ 0.9 ⇒ unstranded nodes default to ~90 % gDNA ⇒ the unstranded
  over/under-call is pulled toward the truth (the +18.7 % shrinks toward the correct split).
- **Stranded (ss=0.99):** `τ = 4N(0.49)² ≈ N ≫ τ₀` ⇒ the prior is irrelevant ⇒ g →ĝ ⇒ **unchanged**
  (the stranded win, and the gdna_none ss0.99 ~0 FP, are preserved).
- **Observable regions:** g₀=1 keeps their correct gDNA-pure default at any SS.

The shrinkage onset is automatic and per-node: a node with abundant well-stranded fragments follows
its own data; a node that is sparse *or* unstranded leans on the library rate, in exact proportion
to `(2κ−1)²·N`.

## 5. Benefits

- **Fixes the unstranded false-positive** (the headline failure: pure-RNA libraries no longer
  manufacture gDNA) and improves the unstranded leak at all contamination levels.
- **Robust across the whole matrix** — contamination level (the target *is* the contamination rate),
  strand specificity (smooth κ→½), and sparsity (smooth N→0). Unlike the global-FL idea (Issue #3,
  which only worked because gDNA dominated), the target here is *derived from* the contamination
  level, so it is correct at low and high gDNA alike.
- **Principled, minimal magic numbers** — `ĝ_global` is data-derived (no constant); the precision
  weight `(2κ−1)²·N` is the Fisher information (no constant); the only knob is τ₀ (the prior
  strength, ≈ "one stranded fragment"), and it only bites where strand is genuinely uninformative.
- **No high-SS regression** (proven by the g₀=1 trial — the shrinkage is inert when strand is
  strong; changing the *target* only changes the low-SS fallback).
- **Cheap** — one extra strand-free global reduction; the per-node formula is vectorized.

## 6. Implementation sketch

1. In `node_gdna_density`, compute `ρ_global` (strand-free, from count-observable regions) → derive
   `ĝ_global` per node as `clip(ρ_global·eff_len/count, 0, 1)` (this is ~1 for observable regions and
   the count-based fraction for exons — node-type falls out, no explicit branch).
2. Pass `g0 = ĝ_global` into `strand_clean_gdna_frac` (currently hard-wired to 1.0); shrink toward it.
3. Carry `ρ_global` to the boundary-side clean (`_compute_side`) so sides use the same target.
4. Keep `τ₀ = strand_clean_prior_weight` (config; trial default 1.0). `τ₀=0` still reproduces legacy.

## 7. Validation plan

Re-run the full 20-scenario benchmark. Acceptance:
- **gdna_none ss0.50 FP: −4.5 % → ~0** (the headline win).
- unstranded+capture leaks (gdna100/400/1000 ss0.50) shrink toward the truth.
- **ss0.99 conditions unchanged** (win preserved); gdna_none ss0.99 still ~0.
- Add a unit test: `g_robust → ĝ_global` as κ→½, `→ ĝ` as τ≫τ₀, finite everywhere.

## 8. Relationship to the larger restructure

This is the *incremental* fix within the current acyclic architecture (pre-clean → joint deconv).
The larger alternative — a single robust **strand deconvolution** step feeding a **count** step (no
pre-clean) — remains on the table as a future simplification; g₀=global is the low-risk move that
makes the current pre-clean robust without that restructure.
