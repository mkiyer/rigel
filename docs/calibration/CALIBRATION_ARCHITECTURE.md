# Calibration Architecture — the information & precision model (AUTHORITATIVE)

**Status:** the foundational architectural principles of Rigel calibration (2026-06-17). This is the
**reference any future work must internalize first.** It states *what information exists, where precision comes
from, and how a node is solved* — the invariants that the mechanism specs sit on:
- `CALIBRATION_PLAN_v6.md` — the bipartite region/boundary node graph (the mechanism).
- `rna_imputation_transcript_structure.md` — the transcript-structure rules for RNA imputation.

If a mechanism spec ever contradicts this document, **this document wins** and the spec is wrong. The bugs that
motivated this doc were all violations of §0.

---

## 0. THE ONE PRINCIPLE EVERYTHING FOLLOWS FROM

> **A fragment COUNT carries ZERO intrinsic information about a node's gDNA/RNA composition.**

A node observes `N` unspliced fragments, split `(n_pos, n_neg)` by genome strand. Considered **alone**, in
isolation, this tells you **nothing** about `{f_pos, f_neg, f_gdna}`. Example: 100 unspliced, 50 genome-`+` /
50 genome-`−`. That 50/50 is *equally consistent* with pure gDNA (genomic strand is symmetric, sense ½) and pure
RNA at an unstranded library (κ=½). **The count alone cannot say "this is gDNA" or "this is RNA."**

The count's *only* intrinsic role is to set the **precision of the STRAND signal** — more fragments ⇒ a tighter
strand posterior. There is no other legitimate way for the raw count to influence a node's solution or its
precision.

**Molecular sampling.** Re-running the experiment yields a different total and split (103; 51/52) — and for
RNA-seq this is **overdispersed**. The strand model **must** include this overdispersion (Beta-Binomial). The
count's contribution to precision is therefore the *overdispersed* Fisher information, not a naïve Poisson.

> **THE INVARIANT (violating it is a bug):** *No code path may let the raw fragment count determine a node's
> gDNA/RNA fraction, or the precision of that fraction, EXCEPT through the strand model's (overdispersed) Fisher
> information.*

Every calibration bug we have chased traces to a violation of this invariant (the count prior voting gDNA at
introns; `mass_u` Poisson caps on the fraction precisions; the DIRECT own-count `var~mean` used as a
fraction-confidence; the `I₀` count→strand suppression).

---

## 1. THE THREE — AND ONLY THREE — INFORMATION SOURCES

A node's `{f_pos, f_neg, f_gdna}` is determined by exactly these, in priority order:

### 1.1 The STRAND signal — the only INTRINSIC per-node signal
The Beta-Binomial tilt of `(n_pos, n_neg)`: **RNA tilts the genome-strand rate toward κ** (the library sense
fraction), while **gDNA stays at ½** (genomic strand symmetry). The mixture plus-strand rate is
`p = ½·f_gdna + κ·f_+ + (1−κ)·f_−`, which depends only on the tilt `t = f_+ − f_−`. So the strand constrains
the gDNA-vs-RNA split *and* the ± split, **with no help from the count except its statistical power**:
- Precision ∝ the overdispersed BB Fisher information `~ N·(2κ−1)²`. High at high N + high strand specificity;
- **Uninformative (near-infinite variance, zero precision) at κ=½ (unstranded) or low N.** Then the node has
  NO intrinsic signal and must be solved from outside (§1.2, §1.3).

This is the **PRIORITY** signal and the **ONLY** legitimate count→precision path.

### 1.2 IMPUTATION across the node chain — OUTSIDE information
A neighbour node with **higher precision** propagates its solution onto a lower-precision node. The
**precision of the propagation is a MODELED RELIABILITY** — the `var~mean` of the prediction error — which
depends on count **density**, *not* the raw count as a confidence. Empirically **imprecise**. It is
transcript-structure-aware (§4): gDNA propagates genomically (smooth); nascent RNA propagates within a
transcript across exon↔intron; mature RNA is anchored, not propagated; nothing propagates across a TSS/TES.

### 1.3 The GLOBAL gDNA hyperprior — OUTSIDE information (the foundation)
The population baseline gDNA density; it sways an *unanchored* node toward gDNA when there is global evidence of
gDNA elsewhere. Precision = the robust **MAD population spread** of the per-node densities (NOT a count). It is
the stabilizing foundation for nodes the strand and imputation both leave undetermined.

**A node with an informative strand is solved by the strand. A node with an uninformative strand (unstranded /
AMBIG / low count) is solved by imputation + the global hyperprior. The COUNT never independently solves a node.**

---

## 2. THE PRECISION MODEL — where variance is allowed to come from

Each node carries a per-component variance `{var_pos, var_neg, var_gdna}` (precision = 1/var). Variance is set
ONLY by:

| source | what it sets | count-dependence |
|---|---|---|
| **Initialization** (§3) | the signature-binary starting variance | **NONE** (signature only) |
| **Strand deconvolution** | updates variance after the BB solve | **YES, and only here** — the overdispersed Fisher info `N·(2κ−1)²` |
| **Imputation** | the propagated prior's precision | the `var~mean` **reliability** (a measured prediction error; density-based, NOT a raw-count confidence) |
| **Global hyperprior** | the foundation precision | the **MAD** population spread (NOT count) |

**FORBIDDEN precision sources (the contamination to delete — §6):**
- the raw count (`mass_u`) **capping** a fraction precision (`min(·, mass_u)` on `τ_count`/`τ_global`/`τ_rna`) —
  the count as a precision *ceiling*;
- **count-derived SOURCE variances/means used as composition votes**: the DIRECT own-count `var~mean` (a node's
  own-density Poisson disagreement) used as a *fraction*-confidence (it measures how precisely you counted, NOT
  how sure you are of the gDNA *fraction*); the count-prior `count_gdna_frac`/`τ_count` vote; the RNA-prior
  count-derived mean/precision (§6.1, the symmetric twin);
- the **`I₀` count→strand suppression** (`w=I/(I+I₀)`) — a hand-tuned weight; deference must be **emergent**
  (a sharp strand likelihood dominates a diffuse prior by their honest precisions, with no tuned constant).

> **NOT a violation — the definitional density↔fraction jacobian.** The change-of-variables `(eff_len/mass_u)²`
> that converts a *density*-variance to a *fraction*-variance is the SAME `mass_u` that **defines** the fraction
> `f_g = M_gdna / mass_u` (mass-over-observed-total). It is the legitimate normalizer, not a count-confidence — so
> the surviving GLOBAL prior keeps the jacobian + its `mu_global = ρ_global·eff/mass_u` normalization, and only
> the **caps** and the count-derived **sources** are removed. The §0 violation is the count as a *ceiling* or a
> *composition vote*, never the count as the *normalizer of a fraction*. (Decision 2026-06-17; resolves the
> two-critic split — the jacobian was wrongly conflated with the cap in an earlier draft.)

---

## 3. INITIALIZATION — signature-binary, derived from first principles

Every node starts at the all-gDNA default with a **signature-determined** (never count-determined) variance.

**Fractions:** `{f_pos=0, f_neg=0, f_gdna=1}` everywhere (the all-gDNA default).

**Variances** (4-bit signature: `ex+=0x2, ex-=0x1, in+=0x8, in-=0x4`; +bits=`ex+|in+`, −bits=`ex-|in-`):
- `var_gdna = MAX` (numerically-stable high) **always**;
- `var_pos = MAX` iff the node has **any +bit** (a +transcript is present), else **0** (locked: no +RNA possible);
- `var_neg = MAX` iff the node has **any −bit**, else **0** (locked).

By strand class (these coincide with `transcript_strand_class` / `allow_pos`/`allow_neg`):

| class | signature | var_pos | var_neg | var_gdna | meaning |
|---|---|---|---|---|---|
| **Intergenic** (TS_NONE) | no bits | 0 | 0 | MAX | `f_pos=f_neg=0` locked ⇒ `f_gdna=1` locked by the simplex constraint. **Pure gDNA, signature-KNOWN, zero variance.** |
| **single-strand +** (TS_POS) | +bits only | MAX | **0** | MAX | `f_neg=0` locked; `f_pos`/`f_gdna` unknown — solved by the strand. |
| **single-strand −** (TS_NEG) | −bits only | **0** | MAX | MAX | symmetric. |
| **AMBIG** (TS_AMBIG) | +bit and −bit | MAX | MAX | MAX | fully unknown — solved by strand (weak) + imputation. |

**The count plays NO role in initialization.** `var_gdna=MAX` always: gDNA is the *residual* of the simplex, so
locking the strand fractions (var=0) locks gDNA by construction (intergenic), while leaving them free (var=MAX)
leaves gDNA unknown. This replaces the current count-based all-gDNA bootstrap (`gdna_c = u_tot`), whose pass-0
precision was count-derived and count-capped — the locus of the bug.

---

## 4. THE SOLVE FLOW (iterated)

```
1. INIT            — signature-binary fractions + variances (§3). NO count.
2. STRAND DECONV   — the BB tilt solver updates {fractions, variances} per node where the strand is informative.
                     COUNT → PRECISION happens ONLY here (overdispersed Fisher info).
3. GLOBAL prior    — compute the baseline gDNA density from the gDNA EVIDENCE (strand-deconvolved gDNA +
                     signature-known intergenic); apply it as a prior to UNANCHORED nodes (MAD precision).
4. IMPUTATION      — propagate higher-precision neighbours' pies to lower-precision nodes along the bipartite
                     chain (var~mean reliability; transcript-structure-aware §5).
5. ITERATE 2–4     — to convergence; the global + imputation sharpen as the gDNA estimate refines.
```

**The global hyperprior is computed and applied AFTER the strand deconvolution, iteratively — NOT at init.** At
initialization there is no gDNA evidence yet (every node is the all-gDNA *unknown* default); the baseline only
becomes meaningful once the strand deconvolution has separated gDNA from RNA where it can. So the global is a
per-iteration foundation derived from the running gDNA estimate, never an initial value.

---

## 5. THE BIPARTITE NODE GRAPH + TRANSCRIPT STRUCTURE (mechanism recap)

Detailed in `CALIBRATION_PLAN_v6.md` (graph) and `rna_imputation_transcript_structure.md` (RNA rules):
- **Linear bipartite chain** `region ↔ boundary ↔ region`; both are first-class nodes owning mass, eff-len, a pie.
- **Boundary nodes own the one-sided, motif-stranded SPLICED mass** = a FIXED mature-RNA floor (never deconvolved;
  it anchors, not propagates).
- **Imputation respects transcript structure:** TSS/TES (intergenic↔exon) = a zero-RNA black hole (no RNA across);
  exon↔intron carries the **unspliced/nascent only** (spliced is one-sided, skips the intron); an exon's contained
  mature is anchored by the spliced floor + strand, not imputed from the thin junction crossing.
- **gDNA is genomically smooth (chains); RNA is spiky** (nascent contiguous *within* a transcript; mature
  contained in exons). This is why gDNA imputation propagates broadly while RNA imputation is transcript-gated.

---

## 6. WHAT IS REMOVED (the contamination) AND WHY

All of the following violate §0 (the count, or a tuned constant, illegitimately determining a fraction or its
precision). They are leftover pre-v6 count-channel / suppression code, NOT part of this architecture:

1. **The count prior AND its symmetric RNA twin** — `count_gdna_frac`/`τ_count` (the gDNA fraction vote) AND the
   **RNA prior** `_rna_prior_fraction` (the `f_pos`/`f_neg` vote: mean built from the raw count, precision via
   `geom_rna=(eff_rna/mass_u)²` + a `mass_u` cap). The dry-run (2026-06-17) found the RNA prior is the *exact*
   symmetric count→fraction-precision twin and was missing from this list; because `f_pos+f_neg+f_g=1`, pinning
   the RNA fractions from the count re-pins `f_g` by the back door. **Decision: remove the RNA prior SYMMETRICALLY
   with the count prior** (not de-count-and-keep) — the strand alone already solves the diagnostic introns to
   `f_g≈0` (the count prior was the pull-*up*, the RNA prior the counter-crutch pull-*down*; remove both and the
   strand wins). Removing the RNA prior orphans **`rna_density_model.py`** (its only production consumer) — delete
   the module; Step 2 rebuilds RNA imputation per the corrected transcript-structure spec. The count's magnitude
   survives ONLY as the **density normalizer** (the jacobian, §2) and as the strand BB Fisher info — never as a
   confident "this is gDNA/RNA" vote. (Intergenic is gDNA by **signature** (§3), not by the count.)
2. **The `mass_u` Poisson CAPS** on `τ_count`/`τ_global`/`τ_rna` (`calibrate.py`): the count as a precision
   *ceiling*. (The density↔fraction jacobian itself is KEPT for the global prior — see §2; only the caps go.)
3. **The DIRECT own-count `var~mean`** used as a fraction-confidence (`variance_model.direct_points` /
   `fit_direct_varmean`): the own-density disagreement is the count precision, not the gDNA-fraction precision.
   Folded INTO the count-prior removal (its only consumer `var_d→sig2_frac→τ_count` dies with the count prior).
4. **`I₀` (`gdna_strand_info_scale`) in the SWEEP ONLY** (PATH A: `w_strand=I/(I+I₀)`, `eff_count_precision`) and
   its `(1−w_strand)` suppression — the hand-tuned count→strand weight; deference is **emergent**. **The config
   field + PATH B (`deconv_sides`/`cleaned_gdna_count`, which re-validate `>0`) stay until Step 3** when the
   boundary nodes retire them — deleting the field in Step 1 is a Step-1/Step-3 ordering conflict.
5. **`q_rna = 0.25`** + the chain coupling machinery (`_sweep_chain`/`_edge_logphi`): a magic edge-coupling
   constant. **Decision: delete in Step 1** (it is not a count→precision violation, but its principled
   replacement — reliability-weighted imputation — lands in Step 2; deleting it now collapses the chain to a
   direct per-node posterior). The complex-loci AMBIG regression is pre-accepted and **measured** on the battery.
6. **The dead density-model variance** (`_count_fraction_variance`, `count_gdna_frac_var`, `count_rel_var` in
   `density_model.py` ONLY): superseded duplicate, off on the default path. Do NOT touch the same-named
   `strand_deconv._SideQuantities.count_gdna_frac_var` (a separate live field — name collision). Confirm the
   `gdna_deconv_quantile` knob's disposition (the non-default-quantile path is the only consumer).

**KEEP:** the strand model (BB, overdispersed — the only count→precision); `count_gdna_frac` **MEAN** as the gDNA
strand-overdispersion fit **seed selector** (`gdna_strand.py`) + `count_observable_masks` as a *signature
partition* (intergenic vs transcript — NOT a gDNA-confidence); the `var~mean` IMPUTATION machinery
(`MonotoneVarMean` / `pair_imputation_points` / `fit_pair_imputation_varmean` — production-uncalled after Step 1
but the Step-2 substrate); the GLOBAL prior (MAD spread + jacobian-normalized mean, §2); the signature-binary init.

**The `U>0` active gate** (`simplex_sweep.py`, forces `f_g=0`/`var=0` at zero-count nodes) is a latent count→
fraction violation but is **inert today** (`gdna_mass=f_g·mass_unspl≈0` at zero count) — NOTED, deferred to Step 3
(the variance-state object), to avoid gratuitous golden churn now.

---

## 7. MAGIC-NUMBER LEDGER (the end state)

After this architecture is implemented, calibration carries **no count→fraction-precision and no hand-tuned
suppression constant**. Every remaining constant is **canonical** (Jeffreys ½, MAD `1.4826`), a **numerical
resolution** (lattice `K`, GCV grid, `_EPS`, the MAX/0 init variances are numerically-stable sentinels), or a
**measured reliability** (the `var~mean`). `I₀` and `q_rna` are deleted.

---

## 8. IMPLEMENTATION ORDER (recommended)

Each step independently testable; the zero-gDNA diagnostic (`scripts/debug/rna_imputation_diagnostic.py`, true
`f_g=0` everywhere) + the net-flow before/after are the gates. A temporary recalibration regression during the
count-prior/`I₀` removal is **expected and pre-accepted** (recover it by improving the reliability layer, never
by re-adding a suppression constant).

- **Step 1 — the precision rebuild (the core fix; resolves the intron bug).** Make the init signature-binary
  (§3); **remove the count prior + its RNA twin** and ALL count→fraction-precision (§6.1–6.3); **delete `q_rna`**
  and `I₀`-in-the-sweep (§6.4–6.5) and the dead density-variance (§6.6). The strand becomes the primary per-node
  solve; the global handles strand-weak nodes (imputation returns in Step 2). *Gate: zero-gDNA introns solve to ~0
  via the strand; net-flow non-regressing on the gDNA-present conditions.* Mostly **removal + the explicit init**
  inside the existing region sweep — the most tractable starting point, and it directly kills the bug.

  **Init scope = MINIMAL (decision 2026-06-17).** §3's per-node `{var_pos,var_neg,var_gdna}` *state object* does
  NOT exist in the current sweep — the sweep is a **stateless sum-product** (ψ rebuilt each pass, static lattice;
  the only inter-pass state is the gDNA *mass*). The locks are already realized by the hard `allow_pos`/`allow_neg`
  forbid mask (which coincides exactly with the §3 table). So Step 1's "signature-binary init" is: **replace the
  count-based seed `gdna_c = u_tot.copy()` with the count-free mass seed `gdna_c = mass_u.copy()`** (so `f_gdna=1`
  means the full unspliced *mass*, not the *count*) and audit the seed's downstream feeders (`ρ_global`,
  `active_mu` gate, `region_fg`) to ensure no `u_tot` count leaks into the refit inputs. §3's MAX/0 variances are
  **conceptual** for the current sum-product; the explicit per-node variance-state object arrives with the
  first-class bipartite nodes in **Step 3**.

  **Substep sequence (green-suite checkpoint after each):** (0) `q_rna` + chain-machinery teardown — the per-ref
  chain loop collapses to a direct per-node posterior from ψ; golden update + measure the complex-loci AMBIG
  regression. (1) signature-binary count-free init — near-no-op (forbid mask unchanged). (2) **one commit**:
  remove the count prior + RNA prior + DIRECT var~mean + de-count the global precision (drop the `mass_u` caps,
  keep the jacobian) + delete `rna_density_model.py`; gate on the zero-gDNA diagnostic (introns→0 via
  strand+Jeffreys+global) + net-flow before/after. (3) dead density-variance cleanup.
- **Step 2 — the imputation source = strand-deconvolved gDNA.** Ensure the gDNA imputation propagates the
  **strand-deconvolved / signature-known** gDNA density (clean), not the raw count density; transcript-structure
  RNA imputation already lands here (the RNA-prior gate fix). *Gate: weak-SS / unstranded introns rescued by
  imputation + global without a count prior.*
- **Step 3 — the bipartite boundary nodes (v6 Phase C).** Promote boundaries to first-class solved nodes with the
  one-sided spliced floor; retire `deconv_sides`/`cleaned_gdna_count`; the per-node pie + the signature-binary
  init + the strand/imputation/global apply uniformly to region and boundary nodes.

---

## 9. THE LITMUS TEST (for any future change)

Before adding or keeping any term, ask: **"Does this let the raw fragment count determine a fraction or its
precision, other than through the strand's overdispersed Fisher information?"** If yes, it is a bug. The only
count→precision is the strand. Everything else is signature (init), measured reliability (imputation), or the
population spread (global).
