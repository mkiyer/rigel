<!-- title: Step 2 — Jacobian → count-space (implementation plan, pending open-question resolution) -->
# Step 2 — get out of `(M/E)²` Jacobian space into count space

**Status:** detailed implementation plan. **DO NOT implement yet** — §6 has blocking open questions (the
global is genuinely hard). Companion: `honest_precision_unified_design.md` (the why). Step 1 (σ²_bio raw fit)
is **dropped** (no `f_g`-independent gDNA density exists; the belief-derived σ²_bio is strand-anchored and
converges — see the conversation). So Step 2 reuses the **existing** σ²_bio fit unchanged.

**The goal, stated narrowly (user):** the message, the global prior, and their delivery currently convert a
density precision to a *fraction* precision by multiplying by **`(M_dst/E_dst)²`**. That over-amplifies a prior's
weight by the destination mass-squared, while the strand likelihood — the honest exemplar — weights by the
count linearly. **Replace every `(M/E)²` precision with a count-space pseudo-count** whose strength is set by
counts and length (pseudo-fragments competing one-for-one with the strand likelihood's real fragments), **and
do it at every site the Jacobian exists.** Nothing else changes.

---

## 1. The complete Jacobian enumeration (grep-verified — this is the whole list)

The `(M/E)²` precision-amplifier exists at **exactly three sites**, all in `bp_solver.py`, landing at three
application sites in `simplex_sweep.py`:

| # | site | code (verbatim) | role |
|---|---|---|---|
| **J1** | `bp_solver.py:575` | `tau_f = tau_rho * (mass_dst / max(eff_dst,_EPS)) ** 2` | message density→fraction precision |
| **J2** | `bp_solver.py:581` | `var_floor = mu_f*(1-mu_f)*(1/mass_src + 1/mass_dst)` | the wall-vanishing partial cap (paired with J1) |
| **J3** | `bp_solver.py:636,710` | `geom2_global=(eff_global/mass_global)**2`; `tau_global=1/(sigma2_n*geom2_global)` | global density→fraction precision |
| applied | `simplex_sweep.py:84` | `psi -= 0.5*ti*(f_g - gdna_imp_frac)**2` | gDNA message Gaussian |
| applied | `simplex_sweep.py:95` | `psi -= 0.5*tr*(f_axis - mu_s)**2` | per-strand RNA message Gaussian |
| applied | `simplex_sweep.py:76` | `gpen = -0.5*gt*(f_g - global_mu)**2` | global Gaussian |

**What is NOT a Jacobian (stays untouched) — verified each `**2`/`mass/eff` hit:**
- `μ = clip(ρ·E/M, 0, 1)` mean-conversions (`bp_solver.py:566,708`; `density_model`, `strand_deconv`, `derive`):
  the **mean** of a density-imputed fraction *should* depend on `M/E`. Only the **precision** `(M/E)²` is wrong.
- σ²_bio fit response `(ρ_dst−ρ_src)²` (`bp_solver.py:504,705`) — the var~mean data; kept (Step 1 dropped).
- strand model variance `(n·f)²·od` (`simplex.py:74-76`, `strand_likelihood.py:72-73`) — honest count-space.
- spliced floor `short_p²` (`simplex_sweep.py:64`) — count-weighted hinge, **not** a `(M/E)²` site; kept.
- EM eff-len IPR `g²/supp` (`capture_eff_length.py:164`, `priors.py:176-177,290`) — the *EM-side* effective
  length, a different track; out of scope.
- `_fg_var` readout `post@f² − mean²` (`simplex_sweep.py:119`); `strand_deconv.py:138` count-clue variance —
  not sweep priors.

So Step 2 touches **J1, J2, J3 + their three application sites + the strand_obs fork**, and nothing else.

---

## 2. The unifying count-space form (the replacement for every Jacobian)

Every prior becomes a **two-component (this-component vs rest) binomial pseudo-count** on its own axis:

```
prior_c(f_c | μ_c, N)  =  α_c · log f_c  +  (N − α_c) · log(1 − f_c) ,   α_c = N · μ_c
```

- **mode at `f_c = μ_c`** (∂=0 ⇒ α(1−f)=(N−α)f ⇒ f=α/N=μ) — two-sided, pulls toward the imputed value.
- **curvature at the mode = −N/(μ(1−μ))** — the binomial Fisher information of **N pseudo-fragments**. The
  weight is **N (a count)**, with **no `(M/E)²`**. It competes one-for-one with the strand likelihood, whose
  info is `∝ N_dst` (the node's real fragments).
- **pins a wall via the complement:** as `μ_c → 0`, `α_c → 0` and the term → `N·log(1−f_c)` → pulls `f_c → 0`
  with strength `N`. (This is the piece the rejected *one-sided* `α_c·log f_c` lacked — it's why that form was
  a monotone push to the vertex. The complement makes it honest and two-sided.)
- **length** enters only through the **mean** `μ_c = ρ·E/M` (the imputed fraction), exactly where it belongs;
  the **count** `N` sets the weight. Counts + length, no Jacobian. ✔ the stated goal.

`f_c` is clipped to `[_PRIOR_EPS, 1−_PRIOR_EPS]` (the existing Jeffreys clip) for the `log`. New helper in
`simplex_sweep.py`:
```python
def _binom_pseudo(f, mu, N):   # N, mu are per-node columns; f is the lattice row
    a = N * mu
    fc = np.clip(f, _PRIOR_EPS, 1.0 - _PRIOR_EPS)
    return a * np.log(fc) + (N - a) * np.log1p(-fc)
```

---

## 3. The MESSAGE transformation (confident — no open questions)

**`_message` returns `(μ_c, N_src)` instead of `(μ_f, τ_f)`; delete J1 + J2.**
```python
def _message(rho_src, eff_src, eff_dst, mass_dst, rho_dst_cur, varmean, spliced_dst=0.0):
    mu_c   = clip((rho_src*eff_dst - spliced_dst)/max(mass_dst,_EPS), 0, 1)   # UNCHANGED — imputed composition
    mu_q   = max(0.5*(rho_src + rho_dst_cur), _EPS)                            # UNCHANGED — σ²_bio query (belief-based; Step 1 dropped)
    sigma2 = varmean.predict([mu_q])
    pois   = (rho_src + _MSG_PSEUDOCOUNT/max(eff_src,_EPS))/max(eff_src,_EPS)  # UNCHANGED — (n_src+1)/E_src²
    tau_rho= 1/max(sigma2 + pois, _EPS)                                        # UNCHANGED — density precision
    N_src  = rho_src**2 * tau_rho                                              # NEW — source effective count (1/CV²_total), M_dst-INDEPENDENT
    return mu_c, N_src
    # DELETED: tau_f = tau_rho*(mass_dst/eff_dst)**2  (J1) ; var_floor + floor compose (J2)
```
**Why `N_src = ρ_src²·τ_ρ` is right:** when sampling-dominated (small σ²_bio), `N_src ≈ ρ_src²·E_src²/(n_src+1)
= n_src²/(n_src+1) ≈ n_src` — the **source's gDNA count**. The message carries ≈ the source's fragments, at
composition μ_c, M_dst-independent. Applied via §2: `psi += _binom_pseudo(f_g, μ_c, N_src)` (gDNA) and the same
on `f_pos`/`f_neg` for the per-strand RNA messages. `mass_src` drops out of the signature (was only in J2).

**Behaviour check (stranded-0):** the AMBIG phantom node has a *decisive* strand likelihood (the strand pins
`f_g→0`), and in zero-gDNA the messages fade (`ρ_src→0 ⇒ N_src→0`). So stranded-0 stays ≈0 via the strand,
independent of the global. ✔

---

## 4. The GLOBAL transformation — derived (length-driven zero-pin), and the coupling it exposes

**Decision (user): derive the length-driven zero-pin first.** Done. Here is the derivation, the exact form, and
the one non-obvious thing it forced into scope.

### 4.1 The model (count space, no fraction anywhere)
The node's gDNA fragments are a Poisson count over the gDNA exposure `E_g`, at a rate that varies between nodes:
```
g_i ~ Poisson(ρ_i · E_g) ,   ρ_i ~ Gamma(mean = ρ_global, var = σ²_g)
```
`ρ_global` = the exposure-pooled population rate (the existing `rho_mean`); `σ²_g` = the **gDNA between-node
density spread** (the var~mean). Marginalising `ρ_i` gives a **Negative-Binomial** prior on the count `g_i`:
```
g_i ~ NB(r, p) ,   m = ρ_global·E_g ,   r = ρ_global²/σ²_g ,   p = ρ_global / (ρ_global + σ²_g·E_g)
```
Applied to the lattice as a log-prior on `f_g` (with `g = f_g·M`, `M` observed):
```
glob(f_g)  =  gammaln(f_g·M + r) − gammaln(f_g·M + 1) + (f_g·M)·log(1 − p)        [+ f_g-independent const]
```
**This is Jacobian-free:** `M` enters **linearly** inside the log (via `g = f_g·M`), never as `(M/E)²`. Length
enters honestly through `m = ρ_global·E_g` and through `p`. The only constant is the Gamma shape pseudo-count
`a = 1` (reused `_MSG_PSEUDOCOUNT`), added to `r` for well-behavedness: `r = a + ρ_global²/σ²_g`. Mode ≈
`μ_global = m/M`; two-sided.

### 4.2 The zero-pin and capture-deference are BOTH governed by `σ²_g` — and that is correct
The down-pull toward `f_g=0` is the `(f_g·M)·log(1−p)` term. Its strength is set by `1−p = σ²_g·E_g /
(ρ_global + σ²_g·E_g)`:
- **zero-gDNA, `σ²_g → 0`:** `p → 1`, `log(1−p) → −∞` ⇒ **strong length-driven pin** `f_g → 0`. ✔ (the pin
  scales with `E_g` — a long zero-count region is confidently gDNA-free: exactly the user's length intuition.)
- **capture, `σ²_g` large:** `p → 0`, `log(1−p) → 0` ⇒ the global **defers** (no down-push), letting the
  messages (real `N_src` from captured crossings) carry the enrichment. ✔

So the NB form does the right thing in both regimes, with the regime selected by `σ²_g`. **No `(M/E)²`, no new
free constant.** This is the length-driven zero-pin, derived.

### 4.3 The coupling this forces into scope, and the seed set that resolves it (user-specified)
The pin **requires both `ρ_global` and `σ²_g` to be honest** (`μ_global = ρ_global·E/M` must → 0 to pin, and
`σ²_g` must → 0 to make `p→1`). Both are currently **circular** — fit on the belief `ρ_g = f_g·M/E`, so a
phantom `f_g` inflates them and the phantom self-sustains. The fix is to estimate the global's mean **and**
spread from **gDNA-clean, non-circular seed nodes**, fit **once before the sweep** (the inputs are
belief-independent ⇒ no per-pass refit ⇒ the circularity is broken at the root, and the global becomes a fixed
population anchor rather than a moving target).

**The seed set (user-specified):**
- **Structural seeds (always — the only path for UNSTRANDED data):**
  - intergenic regions (`M` ≈ pure gDNA) and intron regions (**assumption: nascent RNA is sparse** ⇒ ≈ gDNA);
  - **intergenic-exon and intron-exon boundary crossings.** The unspliced mass spanning these seams is gDNA:
    spliced RNA reads the junction in the *spliced* channel (not the unspliced crossing), and nascent is sparse.
    Crucially, one side is an **exon under the capture probe**, so these boundaries are the **only gDNA-clean
    nodes that sample the capture-ENRICHED (high-μ) density** — intron/intergenic *interiors* are **depleted**
    under hybrid capture, so without the boundaries the high-μ curve would be unobserved. This is what lets
    unstranded data estimate gDNA under capture.
  - gDNA density for these = the **observed** `M/E` (RNA-free by structure ⇒ non-circular, no contamination).
- **Strand seeds (stranded data only):** all single-strand (POS/NEG) nodes, using their **strand-deconvolved**
  gDNA mass `f_g(strand-only)·M` (the init G2 strand solve — strand-anchored, **not** the iterated belief ⇒
  non-circular). These reach the full enriched density range directly. To avoid a hard κ threshold, **weight**
  each strand seed by its strand-deconvolution confidence (`∝ counts·(2κ−1)²`) so it **fades to 0 as κ→½**
  (magic-free): stranded libraries get the extra coverage, unstranded libraries fall back to the structural
  seeds automatically.

**Result on the failure faces:**
- **unstranded-0:** structural seeds are all ≈0 gDNA (no fragments) ⇒ `ρ_global→0`, `σ²_g→0` ⇒ the pin fires ⇒
  the zero line holds. ✔
- **capture-unstranded:** intron/intergenic interiors depleted, but the **boundary crossings** carry the
  enriched-edge gDNA ⇒ `ρ_global`/`σ²_g` reach the high-μ regime ⇒ the global can recover part of the −367K
  siphon (the rest is the FL floor). ✔
- **capture (stranded):** strand seeds span the full enriched range ⇒ honest high-μ `σ²_g`, no extrapolation. ✔
- **residual:** for capture-UNSTRANDED, the boundaries reach the exon-*edge* density (crossing eff-len), which
  may sit below the exon-*interior* density — a **smaller** extrapolation than intron-only, flagged + measured.

### 4.4 Resulting Step-2 scope (three pieces, ONE gDNA currency, one constant `a=1`)
1. **Message** — `_binom_pseudo(f_c, μ_c, N_src)`, `N_src = ρ_src²·τ_ρ` (§3). Solid.
2. **Global** — the NB-count log-prior (§4.1), pinning/deferring via `ρ_global`, `σ²_g` (§4.2).
3. **The honest gDNA `σ²_g` (and `ρ_global`)** — fit **once** on the §4.3 seed set, **non-circular**. Used by
   **both** the global **and** the gDNA messages (one gDNA between-node-spread curve, one currency). The
   **RNA** `σ²_bio` stays the existing belief-fit for now (the RNA observable is gDNA-contaminated — the hard
   channel; deferred, OQ4). C++: none.

---

## 5. The `strand_obs` fork dissolution (the §3.3 refinement — fixes Bug C)

Today (`simplex_sweep.py:66-77`): `strand_obs ⇒ Jeffreys(½)` XOR `~strand_obs ⇒ global`. In one currency there
is **no fork** — every node gets **all** of: strand likelihood + global `_binom_pseudo` + the Jeffreys-½
baseline. Deference is emergent (strand wins where decisive; global wins where strand is flat). This is what
brings **Bug C** (the 133K single-strand-exon Jeffreys vertex-push in unstranded-0) into scope: a flat-κ
single-strand exon now gets the global's `→0` zero-pin instead of the vertex-pushing Jeffreys alone.

**Risk to validate:** applying the global to *informative* single-strand nodes (κ far from ½) must not drag
their strand-determined `f_g` toward the genome-average. The global must be weak there (it is, under capture:
large spread ⇒ small `N_global`) so the strand dominates. **Gate this explicitly** (a stranded single-strand
node with real gDNA must keep its strand answer).

---

## 6. Open questions / problems to resolve BEFORE implementing

- **OQ1 — RESOLVED (derived, §4).** Global = NB-count log-prior; zero-pin/deference via `ρ_global`,`σ²_g`;
  Jacobian-free; one constant `a=1`.
- **OQ1b — RESOLVED (seed set, §4.3, user-specified):** structural seeds (intergenic + intron regions +
  intergenic-exon + intron-exon boundaries) always; single-strand nodes (strand-deconv gDNA) for stranded data,
  confidence-weighted to fade at κ→½. Remaining **implementation details** to confirm during build (not
  blocking the design):
  - (i) **boundary gDNA density** for a seed = crossing mass / boundary eff-len; use the **exon-facing side**
    (the high-μ datum) vs the total crossing — confirm which.
  - (ii) **strand-seed weight** = strand-deconv confidence `∝ n·(2κ−1)²` (the strand Fisher info) — confirm the
    exact expression feeds the `var~mean` IRLS weight.
  - (iii) **nascent-sparse assumption** for intron + boundary seeds — documented; the one stated assumption.
  - (iv) **seed-seed edge sufficiency / sparse fallback** — the `var~mean` fits on adjacent seed pairs; confirm
    there are enough across the μ range (esp. high-μ via boundaries) and the `_constant_offset` fallback covers
    sparse chains.
- **OQ2 — expected direction.** With the seed-based `ρ_global`/`σ²_g`, the intent is to **hold** the unstranded
  zero line (pin fires on honest `≈0`) and **partially recover** capture-unstranded via the boundary seeds.
  Residual = the exon-edge-vs-interior gap under capture-unstranded (smaller than intron-only). Confirm we
  accept "measure it" rather than require a guaranteed capture number.
- **OQ3 — wall/median numerics.** With large `N`, the `α·log(_PRIOR_EPS)` edge term is a big constant at the
  `f=0/1` lattice points. It should be *correct* (a confident message rightly kills the far wall), but verify
  it doesn't distort the `_fg_median` readout near the walls (re-run the +8.7-pt median-vs-mean case; add a
  direct bimodal/wall unit test).
- **OQ4 — σ²_bio's larger role.** After the change, `N_src` depends on σ²_bio *directly* (no Jacobian to swamp
  it). We're keeping the existing belief-derived σ²_bio (Step 1 dropped). Confirm we accept that `N_src` is "as
  honest as the converged σ²_bio" — and that we measure, not assume, it's adequate.
- **OQ5 — state migration (scope).** I propose **NOT** migrating `NodeBelief` to `(n_c, τ_c)` in Step 2 — it's
  not needed to remove the Jacobian (the GS sweep keeps fraction state; only the prior *forms* change). It's
  the FB-relay door, deferrable. Confirm we keep state as fractions for Step 2.
- **OQ6 — the global's two var~means.** `gprior_vm` (its σ²_prior, fit on raw `M/E` index, `bp_solver.py:703`)
  is reused unchanged for `N_global`. Confirm.

---

## 7. Exact edit sites (once OQ1/OQ2/OQ5 are locked)

| file | function | change |
|---|---|---|
| `bp_solver.py` | `_message` (:559) | return `(μ_c, N_src)`; delete J1 (:575) + J2 (:581-582); drop `mass_src` arg |
| `bp_solver.py` | `node_sweep` global block (:632-710) | delete `geom2_global` (:636); compute `N_global` (OQ1) not `tau_global`; pass `(μ_global, N_global)` |
| `bp_solver.py` | `node_sweep` `_solve` (:656) + call sites (:727,735,747,757,766) | rename precision→count args; `_message` now returns `N_src`; **pass global to ALL solvable nodes** (drop the `~strand_obs` restriction → §5) |
| `simplex_sweep.py` | `_local_loglik` (:30) | add `_binom_pseudo`; replace the three Gaussian terms (:76,84,95) with `_binom_pseudo`; apply global + Jeffreys-½ to **all** nodes (dissolve the `strand_obs` fork); rename `*_precision`/`global_tau` → `*_count`/`global_count` |
| `simplex_sweep.py` | `_solve_nodes` (:129) | thread the renamed count args; readout (`_fg_median`/`_axis_mean`/`_fg_var`) **unchanged** |
| `simplex.py` / strand / spliced / eff-len | — | **untouched** |
| C++ | — | **none** |

No `NodeBelief` change (OQ5). No σ²_bio change (OQ4). No new constant if OQ1 resolves to a derived `N_global`.

---

## 8. Validation (the gates — all measured, judged as one landing per the user)

- **Portable pytest:** `factor-1-under-uniform` (`test_gdna_sweep_factor1_uniform` + priors/capture variants) —
  **must hold**; a new **two-sided unit test** (a `μ=0.2` message ⇒ posterior `f_g` mode ≈ 0.2, not 0.77);
  a **deference unit test** (`N_src ≪ N_dst` ⇒ a decisive strand survives the message); a **wall/median** test
  (OQ3); full suite green.
- **Local gate** (`honest_precision_baseline.py`, recorded Step-0 baseline): per-condition signed phantom/siphon.
  Expected: `stranded-0` AMBIG-exon stays ≈0; `capture` net no worse (watch the AMBIG over-call vs a new
  siphon); `unstranded-0`/`capture-unstranded` move per OQ2 (improve if OQ1=(b) holds the zero-pin; otherwise
  regress as the exposed floor — **documented, not hidden**).
- **Secondary:** post-EM `net_gdna_to_rna` via `evaluate_suite.py` (does it propagate through quant).
- **Convergence:** `sweep_deltas` still settle within `sweep_max_passes` (the prior-form change must not
  destabilise the fixed point).
