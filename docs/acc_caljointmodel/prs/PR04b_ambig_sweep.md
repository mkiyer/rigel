# PR 4b — AMBIG boundary↔region sweep imputation (D7 leg 2)

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §4 D7, §5.
**Builds on:** PR 4a ([`PR04a_estep_exposure.md`](PR04a_estep_exposure.md)).
**Type:** Python-only. **Build required:** no.
**Status:** **DESIGN FINAL (rev. 3) — awaiting go-ahead to implement.** Theory
corrected per your note: gDNA effective length comes from the gDNA FL
distribution, *unconstrained* by region size; region size constrains only the
fractional mass. Decisions III′.1–.7 + scope (single PR, §II.6) resolved.

> **Record fix.** The PR04a doc (III.3) and master-plan §5 say the recovered
> `boundary_sweep.py` has a `confidence_floor = 0.6` cliff. The `fc96902` source
> has none — its transfer weight is already smooth. I'll amend those two lines if
> you want the record clean on `main`.

---

# Part I — Theory

## I.1 gDNA sampling and the two effective lengths

gDNA fragments are sheared from the **entire chromosome**: left-endpoints are a
Poisson process at rate `ρ_0·ω(x)` per bp (`ω` = local accessibility/exposure),
lengths drawn from the gDNA fragment-length pmf `f` — **independent of any region
boundary**. Let `μ_FL = Σ_ℓ ℓ·f(ℓ)`.

An element's gDNA **count exposure** is the genomic measure of fragment
start-positions that produce the event in question:

- **Region** `R` of physical length `L` — a fragment is *contained* iff it fits,
  so it needs `ℓ ≤ L` and has `L−ℓ` admissible starts:
  $$ L_{\text{eff}}(R) \;=\; \mathbb{E}_f[\max(0,\,L-\ell)] \;=\; \sum_{\ell \le L}(L-\ell)\,f(\ell). $$
  FL-corrected and region-limited (contained fragments are shorter than the
  region — your III′.6). `→ L` for `L ≫ μ_FL`; small for `L ≲ μ_FL`.

- **Boundary** `B` — a fragment *crosses* a point iff its start lies in the `ℓ`
  bp upstream of it, **regardless of region sizes** (a long fragment simply
  spills onto further regions and still crosses `B`):
  $$ \boxed{\;L_{\text{eff}}(B) \;=\; \mathbb{E}_f[\ell] \;=\; \mu_{\text{FL}}.\;} $$
  Just the mean gDNA fragment length. No region-size term.

> **Mass ≠ effective length (the disentanglement).** The fractional accumulator
> *splits* a crossing fragment's mass across the regions it overlaps, in
> proportion to **overlap in bases**, and per-side overlap **is** capped by the
> region size. That region-limited *mass* is the D1 quantity PR 4a attributes to
> regions; it is **not** the count exposure. The molecule was still drawn from
> the unconstrained gDNA FL, so the crossing exposure stays `μ_FL`. Keeping these
> separate is what was previously obfuscated.

> **Fair to RNA and gDNA — same geometry, different FL.** The contained/crossing
> formulas are **pure geometry**; they apply identically to RNA with the RNA FL
> pmf, so `effective_length` is FL-agnostic. The model *consumes* a **gDNA**
> effective length because gDNA is parameterised as **density × exposure**
> (`ρ_0, ω`) and must predict gDNA counts. RNA is the **observed residual** mass
> (doc 01 §9 — no per-region RNA exposure factor), so its geometric constraints
> are already carried by the observed counts, not a modelled length. This
> self-corrects fairly: in a region shorter than the gDNA FL, `L_eff_gdna(R)` is
> small (gDNA can't be contained → routed to its boundary crossings) and the
> contained count correctly falls to the shorter RNA. RNA's truncation by its
> (unknown) source transcript is absorbed into the empirically observed RNA FL
> (spliced fragments). **Full-symmetry flag (PR 5):** the only place an explicit
> `E_{f_rna}` mean would enter is the *live* count channel — inert on the single
> pass, so deferred to PR 5 rather than introducing an RNA rate now (§III).

Both effective lengths use the **same** gDNA `f`, so every node's gDNA density
estimate `(gDNA count)/L_eff ≈ ρ_0·ω` is on **one comparable scale** — the
prerequisite for sweeping densities between regions and boundaries.

## I.2 The alternating region↔boundary chain

A reference is `… B → R → B → R …` (`k` regions, `k+1` boundaries; ends are
terminals). **Every node carries a gDNA Gamma pair** `(α, β) = (gDNA count,
L_eff)` and estimates the local density `ρ_0ω`:

| Node | α (gDNA count) | β (effective length) |
|---|---|---|
| region `R` | `M_g_contained` = `π_g(R)·n_u^{contained}` (E-step) | `E_f[max(0, L−ℓ)]` |
| boundary `B` | gDNA crossing flux `π_g(B)·u_B` (E-step on the crossing flux) | `μ_FL` |

The asymmetry you flagged drives the sweep:

| Boundary (by flanking signature) | Unspliced crossing | Role |
|---|---|---|
| **exon→intron / exon→intergenic** | contiguous **gDNA** — mature RNA cannot cross | **clean gDNA anchor** for its regions |
| intron/intergenic (no exon) | gDNA (+ nascent) | clean-ish |
| **exon→exon** | gDNA **and** contiguous mature/nascent RNA | **convoluted mixture**; itself needs imputing |

**Worked example (hybrid capture).** `intron region — B:intron→exon — exon
region`. The intron region is **depleted** (uncaptured); the `intron→exon`
boundary carries **enriched, RNA-clean spillover gDNA**; the exon region is a
gDNA+RNA mixture. The clean boundary's density **imputes the exon region**
(boundary → region). Symmetrically a deconvolved region imputes an ambiguous
`exon→exon` boundary (region → boundary). Hence the alternating sweep.

## I.3 Purity (transfer) weight `w_B ∈ [0,1]`

Evidence crossing boundary `B` is attenuated by how reliably `B` conducts a
*clean gDNA* density signal:

```
traffic_B   = u_B / (u_B + 1)          # more crossing fragments ⇒ more reliable  (unit pc)
rna_clean_B = 1 / (1 + s_B)            # spliced crossings s_B = mature-RNA junction (unit pc)
w_B = traffic_B · rna_clean_B · π_g(B)         # π_g(B) ∈ (0,1): the boundary's own gDNA purity
w_B = 0 at reference terminals
```

`π_g(B)` (the boundary's E-step gDNA fraction) is the **smooth** exon-exon
down-weight (III′.2): at an exon→exon mixture the count/flux give `π_g(B) < 1`,
choking the conduit, while a clean exon→intron boundary has `π_g(B) ≈ 1`. No
binary type gate, no `confidence_floor`. Two unit pseudocounts; no other
constants.

## I.4 The sweep — forward / reverse decay scan

Two linear passes per reference over the alternating node sequence. The running
`(α,β)` is decayed by `w_B` when crossing each internal boundary (regions are
contiguous → no decay), and each node injects its own local evidence:

```
# forward (left → right):
run_a = run_b = 0
for node in sequence [R0, B1, R1, B2, …]:
    if node is internal boundary B:  run_a *= w_B;  run_b *= w_B    # decay across the conduit
    from_left_a[node], from_left_b[node] = run_a, run_b             # inflow from the left
    run_a += alpha_local[node];  run_b += beta_local[node]          # inject own evidence
# reverse (right → left): symmetric → from_right_a/b
alpha_swept[node] = alpha_local + from_left_a + from_right_a
beta_swept[node]  = beta_local  + from_left_b + from_right_b
```

Evidence decays as the **product of the `w_B` crossed**, so the **immediately
adjacent** nodes dominate (III′.3, nearest-neighbour, no extra constant): an
AMBIG exon's nearest evidence is its two flanking boundaries (one `w` each) — the
clean `intron→exon` conduits — then the regions beyond (two `w`'s), etc. A choked
boundary smoothly cuts the chain; past it, the region reverts to the global
density (leg 3).

## I.5 Imputed AMBIG exposure (Gamma posterior; AMBIG-only)

The swept `(α, β)` feeds the **same** closed-form posterior as PR 4a — it is just
a neighbour-pooled `M_g` / `L_eff`:

```
for region r with ts_class[r] == TS_AMBIG:
    ω[r]            = (1/φ + alpha_swept[r]) / (1/φ + ρ_0 · beta_swept[r])
    log_omega_var[r]= 1 / (1/φ + alpha_swept[r])
# decodable regions (POS/NEG/NONE): keep PR 4a exposure untouched (III′.4)
```

AMBIG regions contribute their own count-only local evidence (strand-blind —
they are a both-strand-RNA + gDNA mixture; III′.5); NONE/intergenic are strong
clean anchors (III′.5). Where the chain is choked, `α_swept,β_swept → 0 ⇒ ω → 1`
(global fallback).

---

# Part II — Implementation plan

### II.1 `effective_length.py` (NEW) — gDNA FL-corrected exposures

```python
def region_eff_length(L_bp: np.ndarray, fl_pmf: np.ndarray) -> np.ndarray:
    """E_f[max(0, L−ℓ)] per region; vectorised cumulative sum over the FL pmf."""
def boundary_eff_length(fl_pmf: np.ndarray) -> float:
    """μ_FL = Σ ℓ·f(ℓ): the gDNA crossing exposure (region-independent)."""
```

Both are **FL-agnostic** (any `fl_pmf`) — the fairness point: 4b calls them with
the **gDNA** pmf (post-scan gDNA/intergenic histogram via `.pmf()`, III′.7); the
identical functions serve RNA if PR 5's count channel adopts an `E_{f_rna}` mean.
`region_eff_length` is a single `O(R + FL_BINS)` pass (precompute
`Σ_{ℓ≤L}(L−ℓ)f(ℓ)` from cumulative sums of `f(ℓ)` and `ℓ f(ℓ)`). Pure, no tunables.

### II.2 `sweep.py` (NEW) — purity weight + alternating scan + AMBIG exposure

```python
def boundary_purity_weight(payload, region_arrays, pi_g_boundary) -> np.ndarray   # w_B, len = n_internal_boundaries
def sweep_ambig_exposure(substrate, region_arrays, *, alloc_region, alloc_boundary,
                         rho_0, phi, mu_fl, region_eff_len) -> (omega, log_omega_var)
```

`alloc_region`/`alloc_boundary` are PR 4a E-step `Allocation`s (region nodes use
contained mass; boundary nodes use the crossing **flux** — §I.2). The scan is the
two loops above over the per-reference node order from `region_arrays`; it writes
refined `(ω, log_omega_var)` **only** at AMBIG rows.

### II.3 `calibrate.py` (EDIT) — wire FL + replace AMBIG exposure

1. New arg `frag_length_models`; derive `fl_pmf`, `mu_fl`, and
   `region_eff_len = region_eff_length(L_bp, fl_pmf)`.
2. Use `region_eff_len` (not physical length) wherever `L_eff` enters — the
   `ρ_0` seed (`density.py`) and the exposure denominator. **This refines every
   region's exposure** (decodable included), so it is a model upgrade, not just
   an AMBIG change (§II.6).
3. After `exposure_posterior(...)`, call `sweep_ambig_exposure(...)` and overwrite
   the AMBIG rows of `(ω, log_omega_var)`.

### II.4 Constants (Q6)

`_TRAFFIC_PSEUDOCOUNT = 1.0`, `_SPLICED_PSEUDOCOUNT = 1.0` (the two unit
pseudocounts in `w_B`). No locality `λ` (III′.3). No floors beyond doc 03 §8.

### II.5 Tests

- `test_effective_length.py`: `boundary_eff_length == μ_FL`; `region_eff_length`
  vs hand sums incl. `L < μ_FL` (small exposure) and `L ≫ μ_FL` (`≈ L`); your
  300 bp-FL / 400 bp / 200 bp case — boundary exposure `= μ_FL`, **independent**
  of the 200 bp region (the explicit regression against the old tent error).
- `test_sweep.py`: `w_B` high for `intron→exon`, low for `exon→exon` / spliced /
  terminal; an AMBIG exon flanked by a clean boundary + gDNA-rich intron inherits
  that density; behind a choked boundary → `ω ≈ 1`; decodable rows byte-identical
  to PR 4a; reuse `fc96902` `test_boundary_sweep.py` for the scan arithmetic.
- Update PR 4a calibration tests for FL-corrected `L_eff` (they currently assume
  physical length); add a synthetic gDNA FL pmf to the fixtures.

### II.6 Scope — single PR (resolved)

The FL-corrected `L_eff` (II.3.2) changes **all** regions' exposure and the `ρ_0`
seed — a coherent accuracy upgrade bundled with the AMBIG sweep as **one PR 4b**
(your decision). Implementation order within the PR: (1) `effective_length.py` +
wire FL-corrected `L_eff` into the seed/exposure and update PR 4a tests; (2)
`sweep.py` + AMBIG exposure overwrite. Each step keeps the suite green before the
next.

---

# Part III — Decisions (resolved)

III′.1 region `(M_g,L_eff)` + boundary `(flux, μ_FL)` nodes ✓ · III′.2 smooth
purity via `π_g(B)`, no type gate ✓ · III′.3 `w`-product decay, **no** new
constant ✓ · III′.4 AMBIG-only ✓ · III′.5 NONE anchors, AMBIG contributes
count-only ✓ · III′.6 FL-correct region `L_eff = E_f[max(0,L−ℓ)]` ✓ · III′.7 gDNA
FL — **INTERIM**: the pipeline passes the UNSPLICED-category histogram, which is
gDNA **+ nascent RNA** (an over-broad proxy, not the true gDNA FL). The proper
derivation (gDNA-dominated regions/boundaries → mixture-EM subtracting RNA → EB
shrink) was burnt with SRD-v1 and is **resurrected in PR 4c**; only the AMBIG
sweep's effective lengths depend on it. **Scope:** single PR (effective-length
upgrade + sweep), §II.6 ✓.

**Fairness (RNA vs gDNA):** the effective-length geometry is FL-agnostic and
applies to both (same formulas, different FL); the model consumes a gDNA
effective length only because gDNA is density × exposure while RNA is the
observed residual (doc 01 §9). **PR 5 flag:** decide whether the live count
channel expresses the RNA mean via `E_{f_rna}` (full symmetry) — inert on the
single pass, so it does not affect PR 4b.

All decisions resolved; design is final.

## Rollback

Revert the sub-PR. `calibrate()` reverts to the PR 4a single-pass result
(physical `L_eff`, no sweep). No on-disk artifacts change.
