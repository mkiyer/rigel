# Phase 4 — Derive ρ₀ & Exposure: Implementation Plan (v3, all clarifications resolved)

All §7 questions are resolved (recorded in §0/§3/§6 as math, not prose, so they can't drift).
**Two resolutions re-open earlier work** — the boundary unit (Phase 3) and the boundary
effective length (Phase 1) — flagged in §7 as residual.

## 0. Resolved decisions (math, not prose)

- **A — ρ₀ = global average gDNA density;** exposure relative to it (`ω<1` depleted, `=1`
  uniform, `>1` enriched).
- **B — start simple;** no robustness/shrinkage machinery up front (Phase 5 adds it from
  benchmarks).
- **C — locus = its nodes** (regions + boundary-sides, incl. START/END); per-locus attributes
  use all of them (§6).
- **D — refactor, no backward compatibility;** redesign the result type + `assemble_priors`.
- **E — exposure is the pure derived ratio; no separate shrinkage** (the mass is already
  regularized by Phase 3's count prior — shrinking ω too would double-count).
- **F — graceful across the whole zero-gDNA ↔ zero-RNA spectrum.**
- **G — boundary length is per-side and has two distinct meanings** (§3).
- **i/ii — the exposure↔length relationship (settled):**
  ```
  gdna_mass_node = ρ₀ · ω_node · L_node            (exposure is a MASS factor; this is what we keep)
  exposure-weighted eff_len_node = ω_node · L_node  (ω in the numerator)
  abundance = mass / eff_len = ρ₀                   (enrichment carried consistently, cancels to the rate)
  ```
  The component effective length is the **exposure-weighted sum** `Σ_node ω_node · L_node`
  (the per-read normalizer). `length / exposure` was the earlier mistake.
- **v — transcript components reuse the gDNA-derived ω** for their effective lengths. gDNA is
  the *only* viable exposure estimator (its baseline is uniformity; RNA expression is not).
- **vi — keep ρ₀/ω derivation and node→locus assembly as clean re-runnable functions** so a
  future iterative-calibration outer loop wraps them without restructuring.

## 1. ρ₀ (global density) — §1 unchanged: `ρ₀ = Σ_nodes gdna_mass / Σ_nodes L_eff(density)`.

## 2. Per-node exposure ω — pure ratio (decision E): `ω_node = (gdna_mass/L_eff) / ρ₀`. No shrink.

## 3. Boundary length has TWO meanings (decisions G + iv) — both per-side

A boundary is one position with a **left** and a **right** side; each side belongs to one
region (its `left`/`right` substrate view). The two lengths are **not the same quantity**:

1. **Density effective length (for mass→density):** the length the fractional crossing mass
   is divided by. It is region-constrained and is an **integral over the FL distribution**,
   not a single length:
   ```
   L_eff_side  =  E_FL[ min(ℓ, region_side_length) ]        (∫ f(ℓ)·min(ℓ, R_side) dℓ)
   ```
   A long fragment over a short side contributes only `R_side` of length; over a long side,
   up to `ℓ`. This refines the shipped `μ_FL` (§7-ii residual → patch Phase 1).
2. **Statistical power (for the count prior's precision):** the **discrete** crossing-read
   count, drawn from the **full** FL distribution (long fragments cross many regions and are
   *not* region-constrained). This is the `count_evidence` (a fragment count), kept separate
   from the density length per decision iii/§0-i-ii.

## 4. Per-node granularity = boundary-SIDE (decision iii) — **revises shipped Phase 3**

The deconvolution unit is the **boundary-side**, not the whole boundary. Iterating over
regions, each region aggregates three nodes: its **contained** mass, the **right side of its
left boundary** (`substrate.left[r]`), and the **left side of its right boundary**
(`substrate.right[r]`). So the joint decode runs on the substrate's three views per region,
each with its own per-side density length (§3.1) and discrete power (§3.2). **Phase 3's
`decode_boundaries` (which combined the two sides into one node) must be revised to decode
per side** (§7 residual).

## 5. Every-node guarantee + zero spectrum (decision F): finite `(gdna, rna, L_eff, ω)` for
every node; zero-gDNA ⇒ ρ₀=0, ω=1; zero-RNA ⇒ ω≈1; no divide-by-zero.

## 6. Per-locus assembly — specified, built in Phase 6 (decisions C, D, i/ii, v)

- **Locus** = its regions + its boundary-sides (incl. START/END).
- **Locus / component gDNA mass** = `Σ_node gdna_mass` over the component's nodes (gDNA
  component = all locus nodes; transcript component = nodes overlapping its exons).
- **Component effective length** = `Σ_node ω_node · L_eff_node` (the §0-i/ii exposure-weighted
  sum); transcripts reuse the gDNA-derived ω (decision v).
- Refactored result type + `assemble_priors` (decision D).

## 7. Residual issues surfaced by these resolutions

**R1 — decision iii re-opens shipped Phase 3.** The per-boundary-side unit means
`joint_decode.decode_boundaries` must be rewritten to decode each side independently (three
per-region views), each with its §3 per-side length. Plan: a small **Phase 3.1 revision**
(per-side boundary decode) *before* Phase 4, re-validated against the same scripts. (Phase 3's
core — regions, the joint posterior — is unchanged.)

**R2 — decision G/iv re-opens shipped Phase 1.** The density length becomes
`E_FL[min(ℓ, R_side)]` per side, replacing `boundary_eff_length = μ_FL`. This changes the
Phase-1 density numbers. Patch `effective_length.py` (a new `boundary_side_eff_length(fl_pmf,
region_side_len)`) and re-run the Phase-1 validation. I'll do R1 and R2 together since both
touch the boundary representation.

**R3 — the FL-integral effective length needs the exact form.** `E_FL[min(ℓ, R)]` is my read
of "an integral over the compatible FLs, constrained by the region side." Confirm that's the
intended functional (vs e.g. `E_FL[max(0, R − ℓ)]`, the contained-style form, or a
crossing-specific kernel). The two differ for short regions and set the density scale.

**R4 — residual wording vs math on i.** I've committed to `eff_len = ω·L` (so an enriched node
has a *larger* exposure-weighted length, and abundance cancels to ρ₀). Your verbal "larger
exposure → smaller effective length" would instead give abundance `= ρ₀·ω` (a single ω, in the
length). These are different downstream scales. I went with the **math you confirmed**
(`ω·L`); flag if the single-ω/`ρ₀·ω` reading was intended.

**R5 — two lengths, one count prior.** §3 gives a node a density length (region-constrained
FL-integral) *and* a discrete power (full-FL count). In the joint count prior they feed
different slots: the density length sets `π_count = density·L_eff/M`; the discrete count sets
the prior concentration `κ_c`. Confirm this split is what "stay consistent everywhere" means
(density uses constrained length; precision uses unconstrained count).

## 8. Sequenced plan

1. **Phase 3.1** (R1+R2): per-boundary-side decode + per-side density length; re-validate.
2. **Phase 4**: `derive.py` — ρ₀ (§1), per-node ω (§2), per-node lengths (§3); validate
   zero/capture; report `min ω` and component eff-lens (informs Phase 5).
3. **Phase 5**: measure exposure stability; add robustness only if warranted (decision B).
4. **Phase 6**: node→locus assembly (§6), refactored result + `assemble_priors`, wire in,
   tear down `update_rho_0`/`exposure.py`, regenerate goldens.

## 9. Implementation surface

- `effective_length.py`: `boundary_side_eff_length(fl_pmf, region_side_len)` (R2/R3).
- `joint_decode.py`: per-side boundary decode (R1).
- new `derive.py`: ρ₀ + per-node ω + per-node lengths (Phase 4).
- `scripts/debug/proto_derive.py`: validation.
- Phase 6: result refactor + `assemble_priors` + teardown.
