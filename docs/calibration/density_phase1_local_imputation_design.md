# Density Phase 1 — local boundary-anchored count density (design + dry-run)

**Status:** design / dry-run, ready to implement next turn. 2026-06-09. Implements Phase 1 of
`density_sweep_implementation_plan.md`, now that the accumulator span redesign makes the
boundary-crossing flux an **unbiased** density estimator (no `L_cross`, no correction factor).

> **One-line goal.** Replace the global region↔boundary **sweep** in
> `density_model.node_gdna_density` (which collapses every node to the genome-wide average density —
> the root cause of the captured-exon count-prior bug and the gDNA→RNA leak) with a **local**
> estimator: an observable region uses its own contained gDNA density; a non-observable (exon/AMBIG)
> region is anchored from its **observable boundary sides** (`crossing gDNA count / fl_mean`), with
> run interiors filled by a directional inward carry.

---

## 1. What is broken, precisely

`node_gdna_density` runs an alternating region↔boundary sweep with conduit weight
`w = flux/(flux+1) ≈ 1` at real depth → **undamped** → `density[node] ≈ Σ observable gDNA / Σ observable
length` = the **global average** for every node (measured ~0.62/bp uniformly on GENE0037). So a
captured exon (non-count-observable) inherits the global density instead of its own enriched gDNA
density. In `deconv_regions` that yields `count_gdna_frac = density·eff_len/mass ≈ 0.02` with a
fabricated concentration `count_evidence = density·eff_len ≈ 700` that **overrules the correct strand
clue** (gf≈0.8) → exon gDNA under-deconvolved ~2× → low per-locus prior → EM leak.

## 2. Data flow — where `node_density` is consumed (the contract to preserve)

`node_gdna_density → NodeDensity{density, gdna_mass, count_evidence, region_count_observable,
boundary_count_observable, …}`. Consumers (verified):

| Field | Consumer | Use |
|---|---|---|
| `density[r]` | `deconv_regions` | region count-prior **mean** `clip(density·region_eff_len/mass,0,1)` |
| `density[r]` | `deconv_sides` (`region_density`) | **fallback** density for region `r`'s count-**non**-observable boundary sides |
| `count_evidence[r]` | `deconv_regions` | region count-prior **concentration** |
| `count_evidence[r]` | `gdna_strand._region_seeds` | seed **weight** `clip(count_evidence/mass,0,1)` — **only on `region_count_observable` regions** |
| `region_count_observable`, `boundary_count_observable` | `gdna_strand`, `deconv_sides` | seed / fallback masks |
| `gdna_mass[r]` | — | **unused** (derive uses the *deconv* gdna_mass); keep `= density·eff_len` for struct stability |

Two consequences that shape the design:
- A count-observable region's `density`/`count_evidence` already = its **own contained** quantity in
  the new design (the sweep only contaminated it). So the gDNA strand fit's seeds get *cleaner*, not
  broken. The imputed **exon** density is consumed **only** by `deconv_regions` (the leak fix) and as
  the side fallback — never by the strand fit (seeds gate on `region_count_observable`).
- `deconv_sides`/`_compute_side` **already** computes each count-observable side's *own* crossing
  density locally; `node_density.density` is only its fallback for non-observable sides. So Phase 1
  changes the *region* density (and that fallback), not the side-deconv's own-density path.

## 3. The new algorithm (dry-run pseudocode)

```
node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean, rna_sense_frac):
    ts   = strand_class ; ref = ref_id ; R = n_regions ; EPS = 1e-9
    reg_obs, bnd_obs = count_observable_masks(signature, ref)

    # --- strand-cleaned gDNA COUNTS (closed form, shared helper) ---
    # contained (region) and the two boundary sides, each oriented by region r's strand.
    cont = clean_count(contained, ts, rna_sense_frac)      # = clean_gf * (pos+neg)
    lft  = clean_count(left,      ts, rna_sense_frac)
    rgt  = clean_count(right,     ts, rna_sense_frac)
    # clean_gf: (sense_frac − κ)/(½ − κ) clipped [0,1]; 1.0 where κ≈½ (unstranded),
    #           ts==NONE (intergenic, pure gDNA), or ts==AMBIG (no sense → defer, Phase 4).

    # --- per-side boundary observability for region r ---
    left_anchor[r]  = bnd_obs[r-1] and ref[r]==ref[r-1]    # boundary (r-1,r) observable
    right_anchor[r] = bnd_obs[r]   and ref[r]==ref[r+1]    # boundary (r,  r+1) observable

    density = full(R, NaN)
    for r in range(R):
        if reg_obs[r] and region_eff_len[r] > EPS:
            density[r] = cont[r] / region_eff_len[r]                 # own contained density
        else:
            est = []
            if left_anchor[r]:  est.append(lft[r] / fl_mean)         # crossing gDNA → density
            if right_anchor[r]: est.append(rgt[r] / fl_mean)
            if est: density[r] = mean(est)                           # else NaN → run-fill

    # --- run-fill: NaN interiors inherit the nearest anchored neighbour, per ref ---
    # forward carry (L→R) and reverse carry (R→L) within a ref; average where both reach.
    density = directional_fill(density, ref)

    # --- fallback: a region the fill never reached (no observable node in its ref) ---
    density = where(isnan(density), 0.0, density)   # → count_gf≈0, conc≈0 ⇒ Beta(½,½) ⇒ defer to strand

    gdna_mass      = density * region_eff_len
    count_evidence = gdna_mass                       # Phase-1 (Phase 3 revises to observed crossing count)
    return NodeDensity(density, gdna_mass, count_evidence, reg_obs, bnd_obs, …)
```

`directional_fill`: forward pass copies the last non-NaN value across NaNs within a ref; reverse pass
likewise; a NaN reached from both sides takes the mean (the prototype's run-fill; demonstrated on the
antisense `ex+|ex+|ex−|ex−` interior).

## 4. Why the leak is fixed in Phase 1 (before Phase 3's concentration change)

Post-accumulator-fix, a captured exon's flanking intron|exon boundaries are count-observable, and
their crossing gDNA fragments **overlap the exon ⇒ are captured/enriched**, so `side_flux/fl_mean ≈
the exon's enriched gDNA density` (dry-run: GENE0037 ~21 ≈ truth). Then in `deconv_regions`:
`count_gdna_frac = density·eff_len/mass ≈ (enriched gDNA expected)/(total contained) ≈ true exon gDNA
fraction (~0.8)`, with concentration `= expected gDNA count` (a *sensible* confidence). The count
prior now **agrees with** the strand clue instead of fighting it. (Phase 3 sharpens the concentration
to the directly-observed crossing count; Phase 1's `density·eff_len` is already correct-in-mean.)

## 5. Strand-cleaning consolidation (DRY)

The closed form `(sense_frac − κ)/(½ − κ)` clipped to [0,1], with the `κ≈½ / NONE / AMBIG → 1.0`
guards, currently lives in **three** places: `calibrate.py` (region `gdna_frac`), `joint_deconv.
_compute_side` (side `clean_gdna_frac`), and the prototype. Phase 1 factors **one** helper —
`density_model._clean_gdna_frac(sense, total, rna_sense_frac)` (returns the per-element fraction; the
NONE/AMBIG masks applied by the caller using `strand_class`) — and uses it in `node_gdna_density`
(region + both sides). It also:
- **removes** the `gdna_frac` block from `calibrate.py` (its only consumer was `node_gdna_density`,
  which now computes it internally from `rna_sense_frac`); and
- **switches** `_compute_side` to the same helper (formula-identical ⇒ the side-deconv stays
  bit-for-bit, verified on golden), so there is exactly one implementation.

New `node_gdna_density` signature: drop `gdna_frac`, add `rna_sense_frac` →
`node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean, rna_sense_frac)`.

## 6. Dry-run: issues surfaced

1. **Degenerate tiny regions (`region_eff_len ≈ 0`).** A region shorter than the gDNA FL holds no
   contained fragment → `count/eff_len` blows up. Resolved by the `reg_obs[r] and eff_len>EPS`
   guard: such regions fall through to **boundary anchoring** (crossing flux/`fl_mean` is
   region-size-free) and then run-fill. (The prototype's `reg_eff>1.0` guard, generalized.)
2. **`count_evidence` is still `density·eff_len` in Phase 1.** For an imputed exon this is an
   *expected* (not directly observed) count → a possibly-large concentration. It is **correct in
   mean** (§4) so it fixes the leak, but Phase 3 must replace it with the supporting crossing count to
   avoid over-confident priors on thinly-anchored exons. Flagged, not fixed here.
3. **No-anchor fallback = density 0 (⇒ defer to strand), never global.** A region with no observable
   node anywhere in its reference (e.g. a single-exon-only contig) stays NaN → set to 0 → count prior
   Beta(½,½) (Jeffreys) → strand governs. This is the Phase-4 "near-zero concentration" tier; Phase 4
   later upgrades tier-1/2 (neighbour / intron-baseline). Acceptable Phase-1 behaviour.
4. **AMBIG sides/regions.** No defined sense ⇒ `clean_gf = 1.0` (treat crossing/contained unspliced as
   all gDNA). This over-calls gDNA on antisense overlaps at low SS — the genuine floor; Phase 4 owns
   the explicit AMBIG policy. Matches `_compute_side`'s current AMBIG handling, so consistent.
5. **Per-side orientation across antisense seams.** Each side is cleaned by **region `r`'s** strand
   (the region we are anchoring), exactly as `_compute_side` orients its own side — consistent. A
   side whose boundary's *other* region differs in strand is still oriented by `r` (we want `r`'s
   density). Confirmed against the prototype.
6. **Run-fill direction & conflicts.** Forward+reverse carry with mean-on-overlap handles runs of
   length ≥1; for a run reachable from only one side, that side's value carries (no global leak). Unit
   tests must cover a 1-region run, a 3-region antisense run, and a one-sided run.
7. **gDNA strand fit unaffected in kind.** Its seeds gate on `region_count_observable`; their
   `count_evidence` becomes the *own contained* count (cleaner than the swept value). Improvement, not
   regression — but it *will* shift the fitted overdispersion slightly (golden churn, expected).
8. **`deconv_sides` fallback.** Non-count-observable sides borrow `node_density.density` of region
   `r`, now the local exon estimate instead of the global average — strictly more correct.

## 7. Tests & acceptance

- **Unit (`density_model`):** a single exon between two introns recovers the known uniform density
  from either/both sides; a tiny (<FL) observable region anchors from its boundaries, not contained;
  a 3-region antisense run fills the interior; a no-anchor region → density 0.
- **Helper:** `_clean_gdna_frac` table (sense=κ→0; symmetric→1; clipping; κ≈½→1).
- **Bit-for-bit:** the `_compute_side` refactor leaves the side-deconv unchanged on golden.
- **Worked (oracle):** GENE0037 exon density 0.62 → ~20; `count_gdna_frac` ~0.8; joint `gdna_frac`
  ~0.8 (a debug-script check, not necessarily a unit test).
- **Suite + golden:** gDNA/combo + nRNA goldens regenerate (calibration density changes); full suite
  green. Then the benchmark (`evaluate_suite.py` / `bench_calibration`) shows the captured-exon leak
  drop with no rise in `gdna_none` false gDNA.

## 8. Risks

- (a) Golden churn across nearly all scenarios (calibration density feeds every gDNA/nRNA case) —
  expected; regenerate and sanity-check magnitudes (RNA-only scenarios should move little).
- (b) The `_compute_side` refactor must be formula-identical — verify by golden diff isolated to
  the density change (ideally land the helper refactor as a bit-identical commit first).
- (c) Run-fill edge cases (one-sided runs, all-NaN refs) — covered by unit tests above.
- (d) `count_evidence` over-confidence on thin exons (Issue 2) — bounded; Phase 3 addresses.

## 9. Implementation checklist (next turn)

1. `density_model.py`: add `_clean_gdna_frac`; rewrite `node_gdna_density` (local estimator + run-fill;
   new signature); delete the sweep (`from_left/from_right`, `weight`, `_TRAFFIC_PSEUDOCOUNT`,
   `seed_density`).
2. `calibrate.py`: drop the `gdna_frac` block; pass `rna_sense_frac` to `node_gdna_density`.
3. `joint_deconv._compute_side`: use `_clean_gdna_frac` (bit-identical).
4. Tests: new `tests/calibration/test_density_model.py` (unit + helper + run-fill); regenerate golden;
   full suite.
5. Debug check: `impute_prototype.py` / a GENE0037 dive confirming exon density ~20 and the leak drop.
