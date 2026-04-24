# Capture-class density channel: implementation plan

**Date:** 2026-04-23
**Status:** Design approved. Ready to implement.
**Supersedes:** [density_nb_model_plan.md](density_nb_model_plan.md). The
Negative-Binomial / Cox-Reid / φ_G approach is **retired**; the
underlying identifiability problem (which regions are captured?) cannot
be recovered from count moments alone and must be supplied as side
information.
**Scope:** Fix the calibration density channel on hybrid-capture data
by partitioning regions into `{on, off}` capture classes via an
optional targets BED. Both λ_G **and** (μ_R, σ_R²) are stratified.
Strand and FL channels are untouched. Fusion stays unit-weight sum of
LLRs. When no BED is supplied, behavior is bit-identical to today.

---

## 1. Why BED-based, not auto-detection

Within a single library the density channel cannot distinguish:

- **Off-expression enriched background.** An intron or UTR adjacent to
  a captured exon receiving 50–150× gDNA oversampling from probe
  bleed.
- **Lowly-expressed authentic mRNA.** An off-target transcript at low
  abundance producing the same per-bp rate.

Both yield indistinguishable `k/E` moments. Any "intelligent"
auto-detector that relies on neighboring-region rates, GC content, or
probe-density proxies is just rediscovering the BED file from noisy
surrogates. The vendor-supplied BED **is** the ground truth those
proxies approximate; we use it directly.

The previously-considered NB-with-learned-φ_G approach absorbs the
on/off rate heterogeneity by widening variance globally — losing
discriminating power everywhere, not just in off-target zones. It is
retired in favor of the partitioned Poisson described below.

---

## 2. Model

### 2.1 Capture class

Each region *i* is assigned a capture class `c_i ∈ {on, off}` at
calibration-prep time (fresh per `rigel quant` call, not baked into
the index). Rule, with padding `p` (default 150 bp):

```
c_i = on  ⇔  (mappable_bp overlap with any BED interval padded by ±p)
             ≥ 0.5 · mappable_bp_i
```

When `targets_bed=None`, every region is `on`; the partitioned model
collapses exactly to the current single-class model (see §3.3).

### 2.2 Stratified density kernel

For each class `c ∈ {on, off}` the two-component Poisson-LogNormal
mixture runs *independently*:

$$P(k_i \mid G) = \text{Poisson}\!\left(k_i;\; \lambda_G^{(c_i)} E_i\right)$$

$$P(k_i \mid R) = \int_0^\infty \text{Poisson}\!\left(k_i;\; (\lambda_G^{(c_i)} + \mu)\,E_i\right)\,\text{LogNormal}\!\left(\mu;\; \mu_R^{(c_i)}, (\sigma_R^{(c_i)})^2\right)\,d\mu$$

evaluated by the existing 21-node Gauss-HermiteE quadrature. Six
density parameters total:
`(λ_G_on, μ_R_on, σ_R_on, λ_G_off, μ_R_off, σ_R_off)`.

**Why both classes stratify.** Hybrid-capture probes enrich
*everything* in the targeted region — gDNA **and** mature RNA. Expressed
on-target transcripts produce read counts 100–1000× higher than their
off-target counterparts. If we enforced a single `(μ_R, σ_R²)`, the
LogNormal would either have to explode σ_R² to accommodate the
bimodal mass, or μ_R would drift high enough that off-target mRNA
looks like gDNA. Stratifying R parallels G and eliminates the trap.

### 2.3 Strand, FL, and fusion — unchanged

- Strand channel: same BetaBinomial on unique `(k, n)` pairs.
- FL channel: same γ-weighted histogram comparison.
- Fusion: `logit(γ_i) = log(π_soft / (1-π_soft)) + LLR_count_i +
  LLR_strand_i + LLR_fl_i`, unit-weight. Well-calibrated density
  makes this fusion correct by construction; the `(2·SS−1)²`
  weighting question and the temperature-scaling / Fisher-info
  adjustments in [channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md)
  become moot.

---

## 3. M-step: partitioned closed-form fits

Let `m_c = eligible ∧ (c_i = c)` for `c ∈ {on, off}`. Closed-form
updates, evaluated twice:

$$\hat\lambda_G^{(c)} = \frac{\sum_{i \in m_c}\gamma_i k^u_i}{\sum_{i \in m_c}\gamma_i E_i}$$

`(μ_R^{(c)}, σ_R^{(c)})` fit on the hard-spliced R anchors within
class `c` (`n_s > 0 ∧ eligible ∧ c_i = c`), using the existing
log-rate excess formula:

```
μ̂_i = max(k_i / E_i − λ_G^{(c)}, ε(E_i))
log_μ̂_i over anchors in class c
μ_R^{(c)}  = mean(log_μ̂)
σ_R^{(c)}² = max(var(log_μ̂), 1e-4)
```

π, π_soft remain **global** scalars — the soft fraction is a property
of the library's contamination, not the panel.

### 3.1 Fallbacks for sparse classes

If an anchor class has `< 2` R-anchors (possible on highly
unbalanced panels when only a handful of targets contain annotated
transcripts), fall back that class's `(μ_R, σ_R)` to the **other
class's** fit rather than the broad uninformative prior. Rationale:
the expression prior has some biological continuity across capture
classes; it's strictly a worse mis-specification to use a diffuse
prior than to borrow the neighboring partition's fit.

If *no* anchors exist in either class (synthetic/tiny cases), both
classes share the existing `mu_R = log(max(10·λ_G_on, 1e-6))`,
`σ_R = 2.0` fallback.

### 3.2 Fallback for sparse G-weight

If `Σ γ_i · E_i < ε` for a class (can happen for `off` on very clean
WGS, or `on` on unstranded pure-RNA), `λ_G^{(c)} = 0` with the
existing `_EPS` guard. Poisson at `λ = 0` degenerates correctly:
`P(k | G) = 1[k=0]`. Regions with any read in that class get γ → 0
(pure R), which is the right answer.

### 3.3 Reduction to current behavior (bit-identical)

When `targets_bed=None`:

- `c_i ≡ on` for all regions.
- `m_off = ∅`; `λ_G_off` initialized to 0 and never updated.
- `m_on = eligible` exactly.
- Partitioned MLE on `m_on` = current global MLE.
- Anchor fit on `m_on` = current anchor fit.

No numerical path differs from today's code. Existing golden outputs
and all 1065 tests should pass unchanged.

---

## 4. BED annotation

### 4.1 Input validation

- File must exist, be BED3+ (cols: ref, start, end; extra cols
  ignored).
- Reject if: `< 100` intervals (likely user error); total interval
  coverage `> 50%` of genome (not a capture panel — probably a full
  chromosome list supplied by mistake).
- Refs in BED not matching `region_df["ref"]` produce a warning but
  not a hard error (common: `chr1` vs `1` prefix mismatches — but
  we do **not** auto-fix; user must supply a BED matching their
  reference naming).

### 4.2 Algorithm

Use the existing vendored cgranges (`_cgranges_impl`) for the
interval join. Pseudocode:

```python
def annotate_capture_class(
    region_df: pd.DataFrame,
    bed_path: Path,
    pad: int = 150,
    min_overlap_frac: float = 0.5,
) -> np.ndarray:  # bool, len == len(region_df), True = on-target
    intervals = _read_bed3(bed_path, pad=pad)
    cg = cgranges.build(intervals)  # one per chromosome
    on = np.zeros(len(region_df), dtype=bool)
    mappable = region_df["mappable_effective_length"].values  # or "length"
    for chrom, group in region_df.groupby("ref", sort=False):
        overlaps = cg.overlap_sum(chrom, group["start"], group["end"])
        on[group.index] = overlaps >= min_overlap_frac * mappable[group.index]
    return on
```

Single sorted-interval sweep, O((N_regions + N_bed) log) per
chromosome. On the human index (~684k regions, typical panel ~300k
intervals), annotation runs in < 5 s.

### 4.3 Where to wire it in

At `rigel quant` startup, after `BamScanConfig` is loaded and before
`calibrate_gdna` is called (in `src/rigel/pipeline.py`, near line
903). Emit the capture mask as a `bool` array alongside `stats_dict`,
passed into `run_em` as a new `capture_class` kwarg. This keeps the
index panel-agnostic — same index supports any number of capture
panels.

### 4.4 No-BED path

When `targets_bed is None`:

- Skip annotation entirely.
- `capture_class = np.ones(N, dtype=bool)`.
- All subsequent code paths are single-class (see §3.3).

---

## 5. Code changes

### 5.1 Config

`src/rigel/config.py`, `CalibrationConfig`:

```python
#: Path to targets BED (BED3+) for hybrid-capture panels. When set,
#: regions are partitioned into {on, off} capture classes and the
#: density channel fits (λ_G, μ_R, σ_R²) independently per class.
#: When None (default), all regions are on-target and behavior is
#: identical to whole-genome / total-RNA inputs.
targets_bed: Path | None = None

#: Padding (bp) applied to every BED interval before overlap
#: scoring. Compensates for fragment bleed from capture probes:
#: a 300 bp fragment tiled on a 100 bp exon will overhang the
#: exon edges by ~150 bp.
targets_pad: int = 150

#: Minimum fraction of a region's mappable_bp that must lie inside
#: the padded target intervals to be classified on-target.
targets_min_overlap_frac: float = 0.5
```

### 5.2 CLI

`src/rigel/cli.py`, **standard** (not advanced) quant options:

- `--targets PATH` (dest: `targets_bed`)
- `--targets-pad INT` (dest: `targets_pad`, default 150)
- `--targets-min-overlap-frac FLOAT` (dest: `targets_min_overlap_frac`,
  default 0.5)

Plumbed via the existing `_ParamSpec` registry:

```python
_ParamSpec("targets_bed", "calibration.targets_bed", "path_or_none"),
_ParamSpec("targets_pad", "calibration.targets_pad"),
_ParamSpec("targets_min_overlap_frac", "calibration.targets_min_overlap_frac"),
```

### 5.3 Calibration code

`src/rigel/calibration/_annotate.py` (**new**):

- `annotate_capture_class(region_df, bed_path, pad, min_overlap_frac)
  -> np.ndarray[bool]`
- `_read_bed3(bed_path, pad) -> dict[str, list[(start, end)]]`

`src/rigel/calibration/_em.py`:

- `_count_llr_poisson_ln` gains a `capture_class: np.ndarray[bool]`
  kwarg; evaluates the quadrature on-class and off-class with
  different `(λ_G, μ_R, σ_R)` tuples and concatenates into a single
  per-region LLR vector.
- `_m_step` signature gains `capture_class`; returns tuples
  `(λ_G_on, λ_G_off), (μ_R_on, μ_R_off), (σ_R_on, σ_R_off), π, π_soft`.
- `run_em` signature gains `capture_class: np.ndarray[bool] | None`.
  When `None`, defaults to all-True (identical to §3.3). Diagnostics
  dict records all six parameters per iter.
- `EMFit` / `CalibrationResult` gains fields:
  `lam_G_on, lam_G_off, mu_R_on, mu_R_off, sigma_R_on, sigma_R_off,
  n_regions_on, n_regions_off, on_target_enrichment_ratio`
  (last one = `lam_G_on / max(lam_G_off, ε)`).

`src/rigel/calibration/_calibrate.py`:

- `calibrate_gdna` gains a `capture_class` kwarg, plumbs it through
  to `run_em`. No BED I/O here — annotation is upstream.

`src/rigel/pipeline.py`:

- Before `calibrate_gdna`, call `annotate_capture_class` if
  `cal_cfg.targets_bed is not None`. Log: target-region count,
  on-target region count, total mappable_bp on/off.

### 5.4 Summary / QC

`summary.json` gains, under `calibration`:

```
"density_model": {
  "lam_G_on":  7.4e-3,
  "lam_G_off": 8.1e-5,
  "on_target_enrichment_ratio": 91.4,
  "mu_R_on":   -3.1,
  "sigma_R_on": 1.4,
  "mu_R_off":  -4.5,
  "sigma_R_off": 1.8,
  "n_regions_on":  42318,
  "n_regions_off": 641598,
  "targets_bed": "/path/to/panel.bed",
  "targets_pad": 150
}
```

When `targets_bed is None`, the on/off split is omitted and only
`lam_G, mu_R, sigma_R` appear (as today).

---

## 6. Validation

### 6.1 Unit tests (new)

- `test_annotate_capture_class_basic`: synthetic 3-region toy BED,
  verify overlap fraction thresholding.
- `test_annotate_capture_class_padding`: a region that's 0-overlap
  without padding but ≥50% with 150 bp pad flips to on-target.
- `test_annotate_capture_class_mismatched_refs`: BED with `chr1`
  vs region_df with `1` emits warning, classifies those regions as
  off.
- `test_m_step_partitioned`: synthetic stats with bimodal λ_G
  (1e-4 off, 1e-2 on); verify partitioned MLE recovers both rates
  within 5%.
- `test_run_em_no_bed_bit_identical`: with `capture_class=None`,
  assert `EMFit` matches current code to 1e-12.

### 6.2 Regression (existing)

- `tests/test_calibration*.py`: pass unchanged with `targets_bed=None`.
- `tests/test_golden_output.py`: pass unchanged.
- Full suite (1065 tests) green before any default change.

### 6.3 VCaP titration validation

Eight libraries at
`/scratch/.../runs/human/mctp_vcap_rna20m_dna{00,01,02,05,10,20,40,80}m/`,
with MCTP capture panel BED. For each, run calibration with and
without `--targets`.

Primary checks:

1. **Enrichment ratio monotone with gDNA spike.** At dna00m, enrichment
   ratio is undefined (no gDNA); at dna80m, should plateau near the
   panel's true on-target capture efficiency (typically 40–150×).
2. **`quant.gdna_total / strand_inferred_gdna` drops to 1.0–1.1×** across
   the titration (currently 1.1–1.5× without BED).
3. **Per-region γ correlation with strand-only baseline** improves
   from ~0.5 to > 0.9 at dna10m+.
4. **dna00m behavior unchanged** — the strand channel already handles
   pure-RNA cleanly; BED-mode should not harm it.

### 6.4 No-regression on synthetic benchmarks

`scripts/benchmark/configs/locus_simple_*.yaml` — these have no
capture structure. Runs with `targets_bed=None` (the default) should
be bit-identical to current goldens. This is the primary regression
guard for "did we break WGS / total-RNA?"

---

## 7. Rollout

| Phase | Action | Gate |
|---|---|---|
| 1 | Implement annotator + partitioned EM, all behind `targets_bed=None` default. | Full test suite green; bit-identical goldens with no BED. |
| 2 | VCaP titration validation (§6.3). Publish enrichment-ratio, γ-correlation, and gdna-total-ratio plots. | §6.3 predictions 1–4 hold. |
| 3 | Document in [MANUAL.md](../MANUAL.md) and [METHODS.md](../METHODS.md): capture-aware density is available via `--targets`; remove `(2·SS−1)²` section from METHODS since fusion is now well-calibrated. | §6.3 gate passed. |
| 4 | (Optional, future) Extend partitioning to support multi-class panels (e.g., WGS-baseline + targeted). Not in first release. | — |

The `--targets` flag ships as a standard, user-facing option —
**not** opt-in advanced. Users with capture panels simply supply
their panel BED; users without one don't pass the flag and see
identical behavior to today.

---

## 8. What this deletes / simplifies

- **Retires** [density_nb_model_plan.md](density_nb_model_plan.md)
  — no NB kernel, no φ_G, no Cox-Reid, no warm-started golden
  section, no iter-0 skip, no hard-G subset bookkeeping.
- **Kills** the channel-fusion weighting question: unit-weight LLR
  summation is correct once density is correct. No
  `(2·SS−1)²`-style weights needed. No temperature scaling. No
  Fisher-info-based reweighting.
- **Removes** the "strand-anchored local λ_G re-initialization"
  line of thinking (Class C in
  [channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md)).
  Density stops being confidently wrong, so none of the re-ordering
  gymnastics around it are needed.
- **Leaves untouched**: strand channel (BetaBinomial with learned
  κ_G_min floor), FL channel (γ-weighted histogram), π_soft update,
  all E/M ordering.

---

## 9. Open questions (none currently blocking)

- **Multi-panel / tiered panels.** Some vendors ship "dense" and
  "sparse" target regions with different expected enrichment. First
  release treats all BED intervals identically; extending to
  multi-class is a straightforward generalization of §3 (N classes
  instead of 2) but deferred.
- **Probe-density weights.** Some vendors supply per-interval
  hybridization strength. We ignore this in v1; the padded
  intersection is good enough given the Poisson noise floor.
- **Unstranded hybrid capture.** Strand channel provides no signal;
  density + FL is the only path. Capture-class density restores
  density's validity, so unstranded + `--targets` should work as
  well as stranded + `--targets`. Worth confirming empirically in
  Phase 2.

---

## 10. References

- [channel_fusion_hybrid_capture.md](channel_fusion_hybrid_capture.md)
  — parent research doc (Class A.1 "targets BED mode" was flagged
  as the purest solution; this plan realizes it and supersedes the
  NB variance-absorption alternative).
- [density_nb_model_plan.md](density_nb_model_plan.md) — retired.
- [SESSION_HANDOFF_2026-04-23.md](SESSION_HANDOFF_2026-04-23.md)
  — prior session state that decided the NB direction; that
  decision is reversed by this plan based on the
  "captured RNA is enriched too" biological observation, which
  forces stratified `(μ_R, σ_R²)` alongside `λ_G`.
- [src/rigel/calibration/_em.py](../../src/rigel/calibration/_em.py)
  — current single-class density kernel (`_count_llr_poisson_ln`,
  `_m_step`) to be partitioned.
