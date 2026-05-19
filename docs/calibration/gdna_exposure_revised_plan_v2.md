# gDNA Regional Exposure Revised Plan v2

**Status**: implementation-ready after review fixes R1-R6, not implemented.
**Date**: 2026-05-18
**Supersedes**: `gdna_exposure_revised_plan_v1.md`.
**Implements lessons from**: `gdna_exposure_plan_v3.md` (the failed production attempt).

This revision keeps v1's core insight — regional exposure must be applied
symmetrically to every hypothesis whose molecule occupies the exposed
genomic location — and adds (a) an explicit generative model so the
implementation cannot drift from the math, (b) a mitigation for a new
mRNA-vs-gDNA risk that v1 underweighted, (c) quantitative validation
gates derived from the observed regression, and (d) specific rules for
shared nRNA spans, cross-reference candidates, and the v6 gDNA prior
count calibration.

Review fixes folded into this version:

- R1: exposure numerator weighting is restricted to units with trustworthy
  unit-level genomic coordinates; multimappers stay unweighted unless a later
  native change emits per-candidate/per-hit genomic positions.
- R2: nRNA exposure follows the current scorer semantics: v2 weights synthetic
  nRNA rows, while annotated nascent-equivalent single-exon rows remain on the
  mRNA/annotated-transcript path and are monitored by the mRNA guard rail.
- R3: nRNA span ids are derived at index load time from existing transcript
  columns, so the first implementation does not require an index format bump.
- R4: the pipeline order is aligned to the current architecture, where
  `CalibrationResult.regional_exposure` is built by calibration before
  `quant_from_buffer` begins.
- R5: Mode 3b is explicitly a guard-rail branch with additional data
  requirements, not hidden scope inside the symmetric-unspliced patch.
- R6: bit-exactness is required only for uniform fast paths; hand-constructed
  constant non-identity exposure fields use numerical equivalence tolerances.

---

## 1. Executive Summary

The v3 production attempt failed on hybrid-capture data not primarily
because of the cross-class normalization bug, but because the exposure
field `A(x)` was applied only to gDNA. nRNA — a competing unspliced
genomic hypothesis at the same physical location — was left unweighted,
producing a `log(A(x) / Ā_locus)` likelihood-ratio penalty against
gDNA that flipped true gDNA fragments to nRNA.

v2's principle:

> Regional exposure `A(x)` is an assay observation-process term. Apply
> it to every component whose molecule could plausibly emit a fragment
> from genomic position `x`, with that component's own
> fragment-length distribution and span geometry.

Operationally, in v2 we apply `A(x)` symmetrically to **gDNA and
synthetic nRNA**. This matches the current native scorer's pool
separation, where only `is_synthetic` rows are treated as the nRNA
competition pool. Annotated nascent-equivalent single-exon rows
(`is_nrna=True, is_synthetic=False`) remain on the annotated transcript
path in v2 and are covered by the mRNA guard rail (§6). mRNA is treated
specially because its genomic footprint is restricted to exons and the
dominant hybrid-capture exposure regime makes uniform mRNA effective
length an acceptable approximation — *but* we add a Mode-3b test and a
guard rail to catch mRNA-vs-gDNA regressions if that approximation
breaks.

We do **not** ship v3's class-collapsed `rho_ref`. We use class-specific
references as the conservative bridge field, and treat the resulting
`A(x)` as "approximate within each class, inconsistent across classes"
until a unified smoothed exposure field replaces it.

---

## 2. Generative Model

This section is normative. The implementation must agree with these
equations.

### 2.1 Per-component fragment likelihood

For each EM unit `u` at genomic position `x_u` with length `ell_u`:

#### gDNA (component `g`, per-locus `M`)

$$
P(u \mid \text{gDNA}, M) \;=\;
\frac{A(x_u)\;\cdot\;h_G(\ell_u)\;\cdot\;\mathbb{1}[u \text{ valid in } M]}
     {\displaystyle\sum_\ell h_G(\ell)\;\int_{\text{valid}(M, \ell)} A(s)\,ds}
$$

Let $L_g^A(M) \equiv \sum_\ell h_G(\ell)\,\int A(s)\,ds$ — the
exposure-weighted gDNA effective length already implemented as
`weighted_gdna_eff_len_for_loci`.

#### Synthetic nRNA (component `n`, per nRNA-bearing transcript `t`)

$$
P(u \mid \text{nRNA}, t) \;=\;
\frac{A(x_u)\;\cdot\;h_{RNA}(\ell_u)\;\cdot\;\mathbb{1}[u \text{ valid on } t]}
     {\displaystyle\sum_\ell h_{RNA}(\ell)\;\int_{\text{valid}(t, \ell)} A(s)\,ds}
$$

Let $L_n^A(t) \equiv \sum_\ell h_{RNA}(\ell)\,\int A(s)\,ds$. Three
important details:

- The FL distribution is the **RNA** FL ($h_{RNA}$), not the gDNA FL.
- The integration domain is the **unspliced genomic span** of the nRNA
  candidate (i.e., the shared nRNA span keyed by
  `(ref, strand, start, end)`), not the transcript's exon blocks.
- Multiple transcripts may share an nRNA span. $L_n^A$ is computed
  once per span and broadcast to all owning transcripts (mirrors the
  existing unweighted nRNA architecture).

#### Annotated mRNA (component `m`, per transcript `t`)

$$
P(u \mid \text{mRNA}, t) \;=\;
\frac{A(x_u)\;\cdot\;h_{RNA}(\ell_u)\;\cdot\;\mathbb{1}[u \text{ valid on exons of } t]}
     {\displaystyle\sum_\ell h_{RNA}(\ell)\;\int_{\text{exonic valid}(t, \ell)} A(s)\,ds}
$$

Two practical observations:

- mRNA candidates can only be evaluated at positions on transcript
  exons. Outside of capture targets, `A(x)` on exons should
  systematically equal the high-exposure tail of the field.
- Therefore $L_m^A(t) \approx \bar{A}_{\text{exonic}}\,L_m(t)$, where
  $\bar{A}_{\text{exonic}}$ is approximately constant across captured
  transcripts. The per-fragment density on exons becomes
  $A(x_u) / (\bar{A}_{\text{exonic}}\,L_m(t)) \approx 1/L_m(t)$ at any
  exonic $x_u$ in a captured target. Weighting mRNA is approximately
  the identity in the dominant case.

This is the justification for **not** weighting mRNA by default in v2,
but we will measure it (§9.2 Mode 3b) rather than assume it.

### 2.2 EM scoring equation (current native EM, unchanged ABI)

The native EM scores each candidate as:

$$
\text{score}_k(u) = \log\!\text{lik}_k(u) + \log\theta_k - \log L_k
$$

To realize the equations in §2.1 we therefore:

- **Add** `log A(x_u)` to `log lik_k(u)` for $k \in \{\text{gDNA}, \text{nRNA}\}$
  candidates of unspliced units.
- **Replace** $L_g$ with $L_g^A$ in the per-locus gDNA effective length
  array.
- **Replace** $L_n(t)$ with $L_n^A(t)$ in the per-transcript effective
  length array for synthetic nRNA rows only.

The native solver does not need to know about exposure. All math lives
in the Python layer that prepares the numerator log-likelihoods and the
effective length arrays.

### 2.3 Cancellation invariant

For a single component $k$ in a locus with constant $A(x) = A_0$:

$$
\Delta\text{score}_k(u) = \log A_0 - \log(A_0 \cdot L_k / L_k) = 0
$$

The current implementation already satisfies this for gDNA in the
uniform fast path. v2 extends this invariant to nRNA.

Implementation invariant:

- If `exposure.mode == "uniform"`, the path must be bit-exact versus the
  unweighted baseline.
- If a test hand-constructs a constant non-identity exposure field
  (`A(x) = A_0 != 1`), scores and posteriors must be numerically
  equivalent to the unweighted baseline within floating-point tolerance;
  bit-exactness is not required because the weighted-length integration
  and `log(A_0)` additions may round differently.

---

## 3. Why Symmetry Matters (concise re-derivation)

Comparing two components $k_1, k_2$ at fragment $u$, the
exposure-induced shift in their log-odds is:

$$
\Delta\!\log\text{LR}_{k_1/k_2}(u) =
\bigl[\log A(x_u) - \log\bar{A}_{k_1}(u)\bigr] -
\bigl[\log A(x_u) - \log\bar{A}_{k_2}(u)\bigr]
= \log\bar{A}_{k_2}(u) - \log\bar{A}_{k_1}(u)
$$

where $\bar{A}_k(u) \equiv L_k^A / L_k$ is the FL-weighted average
exposure over component $k$'s valid window at $u$.

If both $k_1, k_2$ receive identical exposure treatment **and** their
valid windows have the same shape (gDNA over the locus, nRNA over the
same unspliced span), then $\bar{A}_{k_1} \approx \bar{A}_{k_2}$ and
the exposure terms cancel — *exposure cannot move mass between them*.

If one is weighted and the other is uniform (v3 today), the log-odds
shift becomes $\log A(x_u) - \log\bar{A}_{k_1}$, which is exactly the
$-4$ to $-9$ nat penalty against gDNA at intronic positions that we
observed.

This is the formal reason gDNA-and-nRNA-symmetric weighting eliminates
the regression while preserving the across-locus gDNA boost in
low-exposure loci (because $\bar{A}_n$ and $\bar{A}_g$ shrink together
when the locus is dominated by low-A intron).

---

## 4. Exposure Field `A(x)` for v2

### 4.1 Bridge construction (ship with v2)

- Use the existing `_regional_exposure.py` skeleton: per-region
  conservative-count tables `Y_r`, opportunity `E_r`, EB-shrunk
  `rho_hat_r`.
- **Fix the cross-class normalization bug**: each class normalizes to
  its own reference quantile.

  ```python
  for cname in CLASSES:
      log_rho_ref_k = log(_weighted_quantile(rho_hat_k, E_k, REFERENCE_QUANTILE))
      raw_log_ratio_r = min(log(rho_hat_r) - log_rho_ref_k, 0.0)
      log_A_r = signal_k * raw_log_ratio_r
      log_A_r = max(log_A_r, LOG_A_FLOOR)
  ```

  Summary schema keeps the existing numeric `rho_ref` field as the
  informational maximum for backward compatibility, and adds
  `rho_ref_per_class` as a dictionary keyed by class. The numeric maximum is
  never used to compute weights.

- **Keep auto-uniform signal logic and the deterministic null-spread
  approximation** unchanged.

### 4.2 Acknowledged limitation of the bridge field

After class-specific normalization, `A(x)` is internally consistent
*within* each class but only approximately comparable across classes.
A per-locus integral $\int_M A(s)\,ds$ that crosses class boundaries
mixes three differently-calibrated relative scales.

For v2 this is acceptable because:

1. The symmetric application means class-mismatch errors are absorbed
   identically by both gDNA and nRNA denominators, so they largely
   cancel in the gDNA-vs-nRNA log-odds.
2. The remaining mRNA-vs-{gDNA,nRNA} bias is bounded (§6 and the
   exonic-flat approximation in §2.1).

A unified smoothed field is deferred to v3 of this plan (separate doc).

### 4.3 Per-class normalization invariants enforced in code

- `rho_ref_per_class[k] > 0` for every class with positive total
  exposure. If a class has zero exposure, its regions get `A_r = 1`
  (identity).
- `signal_k = 0 ⇒ A_r = 1` for the class (auto-uniform fallback).
- All-classes `signal_k == 0 ⇒` exposure object reports
  `mode == "uniform"` and the entire downstream path goes through the
  bit-exact fast path.

### 4.4 Sentinel test (lock in the v3 RCA)

Add a regression test that constructs three synthetic classes with
`rho_ref_EXON / rho_ref_INTRON ≥ 100`. Assert:

- An INTRON region whose `rho_hat` equals INTRON's `rho_ref_class`
  receives `A_r` within `1e-9` of `1.0`.
- An EXON region whose `rho_hat` equals EXON's `rho_ref_class` also
  receives `A_r` within `1e-9` of `1.0`.

This locks the class-normalization fix in place permanently.

---

## 5. nRNA Implementation Specifics

### 5.1 Per-candidate nRNA weighting on safe units

For a unit `u` with candidates `j` in CSR positions `offsets[u]:offsets[u+1]`,
weight only candidates whose transcript is a synthetic nRNA row and
whose unit is unspliced:

```python
is_synth_t = index.t_df["is_synthetic"].to_numpy(dtype=bool)
is_nrna_cand = is_synth_t[em_data.t_indices]

# Unit-level safety mask shared with the gDNA numerator path.
# The current CSR stores one genomic_midpoint per EM unit, not one per
# candidate or per multimapping hit. Therefore v2 must only apply local
# numerator exposure to units for which that single midpoint is trustworthy
# for every weighted candidate.
safe_unit = (
    (~em_data.is_spliced.astype(bool))
    & (em_data.genomic_midpoint != INT64_MIN)
    & same_ref_unit(em_data, index)
    & non_multimapper_unit(em_data)
)

for unit_block in iter_unit_blocks(em_data.offsets):
    c0 = int(em_data.offsets[unit_block.start])
    c1 = int(em_data.offsets[unit_block.stop])
    cand_idx = np.arange(c0, c1, dtype=np.int64)
    cand_unit = local_candidate_unit_ids(unit_block, em_data.offsets)
    mask = (
        safe_unit[cand_unit]
        & is_nrna_cand[cand_idx]
        & np.isfinite(em_data.log_liks[cand_idx])
    )
    if not mask.any():
        continue
    selected = cand_idx[mask]
    unit_ids = cand_unit[mask]
    ref_ids = index.t_to_ref_arr[em_data.t_indices[selected]]
    positions = em_data.genomic_midpoint[unit_ids]
    log_w = exposure.log_weights_for_positions(ref_ids, positions)
    em_data.log_liks[selected] += log_w.astype(em_data.log_liks.dtype, copy=False)
```

Safety rules:

- **Per-candidate candidate selection**: only synthetic nRNA candidates
  inside a safe unit are mutated. Other candidates in the unit are left
  unchanged.
- **Shared safe-unit mask with gDNA**: v2 updates `_apply_unit_gdna_weights`
  to use the same safety rule. If gDNA is not safely weightable for a unit,
  nRNA is not locally weighted either. This preserves symmetry.
- **No multimapper local exposure in v2**: multimapper units can contain
  candidates from different alignments/positions but the current CSR keeps only
  one `genomic_midpoint` (`mm_gmid_first`). Applying that coordinate to every
  nRNA candidate would be wrong. Multimapper numerator exposure remains
  identity until a later native change emits per-candidate or per-hit genomic
  positions.
- **Memory bounded**: do not materialize a whole-run `cand_unit = np.repeat(...)`
  vector for large human runs. Process units in blocks and build block-local
  candidate-unit ids.
- **Spliced units are skipped**: nRNA is by definition unspliced;
  spliced units are not nRNA candidates regardless of `t_indices` flag.

Diagnostic counters should distinguish `cross_ref_skipped` and
`multimapper_skipped` so the VCaP report can quantify how much data remains in
the identity path.

### 5.2 Per-span weighted effective length

Add to `src/rigel/calibration/_exposure.py`:

```python
def weighted_nrna_eff_lens(
    nrna_spans: NrnaSpanTable,   # (ref, strand, start, end, span_id)
    rna_fl: FragmentLengthModel,
    exposure: RegionalGdnaExposure,
    *,
    min_value: float = 1.0,
) -> np.ndarray:
    """Per-unique-span weighted effective length under RNA FL."""
```

Algorithm mirrors `weighted_gdna_eff_len_for_loci` but:

- Operates on individual `(ref, start, end)` spans (no `Locus`
  intervals).
- Uses `rna_fl.pmf` instead of `gdna_fl.pmf`.
- Uses `min_value=1.0` for numerical safety.

Then in the pipeline:

```python
nrna_span_weighted_L = weighted_nrna_eff_lens(...)
# Broadcast span -> transcript:
t_eff_lens_weighted = t_eff_lens.copy()
synth_mask = index.t_df["is_synthetic"].to_numpy(dtype=bool)
t_eff_lens_weighted[synth_mask] = nrna_span_weighted_L[
    index.t_to_nrna_span_id[synth_mask]
]
```

`index.t_to_nrna_span_id` is the per-transcript span id for synthetic nRNA
rows. For v2, derive it at `TranscriptIndex.load()` time from existing
`transcripts.feather` columns (`ref`, `strand`, `start`, `end`, `is_synthetic`,
`nrna_t_index`) and expose:

```python
index.nrna_span_df          # span_id, ref_id, strand, start, end, nrna_t_index
index.t_to_nrna_span_id     # int32[n_transcripts], -1 for non-synthetic rows
```

This avoids an index format bump for the first implementation. If a later patch
persists an explicit span table, bump `INDEX_FORMAT_VERSION` and require index
rebuilds in that patch.

When `exposure.mode == "uniform"`, return the existing unweighted nRNA
effective lengths bit-exactly (fast path).

### 5.3 EM hand-off

Currently `AbundanceEstimator` consumes a single per-transcript
`t_eff_lens` for all transcript-like components. Replace it with
`t_eff_lens_weighted` (constructed once per pipeline run, after
calibration). No native ABI changes.

---

## 6. mRNA Treatment and the mRNA→gDNA Risk

### 6.1 Default: leave mRNA unweighted in `auto`

Per §2.1, mRNA's exonic footprint sits in the high-A tail of the field
for hybrid-capture targets, so weighting mRNA is approximately the
identity. Keeping it unweighted limits blast radius.

### 6.2 Quantified risk and detection

The risk is: at exonic positions in a low-`Ā_locus` locus, gDNA
and nRNA denominators shrink (because `L^A` integrates over both
high-A exons and low-A introns), boosting gDNA's per-fragment density.
mRNA's denominator does not shrink. The relative log-odds at an exonic
position becomes:

$$
\Delta\!\log\text{LR}_{g/m}(u_\text{exon}) =
\log A(x_u) - \log\bar{A}_g(M)
\;\approx\;
\log A_\text{exon} - \log\bar{A}_g(M)\,.
$$

If `A_exon ≈ 1` (capture target) and $\bar{A}_g(M) \ll 1$
(intron-dominated locus), this is positive — mRNA loses mass to gDNA
at exonic positions in those loci. Magnitude is bounded by
$-\log\bar{A}_g(M)$, which equals `-log(gdna_eff_len_weight_ratio)` per
locus from the diagnostic schema we already emit.

### 6.3 Guard rail: Mode 3b

Add an internal mode that also weights mRNA candidates and denominators
(same algorithm as nRNA, swapping the integration domain to per-exon
windows of the mature transcript). Run it in the validation sweep
(§9.2). If Mode 3 shows an mRNA→gDNA regression on VCaP truth that
Mode 3b corrects, the production default shifts to Mode 3b.

Mode 3b has extra data requirements that Mode 3 does not: numerator weighting
for mRNA needs an exonic footprint or midpoint for each transcript candidate,
not just the unit-level genomic midpoint. The current global CSR does not carry
that information for every candidate, especially for multimappers and spliced
fragments. Therefore Mode 3b is a guard-rail branch with an explicit decision
point:

1. First implement Mode 3 and the `mrna_gdna_logratio_shift` diagnostic.
2. Run Mode 3 on VCaP.
3. If G4 passes, ship Mode 3 and leave Mode 3b as a follow-up.
4. If G4 fails, do not ship Mode 3. Implement Mode 3b by adding the required
  candidate-level exonic exposure coordinates (or an equivalent candidate
  footprint integral) and rerun the gates.

We commit to **using Mode 3b as the escape hatch if G4 fails**, not to hiding a
native/scoring ABI expansion inside the first symmetric-unspliced patch. This
keeps the common case small while still blocking an mRNA regression before
merge.

### 6.4 Per-locus auditing column

Add to `loci.feather`:

- `mrna_gdna_logratio_shift` = `-log(gdna_eff_len_weight_ratio)` (this
  is exactly the predicted bound on the mRNA→gDNA shift in nats).

This makes it possible to identify loci that are at risk of an mRNA
regression by simple thresholding (e.g., shift > 4 nats).

---

## 7. v6 gDNA Prior Count Calibration

The v6 prior count `eta_g(M)` was derived assuming a uniform exposure
denominator. Once we ship symmetric weighting, two strategies exist:

**A. Keep `eta_g` numerically identical, accept that the prior is
slightly miscalibrated.** Justification: `eta_g` is a Dirichlet
pseudocount, and its role is purely regularization of small-data loci.
In hybrid-capture targets `gdna_eff_len_weight_ratio` is typically
in `[1e-3, 1]`; the implicit shift in posterior gDNA mass is bounded
and the prior is dominated by data for any locus with non-trivial
fragment count.

**B. Renormalize: `eta_g_v2 = eta_g * (L_g^A / L_g)`.** This keeps
`eta_g / L_g` invariant, i.e., the prior rate per unit of *weighted*
opportunity stays equal to the global gDNA density.

**v2 ships Strategy A** for parity simplicity, **and emits Strategy B
as a diagnostic** in `loci.feather` (`gdna_prior_count_renormalized`).
The validation sweep (§9.2) will include a Mode-3-with-Strategy-B run.
If Strategy B improves the gDNA→nRNA leak gate by ≥1 absolute percentage point
and does not worsen any of G1, G3, or G4, flip the default before merge.

---

## 8. Pipeline Order

Calibration already builds `CalibrationResult.regional_exposure` in
`rigel.calibration._orchestrator.calibrate(...)`, before
`quant_from_buffer(...)` starts. v2 keeps that ownership: exposure learning
stays in calibration; quantification consumes the learned table.

Inside `quant_from_buffer`:

1. Read `calibration.regional_exposure` (already class-normalized by §4.1).
2. Build per-transcript weighted effective lengths before constructing
  `AbundanceEstimator`:
  - annotated mature rows keep existing RNA effective lengths;
  - synthetic nRNA rows receive `L_nRNA^A` from §5.2;
  - uniform-mode fast path is bit-exact to the current unweighted array.
3. Construct `TranscriptGeometry` / `AbundanceEstimator` with the weighted
  transcript effective-length array.
4. Score fragments into global `ScoredFragments` (with `genomic_midpoint`
  already populated from v3 step 8).
5. Build `multi_loci`.
6. `assemble_priors(... regional_exposure=...)`. Computes `L_g^A` per-locus
  and diagnostic regional `eta` per-locus.
7. Apply local numerator exposure on the global CSR arrays using the shared
  safe-unit mask from §5.1:
  - mutate gDNA `gdna_log_liks` for safe unspliced units;
  - mutate synthetic nRNA candidate `log_liks` for the same safe units.
8. Partition and free.
9. Run EM with weighted transcript effective lengths and `L_g^A`.

The two weight-application calls must run *after* `assemble_priors`
(single-owner discipline for the global CSR arrays) and *before*
`partition_and_free`.

---

## 9. Validation Plan

### 9.1 Quantitative gates (production criteria for `auto`)

All gates evaluated on the same VCaP `rna20m_gdna20m` flowcell-truth
input with a single fixed seed.

| Gate | Threshold |
|---|---|
| G1: gDNA→RNA (mRNA+nRNA) leak | ≤ same-version uniform baseline + 1 percentage point |
| G2: gDNA→nRNA specifically | ≤ same-version uniform baseline + 0.5 percentage points |
| G3: RNA→gDNA leak | within ±0.5 percentage points of uniform baseline |
| G4: mRNA→gDNA leak | ≤ uniform baseline + 0.5 percentage points |
| G5: synthetic 24-condition mRNA relative error | ≤ uniform baseline + 0.5% |
| G6: synthetic 24-condition nRNA relative error | ≤ uniform baseline + 0.5% |
| G7: `--regional-exposure off` regression suite | bit-exact vs current `off` |

If G1–G4 hold for the symmetric-unspliced mode (Mode 3), ship Mode 3.
If G4 fails for Mode 3, do not ship; implement and test Mode 3b. If Mode 3b
then satisfies G1–G4, ship Mode 3b. If G2 fails for both, iterate on the
field-construction layer rather than tuning constants.

### 9.2 VCaP mode sweep

Run, in order:

| Mode | Description |
|---|---|
| 0 | Uniform (current `off`). Baseline. |
| 1 | v3 gDNA-only (current `auto`). Failed comparator. |
| 2 | gDNA denominator-only (no per-unit numerator on gDNA). |
| 3 | Symmetric unspliced: gDNA + nRNA weighted. Production candidate. |
| 3a | Mode 3 + Strategy B prior renormalization. |
| 3b | Mode 3 + mRNA weighting on exon footprints. Guard-rail branch; required only if Mode 3 fails G4, because it needs candidate-level mRNA footprint data not present in the current CSR. |
| 4 | Mean-centered numerator only (diagnostic, no denominator change). |

For each implemented mode, emit:

- Full truth-by-predicted confusion matrix (3×3).
- Pairwise transition matrix vs Mode 0.
- For false gDNA→RNA fragments only: histogram of
  `log A(x_u)`, `log Ā_locus`, `log(A(x_u)/Ā_locus)`.
- Per-locus `gdna_eff_len_weight_ratio`,
  `mrna_gdna_logratio_shift`.

### 9.3 CI-runnable synthetic locus

Add `tests/test_regional_exposure_locus.py` with one deterministic
synthetic locus:

- 5 kb high-exposure exon (`A_r = 1.0`).
- 95 kb low-exposure intron (`A_r = 0.01`).
- 100 true gDNA fragments uniformly along the genomic span.
- 100 true nRNA fragments uniformly along the same span.
- 50 true mRNA fragments restricted to the exon.

Assertions:

- Mode 3 (symmetric) recovers gDNA / nRNA counts within ±10% of
  truth and does not over-call nRNA by more than 15%.
- Mode 1 (gDNA-only) systematically misroutes intronic gDNA to nRNA
  (reproduces the regression as a *positive* assertion that v2
  preserves the diagnostic power).
- Mode 0 and Mode 3 with `exposure.mode == "uniform"` produce bit-identical
  assignments; a hand-constructed constant non-identity `A` field produces
  equivalent posteriors within floating-point tolerance (cancellation invariant
  from §2.3).

Runtime budget: < 2 s.

### 9.4 Synthetic 24-condition regression

Rerun the full benchmark grid with Mode 3. Apply gates G5–G7.

### 9.5 Oracle confirmation (one-shot, not gated)

For the next available oracle-tagged simulation, confirm that the
implicated gDNA→nRNA transitions in v3 are truly gDNA-origin. Flowcell
truth on VCaP is strong; oracle on synthetic is sanity.

---

## 10. CLI and Config Surface

User-facing flag is unchanged:

```text
--regional-exposure {off,auto}
```

Where `auto` resolves at runtime to the validated production mode
(Mode 3, or Mode 3a/3b if the gates require).

For development we add an internal-only enum (debug logging only, not
documented in `--help`):

```text
RIGEL_REGIONAL_MODE = off | gdna_only | gdna_denominator_only | unspliced_symmetric | unspliced_symmetric_renorm | full_symmetric | mean_centered
```

`full_symmetric` is unavailable until the Mode 3b data requirement in §6.3 is
implemented; selecting it before that should fail fast with a clear error.

Resolved via an env-var or a hidden CLI argument
`--regional-mode-internal` for use by the mode-sweep harness. Production
config files should never reference it.

`CalibrationConfig.regional_exposure_enabled` stays as the
top-level user knob (`auto` ↔ `True`). A second field
`CalibrationConfig.regional_mode_internal: str | None = None`
selects the internal mode when present.

---

## 11. Implementation Checklist

1. **Hygiene fix (small).** Replace global-max `rho_ref` with
   class-specific references in `_regional_exposure.py::build()`. Add
   the §4.4 sentinel test. *Land first.*
2. **Per-span nRNA bookkeeping.** Derive `index.t_to_nrna_span_id` and
  `index.nrna_span_df` at `TranscriptIndex.load()` from existing transcript
  columns. Do not bump the index format unless these are later persisted.
3. **Weighted nRNA effective lengths.** Implement
   `weighted_nrna_eff_lens` in `_exposure.py` with bit-exact uniform
   fast path. Tests in `tests/test_weighted_eff_len.py`.
4. **Shared safe-unit exposure mask.** Refactor the gDNA numerator applier and
   the new nRNA numerator applier to use one safe-unit mask: unspliced, valid
   midpoint, same transcript reference, and non-multimapper. Add counters for
   cross-ref and multimapper skips.
5. **Per-candidate nRNA numerator weights.** Implement
   `_apply_candidate_nrna_weights` in `pipeline.py`, mutating only synthetic
   nRNA candidates inside safe units. Tests in `tests/test_pipeline_routing.py`
   family.
6. **Effective length hand-off.** Construct `t_eff_lens_weighted` in
   the pipeline; pass through to `AbundanceEstimator`. Bit-exact
   parity test when `mode == "uniform"`.
7. **Locus diagnostics.** Add `mrna_gdna_logratio_shift` and
   `gdna_prior_count_renormalized` columns to `loci.feather`.
8. **Internal mode plumbing.** Add the hidden enum and env-var. Wire Modes 0–4
  and 3a. `full_symmetric`/Mode 3b may fail fast until candidate-level mRNA
  footprint data is implemented.
9. **CI synthetic locus test.** §9.3.
10. **VCaP mode sweep.** §9.2. Apply gates §9.1.
11. **Synthetic 24-condition regression.** §9.4.
12. **Decide production mode.** If G1–G4 pass for Mode 3, ship Mode 3
    behind `auto`. If G4 fails, implement and test Mode 3b. If G2 fails for both,
    *do not ship*; iterate on the field-construction layer.
13. **Document.** Update `docs/calibration/gdna_exposure_plan_v3.md`
    status to "superseded by revised plan v2" and link this document.

---

## 12. Open Questions and Risks

### 12.1 What if `Y_r` is dominated by RNA at exon boundaries?

The conservative count `Y_r` in the EXON class is computed from
boundary-crossing fragments; these can include both spliced and
unspliced fragments. If a meaningful fraction is actually RNA, the
EXON `rho_hat` is inflated, which propagates into `Ā_locus` and
modulates Strategy-B prior renormalization. This is bounded by the
existing v6 strand-correction logic but is worth measuring per-mode.
Diagnostic: per-class `Y_r` totals already in `summary.json`.

### 12.2 Footprint vs midpoint integration

We use midpoint sampling of `A(x)` in v2 for parity with v3. If
validation shows boundary artifacts for long fragments spanning the
exon→intron edge, switch the numerator to a 2-point footprint average
(start, end) before considering a full integral. No change to the
denominator is needed (it already integrates over a window).

### 12.3 Loci with mixed-ref or multimapping units

The current CSR has one unit-level genomic midpoint, not a per-candidate or
per-hit coordinate. Therefore both gDNA and nRNA numerator weighting skip units
whose candidates span multiple references or whose fragment class is
multimapper. This is conservative and symmetric: it leaves the local exposure
factor at identity for both unspliced hypotheses. We count
`n_units_cross_ref_skipped` and `n_units_multimapper_skipped` for transparency.

If future profiling shows these skips dominate a real assay, the next step is a
native scoring extension that emits candidate-level or hit-level genomic
positions and applies exposure inside the multimapper gDNA LSE.

### 12.4 The unified smoothed field is still deferred

Class-specific normalization is a bridge. The asymmetric
calibration of EXON vs INTRON `rho_hat` (boundary-crossing density vs
contained-fragment density) means that `∫ A(s) ds` crossing class
boundaries integrates inconsistent quantities. Symmetric application
absorbs most of this in cancellation (§3), but **only most**. The v3
of this document line of work is: replace per-class densities with a
single smoothed `coverage(x)` estimator that converts all evidence
sources onto one target. Tracked in `TODO.md`.

---

## 13. Bottom Line

The fair model is: any hypothesis whose molecule occupies the exposed
genomic location must be penalized by the same `log A(x)`. In v2
that is gDNA and synthetic nRNA. mRNA is approximated as already-at-1
on captured exons, with a measured guard rail (Mode 3b) so we catch
the approximation breaking. Everything else — the class-specific
`rho_ref` fix, the prior-count Strategy A/B choice, the deferred
unified field — is hygiene or follow-up.

The validation gates in §9.1 are the contract. If the symmetric
model cannot meet them on the VCaP truth case, we do not ship; we
revisit the field construction. If it can, we ship behind the
unchanged `--regional-exposure auto` user flag.
