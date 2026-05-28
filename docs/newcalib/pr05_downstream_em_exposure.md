# PR 05 - Downstream EM Exposure Normalization

## Goal

Fix the downstream sign and aggregation of exposure factors in locus EM.

The current code computes a regional exposure average and multiplies gDNA effective length by it.
Native EM then subtracts `log(gdna_eff_len)`, so high exposure lowers gDNA likelihood. This PR makes
high exposure increase gDNA competitiveness.

## Current State

In `src/rigel/calibration/prior.py`, `assemble_priors()` currently does:

```text
exposure_weight = bp_weighted_mean_exposure_over_blocks(...)
gdna_eff_len = max(unweighted * exposure_weight, 1.0)
gdna_em_exposure_weight = exposure_weight
```

In `src/rigel/native/em_solver.cpp`, the gDNA component receives:

```text
sub.log_eff_len[gdna_idx] = log(gdna_eff_len)
```

The EM likelihood contract is therefore:

```text
log P(read | gDNA component) includes -log(gdna_eff_len)
```

Multiplying `gdna_eff_len` by exposure has the wrong sign.

## New Contract

For the first production exposure model, exposure is a denominator-only visibility adjustment:

```text
omega_locus = exposure factor for the gDNA candidate window
gdna_eff_len_em = max(gdna_eff_len_unweighted / omega_locus, 1.0)
```

No per-unit `+log(omega)` term should be applied in the same model.

## Locus Exposure Aggregation

Replace `bp_weighted_mean_exposure_over_blocks()` for EM normalization with an FL-window-aware helper
that mirrors `gdna_eff_len_for_loci()`:

```text
gdna_exposure_factor_for_loci(loci, ref_lengths, gdna_fl, region_arrays, exposure)
```

The helper should compute:

```text
omega_locus = weighted_exposed_opportunity / unweighted_opportunity
```

where the opportunity windows are the same fragment-start windows used by `gdna_eff_len_for_loci()`.
This matters for long loci with small captured subregions.

If exact FL-window aggregation is too large for this PR, use the existing bp-weighted approximation
only behind a clearly named temporary helper and add a failing or skipped test documenting the exact
aggregation requirement. Do not keep the old helper name for the new contract.

## Code Changes

### `src/rigel/calibration/_exposure.py`

- Add the FL-window-aware locus exposure helper.
- Keep `gdna_eff_len_for_loci()` as the unweighted base opportunity helper.
- Delete or deprecate production use of `bp_weighted_mean_exposure_over_blocks()` for EM.

### `src/rigel/calibration/prior.py`

Change `assemble_priors()` to:

```text
unweighted = gdna_eff_len_for_loci(...)
omega_locus = gdna_exposure_factor_for_loci(...)
gdna_eff_len = max(unweighted / omega_locus, 1.0)
```

Rename diagnostics so sign is unambiguous:

```text
gdna_eff_len_unweighted
gdna_eff_len_em
gdna_exposure_factor
gdna_eff_len_adjustment_ratio = gdna_eff_len_em / gdna_eff_len_unweighted
```

No backward-compatible field names are required. If a transitional alias is needed inside one PR,
delete it before merging the full redesign.

### `src/rigel/pipeline.py` and `src/rigel/estimator.py`

- Update locus metadata and output schemas to the new names.
- Remove `gdna_em_exposure_weight` if it is only a duplicate of `omega_locus`.
- Ensure `loci.feather` makes the sign obvious: enriched loci should have high
  `gdna_exposure_factor` and low `gdna_eff_len_adjustment_ratio`.

### Native EM

No native EM math change is expected for the denominator-only model. The native solver already
subtracts log effective length.

## Tests

- Unit test with `unweighted = 100000` and `omega_locus = 50`: EM effective length is `2000`, not
  `5000000`.
- Unit test with `omega_locus = 0.1`: EM effective length is larger than unweighted.
- Locus metadata test checks diagnostic ratio equals `1 / omega_locus` within tolerance.
- Integration smoke test with a mock high-exposure gDNA locus shows gDNA posterior increases versus
  uniform exposure.
- Golden outputs updated only after the sign fix is confirmed intentional.

## Future Model Not in This PR

A fully normalized local-exposure likelihood would use:

```text
log_lik_unit += log(omega_unit)
denominator = integral of omega over candidate windows
```

That is a different representation. It must replace the denominator-only model and come with tests
showing the expected behavior in mixed-exposure loci. Do not stack it on top of this PR.

## Done Means

- High exposure decreases the gDNA denominator in native EM inputs.
- Low exposure increases the gDNA denominator.
- Output diagnostics make the sign inspectable.
- No production path multiplies gDNA effective length by exposure.
