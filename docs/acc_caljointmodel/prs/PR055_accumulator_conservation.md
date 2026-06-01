# PR 5.5 — Accumulator: conserve fragment mass for encompassed regions

A small, dedicated bug fix that lands **before PR 6** (which consumes the
calibration mass + exposure this bug inflates).

## The bug

The fractional accumulator's boundary attribution (spec
`docs/accumulator/00_design.md` §4.3) over-counts the mass of any region a
fragment **encompasses** — i.e. fully traverses, overlapping *both* the region's
boundaries. The interior region's slice mass is deposited at **full value on each
side** instead of split, so total deposited mass exceeds 1.0:

| fragment | region | should deposit | actually deposited |
|---|---|---|---|
| `[50,250)` over `[0,100)[100,200)[200,300)` | `[100,200)` (interior) | 0.25 + 0.25 = 0.50 | 0.50 + 0.50 = **1.00** |

The old tests *asserted* the bug: `test_t6` expected `total_mass == 1.5`,
`test_t5` expected `1.0 + 80/L`. It is consistent across spec + Python reference
(`_accumulator_reference.py`) + C++ (`accumulator.cpp`) + substrate (no downstream
halving), so it is a spec-level oversight, not a regression.

It breaks the exposure's `∫overlap/ℓ\,dx = L` mass-conservation identity and
inflates per-region gDNA/RNA mass for short exons/introns spanned by contiguous
fragments — exactly the gDNA signal calibration estimates.

## The fix (per `D-ACC`)

A fragment encompasses a region ⇔ it overlaps both region boundaries (an
*interior* slice in the per-region slice decomposition). Such a region splits its
mass **50/50**: half to the RIGHT side of its LEFT boundary (`mass_right`), half
to the LEFT side of its RIGHT boundary (`mass_left`). End slices (the fragment's
first/last region) keep full mass on their single crossed side.

**Implementation** — rewrite the crossing path from *per-adjacent-pair* to
*per-slice*: each slice deposits `slice_mass / n_cross` to each side it crosses,
where `n_cross` = (crosses-left? 1:0) + (crosses-right? 1:0) ∈ {1, 2}. This is
exactly the half-split for interior slices and a no-op for end slices / simple
2-region crossings. It uniformly handles the contiguous span *and* a covered
spliced middle exon (both are "encompassed").

- **Flux is unchanged** — `flux_left/flux_right` stay the integer per-side
  crossing-event count (1 per crossing per side). Only `mass_left/mass_right`
  is split. "Power from flux, density from mass."
- **Downstream is unchanged** — same dtypes, same interfaces; only mass *values*
  become correct. No Python calibration code changes.

## Scope of change

| File | Change |
|---|---|
| `src/rigel/native/calibration/accumulator.cpp` | crossing loop → per-slice with `n_cross` division |
| `tests/native/_accumulator_reference.py` | same (the byte-for-byte reference) |
| `docs/accumulator/00_design.md` §4.3 / §6 | document the encompassing/interior-slice 50/50 rule + restate per-fragment mass = 1 |
| `tests/native/test_accumulator_spec.py` | T5 (`b1.mass_right`, `b2.mass_left` → half; total → 1.0), T6 (same; total → 1.0); add a spanning-fragment conservation test |

Unaffected: T1–T4, T7, T8, T10–T12 and the FL-pool tests (no interior slice, or a
different quantity). The FL-pool deposit (per-slice, full mass) already conserves
per fragment and is left as-is.

## Verification

Rebuild C++ (`pip install --no-build-isolation -e .`), run
`pytest tests/native/ -v` (reference == native byte-for-byte; new conservation
test), then the full suite + `ruff`.

## Rollback

Revert the per-slice loop in `accumulator.cpp` + reference, the spec §4.3 text,
and the T5/T6 assertions. Self-contained.
