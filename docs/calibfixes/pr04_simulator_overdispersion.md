# PR 4: Simulator Overdispersion Knobs

## Goal

Make the synthetic suite stress the uncertainty models that PR 1 and PR 2 rely on. The current suite can produce near-binomial gDNA strand balance (`kappa_d` effectively at the hard maximum), which is useful but too clean. This PR adds or wires simulation knobs for gDNA strand-balance overdispersion, RNA minor-orientation overdispersion, and optional region-level strand-bias heterogeneity.

This PR does not change calibration behavior directly. It changes simulation configuration, truth/manifest reporting, and tests.

## Current code status

`src/rigel/sim/reads.py` already has `GDNAConfig.strand_kappa` and region-level gDNA strand overdispersion for the smaller read simulator. The whole-genome simulator and suite YAML parser do not appear to expose an equivalent field yet. Therefore this PR should first wire the existing idea through the whole-genome/suite path, then add the missing RNA and heterogeneity knobs.

## Non-goals

- Do not change default simulator behavior. Existing configs should remain near-binomial unless they opt into overdispersion.
- Do not make calibration tests depend on stochastic high-variance defaults.
- Do not alter post-capture truth semantics.
- Do not use simulator overdispersion to tune calibration constants.

## Files to edit

| File | Purpose of edit |
|---|---|
| `src/rigel/sim/reads.py` | Add RNA minor-orientation and optional region bias knobs to the small simulator; keep existing `GDNAConfig.strand_kappa`. |
| `src/rigel/sim/whole_genome.py` | Add fields to `GDNASimConfig` and `SimulationParams`/strand config; implement overdispersed strand generation in whole-genome simulation. |
| `src/rigel/sim/suite.py` | Parse suite YAML fields and propagate them to whole-genome configs. |
| `src/rigel/sim/manifest.py` or manifest writing code | Record overdispersion settings per condition. |
| `src/rigel/sim/truth.py` | Ensure post-capture truth and FL truth remain unchanged except for metadata if needed. |
| `tests/test_sim.py` | Unit tests for small simulator strand overdispersion. |
| `tests/test_sim_capture.py` | YAML parsing and suite config tests for new fields. |
| `tests/test_sim_truth.py` | Non-regression for post-capture truth. |
| `docs/SIMULATOR.md` or `docs/calibfixes` follow-up note | Document the new knobs and recommended stress matrices. |

## Configuration additions

### gDNA strand-balance overdispersion

Add to whole-genome `GDNASimConfig`:

```python
strand_kappa: float | None = None
strand_region_size: int | None = None
```

Semantics:

- `None` or `<= 0`: current behavior, independent 50/50 gDNA strands per fragment.
- `strand_kappa > 0`: each region draws `p_plus ~ Beta(kappa/2, kappa/2)` and gDNA fragments in that region use that plus-strand probability.
- `strand_region_size is None`: use annotation-derived fine-region boundaries where available, matching the small simulator idea.
- `strand_region_size > 0`: use fixed genomic tiles of this size for references without useful annotation or for controlled stress tests.

YAML:

```yaml
gdna:
  strand_kappa: 50
  strand_region_size: 5000
```

### RNA minor-orientation overdispersion

Add to the RNA simulation config. If `SimulationParams` is the only strand-specificity home, add fields there:

```python
rna_orientation_kappa: float | None = None
rna_orientation_region_size: int | None = None
```

Semantics:

- Existing `strand_specificity` remains the mean probability that RNA preserves the expected read orientation.
- `None` or `<= 0`: current behavior, each RNA fragment flips with probability `1 - strand_specificity` independently.
- `rna_orientation_kappa > 0`: local minor-orientation rate varies by region/gene/transcript around the global mean.

Use a Beta distribution parameterized by mean `q = min(strand_specificity, 1 - strand_specificity)` and concentration `kappa`:

```text
q_region ~ Beta(q * kappa, (1 - q) * kappa)
```

Then for each RNA fragment in that region:

```text
flip ~ Bernoulli(q_region)
```

For unstranded `strand_specificity = 0.5`, this gives `q = 0.5`; the simulator remains unstranded on average but can have local imbalance when `kappa` is small. For perfectly stranded `strand_specificity = 1.0`, `q = 0`; handle this as a point mass at zero flips unless a separate explicit error floor is provided.

YAML:

```yaml
simulation:
  rna_orientation_kappa: 200
  rna_orientation_region_size: 5000
```

### Optional region-level strand-bias heterogeneity

Add only if it can be implemented without disrupting the first two knobs:

```python
strand_bias_mode: Literal["none", "gdna", "rna", "both"] = "none"
```

This is mostly a convenience switch for suite generation. The underlying fields above should remain explicit so benchmark configs are reproducible.

## Region assignment recipe

Use one helper for both gDNA and RNA local strand parameters:

```python
def build_strand_bias_regions(
    *,
    ref_lengths: Mapping[str, int],
    transcripts: Sequence[Transcript],
    region_size: int | None,
) -> dict[str, np.ndarray]:
    ...
```

Rules:

1. If `region_size` is provided, boundaries are `0, region_size, 2 * region_size, ..., ref_len`.
2. Else, use annotation-derived boundaries: transcript starts/ends and exon starts/ends, plus `0` and `ref_len`.
3. Always include reference end.
4. Boundaries must be sorted unique `int64` arrays.
5. Region lookup uses `np.searchsorted(boundaries, starts, side="right") - 1`, then clipping.

Small simulator `ReadSimulator._init_strand_regions()` can either keep its current implementation or call this helper.

## gDNA implementation recipe

### Small simulator

`reads.py` already has:

```python
GDNAConfig.strand_kappa
ReadSimulator._init_strand_regions()
ReadSimulator._gen_reads_from_genome()
```

Tasks:

1. Add `strand_region_size` if needed.
2. Reuse the shared region helper.
3. Add tests proving existing default remains 50/50 independent and `strand_kappa` creates correlated local strand ratios.

### Whole-genome simulator

In `whole_genome.py`, update `GDNASimConfig` and the gDNA write path.

Current `_write_gdna_chunk()` uses:

```python
chrom_strands = rng.integers(0, 2, size=count)
```

Replace with a helper:

```python
chrom_strands = self._sample_gdna_strands(ref, starts, count)
```

The helper should:

- return independent 0/1 strands when `gdna.strand_kappa` is unset;
- lazily initialize per-reference region boundaries and `p_plus` arrays when `strand_kappa > 0`;
- sample `is_negative = rng.random(count) >= p_plus_for_start` so `0=+`, `1=-`, matching existing qname conventions;
- use the same RNG seed path as the rest of the simulator for reproducibility.

## RNA implementation recipe

Find mature and nascent RNA read-generation paths in `whole_genome.py` and the small simulator. In `reads.py`, the current pattern is:

```python
ss = self.config.strand_specificity
if ss < 1.0:
    flip_mask = rng.random(count) >= ss
else:
    flip_mask = None
```

Replace with a helper:

```python
flip_mask = self._sample_rna_flip_mask(
    ref=t.ref,
    starts=genomic_or_transcript_starts,
    count=count,
    strand_specificity=self.config.strand_specificity,
    is_nrna=is_nrna,
)
```

For transcript-coordinate starts, map to genomic starts if possible. If mapping is expensive for spliced mature RNA, first implementation may assign by transcript/gene rather than exact genomic start:

- build a per-transcript `q_minor` from the transcript midpoint region;
- use that `q_minor` for all fragments from the transcript;
- document this approximation and prefer exact genomic mapping in a later cleanup.

Helper semantics:

```text
if rna_orientation_kappa is None or <= 0:
    flip probability = 1 - strand_specificity
else:
    q_mean = min(strand_specificity, 1 - strand_specificity)
    q_local = local Beta(q_mean * kappa, (1 - q_mean) * kappa)
    flip probability = q_local if the protocol's major orientation is preserved by default
```

Because the simulator currently treats `strand_specificity` as probability of preserving expected orientation, the flip probability is `1 - strand_specificity`. For `strand_specificity < 0.5`, either reject the config or define it explicitly. Current suite values are `0.50` and `0.99`, so keep validation simple:

```text
0.5 <= strand_specificity <= 1.0
```

## YAML and suite parsing

Update both config loaders:

- `src/rigel/sim/whole_genome.py::parse_yaml_config()`
- `src/rigel/sim/suite.py::_load_suite_config()`

Accepted YAML:

```yaml
gdna:
  rates: [0.0, 1.0]
  rate_labels: [none, high]
  strand_kappa: 100
  strand_region_size: 5000

simulation:
  n_rna_fragments: 100000
  rna_orientation_kappa: 200
  rna_orientation_region_size: 5000
```

Suite CLI overrides can be added later if not needed immediately. YAML support is the requirement for benchmark reproducibility.

## Manifest and metadata

Each condition manifest entry should include:

```json
{
  "gdna_strand_kappa": 100,
  "gdna_strand_region_size": 5000,
  "rna_orientation_kappa": 200,
  "rna_orientation_region_size": 5000
}
```

These fields do not change truth abundances directly. They explain the generated read-orientation distribution and should be present in the manifest for auditability.

## Recommended stress matrix

Keep the existing near-binomial suite as baseline. Add a second overdispersed suite or config variant:

| Knob | Baseline | Moderate stress | Strong stress |
|---|---:|---:|---:|
| `gdna.strand_kappa` | null | 200 | 20 |
| `simulation.rna_orientation_kappa` | null | 500 | 50 |
| `strand_region_size` | annotation | 10000 | 2000 |

Do not use strong stress as the first acceptance gate for PR 1/PR 2. Use it after the baseline eight-condition suite passes.

## Tests

### Small simulator tests in `tests/test_sim.py`

Add tests for:

1. Default `GDNAConfig.strand_kappa=None` keeps both strand labels present over enough fragments.
2. Small `strand_kappa` creates region-correlated gDNA strand ratios. Use a deterministic seed and assert the per-region variance is larger than the independent baseline.
3. `rna_orientation_kappa=None` reproduces the expected binomial flip rate within a loose tolerance.
4. Small `rna_orientation_kappa` creates transcript/region-level variation in flip rates.
5. `strand_specificity=1.0` with RNA orientation kappa still produces no flips unless an explicit error floor is added.

### Whole-genome and YAML tests in `tests/test_sim_capture.py`

Add tests for:

1. `parse_yaml_config()` reads `gdna.strand_kappa` and `gdna.strand_region_size`.
2. `parse_yaml_config()` reads `simulation.rna_orientation_kappa` and `simulation.rna_orientation_region_size`.
3. `_load_suite_config()` forwards the same fields.
4. Manifest entries include all four overdispersion metadata fields.

### Truth non-regression in `tests/test_sim_truth.py`

Add or update a test proving post-capture abundance truth and FL truth are unchanged by orientation overdispersion when the same molecular/capture config is used. Read orientation changes alignment/calibration stress, not molecular truth.

## Validation commands

```bash
conda activate rigel && ruff check src/rigel/sim tests/test_sim.py tests/test_sim_capture.py tests/test_sim_truth.py
conda activate rigel && pytest tests/test_sim.py tests/test_sim_capture.py tests/test_sim_truth.py -v
```

Run a smoke suite with overdispersion:

```bash
conda activate rigel && python scripts/sim/simulate_suite.py -c scripts/sim/configs/<overdispersed_smoke>.yaml
```

Use an actual config path once the PR adds it. Do not run large full-scale simulations as part of the unit-test gate.

## Benchmark gate

After PR 4, run calibration benchmarks in two phases:

1. Baseline near-binomial suite: must match PR 1 + PR 2 acceptance behavior.
2. Moderate overdispersion suite: calibration diagnostics should report finite `kappa_d`, wider reliability uncertainty, and no cliff-like source mass changes.

Expected diagnostic changes:

- `kappa_d` should no longer saturate at `1e6` in moderately overdispersed high-gDNA conditions.
- PR 1 reliability should decrease smoothly in noisy strand regions.
- No-gDNA scenarios should not learn large exposure from RNA orientation noise.

## Review checklist

- [ ] Defaults reproduce existing simulation behavior.
- [ ] Whole-genome simulator supports gDNA strand overdispersion.
- [ ] RNA minor-orientation overdispersion is independent of molecular truth.
- [ ] YAML and manifest record every new knob.
- [ ] Tests verify deterministic behavior under fixed seeds.
- [ ] Overdispersed stress configs are additive; they do not replace the baseline suite.
