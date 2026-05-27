# Hybrid Capture Simulation Model

## Goal

Hybrid capture changes the fragment sampling distribution. Standard RNA-seq
simulation samples uniformly over legal fragment starts on each transcript or
genomic reference. Capture simulation needs a non-uniform landscape where
fragment starts overlapping bait/probe sequence are much more likely, while
off-target fragments remain possible.

This implementation is intentionally simple and fast enough for whole-genome
benchmarks.

## Probe Inputs

Two probe formats are supported by `rigel.sim.capture`:

1. Transcript-coordinate TSV: `transcript_id, start, end`, with an optional
   header. Coordinates are 0-based half-open transcript coordinates.
2. BED12: genomic 0-based half-open coordinates. Split BED12 blocks represent
   probes that span exon-exon junctions or other multi-block genomic layouts.

Transcript-coordinate probes are projected back to genomic blocks using the
transcript exon map. BED12 probes are projected onto any compatible transcript
whose exons fully contain the probe blocks.

When `scripts/sim/simulate_suite.py` generates random capture probes, it writes
both files:

- `reference/capture_probes_<group>.tsv` for transcript-design provenance.
- `reference/capture_probes_<group>.bed` for simulator input and genomic BED12
   inspection. Probes spanning exon-exon junctions are represented as multi-block
   BED12 records.

Generated probe design remains transcript-targeted, but the generated BED12 is
the capture model input. A probe designed against one isoform can therefore
enrich every compatible transcript that contains the same genomic probe blocks.

The suite can also use a provided probe panel directly. Supplying `probes:` on a
capture scenario bypasses random probe generation for that scenario and passes
the panel straight to the capture sampler:

```yaml
capture:
   configs:
      - label: off
         enabled: false
      - label: panel
         probes: /path/to/probes.bed
         format: bed12
         binding_per_base: 10.0
```

If `capture.probes` is supplied without a `capture.configs` list, the suite runs
a single capture-on scenario using that panel. Provided panels may be BED12 or
transcript-coordinate TSV; use `format: bed12`, `format: transcript`, or leave
`format: auto` for detection.

The random generated-probe designer uses a deterministic greedy heuristic:

1. Group eligible transcripts by gene.
2. Randomly select captured gene groups using `capture_fraction`; every eligible
   isoform in a selected gene enters the captured design pool, and every isoform
   in an unselected gene stays out of that pool.
3. Sort selected transcripts by mature-RNA molecular abundance, high to low.
4. Before designing probes for each transcript, project already-generated
    genomic probes back onto that transcript and mask those transcript-space
    intervals.
5. Tile new probes only in the remaining open transcript-space intervals.

The requested probe density is a soft target under masking. Shared exons or
short open intervals can make the realized number of new probes lower than the
nominal tiling density for a selected transcript.

Probe design uses the pre-capture mature-RNA molecular abundances assigned by
the suite abundance model. Post-capture abundances are not known until reads are
sampled through the full capture landscape.

Within a suite, generated probe sets are shared across compatible capture
scenarios. The probe-set key is the geometry of the design: capture fraction,
probe length, and probe density. Runtime weighting parameters such as
`binding_per_base`, `off_target_weight`, `gdna_split_penalty`, and `min_overlap`
reuse the same probe file so binding-energy sweeps remain head-to-head
comparable. Changing probe density, length, or captured-gene fraction
creates a separate probe set.

## Weight Model

For a fragment start `s` with fragment interval `[s, s + L)`, the sampling
weight is:

```text
w(s) = off_target_weight + binding_per_base * sum_i scale_i * overlap_bp(fragment, probe_i)
```

The overlap term is linear in hybridized bases. This is a coarse proxy for the
fact that longer contiguous probe-template hybridization is more stable. A more
thermodynamic Boltzmann-style model would exponentiate an energy term, but the
linear model has two practical advantages for Rigel simulation:

1. It is easy to reason about and tune.
2. It is exactly sparse-sampleable without materializing chromosome-sized
   probability arrays.

With `off_target_weight=1` and `binding_per_base=10`, a full 120 bp probe hit
has weight `1201` relative to weight `1` off target. A 10 bp overlap has weight
`101`. This is strong enrichment without making off-target sampling impossible.

## Transcript Choice

Capture affects transcript choice through a capture-aware effective length:

```text
Z_t(L) = sum_s w_t,L(s)
P(t | L, pool) proportional to abundance_t * Z_t(L)
```

Without probes, `Z_t(L)` reduces to the ordinary legal start count
`max(length_t - L + 1, 0)` when `off_target_weight=1`.

## Start Sampling

After a template is selected, fragment starts are sampled from `w(s)`. The
sampler uses an exact mixture for the linear model:

1. Baseline component: sample uniformly over all legal starts.
2. Probe component: choose a probe interval proportional to its total overlap
   mass, then sample a local start proportional to overlap with that interval.

The local start window for one probe is only about `fragment_length + probe_len`,
so gDNA capture can be sampled over full chromosomes without allocating a full
per-base or per-start chromosome array.

## Split Probe Behavior

Mature RNA gets the full contiguous transcript-coordinate probe. When the same
probe is projected to genomic DNA or pre-mRNA and splits across exon-exon
junctions, each genomic block is downweighted by `gdna_split_penalty`.

This captures the main biological asymmetry:

- mRNA can bind the full spliced probe sequence contiguously.
- gDNA and pre-mRNA contain intronic sequence between split probe blocks, so a
  junction-spanning probe has much lower effective binding potential.

For example, a transcript probe `[80, 140)` crossing a two-exon junction may map
to genomic blocks `[180, 200)` and `[400, 440)`. With `gdna_split_penalty=0.2`,
a gDNA fragment overlapping all 60 probe bases receives only `12` effective
hybridized bases, while the mature RNA fragment receives `60`.

## YAML Example

```yaml
capture:
  probes: probes.tsv
  format: transcript     # auto, transcript, or bed12
  off_target_weight: 1.0
  binding_per_base: 10.0
  gdna_split_penalty: 0.2
  min_overlap: 1
```

These settings are recorded in `manifest.json` under the `capture` key.
For generated probe suites, each condition also records `capture_probe_tsv` and
`capture_probe_bed` paths. The `capture_config.probes` field points to the BED12
file used by the simulator.

Each condition records both pre-capture and post-capture truth artifacts:

- `pre_capture_abundances` points to the molecular abundance table used before
   capture sampling.
- `post_capture_abundances` points to the empirical observed-fragment abundance
   table after capture sampling.
- `post_capture_fragment_lengths` points to the empirical post-capture fragment
   length distribution.

The post-capture abundance table includes explicit `pre_capture_*` molecular
columns and `post_capture_*` observed-fragment columns. Legacy-compatible
`mrna_abundance`, `nrna_abundance`, and `total_rna` columns contain the
post-capture observed fragment counts used for evaluation.

## Current Scope

The simulator applies capture weighting to mRNA, nRNA/pre-mRNA, and gDNA pools.
It does not yet model sequence-specific GC/Tm effects, probe depletion,
PCR/sequencing bias, or explicit post-capture molecule competition. Those can be
layered on later if benchmarking shows the simple overlap model is not enough.