# The calibration oracle & benchmarking — the ONE trusted path

This is the single, authoritative reference for how we debug, test, and benchmark the gDNA-vs-RNA
calibration and the per-locus EM. There is **one** oracle, **one** metric definition, and a small fixed set
of tools. If a diagnostic disagrees with what is written here, the diagnostic is wrong.

## The oracle — `scripts/debug/oracle.py`

**Principle: the oracle IS the production accumulator, partitioned by true fragment origin.** We split the
sim BAM into `gdna` / `mrna` / `nrna` by read-name origin (`rigel.sim.read_name.parse_origin`), run the
*same* production scanner+accumulator on each partition, and assert the partitions sum to the full payload.
Because the accumulator deposits each fragment independently, this **sum-to-full identity proves** the
partition is the production payload split by origin — no reimplementation, nothing to get subtly wrong.

`OracleTruth.from_bam(...)` runs four hard validation gates (raises on any failure — a silently-wrong
oracle can never be returned):

1. **Sum-to-full (integer channels):** `Σ_origin region_contained == full` exactly.
2. **Sum-to-full (flux/mass):** `boundary_flux_*` exact; `boundary_mass_*` within float32 rounding.
3. **gDNA is never spliced:** the `gdna` partition has zero spliced (ch2/3) contained mass.
4. **Every fragment accounted for:** every read is written to exactly one partition; totals reconcile.

CI: `tests/calibration/test_oracle.py` builds a tiny gDNA+mRNA+nascent scenario and asserts these hold.

### Accumulator channels (the basis everything is measured on)

`region_contained[R,4]` and `boundary_mass_{left,right}[B,4]`:

| ch | meaning |
|----|---------|
| 0 | unspliced, genome + |
| 1 | unspliced, genome − |
| 2 | spliced, sense (motif-relative) |
| 3 | spliced, antisense |

- **The gDNA-vs-RNA deconvolution (calibration) is over the UNSPLICED channels (0,1).**
- **Spliced (2,3) is mature RNA that NEVER competes with gDNA** — the EM separates it, and
  `result.py` withholds `mass_rna_spliced` from the RNA prior (`a_g + a_r` is the *unspliced* competing
  mass). Contained fragments live in `region_contained`; crossing/spliced fragments live in `boundary_mass`.
- **Intergenic gDNA is held out of the EM** (known DNA, no transcript there) — it is not a source of
  competition or error.

### Two truth perspectives (both supported)

- **All fragments accounted for** — every fragment lands in exactly one origin partition and one
  region/boundary node. Use this to know where every fragment goes.
- **Competing basis** — restrict to fragments that actually compete: *unspliced* gДNA-vs-RNA in genic
  regions (exclude intergenic gDNA and spliced mature). This lowers the denominator (e.g. below the raw
  library gDNA count) and is the honest basis for "calibration/EM deconvolution error."

## The retired oracle — do NOT resurrect

The deleted `_metrics.oracle_node_masses` (and its scripts `oracle_3pool.py`, `oracle_levered.py`,
`oracle_rna_dist.py`, `oracle_lever_test.py`, `oracle_calibration.py`, `ambig_locus_em_proof.py`,
`capon_dist.py`) deposited **whole fragments by span, with no intron-cutting** — an **incompatible basis**
with the accumulator the calibration actually consumes (per-base coverage, introns cut). It disagreed with
the accumulator on library-total gDNA (its 1.66M "contained" vs the accumulator's real mass) and wildly per
region (it reported ~0 RNA in high-expression exons that in fact hold real unspliced exon-body mRNA).
Overriding calibration masses with it **starved the RNA prior** via wrong-basis "0 RNA in exons", producing
a spurious 157k "fix" of the capON leak. This confounded earlier conclusions. Any oracle must be validated
by sum-to-full against the production accumulator (as `oracle.py` is) before it is trusted.

## The metric — soft 3-pool net surplus

- **Metric:** per-pool net surplus `assigned − true` for **gDNA / nascent / mature** (the soft EM counts).
  Sensitive to a calibration-prior shift, unlike the hard-label `net_flow` — **never gate on `net_flow`**
  (it is byte-identical under changes that move the soft pools by tens of thousands of fragments).
- **The nascent siphon = `estimator.nrna_em_count`** = EM mass on the **SYNTHETIC** nascent-shadow spans
  only. On an `nrna_none` scenario (true nascent = 0) any nonzero value is pure siphon. Three traps that
  give WRONG siphon numbers:
  1. `get_counts_df(index)["count"]` is `np.where(is_synthetic, 0, ...)` — it ZEROES the synthetic shadows.
  2. Filtering by `is_nrna` conflates real single-exon transcripts (legitimately `is_nrna=True`, both
     mature and nascent) with the synthetic shadows — their mass is real, not siphon.
  3. `net_flow` ZF/ZT hard labels are per-fragment MAP labels, not the EM abundance.
- Single-exon `is_nrna` transcripts count in the **mature** pool (`mrna_total`); the siphon is
  synthetic-only. The production `summary.json` `quantification` block already reports these correctly
  (`nrna_total == nrna_em_count`).

## The tools (the whole set)

| tool | purpose |
|------|---------|
| `scripts/debug/oracle.py` | THE oracle: validated per-origin truth + `override_masses()` (perfect-calibration lever) + per-region true `f_g`. |
| `scripts/debug/oracle_reattribute.py` | baseline vs oracle-override → EM → 3-pool: is a leak calibration error or downstream (EM/eff-length)? |
| `scripts/debug/pass_trace.py` | per-node calibration stage trace (init → strand → +prior → +messages → +KDE) vs oracle, on the **accumulator basis** — correct for localizing *where in calibration* an error enters. |
| `scripts/debug/benchmark_ab_report.py` | end-to-end A/B across the 16-condition suite → soft 3-pool + transcript accuracy JSON (`benchmark_ab_render.py` renders it). Arms may vary env (e.g. `RIGEL_MSG_MODE`). |

Suites: `quick_3to1_5mb` (16 conditions, easy), `ambig_dense_10mb` (AMBIG, hard). Run inside the `rigel`
conda env; `OMP_NUM_THREADS=1` for determinism.

## Worked conclusion (why this doc exists)

The capON gDNA→RNA leak was mis-attributed to calibration error by the retired oracle. On the validated
basis: the calibration's per-region gDNA fraction is accurate to <1% (mass-weighted |Δf_g| ≈ 0.007), and
feeding *perfect* calibration masses (`oracle_reattribute`) barely moves the leak (−112k → −106k). So the
capON leak is **downstream of an accurate calibration** — the EM assignment and/or the gDNA eff-length
translation — not the deconvolution.
