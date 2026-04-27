# Rigel Streamlining & Performance Review

**Sample profiled:** VCAP RNA20M + DNA20M (`mctp_vcap_rna20m_dna20m`),
8 threads, full human index (457k transcripts).
**Pipeline state:** post-SRD v1, post-region removal (≥ v0.5).

---

## §0 Executive Summary

The two-stage architecture (C++ BAM scan → Python/C++ locus EM) is sound and
the calibration overhaul (SRD v1) leaves the codebase in much better shape than
the previous v5/region/Beta-Binomial era. The largest opportunities now are:

1. **Wall-time:** scan_and_buffer (58%) and fragment_router_scan (15%) together
   eat 73% of runtime. Both are I/O-shaped but have measurable Python-side fat
   that can be cut without algorithmic changes.
2. **Memory:** peak RSS 11.9 GB on a 32M-fragment library, with a 7.7 GB
   delta. The buffer (~4 GB in-memory + spilled chunks) and the partition step
   are the dominant consumers.
3. **Correctness/clarity:** several stale modules and v5-era docstrings are
   still present, and a handful of constants ("magic numbers") govern
   calibration / EM behavior without being exposed as config.
4. **EM:** SQUAREM dominates the EM stage at 83% of CPU; the implementation is
   already well-optimized (Kahan, parallel E-step). Wins here are about
   reducing iteration counts on the long-tail of small loci, not micro-opt.

Recommended ordering: **§5 dead-code purge first (zero risk, immediate
clarity)**, then **§4.1 fast-path chunk materialization (highest ROI on
wall-time)**, then **§4.3 fused calibration+scoring pass (deferred until
benchmarked)**.

---

## §1 Empirical Profile (dna20m, 8 threads)

```
Wall time:       346.90 s
Throughput:      190,583 frags/s
Buffered frags:  31,985,707
Multimapping:    275,868 molecules
Peak RSS:        11,947 MB   (Δ 7,669 MB above pre-quant baseline)
```

### Stage breakdown

| Stage                  |   Wall (s) |    %  | Notes |
|------------------------|-----------:|------:|-------|
| scan_and_buffer        |    199.585 | 57.5% | C++ BAM scan + train + spill |
| finalize_models        |      0.000 |  0.0% | no-op (folded into scan) |
| calibration (SRD v1)   |     11.313 |  3.3% | mixture EM over pool histogram |
| compute_geometry       |      0.026 |  0.0% | trivial |
| fragment_router_scan   |     50.558 | 14.6% | likelihoods + routing per chunk |
| build_loci             |      3.938 |  1.1% | UF over multimapper edges |
| partition              |     40.709 | 11.7% | scatter CSRs into per-locus arrays |
| eb_gdna_priors         |      3.181 |  0.9% | `compute_locus_priors_from_partitions` |
| locus_em               |     37.590 | 10.8% | parallel SQUAREM over all loci |
| **TOTAL**              | **346.900**|       | |

### RSS timeline

| Phase | RSS (MB) |
|-------|---------:|
| before                | 4,278 |
| after_scan            | 6,466 |
| after_calibration     | 6,720 |
| after_router_scan     | 9,684 |
| after_locus_em        | 9,997 |

### Locus EM concentration (per-locus stats from `locus_stats_default.json`)

The C++ EM stage is well parallelized: **CPU-time sum = 240 s** vs **wall = 37.6 s**
→ effective parallel speedup ≈ 6.4× on 8 threads (good). Within EM:

| Sub-stage     | CPU (s) |    %  |
|---------------|--------:|------:|
| extract       |   12.97 |  5.4% |
| bias          |    6.73 |  2.8% |
| build_ec      |    9.18 |  3.8% |
| warm_start    |    1.00 |  0.4% |
| **squarem**   | **200.14** | **83.2%** |
| assign        |   10.43 |  4.3% |

Cost concentration across 25,055 loci:

| Top-N loci | Cumulative % of EM CPU |
|-----------:|-----------------------:|
| 1          |  3.2% |
| 10         |  9.1% |
| 100        | 22.2% |
| 1,000      | 55.1% |
| 5,000      | 88.5% |

The single mega-locus (29,713 transcripts × 1.75M units) consumed only 7.7 s
total CPU (44 SQUAREM iters) — the parallel E-step (8 threads) and warm-start
keep it bounded. The next-largest loci are *medium*-sized loci that hit the
SQUAREM iteration cap with 100–200 iters and run single-threaded.
**The long tail is the budget here, not the mega-locus.**

---

## §2 Bottleneck Analysis

### B1. `scan_and_buffer` — 199.6 s (57.5%)

Wall-clock dominator. From log timestamps:

- 31,985,707 fragments scanned in ~140 s ⇒ ~228k frag/s sustained.
- 12 chunks of 1 M frag spilled to `/tmp` (~197 MB each, ~2.4 GB total spill).
- C++ resolver, exon-fit, strand training, and FL training all happen here.

Sub-bottleneck signals:
- Three pauses visible (~50 s, ~50 s, ~50 s gaps between spill bursts) suggest
  the consumer (Python `_on_chunk` → `_FinalizedChunk.from_raw`) is at times
  blocking the C++ producer. The chunk-materialization wrapper (review §1.2)
  is the prime suspect.
- Spill is triggered when in-memory exceeds the buffer cap; on this run we
  spilled 12 / 32 chunks. Spill itself is fast (sub-second per chunk).

### B2. `fragment_router_scan` — 50.6 s (14.6%)

Walks the (now-finalized) buffer chunk-by-chunk, calls
`FragmentScorer.score_chunk` (C++) for each. Throughput ≈ 632k frag/s — much
higher than B1, but still significant in absolute terms because the buffer is
re-walked. Together with B1, **the fragments are touched twice in Python
(chunk wrapping + scoring iteration)**, which is the largest single
streamlining opportunity.

### B3. `partition` — 40.7 s (11.7%)

`partition_and_free()` scatters global CSRs into per-locus partitions and
nulls the originals. RSS grows during scatter (visible in `after_router_scan`
9.7 GB → stays flat through `after_locus_em` 10.0 GB). The scatter is C++ but
with Python orchestration; the time is plausibly dominated by the scatter
itself, not the orchestration. **No obvious algorithmic improvement** without
restructuring how partitions are produced (e.g., emit per-locus during scoring,
not scatter post-hoc).

### B4. `locus_em` — 37.6 s (10.8%)

Already parallelized 6.4×; SQUAREM dominates within (83%). Long-tail
concentration: the bottom 80% of loci (~20k) cost only ~12% of EM CPU. The
top 1000 loci cost 55%. Streamlining within EM itself has diminishing
returns; the bigger lever is **convergence**:
- Mean SQUAREM iters across all loci ≈ 246,167 / 25,055 ≈ **9.8** (low).
- But the largest loci hit 100–200 (capped). A tighter warm-start or better
  step-size policy on those few hundred loci would matter most.

### B5. Memory — peak 11.9 GB, Δ 7.7 GB

| Consumer | Approx |
|----------|-------:|
| FragmentBuffer (32M frag in-memory portion) | ~4.0 GB |
| Per-locus partitions (post-scatter) | ~3.5 GB |
| Index + scorer caches | ~1.5 GB |
| EM scratch (per-thread Kahan, posteriors) | ~0.5 GB |

The 9.7 → 10.0 GB jump during EM shows partition arrays are NOT being freed
as loci finish — this is the largest memory footgun. Per-locus partitions
should be released as their EM completes.

---

## §3 Code Review Findings (from full structural sweep)

A complete fresh-eyes review by a separate exploration pass identified ~30
discrete items. The high-signal ones are summarized below; the full report
covered also stale jargon and minor refactors that are bundled into §5.

### 3.1 Dead modules / docstrings

- **[src/rigel/mappability.py](src/rigel/mappability.py)** — entire ~250-line
  module is dead post-region removal. `grep -r "mappability" src/ tests/`
  shows no production importer. **Delete.**
- Stale v5/Beta-Binomial/regional references in docstrings:
  - [src/rigel/config.py](src/rigel/config.py) `CalibrationConfig` ("No regional priors, no Beta-Binomial…" reads as negation list).
  - [src/rigel/locus.py](src/rigel/locus.py) "matches the v5 `gdna_prior_c_base` default".
  - [src/rigel/estimator.py](src/rigel/estimator.py) "Pass 1" terminology.
  - [src/rigel/index.py](src/rigel/index.py) `intervals.feather` mentioned but unused for calibration now.

### 3.2 Magic constants that should be config

- `_C_BASE_DEFAULT = 5.0` and `_PI_FLOOR = 1e-6` in
  [src/rigel/locus.py](src/rigel/locus.py). The code already comments
  "promote to CalibrationConfig if Phase 5 benchmarks show sensitivity" —
  benchmarks have run; either commit to making them config or delete the
  comment.

### 3.3 Untuned / unverified config knobs

- `BamScanConfig.n_decomp_threads` default 4 — never benchmarked, may be
  redundant if htslib auto-tunes. Either remove or document the tuning.
- The mega-locus heuristic (`work ≥ fair_share`) in
  [src/rigel/pipeline.py](src/rigel/pipeline.py) `_run_locus_em_partitioned`
  is implicit, not exposed, and never falls back to mega-mode in
  single-thread runs. The profile shows it correctly fired on the 1.75M-unit
  locus, but it's worth a `mega_locus_threshold_work` knob with logging.

### 3.4 Refactor opportunities (clarity, not perf)

- `FragmentScorer.from_models()` is a ~100-line factory — split into
  `_build_fl_lut()`, `_build_strand_context()`, `_build_penalty_dict()`
  for testability.
- Two nested closures inside `_run_locus_em_partitioned()` capture outer
  scope; promote to module-level `_build_locus_meta_from_partition` /
  `_call_batch_em` for unit-testability.
- Fragment-class label constants are echoed in C++ and Python; consolidate
  into a single source of truth.

### 3.5 Algorithmic opportunities (deferred)

- **Vectorize the outer loop** of `compute_locus_priors_from_partitions`
  using `np.add.reduceat` over a flattened global posterior array. Stage is
  only 3.2 s today, so this is a clarity win, not a wall-time win.
- **Push categorize_chunk into C++**: only justified if calibration becomes
  a bottleneck (currently 3.3%). Skip.

---

## §4 Prioritized Improvement Roadmap

Ordering combines profile evidence with code-review risk/impact estimates.

### P0 — High impact, low risk (do first)

| # | Item | Files | Expected win |
|---|------|-------|-------------:|
| P0-1 | **Delete `mappability.py`** + remove `intervals.feather` references in `index.py` | [mappability.py](src/rigel/mappability.py), [index.py](src/rigel/index.py) | Code clarity, ~250 LOC removed |
| P0-2 | **Strip stale v5/Beta-Binomial/regional/Pass-1 docstrings** | [config.py](src/rigel/config.py), [locus.py](src/rigel/locus.py), [estimator.py](src/rigel/estimator.py), [buffer.py](src/rigel/buffer.py) | Avoids future-reader confusion |
| P0-3 | **Fast-path chunk materialization**: skip `_FinalizedChunk.from_raw` for in-memory chunks; direct dict→C++ scorer | [buffer.py](src/rigel/buffer.py), [pipeline.py](src/rigel/pipeline.py) `_on_chunk` | **Targeting B1 + B2: ~5–10% wall (~15–30 s)** |
| P0-4 | **Free per-locus partitions as EM completes** (currently held until full EM finishes) | [pipeline.py](src/rigel/pipeline.py) `_run_locus_em_partitioned`, [locus_partition.py](src/rigel/locus_partition.py) | **−1 to −2 GB peak RSS** |

### P1 — Targeted perf (after P0)

| # | Item | Expected win |
|---|------|-------------:|
| P1-1 | **Fuse calibration + scoring buffer passes** so chunks are walked once. Requires honoring the fact that calibration must produce gDNA FL model before scoring; could be done by collecting calibration counts during scoring of the *first* in-memory portion, then finalizing. | **~10–15% wall** |
| P1-2 | **SQUAREM iteration cap tightening** for medium loci that hit 100–200 iters. Investigate convergence criterion (delta vs. ELBO trend) — top-1000 loci own 55% of EM CPU; halving their iter count saves ~10 s wall. | ~3–5% wall |
| P1-3 | **Promote `_C_BASE_DEFAULT` / `_PI_FLOOR` to `CalibrationConfig`** *or* delete the "promote later" comment and lock them. | Hygiene |
| P1-4 | **Audit `BamScanConfig.n_decomp_threads`** — benchmark or remove. | Hygiene |

### P2 — Refactor / testability

| # | Item |
|---|------|
| P2-1 | Split `FragmentScorer.from_models()` into typed helpers. |
| P2-2 | Promote nested closures in `_run_locus_em_partitioned` to module-level testable functions. |
| P2-3 | Single source of truth for fragment-class constants (C++ `constants.h` re-exported to Python). |
| P2-4 | Vectorize outer loop of `compute_locus_priors_from_partitions` (clarity, not perf). |

### P3 — Investigate / defer

| # | Item |
|---|------|
| P3-1 | Mega-locus heuristic: expose `mega_locus_threshold_work`; verify against wall-time on a sweep, not theoretical work. |
| P3-2 | Lazy `category_counts` in `CalibrationResult` (small memory, but cleaner API). |
| P3-3 | C++-side categorize_chunk port — only if calibration becomes >10% of wall time after P1-1 lands. |

---

## §5 Honest Concerns

1. **Partition memory accounting is opaque.** The 9.7 GB → 10.0 GB rise
   during EM tells us partitions are not being released per-locus. We should
   either confirm per-locus release happens via C++ scope and re-measure, or
   add explicit `gc.collect()` at locus boundaries inside
   `_run_locus_em_partitioned`. (Item P0-4.)

2. **Mega-locus heuristic is ad-hoc.** It fires correctly on this sample but
   has no exposed knob, no fallback in single-thread, and is justified by a
   simplistic `n_tx × n_units` cost proxy. Add a knob and log the decision.

3. **Scan stage variability.** Three multi-second gaps in the scan log
   suggest occasional consumer back-pressure (Python wrapping the chunk
   slower than C++ producing). If P0-3 lands and the gaps shrink, we have
   evidence the wrapping was the culprit; if not, it's the BAM I/O itself.

4. **The "promote if benchmarks show sensitivity" comment in `locus.py`
   has been there a while.** Either act on it or delete it; long-lived
   conditional TODOs erode trust in the codebase.

5. **The locus-level priors stage (3.2 s) is well-vectorized within each
   locus** but the per-locus Python loop costs ~125 µs × 25k loci = 3 s,
   essentially the entire stage. A vectorized rewrite would not save wall
   time; the value is in clarity, not perf.

---

## §6 Verification Pointers

- Profile artifacts: `/tmp/rigel_profile/dna20m/`
  - `profile_report.txt` (human-readable)
  - `profile_summary.json` (structured)
  - `locus_stats_default.json` (per-locus EM stats)
  - `memory_timeline_default.csv` (RSS samples, 2613 rows)
- After landing P0-3 + P0-4, re-run with the same command:
  ```bash
  python scripts/profiling/profiler.py \
    --bam /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna20m/rigel/annotated.bam \
    --index /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index \
    --outdir /tmp/rigel_profile/dna20m_postp0 \
    --stages --threads 8 --tmpdir /tmp
  ```
  Targets: scan_and_buffer + fragment_router_scan ≤ 220 s (currently 250 s),
  peak RSS ≤ 10 GB (currently 11.9 GB).
