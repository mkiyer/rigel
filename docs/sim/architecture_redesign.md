# `rigel.sim` architecture redesign

**Status:** design for review (2026-06-07). Supersedes `suite_driver_unification_plan.md` (the
suite-driver unification is now Phase 4 of this larger plan). Goal: a clean, layered simulator
package that removes the real duplication and untangles the conflated modules — **with zero
behaviour change** (the simulator generates our test/benchmark ground truth; output must be
byte-identical across the refactor).

## 1. Current state — 16 modules, 11.3k lines, organically grown

Two simulation lineages grew in parallel and accreted duplication:

- **Scenario lineage** (small, single-condition, in-memory): `scenario.py` → `genome`/`annotation`
  + `reads.ReadSimulator` + `oracle_bam` + `truth`. Heavily used by unit/integration tests.
- **Suite lineage** (large, multi-condition, sharded): `suite.py` + `whole_genome.py`
  (`WholeGenomeSimulator` + `run_simulation`) + `synthetic_genome` + `manifest` + `truth`.

The duplication / tangle map (from the audit):

| # | Duplication / smell | Where |
|---|---|---|
| D1 | **Condition-grid orchestration** (nrna×gdna×od×ss×capture loop + manifest/truth) | `suite.py:1155–1313` ≈ `whole_genome.py:2055–2250` |
| D2 | **Fragment-length sampling** (truncated-normal, identical) | `reads.py:340–362` ≈ `whole_genome.py:1072–1088` |
| D3 | **gDNA strand-overdispersion region** Beta partition + searchsorted | `reads.py:272–294` ≈ `whole_genome.py:1096–1117` |
| D4 | **Sim config dataclasses** (same concepts, two definitions) | `reads.GDNAConfig`/`ReadSimConfig` ≈ `whole_genome.GDNASimConfig`/`SimulationParams` |
| D5 | **Interval math** (`merge_intervals`, `project_genomic_block_to_transcript`) | `capture.py:641,668` ≈ `suite.py:409,442` |
| D6 | **Splice-motif injection** (same logic, str-array vs MutableGenome) | `synthetic_genome.py:104` ≈ `annotation.py:165` |
| T1 | `truth.py` conflates **parse-origin** + **truth-write** + **counting** | `truth.py` |
| T2 | `capture.py` conflates **config** + **sampler** + **probe-loading**; probe **design** is off in `suite.py` | `capture.py`, `suite.py:375–537` |
| S1 | **Stale** `analyze_calibration`/`main` (old `summary.json` schema) amid live confusion helpers | `analysis.py` |

No circular dependencies today; `manifest.py` is the model of cohesion to aim for.

## 2. Design principles

1. **Behaviour-preserving.** Every phase is a pure move/extract; verified by the existing test
   suite **plus** a before/after suite-output diff (manifest + truth + BAM counts, pinned seeds).
2. **One concept, one home.** A shared algorithm (fragment sampling, interval math, splice motif,
   origin parsing) lives in exactly one module that others import.
3. **Separate data from behaviour.** Config dataclasses are tiny, serializable, dependency-free;
   runtime engines/samplers live apart from their config.
4. **Keep two engines, share their guts.** The in-memory `ReadSimulator` and the vectorized
   `WholeGenomeSimulator` serve genuinely different needs (audit: merging is *not* advised — string
   generator vs numpy/pysam/parallel). We do **not** merge them; we extract their duplicated
   internals (D2/D3/D4) into shared modules both consume.
5. **Flat-but-clean over deep nesting.** Award-winning ≠ deeply subpackaged. Prefer focused
   top-level modules with single responsibilities; one small `capture/` subpackage where the
   three-concern split clearly warrants it. Preserve the public `rigel.sim` API via `__init__`
   re-exports so no external caller breaks.

## 3. Target module layout

```
rigel/sim/
  # ── primitives (no sim-internal deps) ─────────────────────────────
  genome.py            MutableGenome, reverse_complement                      (keep)
  intervals.py    NEW  merge_intervals, project_genomic_block_to_transcript,  (D5)
                       clip_interval — the shared interval math
  splice_motif.py NEW  inject_splice_motif(seq|genome, exons, strand)         (D6)
                       (named to avoid collision with top-level rigel/splice.py)
  bam.py               coordinate transforms, CIGAR, segment builders         (keep)
  read_name.py    NEW  Origin, parse_origin  (extracted from truth.py)        (T1)

  # ── config (data only, serializable) ──────────────────────────────
  config.py       NEW  FragmentParams, ReadParams, GdnaParams, NrnaParams,    (D4)
                       AbundanceParams — unified sim params (one definition,
                       both engines consume). Capture config stays in capture/.

  # ── shared engine internals ───────────────────────────────────────
  sampling.py     NEW  truncated_normal_frag_lengths(...);                    (D2,D3)
                       build_gdna_strand_regions(...) (Beta partition) —
                       the algorithms both engines duplicated

  # ── read-generation engines (distinct execution models) ───────────
  reads.py             ReadSimulator (in-memory generator) — uses sampling/   (slim)
                       config; gDNA-region + frag-len dup removed
  wgs_engine.py   NEW  WholeGenomeSimulator (vectorized/parallel/pysam),      (split
                       extracted from whole_genome.py; uses sampling/config    out)

  # ── genome+annotation building ────────────────────────────────────
  annotation.py        GeneBuilder (OOP incremental) — uses splice_motif      (slim)
  synthetic_genome.py  batch functional gen — uses splice_motif, intervals    (slim)
                       (genome-path unification deferred — see §6)

  # ── capture (subpackage: the clean 3-way split) ───────────────────
  capture/
    __init__.py        re-export CaptureConfig, CaptureScenario, CaptureSampler
    config.py     NEW  CaptureConfig, CaptureScenario (data)                   (T2)
    sampler.py    NEW  CaptureSampler + probe loading (runtime)               (T2)
    design.py     NEW  probe DESIGN (moved out of suite.py:375–537)           (T2)

  # ── truth (parsing split out to read_name.py) ─────────────────────
  truth.py             write_post_capture_truth, summarize, counting          (slim)
  manifest.py          condition_dir_name, write/load_manifest                (keep — model)

  # ── orchestration ─────────────────────────────────────────────────
  orchestrator.py NEW  the single condition-grid loop (D1): grid build +      (D1)
                       per-condition simulate + manifest/truth emission
  scenario.py          Scenario/ScenarioResult (single-condition, in-memory)  (keep)
  suite.py             THIN: synth genome + probe design + build config +     (slim)
                       call orchestrator
  whole_genome.py      THIN: parse YAML config + call orchestrator + CLI main (slim)
  locus_sweep.py       YAML parametric locus sweep (keep)

  # ── evaluation ────────────────────────────────────────────────────
  benchmark.py         run_benchmark, BenchmarkResult, TranscriptAccuracy     (keep)
  analysis.py          EXCISE stale analyze_calibration/main; keep the live   (S1)
                       fragment-confusion/parse helpers bench_calibration uses
                       (or fold those into bench_calibration's module)
```

## 4. Phased migration (each phase shippable + verified independently)

Ordered safest-first; every phase keeps `pytest tests/` green and (from Phase 4) passes the
suite-output diff gate.

- **P0 — safety net.** Add a `tests/sim/test_suite_output_golden.py` (or a script) that runs a
  tiny pinned suite and snapshots manifest + truth + per-condition BAM record counts. This is the
  diff gate for D1; build it *before* touching orchestration.
- **P1 — primitives extraction (D5, D6).** Create `intervals.py` + `splice_motif.py`; point
  `capture.py`, `suite.py`, `synthetic_genome.py`, `annotation.py` at them; delete the dup copies.
  Pure move; verified by existing tests. **Lowest risk, immediate readability win.**
- **P2 — read_name split (T1).** Move `Origin`/`parse_origin` to `read_name.py`; `truth.py` and
  `analysis.py`/`benchmark.py` import from there. `__init__` re-exports unchanged.
- **P3 — shared sampling + unified config (D2, D3, D4).** Create `sampling.py` + `config.py`; make
  both engines consume them. The config merge is the fiddly bit (field-name reconciliation:
  `sim_seed`↔`seed`, sweep-axis fields); keep back-compat shims if any caller constructs the old
  dataclasses directly (tests do — `test_oracle_bam.py`, `test_sim_gdna_overdispersion.py`).
- **P4 — orchestrator extraction + driver unification (D1). DONE.** `orchestrator.run_condition_grid`
  is the single condition loop; `suite.main` and `whole_genome.run_simulation` both build the grid
  and call it. Per-condition seeding (Option A, user-approved) is now used by both — preserving
  `suite.main` byte-for-byte (snapshot-verified) and **fixing `run_simulation`'s former
  same-seed-for-all-conditions latent bug** (now distinct per-condition seeds; the manifest gains a
  unified schema incl. `seed`). The seed helpers moved to `orchestrator.py` (suite re-exports for
  tests).
- **P5 — capture subpackage (T2).** Split `capture.py` → `capture/{config,sampler}.py`; move probe
  design `suite.py:375–537` → `capture/design.py`. `capture/__init__.py` re-exports; `test_sim_capture.py`
  imports may need path updates (acceptable — internal test).
- **P6 — engine module split.** Move `WholeGenomeSimulator` out of the 2343-line `whole_genome.py`
  into `wgs_engine.py`; `whole_genome.py` becomes the thin config-parse + CLI frontend.
- **P7 — analysis staleness (S1).** Delete dead `analyze_calibration`/`main`/`evaluate_suite.py`
  path (confirm nothing live calls them); relocate the live confusion/parse helpers next to
  `bench_calibration.py`. (Verify against the calibration-benchmark skill's documented routing.)

## 5. Risk & verification

- **Highest risk:** P4 (orchestration identity) and P3 (config merge — tests construct the configs
  directly). The P0 snapshot + the existing `test_sim_capture.py`/`test_whole_genome_sim_config.py`/
  `test_oracle_bam.py`/`test_sim_gdna_overdispersion.py` are the guardrails.
- **Public API:** `rigel/sim/__init__.py` re-exports (`Scenario`, `run_benchmark`, `GDNAConfig`,
  `CaptureConfig`, …) stay stable; external imports (`scripts/sim/*`, `tests/`) are unaffected
  except internal test imports for P5/P6 (flagged).
- **Each phase is its own commit**, independently green, so we can stop/ship at any phase boundary.

## 6. Explicitly deferred (with rationale)

- **Merging the two read engines** — deliberately NOT done. Different execution models (string
  generator vs numpy/pysam/parallel); the audit found merging would lose vectorization and add a
  genome-IO abstraction for no behavioural gain. We share their internals (P3) and stop there.
- **Unifying the two genome-building paths** (`genome`+`annotation` vs `synthetic_genome`) — leave
  both; only the duplicated splice-motif injection is extracted (P1). Full unification (one gene
  model) is a larger semantic change with test risk and little duplication payoff beyond P1.

## 7. Open questions for review

1. **Scope for now:** execute all of P0–P7, or land P0–P4 (the duplication that bit you — sampling,
   config, orchestration) and treat P5–P7 (capture subpackage, engine split, analysis cleanup) as a
   follow-up? (Recommend: P0–P4 first, as one reviewable arc; then P5–P7.)
2. **`capture/` subpackage vs flat** `capture_config.py`/`capture_sampler.py`/`capture_design.py` —
   subpackage reads cleaner; flat is less import churn. (Lean subpackage.)
3. **`config.py` naming** — there's already a top-level `rigel/config.py`. The sim one is
   `rigel/sim/config.py` (distinct path); acceptable, or name it `sim_params.py` to avoid any
   mental collision? (Lean `sim/config.py` — the package prefix disambiguates.)
