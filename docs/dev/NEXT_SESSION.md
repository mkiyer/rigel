# NEXT SESSION — start here (2026-09-13, after W10; W11 measured, W12 derived; W13 and the port remain)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C, and §8 is the design of the
remaining steps. The θ quadrature's derivation is `docs/dev/THETA_QUADRATURE.md`. This file is only how to
begin, and it carries the decisions the owner still owes from W11.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, the commit after `a5f68a99`). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md`,
> `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` and `docs/dev/THETA_QUADRATURE.md`. The PRE-PORT WORKLIST stands
> at: W1–W10 DONE and committed; W11 MEASURED (the coverage census; my decisions on its tiers are below and
> are the first thing to act on); W12 DERIVED, prototype A A/B'd (the result is in the plan's W12 row);
> W13 (the re-capture of the deep library's sweeps and fresh identity references for the port thread) not
> started. Run `scripts/design/preflight.py` and the suite before touching anything; `CLAUDE.md`'s baseline
> line is the count to reproduce.
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel, the suite, and
> timing on back-to-back profiler pairs; a restructure is proven bit-identical against fresh references
> (`hygiene_identity_*` are current as of `bdc1089c`); a removal on the default rows of both substrates
> (`arms/2026-09-13_preport/oracle_{ladder,test}_default.json`); one mechanism at a time; no magic numbers;
> each step's plan presented and gated; each item its own commit, committed on my go; never `ruff format
> scripts/`; patch `calibrate` through `importlib.import_module`; any config value is a `--set` arm. Before
> migrating or extending an instrument, census it and retire it if its question is closed (2026-09-13).

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## W11 — the decisions owed (the census is in the plan's W11 row; the data is not in the tree)

Every removal is its own commit, byte-identical on the default rows of both substrates; a retired
`scripts/` file moves collected by −4 (design/) or −2 (a test), a `src/` module by −3.

**Tier 1 — dead by construction, trivial, no decision needed beyond "go":** `simplex_logodds._sp` (a local
helper with no call); `region_chain.is_region` (a property with no reader); `substrate.n_objects` and
`splice_graph.n_transcripts` (properties whose only "callers" are same-named members of other classes);
`calibration/strand_summary.py`'s `strand_specificity` and `read1_sense` (duplicates of `StrandModel`'s,
never read); `frag_length_model.observe_batch` (no caller).

**Tier 2 — design sockets nothing feeds; the owner decides whether they stay:**
* `CompositionPriors.rna` and `_rna_arm`'s fitted branch — never fed; the docstring calls it "the repair's
  landing point rather than speculative surface" for the `f_g → 1` vertex asymmetry.
* `boundary_rna_reach` — a keyword parameter of `calibrate` and `build_region_geometry` that nothing sets:
  the RNA reach taper at contiguous boundaries, kept as "one argument so an A/B varies one thing"; the
  unbounded form ships. (`crossing_eff_length`'s reaches stay: the sj route rates use them.)
* the pre-v3 legacy index path — `scan_and_buffer` returns no payload when the index has no
  `regions.feather`, and three `calibration is None` branches downstream carry a no-calibration mode.
* `DrainQC.from_dict` — never reached because the cache refuses a drained payload.

**Tier 3 — dead or entry-point-less simulator features (the suite exercises none; panel builds exercise
the live simulator paths but not these):**
* `src/rigel/sim/locus_sweep.py` (1,020 lines, 0 % executed) + `scripts/sim/locus_sweep.py` (a 21-line
  wrapper); no manual names it; the toy harness covers the question.
* `src/rigel/sim/net_flow.py` (678 lines, 19.5 %) — a CLAUDE.md row and a TESTING.md line, but no script or
  CLI calls it: only `tests/test_net_flow.py` (100 lines). Retire, or give it an entry point.
* `src/rigel/sim/suite.py` (996 lines; `main` and its helpers never executed) with `scripts/sim/simulate_suite.py`
  and `scripts/sim/snapshot_suite.py` (154), and `scripts/sim/build_toy_2exon_reference.py` (113, "for the gDNA
  effective-length study", which is closed) — the synthetic-mini-genome path `panel.py` superseded. No manual
  names any of them. `sim/synthetic_genome.py` (553) is shared with `tests/test_sim.py`.

**Tier 4 — coverage gaps, not dead code (no removal):** the five CLI command bodies; the silent policy through
`calibrate` (only the instruments run it); the simulator's sharded writers; the zarr splice-blacklist loader
and `index.build`'s blacklist branch (a real feature with no test).

## Where the worklist stands

W1–W10 DONE and committed (W10: `0f3ca3e3`…`bdc1089c`, nine commits; `84fa14b9` the plan's tick). W11 MEASURED
(`a5f68a99`, decisions above). W12 DERIVED (`docs/dev/THETA_QUADRATURE.md`, `a5f68a99`); prototype A (the
trapezoid endpoint weights, patched outside `src/`) A/B'd on both panels' `g00` rows at `sweep_n_tilt` 60 and 30
and REFUTED (≤ 0.7 fragments anywhere; the K_t 30 failure untouched): the resolution term is the whole
mechanism. Next after the owner's W11 decisions: the W12 design decision (prototype B, Gauss–Hermite per λ, or
the analytic marginal — both change the cube's shape and must handle the delivered cube rows), its prototype
and A/B, then `src/`; W13; then the port (plan §F).

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee.
* The two lattice knobs are the model of a kept knob: dimensionless, a guarantee, a recorded ladder, not on
  the CLI. The EM's `gdna_em_llr_bias` stays (the owner's).
* `drain`, `row` and `face` stay (2026-09-13). The arcsine coordinate stays refused; the vertex atom is parked.
* The deconvolved arrays are `count_<population>_<axis>` (2026-09-13): counts on their axis, never masses.
* Four instruments retired rather than migrated (2026-09-13, the owner's ruling after a census); the record and
  the reasons are in `ISSUES: hygiene-ledger`. Second-tier candidates the owner has not ruled on:
  `native_parity_on_real_data`, `prior_units_check`, `accumulator_cost`, `verify_index_rebuild`.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
