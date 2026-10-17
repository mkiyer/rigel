# Bare constants: what the last cleanup stage needs from the owner

The constants stage moves every bare number into `config.CONSTANTS` (or `native/constants.h`), each documented
with WHY it has its value (CLAUDE.md: no magic numbers). A read-only inventory (2026-10-16, at `d5a08e46`)
traced each one through its comments, the docs and `git log`. It changes no number: every value stays, and
only its home and its stated reason move.

`tests/test_constants.py` sees module-level UPPER_CASE Python assignments only: inline literals, function
defaults and every C++ number are outside it. That is how these survived.

## Settled without you (done in the stage)

| constant | where | the reason, as found |
|---|---|---|
| `SQUAREM_BUDGET_DIVISOR = 3` | `em_solver.cpp` | each SQUAREM iteration runs three VBEM steps (two plain, one stabilising), so the budget stays at `EMConfig.iterations` |
| `_block_rows(grid.size, 8)` | `landscape.py` | the argument is an itemsize: 8 is `float64` |
| the resolver's `-999` | `resolve_context.h` | any value outside the strand codes does; it becomes a named sentinel |
| `[:3]` examples in two index errors | `index.py` | `CONSTANTS.qc.report_examples` (5) already exists for the same job |
| `min_count = 2` defaults | `splice_blacklist.py` | duplicates of `IndexConfig.splice_blacklist_min_count`, whose reason is recorded |
| `calib_refit_iters = 3` | `config.py` | recorded only in `daa32a13`: refit 3 captured 96 % of the available gain on the July suite |

## Ruled by the owner (2026-10-17)

1. `max_spliced_reads = 1000`: repaired before 0.8.0, not named: the junction-strand tag is read per record, in wave
   1's step 2 (`ISSUES: sj-strand-tag-chosen-from-the-first-reads`, now 0.8.0).
2. K = 3 is a valid tolerance for implicit splicing; it is renamed to say what it is (`intron_overhang_tolerance`,
   in the renames stage), and wave 1's step 3 judges its other uses.
3. The report's tiles: agreed (the strandedness tile reads `strand_live`; the other cuts are display constants).
4. The junction-depth keys: deleted with the schema-4 bump.
5. The summary's rounding: named beside `SUMMARY_DECIMALS`.
6. The resource numbers: agreed, with a gate for the landscape's tiling.

## The questions as put

Each line says what the number controls, what the record holds, and a recommendation.

1. **`max_spliced_reads = 1000`** (the tag detector; `bam_scanner.cpp`). No reason recorded since 2026-02. The
   0.8.1 per-record repair (`ISSUES: sj-strand-tag-chosen-from-the-first-reads`) deletes it.
   *Recommended:* name it in `constants.h` as "a sample, deleted by the per-record repair".
2. **The anchor tolerance K = 3** (`BamScanConfig.splicing_anchor_tolerance`). Its one recorded reason ("mapping
   noise ≤ 2 bp; 3 is a margin") was for a boundary-clearance use that no longer exists. Today K lets an
   annotated intron overhang a mate gap, contradicts a transcript only beyond K bp, and cuts up to K aligned
   bases inside an intron. Wave 1 step 3 audits each use (`ISSUES:
   the-anchor-tolerance-applies-beyond-implicit-splicing`). The resolver's own default is 0, which contradicts
   the config. *Recommended:* restate the reason for the implicit-splice use, and let step 3 decide the others.
3. **The report's tiles** (`report/model.py`): mapping 0.75/0.5, strandedness 0.75/0.6, gDNA 0.1/0.3, genic 0.6,
   the "> 50 % gDNA" KPI, a 20,000-gene embed cap and the rounding digits. No reason recorded for any of them.
   `summary.json` already carries the model's own strand verdict (`strand_live`), and a precedent replaced the
   report's capture thresholds with calibration's own answer. *Recommended:* the strandedness tile reads
   `strand_live`; the other cuts move to `CONSTANTS.qc` as display choices, with your sign-off as the reason.
4. **The summary's junction-depth keys** `n_sj_depth_ge_100` / `_ge_1000` (`strand_model.py`). They were added
   for the RNA strand-overdispersion fit, which is deleted. *Recommended:* delete them with the schema-4 bump in
   the renames stage.
5. **The summary's rounding:** `round(total, 2)` and `round(posterior_variance, 8)` (`cli.py`). No reason recorded;
   the variance is of order p(1−p)/n, about 1e-8 for a million training fragments, so 6 decimals would print 0. *Recommended:* name them beside `SUMMARY_DECIMALS`,
   the variance in significant digits.
6. **Resource numbers that move no result:** the spill writer's queue of 2 chunks (memory can exceed
   `buffer_size_bytes` by that many), capture efficiency's tile `16` (changed from 8 without comment), and the
   simulator's gzip level 3, pgzip threads 4 (1 per shard), the `I` base quality and the oracle BAM's @PG version
   "1.0". *Recommended:* into `CONSTANTS.resources` as resource choices; the @PG version stamps
   `rigel.__version__`; the quality becomes one named constant.

Not a constant, but found on the way: the landscape's tile height changes the order of its float sums
(`sums[b] += …`), so `config.py`'s claim that "the block size cannot reach the arithmetic" is false for it, and
no gate holds the tiling to the bit. *Recommended:* a gate that two tile heights give bit-identical sums, and,
if they do not, a fixed summation order.
