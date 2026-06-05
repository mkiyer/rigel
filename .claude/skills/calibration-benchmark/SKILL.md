---
name: calibration-benchmark
description: Run and evaluate Rigel's calibration accuracy (gDNA-vs-RNA classification) on the synthetic hybrid-capture benchmark suite, with pool-level + per-fragment confusion metrics and root-cause hooks. Use when validating calibration/quant changes, before a release, or when the user asks to "benchmark calibration", "evaluate gDNA/RNA classification", or "run the synthetic capture suite".
---

# Calibration accuracy benchmark (synthetic hybrid-capture)

Validates how well calibration + the per-locus EM separate **gDNA from RNA** against
simulated ground truth, across `gDNA level × strand-specificity × capture on/off`.

## Prerequisites
- Run everything inside the activated `rigel` conda env:
  `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`
- The suite lives at `/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb`
  (a `manifest.json` + per-condition dirs with `sim_oracle.bam`, `truth_*.tsv/json`, `reference/`).
  Each condition dir is `gdna_{none,high}_ss_{0.99,0.50}_nrna_none_capture_{off,on}`.
  Truth is in `manifest.json` (`n_*_observed`) and each `truth_summary.json` (`origin_counts`).

## Run + evaluate (one command)
```bash
SUITE=/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb
python scripts/sim/bench_calibration.py --sim-base "$SUITE" --run --force
```
- `--run` builds the index (`rigel index`, fresh) and runs `rigel quant` on each
  condition's **`sim_oracle.bam`** (oracle alignment — isolates calibration/quant from
  aligner noise). `--force` re-quants even if `rigel_out/` exists (use it after any
  calibration/scoring change). Omit `--run` to evaluate existing `rigel_out/`.
- Writes `bench_calibration_metrics.tsv` + `bench_calibration_report.txt` into the suite.
- NOTE: `scripts/sim/evaluate_suite.py` (the older harness, `rigel.sim.analysis.main`) has
  a stale `analyze_calibration()` that assumes the pre-rebuild `summary.json` schema and
  crashes (`calibration` is now `None`). `bench_calibration.py` reads the current schema and
  reuses that harness's working `build_index`/`run_quant`/`parse_origin`/fragment-decode.

## Two metrics — and why both matter (they can disagree)
1. **POOL (fractional)** — library-wide gDNA fraction from `summary.json` `quantification`
   (`gdna_total + intergenic_total` vs `mrna_total + nrna_total`). This is what calibration
   is *for*. Truth: `gdna_none` → 0%, `gdna_high` → 50%.
2. **PER-FRAGMENT (hard)** — each fragment's hard gDNA/RNA label (annotated-BAM `ZF` bit
   `0x04`) vs true origin (oracle read name via `parse_origin`). Reports:
   - **gDNA→RNA leak** = gDNA called RNA — the **FP-deleterious** direction the user cares
     about most (gDNA inflating RNA → false positives). Split into **in-locus** (on a gene,
     `ZL≥0`) vs **intergenic** (`ZL<0`).
   - **RNA→gDNA siphon** = RNA called gDNA — the FP-safe direction.
   - **FL leak vs caught** — mean template length of leaked vs caught gDNA (gDNA truth ~150,
     RNA ~250); leaked-gDNA FL between them signals FL isn't separating them.

The pool fraction can look perfect while per-fragment leaks badly, because leak and siphon
**cancel** at the pool level. Always report both.

## Interpreting / root-causing
- **gDNA_none** conditions are pure FP tests (any called gDNA is false). SS=0.50 (unstranded)
  is the hard case (no strand clue to clean ρ₀ / the global density).
- If leak concentrates **in-locus** and grows with capture, the cause is capture enriching
  gDNA onto exons where it is co-located with mRNA (and, at low SS, strand-indistinguishable);
  the per-fragment likelihood (location + FL + strand) can't separate them. Localise with a
  per-locus join of annotated-BAM truth (`parse_origin` + `ZL`/`ZF` tags) against
  `rigel_out/loci.tsv` (`mrna`/`nrna`/`gdna`/`gdna_prior_count`).
- Distinguish **calibration prior** error (per-locus `gdna_prior_count` wrong) from
  **EM/assignment** error (prior ~right, fragments mis-labeled by likelihood / hard sampling).

## Baseline (current `main`, June 2026) — flag regressions against this
- Pool-level gDNA fraction: within **±2.2%** of truth everywhere; FP at SS=0.99 ≈ 0;
  FP at SS=0.50 (unstranded, gdna_none) ≈ 1.6–1.7%.
- Per-fragment gDNA→RNA leak: capture-off 2.9% (ss99) / 5.9% (ss50); **capture-on 10% (ss99)
  / 20% (ss50)** — all in-locus. This capture leak is the open production-blocker
  (see the findings memory). Intergenic gDNA is always caught.

## Regenerating the suite (only if needed)
The suite is pre-generated. To rebuild it: `scripts/sim/simulate_suite.py` (see
`rigel.sim.suite` / `rigel.sim.whole_genome`). The synthetic genome/transcriptome +
capture probes are under the suite's `reference/`.
