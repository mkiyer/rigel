---
name: calibration-benchmark
description: Run and evaluate Rigel's gDNA-vs-RNA deconvolution on the synthetic benchmark suite via the net fragment-flow analysis (transcript- and locus-level net leak, with pool + per-fragment confusion and root-cause hooks). Use when validating calibration/quant changes, before a release, or when the user asks to "benchmark calibration", "evaluate gDNA/RNA classification", or "run the synthetic suite".
---

# gDNA-vs-RNA deconvolution benchmark (synthetic)

Validates how well calibration + the per-locus EM separate **gDNA from RNA** against simulated
ground truth, across `gDNA level × strand-specificity × capture on/off`. The canonical workflow
is `rigel.sim.analysis.main` (run via `scripts/sim/evaluate_suite.py`); the **primary** metric is
the **net fragment-flow** deconvolution.

## Prerequisites
- Run everything inside the activated `rigel` conda env:
  `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`
- Release suite: `/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb` (20 conditions:
  gDNA ∈ {none, gdna25, gdna100, gdna400, gdna1000} × ss ∈ {0.99, 0.50} × capture ∈ {off, on},
  zero nRNA). Config: `scripts/sim/configs/gdna_benchmark_5mb.yaml`. Each condition dir has
  `sim_oracle.bam`, `truth_*.tsv/json`; truth is in `manifest.json` (`n_*_observed`) and each
  `truth_summary.json` (`origin_counts`).

## Run + evaluate (one command)
```bash
SUITE=/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb
python scripts/sim/evaluate_suite.py --sim-base "$SUITE"            # build index + quant + analyze
# subset / re-evaluate without re-quant:
python scripts/sim/evaluate_suite.py --sim-base "$SUITE" --conditions <name> --skip-quant
```
- Builds the index (`rigel index`) and runs `rigel quant` on each condition's **`sim_oracle.bam`**
  with `--annotated-bam --emit-locus-stats` (oracle alignment isolates calibration/quant from
  aligner noise). Quant is skipped per-condition when `rigel_out/quant.feather` exists; delete it
  to force a re-quant after a scoring/calibration change.
- `scripts/sim/bench_calibration.py` is a deprecated shim forwarding here (`--run`/`--force` kept).
- Writes into the suite: `analysis_report.txt` (all sections) plus four focused TSVs —
  **`net_flow_per_condition.tsv`** (pool rollup), **`net_flow_per_locus.tsv`**,
  **`net_flow_per_transcript.tsv`**, and `abundance_per_condition.tsv`.

## Primary metric — net fragment flow (why it's the right target)
Hard per-fragment label recovery is the *wrong* target: an unspliced RNA fragment and a gDNA
fragment from the same locus can be **sequence-identical and unrecoverable** — and that's fine.
What matters is the **net** deconvolution error per component. The analysis builds a per-locus
flow matrix `flow[true_origin][assigned]` (true origin from the oracle read name via
`parse_origin`; assigned from the hard MAP label in annotated-BAM `ZT`/`ZF`), then reduces to
**net** flow `net(a→b)=flow[a][b]-flow[b][a]`. Symmetric (unrecoverable) misassignment cancels;
only systematic bias survives. Components per locus = its transcripts + one `gdna@L` (home locus
is each component's modal `ZL`); intergenic gDNA is `gdna@-1`.

Outputs answer the research questions directly:
- **Pool net leak & direction** — signed `net_gdna_to_rna` per condition (+ = gDNA→RNA leak,
  − = RNA→gDNA siphon), split in-locus vs intergenic. Balanced ⇔ ~0.
- **Per-transcript Δ = observed − expected**, decomposed into `net_from_gdna` (contamination)
  vs `net_from_rna_isoforms` (zero-sum RNA reallocation) vs `cross_locus` (small residual).
- **Root cause** — Spearman of `net_from_gdna` vs covariates (depth, n_exons, single-exon,
  spliced_length); **identifiability** diagnostic comparing single-exon (gDNA-identical) vs
  multi-exon transcripts. Per-transcript/per-locus TSVs carry the covariate-joined detail.

## Secondary metrics (still reported)
- **Pool fraction** — library-wide gDNA fraction from `summary.json` `quantification`.
- **Gross per-fragment confusion** (hard `ZF` bit `0x04`) — `gDNA→RNA leak` / `RNA→gDNA siphon`;
  useful but **inflated by unrecoverable identical-sequence cases**, which the net metric removes.
- Transcript abundance accuracy (Spearman/MARD/false positives) and calibration scalars.

## Interpreting / root-causing
- **gdna_none** conditions are pure FP tests (any net gDNA→RNA is false). SS=0.50 (unstranded)
  is the hard case (no strand clue). Net leak should grow with gDNA level and fall with stranding;
  capture-on concentrates gDNA on exons (co-located with mRNA) → expect more in-locus leak.
- Localize with `net_flow_per_locus.tsv` (join to `rigel_out/loci.feather`:
  `gdna_prior_count`/`rna_prior_count`/`gdna_eff_len_em`) to separate **calibration prior** error
  (per-locus prior wrong) from **EM/assignment** error (prior ~right, likelihood mis-labels).
- Single-exon transcripts have no isoforms, so their Δ ≈ `net_from_gdna`; a nonzero class-mean
  signals systematic gDNA→RNA leak concentrated where RNA is sequence-identical to gDNA.

## Regenerating the suite (only if needed)
Pre-generated. Rebuild: `python scripts/sim/simulate_suite.py --config scripts/sim/configs/gdna_benchmark_5mb.yaml`.
The synthetic genome/transcriptome + capture probes are under the suite's `reference/`.
