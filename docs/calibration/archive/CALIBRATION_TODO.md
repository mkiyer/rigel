# CALIBRATION_TODO — open issues toward a production-ready deconvolution

**Status:** living tracker. Opened 2026-06-09; rewritten 2026-06-10 to the **decoupled + count-channel**
architecture. The phased roadmap is `count_channel_capture_design.md` (authoritative); this file tracks
the live open fronts and their state. Prior issue numbering (#1–#5) and the count-overdispersion /
strand-clean-O2 prototypes are **superseded** by the decoupling (`decoupled_calibration_design.md`,
main@8d8b896) and the count-channel redesign — see git history / `archive/` if resurrecting.

## Architecture (shipped)

- **Decoupled strand/count** (main@8d8b896): each node routes to the strand module (Beta-Binomial
  posterior, unbiased) or the count module (raw density ratio, biased under capture) — never a product.
  Joint archived in `archive/joint_deconvolution.md`.
- **gDNA + RNA strand Beta-Binomial overdispersions** fitted symmetrically (`gdna_strand.py`), so
  unstranded data is uninformative.

## Residual the count work must close

The post-decoupling leak lives in **unstranded + capture** (~18–21% gDNA→RNA) and is a count-module
*mean* bias: under capture, exon gDNA imputed from depleted off-target neighbours under-calls ~2×
(0.41 vs 0.91 truth; `diag_imputation_truth.py`). Stranded+capture is fine (flagship ~3.7–4.6%, strand
rescues). Phase 4-mean (splice-junction gDNA-fraction) attacks this directly.

## Open fronts (map to `count_channel_capture_design.md`)

| front | state | next action |
|---|---|---|
| **Phase 1** precision-weighted blend `w=(2κ−1)²` | implemented, **uncommitted** | commit once 4-mean debiases count and the flagship +1pp clears |
| **Phase 2** nascent 8-scenario benchmark | sim committed (branch `phase2-nascent-benchmark`) | **generate + run the suite** — the validation harness for all nascent-aware work |
| **Phase 4-mean.1** eligible 2-term splice fraction | **implemented** (`splice_junction.py`, Tests A/B) | validating on the 20-cond benchmark (net-flow rerun in progress) |
| **Phase 4-mean.2** run-interior sweep (carry-over) | designed (`count_mean_bias_design.md` §5) | implement after .1 validated |
| **Phase 4-mean.3** 3-term + strand-resolved per-strand sweep | designed (§5, §5.1) | needs nascent suite; absorbs the retired Phase 3 (clean→count) |
| **Phase 4-var** count posterior variance (`var∝mean²`) = remaining_phases **Phase 1** | **SHIPPED `7ccc3b0`** | Poisson floor everywhere (regions LOESS + boundary sides) + LOESS variance~mean; feeds the quantile |
| **Phase 5** FP-rate quantile knob = remaining_phases **Phase 2** | **SHIPPED `7ccc3b0`** | `gdna_deconv_quantile` (default 0.5); retired `gdna_strand_llr_bias`; validated on the 3-pool suite (`phase2_design.md` §5) |

## Deferred / future

- **FL-consistency of the splice fraction** (`f_b·M_region` vs the region eff-length ratio) — future
  investigation, matters for short exons. Memory `calibration_gdna_fraction_fl_consistency`,
  `count_mean_bias_design.md` §9.1.
- **RNA strand BB junction double-count** — minor; one side per boundary per fragment. Memory
  `rna_strand_bb_double_count_followup`.
- **DNA-fraction external lever** (old Issue #4) — unscoped; revisit only if a residual remains after
  Phase 4/5.

## Validation harnesses

- 20-condition gDNA suite (`~/Downloads/rigel_runs/gdna_benchmark_5mb`), skill `calibration-benchmark`.
  **Note:** `--skip-frag-analysis` also skips the **net-flow** metric — rerun analysis *without* it (and
  with `--skip-quant` if quant is fresh) to get `net_flow_per_condition.tsv`.
- 8-condition nascent suite (`scripts/sim/configs/nascent_benchmark_5mb.yaml`) — pending generation.
