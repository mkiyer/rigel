# Rigel Workspace Context — Session Transfer

_Exported: 2026-04-18_

## What we accomplished this session

### 1. Mappability-aware calibration: analytical stress test (DONE)

Built and ran a 24-scenario synthetic stress test (`scripts/benchmark/mappability_stress/stress_calibration.py`) that exercises the density calibration pathway with and without the mappable-region mask across a grid of λ_G, strand specificity, and aligner retention values. Uses real human genome region tables from a pre-built rigel index.

**Key findings** (published in `docs/mappability/stress_test_findings.md`):
- Mask recovers λ_G when aligner drops unmappable reads (retention ≤ 0.1): |rel_err| 91% → 48%
- Density pathway blind below ~1e-3 frags/bp (Poisson identifiability wall)
- Length-weighted P₁₀ has structural -20% to -30% bias
- Mask loses 49% of regions due to strict containment filter

### 2. End-to-end pipeline validation on STAR-aligned BAM (DONE)

Simulated 5M RNA + 5M gDNA fragments (unstranded, SS=0.5), aligned with STAR 2.7.11b, ran `rigel quant --em-mode map` with two indexes (mappable vs no_mappable).

**Critical finding:** The mappability mask **collapsed the density pathway to λ_G=0** on real STAR data because contained regions are predominantly exonic (RNA-dominated). Without the mask, density pathway correctly estimated λ_G=9.4e-4 (truth ~1e-3). However, the downstream EM was remarkably robust — mRNA/gDNA totals matched within <1% regardless of calibration.

Results published as "Pipeline validation" section in `docs/mappability/stress_test_findings.md`.

### 3. Unified joint-likelihood calibration proposal (DONE)

Developed and published a comprehensive proposal (`docs/mappability/unified_proposal.md`) to replace the current two-pathway blend architecture with a single joint log-likelihood. Three key innovations:

1. **Replace P₁₀ with Poisson mixture MLE** — zeros become informative; no arbitrary cutoff; continuous mappable fraction $f_i$ replaces hard mask
2. **Unmapped-read pathway (uT:3)** — STAR's "too many multimapping loci" unmapped reads are overwhelmingly gDNA; the unmappable genome span (~300 Mb) provides enormous Fisher information; single pathway alone predicted to achieve sub-percent λ_G precision
3. **Joint likelihood** — all three pathways (density, strand, unmapped) as terms in the same objective, self-weighting by Fisher information

5-phase implementation plan in the doc. Phase 1 (unmapped pathway probe) is the recommended next step.

---

## Key file locations

### On-disk artifacts (scratch)

```
/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/mappability_stress/
├── idx_mappable/          # rigel index with mappable.feather (693k intervals)
├── idx_no_mappable/       # identical index WITHOUT mappable.feather
├── sim/
│   ├── gdna_high_ss_0.50_nrna_none/
│   │   ├── sim_R1.fq.gz   # 10M pairs (5M RNA + 5M gDNA)
│   │   └── sim_R2.fq.gz
│   └── truth_abundances_nrna_none.tsv
├── runs/gdna_high_ss_0.50_nrna_none/
│   ├── star/
│   │   ├── Aligned.out.bam       # ~1 GB
│   │   ├── name_sorted.bam       # ~800 MB (rigel input)
│   │   └── Log.final.out
│   ├── rigel_mappable/
│   │   ├── summary.json          # λ_G=0 (BROKEN — mask collapsed density)
│   │   ├── quant.feather
│   │   └── ...
│   └── rigel_no_mappable/
│       ├── summary.json          # λ_G=9.42e-4 (correct)
│       ├── quant.feather
│       └── ...
└── results/
    ├── results.tsv               # 24-scenario stress test results
    └── metadata.json
```

### Reference data

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/
├── genome_controls.fasta.bgz     # GRCh38 + controls (286 refs)
├── genes_controls.gtf.gz         # GENCODE v46 + controls
├── star_index/                   # STAR 2.7.11b index
├── rigel_index/                  # production rigel index
└── salmon_index/

/scratch/mkiyer_root/mkiyer0/shared_data/alignable/
├── alignable_mappable_rl50_grch38_star.bed      # mappable regions (rl=50)
├── alignable_unmappable_rl125_grch38_star.bed   # unmappable regions (rl=125)
└── ...
```

### Repo files created/modified this session

```
docs/mappability/stress_test_findings.md          # analytical + pipeline findings
docs/mappability/unified_proposal.md              # the joint-likelihood proposal
scripts/benchmark/mappability_stress/
├── stress_calibration.py                         # 24-scenario stress test driver
├── sim_unstranded_gdna.yaml                      # sim config for pipeline validation
└── run_pipeline.sh                               # STAR → sort → quant × 2 pipeline
```

### STAR location

```
/home/mkiyer/sw/miniforge3/envs/star/bin/STAR     # STAR 2.7.11b
# OR: module load Bioinformatics star/2.7.11a-clzmxha
```

---

## Key calibration numbers for reference

### Pipeline validation (unstranded, 50% gDNA)

| Metric | idx_mappable | idx_no_mappable | Truth |
|---|---|---|---|
| Calibration λ_G | **0.00e+00** | 9.42e-4 | ~1e-3 |
| Calibration gDNA_frac | 0.0% | 27.1% | ~46% |
| EM mRNA total | 4.990M | 4.987M | 5.00M |
| EM gDNA total | 4.599M | 4.632M | 4.60M |
| EM nRNA (FP) | 8.4e4 | 5.5e4 | 0 |

### Stress test summary (where mask helps)

- retention=0.0, λ≥1e-3: |rel_err| no_mask=0.910, mask=0.477
- retention=1.0, λ≥1e-3: |rel_err| no_mask=0.146, mask=0.485
- Mask helps when aligner drops unmappable reads; hurts when aligner retains most reads

---

## Recommended next steps (in priority order)

### Step 1: Unmapped pathway probe (highest leverage)

Write `scripts/debug/unmapped_pathway_probe.py`:
- Read the existing name-sorted BAM via pysam
- Count fragments by `uT` tag value (0–4)
- Compute $S_\text{rep}$ from `alignable_unmappable_rl125_grch38_star.bed`
- Estimate $\hat\lambda_G = (1-\pi_\text{RNA}) \cdot N_\text{uT3} / S_\text{rep}$
- Compare to truth (~1e-3) and existing pathway estimates

This is the single highest-value experiment. If uT:3 gives a good λ_G → proceed to implement it in the production code.

### Step 2: Poisson mixture MLE to replace P₁₀

Replace `_density_pathway()` internals in `calibration.py` with a two-component Poisson mixture EM. Use continuous mappable fraction $f_i$ instead of hard containment mask. Zero-count moment estimator as diagnostic bound.

### Step 3: Joint likelihood unification

Combine density + strand + unmapped into single log-likelihood (§5 of unified_proposal.md).

### Step 4: Validation sweep

Re-run analytical stress test + STAR pipeline at SS={0.5, 0.9, 1.0} × gDNA={low, high}.

---

## Build/test reminders

```bash
conda activate rigel
pip install --no-build-isolation -e .   # after ANY C++ change
pytest tests/ -v                         # 1057 tests should pass
ruff check src/ tests/ && ruff format src/ tests/
```

Known pre-existing failure: `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna` — unrelated to this work.

---

## Current code architecture (calibration)

- `src/rigel/calibration.py` — main calibration logic: `calibrate_gdna()` with strand pathway + density pathway + blending
- `src/rigel/calibration.py::_density_pathway()` — **this is where P₁₀ lives** (the thing to replace)
- `src/rigel/pipeline.py` — orchestrator that calls calibration
- `src/rigel/native/bam_scanner.cpp` — C++ BAM scanner (needs uT tag counter for unmapped pathway)
- `src/rigel/index.py` — index builder (needs $f_i$ mappable fraction column + $S_\text{rep}$ scalar)
- `src/rigel/_bed_io.py` — BED loading/merging utilities
- `src/rigel/mappable.py` — `load_mappable_bed()`, `compute_containment_mask()`
