"""Diagnostic script: analyze strand model behavior with gDNA simulation."""
import logging
import tempfile
from pathlib import Path

logging.basicConfig(level=logging.WARNING)

from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from hulkrna.pipeline import run_pipeline

tmp = Path(tempfile.mkdtemp(prefix="hulkrna_diag_"))

sc = Scenario("diag", genome_length=5000, seed=42, work_dir=tmp / "diag")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
])

sim_config = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80,
    frag_max=450, read_length=100, seed=42,
)
gdna_config = GDNAConfig(
    abundance=50, frag_mean=350, frag_std=100,
    frag_min=100, frag_max=1000,
)

result = sc.build(
    n_fragments=500, sim_config=sim_config, gdna_config=gdna_config,
)
pipeline_result = run_pipeline(
    result.bam_path, result.index, sj_strand_tag="ts",
)
bench = run_benchmark(result, pipeline_result, scenario_name="diag_gdna")

print("=== BENCHMARK SUMMARY ===")
print(bench.summary())
print()

sm = pipeline_result.strand_models
print("=== STRAND MODELS ===")
print(f"exonic_spliced: p_r1_sense={sm.exonic_spliced.p_r1_sense:.4f}, n_obs={sm.exonic_spliced.n_observations}")
print(f"exonic:         p_r1_sense={sm.exonic.p_r1_sense:.4f}, n_obs={sm.exonic.n_observations}")
print(f"intronic:       p_r1_sense={sm.intronic.p_r1_sense:.4f}, n_obs={sm.intronic.n_observations}")
print(f"intergenic:     p_r1_sense={sm.intergenic.p_r1_sense:.4f}, n_obs={sm.intergenic.n_observations}")
print()

stats = pipeline_result.stats
print("=== PIPELINE STATS ===")
print(f"n_fragments:         {stats.n_fragments}")
print(f"n_intergenic:        {stats.n_intergenic}")
print(f"n_chimeric:          {stats.n_chimeric}")
print(f"n_with_exon:         {stats.n_with_exon}")
print(f"n_with_annotated_sj: {stats.n_with_annotated_sj}")
print(f"n_strand_trained:    {stats.n_strand_trained}")
print()

rna_gt = result.ground_truth_from_fastq()
gdna_gt = result.ground_truth_gdna_count()
print("=== GROUND TRUTH ===")
print(f"RNA fragments:  {rna_gt}")
print(f"gDNA fragments: {gdna_gt}")
print()

counter = pipeline_result.counter
print("=== COUNTER ===")
print(f"gdna_total:     {counter.gdna_total:.1f}")
t_counts = counter.t_counts
print(f"t_counts sum:   {t_counts.sum():.1f}")
print(f"t_counts cols:  {t_counts.sum(axis=0)}")

# Check what the exonic strand model actually sees
print()
print("=== STRAND LOG-LIKELIHOOD COMPARISON ===")
import numpy as np
from hulkrna.types import Strand

# For a sense-aligned fragment on a + strand gene:
p_exonic_sense = sm.exonic.strand_likelihood(Strand.POS, Strand.POS)
p_spliced_sense = sm.exonic_spliced.strand_likelihood(Strand.POS, Strand.POS)
p_intergenic_sense = sm.intergenic.strand_likelihood(Strand.POS, Strand.POS)

p_exonic_anti = sm.exonic.strand_likelihood(Strand.NEG, Strand.POS)
p_spliced_anti = sm.exonic_spliced.strand_likelihood(Strand.NEG, Strand.POS)
p_intergenic_anti = sm.intergenic.strand_likelihood(Strand.NEG, Strand.POS)

print(f"Sense-aligned fragment on + gene:")
print(f"  exonic_spliced:  P={p_spliced_sense:.4f}  log={np.log(p_spliced_sense):.4f}")
print(f"  exonic (diluted): P={p_exonic_sense:.4f}  log={np.log(p_exonic_sense):.4f}")
print(f"  intergenic (gDNA):P={p_intergenic_sense:.4f}  log={np.log(p_intergenic_sense):.4f}")
print(f"  LL contrast (exonic vs gDNA):  {np.log(p_exonic_sense) - np.log(p_intergenic_sense):.4f}")
print(f"  LL contrast (spliced vs gDNA): {np.log(p_spliced_sense) - np.log(p_intergenic_sense):.4f}")
print()
print(f"Anti-sense-aligned fragment on + gene:")
print(f"  exonic_spliced:  P={p_spliced_anti:.4f}  log={np.log(p_spliced_anti):.4f}")
print(f"  exonic (diluted): P={p_exonic_anti:.4f}  log={np.log(p_exonic_anti):.4f}")
print(f"  intergenic (gDNA):P={p_intergenic_anti:.4f}  log={np.log(p_intergenic_anti):.4f}")
print(f"  LL contrast (exonic vs gDNA):  {np.log(p_exonic_anti) - np.log(p_intergenic_anti):.4f}")
print(f"  LL contrast (spliced vs gDNA): {np.log(p_spliced_anti) - np.log(p_intergenic_anti):.4f}")

sc.cleanup()

# Additional debugging: shadow init and unique counts
print()
print("=== SHADOW INIT / UNIQUE COUNTS ===")
import numpy as np
counter = pipeline_result.counter
print(f"unique_counts shape: {counter.unique_counts.shape}")
print(f"unique_counts:\n{counter.unique_counts}")
print(f"shadow_init: {counter.shadow_init}")
print(f"gdna_locus_counts: {counter.gdna_locus_counts}")
print()

# Show what classify_strand actually does  
from hulkrna.counter import ReadCounter
from hulkrna.types import Strand
cs_sense = ReadCounter.classify_strand(Strand.POS, Strand.POS, sm.exonic)
cs_anti = ReadCounter.classify_strand(Strand.NEG, Strand.POS, sm.exonic)
print(f"classify_strand(POS frag, POS gene, exonic model): {cs_sense}")
print(f"classify_strand(NEG frag, POS gene, exonic model): {cs_anti}")

# What fraction of exonic fragments are sense vs antisense?
t_counts_arr = counter.t_counts
print(f"\nt_counts array:\n{t_counts_arr}")
print(f"unique_counts array:\n{counter.unique_counts}")
print(f"Unique by column: {counter.unique_counts.sum(axis=0)}")

# Check if the exonic model's strand specificity is actually being used
# for EM scoring
print()
print("=== MODEL ROUTING FOR UNSPLICED ===")
from hulkrna.categories import CountCategory
model = sm.model_for_category(CountCategory.UNSPLICED)
print(f"model_for_category(UNSPLICED) → p_r1_sense={model.p_r1_sense:.4f}")
model_sa = sm.model_for_category(CountCategory.SPLICED_ANNOT)
print(f"model_for_category(SPLICED_ANNOT) → p_r1_sense={model_sa.p_r1_sense:.4f}")

# EM convergence details
print()
print("=== EM CONVERGENCE ===")
if hasattr(counter, '_converged_theta') and counter._converged_theta is not None:
    theta = counter._converged_theta
    print(f"converged theta shape: {theta.shape}")
    print(f"theta[transcript] = {theta[0]:.6f}")
    gdna_base = counter.gdna_base_index
    print(f"theta[gDNA shadow] = {theta[gdna_base]:.6f}")
    print(f"gdna_base_index = {gdna_base}")
    print(f"gdna_locus_counts = {counter.gdna_locus_counts}")
    print(f"gdna_em_count = {counter.gdna_em_count:.1f}")
    print(f"gdna_unique_count = {counter.gdna_unique_count:.1f}")
    print(f"gdna_total = {counter.gdna_total:.1f}")

# Show actual EM scoring for a few representative fragments
print()
print("=== EM SCORING: RNA model vs gDNA model ===")
rna_model = sm.model_for_category(CountCategory.UNSPLICED)
gDNA_model = sm.intergenic
# Sense-aligned on + gene
p_rna_sense = rna_model.strand_likelihood(Strand.POS, Strand.POS)
p_gdna_sense = gDNA_model.strand_likelihood(Strand.POS, Strand.POS)
print(f"Sense-aligned:   RNA P={p_rna_sense:.6f} (LL={np.log(p_rna_sense):.4f}) | gDNA P={p_gdna_sense:.6f} (LL={np.log(p_gdna_sense):.4f}) | ΔLL={np.log(p_rna_sense)-np.log(p_gdna_sense):.4f}")
# Antisense-aligned on + gene
p_rna_anti = rna_model.strand_likelihood(Strand.NEG, Strand.POS)
p_gdna_anti = gDNA_model.strand_likelihood(Strand.NEG, Strand.POS)
print(f"Anti-aligned:    RNA P={p_rna_anti:.6f} (LL={np.log(p_rna_anti):.4f}) | gDNA P={p_gdna_anti:.6f} (LL={np.log(p_gdna_anti):.4f}) | ΔLL={np.log(p_rna_anti)-np.log(p_gdna_anti):.4f}")
print(f"\nFor sense-aligned:  gDNA wins by {np.log(p_gdna_sense)-np.log(p_rna_sense):.2f} log units")
print(f"For anti-aligned:   RNA wins by {np.log(p_rna_anti)-np.log(p_gdna_anti):.2f} log units")

# Show the insert size contrast
print()
print("=== INSERT SIZE MODELS ===")
from hulkrna.insert_model import InsertSizeModels
ism = pipeline_result.insert_models
cat_model = ism.category_models.get(int(CountCategory.UNSPLICED), ism.global_model)
ig_model = ism.intergenic
for isize in [150, 200, 250, 300, 350, 400]:
    ll_rna = cat_model.log_likelihood(isize)
    ll_gdna = ig_model.log_likelihood(isize)
    print(f"  isize={isize}: RNA LL={ll_rna:.4f} | gDNA LL={ll_gdna:.4f} | Δ={ll_rna-ll_gdna:.4f}")
