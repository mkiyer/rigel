"""Quantify the impact of ts tag misinterpretation on fragment resolution."""
from hulkrna.sim import Scenario, SimConfig
from hulkrna.pipeline import run_pipeline
import tempfile
from pathlib import Path

tmp = Path(tempfile.mkdtemp())
sc = Scenario("ts_impact", genome_length=5000, seed=42, work_dir=tmp / "ts_impact")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (2000, 2300)],
     "abundance": 100},
])
r = sc.build(n_fragments=500, sim_config=SimConfig(seed=42))

# Run with ts tag (current behavior - alignment-relative, treated as reference)
pr_ts = run_pipeline(r.bam_path, r.index, sj_strand_tag="ts", seed=42)
# Run with auto (should detect ts since no XS)
pr_auto = run_pipeline(r.bam_path, r.index, sj_strand_tag="auto", seed=42)

print("=== Fragment category counts (sj_strand_tag='ts') ===")
print(f"  Total fragments: {pr_ts.stats.n_fragments}")
print(f"  SPLICED_ANNOT: {pr_ts.stats.n_spliced_annot}")
print(f"  SPLICED_UNANNOT: {pr_ts.stats.n_spliced_unannot}")
print(f"  UNSPLICED: {pr_ts.stats.n_unspliced}")
print(f"  INTRON: {pr_ts.stats.n_intron}")
print(f"  Intergenic: {pr_ts.stats.n_intergenic}")
print(f"  Chimeric: {pr_ts.stats.n_chimeric}")

print(f"\n=== Fragment category counts (sj_strand_tag='auto') ===")
print(f"  Total fragments: {pr_auto.stats.n_fragments}")
print(f"  SPLICED_ANNOT: {pr_auto.stats.n_spliced_annot}")
print(f"  SPLICED_UNANNOT: {pr_auto.stats.n_spliced_unannot}")
print(f"  UNSPLICED: {pr_auto.stats.n_unspliced}")
print(f"  INTRON: {pr_auto.stats.n_intron}")
print(f"  Intergenic: {pr_auto.stats.n_intergenic}")
print(f"  Chimeric: {pr_auto.stats.n_chimeric}")

# Strand model diagnostics
sm_ts = pr_ts.strand_models
sm_auto = pr_auto.strand_models
print(f"\n=== Strand models (ts) ===")
if sm_ts.exonic_spliced.n_observations > 0:
    print(f"  exonic_spliced: p_r1_sense={sm_ts.exonic_spliced.p_r1_sense:.4f} n={sm_ts.exonic_spliced.n_observations}")
else:
    print(f"  exonic_spliced: NO OBSERVATIONS")
print(f"  exonic: p_r1_sense={sm_ts.exonic.p_r1_sense:.4f} n={sm_ts.exonic.n_observations}")

print(f"\n=== Strand models (auto) ===")
if sm_auto.exonic_spliced.n_observations > 0:
    print(f"  exonic_spliced: p_r1_sense={sm_auto.exonic_spliced.p_r1_sense:.4f} n={sm_auto.exonic_spliced.n_observations}")
else:
    print(f"  exonic_spliced: NO OBSERVATIONS")
print(f"  exonic: p_r1_sense={sm_auto.exonic.p_r1_sense:.4f} n={sm_auto.exonic.n_observations}")

# Also check: what does "auto" detect?
from hulkrna.bam import detect_sj_strand_tag
detected = detect_sj_strand_tag(r.bam_path)
print(f"\n=== Auto-detection result ===")
print(f"  Detected tags: {detected}")

sc.cleanup()
