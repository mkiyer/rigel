"""Diagnose test_ss_gradient_unspliced failure after nRNA flag refactor."""
import sys
sys.path.insert(0, "src")

from rigel.sim import Scenario, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline
from rigel.config import PipelineConfig, EMConfig, BamScanConfig
import tempfile, pathlib

tmp = pathlib.Path(tempfile.mkdtemp())
sc = Scenario("diag", genome_length=5000, seed=42, work_dir=tmp / "diag")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
])

for ss in [0.9, 1.0]:
    result = sc.build_oracle(
        n_fragments=2000,
        sim_config=SimConfig(strand_specificity=ss),
    )
    gt = result.ground_truth_auto()
    pr = run_pipeline(
        result.bam_path, result.index,
        config=PipelineConfig(
            em=EMConfig(seed=99),
            scan=BamScanConfig(sj_strand_tag="ts"),
        ),
    )
    est = pr.estimator
    idx = result.index

    print(f"\n=== ss={ss} ===")
    print(f"Ground truth: {gt}")
    print(f"t_df columns: {list(idx.t_df.columns)}")
    print(f"t_df is_nrna: {idx.t_df['is_nrna'].tolist()}")
    print(f"t_df is_synthetic: {idx.t_df['is_synthetic'].tolist()}")
    print(f"num_transcripts: {idx.num_transcripts}")
    print(f"Pipeline stats:")
    stats = pr.stats
    print(f"  det_unambig: {stats.deterministic_unambig_units}")
    print(f"  em_routed_unambig: {stats.em_routed_unambig_units}")
    print(f"  em_routed_ambig_same: {stats.em_routed_ambig_same_strand_units}")
    print(f"  em_routed_ambig_opp: {stats.em_routed_ambig_opp_strand_units}")
    print(f"  gated_out: {stats.n_gated_out}")
    print(f"  n_frags: {stats.n_fragments}")
    print(f"  n_intergenic: {stats.n_intergenic}")

    t_counts = est.t_counts.sum(axis=1)
    for i in range(idx.num_transcripts):
        t_id = idx.t_df.iloc[i]["t_id"]
        is_nrna = idx.t_df.iloc[i]["is_nrna"]
        is_syn = idx.t_df.iloc[i]["is_synthetic"]
        unambig = est.unambig_counts[i].sum()
        em = est.em_counts[i].sum()
        total = t_counts[i]
        print(f"  t{i} ({t_id}): is_nrna={is_nrna}, is_synthetic={is_syn}, "
              f"unambig={unambig:.1f}, em={em:.1f}, total={total:.1f}")

    print(f"  gDNA em total: {est.gdna_em_count:.1f}")
    print(f"  nRNA em count (synthetic): {est.nrna_em_count:.1f}")

    bench = run_benchmark(result, pr, scenario_name=f"ss_{ss}")
    for t in bench.transcripts:
        print(f"  Benchmark {t.t_id}: expected={t.expected}, observed={t.observed:.1f}")

sc.cleanup()
