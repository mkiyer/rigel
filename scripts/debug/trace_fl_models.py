"""Trace FL model states through the pipeline for the failing test scenario."""
import tempfile
import logging
from pathlib import Path
from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import Scenario, SimConfig
from rigel.calibration import calibrate_gdna
from rigel.splice import SpliceType

logging.basicConfig(level=logging.INFO)

with tempfile.TemporaryDirectory() as tmp:
    sc = Scenario("test", genome_length=20000, seed=42, work_dir=Path(tmp) / "t")
    sc.add_gene("g1", "+", [{
        "t_id": "t1",
        "exons": [(2000, 4000), (8000, 10000)],
        "abundance": 100,
    }])
    sc.add_gene("g_ctrl", "-", [{
        "t_id": "t_ctrl",
        "exons": [(14000, 16000), (18000, 19000)],
        "abundance": 0,
    }])
    sim_config = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=0.65, seed=42,
    )
    result = sc.build_oracle(
        n_fragments=500, sim_config=sim_config,
        gdna_config=None, nrna_abundance=70,
    )
    config = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    stats, strand_models, fl_models, buffer, region_counts, fl_table = (
        scan_and_buffer(str(result.bam_path), result.index, config.scan)
    )

    strand_models.finalize()
    fl_models.build_scoring_models()
    fl_models.finalize()

    print("\n=== FL Model Summary (post-finalize, pre-calibration) ===")
    print(f"  Global: {fl_models.global_model.total_weight:.0f} obs, "
          f"mean={fl_models.global_model.mean:.1f}")
    print(f"  Intergenic: {fl_models.intergenic.total_weight:.0f} obs")
    for st in SpliceType:
        m = fl_models.category_models[st]
        if m.total_weight > 0:
            print(f"  {st.name}: {m.total_weight:.0f} obs, mean={m.mean:.1f}")
        else:
            print(f"  {st.name}: 0 obs")
    print(f"  RNA model: {fl_models.rna_model.total_weight:.0f} obs, "
          f"mean={fl_models.rna_model.mean:.1f}")
    print(f"  gDNA model (pre-cal): weight={fl_models.gdna_model.total_weight:.0f}, "
          f"_log_prob is None={fl_models.gdna_model._log_prob is None}")

    # Calibrate
    calibration = calibrate_gdna(
        region_counts, fl_table, result.index.region_df,
        strand_models.strand_specificity,
        intergenic_fl_model=fl_models.intergenic,
    )

    print("\n=== Calibration gDNA FL model ===")
    gfl = calibration.gdna_fl_model
    if gfl is None:
        print("  MODEL IS NONE!")
    else:
        print(f"  total_weight={gfl.total_weight:.0f}")
        print(f"  _log_prob is None: {gfl._log_prob is None}")
        if gfl._log_prob is not None:
            print(f"  log_prob[200]={gfl._log_prob[200]:.4f}")

    # After pipeline assignment (with Dirichlet prior re-finalize)
    cal_gdna = calibration.gdna_fl_model
    from rigel.frag_length_model import FragmentLengthModel
    gdna_copy = FragmentLengthModel(max_size=cal_gdna.max_size)
    gdna_copy.counts = cal_gdna.counts.copy()
    gdna_copy._total_weight = cal_gdna._total_weight
    gdna_copy.finalize(prior_counts=fl_models.global_model.counts)
    fl_models.gdna_model = gdna_copy
    print("\n=== After pipeline assignment ===")
    gdna = fl_models.gdna_model
    if gdna is None:
        print("  gdna_model IS NONE — C++ gets no FL array → returns 0.0")
    else:
        print(f"  gdna_model weight={gdna.total_weight:.0f}")
        print(f"  _log_prob is None: {gdna._log_prob is None}")
        if gdna._log_prob is not None:
            print(f"  log_prob[200]={gdna._log_prob[200]:.4f}")
        else:
            print("  _log_prob is None → C++ gets None → returns 0.0")

    # Summary
    print("\n=== ASYMMETRY ANALYSIS ===")
    rna_lp = fl_models.rna_model._log_prob
    gdna_lp = gdna._log_prob if gdna is not None else None
    if rna_lp is not None:
        print(f"  RNA FL log_prob[200] = {rna_lp[200]:.4f}")
    if gdna_lp is not None:
        print(f"  gDNA FL log_prob[200] = {gdna_lp[200]:.4f}")
        print(f"  Difference = {rna_lp[200] - gdna_lp[200]:.4f}")
    else:
        print(f"  gDNA FL: returns 0.0 in C++")
        print(f"  Asymmetry = {rna_lp[200]:.4f} (RNA penalized, gDNA not)")

    sc.cleanup()
