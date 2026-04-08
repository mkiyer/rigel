#!/usr/bin/env python3
"""Deep diagnostic: nRNA-gDNA misclassification root cause analysis.

Reproduces the worst failure: NTA=500, low mRNA, SS=0.5, no gDNA,
5000 RNA fragments → all nRNA misclassified as gDNA.
"""
import logging
import sys
import tempfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rigel.config import PipelineConfig, EMConfig
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline, scan_and_buffer, quant_from_buffer

logging.basicConfig(level=logging.INFO, format="%(levelname)-5s %(name)s — %(message)s")
logger = logging.getLogger("diag")


def build_scenario(n_rna_fragments=5000, ss=0.5, gdna_fraction=0.0,
                   nta=500, mta1=20, mta2=5, mta3=10, mta4=2):
    """Build the pathological scenario."""
    sc = Scenario("diag", genome_length=50000, seed=42,
                  work_dir=Path(tempfile.mkdtemp(prefix="rigel_diag_")))

    sc.add_gene("NTA", "+", [
        {"t_id": "TA1", "exons": [(1000, 2000), (5000, 5500), (7000, 7500), (9000, 10000)],
         "abundance": mta1, "nrna_abundance": nta},
        {"t_id": "TA2", "exons": [(1000, 2000), (9000, 10000)],
         "abundance": mta2, "nrna_abundance": 0},
        {"t_id": "TA3", "exons": [(1000, 2000), (5000, 5500), (9000, 10000)],
         "abundance": mta3, "nrna_abundance": 0},
        {"t_id": "TA4", "exons": [(4500, 5500), (9000, 10000)],
         "abundance": mta4, "nrna_abundance": 0},
    ])

    sim_cfg = SimConfig(
        frag_mean=250, frag_std=50, frag_min=80, frag_max=500,
        read_length=150, strand_specificity=ss, seed=42,
    )

    gdna_cfg = None
    gdna_frac = gdna_fraction
    if gdna_fraction > 0:
        gdna_cfg = GDNAConfig(abundance=1.0, frag_mean=350, frag_std=100,
                              frag_min=100, frag_max=1000)

    result = sc.build_oracle(
        n_fragments=10000,  # ignored in n_rna_fragments mode
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        nrna_abundance=0.0,
        n_rna_fragments=n_rna_fragments,
        gdna_fraction=gdna_frac,
    )
    return sc, result


def diagnose_case(label, n_rna=5000, ss=0.5, gdna_frac=0.0, nta=500,
                  mta1=20, mta2=5, mta3=10, mta4=2):
    """Run one case and inspect every stage."""
    print(f"\n{'='*70}")
    print(f"  CASE: {label}")
    print(f"  n_rna={n_rna}, SS={ss}, gdna_frac={gdna_frac}, NTA={nta}")
    print(f"  mRNA: TA1={mta1}, TA2={mta2}, TA3={mta3}, TA4={mta4}")
    print(f"{'='*70}")

    sc, result = build_scenario(n_rna, ss, gdna_frac, nta, mta1, mta2, mta3, mta4)

    # Ground truth
    gt = result.ground_truth_auto()
    gt_gdna = result.ground_truth_gdna_auto()
    gt_nrna = result.ground_truth_nrna_auto()
    gt_mrna = sum(gt.values())
    print(f"\n  Ground truth: mRNA={gt_mrna}, nRNA={gt_nrna}, gDNA={gt_gdna}")
    print(f"  Total fragments: {gt_mrna + gt_nrna + gt_gdna}")

    # Stage 1: BAM scan
    pipe_cfg = PipelineConfig()
    bam_str = str(result.bam_path)
    stats, strand_models, frag_length_models, buffer, region_counts, fl_table = \
        scan_and_buffer(bam_str, result.index, scan=pipe_cfg.scan)

    frag_length_models.build_scoring_models()

    print(f"\n  --- BAM SCAN ---")
    print(f"  Total fragments scanned: {stats.n_fragments}")
    print(f"  Unique mappers: {stats.unique}")
    print(f"  Multimappers: {stats.multimapping}")

    # Strand model
    sm = strand_models
    print(f"\n  --- STRAND MODEL ---")
    print(f"  Trained SS: {sm.strand_specificity:.4f}")
    print(f"  Spliced obs: {sm.exonic_spliced.n_observations}")
    print(f"  Exonic diag SS: {sm.exonic.strand_specificity:.4f} "
          f"({sm.exonic.n_observations} obs)")
    print(f"  Intergenic diag SS: {sm.intergenic.strand_specificity:.4f} "
          f"({sm.intergenic.n_observations} obs)")

    # Fragment length models
    fl = frag_length_models
    print(f"\n  --- FRAGMENT LENGTH MODELS ---")
    print(f"  RNA model: {fl.rna_model.n_observations} obs, "
          f"mean={fl.rna_model.mean:.1f}")
    print(f"  gDNA model: {fl.gdna_model.n_observations} obs, "
          f"mean={fl.gdna_model.mean:.1f}" if fl.gdna_model.n_observations > 0
          else f"  gDNA model: 0 obs")

    # Region counts — this is what calibration sees
    print(f"\n  --- REGION COUNTS ---")
    print(f"  Regions: {len(region_counts)}")
    print(f"  Columns: {list(region_counts.columns)}")

    # Get region geometry from index
    region_df = result.index.region_df

    # Show per-region breakdown
    for i, (_, row) in enumerate(region_counts.iterrows()):
        n_up = row.get("n_unspliced_pos", 0)
        n_un = row.get("n_unspliced_neg", 0)
        n_sp = row.get("n_spliced_pos", 0)
        n_sn = row.get("n_spliced_neg", 0)
        n_unspliced = n_up + n_un
        n_spliced = n_sp + n_sn
        n_total = n_unspliced + n_spliced

        # Get region geometry
        if region_df is not None and i < len(region_df):
            rr = region_df.iloc[i]
            ref = rr.get("ref", "?")
            start = rr.get("start", 0)
            end = rr.get("end", 0)
            length = end - start
            strand = rr.get("strand", "?")
        else:
            ref, start, end, length, strand = "?", 0, 0, 0, "?"

        density = n_unspliced / length if length > 0 else 0
        spliced_frac = n_spliced / n_total if n_total > 0 else 0
        print(f"    Region {i}: {ref}:{start}-{end}({strand}) len={length} "
              f"total={n_total:.0f} unspliced={n_unspliced:.0f} "
              f"spliced={n_spliced:.0f} ({spliced_frac:.2f}) "
              f"density={density:.4f}")

    # Stage 2: Calibration
    from rigel.calibration import calibrate_gdna

    region_df = result.index.region_df

    print(f"\n  --- CALIBRATION INPUT ---")
    print(f"  Region DF shape: {region_df.shape}")
    print(f"  Region DF columns: {list(region_df.columns)}")

    cal = calibrate_gdna(
        region_counts=region_counts,
        fl_table=fl_table,
        region_df=region_df,
        strand_specificity=sm.strand_specificity,
        density_percentile=pipe_cfg.calibration.density_percentile,
    )

    print(f"\n  --- CALIBRATION OUTPUT ---")
    print(f"  λ_G (global gDNA density): {cal.lambda_gdna:.6e}")
    print(f"  SS (strand specificity): {cal.strand_specificity:.4f}")
    total_e_gdna = float(cal.region_e_gdna.sum())
    total_n = float(cal.region_n_total.sum())
    gamma = total_e_gdna / max(total_n, 1.0)
    print(f"  Total E[gDNA]: {total_e_gdna:.1f}")
    print(f"  Total N: {total_n:.1f}")
    print(f"  γ (gDNA fraction): {gamma:.4f}")

    # Per-region E[gDNA]
    print(f"\n  Per-region E[gDNA]:")
    for i in range(len(cal.region_e_gdna)):
        e_gdna = cal.region_e_gdna[i]
        n_total = cal.region_n_total[i]
        frac = e_gdna / max(n_total, 1.0)
        if n_total > 0:
            print(f"    Region {i}: E[gDNA]={e_gdna:.1f}/{n_total:.0f} "
                  f"({frac:.3f})")

    # Stage 3: Run full pipeline and check EM output
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    bench = run_benchmark(result, pr, scenario_name="diag")

    print(f"\n  --- PIPELINE OUTPUT ---")
    print(f"  mRNA total: exp={bench.total_expected}, obs={bench.total_observed:.0f}")
    print(f"  nRNA: exp={bench.n_nrna_expected}, obs={bench.n_nrna_pipeline:.0f}")
    print(f"  gDNA: exp={bench.n_gdna_expected}, obs={bench.n_gdna_pipeline:.0f}")
    print(f"  mRNA rel_err: {bench.total_expected and abs(bench.total_observed - bench.total_expected) / bench.total_expected:.3f}")

    for ta in bench.transcripts:
        if ta.expected > 0 or ta.observed > 0:
            print(f"    {ta.t_id}: exp={ta.expected}, obs={ta.observed:.1f}, "
                  f"err={ta.rel_error:.3f}")

    sc.cleanup()
    return cal, bench


def main():
    # Case A: The worst case — no gDNA at all, yet calibration creates γ=0.937
    diagnose_case("NO_GDNA_SS05",
                  n_rna=5000, ss=0.5, gdna_frac=0.0, nta=500)

    # Case B: Same but SS=0.7
    diagnose_case("NO_GDNA_SS07",
                  n_rna=5000, ss=0.7, gdna_frac=0.0, nta=500)

    # Case C: Same but SS=0.9 — should work better
    diagnose_case("NO_GDNA_SS09",
                  n_rna=5000, ss=0.9, gdna_frac=0.0, nta=500)

    # Case D: Same but SS=1.0 — should be fine
    diagnose_case("NO_GDNA_SS10",
                  n_rna=5000, ss=1.0, gdna_frac=0.0, nta=500)

    # Case E: Control — no nRNA, no gDNA, SS=0.5. Should be clean.
    diagnose_case("CONTROL_NO_NRNA_SS05",
                  n_rna=5000, ss=0.5, gdna_frac=0.0, nta=0)

    # Case F: SS=0.5, with actual gDNA, no nRNA — calibration should detect gDNA
    diagnose_case("GDNA_ONLY_SS05",
                  n_rna=5000, ss=0.5, gdna_frac=0.3, nta=0)


if __name__ == "__main__":
    main()
