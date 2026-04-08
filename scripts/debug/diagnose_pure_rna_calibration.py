#!/usr/bin/env python3
"""Diagnose calibration behavior in pure RNA scenarios.

Runs the nrna_double_counting scenario with g=0, n=0 at various SS levels,
and dumps region-level calibration internals to understand where false gDNA
signal comes from.
"""
import sys
import logging
import tempfile
import numpy as np

logging.basicConfig(level=logging.WARNING, format="%(name)s %(message)s")
logger = logging.getLogger(__name__)

from rigel.sim.scenario import Scenario
from rigel.sim.reads import SimConfig, GDNAConfig
from rigel.pipeline import run_pipeline
from rigel.config import PipelineConfig, EMConfig, BamScanConfig

SIM_SEED = 42
PIPELINE_SEED = 123


def make_scenario(tmp_dir):
    """Reproduce the nrna_double_counting scenario."""
    sc = Scenario(
        "diag_pure_rna",
        genome_length=20000,
        seed=SIM_SEED,
        work_dir=tmp_dir,
    )
    sc.add_gene("g1", "+", [
        {
            "t_id": "t1",
            "exons": [(2000, 4000), (8000, 10000)],
            "abundance": 100,
        },
    ])
    sc.add_gene("g_ctrl", "-", [
        {
            "t_id": "t_ctrl",
            "exons": [(14000, 16000), (18000, 19000)],
            "abundance": 0,
        },
    ])
    return sc


def run_diagnostic(ss):
    """Run scenario and extract calibration details."""
    tmp_dir = tempfile.mkdtemp(prefix="diag_cal_")
    sc = make_scenario(tmp_dir)

    sim_cfg = SimConfig(
        read_length=150,
        frag_mean=250,
        frag_std=50,
        frag_min=150,
        frag_max=600,
        strand_specificity=ss,
    )

    oracle_result = sc.build_oracle(
        n_fragments=2000,
        sim_config=sim_cfg,
        gdna_config=None,
        nrna_abundance=0,
    )

    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(oracle_result.bam_path, oracle_result.index, config=config)

    cal = pr.calibration
    est = pr.estimator
    if cal is None:
        print(f"  SS={ss}: No calibration result!")
        return

    index = oracle_result.index
    region_df = index.region_df if hasattr(index, 'region_df') else None

    print(f"\n{'='*70}")
    print(f"SS={ss}  |  λ_gDNA={cal.lambda_gdna:.2e}  |  SS_measured={cal.strand_specificity:.4f}")
    print(f"{'='*70}")

    if cal.region_e_gdna is not None:
        n_regions = len(cal.region_e_gdna)
        n_nonzero_e = int(np.sum(cal.region_e_gdna > 0))
        n_nonzero_n = int(np.sum(cal.region_n_total > 0))
        n_zero_n = n_regions - n_nonzero_n
        total_e_gdna = float(cal.region_e_gdna.sum())
        total_n = float(cal.region_n_total.sum())
        gamma_global = total_e_gdna / max(total_n, 1.0)

        print(f"  Regions: {n_regions} total")
        print(f"  Regions with fragments (n_total > 0): {n_nonzero_n}")
        print(f"  Regions with ZERO fragments: {n_zero_n}")
        print(f"  Regions with E[gDNA] > 0: {n_nonzero_e}")
        print(f"  Total E[gDNA]: {total_e_gdna:.2f}")
        print(f"  Total N: {total_n:.0f}")
        print(f"  Global γ = E[gDNA]/N: {gamma_global:.6f}")

        # Show all regions
        if region_df is not None and len(region_df) == n_regions:
            # Print region_df columns for debugging
            print(f"\n  region_df columns: {list(region_df.columns)}")
            print(f"  {'Idx':>4} {'Ref':>4} {'Start':>7} {'End':>7} {'Len':>6} "
                  f"{'N_tot':>6} {'E_gDNA':>9} {'γ':>9} {'Notes'}")
            print(f"  {'-'*75}")
            for i in range(n_regions):
                row = region_df.iloc[i]
                n_tot = cal.region_n_total[i]
                e_g = cal.region_e_gdna[i]
                gamma_r = e_g / max(n_tot, 1e-12) if n_tot > 0 else 0.0
                ref = row.get('ref', row.get('chrom', '?'))
                start = int(row.get('start', 0))
                end = int(row.get('end', 0))
                length = end - start
                tx_p = bool(row.get('tx_pos', False))
                tx_n = bool(row.get('tx_neg', False))
                ex_p = bool(row.get('exon_pos', False))
                ex_n = bool(row.get('exon_neg', False))
                notes = []
                if n_tot == 0:
                    notes.append("ZERO_FRAGS")
                if e_g > 1e-6:
                    notes.append(f"E[gDNA]={e_g:.2f}")
                if tx_p and not tx_n:
                    notes.append("gene(+)")
                elif tx_n and not tx_p:
                    notes.append("gene(-)")
                elif tx_p and tx_n:
                    notes.append("both_strands")
                else:
                    notes.append("intergenic")
                if ex_p:
                    notes.append("exon+")
                if ex_n:
                    notes.append("exon-")
                note_str = ", ".join(notes)
                print(f"  {i:>4} {str(ref):>4} {start:>7} {end:>7} {length:>6} "
                      f"{n_tot:>6.0f} {e_g:>9.4f} {gamma_r:>9.6f} {note_str}")
        else:
            # Just count regions with/without fragments
            print(f"\n  (No region_df available, raw arrays only)")
            for i in range(min(n_regions, 30)):
                print(f"    region {i}: E[gDNA]={cal.region_e_gdna[i]:.4f}, "
                      f"N={cal.region_n_total[i]:.0f}")

    if cal.gdna_fl_model is not None:
        print(f"\n  gDNA FL model: mean={cal.gdna_fl_model.mean:.1f}, "
              f"n_obs={cal.gdna_fl_model.n_observations}")
    else:
        print(f"\n  gDNA FL model: None")

    # Print locus-level results
    if est is not None and hasattr(est, 'locus_results'):
        print(f"\n  Locus EM results ({len(est.locus_results)} loci):")
        for lr in est.locus_results:
            print(f"    locus {lr['locus_id']}: mrna={lr['mrna']:.1f}, "
                  f"gdna={lr['gdna']:.1f}, "
                  f"alpha_gdna={lr['alpha_gdna']:.4f}, "
                  f"alpha_rna={lr['alpha_rna']:.4f}, "
                  f"n_frags={lr['n_em_fragments']}")

    # Print per-transcript counts
    if est is not None:
        t_df = index.t_df
        print(f"\n  Transcript results:")
        for i in range(len(t_df)):
            t_id = t_df.iloc[i]['transcript_id']
            mrna = float(est.mrna_counts[i])
            nrna = float(est.nrna_counts[i]) if hasattr(est, 'nrna_counts') else 0.0
            print(f"    {t_id}: mrna={mrna:.1f}, nrna={nrna:.1f}")

    sc.cleanup()


if __name__ == "__main__":
    for ss in [1.0, 0.95, 0.90, 0.65, 0.50]:
        try:
            run_diagnostic(ss)
        except Exception as e:
            print(f"  SS={ss}: ERROR: {e}")
            import traceback
            traceback.print_exc()
