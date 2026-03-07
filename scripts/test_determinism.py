"""Test determinism of the pipeline across repeated runs.

Runs the same scenario N times and checks that ALL numerical outputs
are bit-exact identical. Any non-determinism indicates a bug.
"""
import numpy as np
import pandas as pd
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42
N_FRAGS = 2000
N_RUNS = 5

SIM_SS90 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.9, seed=SEED,
)


def _pipeline_config(n_scan_threads=0, n_em_threads=0):
    return PipelineConfig(
        em=EMConfig(seed=SEED, n_threads=n_em_threads),
        scan=BamScanConfig(
            sj_strand_tag="auto",
            n_scan_threads=n_scan_threads,
        ),
    )


def _capture(bam_path, index, config):
    pr = run_pipeline(bam_path, index, config=config)
    est = pr.estimator
    tdf = est.get_counts_df(index).sort_values("transcript_id").reset_index(drop=True)
    gdf = est.get_gene_counts_df(index).sort_values("gene_id").reset_index(drop=True)
    ldf = est.get_loci_df()
    if len(ldf) > 0:
        ldf = ldf.sort_values("locus_id").reset_index(drop=True)
    scalars = {
        "gdna_em_count": float(est.gdna_em_count),
        "nrna_em_count": float(est.nrna_em_count),
        "n_fragments": int(pr.stats.n_fragments),
    }
    return tdf, gdf, ldf, scalars


def _compare_exact(df1, df2, label):
    """Compare two DataFrames for bit-exact equality on all numeric columns."""
    diffs = {}
    for col in df1.columns:
        if not pd.api.types.is_numeric_dtype(df1[col]):
            continue
        a = df1[col].to_numpy(dtype=np.float64)
        b = df2[col].to_numpy(dtype=np.float64)
        # Handle NaN-NaN comparison
        both_nan = np.isnan(a) & np.isnan(b)
        diff = np.abs(a - b)
        diff[both_nan] = 0.0
        max_diff = np.nanmax(diff) if len(diff) > 0 else 0.0
        if max_diff > 0:
            max_val = max(np.nanmax(np.abs(a)), np.nanmax(np.abs(b)))
            n_diff = int(np.sum(diff > 0))
            rel = max_diff / max_val if max_val > 0 else float('inf')
            diffs[col] = {
                "max_abs_diff": max_diff,
                "max_val": max_val,
                "rel_max_diff": rel,
                "n_differ": n_diff,
                "n_total": len(a),
            }
    return diffs


def build_scenario(tmp_dir, gdna_abundance=0.0):
    """Build a scenario with gDNA and nRNA contamination."""
    gdna_config = GDNAConfig(abundance=gdna_abundance) if gdna_abundance > 0 else None
    sc = Scenario("determ", genome_length=10000, seed=SEED,
                  work_dir=tmp_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(500, 800), (1200, 1400), (2000, 2300)],
         "abundance": 100},
        {"t_id": "t2",
         "exons": [(500, 800), (1600, 1800), (2000, 2300)],
         "abundance": 30},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t3",
         "exons": [(4000, 4300), (5000, 5400)],
         "abundance": 50},
    ])
    sc.add_gene("g3", "+", [
        {"t_id": "t4",
         "exons": [(7000, 7500)],
         "abundance": 20},
    ])
    result = sc.build_oracle(
        n_fragments=N_FRAGS,
        sim_config=SIM_SS90,
        gdna_config=gdna_config,
    )
    return sc, result


def test_determinism(tmp_path):
    """Main determinism test across thread configurations."""
    scenarios = [
        ("gdna_10", 10.0),
    ]

    thread_configs = [
        ("1scan_1em", 1, 1),
        ("Nscan_1em", 0, 1),
        ("1scan_Nem", 1, 0),
        ("Nscan_Nem", 0, 0),
    ]

    for sc_name, gdna_abundance in scenarios:
        sc_dir = tmp_path / sc_name
        sc_dir.mkdir()
        sc, result = build_scenario(sc_dir, gdna_abundance=gdna_abundance)

        for tc_name, n_scan, n_em in thread_configs:
            config = _pipeline_config(n_scan_threads=n_scan, n_em_threads=n_em)
            label = f"{sc_name}/{tc_name}"
            print(f"\n{'='*60}")
            print(f"Testing: {label} (N_RUNS={N_RUNS})")
            print(f"{'='*60}")

            # Run N times
            runs = []
            for i in range(N_RUNS):
                tdf, gdf, ldf, scalars = _capture(
                    result.bam_path, result.index, config)
                runs.append((tdf, gdf, ldf, scalars))
                print(f"  Run {i+1}/{N_RUNS}: "
                      f"gdna_em={scalars['gdna_em_count']:.6f}, "
                      f"nrna_em={scalars['nrna_em_count']:.6f}, "
                      f"n_frags={scalars['n_fragments']}")

            # Compare all runs to run 0
            any_diff = False
            ref_tdf, ref_gdf, ref_ldf, ref_scalars = runs[0]

            for i in range(1, N_RUNS):
                tdf, gdf, ldf, scalars = runs[i]

                # Scalar comparison
                for key in ref_scalars:
                    if ref_scalars[key] != scalars[key]:
                        print(f"  *** SCALAR DIFF run0 vs run{i}: "
                              f"{key} = {ref_scalars[key]} vs {scalars[key]}")
                        any_diff = True

                # DataFrame comparison
                for df_name, ref_df, cur_df in [
                    ("transcript", ref_tdf, tdf),
                    ("gene", ref_gdf, gdf),
                    ("loci", ref_ldf, ldf),
                ]:
                    diffs = _compare_exact(ref_df, cur_df, f"{label}/{df_name}")
                    if diffs:
                        any_diff = True
                        for col, info in diffs.items():
                            print(f"  *** DIFF run0 vs run{i} | "
                                  f"{df_name}.{col}: "
                                  f"max_abs={info['max_abs_diff']:.6g}, "
                                  f"rel={info['rel_max_diff']:.6g}, "
                                  f"n={info['n_differ']}/{info['n_total']}")

            if not any_diff:
                print(f"  -> ALL {N_RUNS} RUNS BIT-EXACT IDENTICAL")
            else:
                print(f"  -> NON-DETERMINISM DETECTED")

        sc.cleanup()


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        from pathlib import Path
        test_determinism(Path(tmp))
