#!/usr/bin/env python3
"""
A/B test: equiv-class sorting fix in em_solver.cpp.

This test creates a 20kb two-transcript (opposite strand) scenario with 10,000
fragments including gDNA and nRNA, then runs the pipeline on the ORIGINAL BAM
plus N_SHUFFLES shuffled copies.  It measures the maximum absolute and relative
difference across all pairs of runs.

If the equiv-class sorting makes the EM deterministic with respect to fragment
order, ALL runs (original + shuffles) should produce BIT-EXACT results.

Run:  python -m pytest tests/test_equiv_class_sort_ab.py -s
  or: python tests/test_equiv_class_sort_ab.py
"""

import logging
import tempfile
from pathlib import Path

import numpy as np
import pysam

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
N_FRAGMENTS = 10_000
N_SHUFFLES = 5
SIM_SEED = 42
PIPELINE_SEED = 42
GENOME_LENGTH = 20_000

# Test configs: (gdna_abundance, nrna_abundance, strand_specificity, label)
CONFIGS = [
    (0, 0, 1.0, "mrna_only_ss1"),
    (0, 0, 0.9, "mrna_only_ss09"),
    (0, 50, 0.9, "nrna50_ss09"),
    (50, 0, 0.9, "gdna50_ss09"),
    (50, 50, 0.9, "full_stress"),
    (100, 100, 0.8, "high_contam"),
]


def _sim_config(*, strand_specificity: float = 0.9, seed: int = SIM_SEED):
    return SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=seed,
    )


def _gdna_config(abundance: float):
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def _build_scenario(tmp_path, name="ab_test"):
    sc = Scenario(
        name,
        genome_length=GENOME_LENGTH,
        seed=SIM_SEED,
        work_dir=tmp_path / name,
    )
    sc.add_gene(
        "G1",
        "+",
        [
            {
                "t_id": "T1",
                "exons": [(1000, 1500), (3000, 3500), (8000, 10000)],
                "abundance": 100,
            },
        ],
    )
    sc.add_gene(
        "G2",
        "-",
        [
            {
                "t_id": "T2",
                "exons": [(6000, 7000), (12000, 13000), (15000, 15500)],
                "abundance": 100,
            },
        ],
    )
    return sc


def _shuffle_bam(input_bam: Path, output_bam: Path, seed: int):
    """Shuffle read-name groups in a name-sorted BAM."""
    groups: dict[str, list] = {}
    with pysam.AlignmentFile(str(input_bam), "rb") as inf:
        header = inf.header.to_dict()
        for rec in inf:
            qname = rec.query_name
            if qname not in groups:
                groups[qname] = []
            groups[qname].append(rec)

    rng = np.random.default_rng(seed)
    group_keys = list(groups.keys())
    rng.shuffle(group_keys)

    header["HD"]["SO"] = "queryname"
    with pysam.AlignmentFile(str(output_bam), "wb", header=header) as outf:
        for key in group_keys:
            for rec in groups[key]:
                outf.write(rec)


def _extract_results(pr, result):
    """Extract key metrics from a pipeline run."""
    est = pr.estimator

    info = {}
    for t in result.transcripts:
        ti = t.t_index
        tid = t.t_id
        info[f"{tid}_mrna"] = float(est.em_counts[ti].sum())
        info[f"{tid}_nrna"] = float(est.em_counts[ti].sum()) if est._nrna_mask[ti] else 0.0
        info[f"{tid}_unambig"] = float(est.unambig_counts[ti].sum())

    info["gdna_total"] = float(est.gdna_em_count)
    info["gdna_em"] = float(est.gdna_em_count)
    info["n_loci"] = len(est.locus_results)

    for i, lr in enumerate(est.locus_results):
        if isinstance(lr, dict):
            info[f"locus_{i}_n_transcripts"] = lr.get("n_transcripts", 0)
            info[f"locus_{i}_n_em_fragments"] = lr.get("n_em_fragments", 0)
            info[f"locus_{i}_alpha_gdna"] = lr.get("alpha_gdna", 0.0)
        else:
            info[f"locus_{i}_n_transcripts"] = lr.n_transcripts
            info[f"locus_{i}_n_em_fragments"] = lr.n_em_fragments
            info[f"locus_{i}_alpha_gdna"] = lr.alpha_gdna

    return info


def _compare(r1, r2):
    """Return list of (key, v1, v2, abs_diff, rel_diff) for differing keys."""
    all_keys = sorted(set(r1.keys()) | set(r2.keys()))
    diffs = []
    for k in all_keys:
        v1 = r1.get(k, float("nan"))
        v2 = r2.get(k, float("nan"))
        if v1 != v2:
            ad = abs(v1 - v2)
            denom = max(abs(v1), abs(v2), 1e-30)
            rd = ad / denom
            diffs.append((k, v1, v2, ad, rd))
    return diffs


def run_ab_test(
    tmp_path, *, gdna_abundance, nrna_abundance, strand_specificity, label, n_fragments=N_FRAGMENTS
):
    """Run original + N_SHUFFLES shuffled BAMs, return all pairwise diffs."""
    sc = _build_scenario(tmp_path, name=label)
    sim = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)

    result = sc.build_oracle(
        n_fragments=n_fragments,
        sim_config=sim,
        gdna_config=gdna,
        nrna_abundance=nrna_abundance,
    )

    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )

    # Run on original
    runs = {}
    pr0 = run_pipeline(result.bam_path, result.index, config=config)
    runs["original"] = _extract_results(pr0, result)

    # Run on N shuffled copies
    for i in range(N_SHUFFLES):
        seed = 100 + i * 37  # spread seeds
        shuffled = tmp_path / label / f"shuffled_{i}.bam"
        _shuffle_bam(result.bam_path, shuffled, seed=seed)
        pr_s = run_pipeline(shuffled, result.index, config=config)
        runs[f"shuffle_{i}"] = _extract_results(pr_s, result)

    # Pairwise comparison
    run_names = sorted(runs.keys())
    all_diffs = []
    max_abs = 0.0
    max_rel = 0.0
    structural_diffs = []

    for i_r in range(len(run_names)):
        for j_r in range(i_r + 1, len(run_names)):
            n1, n2 = run_names[i_r], run_names[j_r]
            diffs = _compare(runs[n1], runs[n2])
            for k, v1, v2, ad, rd in diffs:
                max_abs = max(max_abs, ad)
                max_rel = max(max_rel, rd)
                if "n_loci" in k or "n_transcripts" in k or "n_em_fragments" in k:
                    structural_diffs.append((n1, n2, k, v1, v2))
            all_diffs.append((n1, n2, diffs))

    sc.cleanup()
    return runs, all_diffs, max_abs, max_rel, structural_diffs


def main():
    logging.basicConfig(level=logging.WARNING)

    print("=" * 78)
    print("  EQUIV-CLASS SORTING A/B TEST")
    print(f"  N_FRAGMENTS={N_FRAGMENTS}, N_SHUFFLES={N_SHUFFLES}")
    print("=" * 78)

    overall_max_abs = 0.0
    overall_max_rel = 0.0
    any_structural = False

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        for gdna, nrna, ss, label in CONFIGS:
            print(f"\n--- Config: {label} (gdna={gdna}, nrna={nrna}, ss={ss}) ---")

            runs, all_diffs, max_abs, max_rel, structural_diffs = run_ab_test(
                tmp,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                label=label,
            )

            # Report
            n_pairs = len(all_diffs)
            n_diff_pairs = sum(1 for _, _, d in all_diffs if d)
            print(f"  Pairs compared: {n_pairs}")
            print(f"  Pairs with ANY difference: {n_diff_pairs}")
            print(f"  Max abs diff: {max_abs:.2e}")
            print(f"  Max rel diff: {max_rel:.2e}")

            if structural_diffs:
                print(f"  ** STRUCTURAL DIFFERENCES: {len(structural_diffs)} **")
                for n1, n2, k, v1, v2 in structural_diffs:
                    print(f"     {n1} vs {n2}: {k} = {v1} vs {v2}")
                any_structural = True
            else:
                print("  No structural differences.")

            # Show worst diffs by abs magnitude
            worst = []
            for n1, n2, diffs in all_diffs:
                for k, v1, v2, ad, rd in diffs:
                    worst.append((ad, rd, k, n1, n2, v1, v2))
            worst.sort(reverse=True)
            if worst:
                print("  Top-5 worst diffs:")
                for ad, rd, k, n1, n2, v1, v2 in worst[:5]:
                    print(
                        f"    {k}: {v1:.10f} vs {v2:.10f} "
                        f"(abs={ad:.2e}, rel={rd:.2e}) [{n1} vs {n2}]"
                    )

            overall_max_abs = max(overall_max_abs, max_abs)
            overall_max_rel = max(overall_max_rel, max_rel)

    print("\n" + "=" * 78)
    print("  OVERALL SUMMARY")
    print(f"  Max absolute diff across all configs: {overall_max_abs:.2e}")
    print(f"  Max relative diff across all configs: {overall_max_rel:.2e}")
    if any_structural:
        print("  ** WARNING: STRUCTURAL DIFFERENCES FOUND **")
    else:
        print("  No structural differences in any config.")

    if overall_max_abs == 0.0:
        print("  VERDICT: PERFECTLY BIT-EXACT across all shuffles and configs.")
    elif overall_max_abs < 1e-10:
        print("  VERDICT: ULP-level noise only. Effectively deterministic.")
    elif overall_max_abs < 1e-3:
        print("  VERDICT: Small but non-trivial differences. FP accumulation noise.")
    else:
        print("  VERDICT: LARGE DIFFERENCES. Fragment order sensitivity detected!")
    print("=" * 78)


if __name__ == "__main__":
    main()
