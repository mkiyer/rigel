#!/usr/bin/env python3
"""Diagnostic: Measure ACTUAL burden from both unambig + ambig spliced components.

Captures the logged burden total from the real code path, and also computes
Component A (unambig only) to quantify the gap between the two.
"""
import sys
import tempfile
import logging
from pathlib import Path
import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from rigel.estimator import compute_unspliced_to_spliced_ratios


class BurdenCapture(logging.Handler):
    """Capture the logged burden total from compute_eb_gdna_priors."""
    def __init__(self):
        super().__init__()
        self.burden = None
    def emit(self, record):
        msg = record.getMessage()
        if "mRNA unspliced burden subtracted" in msg:
            # Extract number from "... subtracted: N total predicted phantom reads"
            import re
            m = re.search(r"subtracted:\s*([\d.]+)", msg)
            if m:
                self.burden = float(m.group(1))


def run_scenario(label, n_fragments, nrna_ab, gdna_ab, mrna_ab=100, ss=1.0, seed=42):
    tmpdir = tempfile.mkdtemp(prefix="diag_ba_")
    sc = Scenario("ba", genome_length=50000, seed=seed,
                   work_dir=Path(tmpdir) / "work")
    sc.add_gene("geneA", "+", [
        {"t_id": "TA1", "exons": [(1000, 2000), (5000, 5500), (7000, 7500), (9000, 10000)],
         "abundance": mrna_ab, "nrna_abundance": nrna_ab},
        {"t_id": "TA2", "exons": [(1000, 2000), (9000, 10000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
        {"t_id": "TA3", "exons": [(1000, 2000), (5000, 5500), (9000, 10000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
        {"t_id": "TA4", "exons": [(4500, 5500), (9000, 10000)],
         "abundance": 3.0, "nrna_abundance": 0.0},
    ])
    sc.add_gene("geneB", "-", [
        {"t_id": "TB", "exons": [(20000, 24000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
    ])
    sc.add_gene("geneC", "-", [
        {"t_id": "TC", "exons": [(30000, 31000), (35000, 36000)],
         "abundance": 0.0, "nrna_abundance": 0.0},
    ])

    gdna_cfg = None
    if gdna_ab > 0:
        gdna_cfg = GDNAConfig(abundance=gdna_ab, frag_mean=350, frag_std=100,
                              frag_min=100, frag_max=1000)

    sim_cfg = SimConfig(frag_mean=200, frag_std=30, frag_min=80,
                       frag_max=450, read_length=100,
                       strand_specificity=ss, seed=seed)

    result = sc.build_oracle(n_fragments=n_fragments,
                            sim_config=sim_cfg,
                            gdna_config=gdna_cfg,
                            nrna_abundance=0.0)

    # Install log handler to capture burden
    cap = BurdenCapture()
    locus_log = logging.getLogger("rigel.locus")
    locus_log.addHandler(cap)
    locus_log.setLevel(logging.DEBUG)

    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=42, strand_symmetry_kappa=4.0),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    bench = run_benchmark(result, pr, scenario_name="ba")
    locus_log.removeHandler(cap)

    # Compute R_t ratios and Component A-only burden
    ratios = compute_unspliced_to_spliced_ratios(
        result.index, pr.frag_length_models.global_model,
    )
    estimator = pr.estimator

    burden_a_only = 0.0
    for t_idx in range(len(ratios)):
        c_spliced = estimator.unambig_counts[int(t_idx), 2:6].sum()
        r_t = ratios[int(t_idx)]
        if np.isfinite(r_t) and c_spliced > 0:
            burden_a_only += c_spliced * r_t

    # Unspliced totals
    total_unspliced = float(estimator.transcript_unspliced_sense.sum() +
                           estimator.transcript_unspliced_antisense.sum())

    # Unambig spliced total
    total_unambig_spliced = float(estimator.unambig_counts[:, 2:6].sum())

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")
    print(f"  Oracle: mRNA={bench.total_expected}, nRNA={bench.n_nrna_expected}, gDNA={bench.n_gdna_expected}")
    print(f"  Pipeline: mRNA={bench.total_observed:.1f}, nRNA={bench.n_nrna_pipeline:.1f}, gDNA={bench.n_gdna_pipeline:.1f}")
    print(f"  mRNA rel err: {abs(bench.total_observed - bench.total_expected)/max(bench.total_expected, 1)*100:.1f}%")
    print()
    print(f"  Total unambig spliced reads: {total_unambig_spliced:.0f}")
    print(f"  Burden Component A (unambig × R_t):  {burden_a_only:.0f}")
    print(f"  Burden A+B (real code, from log):     {cap.burden}")
    if cap.burden is not None:
        print(f"  Ambig contribution (B):              {cap.burden - burden_a_only:.0f}")
    print(f"  Total unspliced pool:                {total_unspliced:.0f}")
    if cap.burden is not None:
        print(f"  Burden / unspliced pool:             {cap.burden / max(total_unspliced, 1) * 100:.1f}%")

    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    print("="*70)
    print("  DIAGNOSTIC: Actual burden from real code (Component A + B)")
    print("="*70)

    run_scenario("CLEAN: mRNA only", 50000, nrna_ab=0, gdna_ab=0)
    run_scenario("MODERATE gDNA=100", 50000, nrna_ab=0, gdna_ab=100)
    run_scenario("EXTREME gDNA=1000", 50000, nrna_ab=0, gdna_ab=1000)
    run_scenario("gDNA=100, nRNA=100", 50000, nrna_ab=100, gdna_ab=100)
    run_scenario("EXTREME gDNA=1000, nRNA=100", 50000, nrna_ab=100, gdna_ab=1000)
