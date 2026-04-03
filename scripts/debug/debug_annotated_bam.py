#!/usr/bin/env python3
"""Debug script to verify annotation table completeness.

Runs a small pipeline with annotation and checks how many
fragments get annotations from each source (det-unambig, chimeric, EM).
"""

import sys
import numpy as np
import pysam
import tempfile
from pathlib import Path

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import SimConfig, Scenario
from rigel.annotate import (
    AF_RESOLVED,
    AF_GDNA,
    AF_NRNA,
    AF_SYNTHETIC,
    AF_TRANSCRIPT,
    AF_GDNA_RESOLVED,
    AF_NRNA_RESOLVED,
    AF_SYNTH_RESOLVED,
    AF_UNRESOLVED,
)


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Build a scenario similar to test_annotate.py
        sc = Scenario(
            "debug_annot", genome_length=8000, seed=42,
            work_dir=tmpdir / "scenario",
        )
        # Spliced gene for strand training + deterministic unambig
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 1000), (1500, 2000)],
             "abundance": 80},
        ])
        # Second gene for EM-routed isoform ambiguity
        sc.add_gene("g2", "+", [
            {"t_id": "t2",
             "exons": [(3000, 3500), (4000, 4500)],
             "abundance": 40},
            {"t_id": "t3",
             "exons": [(3000, 3500), (4000, 4500), (5000, 5500)],
             "abundance": 40},
        ])

        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
            read_length=100, strand_specificity=1.0, seed=42,
        )
        result = sc.build(
            n_fragments=500, sim_config=sim_config,
        )

        annotated_bam = tmpdir / "annotated.bam"
        pr = run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=42),
                scan=BamScanConfig(sj_strand_tag="ts"),
                annotated_bam_path=annotated_bam,
            ),
        )

        # Check annotation table
        print("\n=== Pipeline Stats ===")
        stats = pr.stats
        print(f"  Total fragments: {stats.n_fragments}")
        print(f"  Det-unambig: {stats.deterministic_unambig_units}")
        print(f"  EM-routed unambig: {stats.em_routed_unambig_units}")
        print(f"  EM-routed ambig (same strand): {stats.em_routed_ambig_same_strand_units}")
        print(f"  EM-routed ambig (opp strand): {stats.em_routed_ambig_opp_strand_units}")
        print(f"  Gated: {stats.n_gated_out}")
        print(f"  Multimapper groups: {stats.em_routed_multimapper_units}")

        # Count annotated BAM tags
        print("\n=== Annotated BAM ===")
        bam = pysam.AlignmentFile(str(annotated_bam), "rb")
        records = list(bam.fetch(until_eof=True))
        bam.close()

        flag_labels = {
            AF_UNRESOLVED: "unresolved",
            AF_TRANSCRIPT: "transcript",
            AF_GDNA_RESOLVED: "gDNA",
            AF_NRNA_RESOLVED: "nRNA",
            AF_SYNTH_RESOLVED: "synthetic_nRNA",
        }

        flag_counts = {}
        fclass_counts = {}
        for rec in records:
            flags = rec.get_tag("ZF")
            label = flag_labels.get(flags, f"unknown({flags})")
            flag_counts[label] = flag_counts.get(label, 0) + 1
            zc = rec.get_tag("ZC")
            fclass_counts[zc] = fclass_counts.get(zc, 0) + 1

        print(f"  Total records: {len(records)}")
        print(f"  Assignment flags distribution:")
        for k, v in sorted(flag_counts.items(), key=lambda x: -x[1]):
            print(f"    {k}: {v} ({v/len(records)*100:.1f}%)")
        print(f"  Fragment class distribution:")
        for k, v in sorted(fclass_counts.items(), key=lambda x: -x[1]):
            print(f"    {k}: {v} ({v/len(records)*100:.1f}%)")

        # Count unique qnames per assignment flag category
        qname_flags = {}
        for rec in records:
            qn = rec.query_name
            flags = rec.get_tag("ZF")
            label = flag_labels.get(flags, f"unknown({flags})")
            if qn not in qname_flags:
                qname_flags[qn] = set()
            qname_flags[qn].add(label)

        flag_frag_counts = {}
        for qn, labels in qname_flags.items():
            for lbl in labels:
                flag_frag_counts[lbl] = flag_frag_counts.get(lbl, 0) + 1

        print(f"\n  Fragments (unique qnames) per assignment flag category:")
        total_frags = len(qname_flags)
        for k, v in sorted(flag_frag_counts.items(), key=lambda x: -x[1]):
            print(f"    {k}: {v} ({v/total_frags*100:.1f}%)")

        # Expected: det-unambig + EM-routed + chimeric = resolved
        expected_annotated = (
            stats.deterministic_unambig_units
            + stats.em_routed_unambig_units
            + stats.em_routed_ambig_same_strand_units
            + stats.em_routed_ambig_opp_strand_units
            + stats.em_routed_multimapper_units
        )
        actual_annotated = sum(
            v for k, v in flag_frag_counts.items() if k != "unresolved"
        )
        print(f"\n  Expected annotated fragments: {expected_annotated}")
        print(f"  Actual annotated fragments: {actual_annotated}")

        if actual_annotated < expected_annotated * 0.9:
            print("\n  WARNING: Significant annotation loss detected!")
            print(f"  Missing: {expected_annotated - actual_annotated} fragments")
        else:
            print("\n  OK: Annotation counts match expected.")

        sc.cleanup()


if __name__ == "__main__":
    main()
