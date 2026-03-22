#!/usr/bin/env python3
"""Verify nRNA definition correctness.

Compares old nrna.feather spans against new synthetic transcripts.
With tolerance=0, the set of unique (ref, strand, start, end) spans
should be IDENTICAL to the old nRNA table.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

OLD_INDEX = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune/rigel_index")
NEW_INDEX = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune/phase_e_results/rigel_index")


def main():
    # Load old nRNA spans
    old_nrna = pd.read_feather(OLD_INDEX / "nrna.feather")
    print(f"Old nRNA table: {len(old_nrna)} entries")
    print(f"  Columns: {old_nrna.columns.tolist()}")
    print(f"  Sample:\n{old_nrna.head()}\n")

    old_spans = set()
    for _, row in old_nrna.iterrows():
        old_spans.add((row["ref"], int(row["strand"]), int(row["start"]), int(row["end"])))
    print(f"Old: {len(old_spans)} unique spans")

    # Load new transcripts (synthetics)
    new_t = pd.read_feather(NEW_INDEX / "transcripts.feather")
    synthetics = new_t[new_t["is_synthetic_nrna"]].copy()
    print(f"\nNew synthetic nRNA transcripts: {len(synthetics)}")

    new_spans = set()
    for _, row in synthetics.iterrows():
        new_spans.add((row["ref"], int(row["strand"]), int(row["start"]), int(row["end"])))
    print(f"New: {len(new_spans)} unique spans")

    # Also load annotated single-exon equivalents
    annotated_equiv = new_t[new_t["is_nascent_equiv"]].copy()
    print(f"Annotated equivalents (is_nascent_equiv=True): {len(annotated_equiv)}")

    # Annotated equiv spans
    equiv_spans = set()
    for _, row in annotated_equiv.iterrows():
        equiv_spans.add((row["ref"], int(row["strand"]), int(row["start"]), int(row["end"])))

    # Combined new spans (synthetic + annotated equiv)
    all_new_spans = new_spans | equiv_spans
    print(f"Combined new (synthetic + equiv): {len(all_new_spans)} unique spans")

    # Compare
    in_old_not_new = old_spans - all_new_spans
    in_new_not_old = all_new_spans - old_spans

    print(f"\nIn old but not new: {len(in_old_not_new)}")
    print(f"In new but not old: {len(in_new_not_old)}")
    print(f"In both:            {len(old_spans & all_new_spans)}")

    if in_old_not_new:
        print(f"\nSample of {min(10, len(in_old_not_new))} old-only spans:")
        for s in sorted(in_old_not_new)[:10]:
            print(f"  {s}")

    if in_new_not_old:
        print(f"\nSample of {min(10, len(in_new_not_old))} new-only spans:")
        for s in sorted(in_new_not_old)[:10]:
            print(f"  {s}")

    # With tolerance=0, the old system used exact spans.
    # The new system with tol=20 merges nearby spans.
    # Let's check: how many old spans differ from their nearest new span?
    print(f"\n--- Tolerance analysis ---")
    print(f"The new system uses tolerance=20bp for TSS/TES merging.")
    print(f"Old exact spans: {len(old_spans)}")
    print(f"New merged spans (synthetic only): {len(new_spans)}")
    print(f"Reduction: {len(old_spans) - len(new_spans)} fewer spans "
          f"({100*(1 - len(new_spans)/len(old_spans)):.1f}%)")

    # Now rebuild synthetics with tolerance=0 and compare
    from rigel.index import create_nrna_transcripts
    from rigel.transcript import Transcript

    # Load ALL transcripts from the new index
    all_transcripts = Transcript.read_gtf(
        "/Users/mkiyer/Downloads/rigel_runs/refs/human/genes_controls.gtf.gz",
        parse_mode="warn-skip",
    )
    print(f"\nLoaded {len(all_transcripts)} GTF transcripts")

    # Create with tol=0
    syn_tol0 = create_nrna_transcripts(all_transcripts, tolerance=0)
    tol0_spans = set()
    for syn in syn_tol0:
        tol0_spans.add((syn.ref, int(syn.strand), syn.start, syn.end))

    # Annotated equiv with tol=0
    equiv_tol0_spans = set()
    for t in all_transcripts:
        if t.is_nascent_equiv:
            equiv_tol0_spans.add((t.ref, int(t.strand), t.start, t.end))

    all_tol0 = tol0_spans | equiv_tol0_spans
    print(f"\nWith tolerance=0:")
    print(f"  Synthetic: {len(tol0_spans)}")
    print(f"  Annotated equiv: {len(equiv_tol0_spans)}")
    print(f"  Combined: {len(all_tol0)}")

    in_old_not_tol0 = old_spans - all_tol0
    in_tol0_not_old = all_tol0 - old_spans

    print(f"\n  In old but not tol0: {len(in_old_not_tol0)}")
    print(f"  In tol0 but not old: {len(in_tol0_not_old)}")
    print(f"  In both:             {len(old_spans & all_tol0)}")

    if in_old_not_tol0:
        print(f"\n  Sample old-only (tol=0):")
        for s in sorted(in_old_not_tol0)[:10]:
            print(f"    {s}")

    if in_tol0_not_old:
        print(f"\n  Sample tol0-only:")
        for s in sorted(in_tol0_not_old)[:10]:
            print(f"    {s}")


if __name__ == "__main__":
    main()
