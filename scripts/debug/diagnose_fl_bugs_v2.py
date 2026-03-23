#!/usr/bin/env python3
"""Diagnose fragment length bugs using the actual rigel pipeline.

Runs the pipeline on Oracle and Minimap2 BAMs and extracts FL diagnostics
from the trained fragment length model and the resolved fragments.

Usage:
    conda activate rigel
    python scripts/debug/diagnose_fl_bugs_v2.py
"""
import sys
import logging
import numpy as np
from collections import Counter
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, "src")

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType

# --- Paths ---
BENCH = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3")
INDEX = BENCH / "rigel_index"
COND = BENCH / "gdna_none_ss_0.95_nrna_none"
ORACLE_BAM = COND / "align_oracle" / "reads_namesort.bam"
MM2_BAM = COND / "align_minimap2" / "reads_namesort.bam"

# --- Load index ---
logger.info("Loading index...")
idx = TranscriptIndex.load(str(INDEX))

# Look up nRNA status
tx_df = idx.t_df
is_nrna = tx_df["is_synthetic_nrna"].values
n_nrna = is_nrna.sum()
n_annotated = len(tx_df) - n_nrna
logger.info(f"Index: {len(tx_df)} tx, {n_annotated} annotated, {n_nrna} nRNA synthetics")


def analyze_fl_model(fl_models, label):
    """Print FL model diagnostics."""
    print(f"\n--- {label}: Fragment Length Model (trained from unique mappers) ---")

    for st in SpliceType:
        model = fl_models.category_models.get(st)
        if model and model.total_weight > 0:
            mode_idx = int(np.argmax(model.counts))
            print(f"  {st.name}: n_obs={model.total_weight:.0f}, "
                  f"mean={model.mean:.1f}, mode_FL={mode_idx}, "
                  f"mode_count={model.counts[mode_idx]:.0f}")
            # Check for tail contamination (FL > 500)
            tail = model.counts[500:]
            tail_sum = tail.sum()
            if tail_sum > 0:
                print(f"    Tail FL>500: {tail_sum:.0f} obs "
                      f"({tail_sum / model.total_weight * 100:.2f}%)")
                # Find peaks in tail
                for i in range(len(tail)):
                    if tail[i] > 50:
                        print(f"      FL={i + 500}: {tail[i]:.0f}")

    gm = fl_models.global_model
    if gm.total_weight > 0:
        gmode = int(np.argmax(gm.counts))
        print(f"  GLOBAL: n_obs={gm.total_weight:.0f}, mean={gm.mean:.1f}, "
              f"mode={gmode}")

    # Build scoring models and check rna_model
    fl_models.build_scoring_models()
    rm = fl_models.rna_model
    if rm.total_weight > 0:
        rmode = int(np.argmax(rm.counts))
        print(f"  RNA_MODEL (for scoring): n_obs={rm.total_weight:.0f}, "
              f"mean={rm.mean:.1f}, mode={rmode}")
    else:
        print(f"  RNA_MODEL: NO OBSERVATIONS")


def analyze_buffer(buf, is_nrna, label):
    """Analyze buffer fragments for nRNA FL contamination and FL issues."""
    print(f"\n--- {label}: Buffer Fragment Analysis ---")

    # Counters
    nrna_frag_count = 0
    nrna_only_frag_count = 0
    nrna_fl_disagree = 0
    nrna_fl_inflated = 0
    spliced_annot_with_nrna = 0
    ufl_total = 0          # fragments with unambiguous FL (all cands agree)
    ufl_nrna_only = 0      # ... where only nRNA candidates present
    ufl_nrna_mixed = 0     # ... where nRNA + mRNA both present

    # Per-splice-type FL distribution trackers
    fl_by_type = {st: [] for st in SpliceType}
    fl_by_type_nrna_only = {st: [] for st in SpliceType}

    sample_traces = []
    total_frags = 0

    for chunk in buf.iter_chunks():
        n = chunk.size
        for i in range(n):
            total_frags += 1
            s = chunk.t_offsets[i]
            e = chunk.t_offsets[i + 1]
            t_inds = chunk.t_indices[s:e]
            fls = chunk.frag_lengths[s:e]
            st_val = chunk.splice_type[i]
            nh = chunk.num_hits[i]

            if len(t_inds) == 0:
                continue

            nrna_mask = is_nrna[t_inds]
            has_nrna = nrna_mask.any()
            has_mrna = (~nrna_mask).any()

            st = SpliceType(st_val)

            if has_nrna:
                nrna_frag_count += 1
                if not has_mrna:
                    nrna_only_frag_count += 1
                if st == SpliceType.SPLICED_ANNOT:
                    spliced_annot_with_nrna += 1

                # Compare nRNA vs mRNA FL values
                nrna_fls = fls[nrna_mask]
                mrna_fls = fls[~nrna_mask]
                valid_nrna = nrna_fls[nrna_fls > 0]
                valid_mrna = mrna_fls[mrna_fls > 0]
                if len(valid_nrna) > 0 and len(valid_mrna) > 0:
                    if valid_nrna[0] != valid_mrna[0]:
                        nrna_fl_disagree += 1
                    if valid_nrna[0] > 2 * valid_mrna[0]:
                        nrna_fl_inflated += 1

            # Replicate get_unique_frag_length logic: all valid FLs must agree
            valid_fls = fls[fls > 0]
            if len(valid_fls) > 0 and np.all(valid_fls == valid_fls[0]):
                ufl = int(valid_fls[0])
                # Only unique mappers matter for training (NH==1)
                if nh == 1:
                    ufl_total += 1
                    fl_by_type[st].append(ufl)
                    if has_nrna and not has_mrna:
                        ufl_nrna_only += 1
                        fl_by_type_nrna_only[st].append(ufl)
                    elif has_nrna and has_mrna:
                        ufl_nrna_mixed += 1

            # Collect sample traces
            if has_nrna and len(sample_traces) < 30:
                trace = {
                    "frag_idx": total_frags - 1,
                    "splice_type": st.name,
                    "nh": int(nh),
                    "n_candidates": len(t_inds),
                    "n_nrna": int(nrna_mask.sum()),
                    "n_mrna": int((~nrna_mask).sum()),
                    "nrna_fls": fls[nrna_mask].tolist(),
                    "mrna_fls": fls[~nrna_mask].tolist(),
                    "genomic_fp": int(chunk.genomic_footprint[i]),
                }
                sample_traces.append(trace)

    print(f"  Total fragments scanned: {total_frags}")
    print(f"  Fragments with nRNA in candidates: {nrna_frag_count}")
    print(f"  Fragments with nRNA ONLY:          {nrna_only_frag_count}")
    print(f"  nRNA FL disagrees with mRNA FL:    {nrna_fl_disagree}")
    print(f"  nRNA FL > 2x mRNA FL:              {nrna_fl_inflated}")
    print(f"  SPLICED_ANNOT with nRNA:           {spliced_annot_with_nrna}")
    print(f"  Unambiguous FL (unique mappers):   {ufl_total}")
    print(f"    nRNA-only contributors:          {ufl_nrna_only}")
    print(f"    nRNA+mRNA mixed:                 {ufl_nrna_mixed}")

    # FL distribution per splice type
    print(f"\n--- {label}: Unambiguous FL Distributions (unique mappers) ---")
    for st in SpliceType:
        arr = np.array(fl_by_type[st]) if fl_by_type[st] else np.array([])
        nrna_arr = np.array(fl_by_type_nrna_only[st]) if fl_by_type_nrna_only[st] else np.array([])
        if len(arr) == 0:
            continue
        n_big = int((arr > 500).sum())
        mode_val = Counter(arr).most_common(1)[0][0]
        print(f"  {st.name}: n={len(arr)}, mean={arr.mean():.1f}, "
              f"median={np.median(arr):.1f}, mode={mode_val}, "
              f"FL>500={n_big} ({n_big / len(arr) * 100:.2f}%)")
        if n_big > 0:
            big = arr[arr > 500]
            print(f"    Large FL: min={big.min()}, max={big.max()}, "
                  f"mean={big.mean():.0f}, n={len(big)}")
        if len(nrna_arr) > 0:
            print(f"    nRNA-only subset: n={len(nrna_arr)}, "
                  f"mean={nrna_arr.mean():.1f}, "
                  f"FL>500={int((nrna_arr>500).sum())}")

    # Sample traces
    if sample_traces:
        print(f"\n--- {label}: Sample nRNA Traces ---")
        for tr in sample_traces[:15]:
            print(f"  Frag#{tr['frag_idx']}: {tr['splice_type']} NH={tr['nh']} "
                  f"{tr['n_mrna']}mRNA+{tr['n_nrna']}nRNA "
                  f"mFL={tr['mrna_fls']} nFL={tr['nrna_fls']} "
                  f"gfp={tr['genomic_fp']}")


def run_scan(bam_path, label):
    """Run scan_and_buffer and extract FL diagnostics."""
    logger.info(f"\n{'=' * 70}")
    logger.info(f"Scanning {label}: {bam_path}")
    logger.info(f"{'=' * 70}")

    scan_cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True)
    stats, strand_models, fl_models, buf, region_counts, fl_table = \
        scan_and_buffer(str(bam_path), idx, scan=scan_cfg)

    logger.info(f"{label}: {buf.total_fragments} fragments buffered")

    analyze_fl_model(fl_models, label)
    analyze_buffer(buf, is_nrna, label)

    return stats, fl_models, buf


# --- Run ---
orc_stats, orc_fl, orc_buf = run_scan(ORACLE_BAM, "Oracle")
mm2_stats, mm2_fl, mm2_buf = run_scan(MM2_BAM, "Minimap2")

# --- Cross-comparison ---
print("\n" + "=" * 80)
print("CROSS-COMPARISON: Oracle vs Minimap2")
print("=" * 80)

for st in SpliceType:
    orc_m = orc_fl.category_models.get(st)
    mm2_m = mm2_fl.category_models.get(st)
    if orc_m and mm2_m and orc_m.total_weight > 0 and mm2_m.total_weight > 0:
        # Compare distributions
        orc_counts = orc_m.counts / orc_m.total_weight
        mm2_counts = mm2_m.counts / mm2_m.total_weight
        # KL divergence (where both > 0)
        mask = (orc_counts > 0) & (mm2_counts > 0)
        kl = np.sum(mm2_counts[mask] * np.log(mm2_counts[mask] / orc_counts[mask]))
        print(f"  {st.name}: KL(mm2||orc) = {kl:.4f}, "
              f"orc_mean={orc_m.mean:.1f} mm2_mean={mm2_m.mean:.1f} "
              f"orc_n={orc_m.total_weight:.0f} mm2_n={mm2_m.total_weight:.0f}")

print("\nDone.")
