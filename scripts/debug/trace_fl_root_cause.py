#!/usr/bin/env python3
"""Trace the exact root cause of extreme FL outliers in Oracle data.

For each Oracle SPLICED_ANNOT fragment with FL > FL_THRESHOLD, this script:
1. Extracts the fragment's read name (ground truth) from the BAM
2. Determines the true source transcript and true FL
3. Compares true FL to computed FL for each candidate transcript
4. Identifies why compute_frag_lengths produces the wrong value

Key insight: simulated fragments have frag_max=1000 and are drawn from
N(300,50) truncated at [50,1000]. Any computed FL > 1000 is definitively
a bug in compute_frag_lengths(). Even FL > 500 in Oracle data is extreme
(>4 sigma) and warrants investigation.

Usage:
    conda activate rigel
    python scripts/debug/trace_fl_root_cause.py
"""
import sys
import logging
import pysam
import numpy as np
from collections import defaultdict
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

FL_THRESHOLD = 500  # investigate Oracle fragments with FL above this

# --- Load index ---
logger.info("Loading index...")
idx = TranscriptIndex.load(str(INDEX))
tx_df = idx.t_df
is_nrna = tx_df["is_synthetic_nrna"].values
t_id_col = tx_df["t_id"].values
t_strand_col = tx_df["strand"].values

# Build t_id -> t_index lookup
t_id_to_idx = {tid: i for i, tid in enumerate(t_id_col)}

# Load SJ data for intron info
import pandas as pd
sj_df = pd.read_feather(INDEX / "sj.feather")
# Build per-transcript intron list from SJ DataFrame
tx_introns = defaultdict(list)
for _, row in sj_df.iterrows():
    tx_introns[row["t_index"]].append((row["start"], row["end"]))


def parse_oracle_readname(qname):
    """Parse ground truth from Oracle BAM read name.
    
    Format: {t_id}:{frag_start}-{frag_end}:{strand}:{idx}
    Returns: (t_id, frag_start, frag_end, strand_str, frag_idx) or None
    """
    parts = qname.split(":")
    if len(parts) < 4:
        return None
    if parts[0].startswith("gdna") or parts[0].startswith("nrna"):
        return None
    t_id = parts[0]
    try:
        coords = parts[1].split("-")
        frag_start = int(coords[0])
        frag_end = int(coords[1])
        strand_str = parts[2]
        frag_idx = int(parts[3]) if len(parts) > 3 else -1
    except (ValueError, IndexError):
        return None
    return (t_id, frag_start, frag_end, strand_str, frag_idx)


# --- Phase 1: Scan and buffer ---
logger.info("Scanning Oracle BAM...")
scan_cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True)
stats, strand_models, fl_models, buf, region_counts, fl_table = \
    scan_and_buffer(str(ORACLE_BAM), idx, scan=scan_cfg)

logger.info(f"Buffered {buf.total_fragments} fragments")

# --- Phase 2: Find outlier fragments ---
logger.info(f"Scanning buffer for FL outliers (FL > {FL_THRESHOLD})...")

outlier_frag_ids = []  # (chunk_idx, frag_in_chunk, max_fl, fls, t_inds, footprint, splice_type, nh)

chunk_idx = 0
for chunk in buf.iter_chunks():
    for i in range(chunk.size):
        s = chunk.t_offsets[i]
        e = chunk.t_offsets[i + 1]
        t_inds = chunk.t_indices[s:e]
        fls = chunk.frag_lengths[s:e]
        st = chunk.splice_type[i]
        nh = chunk.num_hits[i]
        fp = chunk.genomic_footprint[i]
        gstart = chunk.genomic_start[i]

        if len(t_inds) == 0:
            continue

        valid_fls = fls[fls > 0]
        if len(valid_fls) == 0:
            continue

        max_fl = int(valid_fls.max())
        if max_fl > FL_THRESHOLD:
            # Check if this would have entered FL training
            # (get_unique_frag_length_mrna: all non-nRNA candidates agree)
            mrna_fls = []
            for j in range(len(t_inds)):
                if not is_nrna[t_inds[j]] and fls[j] > 0:
                    mrna_fls.append(fls[j])

            ufl = -1
            if mrna_fls:
                if all(f == mrna_fls[0] for f in mrna_fls):
                    ufl = mrna_fls[0]

            outlier_frag_ids.append({
                "chunk_idx": chunk_idx,
                "frag_in_chunk": i,
                "frag_id": int(chunk.frag_id[i]),
                "max_fl": max_fl,
                "ufl": int(ufl),
                "fls": fls.tolist(),
                "t_inds": t_inds.tolist(),
                "footprint": int(fp),
                "genomic_start": int(gstart),
                "splice_type": SpliceType(st).name,
                "nh": int(nh),
                "enters_training": ufl > 0 and nh == 1,
            })
    chunk_idx += 1

logger.info(f"Found {len(outlier_frag_ids)} fragments with FL > {FL_THRESHOLD}")

# How many would enter training?
n_training = sum(1 for o in outlier_frag_ids if o["enters_training"])
logger.info(f"  Of these, {n_training} would enter FL training (unambiguous, NH=1)")

# --- Phase 3: Match outliers to BAM read names ---
# We need to scan the BAM to get read names for the outlier fragments.
# The buffer stores frag_id which is sequential. We need to map frag_id -> read name.

# Get the set of frag_ids we need
target_frag_ids = {o["frag_id"] for o in outlier_frag_ids}
logger.info(f"Looking up {len(target_frag_ids)} unique frag_ids in BAM...")

frag_id_to_qname = {}
with pysam.AlignmentFile(str(ORACLE_BAM), "rb") as bam:
    current_frag_id = 0
    last_qname = None
    for read in bam:
        qname = read.query_name
        if qname != last_qname:
            if last_qname is not None:
                current_frag_id += 1
            last_qname = qname

        if current_frag_id in target_frag_ids:
            frag_id_to_qname[current_frag_id] = qname

        # Early exit once we have all needed
        if len(frag_id_to_qname) >= len(target_frag_ids):
            break

logger.info(f"Matched {len(frag_id_to_qname)} / {len(target_frag_ids)} frag_ids")

# --- Phase 4: Analyze each outlier ---
print(f"\n{'='*90}")
print(f"ROOT CAUSE ANALYSIS: Oracle FL Outliers (FL > {FL_THRESHOLD})")
print(f"{'='*90}")
print(f"Total outliers: {len(outlier_frag_ids)}")
print(f"Would enter FL training: {n_training}")

# Categorize outliers
categories = defaultdict(list)
for o in outlier_frag_ids:
    frag_id = o["frag_id"]
    qname = frag_id_to_qname.get(frag_id, "UNKNOWN")
    truth = parse_oracle_readname(qname)

    if truth is None:
        categories["unmatched"].append(o)
        continue

    true_tid, true_start, true_end, true_strand, _ = truth
    true_fl = true_end - true_start
    true_t_idx = t_id_to_idx.get(true_tid, -1)

    o["qname"] = qname
    o["true_tid"] = true_tid
    o["true_t_idx"] = true_t_idx
    o["true_fl"] = true_fl
    o["true_start"] = true_start
    o["true_end"] = true_end

    # Check if true source is in candidate set
    if true_t_idx in o["t_inds"]:
        # True source IS a candidate. What FL did it get?
        pos = o["t_inds"].index(true_t_idx)
        computed_fl_for_true = o["fls"][pos]
        o["computed_fl_for_true"] = computed_fl_for_true

        if computed_fl_for_true > FL_THRESHOLD:
            categories["true_source_wrong_fl"].append(o)
        else:
            categories["true_source_correct_fl"].append(o)
    else:
        categories["true_source_missing"].append(o)

print(f"\n--- Category Breakdown ---")
for cat, items in sorted(categories.items()):
    print(f"  {cat}: {len(items)}")

# Analyze each category
for cat, items in sorted(categories.items()):
    if not items:
        continue
    print(f"\n{'='*70}")
    print(f"Category: {cat} ({len(items)} fragments)")
    print(f"{'='*70}")

    # Sort by max_fl descending
    items.sort(key=lambda x: x["max_fl"], reverse=True)

    for o in items[:20]:  # top 20
        print(f"\n  frag_id={o['frag_id']}, {o['splice_type']}, NH={o['nh']}, "
              f"fp={o['footprint']}, max_FL={o['max_fl']}")
        print(f"    Read name: {o.get('qname', 'UNKNOWN')}")
        if "true_tid" in o:
            print(f"    True source: {o['true_tid']} (t_idx={o['true_t_idx']})")
            print(f"    True transcript-coordinate FL: {o['true_fl']} "
                  f"(positions {o['true_start']}-{o['true_end']})")
            if "computed_fl_for_true" in o:
                print(f"    Computed FL for true source: {o['computed_fl_for_true']}")

        # Show per-candidate FL
        for j, t_idx in enumerate(o["t_inds"]):
            fl = o["fls"][j]
            tid = t_id_col[t_idx]
            is_true = t_idx == o.get("true_t_idx", -1)
            is_nr = bool(is_nrna[t_idx])
            marker = " <-- TRUE SOURCE" if is_true else ""
            nr_marker = " [nRNA]" if is_nr else ""
            n_introns = len(tx_introns.get(t_idx, []))
            print(f"      t_idx={t_idx} ({tid}): FL={fl}, "
                  f"n_introns={n_introns}{nr_marker}{marker}")

        # Show intron structure of candidates
        if "true_t_idx" in o and cat == "true_source_wrong_fl":
            true_t = o["true_t_idx"]
            print(f"    True source introns ({len(tx_introns.get(true_t, []))}):")
            for start, end in tx_introns.get(true_t, [])[:10]:
                print(f"      [{start}, {end}) = {end-start}bp")
            # Show gap region  
            gstart = o["genomic_start"]
            gend = gstart + o["footprint"]
            print(f"    Fragment genomic span: [{gstart}, {gend})")

# --- Phase 5: Analyze the training-entering outliers separately ---
training_outliers = [o for o in outlier_frag_ids if o["enters_training"] and "true_fl" in o]
if training_outliers:
    print(f"\n{'='*90}")
    print(f"CRITICAL: FL outliers that ENTER training ({len(training_outliers)})")
    print(f"{'='*90}")

    for o in sorted(training_outliers, key=lambda x: x["ufl"], reverse=True)[:30]:
        true_fl = o.get("true_fl", "?")
        computed_ufl = o["ufl"]
        error = int(computed_ufl) - int(true_fl) if isinstance(true_fl, int) else "?"
        print(f"  frag_id={o['frag_id']}: true_FL={true_fl}, computed_UFL={computed_ufl}, "
              f"error={error}, fp={o['footprint']}, {o['splice_type']}")
        print(f"    Source: {o.get('true_tid', '?')}")
        for j, t_idx in enumerate(o["t_inds"]):
            tid = t_id_col[t_idx]
            fl = o["fls"][j]
            is_true = t_idx == o.get("true_t_idx", -1)
            marker = " <-- TRUE" if is_true else ""
            print(f"      {tid}: FL={fl}{marker}")

print("\nDone.")
