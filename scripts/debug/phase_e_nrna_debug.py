#!/usr/bin/env python3
"""Debug nRNA creation to understand why only 62 synthetics with tol=0."""
import collections
import numpy as np

from rigel.transcript import Transcript
from rigel.index import create_nrna_transcripts, _cluster_coordinates

GTF = "/Users/mkiyer/Downloads/rigel_runs/refs/human/genes_controls.gtf.gz"

transcripts = Transcript.read_gtf(GTF, parse_mode="warn-skip")
for i, t in enumerate(transcripts):
    t.t_index = i

multi = [t for t in transcripts if len(t.exons) >= 2]
single = [t for t in transcripts if len(t.exons) == 1]
print(f"Total: {len(transcripts)}")
print(f"Multi-exon: {len(multi)}")
print(f"Single-exon: {len(single)}")

# Manually compute merged spans with tol=0
groups = collections.defaultdict(list)
for t in multi:
    groups[(t.ref, int(t.strand))].append(t)

t_to_merged = {}
for (ref, strand), grp in groups.items():
    starts = np.array([t.start for t in grp], dtype=np.int64)
    order_s = np.argsort(starts, kind="mergesort")
    sorted_starts = starts[order_s]
    cids_s = _cluster_coordinates(sorted_starts, 0)
    n_cs = cids_s[-1] + 1 if len(cids_s) else 0
    cid_min = np.full(n_cs, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(cid_min, cids_s, sorted_starts)
    rep_starts = np.empty(len(grp), dtype=np.int64)
    for i, oi in enumerate(order_s):
        rep_starts[oi] = cid_min[cids_s[i]]

    ends = np.array([t.end for t in grp], dtype=np.int64)
    order_e = np.argsort(ends, kind="mergesort")
    sorted_ends = ends[order_e]
    cids_e = _cluster_coordinates(sorted_ends, 0)
    n_ce = cids_e[-1] + 1 if len(cids_e) else 0
    cid_max = np.full(n_ce, np.iinfo(np.int64).min, dtype=np.int64)
    np.maximum.at(cid_max, cids_e, sorted_ends)
    rep_ends = np.empty(len(grp), dtype=np.int64)
    for i, oi in enumerate(order_e):
        rep_ends[oi] = cid_max[cids_e[i]]

    for j, t in enumerate(grp):
        t_to_merged[t.t_index] = (int(rep_starts[j]), int(rep_ends[j]))

merged = {}
for t in multi:
    ms, me = t_to_merged[t.t_index]
    key = (t.ref, int(t.strand), ms, me)
    if key not in merged:
        merged[key] = t
print(f"Merged spans (tol=0): {len(merged)}")

# Annotated equivalents
se_by_loc = collections.defaultdict(list)
for t in transcripts:
    if len(t.exons) == 1:
        se_by_loc[(t.ref, int(t.strand))].append((t.start, t.end, t))
for key in se_by_loc:
    se_by_loc[key].sort(key=lambda x: (x[0], x[1]))

covered = set()
n_covered_breakdown = collections.Counter()
for span_key in merged:
    ref, strand, m_start, m_end = span_key
    candidates = se_by_loc.get((ref, strand), [])
    for s_start, s_end, s_tx in candidates:
        if s_start > m_start:
            break
        if s_start <= m_start and m_end <= s_end:
            covered.add(span_key)
            n_covered_breakdown["exact" if (s_start == m_start and s_end == m_end) else "contains"] += 1
            break

uncov_count = len(merged) - len(covered)
print(f"Covered by annotated equiv: {len(covered)}")
print(f"  Coverage breakdown: {dict(n_covered_breakdown)}")
print(f"Uncovered (need synthetic): {uncov_count}")

# Now compare to old nRNA table
import pandas as pd
old_nrna = pd.read_feather("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune/rigel_index/nrna.feather")
old_spans = set()
for _, row in old_nrna.iterrows():
    old_spans.add((row["ref"], int(row["strand"]), int(row["start"]), int(row["end"])))
print(f"\nOld nRNA: {len(old_spans)} unique spans")

# The old system creates an nRNA for every multi-exon transcript with
# t_span > t_exonic (i.e., transcript has intronic space).
# With tol=0, merged_spans should match old_spans.
merged_span_set = set(merged.keys())
in_old_not_merged = old_spans - merged_span_set
in_merged_not_old = merged_span_set - old_spans

print(f"\nIn old but NOT in merged_spans: {len(in_old_not_merged)}")
print(f"In merged_spans but NOT old: {len(in_merged_not_old)}")
print(f"In both: {len(old_spans & merged_span_set)}")

if in_old_not_merged:
    sample = sorted(in_old_not_merged)[:5]
    print(f"\nSample from old but not merged:")
    for s in sample:
        # Find this span in old
        ref, strand, start, end = s
        # Is there a multi-exon transcript matching this?
        matching = [t for t in multi
                    if t.ref == ref and int(t.strand) == strand
                    and t.start == start and t.end == end]
        print(f"  {s} - multi-exon matches: {len(matching)}", end="")
        matching_any = [t for t in transcripts
                        if t.ref == ref and int(t.strand) == strand
                        and t.start == start and t.end == end]
        print(f" - any transcript matches: {len(matching_any)}")
        if matching_any:
            t = matching_any[0]
            print(f"    {t.t_id}: {len(t.exons)} exons, start={t.start}, end={t.end}")
