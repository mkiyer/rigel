"""Fragment span / path-length census on a real BAM -- the evidence for the length limit (design 3.2).

Groups a name-collated BAM by query name, builds each group's genomic span on its dominant reference,
subtracts CIGAR-N introns to get L, and reports how concentrated the edge-crossing load is with and
without a limit on L.

    OMP_NUM_THREADS=1 python scripts/design/bam_spans.py <name-collated.bam> [median_node_bp]
"""

import sys
from collections import defaultdict

import numpy as np
import pysam

BAM = sys.argv[1]
NODE = float(sys.argv[2]) if len(sys.argv) > 2 else 151.0

rows = []
cur, recs = None, []


def flush(recs):
    if not recs:
        return
    byref, nintr = defaultdict(list), defaultdict(int)
    supp = False
    for r in recs:
        if r.is_unmapped or r.is_secondary:
            continue
        supp |= r.is_supplementary
        byref[r.reference_id].append((r.reference_start, r.reference_end))
        for op, ln in r.cigartuples or []:
            if op == 3:
                nintr[r.reference_id] += ln
    if not byref:
        return
    rid = max(byref, key=lambda k: len(byref[k]))
    st = min(a for a, _ in byref[rid])
    en = max(b for _, b in byref[rid])
    rows.append((en - st, (en - st) - nintr[rid], supp, len(byref) > 1))


af = pysam.AlignmentFile(BAM, "rb")
for r in af.fetch(until_eof=True):
    if r.query_name != cur:
        flush(recs)
        recs, cur = [], r.query_name
    recs.append(r)
flush(recs)

span = np.array([x[0] for x in rows])
L = np.array([x[1] for x in rows])
supp = np.array([x[2] for x in rows])
multi = np.array([x[3] for x in rows])
print(f"groups with a mapped block: {L.size:,}")
print(
    f"L percentiles: p50={np.percentile(L, 50):.0f} p90={np.percentile(L, 90):.0f} "
    f"p99={np.percentile(L, 99):.0f} p99.9={np.percentile(L, 99.9):.0f} max={L.max():,}"
)
print(
    f"span > 1 Mb: {(span > 1e6).sum():,} groups; of those {100 * supp[span > 1e6].mean():.1f}% "
    f"carry a supplementary record; blocks on >1 reference: {multi.sum():,}"
)
print(f"\n{'limit on L':>12} {'dropped':>10} {'%':>9}   (compare: limit on SPAN)")
for lim in (500, 1000, 2000, 5000):
    print(
        f"{lim:>12,} {(L > lim).sum():>10,} {100 * (L > lim).mean():>8.3f}%   "
        f"span: {100 * (span > lim).mean():.3f}%"
    )
cross = np.maximum(L / NODE, 0)
o = np.argsort(-cross)
print(
    f"\nunbounded: top 1,000 groups carry {100 * cross[o[:1000]].sum() / cross.sum():.1f}% of crossings"
)
for lim in (1000, 2000):
    k = L <= lim
    c = cross[k]
    print(
        f"L<={lim:,}: top 1,000 carry {100 * c[np.argsort(-c)[:1000]].sum() / c.sum():.2f}%  "
        f"| mean crossings/fragment {c.mean():.2f}  (kept {k.sum():,})"
    )
