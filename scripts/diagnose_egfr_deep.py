#!/usr/bin/env python3
"""Deep-dive into EGFR overlap scoring: how do fragments discriminate T1 vs T2?

ENST00000275493.7 (T1): 28 exons, 9905bp, span=192612, starts at 0
ENST00000450046.2 (T2): 28 exons, 9700bp, span=101814, starts at 90706

They share 9464bp (95.5% of T1, 97.6% of T2).

For a fragment falling in the shared region, both transcripts will have
high exon_frac. The question is: what discriminates them in the EM?
"""
import os

EGFR_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_2026_02_19/seed_101/EGFR"
GTF = os.path.join(EGFR_DIR, "region.gtf")

# Parse all exons
tx_exons = {}
with open(GTF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9 or parts[2] != "exon":
            continue
        start = int(parts[3]) - 1
        end = int(parts[4])
        attrs = parts[8]
        tid = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                tid = attr.split('"')[1]
                break
        if tid:
            tx_exons.setdefault(tid, []).append((start, end))

# Sort exons
for tid in tx_exons:
    tx_exons[tid].sort()

T1 = "ENST00000275493.7"
T2 = "ENST00000450046.2"

def exon_bases(tid):
    return set()
def exon_intervals(tid):
    return sorted(tx_exons.get(tid, []))

print(f"T1 ({T1}): {len(tx_exons[T1])} exons")
for s, e in tx_exons[T1]:
    print(f"  [{s:>7}, {e:>7}) size={e-s}")

print(f"\nT2 ({T2}): {len(tx_exons[T2])} exons")
for s, e in tx_exons[T2]:
    print(f"  [{s:>7}, {e:>7}) size={e-s}")

# Find exons unique to T1 (not in T2) and vice versa
t1_bases = set()
for s, e in tx_exons[T1]:
    t1_bases.update(range(s, e))

t2_bases = set()
for s, e in tx_exons[T2]:
    t2_bases.update(range(s, e))

shared = t1_bases & t2_bases
t1_only = t1_bases - t2_bases
t2_only = t2_bases - t1_bases

print(f"\nShared bases: {len(shared)}")
print(f"T1-only bases: {len(t1_only)}")
print(f"T2-only bases: {len(t2_only)}")
print(f"T1 effective length: {len(t1_bases)}")
print(f"T2 effective length: {len(t2_bases)}")

# Where are the T1-only bases?
if t1_only:
    t1_only_sorted = sorted(t1_only)
    # Find contiguous runs
    runs = []
    start = t1_only_sorted[0]
    end = t1_only_sorted[0] + 1
    for b in t1_only_sorted[1:]:
        if b == end:
            end = b + 1
        else:
            runs.append((start, end))
            start = b
            end = b + 1
    runs.append((start, end))
    print(f"\nT1-only regions:")
    for s, e in runs:
        print(f"  [{s:>7}, {e:>7}) size={e-s}")

if t2_only:
    t2_only_sorted = sorted(t2_only)
    runs = []
    start = t2_only_sorted[0]
    end = t2_only_sorted[0] + 1
    for b in t2_only_sorted[1:]:
        if b == end:
            end = b + 1
        else:
            runs.append((start, end))
            start = b
            end = b + 1
    runs.append((start, end))
    print(f"\nT2-only regions:")
    for s, e in runs:
        print(f"  [{s:>7}, {e:>7}) size={e-s}")

# Simulate the overlap scoring
print("\n" + "=" * 80)
print("SIMULATED OVERLAP SCORING")
print("=" * 80)

# For a hypothetical 150bp fragment in different positions:
import numpy as np

FRAG_LEN = 150

# 1) Fragment in T1-only region (early exons, before T2 starts)
# T1 starts at position 0; T2 starts at position 90706
# So positions 0~90705 are T1-only exonic territory (minus introns)
t1_exon_set = set()
for s, e in tx_exons[T1]:
    for b in range(s, e):
        t1_exon_set.add(b)

t2_exon_set = set()
for s, e in tx_exons[T2]:
    for b in range(s, e):
        t2_exon_set.add(b)

# Sample positions within T1 exons
t1_exon_bases_sorted = sorted(t1_exon_set)
t2_exon_bases_sorted = sorted(t2_exon_set)

# Get T1 exon blocks
t1_blocks = tx_exons[T1]
t2_blocks = tx_exons[T2]

# Check: fragments that land in the first few exons of T1 (before T2 starts)
print("\nFragment scenarios for T1/T2 discrimination:")
print(f"T1 first exon: {t1_blocks[0]}")
print(f"T2 first exon: {t2_blocks[0]}")
print(f"T1 exon region before T2 starts ({t2_blocks[0][0]}):")
early_t1_exon_bp = sum(min(e, t2_blocks[0][0]) - s for s, e in t1_blocks if s < t2_blocks[0][0])
print(f"  T1 exonic bases before T2 start: {early_t1_exon_bp}")

print(f"\nT1 last exon: {t1_blocks[-1]}")
print(f"T2 last exon: {t2_blocks[-1]}")
late_t1_only_bp = sum(e - max(s, t2_blocks[-1][1]) for s, e in t1_blocks if e > t2_blocks[-1][1])
print(f"  T1 exonic bases after T2 end: {late_t1_only_bp}")

# KEY INSIGHT: How does the overlap scoring handle a fragment in shared region?
# Consider a fragment at position ~100000 (in the shared region)
# Both T1 and T2 have exonic coverage there.
# Fragment exon_bp for T1 and T2 would both be ~150bp (or close).
# exon_frac for both ≈ 1.0
# So the EM sees NO discrimination between T1 and T2 for shared-region fragments.
#
# The ONLY way to discriminate is:
# 1. Fragments in T1-only regions (before 90706 or after 192520)
# 2. Fragments in T2-only regions (if any)
# 3. Splice junctions unique to one isoform
# 4. Insert size differences (different intron lengths)

# Check splice junctions
print("\n" + "=" * 80)
print("SPLICE JUNCTION ANALYSIS")
print("=" * 80)

def get_introns(tid):
    exons = sorted(tx_exons[tid])
    introns = []
    for i in range(1, len(exons)):
        introns.append((exons[i-1][1], exons[i][0]))
    return introns

t1_introns = set((s, e) for s, e in get_introns(T1))
t2_introns = set((s, e) for s, e in get_introns(T2))

shared_introns = t1_introns & t2_introns
t1_only_introns = t1_introns - t2_introns
t2_only_introns = t2_introns - t1_introns

print(f"T1 introns: {len(t1_introns)}")
print(f"T2 introns: {len(t2_introns)}")
print(f"Shared introns: {len(shared_introns)}")
print(f"T1-only introns: {len(t1_only_introns)}")
print(f"T2-only introns: {len(t2_only_introns)}")

if t1_only_introns:
    print("\nT1-only splice junctions:")
    for s, e in sorted(t1_only_introns):
        print(f"  [{s:>7}, {e:>7}) size={e-s}")

if t2_only_introns:
    print("\nT2-only splice junctions:")
    for s, e in sorted(t2_only_introns):
        print(f"  [{s:>7}, {e:>7}) size={e-s}")

# Now analyze: what fraction of simulated fragments would be
# discriminating (unique to T1 or T2 based on SJ or exon-only regions)?
print("\n" + "=" * 80)
print("FRAGMENT DISCRIMINATION POWER (approximate)")
print("=" * 80)

# For a 150bp fragment from T1:
# - If it spans a SJ unique to T1 → definitively T1
# - If it lands entirely in T1-only exonic region → T1 favored by overlap
# - If it lands in shared exonic region (no unique SJ) → ambiguous
# Fraction of T1's length in discriminating regions:
t1_len = len(t1_exon_set)
t2_len = len(t2_exon_set)

# T1-only exonic bases give overlap discrimination
# T1-only SJs give SJ discrimination (but only for spliced fragments)
print(f"T1 total exonic length: {t1_len}")
print(f"  T1-only exonic bases: {len(t1_only)} ({100*len(t1_only)/t1_len:.1f}%)")
print(f"  Shared exonic bases:  {len(shared)} ({100*len(shared)/t1_len:.1f}%)")
print(f"  T1 SJ discriminating: {len(t1_only_introns)}/{len(t1_introns)} introns")

print(f"\nT2 total exonic length: {t2_len}")
print(f"  T2-only exonic bases: {len(t2_only)} ({100*len(t2_only)/t2_len:.1f}%)")
print(f"  Shared exonic bases:  {len(shared)} ({100*len(shared)/t2_len:.1f}%)")
print(f"  T2 SJ discriminating: {len(t2_only_introns)}/{len(t2_introns)} introns")

# KEY: For unspliced fragments in the shared region (e.g. first/last exon):
# overlap fraction for BOTH T1 and T2 ≈ 1.0
# → log(overlap_frac) ≈ 0 for both → ZERO discrimination from overlap
# → EM distributes proportional to theta (abundance) × eff_length
#
# For spliced fragments spanning a shared SJ: same problem
# For spliced fragments spanning a unique SJ: perfectly discriminating
#
# The "leak" is that ~95% of T1's bases are shared with T2,
# and for fragments landing in those shared bases, the overlap
# penalty provides NO discrimination between T1 and T2.

print("\n" + "=" * 80)
print("ROOT CAUSE ANALYSIS")
print("=" * 80)
print(f"""
PROBLEM: T1 ({T1}) is underestimated by ~2031 counts.
         T2 ({T2}) is overestimated by ~1656 counts.
         Net: ~375 counts leak elsewhere.

95.5% of T1's exonic bases are shared with T2.
Only 4.5% of T1's bases ({len(t1_only)} bp) are unique to T1.

For a fragment landing in the shared region:
  - Exon overlap fraction ≈ 1.0 for BOTH T1 and T2
  - log(overlap_frac) ≈ 0 → NO overlap penalty for either
  - The EM distributes based on theta/eff_len alone
  - This means the EM's only discriminator is:
    1. Unique fragments (from 4.5% T1-only region)
    2. SJ-specific fragments (unique introns)
    3. Insert size (different intron lengths between SJ pairs)

For fragments in the T1-only region (first ~441 bp of exon 1):
  - These land in T1 EXON but T2 INTERGENIC
  - overlap_frac for T1 = 1.0 (full exon overlap)
  - overlap_frac for T2 = depends...
  
  CRITICAL QUESTION: Is T2 even in the candidate set for these fragments?
  - If the fragment is at position 50 (T1 exon 1), T2's first exon
    starts at 90706 — so T2 has ZERO exon/intron overlap
  - The fragment resolves to T1 (and maybe other early isoforms)
  - These are the ONLY fragments that definitively count for T1

  But for fragments in the SHARED region (positions 90706+):
  - Both T1 and T2 are candidates with nearly identical overlap
  - The EM must split these fragments between T1 and T2
  - If the EM starts with roughly equal priors, each gets ~50%
  - T1 should get MORE (it has higher abundance), but the EM
    can only distinguish them through the small unique evidence

SOLUTION ANALYSIS:
  The current overlap scoring (exon_frac = exon_bp / frag_length)
  is at the fragment level. If a fragment fully overlaps a transcript
  exon, exon_frac = 1.0 regardless of how much of the TRANSCRIPT
  the fragment covers.

  What's MISSING is a transcript-level compatibility signal: 
  "how well does the fragment's position/structure match what we'd
  expect from this specific transcript?"

  Salmon/kallisto handle this implicitly through their equivalence
  class model — fragments are assigned to equivalence classes of
  transcripts, and the EM distributes based on transcript length
  and abundance. The key difference: salmon/kallisto use k-mer
  based pseudo-alignment which naturally creates more specific
  equivalence classes that capture positional information.
  
  hulkrna's overlap scoring sees all shared-region fragments as
  equally compatible with both T1 and T2, losing the positional
  information that a read coming from position 100000 has a different
  probability of originating from T1 vs T2 based on their relative
  transcript structures.
""")
