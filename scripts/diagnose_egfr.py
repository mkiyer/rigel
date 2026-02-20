#!/usr/bin/env python3
"""Diagnose EGFR per-isoform errors across all conditions."""
import csv
import os
import sys

EGFR_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_2026_02_19/seed_101/EGFR"
TOOLS = ["hulkrna", "hulkrna_mm", "salmon", "kallisto"]

# Read the GTF to get transcript structures
GTF = os.path.join(EGFR_DIR, "region.gtf")
tx_exons = {}
with open(GTF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        feat = parts[2]
        if feat != "exon":
            continue
        chrom = parts[0]
        start = int(parts[3]) - 1  # convert to 0-based
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        tid = None
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                tid = attr.split('"')[1]
                break
        if tid:
            tx_exons.setdefault(tid, []).append((start, end))

# Sort exons and compute lengths
print("=" * 100)
print("EGFR TRANSCRIPT STRUCTURES")
print("=" * 100)
print(f"{'Transcript':<25} {'#Exons':>6} {'Length':>8} {'Span':>10} {'Exon ranges (first/last)'}")
print("-" * 100)
for tid in sorted(tx_exons):
    exons = sorted(tx_exons[tid])
    length = sum(e - s for s, e in exons)
    span = exons[-1][1] - exons[0][0]
    first = f"{exons[0][0]}-{exons[0][1]}"
    last = f"{exons[-1][0]}-{exons[-1][1]}" if len(exons) > 1 else ""
    print(f"{tid:<25} {len(exons):>6} {length:>8} {span:>10}   {first} ... {last}")

# Now analyze per-isoform errors across conditions
print("\n" + "=" * 100)
print("PER-ISOFORM ERRORS ACROSS CONDITIONS")
print("=" * 100)

conditions = sorted([
    d for d in os.listdir(EGFR_DIR) if d.startswith("gdna_")
])

# Accumulate per-transcript errors
tx_errors = {}  # tid -> {tool: [errors]}
tx_truths = {}  # tid -> [truths]
tx_abundances = {}  # tid -> [abundances]

for cond in conditions:
    csv_path = os.path.join(EGFR_DIR, cond, "per_transcript_counts.csv")
    if not os.path.exists(csv_path):
        continue
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tid = row["transcript_id"]
            truth = float(row["truth"])
            abundance = float(row["abundance"])
            tx_truths.setdefault(tid, []).append(truth)
            tx_abundances.setdefault(tid, []).append(abundance)
            for tool in TOOLS:
                est = float(row[tool])
                err = est - truth
                tx_errors.setdefault(tid, {}).setdefault(tool, []).append(err)

# Print per-isoform summary
print(f"\n{'Transcript':<25} {'#Ex':>4} {'Len':>6} {'Abund':>10} "
      f"{'Truth':>8} {'hulkrna':>10} {'salmon':>10} {'kallisto':>10} "
      f"{'h_err':>10} {'s_err':>10} {'k_err':>10}")
print("-" * 140)

tx_list = sorted(tx_truths.keys(), key=lambda t: -sum(tx_truths[t]) / len(tx_truths[t]))

for tid in tx_list:
    avg_truth = sum(tx_truths[tid]) / len(tx_truths[tid])
    avg_abund = sum(tx_abundances[tid]) / len(tx_abundances[tid])
    n_exons = len(tx_exons.get(tid, []))
    length = sum(e - s for s, e in tx_exons.get(tid, []))

    avg_h = avg_truth + sum(tx_errors[tid]["hulkrna"]) / len(tx_errors[tid]["hulkrna"])
    avg_s = avg_truth + sum(tx_errors[tid]["salmon"]) / len(tx_errors[tid]["salmon"])
    avg_k = avg_truth + sum(tx_errors[tid]["kallisto"]) / len(tx_errors[tid]["kallisto"])

    h_mae = sum(abs(e) for e in tx_errors[tid]["hulkrna"]) / len(tx_errors[tid]["hulkrna"])
    s_mae = sum(abs(e) for e in tx_errors[tid]["salmon"]) / len(tx_errors[tid]["salmon"])
    k_mae = sum(abs(e) for e in tx_errors[tid]["kallisto"]) / len(tx_errors[tid]["kallisto"])

    h_bias = sum(tx_errors[tid]["hulkrna"]) / len(tx_errors[tid]["hulkrna"])
    s_bias = sum(tx_errors[tid]["salmon"]) / len(tx_errors[tid]["salmon"])
    k_bias = sum(tx_errors[tid]["kallisto"]) / len(tx_errors[tid]["kallisto"])

    print(f"{tid:<25} {n_exons:>4} {length:>6} {avg_abund:>10.1f} "
          f"{avg_truth:>8.0f} {avg_h:>10.1f} {avg_s:>10.1f} {avg_k:>10.1f} "
          f"{h_bias:>+10.1f} {s_bias:>+10.1f} {k_bias:>+10.1f}")

# Detailed breakdown: which pair of isoforms have the biggest confusion?
print("\n" + "=" * 100)
print("OVERLAP ANALYSIS: PAIRWISE EXON SHARING")
print("=" * 100)

# Compute pairwise exon overlap
def exon_set(tid):
    """Return set of (start,end) tuples for transcript."""
    return set((s, e) for s, e in tx_exons.get(tid, []))

def exon_bases(tid):
    """Return set of base positions covered by exons."""
    bases = set()
    for s, e in tx_exons.get(tid, []):
        bases.update(range(s, e))
    return bases

tids = sorted(tx_exons.keys())
print(f"\n{'T1':<25} {'T2':<25} {'T1_len':>8} {'T2_len':>8} {'Shared_bp':>10} {'T1_frac':>8} {'T2_frac':>8}")
print("-" * 100)

for i, t1 in enumerate(tids):
    b1 = exon_bases(t1)
    l1 = len(b1)
    for j, t2 in enumerate(tids):
        if j <= i:
            continue
        b2 = exon_bases(t2)
        l2 = len(b2)
        shared = len(b1 & b2)
        if shared == 0:
            continue
        f1 = shared / l1 if l1 > 0 else 0
        f2 = shared / l2 if l2 > 0 else 0
        print(f"{t1:<25} {t2:<25} {l1:>8} {l2:>8} {shared:>10} {f1:>8.3f} {f2:>8.3f}")

# Per-condition detail for worst isoform(s)
print("\n" + "=" * 100)
print("PER-CONDITION DETAIL FOR LARGEST-ERROR ISOFORMS")
print("=" * 100)

# Find top 3 isoforms by hulkrna MAE
top3 = sorted(tx_list, key=lambda t: -sum(abs(e) for e in tx_errors[t]["hulkrna"]) / len(tx_errors[t]["hulkrna"]))[:3]

for tid in top3:
    print(f"\n--- {tid} (length={sum(e-s for s,e in tx_exons.get(tid,[])):,}, "
          f"exons={len(tx_exons.get(tid,[]))},"
          f" mean_truth={sum(tx_truths[tid])/len(tx_truths[tid]):.0f}) ---")
    print(f"{'Condition':<30} {'Truth':>8} {'hulkrna':>10} {'h_err':>10} {'salmon':>10} {'s_err':>10} {'kallisto':>10} {'k_err':>10}")
    print("-" * 108)
    for ci, cond in enumerate(conditions):
        truth = tx_truths[tid][ci]
        h_err = tx_errors[tid]["hulkrna"][ci]
        s_err = tx_errors[tid]["salmon"][ci]
        k_err = tx_errors[tid]["kallisto"][ci]
        print(f"{cond:<30} {truth:>8.0f} {truth+h_err:>10.1f} {h_err:>+10.1f} "
              f"{truth+s_err:>10.1f} {s_err:>+10.1f} "
              f"{truth+k_err:>10.1f} {k_err:>+10.1f}")

# Also check hulkrna vs hulkrna_mm differences
print("\n" + "=" * 100)
print("HULKRNA UNIQUE vs MULTIMAP MODE (all isoforms, all conditions)")
print("=" * 100)
diff_sum = 0
diff_count = 0
for tid in tx_list:
    h_errs = tx_errors[tid]["hulkrna"]
    hmm_errs = tx_errors[tid]["hulkrna_mm"]
    h_mae = sum(abs(e) for e in h_errs) / len(h_errs)
    hmm_mae = sum(abs(e) for e in hmm_errs) / len(hmm_errs)
    if abs(h_mae - hmm_mae) > 0.1:
        n_exons = len(tx_exons.get(tid, []))
        avg_truth = sum(tx_truths[tid]) / len(tx_truths[tid])
        print(f"  {tid:<25} truth={avg_truth:>8.0f}  hulkrna_MAE={h_mae:>8.1f}  hulkrna_mm_MAE={hmm_mae:>8.1f}  diff={h_mae-hmm_mae:>+8.1f}")
    diff_sum += abs(h_mae - hmm_mae)
    diff_count += 1
if diff_count:
    print(f"\n  Average |unique - multimap| MAE difference: {diff_sum/diff_count:.2f}")
    print(f"  => Multimap mode {'helps' if diff_sum > 0 else 'makes no difference'}")

print("\nDone.")
