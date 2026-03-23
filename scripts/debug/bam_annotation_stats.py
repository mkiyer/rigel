"""Quick BAM annotation stats and GAPDH-specific sanity checks."""
import pysam
import sys
from collections import Counter

bam_path = (
    "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/"
    "gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"
)

# Sample first 1M records for overall stats
zn_counter = Counter()
zp_counter = Counter()
n_records = 0
n_gapdh_truth = 0
n_gapdh_resolved = 0

gapdh_tx = "ENST00000229239.10"

print("Scanning annotated BAM (first 2M records for stats)...")
with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
    for read in bam:
        n_records += 1
        tags = dict(read.get_tags())
        
        # Only count R1 primary (to avoid double-counting)
        if not (read.flag & 0x40) and (read.flag & 0x1):  # PE but not R1
            continue
        if read.flag & 0x100:  # secondary
            continue
            
        zn = tags.get("ZN", -1)
        zp = tags.get("ZP", "?")
        zn_counter[zn] += 1
        zp_counter[zp] += 1
        
        # Check GAPDH truth reads  
        if gapdh_tx in read.query_name:
            n_gapdh_truth += 1
            if zn > 0:
                n_gapdh_resolved += 1
        
        if n_records >= 4_000_000:  # ~2M PE fragments
            break

print(f"\nRecords scanned: {n_records:,}")
print(f"\nZN distribution (R1 primary only):")
total_r1_pri = sum(zn_counter.values())
print(f"  Total R1 primary records: {total_r1_pri:,}")
for zn_val in sorted(zn_counter.keys()):
    pct = 100 * zn_counter[zn_val] / total_r1_pri
    print(f"  ZN={zn_val}: {zn_counter[zn_val]:,} ({pct:.1f}%)")

print(f"\nZP distribution (R1 primary only):")
for zp_val, cnt in zp_counter.most_common():
    pct = 100 * cnt / total_r1_pri
    print(f"  ZP={zp_val}: {cnt:,} ({pct:.1f}%)")

print(f"\nGAPDH truth reads (R1 primary): {n_gapdh_truth}")
print(f"GAPDH resolved (ZN>0): {n_gapdh_resolved}")
