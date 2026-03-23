"""Check all BAM records for a specific GAPDH fragment."""
import pysam

bam_path = (
    "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/"
    "gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"
)
target_qname = "ENST00000229239.10:0-200:f:913"

with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
    for read in bam:
        if read.query_name == target_qname:
            tags = dict(read.get_tags())
            zn = tags.get("ZN", "N/A")
            zt = tags.get("ZT", "N/A")
            zp = tags.get("ZP", "N/A")
            zw = tags.get("ZW", "N/A")
            zs = tags.get("ZS", "N/A")
            zl = tags.get("ZL", "N/A")
            flag = read.flag
            chrom = read.reference_name
            pos = read.reference_start
            cigar = read.cigarstring
            mq = read.mapping_quality
            sec = "sec" if flag & 0x100 else "pri"
            r = "R1" if flag & 0x40 else ("R2" if flag & 0x80 else "SE")
            print(
                f"  {r} {sec} {chrom}:{pos} MQ={mq} CIGAR={cigar} "
                f"ZT={zt} ZP={zp} ZN={zn} ZW={zw} ZS={zs} ZL={zl}"
            )

print("\n--- Checking a few more GAPDH qnames ---")
count = 0
seen_qnames = set()
with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
    for read in bam:
        qname = read.query_name
        if "ENST00000229239.10:" in qname and qname not in seen_qnames:
            seen_qnames.add(qname)
            tags = dict(read.get_tags())
            zn = tags.get("ZN", "N/A")
            zt = tags.get("ZT", "N/A")
            zp = tags.get("ZP", "N/A")
            chrom = read.reference_name
            pos = read.reference_start
            sec = "sec" if read.flag & 0x100 else "pri"
            r = "R1" if read.flag & 0x40 else ("R2" if read.flag & 0x80 else "SE")
            print(f"  {qname}: {r} {sec} {chrom}:{pos} ZT={zt} ZP={zp} ZN={zn}")
            count += 1
            if count >= 20:
                break
