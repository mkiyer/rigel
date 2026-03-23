#!/usr/bin/env python3
"""Compare the actual per-transcript FL values from v3 and v4 using the 
BAM annotated files. The 'annotated_default.bam' files should contain 
the resolved fragment results including FL values.

Actually, let's approach this differently. Let me check if there are any
transcripts where v4 now assigns FL=0 (missing) where v3 had FL>0 or 
vice versa. This would directly change which candidate transcripts get 
scored and participate in the EM.

If a candidate transcript gets FL=0 in v4 but FL>0 in v3, it would be 
EXCLUDED from FL scoring, which changes its likelihood, which changes 
the EM solution.

Actually, let me just check the annotated BAM to see what data is stored.
"""

import subprocess
import sys

bam_v3 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"
bam_v4 = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none/align_minimap2/annotated_default.bam"

# Check first few reads from each BAM
for label, bam in [("v3", bam_v3), ("v4", bam_v4)]:
    print(f"\n{'='*80}")
    print(f"  {label} BAM header tags and first 3 reads:")
    print(f"{'='*80}")
    result = subprocess.run(
        ["samtools", "view", "-h", bam],
        capture_output=True, text=True, timeout=10
    )
    lines = result.stdout.split("\n")
    # Print @CO lines (comments often describe tags)
    for line in lines:
        if line.startswith("@CO"):
            print(f"  {line}")
    
    # Print first 3 alignment lines
    count = 0
    for line in lines:
        if not line.startswith("@") and line.strip():
            print(f"  {line[:200]}")
            count += 1
            if count >= 3:
                break

if __name__ == "__main__":
    pass
