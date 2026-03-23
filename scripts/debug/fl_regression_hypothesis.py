#!/usr/bin/env python3
"""Check how minimap2 aligns reads to pseudogenes vs real genes.

For a fragment that maps to both a multi-exon gene and its intronless pseudogene,
minimap2 will produce:
1. A spliced alignment to the real gene (with N/D operation at the intron)
2. A continuous alignment to the pseudogene (no intron, single block)

For the SAME fragment:
- v3 multi-exon alignment: FL computed via gap-SJ correction
- v3 pseudogene alignment: SINGLE BLOCK → FL = block_end - block_start = genomic span
- v4 multi-exon alignment: FL via transcript-space → correct insert size
- v4 pseudogene alignment: SINGLE BLOCK → FL = block_end - block_start = same as v3

So for the pseudogene, FL is IDENTICAL in v3 and v4 (single-block shortcut).
For the real gene multi-exon alignment, FL changed from gap-SJ-corrected to 
transcript-space-projected.

The question is: did v3 compute a DIFFERENT FL for the real gene compared to v4?

If the gap SJ correction worked correctly in v3, FL should be the same.
If the gap SJ correction FAILED (near-miss SJ), FL would have been inflated in v3.

For the INFLATED-FL case in v3:
- Real gene FL = 500 (inflated, should be 300)  
- Pseudogene FL = 300 (correct, single block)
- FL model distribution peaks at ~294
- P(FL=500) << P(FL=300)
- So real gene is PENALIZED relative to pseudogene
- EM pushes counts AWAY from real gene, TOWARD pseudogene
- Net effect: MORE pseudogene leakage in v3 with inflated FL

But the data shows the OPPOSITE: v4 has MORE pseudogene leakage.
This means the gap-SJ correction was NOT inflating FL for these genes.
OR there's another mechanism at play.

Let me check: did the OBSERVED minimap2 alignments for these genes actually 
produce multi-block (spliced) alignments? If not, ALL alignments would be 
single-block, and FL computation is identical in v3 and v4.

WAIT — if ALL alignments to these genes are single-block (e.g., read pairs where 
both reads land within exons), then the FL did NOT change between v3 and v4.
In that case, the regression must be due to FL model training, not FL computation.

The FL model is trained from UNIQUE mappers (UMQ fragments). If the transcript-space 
projection changed the FL values for unique mappers → the trained FL distribution 
is different → ALL fragments get different FL probabilities → EM converges to 
different solution.

THIS IS THE MOST LIKELY MECHANISM.

Let me verify by checking whether the v4 FL model is narrower/shifted vs v3.
"""

import sys
print("This is a hypothesis document — see comments above.")
print()
print("HYPOTHESIS: The minimap2 regression is caused by the FL model training")
print("being different between v3 and v4, NOT by per-fragment FL computation.")
print()
print("The FL model is trained from unique mappers during the BAM scan phase.")
print("In v3, some unique mapper multi-exon fragments got inflated FL values,")  
print("which broadened the trained distribution. In v4, ALL fragments get correct")
print("FL values, producing a narrower distribution centered at ~294bp.")
print()
print("The narrower distribution changes P(FL) for ALL fragments, not just the")
print("multi-exon ones that changed. This global shift changes the EM solution,")
print("and apparently leads to more pseudogene leakage with minimap2.")
print()
print("To verify: compare the FL distributions trained by v3 vs v4 BAM scans.")
print("Or: run v4 with v3's FL model and see if the regression disappears.")
