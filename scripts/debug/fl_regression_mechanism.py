#!/usr/bin/env python3
"""Deep investigation: How does transcript-space FL projection affect
multi-block fragments that map to BOTH a multi-exon transcript AND 
an intronless pseudogene?

Key scenario for understanding the regression:
A paired-end fragment from GAPDH (12 exons) generates a minimap2 alignment 
that has two blocks (R1 and R2 in different exons). This produces:

Alignment 1 (to real GAPDH): 
  Two exon blocks with an intron gap.
  v3: FL = footprint - observed_introns - gap_correction
  v4: FL = |tx_pos(gend) - tx_pos(gstart)| via transcript-space

Alignment 2 (to GAPDH pseudogene, intronless):
  Single continuous block (because the pseudogene has no intron).
  v3: FL = block_end - block_start = read_span (NO splicing)
  v4: FL = block_end - block_start = same (single-block shortcut)

WAIT - but does minimap2 align a spliced read to a pseudogene as a 
single block? YES - because the pseudogene has no intron, minimap2 
would align R1 and R2 as a single continuous alignment (or two blocks 
with a small insert gap, but no N operation). 

Actually, for a PE fragment where R1 maps to positions 1-150 and R2 
maps to positions 250-400 of the pseudogene, minimap2 would report 
this as... hmm, it depends on the insert size. For a concordant pair 
with inferred insert, the BAM scanner computes fragment from the 
two mates' alignment blocks.

Let me think about this differently. The BAM scanner:
1. Groups reads by query name
2. For a pair where R1 maps to positions A-B and R2 maps to C-D:
   - exon blocks = [(A,B), (C,D)] if they're on the same reference
   - If aligned to pseudogene: both blocks are within the single exon
   - If aligned to real gene: blocks may span an intron

For the pseudogene case:
- exon blocks = [(A,B), (C,D)] with A < B < C < D and the gap [B,C) 
  is within the single pseudogene exon (not an intron)
- This is a MULTI-BLOCK fragment (exons.size() == 2)
- v3: footprint = D-A, observed_introns = 0, upper = D-A
  Then checks gaps: gap [B,C) — not in observed introns
  SJ gap query: looks for annotated introns overlapping [B,C)
  For pseudogene: NO annotated introns → gap_correction = 0
  FL = upper - 0 = D-A = genomic span (INCLUDES the inner distance)
  
- v4: gstart = A, gend = D
  For pseudogene (single exon [start, end)):
    tx_pos(A) = A - start
    tx_pos(D) = D - start  
    FL = |tx_pos(D) - tx_pos(A)| = D - A = same as v3!

So for the pseudogene, FL is IDENTICAL in v3 and v4.

For the real multi-exon gene:
- v3: gap [B,C) is queried against SJ index for the real gene
  If an annotated intron MATCHES: gap_correction = overlap → FL correct
  If NO annotated intron matches (off-by-1, shifted, etc.): FL = genomic span (inflated)
  
- v4: transcript-space projection always gives correct FL
  gstart maps to exon in transcript space
  gend maps to another exon in transcript space
  FL = correct intronic-gap-subtracted value

So the ONLY difference is for real-gene multi-block fragments where v3's 
SJ gap query FAILED. If SJ gap query succeeded, v3 = v4. If it failed, 
v3 had inflated FL for the real gene.

With an inflated FL for the real gene:
- P(FL_inflated) << P(FL_correct)
- The real gene gets PENALIZED → EM pushes counts to pseudogene
- This means v3 had MORE leakage, not less!

But the data shows v4 has more leakage. CONTRADICTION.

Unless the SJ gap query bug worked in the OPPOSITE direction for some genes:
it UNDER-counted the gap correction, leading to a SMALLER FL for the real gene.
Is this possible? Let me check the overlap-based SJ matching...

Actually there's another possibility: the v3 SJ gap correction overshot.
If the annotated intron was LONGER than the actual gap (the alignment block 
endpoints are slightly inside the intron), the overlap could be larger than 
the gap size. In v3:
  t_gap_size[ti] += overlap  (where overlap = min(he, ge) - max(hs, gs))
  sz = upper - corr
If the overlap exceeds the gap size, the FL would be UNDER-corrected 
(but each gap is capped at its own size since overlap <= ge-gs).

Hmm, overlap = min(he, ge) - max(hs, gs). If he > ge, then min = ge.
If hs < gs, then max = gs. So overlap = ge - gs = gap_size. 
The overlap can't exceed the gap size. So no over-correction.

The overlap CAN be less than the gap size if the SJ doesn't fully span 
the gap. This happens when:
- The annotated intron starts after the gap start (hs > gs)
- The annotated intron ends before the gap end (he < ge)
In these cases, the correction is partial → FL is slightly inflated.

So in v3, FL was sometimes SLIGHTLY inflated for real genes when the SJ 
didn't perfectly match the gap. In v4, FL is always correct.

But slight inflation PENALIZES the real gene (lower P(FL)), which should 
push counts TOWARD pseudogenes. So fixing the inflation in v4 should 
REDUCE pseudogene leakage, not increase it.

YET THE DATA SHOWS INCREASED LEAKAGE IN V4.

Something else must be going on. Let me reconsider the FL model training 
hypothesis more carefully.
"""

print("See the extensive analysis in the comments above.")
print()
print("CONCLUSION: The per-fragment FL computation change alone cannot explain")
print("why v4 has MORE pseudogene leakage than v3.")
print()
print("The mechanism must involve the FL MODEL TRAINING or an INDIRECT EM effect.")
print("Need to compare the actual FL distributions from v3 and v4 BAM scans.")
