#!/usr/bin/env python3
"""Check if FL values differ between transcript-space projection and
the simple block-based approach for multi-block fragments.

Specifically: for a PE fragment where gstart and gend fall within exons,
the transcript-space FL should equal the spliced distance. But if gstart 
or gend falls in an INTRON (e.g., soft-clipped R1 starts before an exon
boundary), the overhang is added.

Question: Does this overhang addition change FL in a way that affects 
the EM differently?

Let me also consider: for a fragment mapped to a PSEUDOGENE where the 
alignment has two blocks (inner fragment gap from PE):
- ExonBlocks: [(start1,end1), (start2,end2)]
- Multi-block path: gstart=start1, gend=end2  
- Pseudogene is single exon: genomic_to_tx_pos gives gend-gstart

But wait: the inner fragment gap between the two read pairs (the insert) 
is NOT an intron in the pseudogene case. It's just unsequenced bases in 
the middle of the pseudogene's single exon. The FL = gend - gstart = 
full insert size INCLUDING the inner gap. This is correct.

For the REAL multi-exon gene:
- If the inner PE gap aligns exactly with an intron: FL = spliced distance (correct)
- If the inner PE gap is just within one exon: FL = insert size (correct)
- If R1 ends in exon1, R2 starts in exon2, with an intron between:
  gstart = R1_start, gend = R2_end
  tx_pos(gstart) → position in exon1
  tx_pos(gend) → position in exon2
  FL = tx_pos(gend) - tx_pos(gstart) = spliced distance (correct)

Everything seems correct. Let me look at this from the EM perspective.

The EM for a locus with transcript A (real gene) and B (pseudogene):
- pi_A, pi_B = current abundance estimates
- For fragment f that maps to both A and B:
  - lik_A = P(f|A) = strand_A * fl_A * overhang_A * nm_A
  - lik_B = P(f|B) = strand_B * fl_B * overhang_B * nm_B
  - responsibility_A = pi_A * lik_A / (pi_A * lik_A + pi_B * lik_B)

If FL is the same for both, and strand/overhang/NM are the same, then
responsibility = pi/(pi_A + pi_B). The EM converges to abundance ratio
proportional to prior + initial counts.

So the pseudogene leakage is fundamentally about:
1. How many fragments map to both real gene and pseudogene
2. What the initial abundance estimates are
3. What the prior pseudocount is
4. Whether any scoring feature discriminates (splice type!)

The KEY discriminating feature should be SPLICE TYPE. For spliced fragments,
the real gene is strongly favored. For unspliced fragments, there's no 
discrimination.

But wait -- does the scoring code actually use splice type to discriminate?
Let me check if spliced fragments that map to the real gene AND the pseudogene
get different likelihoods.
"""

print("Need to investigate: how does splice type affect mRNA likelihood scoring")
print("for fragments that map to both a real gene and its pseudogene.")
print()
print("If splice_type = SPLICE_ANNOT for the real gene candidate but")
print("splice_type = UNSPLICED for the pseudogene candidate, does the")
print("scoring code assign different mRNA log-likelihoods?")
print()
print("If YES: splice type provides discrimination → less leakage")
print("If NO: splice type doesn't discriminate → leakage depends on EM prior")
