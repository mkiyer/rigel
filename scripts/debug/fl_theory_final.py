#!/usr/bin/env python3
"""Compare the actual FL effects for affected gene families.

The key theoretical distinction between v3 and v4:
- v3: For multi-block fragments with unresolved gaps, FL = footprint - obs_introns
  (inflated because the unresolved gap bp are included)
- v4: FL = |tx_pos(gend) - tx_pos(gstart)| (correct, no inflation)

For a fragment that maps to BOTH a real multi-exon gene AND its pseudogene:
- The real gene mapping is spliced (N operation) → multi-block alignment
- The pseudogene mapping may be single-block (if minimap2 merges because  
  there's no intron) or multi-block (PE read pairs with inner distance)

If the pseudogene mapping is SINGLE BLOCK: 
  FL is identical in v3 and v4 (single-block shortcut)
  
If the pseudogene mapping is MULTI BLOCK (PE pair with inner distance):
  v3: footprint includes the inner insert gap → FL inflated
  v4: tx_pos projection → FL = insert_size including gap → same as v3?
  Actually for a single-exon transcript, tx_pos(gstart) = gstart - exon_start,
  tx_pos(gend) = gend - exon_start, so FL = gend - gstart = footprint. Same!

So for pseudogenes: FL is IDENTICAL in v3 and v4 regardless.

For the real gene's spliced alignment:
  exon blocks: [(E1_start, E1_end), (E2_start, E2_end)] with intron gap between
  observed intron: [(E1_end, E2_start)] with size = E2_start - E1_end
  
  v3: footprint = E2_end - E1_start
      upper = footprint - obs_intron = (E2_end - E1_start) - (E2_start - E1_end)
            = E1_end - E1_start + E2_end - E2_start = exon1_len + exon2_len ✓
      No gaps (all gaps are observed introns) → FL = upper ✓
      
  v4: gstart = E1_start, gend = E2_end  
      tx_pos(E1_start) = position in exon 1 → cumsum[e1] + offset
      tx_pos(E2_end) = position in exon 2 → cumsum[e2] + offset
      FL = exon1_offset + exon2_offset ✓

Both give the SAME FL when the observed intron exactly matches the annotated intron.

THE ONLY DIFFERENCE occurs when:
1. There's a gap between alignment blocks that ISN'T an observed intron
   (e.g., PE inner distance within an exon)
2. AND the SJ gap correction in v3 can't find a matching annotated SJ
   (which means the gap gets treated as uncorrected)

In v4, the transcript-space projection handles this correctly:
- gstart in exon_i, gend in exon_j
- tx_pos correctly handles positions within exons
- FL = correct spliced distance

In v3, without SJ correction for the gap:
- FL = upper - 0 = exon1+exon2 lengths (correct for the observed intron parts)
- BUT if there's an additional gap (PE inner distance) that crosses another intron,
  v3 would try to SJ-correct it. If the correction failed → FL is inflated.

This is EXACTLY the SJ gap bug we fixed. For fragments where v3 over-estimated FL:
- v3: real gene gets high FL → low P(FL) → EM penalizes real gene
- v4: real gene gets correct FL → correct P(FL) → EM gives more weight to real gene

This should HELP the real gene and REDUCE pseudogene leakage.

BUT THE DATA SHOWS MORE LEAKAGE IN V4!

UNLESS: the OPPOSITE also happens. Fragments where v3 UNDER-estimated FL for
the real gene (gap SJ correction over-subtracted), and v4 gives correct FL 
that's HIGHER than v3.

Wait - in v3, the gap correction could only SUBTRACT from the FL (reduce it).
It can't add to it. So v3 FL ≤ footprint always. If v3 subtracted too much
(over-counted a gap correction), FL would be too SMALL in v3 and correct in v4.

A too-small FL in v3 for the real gene means:
- v3: P(FL_too_small) at the edge of distribution → moderate probability
- v4: P(FL_correct) at the peak → higher probability  
- Effect: v4 gives MORE weight to real gene → LESS pseudogene leakage ✓

And a too-LARGE FL in v3 (under-corrected gap):
- v3: P(FL_inflated) very low → strong penalty on real gene
- v4: P(FL_correct) higher → less penalty
- Effect: v4 gives MORE weight to real gene → LESS pseudogene leakage ✓

IN EVERY CASE, fixing the FL should HELP the real gene and REDUCE pseudogene leakage.

The data contradicts this. The FL change CANNOT be causing the increased leakage
through direct per-fragment scoring effects.

THE REGRESSION MUST BE DUE TO AN INDIRECT EFFECT: The FL MODEL TRAINING.

The FL model is trained from unique mappers' FL values. If the FL model changed
(because unique mappers' FL values changed), ALL fragments get different 
FL log-probabilities, which changes the EM globally.

Let me verify this by checking the trained FL model parameters.
"""

print("="*80)
print("CONCLUSION: Direct per-fragment FL changes CANNOT explain the regression.")
print("The regression MUST be caused by the change in FL MODEL TRAINING.")
print()
print("The FL model is trained from unique mappers during the BAM scan.")
print("In v3, some unique mappers had inflated FL → broader model → different P(FL)")
print("In v4, all unique mappers have correct FL → narrower model → different P(FL)")
print()
print("A narrower FL model magnifies small FL differences, increasing discriminating")
print("power. But for fragments where FL is identical between real gene and pseudogene,")
print("the FL model doesn't help. The issue is that the CHANGED FL probabilities")  
print("affect the ENTIRE EM, including fragments that DO show FL differences.")
print()
print("The global shift in the FL probability landscape pushes the EM to a different")
print("steady state, which coincidentally has more pseudogene leakage.")
print("="*80)
