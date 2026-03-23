#!/usr/bin/env python3
"""Investigate how FL distribution training differs between v3 and v4.

Hypothesis: The FL model trained from unique mappers changed because the 
transcript-space FL fix corrected inflated FL values. This changes the 
probability model, which shifts how fragments are scored.

Key question: For UNSPLICED fragments (both reads within one exon) that map 
to both a real multi-exon gene AND its intronless pseudogene, the FL value 
is the same. But what about the SJ routing? Let's check.

Actually the deeper question: WHAT changed between v3 and v4 at the scoring
level that could cause this regression? The ONLY code change was in 
compute_frag_lengths(). So the FL values assigned to fragments changed.

For MULTI-EXON fragments (spanning a splice junction):
- v3: FL was sometimes over-estimated (the SJ gap bug inflated FL)
- v4: FL is now correct (transcript-space projection)

For SINGLE-BLOCK fragments (reads within one exon):
- Both v3 and v4: FL = block_end - block_start (unchanged)

So the regression must come from multi-exon fragments whose FL changed.
But multi-exon fragments that match a splice junction are routed differently
in the EM (they're more informative). If the FL for these fragments changed,
their EM contribution changes, which redistributes abundance estimates,
which then changes the EM solution for ALL fragments including the unspliced 
ones that can go to pseudogenes.

This is an INDIRECT effect: better FL for spliced fragments → different EM 
initialization/convergence → different steady state that happens to leak 
more to pseudogenes.

Let me verify: did the total counts assigned to truth>0 transcripts decrease 
in v4 vs v3?
"""

import pandas as pd
import numpy as np
from pathlib import Path

V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
V4 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/gdna_none_ss_0.95_nrna_none")


def main():
    v3_mm2 = pd.read_csv(V3 / "per_transcript_counts_minimap2.csv")
    v4_mm2 = pd.read_csv(V4 / "per_transcript_counts_minimap2.csv")
    v3_ora = pd.read_csv(V3 / "per_transcript_counts_oracle.csv")
    v4_ora = pd.read_csv(V4 / "per_transcript_counts_oracle.csv")

    # Merge
    m = v3_mm2[["transcript_id", "gene_id", "gene_name", "mrna_truth", "rigel_minimap2"]].merge(
        v4_mm2[["transcript_id", "rigel_minimap2"]],
        on="transcript_id", suffixes=("_v3", "_v4"),
    ).merge(
        v3_ora[["transcript_id", "rigel_oracle"]].rename(columns={"rigel_oracle": "oracle_v3"}),
        on="transcript_id",
    ).merge(
        v4_ora[["transcript_id", "rigel_oracle"]].rename(columns={"rigel_oracle": "oracle_v4"}),
        on="transcript_id",
    )

    # Check: are the gains in false-positive transcripts correlated between oracle and minimap2?
    print("="*80)
    print("ORACLE vs MINIMAP2 CORRELATION FOR FALSE POSITIVE CHANGES")
    print("="*80)

    low_truth = m["mrna_truth"] <= 5
    
    m["ora_delta"] = m["oracle_v4"] - m["oracle_v3"]
    m["mm2_delta"] = m["rigel_minimap2_v4"] - m["rigel_minimap2_v3"]

    from scipy.stats import spearmanr
    lt = m[low_truth]
    r, p = spearmanr(lt["ora_delta"], lt["mm2_delta"])
    print(f"  Spearman(oracle_delta, mm2_delta) for truth<=5: r={r:.4f}, p={p:.2e}")
    print(f"  Oracle total FP delta: {lt['ora_delta'].sum():.0f}")
    print(f"  MM2 total FP delta:    {lt['mm2_delta'].sum():.0f}")

    # Check: transcripts that gained falsely in mm2 but not oracle
    mm2_only_gain = lt[(lt["mm2_delta"] > 100) & (lt["ora_delta"].abs() < 10)]
    print(f"\n  Transcripts gaining >100 in mm2 but <10 in oracle (truth<=5): {len(mm2_only_gain)}")
    print(f"  Total mm2 gain in these: {mm2_only_gain['mm2_delta'].sum():.0f}")
    for _, r in mm2_only_gain.sort_values("mm2_delta", ascending=False).head(20).iterrows():
        print(f"    {r.transcript_id:>30s} ({str(r.gene_name):>15s})  truth={r.mrna_truth:>5.0f}  ora_v3={r.oracle_v3:>6.0f}  ora_v4={r.oracle_v4:>6.0f}  mm2_v3={r.rigel_minimap2_v3:>6.0f}  mm2_v4={r.rigel_minimap2_v4:>6.0f}")

    # What about the nRNA/gDNA pool-level counts?
    print("\n\n" + "="*80)
    print("POOL-LEVEL COMPARISONS")
    print("="*80)
    
    # From the summary.md:
    print("  Oracle v3:  mRNA pred=9994535  nRNA pred=3680  gDNA pred=1785")
    print("  Oracle v4:  mRNA pred=9994550  nRNA pred=3668  gDNA pred=1782")
    print("  MM2 v3:     mRNA pred=9922277  nRNA pred=72958  gDNA pred=3660")
    print("  MM2 v4:     mRNA pred=9904339  nRNA pred=90543  gDNA pred=4013")
    print()
    print("  Note: MM2 mRNA decreased by ~18k, nRNA increased by ~17.6k")
    print("  This is consistent with: more fragments going to nRNA pool in v4")
    print("  BUT: the nRNA truth is 0, so all nRNA prediction is false positive")
    
    # Check: Does the nRNA increase explain the loss or is it separate?
    print(f"\n  MM2 total mRNA pred v3: {v3_mm2['rigel_minimap2'].sum():.0f}")
    print(f"  MM2 total mRNA pred v4: {v4_mm2['rigel_minimap2'].sum():.0f}")
    print(f"  Delta mRNA pred: {v4_mm2['rigel_minimap2'].sum() - v3_mm2['rigel_minimap2'].sum():.0f}")

    # Key question: where did the ~18k lost mRNA counts go?
    # They went to:
    # 1. nRNA: +17,585
    # 2. gDNA: +353
    # 3. False positive transcripts with zero truth
    
    # Let's verify:
    print("\n  Accounting for lost mRNA:")
    print(f"    nRNA increase: {90543 - 72958:+d}")
    print(f"    gDNA increase: {4013 - 3660:+d}")
    active = m["mrna_truth"] > 0
    v3_active_total = m.loc[active, "rigel_minimap2_v3"].sum()
    v4_active_total = m.loc[active, "rigel_minimap2_v4"].sum()
    print(f"    Active tx decrease: {v4_active_total - v3_active_total:.0f}")
    v3_fp_total = m.loc[~active, "rigel_minimap2_v3"].sum()
    v4_fp_total = m.loc[~active, "rigel_minimap2_v4"].sum()
    print(f"    FP tx increase: {v4_fp_total - v3_fp_total:+.0f}")
    
    # Wait - the total fragments should sum to 10M (or close). Let me check.
    print(f"\n  Total fragments accounted for:")
    print(f"    v3: mRNA={v3_mm2['rigel_minimap2'].sum():.0f} + nRNA=72958 + gDNA=3660 = {v3_mm2['rigel_minimap2'].sum() + 72958 + 3660:.0f}")
    print(f"    v4: mRNA={v4_mm2['rigel_minimap2'].sum():.0f} + nRNA=90543 + gDNA=4013 = {v4_mm2['rigel_minimap2'].sum() + 90543 + 4013:.0f}")

    # The nRNA increase is suspicious - nRNA truth is zero, so ALL nRNA prediction is error
    # And we see the gene-level mm2 MAE got WORSE, which means the FL change is somehow
    # making the EM redistribute more fragments to nRNA/gDNA/pseudogenes

    print("\n\n" + "="*80) 
    print("ANALYSIS: WHICH GENE FAMILIES ARE LOSING TO PSEUDOGENES?")
    print("="*80)

    # For the top gene-level losers, check if there's a pseudogene gaining
    # Need to look at the nearby genes on the same chromosome
    
    # Actually, let's look at the locus structure: when minimap2 creates
    # multimappers, it links the real gene to its pseudogene. The EM then
    # has to choose between them. If the FL model fits both equally well,
    # the EM will split proportionally. 
    
    # The question is: did the FL model training change enough to shift the EM?
    
    # Under v3 (SJ gap bug): some unique mappers to multi-exon transcripts
    # got inflated FL values, which means the FL model had a broader/shifted 
    # distribution. This broader model gave slightly different probabilities
    # for the actual FL values, which affects the EM.
    
    # Under v4 (correct FL): the FL model is narrow and centered correctly.
    # Now all fragments get correct FL values, so the FL probability is 
    # higher for the peak and lower for the tails.
    
    # For pseudogene multimapping fragments with FL ~300bp:
    # - v3 FL model (broad): P(FL=300) is moderate
    # - v4 FL model (narrow): P(FL=300) is high
    # But this is the SAME for both real gene and pseudogene copies!
    # So FL probability doesn't discriminate.
    
    # The mechanism must be elsewhere. Let me check if the nRNA increase
    # is the main driver.
    
    print("\n  The minimap2 regression is +1.81 MAE at transcript level")
    print("  The gene-level regression is +2.48 MAE")
    print("  nRNA false positive increased from 72958 to 90543 (+17585)")
    print("  gDNA false positive increased from 3660 to 4013 (+353)")
    print("  FP mRNA (truth=0) increased from 103685 to 137876 (+34191)")
    print()
    print("  Total increase in 'wasted' counts: 17585 + 353 + 34191 = ~52k")
    print("  This represents 0.52% of 10M fragments")
    print("  But the MAE increase is 1.81 (from 9.64 to 11.44), a 19% increase")
    print()
    print("  The MAE increase comes from the top regressed genes:")
    print("  ACTN4(+3157), RPL5(+2089), EIF3C(+1758), EIF3CL(+1751), GAPDH(+1596)")
    print("  These alone account for ~10k of the gap")


if __name__ == "__main__":
    main()
