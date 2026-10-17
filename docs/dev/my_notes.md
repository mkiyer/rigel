Regarding your comments:

# 1) You wrote: "Fixing an alignment is not proving it. When Y overhangs into X's intron, moving Y's bases across the junction (your fix) gives a sensible reading, S. Often there is a second sensible reading, C: X's short anchor is the artifact and the molecule is contiguous. S can be built now; choosing between S and C needs the likelihood."

Our method should admit BOTH "S" (the RNA hypothesis) and "C" (the gDNA hypothesis). Both hypotheses should be allowed to propagate through and ultimately compete in the EM.

# 2) Multimapping - "Today's multimapper arithmetic can't weigh alternatives." -- your explanation of "best hit per transcript" and "averages gDNA over hits" is murky and obscure. I don't understand what is happening. It certainly sounds WRONG to me but I don't understand the details. I agree we need to finish the multimapping design and ensure it makes sense.

# 3) Calibration vetoing RNA readings -- remember that we have a circular problem. We need to run the EM to know how to assign ambiguous fragments, but we need calibration to run the EM. To circumvent the circularity, the accumulator uses an approximation. It uses a two-pass method -- deposit unique unambiguous fragments first and then deposit ambiguous fragments probabilistically using the first-pass weights. This of course is vulnerable to errors but it is also better than excluding ambiguous fragments from calibration entirely.

Remember the goal of calibration is to estimate the gDNA vs RNA content of the library. The EM may be able to tolerate some calibration error. In the EM the fragments get a chance to compete, albeit with a gDNA prior in place.

# 4) agree that clips can hide a junction and it would be helpful to know this

=====

WAVE 1 COMMENTS

A new design idea --- 

The "purest" way to run Rigel would be to align every fragment twice. Once with a spliced aligner (STAR for example) and once with ONLY the genome (minimap2). The two BAM files would be combined and sorted/collated by read name. Then Rigel would directly see all of the gDNA alternative hypotheses measured empirically and directly with a genome alignment tool, along with the spliced aligner results. During the scan stage, Rigel would then need logic to reconcile which fragments keep the RNA hypothesis (spliced or unspliced), and the gDNA hypothesis (unspliced only), or both.

Of course, there will be MANY alignments where the two aligners agree, and in that case, the scan stage would only keep one alignment (no need for both). The scan stage would output two alignment hypotheses only for alignments that conflict.

The cost is upstream -- two runs of an aligner instead of one, and then the cost of scanning and buffering additional fragment hypotheses. The EM would see more competition between RNA and gDNA than it currently sees. There would be more multimapping hits (fragments with gDNA alignment and an RNA alignment that conflict) and more to be reconciled.

To me, this might be the "pure" way to do this. However, this may be overkill and it has its own potential problems.

======

Using minimap2 during index build

A potentially cheaper option is to extract sequences from the GTF REFERENCE itself and align those back to the GENOME ONLY using a genome aligner mode (minimap2 in non-spliced mode).

This can potentially identify some of the splice junctions that are vulnerable to misalignment. It doesn't account for mismatches as well. It needs to somehow generalize to different read lengths.

These are the same challenges that the "alignable" mappability tool has in generating a splice junction blacklist. This design would supersede the alignable tool's splice blacklist and essentially generate it at index build time.

======

If we continue down our current path, the goal is to scrutinize each fragment alignment. Rigel can do several things: 1) "fix" the alignment (correct an intronic overhang, for example), and 2) generate an alternative gDNA hypothesis for the alignment.. for example a unique spliced alignment that has a valid unspliced alternative.

This path requires putting the genome into the index and incorporating genomic lookups into fragment scrutiny.

We have to answer the question "does this sequence have an alternate alignment?" to the best of our ability, without directly aligning it to the genome with a different aligner. I think quite often we can succesfully check alternative hypotheses locally at the exon-intron boundaries and fix these read pairs.

## Your pairing example -- this example shows "outward" facing reads. Illumina paired end fragments face inward, so these are super super weird corner cases that you are considering. Typically outward facing alignments could indicate circular RNA or a form or chimera (tandem duplication, etc). This is tricky.

## Mate fixing -- for unspliced alignments overhanging an intron, there may be an alternative spliced (rna) alignment. The opposite possible as well -- a spliced alignment with a very short anchor may have a plausible unspliced alignment that overhangs the intron. The goal is for Rigel to consider both alignments. Some of this is handled during scoring.
- Overhang bases are penalized
- Mismatches are penalized


## Unpaired mates

All of this work started when we realized that Rigel was corrupting the fragment length distribution by outputting unpaired reads as "fragments". The unpaired reads are missing their mate and so fragment length cannot be computed. Without a fragment length, a fragment cannot be scored (Rigel uses fragment length to differentiate gdna from rna).

So unpaired reads cannot be estimated by Rigel's current EM. There is no "fix" when we cannot compute the fragment length.

ANY situation where we cannot estimate the fragment length means the read pair cannot be used by the EM -- UNLESS we come up with some sort of rescue scoring mechanism where an unpaired read can still be scored by can't use the fragment length likelihood term and doesn't train the fragment length distribution. This is a separate issue.


## Unannotated events - read-through chimeras, any kind of chimeras, unannotated splice junctions, etc

Rigel has not fully made up its mind on these. The strategies are inconsistent. We need to decide on a concrete strategy that applies to unannotated events. This will require some though.

Ideas:
1) completely withhold fragments with unannotated spliced alignments. they are held out and lost. no contribution to transcript or gene abundance

2) ?allow spliced alignments whose individual aligned blocks are compatible with a single transcript and whose 5' and 3' splice sites are annotated. this is the easiest first step to "rescuing" unannotated spliced alignments. If both of the splice sites are annotated, then this is likely just an unannotated isoform of a known transcript. the problem is that the transcript is not part of the reference -- so rigel would need to add new "implied" transcripts to the reference for EM scoring -- at that point, Rigel becomes a transcriptome assembler, not just an abundance estimation tool. So, i think we have to keep our current scope. These reads cannot be rescued currently. this is my current logic -- do you agree? this applies to any unannotated spliced alignments. they are held out.

=====

WAVE 2 COMMENTS

Rule A - I agree that 'S' and 'C' need to compete in the EM. Why is this design parked? Rigel supports multimapping everywhere now (scan phase, calibration phase, EM phase). Correct?

Rule B - A fragment may have multiple alignments and with the new changes COULD have an RNA and gDNA hypothesis. Each possible alignment is accompanied by a mismatch count and each alignment gets a score that penalizes mismatches.

Rule C - what does this pertain to? what is "support that's held"

=====

WAVE 3 COMMENTS

Yes -- we have to fix all the downstream problems AFTER we make the BAM scan phase correctly handle paired end reads and include the RNA and gDNA hypotheses of read pairs with plausible spliced AND unspliced alignments.

====

Decisions:
1) Don't start wave 1 yet because the design is not done. We need a confident design first.
2) My comments above suggest a strategy for read-throughs. What is "the flank rule"? What transcript will the read-throughs be assigned to?
3) The 'S' and 'C' alternative should compete in the EM
4) Currently the sequence evidence is 'overhang' and 'mismatches'. What else? These can be computed for every paired-end alignment possibility.
5) We progress slowly through this. Wave 1, and then reassess.



































