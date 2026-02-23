Plan: Unambiguous Intronic Overlap for nRNA Initialization
Problem Summary
Currently, intron_bp[k] counts bases of a fragment that fall in the INTRON intervals of candidate transcript k. But transcript k's introns can overlap exons of other transcripts. This means a fragment in "T1's middle exon" gets counted as intronic evidence for T2's nRNA — even though it could equally be T1 mRNA. This ambiguous signal inflates _compute_nrna_init, creating phantom nRNA.

Core Concept: unambig_intron_bp
A new per-candidate metric: bases of the fragment that overlap this candidate's introns AND overlap NO exon from ANY transcript (same gene or different gene). Only these bases are "indisputable" evidence of nascent RNA.

Example (Gene G, T1: 3-exon, T2: 2-exon):

Fragment (2000, 2500): intron_bp[T2]=500, but unambig_intron_bp[T2]=0 (entirely within T1's middle exon)
Fragment (1500, 1700): intron_bp[T1]=200, unambig_intron_bp[T1]=200 (purely intronic, no exon overlap)
Fragment (1900, 2100): intron_bp[T2]=200, unambig_intron_bp[T2]=100 (only the 1900–2000 portion is exon-free)
Step 0: Interval subtraction utility (resolution.py)
Add a helper function:

Given a range [[a, b) and a sorted list of (hole_start, hole_end) intervals to remove, returns the remaining base-pair count. This is a simple O(n) sweep.

Step 1: Modify compute_overlap_profile() (resolution.py)
Return type changes from dict[int, (exon_bp, intron_bp)] to dict[int, (exon_bp, intron_bp, unambig_intron_bp)].

New logic per fragment block:

Query cgranges → get all hits (EXON, INTRON, INTERGENIC)
Collect all EXON hits clipped to block boundaries → merge into sorted non-overlapping union (the "global exonic footprint" within this block)
For each INTRON hit: clip to block → subtract the exon union → remaining length = unambig_intron_bp contribution for that transcript
Accumulate all three values per transcript across blocks
Cost: minimal — we already iterate all hits; adding the merge + subtraction is O(k) per block where k = number of overlapping exons (typically 1–3).

Step 2: Update ResolvedFragment (resolution.py)
overlap_bp dict values become 3-tuples. All consumers of overlap_bp that destructure as (exon_bp, intron_bp) need updating — this includes filter_by_overlap() and any other code unpacking the tuple.

Step 3: Thread through the buffer (buffer.py)
Component	Change
_AccumulatorChunk	Add unambig_intron_bp_flat: list[int]
_FinalizedChunk	Add unambig_intron_bp: np.ndarray  # int32[M_t] parallel to t_indices
BufferedFragment	Add unambig_intron_bp: np.ndarray | None = None
append()	Extract 3rd tuple element
_finalize()	Convert to numpy array
__getitem__()	Slice into BufferedFragment
_spill_chunk() / _load_chunk()	Serialize/deserialize the new array
memory_bytes	Include in sum
Step 4: Use for nRNA initialization only (pipeline.py)
In the pre-EM intronic accumulation loop (~lines 1041–1058), change:

And use unambig_intron_bp[_k] (or just the boolean signal) for the weighted accumulation into transcript_intronic_sense / transcript_intronic_antisense.

The EM itself is unchanged — the assignment matrix and likelihood still use the full intron_bp for the nRNA shadow pool. Ambiguous fragments can still compete in the EM; they just don't seed the nRNA prior. This is the distinction between "what evidence initializes the prior" vs "what fragments the EM can assign."

Step 5: Update tests
test_buffer.py, test_resolution.py, test_fragment.py: update for 3-tuple overlap profile
Re-run the nRNA sweep to verify phantom leakage is eliminated
All 673 existing tests must still pass
What this does NOT change (by design)
gDNA handling: untouched. gDNA pool is separate.
EM assignment matrix: still uses full intron_bp. The EM can assign ambiguous fragments to nRNA if the evidence supports it.
Strand model: untouched.
Index on-disk format: no change to Feather files. The exon union is computed on-the-fly during fragment resolution.
Expected outcome
The nRNA sweep column nRNA=0 (pure mRNA) should show 0 phantom nRNA instead of 355, because the only fragments that can seed nRNA init will be in truly intron-only regions, and with gDNA=0 / mRNA-only simulation, no fragments land in purely intronic space.