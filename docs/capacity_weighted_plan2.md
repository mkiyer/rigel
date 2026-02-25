Plan: Capacity-Weighted EM with Per-Fragment Effective Lengths
The EM currently ignores where a fragment falls on a transcript, causing systematic bias when isoforms share a large exon but differ at a small terminal exon. We fix this in two complementary ways: (1) per-fragment effective lengths in the E-step (Bayesian-correct), and (2) a capacity-weighted warm start that distributes ALL fragments (unique and ambiguous) according to geometric difficulty, eliminating reliance on unique reads. The trapezoid coverage model defines "difficulty": fragments near transcript edges have low coverage capacity, so observing them is stronger evidence for that transcript. This resolves your T1/T2/T3 example even when zero unique reads exist.

Steps

Add t_length_arr to ScoringContext (scoring.py:96-134)

New field: t_length_arr: np.ndarray (int32, transcript exonic lengths from index.t_df["length"])
Set in ScoringContext.from_models() alongside existing t_strand_arr
Per-fragment effective length in _score_wta_mrna() (scan.py:109-170)

After computing log_lik = log_strand + log_fl + oh * penalty + log_nm, subtract math.log(max(1, L_t - frag_length + 1)) where L_t = ctx.t_length_arr[t_idx_int] and frag_length = frag_lengths[k]
This replaces the global effective length with a per-fragment version
Neutralize global mRNA effective lengths (locus.py:301-302)

Set eff_len[:n_t] = 1.0 (instead of counter.effective_lengths[t_arr]) since length normalization is now handled per-fragment in each log-likelihood
nRNA and gDNA effective lengths remain unchanged
Cache per-transcript exon intervals on HulkIndex (index.py:429-480)

During HulkIndex.load(), load exon intervals from intervals.feather, group by t_index, store as self._t_exon_intervals: dict[int, np.ndarray] (Nx2 int32 array of sorted genomic [start, end] pairs per transcript)
Mirror the pattern used by exon_lengths_per_gene() at index.py:624 but group by transcript instead of gene
Add accessor get_exon_intervals(t_idx: int) -> np.ndarray
Store genomic_start in the fragment buffer (buffer.py:54-230)

Add genomic_start: int = -1 to BufferedFragment (one int32 per fragment)
Add genomic_start: list[int] to _AccumulatorChunk.__slots__ and its append() method
Add genomic_start: np.ndarray (int32[N]) to _FinalizedChunk
Extract from ResolvedFragment — set during pipeline.py:220-253 as min(exon.start for exon in frag.exons) (the leftmost genomic position of the fragment)
Combined with existing genomic_footprint, this gives both endpoints: genomic_end = genomic_start + genomic_footprint
Implement genomic_to_transcript_pos() (new utility in resolution.py)

Maps a genomic position to a transcript-relative position (0 = transcript 5' end) given sorted exon intervals and strand
For + strand: walk exons left-to-right, accumulate exonic bases up to the position
For - strand: compute the genomic-order position then mirror (tx_length - pos - 1)
~25 LOC, no external dependencies beyond the exon array
Implement compute_fragment_interval_weight() (new function, likely in scoring.py alongside ScoringContext)

Port the user's draft function: takes (frag_start, frag_end, transcript_length, mean_insert) → float
Computes the trapezoid with w_max = min(mean_insert, transcript_length / 2), intersection areas with left ramp / plateau / right ramp, returns w_max / mean_capacity
~30 LOC
Compute capacity weights during scan (scan.py:109-170)

In _score_wta_mrna(), after scoring: for each WTA-winning candidate, compute transcript-relative frag_start via genomic_to_transcript_pos(bf.genomic_start, index.get_exon_intervals(t_idx), strand), then frag_end = frag_start + frag_length
Call compute_fragment_interval_weight(frag_start, frag_end, L_t, mean_insert) to get the weight
Return the weight alongside (overhang, log_lik, count_col) — extend the return dict value to include a 4th element
The ScoringContext will need mean_insert: float added (from finalized FragmentLengthModel.mean)
Thread capacity weights through ScanData and LocusEMInput (scan.py:597-635, estimator.py:264-335)

Add capacity_weights: np.ndarray (float64) to ScanData — parallel to t_indices and log_liks
Add capacity_weights: np.ndarray to LocusEMInput — parallel to t_indices and log_liks
During build_locus_em_data(), extract the locus slice of capacity weights alongside existing arrays
Capacity-weighted initialization in run_locus_em() (estimator.py:625-658)

Replace the current uniform warm start:
With a capacity-weighted distribution: iterate over all units in the locus CSR, read capacity_weights[start:end] for each unit's candidates, mask by prior > 0, normalize weights to sum to 1, accumulate into theta_init[candidates]
Units where all weights are equal (deep plateau) degenerate to the current 1/k behavior — so this is strictly backward-compatible
This gives transcripts with edge-proximal fragments a head start proportional to geometric difficulty


Optional: Capacity-weighted persistent prior (locus.py:322-335)

Controlled by a new capacity_prior_gamma: float parameter (default 0.0 = disabled, suggested 0.01–0.10)
During build_locus_em_data(), after computing capacity weights, accumulate a per-component capacity prior: for each unit, distribute 1 pseudocount proportional to capacity weights across candidates, sum across all units, scale by gamma
Add to existing prior[:n_t], so it persists across EM iterations
This prevents the EM from "forgetting" the geometric evidence as iterations progress
Add regression test (new file: tests/scenarios/test_short_exon_isoform.py)

Scenario A (the pathological case): T1 exons [(1000,1050), (2000,22000)], T2 exons [(1200,1500), (2000,22000)], abundance 10:1. Assert T2 relative error < 20% (currently 53%)
Scenario B (no unique reads): Your T1/T2/T3 example — three overlapping isoforms where no fragment maps uniquely. Assert reasonable accuracy for T2 and T3 (which have short terminal edges)
Scenario C (single exon control): Verify no regression on simple cases
Cleanup: Remove _diag_short_exon.py (the temporary diagnostic script)

Verification

Run full test suite: PYTHONPATH=src conda run -n hulkrna python -m pytest tests/ -x -q — all 583+ tests pass
Run the diagnostic scenario manually to confirm T2 error drops from 53% to <20%
Run the no-unique-reads T1/T2/T3 scenario to verify capacity weighting resolves ambiguity without unique anchors
Benchmark on the existing 10-region benchmark to verify no aggregate regression
Decisions

Capacity weight as initialization, not E-step penalty: The reviewer correctly noted that adding a position penalty to log-likelihood pushes in the wrong direction (penalizes the transcript you're trying to rescue). The weight is used for warm start and optional prior, keeping the E-step generatively pure.
Per-fragment effective length is separate and complementary: Addresses the mathematical correctness of the denominator. Impact is small for the pathological case (~1.3%) but correct in general.
Fragment genomic_start stored per-fragment (not per-candidate): One int32 per fragment is cheap. Transcript-relative mapping is per-candidate at scan time, using cached exon intervals.
No unique-read dependency: The capacity weight is computed for EVERY fragment-candidate pair, so the approach works for complex genes with 15+ isoforms and zero uniquely-mapping fragments.
Persistent prior is optional (gamma=0 by default): Can be tuned empirically. The initialization alone may be sufficient if the EM converges near its starting point in flat-likelihood regions.