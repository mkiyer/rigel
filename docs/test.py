

def compute_fragment_interval_weight(
    frag_start: int,
    frag_end: int,
    transcript_length: int
) -> float:
    """
    Computes the inverse probability weight of a fragment by analytically integrating 
    the coverage trapezoid over the fragment's exact genomic span.
    """
    f_len = frag_end - frag_start
    if f_len <= 0:
        return 1.0

    # 1. Define the maximum height of the trapezoid (Handles short transcripts natively)
    w_max = min(f_len, transcript_length / 2.0)
    
    total_area = 0.0
    
    # 2. Area intersection with the Left Ramp: x in [0, w_max]
    s_left = max(0, min(frag_start, w_max))
    e_left = max(0, min(frag_end, w_max))
    if s_left < e_left:
        # Area of a trapezoid: width * average height
        total_area += (e_left - s_left) * (s_left + e_left) / 2.0
        
    # 3. Area intersection with the Middle Flat Region: x in [w_max, L - w_max]
    s_mid = max(w_max, min(frag_start, transcript_length - w_max))
    e_mid = max(w_max, min(frag_end, transcript_length - w_max))
    if s_mid < e_mid:
        # Area of a rectangle: width * height
        total_area += (e_mid - s_mid) * w_max
        
    # 4. Area intersection with the Right Ramp: x in [L - w_max, L]
    s_right = max(transcript_length - w_max, min(frag_start, transcript_length))
    e_right = max(transcript_length - w_max, min(frag_end, transcript_length))
    if s_right < e_right:
        # Height at any point x on the right ramp is (L - x)
        h_start = transcript_length - s_right
        h_end = transcript_length - e_right
        total_area += (e_right - s_right) * (h_start + h_end) / 2.0
        
    # 5. Calculate area per base across the fragment span
    area_per_base = total_area / f_len    
    # Safe floor to prevent division by zero in edge cases
    area_per_base = max(area_per_base, 1.0)
    
    # 6. Return the Inverse Weight
    # Normalized by w_max so perfectly centered fragments yield exactly 1.0
    return w_max / area_per_base




# Initialize a blank prior array
alpha_prior = np.full(num_transcripts, em_pseudocount)

# Loop over all fragments (both unique and ambiguous)
for frag_idx, candidates in enumerate(fragment_candidates):
    
    # Calculate geometric weights for this fragment against all compatible transcripts
    weights = []
    for t_idx in candidates:
        transcript = transcripts[t_idx]
        weight = compute_fragment_interval_weight(frag[frag_idx], transcript)
        weights.append(weight)
        
    weights = np.array(weights)
    
    # Normalize weights so they sum to 1.0 (fractional distribution)
    fractional_shares = weights / weights.sum()
    
    # Add the fractional shares to the prior
    for i, t_idx in enumerate(candidates):
        alpha_prior[t_idx] += fractional_shares[i]



# Initialize the prior array with standard machine-epsilon pseudocounts
# (em_pseudocount is usually your existing 1e-8 to prevent log(0))
alpha_prior = np.full(num_transcripts, em_pseudocount)

# We only care about the N fragments that are actually ambiguous
n_ambiguous = len(fragment_candidates) 

# If there are no fragments, do nothing
if n_ambiguous > 0:
    # 1. Accumulate the raw geometric shares
    raw_shares = np.zeros(num_transcripts)
    
    for frag_idx, candidates in enumerate(fragment_candidates):
        frag_tx_start = get_transcript_coordinate(frag[frag_idx].start, transcript)
        frag_tx_end = get_transcript_coordinate(frag[frag_idx].end, transcript)
        
        weights = []
        for t_idx in candidates:
            transcript = transcripts[t_idx]
            weight = compute_fragment_weight(
                frag_start=frag_tx_start, 
                frag_end=frag_tx_end, 
                transcript_length=transcript.length
            )
            weights.append(weight)
            
        weights = np.array(weights)
        fractional_shares = weights / weights.sum()
        
        # Add to raw shares
        for i, t_idx in enumerate(candidates):
            raw_shares[t_idx] += fractional_shares[i]
            
    # 2. The Parameter-Free Scaling: Divide by N 
    # This guarantees sum(scaled_shares) == 1.0 virtual read
    scaled_shares = raw_shares / n_ambiguous
    
    # 3. Add exactly 1.0 total virtual read of structural opinion to the prior
    alpha_prior += scaled_shares