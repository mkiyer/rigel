"""Which per-region enrichment readout e_r should drive the capture eff-length contraction?

Ground truth: gDNA is uniform genomic, so the oracle gDNA per-region density IS the probe enrichment.
Test whether the robust proxy (total unspliced density — gDNA+nascent+exon-body, all probe-enriched)
recovers it, including in the no-gDNA condition where the proxy must stand alone (same genome/probes
⇒ the enrichment pattern is identical across conditions).
"""
import numpy as np, pysam
from scipy.stats import spearmanr
from rigel.index import TranscriptIndex

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
df = idx.region_df
length = (df["end"] - df["start"]).to_numpy().astype(float)
starts, ids = {}, {}
for ref, g in df.groupby("ref_name"):
    starts[ref] = g["start"].to_numpy(); ids[ref] = g["region_id"].to_numpy()
R = len(df)


def per_region(cond):
    """Return (gdna_density, unspliced_density) per region from the oracle BAM."""
    ug = np.zeros(R); un = np.zeros(R)
    for r in pysam.AlignmentFile(f"{SUITE}/{cond}/sim_oracle.bam", "rb").fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
            continue
        ref = r.reference_name
        if ref not in starts:
            continue
        i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
        if i < 0:
            continue
        rid = int(ids[ref][i])
        spliced = bool(r.cigartuples) and any(op == 3 for op, _ in r.cigartuples)
        if not spliced:
            un[rid] += 1.0
            if r.query_name.startswith("gdna"):
                ug[rid] += 1.0
    L = np.maximum(length, 1.0)
    return ug / L, un / L


gd_g, un_g = per_region("gdna_gdna300_ss_0.99_nrna_rnd_capture_on")  # has gDNA → true enrichment
_, un_n = per_region("gdna_none_ss_0.99_nrna_rnd_capture_on")        # NO gDNA → proxy stands alone

# Restrict to regions with some coverage (the enrichment is only defined where there are reads).
cov = (gd_g > 0) & (un_g > 0) & (un_n > 0)
print(f"regions with coverage: {cov.sum()} / {R}")
print(f"true enrichment = oracle gDNA density (gdna300 cap_on)")
print(f"  Spearman(true_enrichment, unspliced_density SAME cond)   = "
      f"{spearmanr(gd_g[cov], un_g[cov])[0]:+.3f}  (proxy tracks enrichment when gDNA present)")
print(f"  Spearman(true_enrichment, unspliced_density NO-gDNA cond)= "
      f"{spearmanr(gd_g[cov], un_n[cov])[0]:+.3f}  (proxy recovers enrichment when gDNA ABSENT)")
# Dynamic range check: does the proxy capture the on/off-target chasm magnitude?
def chasm(x):
    p = x[x > 0]
    return np.percentile(p, 95) / max(np.percentile(p, 50), 1e-12)
print(f"\non/off-target dynamic range (p95/p50 density):")
print(f"  true gDNA enrichment:        {chasm(gd_g):.1f}x")
print(f"  unspliced proxy (gdna300):   {chasm(un_g):.1f}x")
print(f"  unspliced proxy (no-gDNA):   {chasm(un_n):.1f}x")
