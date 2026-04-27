"""Diagnose SRD Pass-2 mixture EM convergence on a fallback library.

Re-runs the mixture EM on dna10m's calibration inputs with max_iter=2000
and prints the pi trajectory to determine whether (a) it truly converges
just past 200, (b) drifts slowly (need wider tol or more iters), or
(c) oscillates (need different stopping rule).
"""
from __future__ import annotations

import numpy as np

from rigel.calibration._categorize import categorize_chunk, FragmentCategory, build_t_to_ref_id
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer

BAM = (
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/"
    "mctp_vcap_rna20m_dna10m/rigel/annotated.bam"
)
INDEX = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"

cfg = BamScanConfig()
idx = TranscriptIndex.load(INDEX)
_stats, _strand, fl_models, buf = scan_and_buffer(BAM, idx, cfg)
print(f"Buffer: {len(buf):,} fragments")

# Build the same pool histogram SRD uses.
n_bins = fl_models.global_model.max_size + 1
pool_hist = np.zeros(n_bins, dtype=np.float64)
from rigel.splice import SpliceType
spliced_counts = np.asarray(
    fl_models.category_models[SpliceType.SPLICED_ANNOT].counts, dtype=np.float64
)
global_counts = np.asarray(fl_models.global_model.counts, dtype=np.float64)

POOL_CATS = {
    int(FragmentCategory.EXON_INCOMPATIBLE),
    int(FragmentCategory.INTRONIC),
    int(FragmentCategory.INTERGENIC),
}
t_strand = np.asarray(idx.t_to_strand_arr, dtype=np.int8)
t_ref = build_t_to_ref_id(idx.t_df)
for chunk in buf.iter_chunks():
    unique = chunk.num_hits == 1
    if not unique.any():
        continue
    cat = categorize_chunk(
        chunk, t_to_strand=t_strand, t_to_ref_id=t_ref, exon_fit_tolerance_bp=5,
    )
    c = cat.category[unique]
    fl = cat.frag_length[unique]
    mask = np.isin(c, list(POOL_CATS))
    fl = fl[mask]
    fl = fl[(fl > 0) & (fl < n_bins)]
    np.add.at(pool_hist, fl, 1.0)

# Add the C++ intergenic FL counts.
ig = np.asarray(fl_models.intergenic.counts, dtype=np.float64)
m = min(ig.size, n_bins)
pool_hist[:m] += ig[:m]

# RNA prob = spliced FL normalized.
if spliced_counts.size != n_bins:
    sc = np.zeros(n_bins, dtype=np.float64)
    sc[: min(spliced_counts.size, n_bins)] = spliced_counts[: min(spliced_counts.size, n_bins)]
    spliced_counts = sc
rna = spliced_counts / max(spliced_counts.sum(), 1.0)

# Run EM with logging at fine granularity.
import logging
logging.basicConfig(level=logging.INFO)

pi_traj = []
pi = 0.1
gdna = pool_hist + 1.0
gdna /= gdna.sum()
for it in range(1, 2001):
    num = pi * gdna
    den = num + (1.0 - pi) * rna
    with np.errstate(divide="ignore", invalid="ignore"):
        r_g = np.where(den > 0, num / den, 0.0)
    e_g = pool_hist * r_g
    new_pi = float(e_g.sum() / max(pool_hist.sum(), 1e-300))
    new_gdna = e_g + 1.0
    new_gdna /= new_gdna.sum()
    delta = abs(new_pi - pi)
    pi_traj.append((it, new_pi, delta))
    if it in (1, 10, 50, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000) or delta < 1e-4:
        print(f"iter {it:5d}: pi={new_pi:.7f}  |dpi|={delta:.2e}")
    pi = new_pi
    gdna = new_gdna
    if delta < 1e-4:
        print(f"--> converged at iter {it}")
        break

print(f"\nFinal pi={pi:.6f}, n_pool={int(pool_hist.sum()):,}")
