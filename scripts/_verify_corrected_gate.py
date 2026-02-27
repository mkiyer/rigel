"""Verify: single-exon blocked, multi-exon nRNA candidates allowed."""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hulkrna.index import HulkIndex
from hulkrna.config import PipelineConfig
from hulkrna.pipeline import run_pipeline

INDEX_DIR = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/"
    "PVT1_MYC/hulkrna_index"
)
BAM_PATH = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/"
    "PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/"
    "reads_namesort.bam"
)

pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=PipelineConfig())
est = pipe.estimator

# 1) LINC02912 (single-exon, t_idx=290) — must have zero nRNA
print(f"nrna_em_counts[290] (single-exon LINC02912) = {est.nrna_em_counts[290]:.1f}")
print(f"t_counts[290] sum = {est.t_counts[290].sum():.1f}")

# 2) Total nRNA
print(f"Total nrna_em_counts = {est.nrna_em_counts.sum():.1f}")

# 3) Multi vs single exon counts
idx = HulkIndex.load(INDEX_DIR)
spans = idx.t_df["end"].values - idx.t_df["start"].values
lengths = idx.t_df["length"].values
multi_exon = spans > lengths
n_multi = int(multi_exon.sum())
n_single = int((~multi_exon).sum())
print(f"Multi-exon: {n_multi}, Single-exon: {n_single}")

# 4) Top multi-exon nRNA
multi_idx = np.where(multi_exon)[0]
nrna_multi = est.nrna_em_counts[multi_idx]
top = np.argsort(nrna_multi)[::-1][:5]
for i, k in enumerate(top):
    gt = multi_idx[k]
    print(
        f"  Multi-exon top nRNA [{i}]: t_idx={gt}, "
        f"t_id={idx.t_df.loc[gt, 't_id']}, "
        f"nrna_em={est.nrna_em_counts[gt]:.1f}, "
        f"t_counts={est.t_counts[gt].sum():.1f}"
    )

# 5) Check no single-exon transcripts have nRNA
single_idx = np.where(~multi_exon)[0]
nrna_single = est.nrna_em_counts[single_idx]
if nrna_single.sum() > 0:
    print(f"WARNING: single-exon nRNA total = {nrna_single.sum():.2f}")
    bad = np.where(nrna_single > 0)[0]
    for b in bad[:5]:
        gt = single_idx[b]
        print(
            f"  Leak: t_idx={gt}, t_id={idx.t_df.loc[gt, 't_id']}, "
            f"nrna_em={est.nrna_em_counts[gt]:.2f}"
        )
else:
    print("OK: no single-exon transcripts have nRNA counts")
