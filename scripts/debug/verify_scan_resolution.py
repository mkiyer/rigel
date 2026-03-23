"""Run BamScanner directly to verify pass 1 resolution rate."""
import sys
sys.path.insert(0, "src")

from rigel.index import TranscriptIndex
from rigel.native import BamScanner

idx_dir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
bam_path = (
    "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/"
    "gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam"
)

# Load index
print("Loading index...")
tx_index = TranscriptIndex.load(idx_dir)
resolver = tx_index.resolver

# Create scanner with same settings as benchmark
scanner = BamScanner(
    resolver,
    "auto",          # sj_strand_tag
    True,            # skip_duplicates
    True,            # include_multimap
)

# Scan (single worker for simplicity)
print("Scanning BAM (single-threaded)...")
result = scanner.scan(bam_path, n_workers=1, n_decomp_threads=2)

# Print stats
stats = result["stats"]
print("\n=== BamScanner Stats ===")
for k, v in sorted(stats.items()):
    print(f"  {k}: {v}")

# Check buffer size
accum = result.get("accumulator", {})
if isinstance(accum, dict):
    n_buffered = accum.get("size", "N/A")
    print(f"\n  Buffer size: {n_buffered}")
