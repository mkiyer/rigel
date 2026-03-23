"""Check buffer frag_id distribution from pass 1 scanner."""
import sys
sys.path.insert(0, "src")

import numpy as np
from rigel.index import TranscriptIndex
from rigel.native import BamScanner

idx_dir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
bam_path = (
    "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/"
    "gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam"
)

print("Loading index...")
tx_index = TranscriptIndex.load(idx_dir)

scanner = BamScanner(tx_index.resolver, "auto", True, True)
print("Scanning...")
result = scanner.scan(bam_path, n_workers=1, n_decomp_threads=2)

buf = result["buffer"] if "buffer" in result else None
print(f"\nResult keys: {list(result.keys())}")
if buf is None:
    print("No 'buffer' key! Looking at what's available...")
    for k, v in result.items():
        if isinstance(v, dict):
            print(f"  {k}: dict with keys {list(v.keys())[:10]}")
        elif isinstance(v, (bytes, bytearray)):
            print(f"  {k}: bytes, len={len(v)}")
        elif hasattr(v, '__len__'):
            print(f"  {k}: {type(v).__name__}, len={len(v)}")
        else:
            print(f"  {k}: {type(v).__name__} = {v}")
    sys.exit(0)

frag_ids = np.frombuffer(buf.get("frag_id", b""), dtype=np.int64)
print(f"Frag ID array length: {len(frag_ids)}")
if len(frag_ids) > 0:
    print(f"  Min frag_id: {frag_ids.min()}")
    print(f"  Max frag_id: {frag_ids.max()}")
    unique_fids = np.unique(frag_ids)
    n_unique = len(unique_fids)
    print(f"  Unique frag_ids: {n_unique}")
    expected_max = int(unique_fids[-1])
    print(f"  Expected range: [0, {expected_max}] = {expected_max + 1} values")
    print(f"  Coverage: {n_unique / (expected_max + 1) * 100:.1f}%")

    # Check for duplicate frag_ids (multimappers have multiple entries)
    fid_counts = np.bincount(frag_ids.astype(np.int64))
    multi = np.sum(fid_counts > 1)
    print(f"  Frag_ids with >1 entry: {multi}")
    print(f"  Frag_ids with 0 entries (gaps): {expected_max + 1 - n_unique}")
