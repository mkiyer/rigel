#!/usr/bin/env python3
"""Symbolicate _em_impl addresses from native sample output."""
import re
import subprocess
import sys
from collections import Counter

sample_file = sys.argv[1]
so_path = "/Users/mkiyer/sw/miniforge3/envs/rigel/lib/python3.12/site-packages/rigel/_em_impl.abi3.so"

# Parse all _em_impl addresses (leaf samples = lines with count at start)
# Format: "1757       ???  (in _em_impl.abi3.so)  load address 0x... + 0x...  [0x...]"
pattern = re.compile(
    r'(\d+)\s+\?\?\?\s+\(in _em_impl\.abi3\.so\)\s+load address (0x[0-9a-f]+) \+ (0x[0-9a-f]+)\s+\[0x[0-9a-f]+\]'
)

# Also match inline format: "| 28 ???  (in _em_impl.abi3.so)  load address ..."
pattern2 = re.compile(
    r'\|\s+(\d+)\s+\?\?\?\s+\(in _em_impl\.abi3\.so\)\s+load address (0x[0-9a-f]+) \+ (0x[0-9a-f]+)'
)

offset_counts = Counter()
load_addr = None

with open(sample_file) as f:
    for line in f:
        for pat in [pattern, pattern2]:
            m = pat.search(line)
            if m:
                count = int(m.group(1))
                load_addr = m.group(2)
                offset = m.group(3)
                offset_counts[offset] += count

if not offset_counts:
    print("No _em_impl addresses found")
    sys.exit(1)

print(f"Load address: {load_addr}")
print(f"Total _em_impl leaf samples: {sum(offset_counts.values())}")
print()

# Symbolicate top offsets
top_offsets = offset_counts.most_common(20)
total = sum(offset_counts.values())

print(f"{'Samples':>8}  {'%':>6}  {'Offset':>10}  Function")
print("-" * 80)

for offset, count in top_offsets:
    addr = hex(int(load_addr, 16) + int(offset, 16))
    try:
        result = subprocess.run(
            ["atos", "-o", so_path, "-l", load_addr, addr],
            capture_output=True, text=True, timeout=5
        )
        sym = result.stdout.strip()
    except Exception:
        sym = "???"
    pct = 100.0 * count / total
    print(f"{count:8d}  {pct:5.1f}%  {offset:>10}  {sym}")
