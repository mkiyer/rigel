#!/usr/bin/env python3
"""Debug tree parsing for one EM worker thread."""
import re
import sys

sample_file = "/Users/mkiyer/Downloads/rigel_runs/profile_phase5_contaminated_sample_real.txt"

with open(sample_file) as f:
    lines = f.readlines()

LINE_RE = re.compile(
    r"^([\s+!:|]*?)\s*(\d+)\s+(\S.+?)\s+\(in ([^)]+)\)"
)

# Find the EM worker thread Thread_9693565
in_thread = False
thread_lines = []
for line in lines:
    if "Thread_9693565" in line and not in_thread:
        in_thread = True
        thread_lines.append(line)
        continue
    if in_thread:
        # Check if we hit the next thread
        if re.match(r"^\s+\d+\s+Thread_\d+", line):
            break
        thread_lines.append(line)

print(f"Thread section: {len(thread_lines)} lines")
print()

# Parse and show depth for first 30 frames
frames = []
for line in thread_lines:
    m = LINE_RE.match(line)
    if not m:
        continue
    depth = len(m.group(1))
    count = int(m.group(2))
    func = m.group(3).strip()
    lib = m.group(4).strip()
    frames.append((depth, count, func, lib, line.rstrip()))

print("First 30 parsed frames (depth, count, func, lib):")
for i, (d, c, f, l, raw) in enumerate(frames[:30]):
    print(f"  [{i:3d}] depth={d:3d}  count={c:>8,}  {f[:50]}  (in {l})")

# Now compute self-time for first 30
print("\nSelf-time computation for first 30 frames:")
for i, (depth, count, func, lib, raw) in enumerate(frames[:30]):
    child_total = 0
    first_child_depth = None
    for j in range(i + 1, len(frames)):
        if frames[j][0] <= depth:
            break
        if first_child_depth is None:
            first_child_depth = frames[j][0]
        if frames[j][0] == first_child_depth:
            child_total += frames[j][1]
        elif frames[j][0] < first_child_depth:
            break
    self_time = count - child_total
    if self_time != 0:
        print(f"  [{i:3d}] depth={depth:3d}  count={count:>8,}  children_sum={child_total:>8,}  self={self_time:>8,}  {func[:40]}  (in {lib})")

# Show raw lines for first 15 frames to check depth
print("\nRaw first 15 frame lines:")
for i, (d, c, f, l, raw) in enumerate(frames[:15]):
    print(f"  [{i:3d}] depth={d:3d} | {raw[:120]}")
