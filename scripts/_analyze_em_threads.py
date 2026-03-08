#!/usr/bin/env python3
"""
Analyze macOS `sample` output using proper self-time (leaf node) extraction.

macOS sample output uses inclusive counts in a tree. Self-time for a frame
is its count minus the sum of its direct children's counts. Only self-time
(where CPU was actually executing) gives accurate attribution.
"""
import re
import sys
from collections import defaultdict

sample_file = (
    sys.argv[1]
    if len(sys.argv) > 1
    else "/Users/mkiyer/Downloads/rigel_runs/profile_phase5_contaminated_sample_real.txt"
)

with open(sample_file) as f:
    lines = f.readlines()

# --------------------------------------------------------------------------
# Step 1: Parse the tree into (depth, count, label, library) tuples
# --------------------------------------------------------------------------
# Each tree line has: prefix-drawing-chars COUNT FUNCNAME (in LIBRARY) ...
# The prefix uses +, !, :, | and spaces as tree connectors.
# Depth = column position of the first digit of COUNT.

LINE_RE = re.compile(
    r"^([\s+!:|]*?)\s*(\d+)\s+(\S.+?)\s+\(in ([^)]+)\)"
)

frames = []  # list of (depth, inclusive_count, func_name, library)
for line in lines:
    m = LINE_RE.match(line)
    if not m:
        continue
    # Depth = column position of count number (unambiguous)
    depth = m.start(2)
    count = int(m.group(2))
    func = m.group(3).strip()
    lib = m.group(4).strip()
    frames.append((depth, count, func, lib))

print(f"Parsed {len(frames):,} call-tree frames")

# --------------------------------------------------------------------------
# Step 2: Compute self-time for each frame
# A frame's self-time = its count - sum of direct children counts
# A "direct child" is the next frame at depth > current, before another
# frame at depth <= current.
# --------------------------------------------------------------------------
self_times = []  # (self_count, func, lib) for frames with self > 0

for i, (depth, count, func, lib) in enumerate(frames):
    children_sum = 0
    for j in range(i + 1, len(frames)):
        child_depth = frames[j][0]
        if child_depth <= depth:
            break  # sibling or parent, stop
        # Direct child: depth exactly one "level" deeper
        # In practice, find the first set of frames deeper than us
        # A direct child is one whose depth is the minimum depth > current
        pass

    # Simpler approach: find all immediate children
    # An immediate child is the next deeper frame(s) at the first deeper level
    child_total = 0
    if i + 1 < len(frames):
        # Find the child depth (first frame deeper than us)
        first_child_depth = None
        for j in range(i + 1, len(frames)):
            if frames[j][0] <= depth:
                break
            if first_child_depth is None:
                first_child_depth = frames[j][0]
            if frames[j][0] == first_child_depth:
                child_total += frames[j][1]
            elif frames[j][0] < first_child_depth:
                break  # back to parent level

    self_count = count - child_total
    if self_count > 0:
        self_times.append((self_count, func, lib))

total_self = sum(s for s, _, _ in self_times)
print(f"Total self-time samples: {total_self:,}")

# --------------------------------------------------------------------------
# Step 3: Aggregate by library
# --------------------------------------------------------------------------
print("\n" + "=" * 72)
print("SELF-TIME BY LIBRARY (top 20)")
print("=" * 72)

lib_self = defaultdict(int)
for s, func, lib in self_times:
    lib_self[lib] += s

for lib, total in sorted(lib_self.items(), key=lambda x: -x[1])[:20]:
    print(f"  {total:>12,}  ({100*total/total_self:5.1f}%)  {lib}")

# --------------------------------------------------------------------------
# Step 4: Aggregate by function (top 30)
# --------------------------------------------------------------------------
print("\n" + "=" * 72)
print("SELF-TIME BY FUNCTION (top 30)")
print("=" * 72)

func_self = defaultdict(int)
for s, func, lib in self_times:
    key = f"{func}  (in {lib})"
    func_self[key] += s

for key, total in sorted(func_self.items(), key=lambda x: -x[1])[:30]:
    print(f"  {total:>12,}  ({100*total/total_self:5.1f}%)  {key}")

# --------------------------------------------------------------------------
# Step 5: Per-thread self-time breakdown
# --------------------------------------------------------------------------
# Re-parse with thread boundaries
print("\n" + "=" * 72)
print("PER-THREAD SELF-TIME BREAKDOWN")
print("=" * 72)

THREAD_RE = re.compile(r"^\s+(\d+)\s+(Thread_\d+)(.*)?$")

# Find thread boundaries in the raw lines
thread_ranges = []  # (start_line, end_line, tid, count, desc)
for i, line in enumerate(lines):
    m = THREAD_RE.match(line)
    if m:
        count = int(m.group(1))
        tid = m.group(2)
        desc = (m.group(3) or "").strip()
        thread_ranges.append((i, None, tid, count, desc))
        if len(thread_ranges) > 1:
            thread_ranges[-2] = (
                thread_ranges[-2][0],
                i,
                thread_ranges[-2][2],
                thread_ranges[-2][3],
                thread_ranges[-2][4],
            )
if thread_ranges:
    thread_ranges[-1] = (
        thread_ranges[-1][0],
        len(lines),
        thread_ranges[-1][2],
        thread_ranges[-1][3],
        thread_ranges[-1][4],
    )


def parse_section_self_times(section_lines):
    """Parse a section of lines and return self-time by (func, lib)."""
    section_frames = []
    for line in section_lines:
        m = LINE_RE.match(line)
        if not m:
            continue
        depth = m.start(2)
        count_val = int(m.group(2))
        func_name = m.group(3).strip()
        lib_name = m.group(4).strip()
        section_frames.append((depth, count_val, func_name, lib_name))

    result = defaultdict(int)
    for i, (depth, count_val, func_name, lib_name) in enumerate(section_frames):
        child_total = 0
        if i + 1 < len(section_frames):
            first_child_depth = None
            for j in range(i + 1, len(section_frames)):
                if section_frames[j][0] <= depth:
                    break
                if first_child_depth is None:
                    first_child_depth = section_frames[j][0]
                if section_frames[j][0] == first_child_depth:
                    child_total += section_frames[j][1]
                elif section_frames[j][0] < first_child_depth:
                    break
        self_val = count_val - child_total
        if self_val > 0:
            result[(func_name, lib_name)] += self_val
    return result


# Classify and analyze threads
def classify_thread(self_by_func):
    """Classify a thread based on its self-time profile."""
    em_self = sum(v for (f, l), v in self_by_func.items() if "_em_impl" in l)
    bam_self = sum(v for (f, l), v in self_by_func.items() if "_bam_impl" in l)
    kmp_self = sum(
        v for (f, l), v in self_by_func.items() if "__kmp" in f or "libomp" in l
    )
    idle_self = sum(
        v
        for (f, l), v in self_by_func.items()
        if "__psynch_cvwait" in f or "__ulock_wait" in f
    )
    total = sum(self_by_func.values())

    if total == 0:
        return "IDLE", 0
    if kmp_self / total > 0.5:
        return "OMP_IDLE", total
    if em_self / total > 0.3:
        return "EM_WORKER", total
    if bam_self / total > 0.2:
        return "BAM_WORKER", total
    if idle_self / total > 0.5:
        return "IDLE", total
    return "OTHER", total


# Categorized totals
cat_threads = defaultdict(list)  # category -> [(tid, count, self_by_func)]

for start, end, tid, count, desc in thread_ranges:
    section_lines = lines[start:end]
    sbf = parse_section_self_times(section_lines)
    cat, active = classify_thread(sbf)

    if "com.apple.main-thread" in desc or "DispatchQueue_1" in desc:
        cat = "MAIN"

    cat_threads[cat].append((tid, count, sbf))

# Print by category
for cat in ["MAIN", "EM_WORKER", "BAM_WORKER", "OMP_IDLE", "IDLE", "OTHER"]:
    entries = cat_threads.get(cat, [])
    if not entries:
        continue

    cat_total_self = sum(sum(sbf.values()) for _, _, sbf in entries)
    print(f"\n--- {cat} ({len(entries)} threads, {cat_total_self:,} self-time samples) ---")

    # Aggregate self-time across all threads in this category
    agg = defaultdict(int)
    for _, _, sbf in entries:
        for k, v in sbf.items():
            agg[k] += v

    for (func, lib), total in sorted(agg.items(), key=lambda x: -x[1])[:15]:
        pct = 100 * total / cat_total_self if cat_total_self else 0
        label = f"{func}  (in {lib})"
        if len(label) > 70:
            label = label[:67] + "..."
        print(f"    {total:>10,}  ({pct:5.1f}%)  {label}")

# --------------------------------------------------------------------------
# Step 6: Wall-time projections
# --------------------------------------------------------------------------
print("\n" + "=" * 72)
print("WALL-TIME PROJECTIONS (contaminated dataset, 191.2s total)")
print("=" * 72)

runtime_s = 191.2
em_runtime_s = 138.8
scan_runtime_s = 44.4

# EM worker self-time breakdown
em_entries = cat_threads.get("EM_WORKER", [])
em_agg = defaultdict(int)
for _, _, sbf in em_entries:
    for k, v in sbf.items():
        em_agg[k] += v

em_total_self = sum(em_agg.values())
if em_total_self > 0:
    # exp() self-time
    exp_self = sum(v for (f, l), v in em_agg.items() if f == "exp" or "STUB$$exp" in f)
    em_code_self = sum(v for (f, l), v in em_agg.items() if "_em_impl" in l)
    other_self = em_total_self - exp_self - em_code_self

    NTHREADS = len(em_entries)
    # EM wall time is divided by thread count for per-thread wall time
    # But self-time is per-thread, so: projected = (self_frac) * em_runtime
    print(f"\nEM phase ({em_runtime_s}s wall, {NTHREADS} threads):")
    print(f"  Total EM self-time samples: {em_total_self:,}")
    print(
        f"  exp() self-time:   {exp_self:>10,}  ({100*exp_self/em_total_self:5.1f}%)"
        f"  => ~{exp_self/em_total_self * em_runtime_s:.1f}s wall"
    )
    print(
        f"  _em_impl code:     {em_code_self:>10,}  ({100*em_code_self/em_total_self:5.1f}%)"
        f"  => ~{em_code_self/em_total_self * em_runtime_s:.1f}s wall"
    )
    print(
        f"  other (sync/alloc):{other_self:>10,}  ({100*other_self/em_total_self:5.1f}%)"
        f"  => ~{other_self/em_total_self * em_runtime_s:.1f}s wall"
    )
