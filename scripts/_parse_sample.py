#!/usr/bin/env python3
"""Parse macOS `sample` output to extract per-thread and per-function breakdown."""
import re
import sys

def parse_sample(path):
    with open(path) as f:
        text = f.read()

    # 1. Overall library-level sample counts
    # The sample file has lines like: "  1234 ???  (in _em_impl.abi3.so) ..."
    # or "  567 std::exp  (in libsystem_m.dylib) ..."
    # Extract all (count, library) pairs from leaf-level samples
    print("=" * 72)
    print("LIBRARY-LEVEL BREAKDOWN (leaf sample counts)")
    print("=" * 72)

    lib_samples = {}
    for m in re.finditer(r'\b(\d+)\s+\S+\s+\(in\s+([^)]+)\)', text):
        count = int(m.group(1))
        lib = m.group(2)
        lib_samples[lib] = lib_samples.get(lib, 0) + count

    sorted_libs = sorted(lib_samples.items(), key=lambda x: -x[1])
    for lib, count in sorted_libs[:30]:
        print(f"  {count:>12,}  {lib}")

    # 2. Thread-level breakdown
    print()
    print("=" * 72)
    print("THREAD-LEVEL BREAKDOWN")
    print("=" * 72)

    # Split by thread headers
    thread_pattern = re.compile(r'^(\s+)(\d+)\s+(Thread_\d+)(?:\s+(.*))?$', re.MULTILINE)
    thread_entries = []
    for m in thread_pattern.finditer(text):
        count = int(m.group(2))
        tid = m.group(3)
        desc = m.group(4) or ""
        thread_entries.append((count, tid, desc.strip()))

    # Classify threads
    # Find position of each thread section header to extract content
    thread_starts = [(m.start(), int(m.group(2)), m.group(3))
                     for m in thread_pattern.finditer(text)]

    for i, (pos, count, tid) in enumerate(thread_starts):
        end_pos = thread_starts[i+1][0] if i+1 < len(thread_starts) else len(text)
        section = text[pos:end_pos]

        # Identify thread type
        if 'com.apple.main-thread' in section:
            ttype = "MAIN"
        elif '_em_impl' in section and '__ulock_wait' not in section[:500]:
            ttype = "EM_WORKER"
        elif '_bam_impl' in section:
            if 'bgzf' in section or 'tpool' in section:
                ttype = "HTS_DECOMP"
            else:
                ttype = "BAM_WORKER"
        elif 'bgzf' in section or 'tpool_worker' in section:
            ttype = "HTS_DECOMP"
        elif 'condition_variable::wait' in section and count > 100000:
            ttype = "IDLE_WORKER"
        elif count < 1000:
            ttype = "IDLE"
        else:
            ttype = "OTHER"

        thread_entries[i] = (count, tid, ttype)

    # Group by type
    type_groups = {}
    for count, tid, ttype in thread_entries:
        if ttype not in type_groups:
            type_groups[ttype] = []
        type_groups[ttype].append((count, tid))

    for ttype in ["MAIN", "EM_WORKER", "BAM_WORKER", "HTS_DECOMP", "IDLE_WORKER", "OTHER", "IDLE"]:
        items = type_groups.get(ttype, [])
        if not items:
            continue
        total = sum(c for c, _ in items)
        print(f"\n  {ttype} ({len(items)} threads, {total:,} total samples):")
        for count, tid in sorted(items, key=lambda x: -x[0])[:5]:
            print(f"    {count:>10,}  {tid}")
        if len(items) > 5:
            print(f"    ... +{len(items)-5} more")

    # 3. Deep dive into main thread
    print()
    print("=" * 72)
    print("MAIN THREAD DEEP DIVE")
    print("=" * 72)

    # Find the main thread section
    for i, (pos, count, tid) in enumerate(thread_starts):
        end_pos = thread_starts[i+1][0] if i+1 < len(thread_starts) else len(text)
        section = text[pos:end_pos]
        if 'com.apple.main-thread' not in section:
            continue

        # Extract top-level call breakdown within _em_impl, _bam_impl, Python
        for lib_name in ['_em_impl.abi3.so', '_bam_impl.abi3.so']:
            samples_in_lib = re.findall(r'(\d+)\s+\S+\s+\(in\s+' + re.escape(lib_name) + r'\)', section)
            total = sum(int(s) for s in samples_in_lib)
            if total > 0:
                print(f"\n  Main thread in {lib_name}: {total:,} samples")

        # Key operations in main thread
        for pattern_name, pat in [
            ("_em_impl calls", r'(\d+)\s+\?\?\?\s+\(in _em_impl'),
            ("_bam_impl calls", r'(\d+)\s+\?\?\?\s+\(in _bam_impl'),
            ("thread::join", r'(\d+)\s+std::thread::join'),
            ("__ulock_wait", r'(\d+)\s+__ulock_wait'),
            ("_pthread_join", r'(\d+)\s+_pthread_join'),
            ("condition_variable::wait", r'(\d+)\s+std::condition_variable::wait'),
            ("mutex_lock", r'(\d+)\s+std::mutex::lock'),
            ("cond_notify", r'(\d+)\s+std::condition_variable::notify'),
            ("std::exp", r'(\d+)\s+exp\s+\(in libsystem_m'),
            ("std::log", r'(\d+)\s+log\s+\(in libsystem_m'),
        ]:
            matches = re.findall(pat, section)
            total_s = sum(int(m) for m in matches)
            if total_s > 0:
                print(f"  {pattern_name:35s}: {total_s:>10,}")
        break

    # 4. Deep dive into EM worker threads (pick one with most samples)
    print()
    print("=" * 72)
    print("EM WORKER THREAD DEEP DIVE (highest-sample worker)")
    print("=" * 72)

    em_workers = [(i, pos, count, tid) for i, (pos, count, tid) in enumerate(thread_starts)
                  if thread_entries[i][2] == "EM_WORKER"]

    if em_workers:
        # Pick the one with most samples
        em_workers.sort(key=lambda x: -x[2])
        i, pos, count, tid = em_workers[0]
        end_pos = thread_starts[i+1][0] if i+1 < len(thread_starts) else len(text)
        section = text[pos:end_pos]
        print(f"  {tid}: {count:,} samples")

        for pattern_name, pat in [
            ("_em_impl code", r'(\d+)\s+\?\?\?\s+\(in _em_impl'),
            ("std::exp", r'(\d+)\s+exp\s+\(in libsystem_m'),
            ("std::log", r'(\d+)\s+log\s+\(in libsystem_m'),
            ("thread::join", r'(\d+)\s+std::thread::join'),
            ("__ulock_wait", r'(\d+)\s+__ulock_wait'),
            ("_pthread_join", r'(\d+)\s+_pthread_join'),
            ("condition_variable::wait", r'(\d+)\s+std::condition_variable::wait'),
            ("mutex contention", r'(\d+)\s+(?:mutex.*lock|_pthread_mutex)'),
            ("cond_notify", r'(\d+)\s+std::condition_variable::notify'),
            ("malloc/free", r'(\d+)\s+(?:malloc|free|realloc|_platform_m)'),
            ("memmove/memcpy", r'(\d+)\s+(?:memmove|memcpy|_platform_mem)'),
            ("libdeflate", r'(\d+)\s+\S+\s+\(in.*libdeflate'),
        ]:
            matches = re.findall(pat, section)
            total_s = sum(int(m) for m in matches)
            if total_s > 0:
                print(f"  {pattern_name:35s}: {total_s:>10,}")
    else:
        print("  (no EM worker threads identified)")

    # 5. BAM worker section
    print()
    print("=" * 72)
    print("BAM/SCAN THREAD DEEP DIVE")
    print("=" * 72)

    for ttype in ["BAM_WORKER"]:
        for i, (pos, count, tid) in enumerate(thread_starts):
            if thread_entries[i][2] != ttype:
                continue
            end_pos = thread_starts[i+1][0] if i+1 < len(thread_starts) else len(text)
            section = text[pos:end_pos]

            print(f"\n  {tid} ({ttype}): {count:,} samples")
            for pattern_name, pat in [
                ("_bam_impl code", r'(\d+)\s+\?\?\?\s+\(in _bam_impl'),
                ("condition_variable::wait", r'(\d+)\s+std::condition_variable::wait'),
                ("mutex contention", r'(\d+)\s+(?:std::mutex::lock|_pthread_mutex)'),
                ("std::condition_variable::notify", r'(\d+)\s+std::condition_variable::notify'),
                ("htslib", r'(\d+)\s+\S+\s+\(in.*libhts'),
                ("sam_read1", r'(\d+)\s+sam_read1'),
                ("bgzf", r'(\d+)\s+bgzf'),
                ("__ulock_wait", r'(\d+)\s+__ulock_wait'),
            ]:
                matches = re.findall(pat, section)
                total_s = sum(int(m) for m in matches)
                if total_s > 0:
                    print(f"    {pattern_name:35s}: {total_s:>10,}")
            break  # just first BAM worker

    # 6. Look at the HTS decompression threads
    hts_items = type_groups.get("HTS_DECOMP", [])
    if hts_items:
        print()
        print("=" * 72)
        print(f"HTS DECOMPRESSION ({len(hts_items)} threads, {sum(c for c,_ in hts_items):,} total)")
        print("=" * 72)
        # Look at first one
        for i, (pos, count, tid) in enumerate(thread_starts):
            if thread_entries[i][2] != "HTS_DECOMP":
                continue
            end_pos = thread_starts[i+1][0] if i+1 < len(thread_starts) else len(text)
            section = text[pos:end_pos]
            for pattern_name, pat in [
                ("cond_wait (idle)", r'(\d+)\s+_pthread_cond_wait'),
                ("bgzf_decode", r'(\d+)\s+bgzf_decode'),
                ("libdeflate", r'(\d+)\s+\S+\s+\(in.*libdeflate'),
            ]:
                matches = re.findall(pat, section)
                total_s = sum(int(m) for m in matches)
                if total_s > 0:
                    print(f"  {pattern_name:35s}: {total_s:>10,}")
            break

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else '/Users/mkiyer/Downloads/rigel_runs/profile_phase5_contaminated_sample_real.txt'
    parse_sample(path)
