#!/usr/bin/env python3
"""Run profiler and capture native sample during EM phase.

Usage: python scripts/_profile_with_sample.py <bam> <index> <outdir> <threads>

Forks the profiler, monitors log output for the EM phase start marker,
then runs macOS `sample` for 30s to capture the C++ EM hot path.
"""
import os
import subprocess
import sys
import time
import signal

def main():
    bam = sys.argv[1]
    index = sys.argv[2]
    outdir = sys.argv[3]
    threads = sys.argv[4]

    os.makedirs(outdir, exist_ok=True)

    # Start the profiler as a subprocess
    cmd = [
        sys.executable, "scripts/profiler.py",
        "--bam", bam,
        "--index", index,
        "--outdir", outdir,
        "--stages",
        "--threads", threads,
        "--verbose",
    ]
    print(f"Starting profiler: {' '.join(cmd)}")
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    pid = proc.pid
    print(f"Profiler PID: {pid}")

    # Monitor output for EM phase start markers
    em_started = False
    sample_file = os.path.join(outdir, "native_sample.txt")

    for line in proc.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()

        # Detect when locus EM begins (after fragment_router_scan completes)
        if "eb_gdna_priors" in line.lower() or "nrna_frac priors" in line.lower():
            if not em_started:
                em_started = True
                print(f"\n>>> EM phase detected! Starting 30s native sample of PID {pid}")
                # Run sample in background
                sample_proc = subprocess.Popen(
                    ["sample", str(pid), "30", "-f", sample_file],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )

    proc.wait()
    print(f"\nProfiler exited with code {proc.returncode}")

    if em_started:
        sample_proc.wait()
        print(f"Native sample written to {sample_file}")
    else:
        print("WARNING: EM phase marker not detected")

if __name__ == "__main__":
    main()
