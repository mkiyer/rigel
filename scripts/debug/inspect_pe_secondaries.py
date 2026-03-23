#!/usr/bin/env python3
"""Inspect PE secondary alignment structure for one GAPDH fragment."""
import pysam

bam = '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file/sim_minimap2.bam'

# Pick a specific fragment to examine
target_qname = 'ENST00000229239.10:0-217:f:1698'

with pysam.AlignmentFile(bam, 'rb', check_sq=False) as f:
    n = 0
    for r in f:
        if r.query_name == target_qname:
            flag = r.flag
            tag = 'R1' if flag & 0x40 else ('R2' if flag & 0x80 else 'SE')
            sec = 'sec' if flag & 0x100 else 'pri'
            rev = '-' if flag & 0x10 else '+'
            mrev = '-' if flag & 0x20 else '+'
            pp = 'proper' if flag & 0x2 else 'improper'
            hi = r.get_tag('HI') if r.has_tag('HI') else -1
            mate_chr = r.next_reference_name or '*'
            mate_pos = r.next_reference_start
            nm = r.get_tag('NM') if r.has_tag('NM') else -1
            cigar = r.cigarstring
            print(f"  {tag} {sec:3s} {r.reference_name}:{r.reference_start:>10d} {rev} "
                  f"mate={mate_chr}:{mate_pos:>10d} {mrev} {pp:8s} "
                  f"flag={flag:>5d} HI={hi} cigar={cigar} NM={nm}")
            n += 1
    print(f"\nTotal records for this qname: {n}")
