"""Localize the gDNA->nascent siphon genomically. For every TRUE-gDNA fragment (read-name origin) that the
EM assigned to NASCENT (annotated-BAM ZF bit 0x08), classify its genomic location by the calibration region
signature: exon/intron x strand-composition (POS / NEG / AMBIG=both-strand overlap). Compare the SIPHONED
distribution to the BACKGROUND (all gDNA) to test: is the siphon enriched in AMBIG introns (opposite-strand
overlap), as hypothesized?
"""
from collections import Counter
import numpy as np
import pysam
from rigel.index import TranscriptIndex
from rigel.calibration.region_arrays import RegionArrays
from rigel.sim.read_name import parse_origin
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/annotated.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
starts = np.asarray(ra.start, np.int64); ends = np.asarray(ra.end, np.int64)
ref_off = np.asarray(ra.ref_offsets, np.int64); sig = np.asarray(ra.signature).astype(np.int64)
name2id = index.ref_name_to_id
POSb, NEGb = BIT_EXON_POS | BIT_INTRON_POS, BIT_EXON_NEG | BIT_INTRON_NEG
EXb, INb = BIT_EXON_POS | BIT_EXON_NEG, BIT_INTRON_POS | BIT_INTRON_NEG


def classify(ref, mid):
    rid = name2id.get(str(ref))
    if rid is None:
        return ("intergenic", "none")
    lo, hi = int(ref_off[rid]), int(ref_off[rid + 1])
    j = lo + int(np.searchsorted(ends[lo:hi], mid, side="right"))
    if not (lo <= j < hi) or not (starts[j] <= mid < ends[j]):
        return ("intergenic", "none")
    s = int(sig[j])
    kind = "exon" if (s & EXb) else ("intron" if (s & INb) else "intergenic")
    p, n = bool(s & POSb), bool(s & NEGb)
    strand = "AMBIG" if (p and n) else ("POS" if p else ("NEG" if n else "none"))
    return (kind, strand)


siphon = Counter(); background = Counter()
n_reads = 0
with pysam.AlignmentFile(bam, "rb") as f:
    default = f.references[0] if f.references else None
    seen = set()
    for r in f:
        n_reads += 1
        q = r.query_name
        if q in seen:
            continue
        seen.add(q)
        o = parse_origin(q)
        if o.kind != "gdna" or o.start is None:
            continue
        mid = (int(o.start) + int(o.end)) // 2
        cls = classify(o.ref if o.ref is not None else (r.reference_name or default), mid)
        background[cls] += 1
        zf = int(r.get_tag("ZF")) if r.has_tag("ZF") else 0
        if zf & 0x08:  # assigned to NASCENT
            siphon[cls] += 1

tot_s = sum(siphon.values()); tot_b = sum(background.values())
print(f"reads scanned={n_reads:,}  gDNA fragments={tot_b:,}  siphoned(gDNA->nascent)={tot_s:,} "
      f"({100*tot_s/max(tot_b,1):.1f}% of gDNA)")
print(f"\n{'(kind, strand)':24} {'siphon':>9} {'sip%':>6} {'background':>11} {'bg%':>6} {'enrich':>7}")
keys = sorted(set(siphon) | set(background), key=lambda k: -siphon.get(k, 0))
for k in keys:
    s, b = siphon.get(k, 0), background.get(k, 0)
    sp, bp = 100 * s / max(tot_s, 1), 100 * b / max(tot_b, 1)
    enr = (sp / bp) if bp > 0 else float("inf")
    print(f"{str(k):24} {s:>9,} {sp:>5.1f}% {b:>11,} {bp:>5.1f}% {enr:>6.2f}x")
# rollups
def roll(pred, label):
    s = sum(v for k, v in siphon.items() if pred(k)); b = sum(v for k, v in background.items() if pred(k))
    print(f"  {label:34} siphon {s:>8,} ({100*s/max(tot_s,1):4.1f}%)  bg {b:>9,} ({100*b/max(tot_b,1):4.1f}%)  "
          f"enrich {(100*s/max(tot_s,1))/(100*b/max(tot_b,1)) if b else float('nan'):.2f}x")
print("\nrollups:")
roll(lambda k: k[0] == "intron", "INTRON (any strand)")
roll(lambda k: k[0] == "exon", "EXON (any strand)")
roll(lambda k: k[1] == "AMBIG", "AMBIG (any kind)")
roll(lambda k: k[0] == "intron" and k[1] == "AMBIG", "AMBIG INTRON (the hypothesis)")
roll(lambda k: k[0] == "exon" and k[1] == "AMBIG", "AMBIG EXON")
