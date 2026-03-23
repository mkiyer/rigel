"""Analyze per-candidate resolution for multimapping GAPDH fragments.

Key question: When a fragment has primary (chr12, NM=0) and secondary (chrX, NM=5),
do the chr12 candidates include non-GAPDH transcripts? What fraction of candidate
transcripts are from pseudogenes?
"""
import sys
sys.path.insert(0, "src")

import pandas as pd
import pysam
from collections import defaultdict
from rigel.index import TranscriptIndex

base = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
idx_dir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
mm2_bam = f"{base}/sim_minimap2.bam"

# Load index
print("Loading index...")
tx_index = TranscriptIndex.load(idx_dir)
resolver = tx_index.resolver
ref_to_id = resolver.get_ref_to_id()

# Load transcript info
tx_df = pd.read_feather(f"{idx_dir}/transcripts.feather")
t_idx_to_info = {}
for _, r in tx_df.iterrows():
    t_idx_to_info[r["t_index"]] = {
        "t_id": r["t_id"], "g_name": r["g_name"],
        "ref": r["ref"], "start": r["start"], "end": r["end"],
    }

# Parse minimap2 BAM: group by qname, get all alignments
print("Parsing minimap2 BAM...")
qname_groups = defaultdict(list)
with pysam.AlignmentFile(mm2_bam, "rb", check_sq=False) as bam:
    for read in bam:
        nm = read.get_tag("NM") if read.has_tag("NM") else 0
        qname_groups[read.query_name].append({
            "chrom": read.reference_name,
            "pos": read.reference_start,
            "end": read.reference_end,
            "cigar": read.cigarstring,
            "flag": read.flag,
            "mapq": read.mapping_quality,
            "nm": nm,
            "is_r1": bool(read.flag & 0x40),
            "is_r2": bool(read.flag & 0x80),
            "is_sec": bool(read.flag & 0x100),
            "is_rev": bool(read.flag & 0x10),
        })

print(f"  {len(qname_groups)} fragments")

# For each multimapper, resolve candidates for primary and each secondary
# Then check if the EM can distinguish

STRAND_NONE = 0
STRAND_POS = 1
STRAND_NEG = 2

def cigar_to_exon_blocks(chrom, pos, cigar_str, is_rev):
    """Convert CIGAR to exon blocks and intron blocks."""
    import re
    exons = []
    introns = []
    p = pos
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str):
        length = int(length)
        if op in ('M', '=', 'X'):
            exons.append((p, p + length))
            p += length
        elif op == 'N':
            introns.append((p, p + length))
            p += length
        elif op == 'D':
            p += length
        elif op in ('I', 'S', 'H', 'P'):
            pass
    return exons, introns

# Analyze a sample of multimappers
print("\nResolving multimapper candidates...")

# Stats
n_analyzed = 0
n_unique_candidates = 0
n_multi_candidates = 0
gapdh_genes = {"GAPDH"}

# Per-fragment stats
frag_stats = []

for qname, reads in sorted(qname_groups.items()):
    primary = [r for r in reads if not r["is_sec"]]
    secondary = [r for r in reads if r["is_sec"]]
    if len(secondary) == 0:
        continue  # skip unique mappers

    truth_tx = qname.split(":")[0]
    
    # Group by alignment hit (R1+R2 pairs per alignment)
    # Primary is one hit, each secondary pair is another hit
    # But secondaries in minimap2 are individual reads, not paired
    
    # Resolve just the primary pair
    pri_r1 = [r for r in primary if r["is_r1"]]
    pri_r2 = [r for r in primary if r["is_r2"]]
    
    if not pri_r1 or not pri_r2:
        continue
    
    r1 = pri_r1[0]
    r2 = pri_r2[0]
    
    # Build combined exon/intron blocks for the fragment
    all_exons = []
    all_introns = []
    
    for read in [r1, r2]:
        chrom = read["chrom"]
        ref_id = ref_to_id.get(chrom, -1)
        if ref_id < 0:
            continue
        exons, introns = cigar_to_exon_blocks(chrom, read["pos"], read["cigar"], read["is_rev"])
        strand = STRAND_NONE
        for s, e in exons:
            all_exons.append((ref_id, s, e, strand))
        for s, e in introns:
            all_introns.append((ref_id, s, e, STRAND_POS))  # assume + strand for GAPDH

    # Resolve
    exon_ref_ids = [e[0] for e in all_exons]
    exon_starts = [e[1] for e in all_exons]
    exon_ends = [e[2] for e in all_exons]
    exon_strands = [e[3] for e in all_exons]
    intron_ref_ids = [i[0] for i in all_introns]
    intron_starts = [i[1] for i in all_introns]
    intron_ends = [i[2] for i in all_introns]
    intron_strands = [i[3] for i in all_introns]
    
    gfp = max(e[2] for e in all_exons) - min(e[1] for e in all_exons) if all_exons else 0
    
    result = resolver.resolve(
        exon_ref_ids=exon_ref_ids,
        exon_starts=exon_starts,
        exon_ends=exon_ends,
        exon_strands=exon_strands,
        intron_ref_ids=intron_ref_ids,
        intron_starts=intron_starts,
        intron_ends=intron_ends,
        intron_strands=intron_strands,
        genomic_footprint=gfp,
    )
    
    if result is None:
        pri_t_inds = []
    else:
        pri_t_inds = list(result[0])  # t_inds is element 0
    
    # Now resolve each secondary read individually
    sec_t_inds_by_read = []
    for read in secondary:
        chrom = read["chrom"]
        ref_id = ref_to_id.get(chrom, -1)
        if ref_id < 0:
            sec_t_inds_by_read.append((read, []))
            continue
        exons, introns = cigar_to_exon_blocks(chrom, read["pos"], read["cigar"], read["is_rev"])
        sec_exon_ref_ids = [ref_id] * len(exons)
        sec_exon_starts = [e[0] for e in exons]
        sec_exon_ends = [e[1] for e in exons]
        sec_exon_strands = [STRAND_NONE] * len(exons)
        sec_intron_ref_ids = [ref_id] * len(introns)
        sec_intron_starts = [i[0] for i in introns]
        sec_intron_ends = [i[1] for i in introns]
        sec_intron_strands = [STRAND_POS] * len(introns)
        
        sec_gfp = max(e[1] for e in exons) - min(e[0] for e in exons) if exons else 0
        
        sec_result = resolver.resolve(
            exon_ref_ids=sec_exon_ref_ids,
            exon_starts=sec_exon_starts,
            exon_ends=sec_exon_ends,
            exon_strands=sec_exon_strands,
            intron_ref_ids=sec_intron_ref_ids,
            intron_starts=sec_intron_starts,
            intron_ends=sec_intron_ends,
            intron_strands=sec_intron_strands,
            genomic_footprint=sec_gfp,
        )
        
        if sec_result is None:
            sec_t_inds_by_read.append((read, []))
        else:
            sec_t_inds_by_read.append((read, list(sec_result[0])))

    # Categorize candidates
    pri_genes = set()
    for ti in pri_t_inds:
        info = t_idx_to_info.get(ti, {})
        pri_genes.add(info.get("g_name", "?"))
    
    sec_genes = set()
    sec_all_t = set()
    for read, t_inds in sec_t_inds_by_read:
        for ti in t_inds:
            sec_all_t.add(ti)
            info = t_idx_to_info.get(ti, {})
            sec_genes.add(info.get("g_name", "?"))
    
    # Most secondaries have NM>0, primaries NM=0
    pri_nm = max(r["nm"] for r in primary)
    sec_nms = [r["nm"] for r in secondary]
    
    frag_stats.append({
        "qname": qname,
        "truth_tx": truth_tx,
        "n_pri_candidates": len(pri_t_inds),
        "n_sec_reads": len(secondary),
        "n_sec_candidates": len(sec_all_t),
        "pri_genes": pri_genes,
        "sec_genes": sec_genes,
        "pri_nm": pri_nm,
        "sec_nms": sec_nms,
        "pri_only_gapdh": pri_genes <= gapdh_genes,
    })
    
    n_analyzed += 1
    if n_analyzed >= 500:
        break

print(f"\nAnalyzed {n_analyzed} multimapping fragments")

# Summary
print(f"\n{'='*80}")
print(f"Primary Resolution Summary")
print(f"{'='*80}")

n_cand_hist = defaultdict(int)
for fs in frag_stats:
    n_cand_hist[fs["n_pri_candidates"]] += 1

print("\nPrimary candidates per fragment:")
for n, cnt in sorted(n_cand_hist.items()):
    print(f"  {n} candidates: {cnt} fragments")

# How many have ONLY GAPDH in primary?
n_pri_only_gapdh = sum(1 for fs in frag_stats if fs["pri_only_gapdh"])
n_pri_mixed = sum(1 for fs in frag_stats if not fs["pri_only_gapdh"])
print(f"\nPrimary candidates gene composition:")
print(f"  Only GAPDH: {n_pri_only_gapdh}/{n_analyzed}")
print(f"  Mixed (includes non-GAPDH): {n_pri_mixed}/{n_analyzed}")

# Show what non-GAPDH genes appear
non_gapdh_genes = defaultdict(int)
for fs in frag_stats:
    for g in fs["pri_genes"]:
        if g not in gapdh_genes:
            non_gapdh_genes[g] += 1

if non_gapdh_genes:
    print("\n  Non-GAPDH genes in primary candidates:")
    for g, cnt in sorted(non_gapdh_genes.items(), key=lambda x: -x[1]):
        print(f"    {g}: {cnt} fragments")

# Secondary resolution
print(f"\n  Secondary resolution genes:")
sec_gene_counts = defaultdict(int)
for fs in frag_stats:
    for g in fs["sec_genes"]:
        sec_gene_counts[g] += 1
for g, cnt in sorted(sec_gene_counts.items(), key=lambda x: -x[1])[:20]:
    print(f"    {g}: {cnt} fragments")

# Detailed examples
print(f"\n{'='*80}")
print(f"Example: Primary has both GAPDH and non-GAPDH candidates")
print(f"{'='*80}")

for fs in frag_stats:
    if not fs["pri_only_gapdh"] and len(fs["pri_genes"]) > 1:
        print(f"\n  {fs['qname']}")
        print(f"    Truth: {fs['truth_tx']}")
        print(f"    Primary NM={fs['pri_nm']}, {fs['n_pri_candidates']} candidates")
        print(f"    Primary genes: {fs['pri_genes']}")
        print(f"    Secondary: {fs['n_sec_reads']} reads, {fs['n_sec_candidates']} candidates")
        print(f"    Secondary genes: {fs['sec_genes']}")
        print(f"    Secondary NMs: {fs['sec_nms']}")
        break
