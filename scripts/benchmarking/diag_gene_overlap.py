#!/usr/bin/env python3
"""Diagnose gene overlap at the LINC02912 locus."""

GTF = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/"
    "PVT1_MYC/region.gtf"
)

# Parse GTF to find genes and transcripts overlapping LINC02912's exon
genes = {}
transcripts = {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    parts = line.strip().split("\t")
    if len(parts) < 9:
        continue
    feature = parts[2]
    start, end, strand = int(parts[3]), int(parts[4]), parts[6]
    attrs = parts[8]

    gid = attrs.split('gene_id "')[1].split('"')[0]
    gname = ""
    if 'gene_name "' in attrs:
        gname = attrs.split('gene_name "')[1].split('"')[0]

    if gid not in genes:
        genes[gid] = {"start": start, "end": end, "strand": strand, "name": gname}
    else:
        genes[gid]["start"] = min(genes[gid]["start"], start)
        genes[gid]["end"] = max(genes[gid]["end"], end)
        if gname:
            genes[gid]["name"] = gname

    if 'transcript_id "' in attrs:
        tid = attrs.split('transcript_id "')[1].split('"')[0]
        if feature == "exon":
            if tid not in transcripts:
                transcripts[tid] = {
                    "gene_id": gid, "gene_name": gname,
                    "strand": strand, "exons": [],
                }
            transcripts[tid]["exons"].append((start, end))

# LINC02912 region
region_start, region_end = 1501119, 1503283

print("=== Genes overlapping LINC02912 region (1501119-1503283) ===")
overlapping_genes = []
for gid, info in sorted(genes.items(), key=lambda x: x[1]["start"]):
    if info["start"] <= region_end and info["end"] >= region_start:
        overlapping_genes.append(gid)
        print(
            f"  {gid:30s} {info['name']:15s} {info['strand']} "
            f"{info['start']:>10d}-{info['end']:>10d}"
        )

print()
print("=== Transcripts with exons overlapping LINC02912 region ===")
for tid, info in sorted(transcripts.items(), key=lambda x: x[1]["exons"][0][0]):
    for (es, ee) in info["exons"]:
        if es <= region_end and ee >= region_start:
            n_exons = len(info["exons"])
            print(
                f"  {tid:30s} {info['gene_name']:15s} ({info['gene_id']}) "
                f"{info['strand']}  exon_at={es}-{ee}  n_exons={n_exons}"
            )
            break

print()
print("=== PCAT1 (ENSG00000253438.5) transcript details ===")
for tid, info in sorted(transcripts.items(), key=lambda x: x[1]["exons"][0][0]):
    if info["gene_id"] == "ENSG00000253438.5":
        exon_str = "; ".join(f"{s}-{e}" for s, e in info["exons"])
        print(f"  {tid:30s} n_exons={len(info['exons'])}  {exon_str}")

print()
print("=== LINC02912 (ENSG00000280055.2) transcript details ===")
for tid, info in sorted(transcripts.items(), key=lambda x: x[1]["exons"][0][0]):
    if info["gene_id"] == "ENSG00000280055.2":
        exon_str = "; ".join(f"{s}-{e}" for s, e in info["exons"])
        print(f"  {tid:30s} n_exons={len(info['exons'])}  {exon_str}")

print()
print("=== PVT1 (ENSG00000249859.14) details ===")
pvt1_info = genes.get("ENSG00000249859.14", {})
print(f"  PVT1 span: {pvt1_info.get('start')}-{pvt1_info.get('end')} {pvt1_info.get('strand')}")
