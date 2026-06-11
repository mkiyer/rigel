"""3-pool net fragment flow: gDNA <-> nRNA <-> mRNA (not the old 2-pool gDNA<->RNA).

The capture eff-length fix resurrected the nRNA component; the 2-pool (gDNA vs RNA) leak conflates nRNA
with mRNA. This reports the full 3-pool picture per condition: the 3x3 gross flow (true origin ×
assigned pool), the three NET flows (symmetric misassignment cancels — only systematic bias survives),
and the key release metric — the **mature-RNA false-positive rate** (gDNA or nRNA mislabeled as mRNA).

True origin: parse_origin(read name) → gdna/nrna/mrna. Assigned pool: ZF bit 0x04 (or ZC=intergenic) =
gDNA, ZF bit 0x08 = nRNA, else mRNA — the same logic as rigel.sim.analysis.collect_fragment_flows.
"""
import sys
from collections import Counter
import pysam
from rigel.sim.read_name import parse_origin

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
POOLS = ["gdna", "nrna", "mrna"]
conds = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_rnd_capture_on", "gdna_gdna300_ss_0.50_nrna_rnd_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on", "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.99_nrna_rnd_capture_on", "gdna_none_ss_0.50_nrna_rnd_capture_on",
]


def pool_flow(cond):
    f = Counter()
    with pysam.AlignmentFile(f"{SUITE}/{cond}/annotated.bam", "rb") as bam:
        for r in bam:
            if r.is_read2 or r.is_secondary or r.is_supplementary:
                continue
            true = parse_origin(r.query_name).kind
            zf = int(r.get_tag("ZF")) if r.has_tag("ZF") else 0
            zc = str(r.get_tag("ZC")) if r.has_tag("ZC") else ""
            if (zf & 0x04) or zc == "intergenic":
                ap = "gdna"
            elif zf & 0x08:
                ap = "nrna"
            else:
                ap = "mrna"
            f[(true, ap)] += 1
    return f


for cond in conds:
    try:
        f = pool_flow(cond)
    except FileNotFoundError:
        print(f"{cond}: no annotated.bam")
        continue
    tot = sum(f.values())
    exp = {p: sum(f[(p, b)] for b in POOLS) for p in POOLS}  # true totals
    obs = {p: sum(f[(a, p)] for a in POOLS) for p in POOLS}  # assigned totals
    print(f"\n=== {cond}  (n={tot:,}) ===")
    print(f"{'true\\assigned':>14} " + "".join(f"{p:>12}" for p in POOLS) + f"{'TRUE_TOTAL':>12}")
    for tp in POOLS:
        row = "".join(f"{f[(tp,ap)]:>12,}" for ap in POOLS)
        print(f"{tp:>14} {row}{exp[tp]:>12,}")
    print(f"{'ASSIGNED_TOT':>14} " + "".join(f"{obs[p]:>12,}" for p in POOLS))
    # net flows (a->b) - (b->a)
    def net(a, b):
        return f[(a, b)] - f[(b, a)]
    # NET flux per pool pair (symmetric, truly-unidentifiable misassignment cancels — only systematic
    # bias survives). Sign: + means net flow from the first pool to the second.
    print(f"  NET flux:  gDNA->nRNA={net('gdna','nrna'):>+9,}   gDNA->mRNA={net('gdna','mrna'):>+9,}   "
          f"nRNA->mRNA={net('nrna','mrna'):>+9,}")
    # Per-pool NET surplus = assigned - true. + = the pool is net-inflated (false gain); - = net-deficit.
    # The mature-RNA "false positive" is the NET mRNA surplus (not the gross sum) — symmetric exchange
    # with sequence-identical fragments cancels.
    print(f"  NET surplus (assigned-true):  gDNA={obs['gdna']-exp['gdna']:>+9,}   "
          f"nRNA={obs['nrna']-exp['nrna']:>+9,}   mRNA={obs['mrna']-exp['mrna']:>+9,}"
          + (f"   (mRNA {100*(obs['mrna']-exp['mrna'])/max(exp['mrna'],1):+.1f}% of true mRNA)"
             if exp['mrna'] else ""))
