"""Phase-2 validation: sweep the FP-rate quantile q and watch the 3-pool NET flow trade.

For each (condition, q) it runs `rigel quant` at `--gdna-deconv-quantile q` into a scratch dir,
computes the 3-pool NET surplus (gDNA / nRNA / mRNA assigned − true) from the annotated BAM, then
DELETES the scratch dir (the annotated BAMs are ~0.5 GB each — stream, don't accumulate).

Claim under test (phase2_design.md §5): as q rises the gDNA→RNA *leak* falls and the RNA→gDNA
*siphon* rises, **monotonically**, and the trade concentrates on strand-observable nodes (the stranded
condition moves more than the count-routed unstranded one). q=0.5 is the no-op baseline.
"""
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pysam

from rigel.sim.read_name import parse_origin

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb")
INDEX = SUITE / "rigel_index"
SCRATCH = Path("/tmp/qsweep")
POOLS = ["gdna", "nrna", "mrna"]
QS = [0.5, 0.7, 0.9, 0.95]
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_rnd_capture_on",  # stranded   — quantile bites via strand posterior
    "gdna_gdna300_ss_0.50_nrna_rnd_capture_on",  # unstranded — count-routed (the hard floor case)
]


def run_quant(cond: str, q: float) -> Path:
    out = SCRATCH / f"{cond}__q{q}"
    out.mkdir(parents=True, exist_ok=True)
    ann = out / "annotated.bam"
    cmd = [
        "rigel", "quant",
        "--bam", str(SUITE / cond / "sim_oracle.bam"),
        "--index", str(INDEX),
        "-o", str(out),
        "--annotated-bam", str(ann),
        "--sj-strand-tag", "auto",
        "--gdna-deconv-quantile", str(q),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"quant failed for {cond} q={q}:\n{r.stderr[-2000:]}")
    return ann


def pool_net(ann: Path) -> dict:
    f = Counter()
    with pysam.AlignmentFile(str(ann), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            true = parse_origin(read.query_name).kind
            zf = int(read.get_tag("ZF")) if read.has_tag("ZF") else 0
            zc = str(read.get_tag("ZC")) if read.has_tag("ZC") else ""
            ap = "gdna" if ((zf & 0x04) or zc == "intergenic") else ("nrna" if zf & 0x08 else "mrna")
            f[(true, ap)] += 1
    exp = {p: sum(f[(p, b)] for b in POOLS) for p in POOLS}
    obs = {p: sum(f[(a, p)] for a in POOLS) for p in POOLS}

    def net(a, b):
        return f[(a, b)] - f[(b, a)]

    return {
        "leak_gdna_to_rna": net("gdna", "mrna") + net("gdna", "nrna"),  # gDNA mislabeled as any RNA
        "net_gdna_to_nrna": net("gdna", "nrna"),
        "net_gdna_to_mrna": net("gdna", "mrna"),
        "net_nrna_to_mrna": net("nrna", "mrna"),
        "surplus_gdna": obs["gdna"] - exp["gdna"],
        "surplus_nrna": obs["nrna"] - exp["nrna"],
        "surplus_mrna": obs["mrna"] - exp["mrna"],
        "exp_mrna": exp["mrna"],
    }


for cond in CONDS:
    print(f"\n=== {cond} ===")
    print(f"{'q':>6} {'leak gDNA→RNA':>15} {'net g→nRNA':>12} {'net g→mRNA':>12} "
          f"{'surplus gDNA':>13} {'surplus nRNA':>13} {'surplus mRNA':>13}")
    for q in QS:
        ann = run_quant(cond, q)
        m = pool_net(ann)
        shutil.rmtree(ann.parent, ignore_errors=True)  # stream: reclaim ~0.5 GB immediately
        mrna_pct = f"({100*m['surplus_mrna']/max(m['exp_mrna'],1):+.1f}%)"
        print(f"{q:>6} {m['leak_gdna_to_rna']:>+15,} {m['net_gdna_to_nrna']:>+12,} "
              f"{m['net_gdna_to_mrna']:>+12,} {m['surplus_gdna']:>+13,} {m['surplus_nrna']:>+13,} "
              f"{m['surplus_mrna']:>+9,} {mrna_pct:>6}")

shutil.rmtree(SCRATCH, ignore_errors=True)
print("\n(scratch cleaned)")
