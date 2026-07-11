"""Verify (or refute) the capture gDNA fragment-length effect — the user's two controlled tests.

The earlier "on-target gDNA is SHORTER (345) than off-target (389)" claim is suspected to be a confound:
(a) capture biophysics + the simulator's own model (weight = binding_per_base × overlapping probe bases)
FAVOR LONGER fragments, so on-target should be LONGER; (b) classifying by fragment MIDPOINT-region is
length-biased (a long fragment over a short exon has its midpoint in the flanking intron ⇒ mislabeled
off-target). Here we measure true-origin gDNA FL classified by EXON OVERLAP (not midpoint), under:

  TEST 1 (no-RNA):  gDNA-dominant, capture OFF vs ON, input gDNA FL=250 — does capture shift gDNA FL UP?
  TEST 2 (flip):    input gDNA FL=100, RNA FL=200, capture ON — does captured gDNA FL rise ABOVE 100?

If capture favors LONGER gDNA in both, the "on-target shorter" signal was an artifact and the FL hypothesis
is refuted.

  OMP_NUM_THREADS=1 python scripts/debug/fl_capture_verify.py
"""
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.sim.read_name import parse_origin

WD = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/abe830c0-6786-4484-8f6b-96f6a75b0c35/scratchpad")
# genes define the capture probe targets (exons); long introns so the exon/intron gDNA split is meaningful.
GENES = [
    ("GA", "+", [("GA.1", [(1000, 4000), (15000, 18000), (28000, 31000)])]),
    ("GB", "-", [("GB.1", [(2000, 5000), (20000, 22000), (27000, 30000)])]),
]
EXONS = [(s, e) for _g, _st, isos in GENES for _t, exs in isos for (s, e) in exs]


def _overlaps_exon(a, b):
    return any(a < ee and b > es for (es, ee) in EXONS)


def run(name, capture, gdna_fm, rna_fm, rna_abundance, gdna_fraction):
    wd = WD / f"flcap_{name}"
    sc = Scenario(name, genome_length=33000, seed=7, work_dir=wd, ref_name="chr1")
    for gid, strand, isos in GENES:
        sc.add_gene(gid, strand, [{"t_id": t, "exons": ex, "abundance": rna_abundance} for t, ex in isos])
    cap_cfg = None
    if capture:
        wd.mkdir(parents=True, exist_ok=True)
        probes = wd / "probes.bed"
        with open(probes, "w") as fh:
            for gid, _s, isos in GENES:
                for t, ex in isos:
                    for i, (s, e) in enumerate(ex):
                        fh.write(f"chr1\t{s}\t{e}\t{t}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=20.0)
    result = sc.build_oracle(
        n_rna_fragments=6000, gdna_fraction=gdna_fraction,
        sim_config=ReadSimConfig(frag_mean=rna_fm, frag_std=50, frag_min=60, frag_max=700,
                                 read_length=100, strand_specificity=0.99, seed=7),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=gdna_fm, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=0.0)
    on, off, rna_n = [], [], 0
    with pysam.AlignmentFile(str(result.bam_path), "rb") as f:
        seen = set()
        for r in f:
            q = r.query_name
            if q in seen:
                continue
            seen.add(q)
            o = parse_origin(q)
            if o.start is None:
                continue
            if o.kind != "gdna":
                rna_n += 1
                continue
            L = int(o.end) - int(o.start)
            (on if _overlaps_exon(int(o.start), int(o.end)) else off).append(L)
    on, off = np.array(on), np.array(off)
    allg = np.concatenate([on, off]) if on.size + off.size else np.array([0.0])
    print(f"  {name:22} cap={'ON ' if capture else 'OFF'}  input gDNA FL={gdna_fm}  "
          f"| n_gdna={allg.size:>6} n_rna={rna_n:>5} | mean gDNA FL ALL={allg.mean():6.1f}  "
          f"exon-overlap={on.mean() if on.size else float('nan'):6.1f} (n={on.size})  "
          f"no-overlap={off.mean() if off.size else float('nan'):6.1f} (n={off.size})")


if __name__ == "__main__":
    print("=== TEST 1: no-RNA (gDNA-dominant), input gDNA FL=250 — does capture shift gDNA FL UP (favor long)? ===")
    run("t1_capOFF", False, gdna_fm=250, rna_fm=250, rna_abundance=0.0, gdna_fraction=0.98)
    run("t1_capON", True, gdna_fm=250, rna_fm=250, rna_abundance=0.0, gdna_fraction=0.98)
    print("\n=== TEST 2: flip (input gDNA FL=100, RNA FL=200), capture ON — captured gDNA FL should rise ABOVE 100 ===")
    run("t2_capOFF", False, gdna_fm=100, rna_fm=200, rna_abundance=100.0, gdna_fraction=0.5)
    run("t2_capON", True, gdna_fm=100, rna_fm=200, rna_abundance=100.0, gdna_fraction=0.5)
