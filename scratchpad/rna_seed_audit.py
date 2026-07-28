"""EXAMINE THE ACTUAL SPLICED FRAGMENTS — which splice types feed the RNA strand training set?

Owner (2026-07-28): *"I cannot believe you have splice junctions with sense/antisense ratios of 0.30-0.49.
That is ludicrous. You need to dig deeply. You'll need to examine some of these spliced fragments."*
and *"we should only be training on annotated splice junctions."*

THE CODE AUDIT ALREADY FOUND TWO THINGS (`bam_scanner.cpp:1433-1448`):

  ✅ SPLICE_ARTIFACT is correctly and completely excluded — `if (st == SPLICE_ARTIFACT) return;`,
     held out of the accumulator entirely, no deposit and no FL pool.
  ⚠  but the spliced channel is `SPLICED_ANNOT || SPLICED_UNANNOT || SPLICED_IMPLICIT` — so it includes
     UNANNOTATED junctions, which the owner says it should not.
  ⚠⚠ and SPLICED_IMPLICIT has NO SEQUENCED MOTIF at all. It is a paired-end fragment whose MATE GAP
     happens to span an annotated intron; the scanner then orients it by that intron's ANNOTATED
     transcript strand (`motif_strand |= isj.strand`). For a true RNA fragment that is ~right; for a gDNA
     fragment whose mate gap merely straddles an intron the orientation is arbitrary, i.e. ~50/50 — and a
     50/50 population is exactly what was measured (MO_3021 at depth >= 10 reads 0.4876).

THIS SCRIPT MEASURES IT PER FRAGMENT, from the buffer the scan produces: for each splice type, how many
fragments, and what fraction have `align_strand == sj_strand` (the "sense" bit the accumulator deposits).

PREDICTION: SPLICED_ANNOT sits at kappa (i.e. ~0.00 or ~1.00 on a stranded library); IMPLICIT and/or
UNANNOT sit near 0.5 and carry enough fragments to drag the pooled aggregate to the observed 0.30-0.49.
FALSIFIER: if every splice type sits at kappa, the pollution is not splice-type-related and the 0.30-0.49
comes from the accumulator's boundary bookkeeping instead.

Run: OMP_NUM_THREADS=1 python scratchpad/rna_seed_audit.py [--sample LBX0190]
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np

from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
_EPS = 1e-12
NAMES = {int(s): s.name for s in SpliceType}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="LBX0190")
    a = ap.parse_args()
    meta = pickle.load(open(CF / f"{a.sample}.pkl", "rb"))
    bam, index_dir = meta["bam"], meta["index_dir"]
    print(f"# sample={a.sample}\n# bam={bam}\n# index={index_dir}\n")
    index = TranscriptIndex.load(str(index_dir))
    cfg = PipelineConfig()
    stats, sm, flm, buf, _payload = scan_and_buffer(str(bam), index, cfg.scan)

    st = np.asarray(buf.splice_type, np.int64)
    al = np.asarray(buf.align_strand, np.int64)
    sj = np.asarray(buf.sj_strand, np.int64)

    print("=== A1. EVERY FRAGMENT, BY SPLICE TYPE — is 'sense' well defined, and what is it? ===")
    print("    sense = (align_strand == sj_strand), the exact bit the accumulator deposits.")
    print("    'both defined' = the deposit's own gate (align_ok && motif_ok).\n")
    print(f"{'splice type':22s} {'fragments':>12s} {'share':>7s} {'both defined':>13s} "
          f"{'SENSE FRACTION':>15s}")
    for t in sorted(set(st.tolist())):
        m = st == t
        ok = m & np.isin(al, (1, 2)) & np.isin(sj, (1, 2))
        sf = float(np.mean(al[ok] == sj[ok])) if ok.any() else float("nan")
        print(f"{NAMES.get(t, str(t)):22s} {int(m.sum()):12,d} {m.mean():6.1%} "
              f"{ok.sum() / max(int(m.sum()), 1):12.1%} {sf:15.4f}")

    print("\n\n=== A2. WHAT THE ACCUMULATOR ACTUALLY DEPOSITS AS 'SPLICED' ===")
    dep = np.isin(st, (int(SpliceType.SPLICED_ANNOT), int(SpliceType.SPLICED_UNANNOT),
                       int(SpliceType.SPLICED_IMPLICIT)))
    ok = dep & np.isin(al, (1, 2)) & np.isin(sj, (1, 2))
    print(f"  deposited-spliced fragments : {int(ok.sum()):,}")
    print(f"  pooled SENSE FRACTION       : {float(np.mean(al[ok] == sj[ok])):.4f}   "
          f"<- this is what trains od_r")
    ann = ok & (st == int(SpliceType.SPLICED_ANNOT))
    print(f"  ANNOTATED only              : {int(ann.sum()):,} fragments, "
          f"sense fraction {float(np.mean(al[ann] == sj[ann])):.4f}   <- what it SHOULD be")

    print("\n\n=== A3. THE COUNTERFACTUAL — annotated-only vs what ships ===")
    for lab, m in (("ships (ANNOT+UNANNOT+IMPLICIT)", ok),
                   ("ANNOT only", ann),
                   ("ANNOT + UNANNOT", ok & (st != int(SpliceType.SPLICED_IMPLICIT))),
                   ("IMPLICIT only", ok & (st == int(SpliceType.SPLICED_IMPLICIT))),
                   ("UNANNOT only", ok & (st == int(SpliceType.SPLICED_UNANNOT)))):
        if not m.any():
            print(f"  {lab:32s} (none)")
            continue
        print(f"  {lab:32s} n={int(m.sum()):10,d}  sense fraction {float(np.mean(al[m] == sj[m])):.4f}")

    print("\n  A splice motif is GT/AG — it is not ambiguous. Any population sitting near 0.5 is not")
    print("  measuring library-prep efficiency; it is measuring something whose orientation is arbitrary.")


if __name__ == "__main__":
    main()
