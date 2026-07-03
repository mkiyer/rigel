"""Audit the spliced-mass signal + its propagation vs oracle (dedicated-session step 1).

Uses the by-origin split scans (gdna-only / rna-only) from dissect_regions to get the oracle RNA/gDNA,
and the production calibrate() to get the solved values. Reports, in COUNTS/DENSITIES:
  (A) overall spliced availability: how much spliced vs unspliced mass, how many boundaries are junctions.
  (B) for AMBIG-exon regions: adjacent-boundary spliced mass (the RNA measurement that SHOULD anchor them),
      vs the region's oracle RNA and solved RNA.
  (C) for boundary nodes: spliced (motif RNA) vs oracle boundary RNA vs solved.

    OMP_NUM_THREADS=1 python scripts/debug/spliced_mass_audit.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402

COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
index, blob = build_or_load_cache(COND, False)
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
rdf = index.region_df.reset_index(drop=True)
sig = rdf["signature"].to_numpy()
scls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])

# oracle (by-origin) region masses
sg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
sr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
sf = CalibrationSubstrate.from_payload(blob["payload_full"], ra)
g_or = np.asarray(sg.contained.mass_unspliced, float)
r_or_uns = np.asarray(sr.contained.mass_unspliced, float)   # oracle RNA unspliced (nascent) in region interior
r_or_spl = np.asarray(sr.contained.mass_spliced, float)     # oracle RNA spliced in region interior (should be ~0)

# boundary substrates (full + oracle)
bf = BoundarySubstrate.from_payload(blob["payload_full"])
bg = BoundarySubstrate.from_payload(blob["payload_gdna"])
brn = BoundarySubstrate.from_payload(blob["payload_rna"])

# ---------- (A) overall spliced availability ----------
print("=" * 70)
print("(A) SPLICED SIGNAL AVAILABILITY (mass)")
tot_uns_reg = np.asarray(sf.contained.mass_unspliced, float).sum()
tot_spl_reg = np.asarray(sf.contained.mass_spliced, float).sum()
bspl = np.asarray(bf.left.mass_spliced, float) + np.asarray(bf.right.mass_spliced, float)
buns = np.asarray(bf.left.mass_unspliced, float) + np.asarray(bf.right.mass_unspliced, float)
print(f"  region interiors: unspliced={tot_uns_reg:,.0f}  spliced={tot_spl_reg:,.0f}")
print(f"  boundaries:       unspliced={buns.sum():,.0f}  spliced={bspl.sum():,.0f}")
nB = bspl.shape[0]
print(f"  boundaries total={nB}  with spliced>0 (junctions)={int((bspl>0).sum())}  "
      f"spliced mass frac of all fragments={bspl.sum()/(tot_uns_reg+tot_spl_reg+buns.sum()+bspl.sum()):.3f}")

# ---------- (B) AMBIG exons: adjacent-boundary spliced vs region oracle RNA ----------
print("=" * 70)
print("(B) AMBIG-EXON regions: adjacent-boundary SPLICED (motif RNA) vs region oracle")
lr = np.asarray(bf.left_region); rr = np.asarray(bf.right_region)
bspl_sense = np.asarray(bf.left.n_spliced_sense, float) + np.asarray(bf.right.n_spliced_sense, float)
bspl_anti = np.asarray(bf.left.n_spliced_antisense, float) + np.asarray(bf.right.n_spliced_antisense, float)
# map region -> adjacent boundaries
from collections import defaultdict
reg_bnds = defaultdict(list)
for b in range(nB):
    if lr[b] >= 0: reg_bnds[int(lr[b])].append(b)
    if rr[b] >= 0: reg_bnds[int(rr[b])].append(b)
ae = [r for r in range(len(sig)) if scls[r] == 3 and tcls[r] == 2]
print(f"  {'reg':>5} {'orc_gDNA':>9} {'orc_RNA':>8} {'orcRNAspl':>9} | {'adjB':>4} {'adjBspl_sense':>13} {'adjBspl_anti':>12}")
for r in sorted(ae, key=lambda r: -(g_or[r]+r_or_uns[r]))[:10]:
    bs = reg_bnds.get(r, [])
    ss = sum(bspl_sense[b] for b in bs); aa = sum(bspl_anti[b] for b in bs)
    print(f"  {r:>5} {g_or[r]:>9,.0f} {r_or_uns[r]:>8,.0f} {r_or_spl[r]:>9,.0f} | {len(bs):>4} {ss:>13,.0f} {aa:>12,.0f}")
# aggregate: do AMBIG exons have adjacent spliced at all?
adj_spl = np.array([sum(bspl[b] for b in reg_bnds.get(r, [])) for r in ae])
print(f"  AMBIG exons with ANY adjacent spliced: {(adj_spl>0).sum()}/{len(ae)}   total adj spliced mass={adj_spl.sum():,.0f}")
print(f"  AMBIG-exon oracle: gDNA={g_or[ae].sum():,.0f}  RNA_unspliced={r_or_uns[ae].sum():,.0f}  (RNA is {'MINOR' if r_or_uns[ae].sum()<g_or[ae].sum() else 'MAJOR'})")

# ---------- (C) boundary nodes: spliced vs oracle boundary composition ----------
print("=" * 70)
print("(C) BOUNDARY nodes: oracle gDNA/RNA(unspl)/RNA(spl) on the crossing")
bg_gdna = np.asarray(bg.left.mass_unspliced, float) + np.asarray(bg.right.mass_unspliced, float)
brn_uns = np.asarray(brn.left.mass_unspliced, float) + np.asarray(brn.right.mass_unspliced, float)
brn_spl = np.asarray(brn.left.mass_spliced, float) + np.asarray(brn.right.mass_spliced, float)
print(f"  boundary oracle totals: gDNA(unspl xing)={bg_gdna.sum():,.0f}  RNA(unspl nascent)={brn_uns.sum():,.0f}  RNA(spliced)={brn_spl.sum():,.0f}")
# are boundaries with spliced ALSO high RNA? and boundaries WITHOUT spliced mostly gDNA?
has_spl = bspl > 0
print(f"  boundaries WITH spliced (n={int(has_spl.sum())}):  oracle gDNA={bg_gdna[has_spl].sum():,.0f}  RNA_uns={brn_uns[has_spl].sum():,.0f}")
print(f"  boundaries NO spliced  (n={int((~has_spl).sum())}):  oracle gDNA={bg_gdna[~has_spl].sum():,.0f}  RNA_uns={brn_uns[~has_spl].sum():,.0f}")
