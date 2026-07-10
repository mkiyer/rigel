"""Quantify spliced-mass ownership: region_contained (c2/c3) vs boundary objects (c2/c3)."""
import dataclasses, sys
import numpy as np
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_rnd_capture_on"
bam = f"{SUITE}/{cond}/sim_oracle.bam"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
cfg = PipelineConfig()
scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
_st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)

rc = pl.region_contained.astype(np.float64)          # (R,4) counts (==mass)
bml = pl.boundary_mass_left.astype(np.float64)        # (B,4) float mass
bmr = pl.boundary_mass_right.astype(np.float64)

# Channels: 0,1 = unspliced +/-, 2,3 = spliced sense/antisense
region_unspl = rc[:, 0:2].sum()
region_spl   = rc[:, 2:4].sum()
bnd_unspl    = bml[:, 0:2].sum() + bmr[:, 0:2].sum()
bnd_spl      = bml[:, 2:4].sum() + bmr[:, 2:4].sum()

total_spl   = region_spl + bnd_spl
total_unspl = region_unspl + bnd_unspl
total_mass  = total_spl + total_unspl

print(f"=== {cond} ===")
print(f"  total deposited mass        : {total_mass:,.1f}")
print(f"  SPLICED  total              : {total_spl:,.1f}")
print(f"     in region_contained      : {region_spl:,.1f}  ({100*region_spl/max(total_spl,1):.3f}% of spliced)")
print(f"     in boundary objects      : {bnd_spl:,.1f}  ({100*bnd_spl/max(total_spl,1):.3f}% of spliced)")
print(f"  UNSPLICED total             : {total_unspl:,.1f}")
print(f"     in region_contained      : {region_unspl:,.1f}  ({100*region_unspl/max(total_unspl,1):.3f}%)")
print(f"     in boundary objects      : {bnd_unspl:,.1f}  ({100*bnd_unspl/max(total_unspl,1):.3f}%)")
# how many regions carry ANY contained spliced?
n_reg_with_spl = int((rc[:, 2:4].sum(1) > 0).sum())
print(f"  regions with contained spliced > 0 : {n_reg_with_spl} / {rc.shape[0]}")
# boundary one-sidedness: of boundaries with spliced, what fraction is one-sided?
bspl_l = bml[:, 2:4].sum(1); bspl_r = bmr[:, 2:4].sum(1)
has = (bspl_l + bspl_r) > 0
onesided = np.where(has, np.maximum(bspl_l, bspl_r) / np.maximum(bspl_l + bspl_r, 1e-9), np.nan)
print(f"  boundaries with spliced: {int(has.sum())}; mean one-sidedness = {np.nanmean(onesided):.4f} (1.0=fully one-sided)")
