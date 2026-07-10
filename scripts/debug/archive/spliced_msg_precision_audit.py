"""Step 2: audit the spliced/RNA message PRECISION — is it honest, tempered, and unable to overwhelm strand?

Captures per-node the combined RNA message (mode_p/prec_p, mode_n/prec_n) vs the local message-free strand
belief (fp_loc/vp_loc). Flags nodes where an RNA message with HIGH precision DISAGREES with a CONFIDENT
local strand belief and pulls it toward RNA (the over-injection). The disagreement-aware gate should make
prec small when the message disagrees with a confident local belief; verify it does.

    OMP_NUM_THREADS=1 python scripts/debug/spliced_msg_precision_audit.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402
import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402

COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
index, blob = build_or_load_cache(COND, False)
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

CAP = {}; CH = {}
_o = bp.node_sweep
def w(*a, **k):
    k["_capture"] = CAP; CH["c"] = a[0]; return _o(*a, **k)
bp.node_sweep = w
sys.modules["rigel.calibration.calibrate"].node_sweep = w
from rigel.calibration.calibrate import calibrate
from rigel.config import CalibrationConfig
cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
ch = CH["c"]; kind = np.asarray(ch.kind); rref = np.asarray(ch.ref_idx)
r2n = {int(rref[i]): i for i in range(len(kind)) if kind[i] == REGION}

sig = index.region_df["signature"].to_numpy()
scls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])
sg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
sr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
go = np.asarray(sg.contained.mass_unspliced, float); ro = np.asarray(sr.contained.mass_unspliced, float)

def arr(k): return np.asarray(CAP[k])
mode_p, prec_p = arr("mode_p"), arr("prec_p")
mode_n, prec_n = arr("mode_n"), arr("prec_n")
fp_loc, fn_loc = arr("fp_loc"), arr("fn_loc")
vp_loc, vn_loc = arr("vp_loc"), arr("vn_loc")
pp_loc = 1.0 / np.maximum(vp_loc, 1e-9); pn_loc = 1.0 / np.maximum(vn_loc, 1e-9)

# For each REGION node: does an RNA message overwhelm the local strand belief?
# message pulls toward RNA if mode (log f) > log(local f). "overwhelm" = prec_msg > prec_loc AND disagree.
print("=== spliced/RNA message vs local strand precision, by node class ===")
print("  ratio = prec_msg / prec_loc (per component, summed pos+neg); >1 means message can move the node")
for cname, mask in [("AMBIG exon", (scls == 3) & (tcls == 2)),
                    ("POS/NEG exon", np.isin(scls, [1, 2]) & (tcls == 2)),
                    ("POS/NEG intron", np.isin(scls, [1, 2]) & (tcls == 1))]:
    regs = [r for r in np.where(mask)[0] if r in r2n]
    nodes = np.array([r2n[r] for r in regs])
    if len(nodes) == 0:
        continue
    # combined RNA message precision vs combined local precision
    pm = prec_p[nodes] + prec_n[nodes]
    pl = pp_loc[nodes] + pn_loc[nodes]
    ratio = pm / np.maximum(pl, 1e-9)
    # oracle: is this node RNA-minor (mostly gDNA)? over-injection risk is high on RNA-minor nodes
    orc_fg = go[np.array(regs)] / (go[np.array(regs)] + ro[np.array(regs)] + 1e-9)
    rmin = orc_fg > 0.6
    print(f"  {cname:16s} n={len(nodes):>4}  median prec_msg/prec_loc={np.median(ratio):.2f}  "
          f"p90={np.percentile(ratio,90):.2f}  frac(ratio>1)={np.mean(ratio>1):.2f}  "
          f"| RNA-minor nodes(n={int(rmin.sum())}) median ratio={np.median(ratio[rmin]) if rmin.any() else float('nan'):.2f}")

# drill: RNA-minor AMBIG/POSNEG exons where the message is strong (potential over-inject)
print("\n=== nodes where RNA message is strong on an oracle-gDNA-rich exon (over-inject suspects) ===")
exon = (tcls == 2)
regs = [r for r in np.where(exon)[0] if r in r2n]
rows = []
for r in regs:
    n = r2n[r]; ofg = go[r] / (go[r] + ro[r] + 1e-9)
    if ofg < 0.6:  # only gDNA-rich (RNA-minor) exons
        continue
    pm = prec_p[n] + prec_n[n]; pl = pp_loc[n] + pn_loc[n]
    # message pushes RNA up if mode_p/mode_n (log f) exceeds local log f
    push = (mode_p[n] - np.log(max(fp_loc[n], 1e-9))) * (prec_p[n]) + (mode_n[n] - np.log(max(fn_loc[n], 1e-9))) * (prec_n[n])
    rows.append((r, ofg, pm, pl, pm / max(pl, 1e-9), push))
rows.sort(key=lambda x: -x[4])
print(f"  {'reg':>5} {'orc_fg':>7} {'prec_msg':>9} {'prec_loc':>9} {'ratio':>7} {'rna_push':>9}")
for r, ofg, pm, pl, ra_, pu in rows[:12]:
    print(f"  {r:>5} {ofg:>7.3f} {pm:>9.2f} {pl:>9.2f} {ra_:>7.2f} {pu:>+9.2f}")
