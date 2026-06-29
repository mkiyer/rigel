"""Solidify the hump-shaped-dispersion finding: (1) on the CONVERGED belief (what actually drives sweep
messages), not just init; (2) rule out a depleted-FLOOR artifact (edges touching the 1/E_gdna clamp);
(3) quantify the monotone-increasing mis-fit precisely vs a free (non-monotone) smoother."""
from pathlib import Path
import numpy as np
from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import build_node_chain
from rigel.calibration import bp_solver as bp
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.splice import SpliceType

ROOT = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
BAM = ROOT / "gdna_gdna300_ss_0.99_nrna_none_capture_on/sim_oracle.bam"
idx = TranscriptIndex.load(str(ROOT / "rigel_index"))
st, sm, flm, buf, pl = scan_and_buffer(str(BAM), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
statics = bp.build_node_statics(chain, sub, bsub, ra)
kappa = float(fit_strand_balance(sm).rna_sense_frac)
belief0 = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
belief, _ = bp.node_sweep(chain, statics, geom, belief0, ra, bsub, rna_sense_frac=kappa,
                          n_grid=60, max_outer=CalibrationConfig().sweep_max_outer)
dens = bp.node_densities(belief, geom)
gdna_vm = bp.fit_gdna_varmean(chain, dens, geom, statics)

left = np.asarray(chain.left); right = np.asarray(chain.right)
live = np.asarray(statics.strand_obs, bool)
rgl = np.asarray(dens.rho_g_left); rgr = np.asarray(dens.rho_g_right)
EGl = np.asarray(geom.eff_gdna_left); EGr = np.asarray(geom.eff_gdna_right)
floor = 2.0 / np.maximum(np.where(EGl > 0, EGl, EGr), 1e-9)  # ~2× the 1/E_gdna min-observable density
rows = []
for nbr, d_rho, s_rho, fl_d in ((left, rgl, rgr, floor), (right, rgr, rgl, floor)):
    for i in np.where((nbr >= 0) & live)[0]:
        s = int(nbr[i])
        if not live[s]:
            continue
        dr, sr = d_rho[i], s_rho[s]
        if dr <= 0 or sr <= 0:
            continue
        at_floor = (dr < fl_d[i]) or (sr < fl_d[s])
        rows.append((0.5 * (dr + sr), (np.log(dr) - np.log(sr)) ** 2, dr, sr, int(at_floor)))
rows = np.array(rows)
mu, raw, dr, sr, atf = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4].astype(bool)
print(f"=== CONVERGED belief: {mu.size:,} gDNA edges  (at-floor={int(atf.sum()):,}, off-floor={int((~atf).sum()):,}) ===")

print("\n(1) raw=(Δlogρ)² vs μ  — ALL edges vs OFF-FLOOR only (rule out clamp artifact)")
qs = np.quantile(mu, np.linspace(0, 1, 9))
for a, b in zip(qs[:-1], qs[1:]):
    m = (mu >= a) & (mu <= b) if b == qs[-1] else (mu >= a) & (mu < b)
    off = m & ~atf
    s2 = float(np.median(gdna_vm.predict(mu[m]))) if m.sum() else 0.0
    print(f"  μ∈[{a:8.4f},{b:8.4f})  n={int(m.sum()):>5}  median raw(all)={np.median(raw[m]):6.2f}  "
          f"median raw(off-floor n={int(off.sum()):>4})={np.median(raw[off]) if off.sum() else float('nan'):6.2f}  "
          f"fitted σ²_bio={s2:5.2f}")

# (2) monotone vs FREE smoother: bin raw by μ, is the per-bin median non-monotone (hump)?
print("\n(2) free per-bin median raw (off-floor) — the shape a monotone-INCREASING spline must distort")
off = ~atf
qo = np.quantile(mu[off], np.linspace(0, 1, 11))
meds = []
for a, b in zip(qo[:-1], qo[1:]):
    m = off & ((mu >= a) & (mu <= b) if b == qo[-1] else (mu >= a) & (mu < b))
    meds.append(np.median(raw[m]) if m.sum() else np.nan)
meds = np.array(meds)
print("  per-decile median raw:", " ".join(f"{v:5.1f}" for v in meds))
peak = int(np.nanargmax(meds))
print(f"  PEAK at decile {peak}/10 (μ≈{0.5*(qo[peak]+qo[peak+1]):.3f}); "
      f"{'HUMP (non-monotone) — monotone-increasing CANNOT fit' if 0 < peak < len(meds)-1 else 'monotone-ish'}")
print(f"\n  fit: edf={gdna_vm.edf:.1f} lam={gdna_vm.lam:.3g} σ²_bio∈[{gdna_vm.predict(mu).min():.3f},{gdna_vm.predict(mu).max():.3f}]"
      f"  peak true raw={np.nanmax(meds):.1f} → under-fit {np.nanmax(meds)/max(gdna_vm.predict(mu).max(),1e-9):.1f}×")
