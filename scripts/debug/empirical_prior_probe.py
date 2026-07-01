"""Confirm the unstranded-intron-leak mechanism + preview the intergenic+intron-REGION floor.

Builds the real calibration objects on a toy, then for each REGION reports its observed gDNA
density M/E_g, and computes:
  * CURRENT global: rho_global + sigma2_g (the current seed set) → n_glob, and the implied f_g
    target for each region (what the global prior pins it toward).
  * PROPOSED floor: rho_global + spread from intergenic+intron REGIONS only (the user's directive)
    → n_glob_floor, implied f_g target.
Shows whether a depleted intron's leak is a global-precision (spread) problem and whether the
region-floor pins it to f_g≈1.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.substrate import CalibrationSubstrate, BoundarySubstrate
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration.bp_solver import (
    build_node_geometry, build_node_statics, init_beliefs, _gdna_seed_estimate,
    _node_region_type,
)
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.effective_length import region_eff_length
from rigel.splice import SpliceType

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/e1a42517-d9da-4b58-b12a-0b84f986ddef/scratchpad")


def _write_probes(path, ref, genes):
    with open(path, "w") as f:
        for gid, _s, exons in genes:
            for i, (s, e) in enumerate(exons):
                f.write(f"{ref}\t{s}\t{e}\t{gid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")


def probe(name, genes, *, kappa, capture, n_rna=20000, gdna_fraction=0.5, capture_strength=20.0,
          genome_length=12000, seed=7):
    wd = SCRATCH / f"eprior_{name}"
    sc = Scenario(name, genome_length=genome_length, seed=seed, work_dir=wd, ref_name="chr1")
    for gid, strand, exons in genes:
        sc.add_gene(gid, strand, [{"t_id": gid, "exons": exons, "abundance": 100.0}])
    cap_cfg = None
    if capture:
        probes = wd / "probes.bed"; wd.mkdir(parents=True, exist_ok=True)
        _write_probes(probes, "chr1", genes)
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=capture_strength)
    res = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=gdna_fraction,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=120, strand_specificity=kappa, seed=seed),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=250, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=0.0,
    )
    idx = res.index
    _st, sm, flm, _buf, pl = scan_and_buffer(str(res.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geo = build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
    st = build_node_statics(chain, sub, bsub, ra)
    kfit = float(fit_strand_balance(sm).rna_sense_frac)
    belief = init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kfit, n_grid=60, statics=st)

    rho_g_cur, gdna_vm, var_mean = _gdna_seed_estimate(chain, st, geo, ra, bsub, belief.f_g, kfit)
    s2_between_cur = max(float(gdna_vm.predict(np.array([max(rho_g_cur, 1e-9)]))[0]), 0.0)
    n_glob_cur = 1.0 / max(var_mean + s2_between_cur, 1e-9)

    node_rtype, rtype = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind); is_reg = kind == REGION
    EGl = geo.eff_gdna_left; Ml = geo.mass_left
    rho_obs = Ml / np.maximum(EGl, 1e-9)

    # PROPOSED floor: intergenic (rtype 0) + intron (rtype 1) REGIONS, NOT boundaries/exons.
    floor = is_reg & ((node_rtype == 0) | (node_rtype == 1))
    fm = floor & (Ml > 0)
    eff_f = EGl[floor]; mass_f = Ml[floor]
    G_f = float(mass_f.sum()); E_f = max(float(eff_f.sum()), 1e-9)
    rho_floor = (1.0 + G_f) / E_f
    var_mean_f = 1.0 / (1.0 + G_f)
    # population spread of log-density over floor regions with mass>0 (eff-weighted), minus Poisson floor.
    if fm.sum() >= 2:
        lr = np.log(np.maximum(rho_obs[fm], 1e-9)); w = EGl[fm]
        mu = np.sum(w * lr) / np.sum(w)
        s2_raw = np.sum(w * (lr - mu) ** 2) / np.sum(w)
        pois = float(np.mean(1.0 / (rho_obs[fm] * EGl[fm] + 1.0)))
        s2_floor = max(s2_raw - pois, 0.0)
    else:
        s2_floor = 0.0
    n_glob_floor = 1.0 / max(var_mean_f + s2_floor, 1e-9)

    print(f"\n===== {name} kappa={kappa}(fit {kfit:.3f}) cap={capture} =====")
    print(f"  CURRENT  rho_global={rho_g_cur:.4f}  s2_between={s2_between_cur:.4f} var_mean={var_mean:.4f}"
          f"  -> n_glob={n_glob_cur:.3f}")
    print(f"  FLOOR    rho_floor ={rho_floor:.4f}  s2_floor  ={s2_floor:.4f} var_mean={var_mean_f:.4f}"
          f"  -> n_glob={n_glob_floor:.3f}   (floor regions: {int(floor.sum())}, with mass: {int(fm.sum())})")
    rdf = idx.region_df
    print(f"  {'reg':>3} {'span':>12} {'rtype':>10} | {'M/E_g':>8} | "
          f"{'cur_tgt_fg':>10} {'floor_tgt_fg':>12}")
    for n in np.where(is_reg)[0]:
        ri = int(np.asarray(chain.ref_idx)[n])
        rt = {0: 'intergenic', 1: 'intron', 2: 'exon'}[int(node_rtype[n])]
        M = Ml[n]; E = EGl[n]
        cur_t = np.clip(rho_g_cur * E / max(M, 1e-9), 1e-9, 1.0)
        flo_t = np.clip(rho_floor * E / max(M, 1e-9), 1e-9, 1.0)
        print(f"  {ri:>3} {f'{rdf.start.iloc[ri]}-{rdf.end.iloc[ri]}':>12} {rt:>10} | "
              f"{rho_obs[n]:>8.4f} | {cur_t:>10.3f} {flo_t:>12.3f}")


if __name__ == "__main__":
    TA = [("TA", "+", [(1000, 2000), (5000, 10000)])]
    probe("S2_unstr_nocap", TA, kappa=0.5, capture=False)
    probe("S4_unstr_cap", TA, kappa=0.5, capture=True)
