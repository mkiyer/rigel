"""Geometry probe: measure the spliced-mass ↔ mature-density relationship directly.

For a toy scenario with nascent=0 (so ALL RNA is mature), build the calibration geometry and
compare, per region/boundary:
  * TRUE RNA density in an exon   = (oracle RNA contained count) / E_rna_contained(L)
  * spliced-implied mature density = (boundary one-sided spliced mass) / E_spl(L_exon_flank)
  * boundary two-sided crossing densities (gDNA / RNA) vs truth

This isolates whether E_spl recovers ρ_m correctly and where any length-dependent factor lives.
"""
from __future__ import annotations

import collections
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.substrate import CalibrationSubstrate, BoundarySubstrate
from rigel.calibration.node_chain import build_node_chain, REGION, BOUNDARY
from rigel.calibration.bp_solver import build_node_geometry, build_node_statics
from rigel.calibration.effective_length import (
    region_eff_length, boundary_side_eff_length, spliced_side_eff_length, fl_mean,
)
from rigel.splice import SpliceType

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/e1a42517-d9da-4b58-b12a-0b84f986ddef/scratchpad")


def _write_probes(path, ref, genes):
    with open(path, "w") as f:
        for gid, _strand, exons in genes:
            for i, (s, e) in enumerate(exons):
                f.write(f"{ref}\t{s}\t{e}\t{gid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")


def _oracle_contained(bam_path, ra, ref_names):
    """Per-region oracle contained counts (gDNA, RNA) via mate+TLEN span fully inside the region."""
    starts = ra.start; ends = ra.end; rid = ra.ref_id; R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]),
                     np.array(arr))
    name2ref = {n: i for i, n in enumerate(ref_names)}
    g = np.zeros(R); r = np.zeros(R)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if (rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1
                or rd.mate_is_unmapped):
            continue
        is_g = "gdna" in rd.query_name.lower()
        rf = name2ref.get(bam.get_reference_name(rd.reference_id))
        if rf is None or rf not in by_ref:
            continue
        tl = abs(rd.template_length)
        if tl == 0:
            continue
        lo = min(rd.reference_start, rd.next_reference_start); hi = lo + tl
        s, e, ii = by_ref[rf]; j = np.searchsorted(s, lo, side="right") - 1
        if 0 <= j < R and s[j] <= lo and hi <= e[j]:
            (g if is_g else r)[ii[j]] += 1.0
    bam.close()
    return g, r


def probe(name, genes, *, kappa=0.5, n_rna=4000, gdna_fraction=0.5, capture=True,
          capture_strength=20.0, nascent=0.0, genome_length=12000, seed=7):
    wd = SCRATCH / f"probe_{name}"
    sc = Scenario(name, genome_length=genome_length, seed=seed, work_dir=wd, ref_name="chr1")
    for gid, strand, exons in genes:
        sc.add_gene(gid, strand, [{"t_id": gid, "exons": exons, "abundance": 100.0}])
    cap_cfg = None
    if capture:
        probes = wd / "probes.bed"; wd.mkdir(parents=True, exist_ok=True)
        _write_probes(probes, "chr1", genes)
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=capture_strength)
    result = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=gdna_fraction,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=120, strand_specificity=kappa, seed=seed),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=250, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=nascent,
    )
    idx = result.index
    _st, sm, flm, _buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geo = build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)

    g, r = _oracle_contained(result.bam_path, ra, idx.ref_names)
    L = np.asarray(ra.region_size_bp, float)
    E_r = region_eff_length(L, fl.rna_pmf)
    E_g = region_eff_length(L, fl.gdna_pmf)
    print(f"\n===== {name} kappa={kappa} cap={capture} nascent={nascent} =====")
    print(f"  FL means: rna={fl_mean(fl.rna_pmf):.1f}  gdna={fl_mean(fl.gdna_pmf):.1f}")
    rdf = idx.region_df
    print("  --- REGIONS: true contained densities (RNA=mature since nascent=0) ---")
    print(f"  {'reg':>3} {'span':>12} {'L':>6} | {'gTrue':>6} {'rTrue':>6} | "
          f"{'E_r':>7} {'E_g':>7} | {'rho_RNA':>8} {'rho_g':>8}")
    rho_rna_true = np.full(len(rdf), np.nan)
    for i in range(len(rdf)):
        st = rdf['start'].iloc[i]; en = rdf['end'].iloc[i]
        rr = r[i] / E_r[i] if E_r[i] > 0 else np.nan
        gg = g[i] / E_g[i] if E_g[i] > 0 else np.nan
        rho_rna_true[i] = rr
        print(f"  {i:>3} {f'{st}-{en}':>12} {L[i]:>6.0f} | {g[i]:>6.0f} {r[i]:>6.0f} | "
              f"{E_r[i]:>7.1f} {E_g[i]:>7.1f} | {rr:>8.4f} {gg:>8.4f}")

    print("  --- BOUNDARIES: spliced-implied mature density vs adjacent exon true ρ_RNA ---")
    side_eff_spl = spliced_side_eff_length(fl.rna_pmf, L)
    side_eff_r = boundary_side_eff_length(fl.rna_pmf, L)
    side_eff_g = boundary_side_eff_length(fl.gdna_pmf, L)
    blr = np.asarray(bsub.left_region); brr = np.asarray(bsub.right_region)
    for b in range(bsub.n_boundaries):
        lr = int(blr[b]); rr_ = int(brr[b])
        spl_l = float(bsub.left.mass_spliced[b]); spl_r = float(bsub.right.mass_spliced[b])
        mu_l = float(bsub.left.mass_unspliced[b]); mu_r = float(bsub.right.mass_unspliced[b])
        # implied mature density on each side (one-sided half-triangle eff-len of the flank exon)
        esp_l = side_eff_spl[lr] if lr >= 0 else np.nan
        esp_r = side_eff_spl[rr_] if rr_ >= 0 else np.nan
        rho_m_l = spl_l / esp_l if (lr >= 0 and esp_l > 0) else np.nan
        rho_m_r = spl_r / esp_r if (rr_ >= 0 and esp_r > 0) else np.nan
        # true RNA density of the flank exon (the target the mature should recover)
        tl = rho_rna_true[lr] if lr >= 0 else np.nan
        tr = rho_rna_true[rr_] if rr_ >= 0 else np.nan
        print(f"  B{b} lr={lr} rr={rr_} | splL={spl_l:7.2f} splR={spl_r:7.2f} | "
              f"E_spl(L)={esp_l if lr>=0 else 0:6.1f}/{esp_r if rr_>=0 else 0:6.1f} | "
              f"rho_m_L={rho_m_l:7.4f}(exon true {tl:.4f}) "
              f"rho_m_R={rho_m_r:7.4f}(exon true {tr:.4f})")
    return rho_rna_true


if __name__ == "__main__":
    TA = [("TA", "+", [(1000, 2000), (5000, 10000)])]
    probe("TA_unstr_cap", TA, kappa=0.5, capture=True)
    probe("TA_str_cap", TA, kappa=1.0, capture=True)
