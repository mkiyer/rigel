#!/usr/bin/env python
"""Trace nRNA→gDNA mis-routing on the antisense_intronic scenario (PR 8 dive).

Mirrors tests/scenarios/test_antisense_intronic.py::TestAntisenseIntronicOverlap
(g1(+) 3-exon, g2(−) single-exon inside g1's intron, g_ctrl(+)), runs the full
pipeline, and dumps the per-region E-step strand channel so we can see *why*
intronic nascent reads route to gDNA:

  - the strand model (p_r1_sense → κ_rna, ρ_r_bb, the boundary double-count)
  - per active region: ts_class, n_u, the +/− unspliced flux, the oriented
    k_sense, ω, μ_g, the count vs strand log-Bayes-factors, π_g, m_g/m_d
  - the final per-transcript / nRNA / gDNA counts

Usage:  python scripts/debug/trace_antisense_nrna.py [--nrna 30] [--ss 1.0]
"""

from __future__ import annotations

import argparse
import tempfile

import numpy as np

from rigel.calibration.effective_length import region_eff_length
from rigel.calibration.estep import _llr_count, _llr_strand, _sense_flux
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType

_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}
_TS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--nrna", type=float, default=30.0)
    ap.add_argument("--gdna", type=float, default=0.0)
    ap.add_argument("--ss", type=float, default=1.0)
    ap.add_argument("--nfrag", type=int, default=2000)
    args = ap.parse_args()

    tmp = tempfile.mkdtemp(prefix="anti_nrna_")
    sc = Scenario("anti_intron", genome_length=8000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(1000, 1200), (2000, 2200), (5000, 5600)], "abundance": 100},
    ])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0}])
    sc.add_gene("g_ctrl", "+", [{"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0}])
    cfg = ReadSimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=args.ss, seed=42,
    )
    gdna = None if args.gdna == 0 else GDNAConfig(
        abundance=args.gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000
    )
    result = sc.build_oracle(
        n_fragments=args.nfrag, sim_config=cfg, gdna_config=gdna, nrna_abundance=args.nrna
    )
    index = result.index
    pr = run_pipeline(
        result.bam_path, index,
        config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto")),
    )
    cal = pr.calibration
    sm = pr.strand_models

    print("\n=== strand model (from spliced channel) ===")
    print(f"  p_r1_sense={sm.p_r1_sense:.4g}  strand_specificity={sm.strand_specificity:.4f}  "
          f"read1_sense={sm.read1_sense}  n_obs={sm.n_observations}")
    es = sm.exonic_spliced
    print(f"  2x2: pos_pos={es.pos_pos} pos_neg={es.pos_neg} neg_pos={es.neg_pos} neg_neg={es.neg_neg}")
    print("\n=== calibration hyperparameters ===")
    print(f"  iters={cal.n_iterations} converged={cal.converged}")
    print(f"  rho_0={cal.rho_0:.5g}  exposure_dispersion={cal.exposure_dispersion:.4g}")
    print(f"  kappa_rna={cal.kappa_rna:.4g}  rho_r_bb={cal.rho_r_bb:.4g}  rho_d_bb={cal.rho_d_bb:.4g}")

    # Re-run a single E-step with the converged params to expose the channels.
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pr.calibration_payload, ra)
    rtype = coarse_type_array(ra.signature)
    cv = sub.contained
    nu = cv.n_unspliced
    ns = cv.n_spliced
    l_phys = sub.L_eff
    flm = pr.frag_length_models
    gdna_fl_pmf = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pr.calibration_payload),
        max_size=flm.max_size,
    ).gdna_pmf
    region_eff_len = region_eff_length(l_phys, gdna_fl_pmf)
    omega = cal.omega
    mu_g = omega * cal.rho_0 * region_eff_len
    m_d_prev = cal.mass_d_contained
    k_sense = _sense_flux(cv, sub.ts_class)
    llr_c = _llr_count(nu.astype(np.float64), mu_g, m_d_prev, 1.0 / cal.exposure_dispersion)
    llr_s = _llr_strand(k_sense, nu, sub.ts_class, kappa_rna=cal.kappa_rna,
                        rho_r_bb=cal.rho_r_bb, rho_d_bb=cal.rho_d_bb)

    m_g = cal.mass_g_contained + cal.mass_g_left + cal.mass_g_right
    m_d = cal.mass_d_contained + cal.mass_d_left + cal.mass_d_right
    tot = m_g + m_d

    print("\n=== per-region contained channels (active regions) ===")
    print("  idx ts    rtype     start  end   n_u  n+   n-   ksen n_s   omega   mu_g  llr_c  llr_s   m_g   m_d  gfrac")
    for i in range(ra.n_regions):
        if nu[i] == 0 and ns[i] == 0 and tot[i] < 1e-9:
            continue
        gfrac = m_g[i] / tot[i] if tot[i] > 1e-12 else float("nan")
        print(f"  {i:3d} {_TS.get(int(sub.ts_class[i]),'?'):>5s} {_RTYPE.get(int(rtype[i]),'?'):>9s}"
              f" {int(ra.start[i]):6d}{int(ra.end[i]):6d} {int(nu[i]):4d} {int(cv.n_unspliced_pos[i]):4d}"
              f" {int(cv.n_unspliced_neg[i]):4d} {int(k_sense[i]):4d} {int(ns[i]):3d} {omega[i]:7.3f}"
              f" {mu_g[i]:6.2f} {llr_c[i]:+.2f} {llr_s[i]:+.2f} {m_g[i]:6.1f}{m_d[i]:6.1f} {gfrac:5.2f}")

    print("\n=== final counts (expected vs observed) ===")
    counts = pr.estimator.get_counts_df(index)
    for _, row in counts.iterrows():
        if row["transcript_id"] in ("t1", "t2", "t_ctrl"):
            print(f"  {row['transcript_id']:8s} count={row['count']:.1f}")
    print(f"  gdna_em={pr.estimator.gdna_em_count:.1f}  nRNA(synth)={_nrna_total(pr, index):.1f}")


def _nrna_total(pr, index) -> float:
    counts = pr.estimator.get_counts_df(index)
    tdf = index.t_df
    syn = tdf["is_nrna"] if "is_nrna" in tdf.columns else tdf.get("is_synthetic")
    if syn is None:
        return float("nan")
    ids = set(tdf.loc[syn.astype(bool), "transcript_id"]) if syn is not None else set()
    return float(counts[counts["transcript_id"].isin(ids)]["count"].sum())


if __name__ == "__main__":
    main()
