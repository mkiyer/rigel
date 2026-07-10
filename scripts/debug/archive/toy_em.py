"""Production-faithful toy driver WITH the EM — sim → oracle BAM → calibrate() → quant_from_buffer().

Extends toy_prod.py past calibration into the per-locus EM, so we can measure the END-TO-END gDNA leak
(soft-EM gDNA fraction vs oracle) on tiny hand-specified transcripts — the controlled bed for prototyping
effective-length / prior-assembly strategies. NO hand-rolled EM; this IS the production quant path.

API:  run_em(name, genes, kappa=, n_rna=, gdna_fraction=, capture=, capture_strength=, nascent=)
      genes = [("GA","+",[(1000,1500),(3000,3500),(5000,5500)], 100.0), ...]  (last elem = mature abundance)
Reports: oracle gDNA fraction, calibration prior gDNA fraction, soft-EM gDNA fraction, gDNA eff-len.
"""
from __future__ import annotations

import collections
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig, CalibrationConfig, EMConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.splice import SpliceType

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/e1a42517-d9da-4b58-b12a-0b84f986ddef/scratchpad")


def _write_probes(path, ref, genes):
    with open(path, "w") as f:
        for gid, _strand, exons, *_rest in genes:
            for i, (s, e) in enumerate(exons):
                f.write(f"{ref}\t{s}\t{e}\t{gid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")


def _oracle_gdna_fraction(bam_path):
    """Global oracle gDNA fraction over read1 fragments (from read-name origin)."""
    g = r = 0
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1 or rd.mate_is_unmapped:
            continue
        if abs(rd.template_length) == 0:
            continue
        if "gdna" in rd.query_name.lower():
            g += 1
        else:
            r += 1
    bam.close()
    return g, r


def run_em(name, genes, *, kappa=0.99, n_rna=8000, gdna_fraction=0.75, capture=True,
           capture_strength=20.0, nascent=0.0, genome_length=14000, seed=7, config=None):
    wd = SCRATCH / f"toyem_{name}"
    sc = Scenario(name, genome_length=genome_length, seed=seed, work_dir=wd, ref_name="chr1")
    for gid, strand, exons, *rest in genes:
        ab = float(rest[0]) if rest else 100.0
        sc.add_gene(gid, strand, [{"t_id": gid, "exons": exons, "abundance": ab}])
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
    st, sm, flm, buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, config or CalibrationConfig())
    # calibration prior (deconvolved) gDNA fraction over contained+sides
    gcont = np.asarray(cal.mass_gdna_contained).sum()
    gsides = np.asarray(cal.mass_gdna_left).sum() + np.asarray(cal.mass_gdna_right).sum()
    rcont = np.asarray(cal.mass_rna_contained).sum()
    rsides = np.asarray(cal.mass_rna_left).sum() + np.asarray(cal.mass_rna_right).sum()
    calib_gdna = gcont + gsides
    calib_frac = calib_gdna / max(calib_gdna + rcont + rsides, 1e-9)
    # EM
    est, _ = quant_from_buffer(buf, idx, sm, fl, ra, st, cal, pl, em_config=EMConfig(), emit_locus_stats=True)
    lr = est.locus_results
    ldf = est.get_loci_df(idx)  # mrna, nrna, gdna, total (same schema as the benchmark loci.feather)
    em_gdna = float(ldf["gdna"].sum())
    em_rna = float(ldf["mrna"].sum() + ldf["nrna"].sum())
    em_frac = em_gdna / max(em_gdna + em_rna, 1e-9)
    gdna_efflens = [float(x["gdna_eff_len_em"]) for x in lr]
    # oracle
    og, orr = _oracle_gdna_fraction(result.bam_path)
    orc_frac = og / max(og + orr, 1e-9)
    print(f"\n===== {name}: kappa={kappa} capture={'ON x'+str(capture_strength) if capture else 'OFF'} "
          f"gdna_frac={gdna_fraction} nascent={nascent} =====")
    print(f"  oracle gDNA count={og} RNA={orr}  -> oracle gDNA frac = {orc_frac:.3f}")
    print(f"  calibration prior gDNA frac = {calib_frac:.3f}  (gDNA mass {calib_gdna:.0f})")
    print(f"  SOFT-EM gDNA frac           = {em_frac:.3f}  (gDNA {em_gdna:.0f} / RNA {em_rna:.0f})")
    print(f"  gDNA eff-len per locus      = {[round(e,0) for e in gdna_efflens]}")
    for x in lr:
        print(f"    locus {x['locus_id']}: n_em={x['n_em_fragments']} n_t={x['n_transcripts']} "
              f"gdna_prior={x['gdna_prior_count']:.0f} rna_prior={x['rna_prior_count']:.0f} "
              f"| EM gdna={x['gdna']:.0f} rna={x['rna_total']:.0f} | gdna_efflen={x['gdna_eff_len_em']:.0f}")
    # transcript eff-lens (EM)
    tel = est._t_eff_len_em if hasattr(est, '_t_eff_len_em') else None
    if tel is not None:
        import numpy as _np
        print(f"    transcript em_eff_len: min={_np.min(tel):.0f} median={_np.median(tel):.0f} max={_np.max(tel):.0f}")
    print(f"  ---> LEAK (oracle - soft-EM) = {orc_frac - em_frac:+.3f}")
    return dict(oracle=orc_frac, calib=calib_frac, soft_em=em_frac, gdna_efflen=gdna_efflens, n_loci=len(lr))


if __name__ == "__main__":
    import os
    # 3-exon gene, capture ON, stranded, gDNA 0.75 + expressed mature — should reproduce the eff-len leak
    GENES = [("GA", "+", [(1000, 1500), (3000, 3500), (5000, 5500)], 100.0)]
    scale = os.environ.get("RIGEL_DBG_GDNA_EFFLEN_SCALE", "1.0")
    print(f"[gdna_efflen_scale={scale}]")
    run_em("multi_exon_capture", GENES, kappa=0.99, n_rna=8000, gdna_fraction=0.75, capture=True, capture_strength=20.0)
