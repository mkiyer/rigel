"""PROVE the isoform-multiplicity siphon on the REAL mechanism: multi-exon isoforms (each spawns a nascent
shadow) + gDNA, ORACLE-calibrated so it is PURELY the EM. Vary the number of isoforms N; measure how much
true gDNA is siphoned into RNA (mostly the nascent shadows). Single-exon can't show this (no shadows) — the
real leak went to shadows, so we need multi-exon.

Design per N: one +strand gene with N 2-exon isoforms sharing exon1 but with DISTINCT exon2 (⇒ N distinct
unspliced spans ⇒ N nascent shadows). Equal RNA abundance split across isoforms (total fixed). gDNA 3:1 with
capture ON (concentrates gDNA onto the exons, into the locus). Calibration replaced by the validated oracle.
"""
import os
import sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
from pathlib import Path
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from oracle import OracleTruth
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType
from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig

WD = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "em_mult_sim"
cfg = PipelineConfig()


def build(N, total_abund=120.0):
    wd = WD / f"n{N}"
    sc = Scenario(f"mult{N}", genome_length=60000, seed=13, work_dir=wd, ref_name="chr1")
    isos = []
    for i in range(N):
        e2 = 3000 + i * 800
        isos.append({"t_id": f"iso{i}", "exons": [(1000, 1400), (e2, e2 + 400)],
                     "abundance": total_abund / N})
    sc.add_gene("G", "+", isos)
    probes = wd / "probes.bed"
    wd.mkdir(parents=True, exist_ok=True)
    with open(probes, "w") as fh:
        for i in range(N):
            for (s, e) in [(1000, 1400), (3000 + i * 800, 3000 + i * 800 + 400)]:
                fh.write(f"chr1\t{s}\t{e}\tp\t0\t+\t{s}\t{e}\t0\t1\t{e-s}\t0\n")
    res = sc.build_oracle(
        n_rna_fragments=4000, gdna_fraction=3.0, nrna_abundance=0.0,
        sim_config=ReadSimConfig(frag_mean=180, frag_std=30, frag_min=80, frag_max=400,
                                 read_length=90, strand_specificity=0.99, seed=13),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=200, frag_std=40),
        capture_config=CaptureConfig(probes=str(probes), binding_per_base=20.0))
    return res


def run(res, N):
    bam = str(res.bam_path)
    index = res.index
    ra = RegionArrays.from_index(index)
    orc = OracleTruth.from_bam(bam, index, cfg, WD, f"n{N}")
    override = orc.override_masses(ra)
    # true observed gDNA fragment count from the oracle read names (nrna_abundance=0 ⇒ true nascent=0,
    # so nrna_em_count IS the pure siphon — the direct proof).
    import pysam
    from rigel.sim.read_name import parse_origin
    true_g = 0
    with pysam.AlignmentFile(bam, "rb") as f:
        seen = set()
        for rd in f:
            if rd.query_name in seen:
                continue
            seen.add(rd.query_name)
            if parse_origin(rd.query_name).kind == "gdna":
                true_g += 1
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, sc)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    cal = dataclasses.replace(cal, **override)  # ORACLE calibration
    est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                               em_config=cfg.em, scoring=cfg.scoring)
    g = float(est.gdna_em_count) + float(stats.n_intergenic)
    return dict(true_g=true_g, assigned_g=g, siphon=float(est.nrna_em_count),
                mature=float(est.get_counts_df(index)["count"].sum()))


print(f"{'N_iso':>6} {'true_gdna':>10} {'assigned_g':>11} {'gdna_leak':>10} {'siphon(nrna)':>13} {'mature':>10}")
for N in [1, 2, 4, 8, 16, 32]:
    res = build(N)
    r = run(res, N)
    print(f"{N:>6} {r['true_g']:>10.0f} {r['assigned_g']:>11.0f} {r['true_g']-r['assigned_g']:>10.0f} "
          f"{r['siphon']:>13.1f} {r['mature']:>10.0f}", flush=True)
