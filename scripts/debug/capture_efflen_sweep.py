"""Capture binding-energy × DNA:RNA diagnostic sweep for the effective-length leak.

Hypothesis: hybrid-capture enrichment is the ROOT of the gDNA->RNA leak. Under capture-OFF the tool
should be pristine (eff-lengths correct, ~zero leak). As the binding-energy knob rises, gDNA concentrates
on the probes (exons); at the extreme, ALL fragments sit over probes and gDNA occupies only the exon
footprint — so the gDNA effective length should CONTRACT to the transcript (exon) footprint. If it does
not, the EM under-weights gDNA at exons and it leaks. We test the FULL tool end-to-end (scan -> calibrate
-> quant), expose every component's effective length (mature / nascent / gDNA, base vs capture-contracted),
and the leak split by pool (-> mature, -> nascent), across a grid of (binding_per_base, DNA:RNA ratio).

    OMP_NUM_THREADS=4 python scripts/debug/capture_efflen_sweep.py
"""
from __future__ import annotations
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

# Controlled gene set: multi-exon at HIGH and LOW abundance (the abundance-bias axis) + a single-exon gene.
# exons doubled as capture probes (see _write_probes), so capture enriches exactly the exon footprint.
GENES = [
    ("HI",  "+", [(1000, 1500), (3000, 3500), (5000, 5500)], 2000.0),   # 3-exon, HIGH abundance
    ("LO",  "+", [(16000, 16500), (18000, 18500), (20000, 20500)], 100.0),  # 3-exon, LOW abundance
    ("SE",  "+", [(31000, 33000)], 500.0),                              # single-exon
]
GENOME = 40000


def _write_probes(path, ref, genes):
    with open(path, "w") as f:
        for gid, _st, exons, *_ in genes:
            for i, (s, e) in enumerate(exons):
                f.write(f"{ref}\t{s}\t{e}\t{gid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")


def _oracle(bam):
    g = r = 0
    b = pysam.AlignmentFile(str(bam), "rb")
    for rd in b.fetch(until_eof=True):
        if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1 or rd.mate_is_unmapped:
            continue
        if abs(rd.template_length) == 0:
            continue
        (g := g + 1) if "gdna" in rd.query_name.lower() else (r := r + 1)  # noqa
    b.close()
    return g, r


def run_cell(binding, gdna_ratio, *, kappa=0.99, n_rna=12000, seed=7):
    tag = f"b{binding if binding is not None else 'off'}_g{gdna_ratio}"
    wd = SCRATCH / f"efflen_{tag}"
    sc = Scenario(tag, genome_length=GENOME, seed=seed, work_dir=wd, ref_name="chr1")
    for gid, strand, exons, ab in GENES:
        sc.add_gene(gid, strand, [{"t_id": gid, "exons": exons, "abundance": ab}])
    cap_cfg = None
    if binding is not None:
        wd.mkdir(parents=True, exist_ok=True)
        probes = wd / "probes.bed"
        _write_probes(probes, "chr1", GENES)
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=binding)
    result = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=gdna_ratio,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=120, strand_specificity=kappa, seed=seed),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=250, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=0.0,
    )
    idx = result.index
    st, sm, flm, buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
    est, _ = quant_from_buffer(buf, idx, sm, fl, ra, st, cal, pl, em_config=EMConfig(), emit_locus_stats=True)
    lr = est.locus_results
    ldf = est.get_loci_df(idx)
    og, orr = _oracle(result.bam_path)
    em_g = float(ldf["gdna"].sum()); em_m = float(ldf["mrna"].sum()); em_n = float(ldf["nrna"].sum())
    # per-component eff-lengths (EM). t_df aligns with _t_eff_len_em by row.
    tdf = idx.t_df.reset_index(drop=True)
    is_nrna = tdf["is_nrna"].to_numpy()
    tel = np.asarray(est._t_eff_len_em, float)       # capture-contracted transcript eff-len
    tout = np.asarray(est._t_eff_len_output, float)  # base (FL-marginal) transcript eff-len
    exonic = np.asarray(est._exonic_lengths, float)
    mat = ~is_nrna
    gdna_eff = np.array([float(x["gdna_eff_len_em"]) for x in lr])
    span = np.array([float(x.get("locus_span_bp", np.nan)) for x in lr])
    return dict(
        binding=binding, gdna_ratio=gdna_ratio,
        oracle_g=og, oracle_r=orr, oracle_frac=og / max(og + orr, 1),
        em_g=em_g, em_m=em_m, em_n=em_n, em_frac=em_g / max(em_g + em_m + em_n, 1),
        total_leak=og - em_g, mature_leak=em_m - orr, nascent_leak=em_n,
        mat_eff=float(np.median(tel[mat])), mat_eff_base=float(np.median(tout[mat])),
        nasc_eff=float(np.median(tel[is_nrna])) if is_nrna.any() else float("nan"),
        exon_foot=float(np.median(exonic[mat])),
        gdna_eff=float(np.median(gdna_eff)), span=float(np.nanmedian(span)),
    )


if __name__ == "__main__":
    bindings = [None, 2.0, 10.0, 50.0, 200.0]     # off -> extreme
    ratios = [0.1, 0.5, 1.0, 3.0]                  # DNA:RNA low -> DNA-dominant
    print(f"{'bind':>5} {'DNA:RNA':>7} | {'orcF':>5} {'emF':>5} {'leak':>7} {'->mat':>7} {'->nasc':>7} "
          f"| {'matEff':>7} {'gdnaEff':>7} {'gE/mE':>6} {'gE/span':>7} {'matContr':>8}")
    for b in bindings:
        for gr in ratios:
            try:
                r = run_cell(b, gr)
            except Exception as e:
                print(f"{str(b):>5} {gr:>7} | FAILED: {e}")
                continue
            print(f"{str(b):>5} {gr:>7.2g} | {r['oracle_frac']:>5.2f} {r['em_frac']:>5.2f} "
                  f"{r['total_leak']:>7.0f} {r['mature_leak']:>7.0f} {r['nascent_leak']:>7.0f} | "
                  f"{r['mat_eff']:>7.0f} {r['gdna_eff']:>7.0f} {r['gdna_eff']/max(r['mat_eff'],1):>6.2f} "
                  f"{r['gdna_eff']/max(r['span'],1):>7.2f} {r['mat_eff']/max(r['mat_eff_base'],1):>8.2f}")
