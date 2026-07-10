"""Phase-3 development scenario: an AMBIG node fed by splice carry-over.

A small toy genome (controlled, no large-suite confounds) to develop + validate the Phase-3
strand-resolved nascent removal / carry-over sweep:

  - several standard multi-exon transcripts (both strands, well separated) so the strand model trains
    and there are ample intergenic/intronic count-observable regions;
  - one overlapping opposite-strand pair → an AMBIG exon node with single-strand splice-junction
    neighbours that can carry a per-strand {gDNA, RNA⁺, RNA⁻} split into it:
        TA+ exons (1000,1500),(4000,6000)   TB− exons (5000,7000),(10000,10500)
    overlap 5000-6000 = E+ & E− (AMBIG); 4000-5000 = E+ only (TA), 6000-7000 = E− only (TB).

Sweeps gDNA × nascent and reports, per AMBIG region, the calibration's estimated gDNA fraction of the
contained-unspliced mass vs the ORACLE (gdna / (gdna+nrna+mature)). The 2-term target is
(gdna+nrna)/total; the 3-term target is the pure-gDNA fraction. Run on HEAD (no Phase 3) vs the working
tree (3-term) to see whether the AMBIG node is solved.
"""
import sys
import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.config import PipelineConfig
from rigel.calibration.calibrate import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import TS_AMBIG, TS_POS, TS_NEG, TS_NONE, BIT_EXON_POS, BIT_EXON_NEG
from rigel.splice import SpliceType

EXON = BIT_EXON_POS | BIT_EXON_NEG
SC_LABEL = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}


def make_scenario(work):
    sc = Scenario("p3ambig", genome_length=30000, seed=7, work_dir=work)
    # The overlapping opposite-strand pair → AMBIG node + single-strand carry-over neighbours.
    sc.add_gene("gA", "+", [{"t_id": "TA", "exons": [(1000, 1500), (4000, 6000)], "abundance": 100}])
    sc.add_gene("gB", "-", [{"t_id": "TB", "exons": [(5000, 7000), (10000, 10500)], "abundance": 100}])
    # Standard multi-exon transcripts (both strands) for strand training + intergenic/intronic regions.
    sc.add_gene("s1", "+", [{"t_id": "S1", "exons": [(12000, 12500), (13500, 14000), (15000, 15500)],
                             "abundance": 120}])
    sc.add_gene("s2", "-", [{"t_id": "S2", "exons": [(17000, 17500), (18500, 19000), (20000, 20500)],
                             "abundance": 120}])
    sc.add_gene("s3", "+", [{"t_id": "S3", "exons": [(22000, 22600), (24000, 24500)], "abundance": 90}])
    sc.add_gene("s4", "-", [{"t_id": "S4", "exons": [(26000, 26500), (27500, 28000)], "abundance": 90}])
    return sc


def oracle_region_unspliced(bam_path, starts, ids, R):
    """Per region: (gdna, nrna, mature) counts of CONTAINED-UNSPLICED fragments, from origin names."""
    g, n, m = np.zeros(R), np.zeros(R), np.zeros(R)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue  # spliced — not contained-unspliced
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i])
            nm = r.query_name
            (g if nm.startswith("gdna") else n if nm.startswith("nrna") else m)[rid] += 1.0
    return g, n, m


def run_condition(gdna_abund, nrna_abund, work):
    sc = make_scenario(work)
    gd = GDNAConfig(abundance=gdna_abund, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000) \
        if gdna_abund > 0 else None
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=7),
        gdna_config=gd, nrna_abundance=float(nrna_abund),
    )
    idx = res.index
    bam = str(res.bam_path)
    df = idx.region_df
    starts = {ref: gg["start"].to_numpy() for ref, gg in df.groupby("ref_name")}
    ids = {ref: gg["region_id"].to_numpy() for ref, gg in df.groupby("ref_name")}
    R = len(df)

    import dataclasses
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    result = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)

    sc_arr = np.asarray(ra.strand_class)
    start = df["start"].to_numpy()
    end = df["end"].to_numpy()
    og, on, om = oracle_region_unspliced(bam, starts, ids, R)
    gm = np.asarray(result.mass_gdna_contained)
    rm = np.asarray(result.mass_rna_contained)
    est = np.where(gm + rm > 0, gm / np.maximum(gm + rm, 1e-9), np.nan)

    print(f"\n===== gDNA={gdna_abund}  nRNA={nrna_abund}  (κ={result.rna_sense_frac:.3f}) =====")
    ambig = np.flatnonzero(sc_arr == TS_AMBIG)
    if len(ambig) == 0:
        print("  (no AMBIG region formed!)")
    print(f"  {'region':>26} {'cls':>5} {'gdna':>6} {'nrna':>6} {'mat':>6} "
          f"{'orcl_g/tot':>10} {'2term':>7} {'est':>7}")
    # Show AMBIG regions + their immediate neighbours (the carry-over context)
    show = set()
    for a in ambig:
        for k in (a - 1, a, a + 1):
            if 0 <= k < R:
                show.add(k)
    for i in sorted(show):
        tot = og[i] + on[i] + om[i]
        orc = og[i] / tot if tot > 0 else float("nan")
        twoterm = (og[i] + on[i]) / tot if tot > 0 else float("nan")
        mark = " <AMBIG" if sc_arr[i] == TS_AMBIG else ""
        print(f"  {f'{start[i]}-{end[i]}':>26} {SC_LABEL[sc_arr[i]]:>5} "
              f"{og[i]:>6.0f} {on[i]:>6.0f} {om[i]:>6.0f} "
              f"{orc:>10.3f} {twoterm:>7.3f} {est[i]:>7.3f}{mark}")
    sc.cleanup()


if __name__ == "__main__":
    import tempfile
    sweep = [(0, 0), (20, 0), (100, 0), (0, 30), (100, 30)] if len(sys.argv) < 2 else [
        tuple(int(x) for x in sys.argv[1].split(","))]
    with tempfile.TemporaryDirectory() as work:
        for g, n in sweep:
            run_condition(g, n, work)
