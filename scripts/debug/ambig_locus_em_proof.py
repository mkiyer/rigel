"""PROOF experiment: is the nascent siphon an EM bias, or calibration error that snowballs?

A controlled ambiguous locus: two overlapping transcripts on OPPOSITE strands, each 2 exons, the SAME
spliced length, splice junctions OFFSET so there is a single intron region that is intronic for BOTH (the
ambiguous intron where gDNA and the two nascent spans compete). NASCENT abundance = 0, so any nascent count
is a pure siphon. We run the per-locus EM with (a) the FITTED calibration and (b) the ORACLE calibration
(true per-node gDNA/RNA masses from the read-name origins) — same fragments, same scoring, only the prior
differs. If the nascent siphons under the ORACLE prior, the EM is intrinsically biased. If it does NOT (and
only siphons under FITTED), the siphon is calibration error — possibly snowballed by the EM from a few
mis-split fragments. Capture OFF and ON both tested.

  OMP_NUM_THREADS=1 python scripts/debug/ambig_locus_em_proof.py
"""
import dataclasses
import os
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig, CalibrationConfig, EMConfig, FragmentScoringConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType
from _metrics import oracle_node_masses, rna_component_breakdown

WD = Path(os.environ["SCRATCH_DIR"]) if "SCRATCH_DIR" in os.environ else Path(
    "/private/tmp/claude-503/-Users-mkiyer-proj-rigel/abe830c0-6786-4484-8f6b-96f6a75b0c35/scratchpad")

# genes spec: list of (gene_id, strand, [(transcript_id, [(exon_start, exon_end), ...]), ...]).
# SINGLE: minimal — 1 isoform/gene, 2 exons, short introns, opposite strands (the clean control).
GENES_SINGLE = [
    ("TA", "+", [("TA.1", [(1000, 4000), (6000, 9000)])]),
    ("TB", "-", [("TB.1", [(1000, 3000), (5000, 9000)])]),
]
# MULTI: flagship-like — 3 isoforms/gene, 3 exons, LONG introns, varied splicing, opposite strands.
GENES_MULTI = [
    ("GA", "+", [("GA.1", [(1000, 4000), (15000, 18000), (28000, 31000)]),
                 ("GA.2", [(1000, 4000), (12000, 14000), (28000, 31000)]),
                 ("GA.3", [(1000, 3000), (28000, 31000)])]),
    ("GB", "-", [("GB.1", [(2000, 5000), (15000, 18000), (27000, 30000)]),
                 ("GB.2", [(2000, 5000), (20000, 22000), (27000, 30000)]),
                 ("GB.3", [(2000, 4000), (27000, 30000)])]),
]
GENES_BY_MODE = {"single": GENES_SINGLE, "multi": GENES_MULTI}
GENOME_LEN = {"single": 11000, "multi": 33000}
N_RNA = int(os.environ.get("N_RNA", "8000"))  # depth (snowball test): raise to flagship per-locus scale
GDNA_FRACTION = 0.5
RNA_FRAG_MEAN = 250
# FL interrogation: the gDNA sim fragment length. Default 250 == RNA (EQUAL FL, no FL signal — the
# bias-ruling-out baseline). Set e.g. 380 to mimic the real flagship's gDNA>RNA length separation.
GDNA_FRAG_MEAN = int(os.environ.get("GDNA_FRAG_MEAN", "250"))


def _pmf_mean(pmf):
    p = np.asarray(pmf, dtype=np.float64)
    return float(np.dot(np.arange(p.size), p) / max(p.sum(), 1e-12))


# oracle_node_masses + rna_component_breakdown live in _metrics.py (single source of truth for the siphon
# metric + oracle calibration — see docs/calibration/siphon_measurement.md).


def truth_counts(bam_path):
    c = {"gdna": 0, "nrna": 0, "mrna": 0}
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        seen = set()
        for rd in f:
            if rd.query_name in seen:
                continue
            seen.add(rd.query_name)
            c[parse_origin(rd.query_name).kind] = c.get(parse_origin(rd.query_name).kind, 0) + 1
    return c


def run_one(name, capture, genes_mode):
    genes = GENES_BY_MODE[genes_mode]
    wd = WD / f"ambigem_{name}"
    sc = Scenario(name, genome_length=GENOME_LEN[genes_mode], seed=7, work_dir=wd, ref_name="chr1")
    for gid, strand, isoforms in genes:
        sc.add_gene(gid, strand, [{"t_id": tid, "exons": exons, "abundance": 100.0}
                                  for tid, exons in isoforms])
    cap_cfg = None
    if capture:
        wd.mkdir(parents=True, exist_ok=True)
        probes = wd / "probes.bed"
        with open(probes, "w") as fh:
            for gid, _s, isoforms in genes:
                for tid, exons in isoforms:
                    for i, (s, e) in enumerate(exons):
                        fh.write(f"chr1\t{s}\t{e}\t{tid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=20.0)
    result = sc.build_oracle(
        n_rna_fragments=N_RNA, gdna_fraction=GDNA_FRACTION,
        sim_config=ReadSimConfig(frag_mean=RNA_FRAG_MEAN, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=120, strand_specificity=0.99, seed=7),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=GDNA_FRAG_MEAN, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=0.0)  # NASCENT = 0
    idx = result.index
    tc = truth_counts(result.bam_path)
    n_nrna_rows = int(idx.t_df["is_nrna"].to_numpy(bool).sum())

    def quant(lever):
        _st, sm, flm, buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
        sm.finalize()
        fl = build_fl_models(global_counts=flm.global_model.counts,
                             rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                             gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
        gdna_fl_mean, rna_fl_mean = _pmf_mean(fl.gdna_pmf), _pmf_mean(fl.rna_pmf)
        ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
        cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
        if lever == "oracle":  # perfect per-node masses (isolates the EM/scoring from calibration error)
            cal = dataclasses.replace(cal, **oracle_node_masses(result.bam_path, ra, idx))
        if lever == "fl_neutral":  # force the scorer's gDNA FL == RNA FL ⇒ ZERO FL gDNA-vs-RNA discrimination
            fl = dataclasses.replace(fl, gdna_pmf=fl.rna_pmf.copy())
        est, _ = quant_from_buffer(buf, idx, sm, fl, ra, _st, cal, pl,
                                   em_config=EMConfig(), scoring=FragmentScoringConfig())
        # SIPHON = EM mass on the SYNTHETIC shadows (see _metrics / docs/calibration/siphon_measurement.md).
        siphon, _single_exon, mature_em = rna_component_breakdown(est, idx)
        return siphon, mature_em, gdna_fl_mean, rna_fl_mean

    print(f"\n===== {name}  (genes={genes_mode}  capture={'ON' if capture else 'OFF'}  "
          f"sim FL: gDNA={GDNA_FRAG_MEAN} RNA={RNA_FRAG_MEAN}) =====")
    print(f"  TRUTH: gdna={tc['gdna']}  mrna={tc['mrna']}  nrna={tc['nrna']} (=0)  |  {n_nrna_rows} nascent EM rows")
    print(f"  {'lever':12} {'nascent_EM(siphon)':>18} {'mature_EM':>10} {'est.FL gDNA/RNA':>18}")
    for lever in ("fitted", "fl_neutral", "oracle"):
        siphon, mature, gfm, rfm = quant(lever)
        print(f"  {lever:12} {siphon:>18,.1f} {mature:>10,.1f} {f'{gfm:.1f}/{rfm:.1f}':>18}")


if __name__ == "__main__":
    gm = os.environ.get("MODE", "multi")
    for cap in (False, True):
        run_one(f"{gm}_{'cap' if cap else 'nocap'}", cap, gm)
