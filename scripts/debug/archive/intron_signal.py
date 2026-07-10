"""Is there an intron-coverage signal distinguishing real nascent from gDNA, under capture?

For each condition, measures the unspliced read DENSITY in intron-class regions (nascent + gDNA there)
vs exon-class. The nascent signal = (nrna_rnd intron density) - (nrna_none intron density) at matched
gDNA. If capture collapses the intron signal, nascent vs gDNA is unidentifiable under capture (no prior
can recover it); if intron density still rises with real nascent under capture, a prior keyed on intron
coverage can break the eff-len tradeoff.
"""
from dataclasses import replace as _dc_replace
import numpy as np
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import RegionType, coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
idx = TranscriptIndex.load(f"{S}/rigel_index")
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
coarse = coarse_type_array(np.asarray(ra.signature))
L = np.asarray(ra.region_size_bp, float)
intron = coarse == int(RegionType.INTRON)
exon = coarse == int(RegionType.EXON)

print(f"{'condition':50s} {'intron_dens':>11} {'exon_dens':>10} {'intron/exon':>11}")
for C in ["gdna_none_ss_0.99_nrna_none_capture_off","gdna_none_ss_0.99_nrna_rnd_capture_off",
          "gdna_none_ss_0.99_nrna_none_capture_on","gdna_none_ss_0.99_nrna_rnd_capture_on"]:
    BAM=f"{S}/{C}/sim_oracle.bam"
    cfg=PipelineConfig(); scan=_dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(BAM))
    stats,sm,flm,buf,pl=scan_and_buffer(BAM,idx,scan); buf.cleanup()
    sub=CalibrationSubstrate.from_payload(pl, ra)
    u = sub.contained.n_unspliced_pos.astype(float)+sub.contained.n_unspliced_neg.astype(float)
    idens = u[intron].sum()/L[intron].sum()
    edens = u[exon].sum()/L[exon].sum()
    print(f"{C:50s} {idens:>11.4f} {edens:>10.4f} {idens/max(edens,1e-9):>11.4f}")
