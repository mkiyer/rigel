"""Oracle diagnostic: does fixing the gDNA FL model and/or the gDNA eff_len shrink the leak?

Each config gets a FRESH scan + calibrate + per-locus EM (the buffer is consumed by one quant
pass, so we re-scan). Calibration is deterministic/identical; only the EM input is changed, to
isolate the EM stage:
  (a) baseline,
  (b) FL test: shift the SCORING gDNA FL pmf to the on-target mean (~343bp) instead of off-target ~386,
  (c) eff_len test: scale the per-locus gDNA eff_len ×0.8 (remove the ~25% over-estimate),
  (d) both.
"""

import sys
import numpy as np
from dataclasses import replace as _dc_replace

import rigel.calibration.priors as priors_mod
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

import os
SUITE = os.environ.get("DIAG_SUITE", "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb")
COND = os.environ.get("DIAG_COND", "gdna_gdna400_ss_0.99_nrna_none_capture_on")
BAM = f"{SUITE}/{COND}/sim_oracle.bam"
import json
_ts=json.load(open(f"{SUITE}/{COND}/truth_summary.json"))
_oc=_ts.get("origin_counts",_ts)
TRUE=float(_oc.get("gdna", _oc.get("gDNA", 4_000_000)))
index = TranscriptIndex.load(f"{SUITE}/rigel_index")


def full_run(shift_fl_to=None, eff_scale=None, gdna_pool_idxs=None, label=""):
    cfg = PipelineConfig()
    scan = _dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(BAM))
    stats, sm, flm, buffer, payload = scan_and_buffer(BAM, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size,
    )
    bins = np.arange(fl.gdna_pmf.size, dtype=np.float64)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    if gdna_pool_idxs is not None:  # build scoring gDNA FL from specific accumulator pools
        from rigel.calibration.fl import _aligned, _normalized
        fp = np.asarray(payload.fl_pool_mass, dtype=np.float64)
        raw = fp[list(gdna_pool_idxs)].sum(axis=0)
        pmf = _normalized(_aligned(raw, fl.max_size))
        mn = float(np.dot(np.arange(pmf.size), pmf))
        print(f"  [scoring gDNA FL from pools {gdna_pool_idxs} → mean {mn:.1f}]", flush=True)
        fl = _dc_replace(fl, gdna_pmf=pmf)
    if shift_fl_to == "rna":      # neutral: gDNA FL == RNA FL ⇒ FL gives ZERO gDNA-vs-RNA discrimination
        fl = _dc_replace(fl, gdna_pmf=fl.rna_pmf.copy())
    elif shift_fl_to == "global":  # neutral: gDNA FL == global FL
        fl = _dc_replace(fl, gdna_pmf=fl.global_pmf.copy())
    elif shift_fl_to is not None:
        mean = float(np.dot(bins, fl.gdna_pmf))
        shift = int(round(mean - shift_fl_to))
        g = fl.gdna_pmf; s = np.zeros_like(g)
        if shift > 0:
            s[:-shift] = g[shift:]; s[-1] += g[:shift].sum()
        else:
            s = g.copy()
        fl = _dc_replace(fl, gdna_pmf=s / s.sum())
    orig = priors_mod.assemble_priors
    if eff_scale is not None:
        priors_mod.assemble_priors = lambda c, r, m: _dc_replace(
            orig(c, r, m), gdna_eff_len=np.maximum(orig(c, r, m).gdna_eff_len * eff_scale, 1.0))
    try:
        est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                                   em_config=cfg.em, scoring=cfg.scoring)
    finally:
        priors_mod.assemble_priors = orig
        buffer.cleanup()
    g = est.gdna_em_count
    print(f"  {label:38s} gDNA_em={g:11.0f}  soft leak={100*(TRUE-g)/TRUE:6.2f}%", flush=True)


mode = sys.argv[1] if len(sys.argv) > 1 else "all"
print(f"[suite={SUITE.split('/')[-1]} cond={COND} TRUE_gdna={TRUE:.0f}]", flush=True)
if mode in ("baseline",):
    full_run(label="(a) baseline (current off-target FL)")
if mode in ("rna", "neutral"):
    full_run(shift_fl_to="rna", label="(e) gDNA FL = RNA FL (no FL discrim)")
if mode in ("global", "neutral"):
    full_run(shift_fl_to="global", label="(f) gDNA FL = global FL (neutral)")
# pool-sourced scoring gDNA FL (test the contained/boundary split idea)
if mode in ("contained", "pools"):
    full_run(gdna_pool_idxs=[0, 2], label="(g) gDNA FL = contained off-target (intgen+intron)")
if mode in ("boundary", "pools"):
    full_run(gdna_pool_idxs=[1, 3, 5], label="(h) gDNA FL = ALL boundary pools")
if mode in ("exonbound", "pools"):
    full_run(gdna_pool_idxs=[5], label="(i) gDNA FL = exonic_boundary only")
