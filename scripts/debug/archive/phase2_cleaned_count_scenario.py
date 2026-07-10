"""Phase-2 validation: the count field on strand-CLEANED counts (the AMBIG fix + robustness).

Isolates the Phase-2 deliverable — ``node_gdna_density(gdna_counts=...)`` fed strand-cleaned counts
(``strand_deconv.cleaned_gdna_count``) — on the toy AMBIG scenario, **without** the Phase-3 combine
(``calibrate`` still runs on raw counts). For each true strand specificity we:

  1. scan the toy genome (gDNA=0 + nascent — the AMBIG over-call condition; oracle gDNA fraction = 0);
  2. fit κ, strand-deconvolve every node at the fitted κ, and clean each view's count
     (``cleaned = (w·g_strand + (1−w)·1)·raw``, ``w = info/(info+I₀)``);
  3. build the count field twice — RAW (``gdna_counts=None``) vs CLEANED — and report the AMBIG
     region's ``count_gdna_frac`` against the oracle (0).

Expectation (the two things Phase 2 must show):
  - **AMBIG fix:** at high ss the cleaned AMBIG ``count_gdna_frac`` collapses toward 0 (the raw field
    over-states it because the imputation anchors carry nascent the count clue can't see).
  - **Robustness / graceful degradation:** as the true ss → ½, the fitted κ → ½, info → 0, and the
    cleaning degrades *continuously* to a no-op (cleaned → raw) — no cliff, no garbage at κ≈½.

Usage:  python scripts/debug/phase2_cleaned_count_scenario.py            # ss sweep
        python scripts/debug/phase2_cleaned_count_scenario.py 0.99       # one ss
"""
import dataclasses
import sys
import tempfile

import numpy as np

from rigel.calibration.calibrate import CalibrationSubstrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_deconv import cleaned_gdna_count, strand_deconvolve
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig
from rigel.splice import SpliceType

# reuse the toy genome + oracle from the Phase-3 scenario (identical geometry)
from phase3_ambig_scenario import make_scenario, oracle_region_unspliced  # noqa: E402

SC_LABEL = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}
INFO_SCALE = 1.0  # I₀ — the strand-information half-trust scale (Phase-3 combine weight scale)


def _raw_counts(substrate):
    return tuple(
        (v.n_unspliced_pos + v.n_unspliced_neg).astype(np.float64)
        for v in (substrate.contained, substrate.left, substrate.right)
    )


def run_condition(strand_specificity, work, nrna_abund=30):
    sc = make_scenario(work)
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(
            frag_mean=250, frag_std=50, frag_min=80, frag_max=600, read_length=100,
            strand_specificity=strand_specificity, seed=7,
        ),
        gdna_config=None, nrna_abundance=float(nrna_abund),  # gDNA=0 ⇒ oracle gDNA fraction = 0
    )
    idx, bam = res.index, str(res.bam_path)
    df = idx.region_df
    starts = {ref: gg["start"].to_numpy() for ref, gg in df.groupby("ref_name")}
    ids = {ref: gg["region_id"].to_numpy() for ref, gg in df.groupby("ref_name")}
    R = len(df)

    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size,
    )

    substrate = CalibrationSubstrate.from_payload(pl, ra)
    region_eff_len = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    fl_mean = boundary_eff_length(fl.gdna_pmf)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)

    # Strand-deconvolve every node at the fitted κ, then clean each view's count.
    contained, left, right = strand_deconvolve(
        substrate, ra, rna_sense_frac=kappa,
        gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0, n_grid=cfg.calibration.n_grid,
    )
    c_raw, l_raw, r_raw = _raw_counts(substrate)
    cleaned = (
        cleaned_gdna_count(contained, c_raw, INFO_SCALE),
        cleaned_gdna_count(left, l_raw, INFO_SCALE),
        cleaned_gdna_count(right, r_raw, INFO_SCALE),
    )

    nd_raw = node_gdna_density(substrate, ra, region_eff_len, fl_mean)
    nd_clean = node_gdna_density(
        substrate, ra, region_eff_len, fl_mean, gdna_counts=cleaned
    )

    sc_arr = np.asarray(ra.strand_class)
    start, end = df["start"].to_numpy(), df["end"].to_numpy()
    og, on, om = oracle_region_unspliced(bam, starts, ids, R)
    ambig = np.flatnonzero(sc_arr == TS_AMBIG)
    show = sorted({k for a in ambig for k in (a - 1, a, a + 1) if 0 <= k < R})

    print(f"\n===== true ss={strand_specificity:.2f}  (fitted κ={kappa:.3f})  gDNA=0 nRNA={nrna_abund} =====")
    print(f"  {'region':>14} {'cls':>5} {'gdna':>5} {'nrna':>5} {'mat':>5} "
          f"{'orcl_g':>7} {'raw_cf':>7} {'cln_cf':>7}")
    for i in show:
        tot = og[i] + on[i] + om[i]
        orc = og[i] / tot if tot > 0 else float("nan")
        mark = " <AMBIG" if sc_arr[i] == TS_AMBIG else ""
        print(f"  {f'{start[i]}-{end[i]}':>14} {SC_LABEL[sc_arr[i]]:>5} "
              f"{og[i]:>5.0f} {on[i]:>5.0f} {om[i]:>5.0f} "
              f"{orc:>7.3f} {nd_raw.count_gdna_frac[i]:>7.3f} {nd_clean.count_gdna_frac[i]:>7.3f}{mark}")
    sc.cleanup()


if __name__ == "__main__":
    sweep = [0.99, 0.80, 0.60, 0.54, 0.50] if len(sys.argv) < 2 else [float(sys.argv[1])]
    with tempfile.TemporaryDirectory() as work:
        for ss in sweep:
            run_condition(ss, work)
