"""Phase A / Issue D: does (spliced-flux / fl_mean_RNA)·E_mu recover the contained-unspliced MATURE?

Clean controlled scenario — ONE multi-exon + gene (internal exons of varying length: short/long/medium),
NO gDNA, NO nascent, high ss. Then contained-unspliced fragments == mature only, so we can validate the
mature-density imputation geometry directly:

    m̂  = mean over bounding eligible junctions of (spliced_flux_side / fl_mean_RNA)
    M̂  = m̂ · E_mu,    E_mu = region_eff_length(L_region, rna_pmf)
    vs  actual contained unspliced count per region (= true mature, since no gDNA/nascent)

If M̂ ≈ actual across exon sizes, the geometry + the accumulator spliced-flux semantics are sound (Issue D
resolved). If a systematic bias appears (esp. vs exon length → an E_ms/E_mu FL-geometry error), quantify it
and propose the eff-length correction.

Usage:  python scripts/debug/phaseA_issueD_mature_density.py
"""
import dataclasses

import numpy as np

from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS,
)
from rigel.calibration.splice_junction import splice_junction_eligibility
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType


_SIG_BITS = ((BIT_EXON_POS, "E+"), (BIT_EXON_NEG, "E-"), (BIT_INTRON_POS, "I+"), (BIT_INTRON_NEG, "I-"))


def decode(s):
    return "|".join(tag for bit, tag in _SIG_BITS if s & bit) or "."


def main():
    import tempfile
    import pathlib
    work = pathlib.Path(tempfile.mkdtemp())
    sc = Scenario("issueD", genome_length=40000, seed=11, work_dir=work)
    # ONE + gene, internal exons of varying length (the geometry stress is exon length vs RNA FL ~250).
    sc.add_gene("gP", "+", [{"t_id": "TP", "exons": [
        (2000, 2500),    # e1 terminal 500
        (4000, 4150),    # e2 SHORT internal 150 (< FL)
        (6000, 7500),    # e3 LONG internal 1500
        (9000, 9500),    # e4 internal 500
        (11000, 11400),  # e5 internal 400
        (13000, 13500),  # e6 terminal 500
    ], "abundance": 300}])
    # a − gene so the strand model + balance train (needs spliced reads on both strands).
    sc.add_gene("gN", "-", [{"t_id": "TN", "exons": [(20000, 20500), (22000, 22500), (24000, 24500)],
                             "abundance": 120}])
    res = sc.build_oracle(
        n_fragments=20000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=None, nrna_abundance=0.0,
    )
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, _sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)

    fl_mean_rna = boundary_eff_length(fl.rna_pmf)
    E_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    sig = np.asarray(ra.signature)
    ref_id = np.asarray(ra.ref_id)
    R = len(sig)

    df = idx.region_df
    rid = df["region_id"].to_numpy()
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()

    L = sub.left.n_spliced_sense.astype(float) + sub.left.n_spliced_antisense.astype(float)
    Rt = sub.right.n_spliced_sense.astype(float) + sub.right.n_spliced_antisense.astype(float)
    contained_unspl = (sub.contained.n_unspliced_pos + sub.contained.n_unspliced_neg).astype(float)

    # per-region: average the bounding ELIGIBLE-junction flux densities (density_model-style mean),
    # where the junction names this region as the exon side.
    print(f"fl_mean_RNA={fl_mean_rna:.1f}")
    print(f"{'rid':>4}{'ref':>6}{'start':>7}{'end':>7}{'sig':>6}{'L_bp':>6}{'Lspl':>6}{'Rspl':>6}"
          f"{'E_mu':>7}{'m_hat':>7}{'M_pred':>8}{'actual':>8}{'ratio':>6}")
    for i in range(R):
        if not (sig[i] & (BIT_EXON_POS | BIT_EXON_NEG)):
            continue  # only exon regions carry mature
        anchors = []
        if i > 0 and ref_id[i] == ref_id[i - 1]:
            el = splice_junction_eligibility(int(sig[i - 1]), int(sig[i]))
            if el is not None and el.exon_side == "R":
                anchors.append(L[i] / fl_mean_rna)
        if i < R - 1 and ref_id[i] == ref_id[i + 1]:
            el = splice_junction_eligibility(int(sig[i]), int(sig[i + 1]))
            if el is not None and el.exon_side == "L":
                anchors.append(Rt[i] / fl_mean_rna)
        if not anchors:
            continue
        m_hat = float(np.mean(anchors))
        M_pred = m_hat * E_mu[i]
        actual = contained_unspl[i]
        ratio = M_pred / actual if actual > 0 else float("nan")
        print(f"{rid[i]:>4}{refn[i]:>6}{s[i]:>7}{e[i]:>7}{decode(int(sig[i])):>6}{e[i]-s[i]:>6}"
              f"{L[i]:>6.0f}{Rt[i]:>6.0f}{E_mu[i]:>7.0f}{m_hat:>7.2f}{M_pred:>8.0f}{actual:>8.0f}{ratio:>6.2f}")
    sc.cleanup()
    print("\nratio = M_pred/actual; want ~1.0 across exon sizes. Systematic deviation vs exon length"
          " ⇒ E_ms/E_mu FL-geometry correction needed.")


if __name__ == "__main__":
    main()
