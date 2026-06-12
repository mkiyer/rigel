"""Phase A validation: mature_density.mature_density vs ORACLE per-strand contained-unspliced mature.

Runs the module on the flagship complex locus (locus 21, gdna 3:1 capture-ON ss0.99) and compares the
predicted per-strand contained-unspliced mature (M_pos / M_neg) to oracle truth (mrna-origin
contained-unspliced read1 fragments, split by transcript strand). Also reports the PAYOFF diagnostic:
the residual fraction (U − M_pos − M_neg)/U vs the oracle gDNA fraction — the number Phase B/C will use
to stop the leak (must move ~0.40 → ~0.56 on the failing regions).

Usage:  python scripts/debug/phaseA_mature_imputation.py [condition]
"""
import dataclasses
import sys

import numpy as np
import pysam

from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.mature_density import mature_density
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"


_SIG_BITS = ((BIT_EXON_POS, "E+"), (BIT_EXON_NEG, "E-"), (BIT_INTRON_POS, "I+"), (BIT_INTRON_NEG, "I-"))


def decode(s):
    return "|".join(tag for bit, tag in _SIG_BITS if s & bit) or "."


def oracle_mature_by_strand(bam_path, starts, ids, R, tid2strand):
    """Per region: contained-unspliced mrna-origin read1 counts split by transcript strand (+,-) and gDNA."""
    mpos, mneg, gdna = np.zeros(R), np.zeros(R), np.zeros(R)
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
            o = parse_origin(r.query_name)
            if o.kind == "gdna":
                gdna[rid] += 1.0
            elif o.kind == "mrna":
                tid = r.query_name.split(":")[0]  # read-name prefix is the origin transcript id
                st = tid2strand.get(tid)
                (mpos if st == "+" else mneg)[rid] += 1.0
    return mpos, mneg, gdna


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
    # strand map: t_df strand is encoded 1=+,2=- (or '+'/'-'); normalize
    st_raw = dict(zip(idx.t_df["t_id"].astype(str), idx.t_df["strand"].astype(str)))
    tid2strand = {k: ("+" if v in ("+", "1") else "-") for k, v in st_raw.items()}

    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, _sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(pl, ra)
    e_mu = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    fl_mean_rna = boundary_eff_length(fl.rna_pmf)

    md = mature_density(sub, ra, e_mu, fl_mean_rna)

    df = idx.region_df
    starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
    ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
    mpos, mneg, gdna = oracle_mature_by_strand(bam, starts, ids, len(df), tid2strand)

    sig = np.asarray(ra.signature)
    rid = df["region_id"].to_numpy()
    s = df["start"].to_numpy()
    e = df["end"].to_numpy()
    refn = df["ref_name"].to_numpy()
    U = (sub.contained.n_unspliced_pos + sub.contained.n_unspliced_neg).astype(float)

    sel = np.flatnonzero((refn == "chr_syn") & (e > 964416) & (s < 1004165))
    print(f"fl_mean_RNA={fl_mean_rna:.1f}  (locus 21)")
    print(f"{'rid':>4}{'start':>8}{'sig':>6}{'U':>7}{'o_m+':>7}{'o_m-':>7}{'o_gd':>7}"
          f"{'M+':>7}{'M-':>7}{'oFrac':>6}{'resid':>6}")
    for i in sel:
        if not (sig[i] & (BIT_EXON_POS | BIT_EXON_NEG)):
            continue
        omat = mpos[i] + mneg[i]
        ogf = gdna[i] / max(gdna[i] + omat, 1e-9)
        resid = (U[i] - md.mature_pos[i] - md.mature_neg[i]) / max(U[i], 1e-9)
        print(f"{rid[i]:>4}{s[i]:>8}{decode(int(sig[i])):>6}{U[i]:>7.0f}{mpos[i]:>7.0f}{mneg[i]:>7.0f}"
              f"{gdna[i]:>7.0f}{md.mature_pos[i]:>7.0f}{md.mature_neg[i]:>7.0f}{ogf:>6.2f}{resid:>6.2f}")
    print("\n  o_m+/o_m- = oracle contained-unspliced mature by strand; M+/M- = predicted.")
    print("  oFrac = oracle gDNA fraction of contained-unspliced; resid = (U-M+-M-)/U (Phase-B target ≈ oFrac).")


if __name__ == "__main__":
    main()
