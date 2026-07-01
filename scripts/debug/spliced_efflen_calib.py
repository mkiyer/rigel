"""Determine E_spl_effective empirically: is the spliced→mature bias a constant or length-dependent?

Build genes with an INTERNAL exon (both junctions real) of varying length, very high depth (low
noise), nascent=0 (RNA=mature). For the internal exon measure:
    E_spl_eff(donor) = spliced_mass_donor_side / rho_RNA_true
    E_spl_eff(accept)= spliced_mass_accept_side / rho_RNA_true
and compare to the coded half-triangle E_spl = fl_mean/2.  Also report the contained-vs-spliced
internal consistency (within-exon mature count / E_r_contained == rho_RNA_true by construction).
"""
from __future__ import annotations

import collections
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.substrate import CalibrationSubstrate, BoundarySubstrate
from rigel.calibration.effective_length import (
    region_eff_length, spliced_side_eff_length, boundary_side_eff_length, fl_mean,
)
from rigel.splice import SpliceType

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/e1a42517-d9da-4b58-b12a-0b84f986ddef/scratchpad")


def _oracle_contained(bam_path, ra, ref_names):
    starts = ra.start; ends = ra.end; rid = ra.ref_id; R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]),
                     np.array(arr))
    name2ref = {n: i for i, n in enumerate(ref_names)}
    g = np.zeros(R); r = np.zeros(R)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if (rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1
                or rd.mate_is_unmapped):
            continue
        is_g = "gdna" in rd.query_name.lower()
        rf = name2ref.get(bam.get_reference_name(rd.reference_id))
        if rf is None or rf not in by_ref:
            continue
        tl = abs(rd.template_length)
        if tl == 0:
            continue
        lo = min(rd.reference_start, rd.next_reference_start); hi = lo + tl
        s, e, ii = by_ref[rf]; j = np.searchsorted(s, lo, side="right") - 1
        if 0 <= j < R and s[j] <= lo and hi <= e[j]:
            (g if is_g else r)[ii[j]] += 1.0
    bam.close()
    return g, r


def run_len(mid_len, *, n_rna=200000, seed=7, frag_mean=250, read_length=120):
    # 3-exon + gene: small flank exons, internal exon of length mid_len.
    e1 = (2000, 2300)
    e2 = (10000, 10000 + mid_len)
    e3 = (e2[1] + 5000, e2[1] + 5300)
    exons = [e1, e2, e3]
    genes = [("T", "+", exons)]
    glen = e3[1] + 3000
    name = f"mid{mid_len}"
    wd = SCRATCH / f"efl_{name}"
    sc = Scenario(name, genome_length=glen, seed=seed, work_dir=wd, ref_name="chr1")
    sc.add_gene("T", "+", [{"t_id": "T", "exons": exons, "abundance": 100.0}])
    result = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=0.05,
        sim_config=ReadSimConfig(frag_mean=frag_mean, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=read_length, strand_specificity=1.0, seed=seed),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=frag_mean, frag_std=50),
        capture_config=None, nrna_abundance=0.0,
    )
    idx = result.index
    _st, sm, flm, _buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    bsub = BoundarySubstrate.from_payload(pl)
    g, r = _oracle_contained(result.bam_path, ra, idx.ref_names)
    L = np.asarray(ra.region_size_bp, float)
    E_r = region_eff_length(L, fl.rna_pmf)
    flm_rna = fl_mean(fl.rna_pmf)
    # find the internal exon region (the one spanning e2)
    rdf = idx.region_df
    mid = None
    for i in range(len(rdf)):
        if int(rdf['start'].iloc[i]) == e2[0] and int(rdf['end'].iloc[i]) == e2[1]:
            mid = i; break
    rho_true = r[mid] / E_r[mid]
    # its flanking boundaries: left boundary's right side (= mid as right_region) carries acceptor
    # spliced; right boundary's left side (= mid as left_region) carries donor spliced.
    blr = np.asarray(bsub.left_region); brr = np.asarray(bsub.right_region)
    # acceptor side: boundary b with right_region==mid → its right side spliced (mid is acceptor)
    acc_b = np.where(brr == mid)[0]
    don_b = np.where(blr == mid)[0]
    spl_acc = float(bsub.right.mass_spliced[acc_b[0]]) if acc_b.size else 0.0
    spl_don = float(bsub.left.mass_spliced[don_b[0]]) if don_b.size else 0.0
    E_spl_coded = float(spliced_side_eff_length(fl.rna_pmf, np.array([L[mid]]))[0])
    eff_acc = spl_acc / rho_true if rho_true > 0 else np.nan
    eff_don = spl_don / rho_true if rho_true > 0 else np.nan
    print(f" mid_len={mid_len:5d} L={L[mid]:6.0f} rho_true={rho_true:.4f} | "
          f"E_spl_coded={E_spl_coded:6.1f} (flmean/2={flm_rna/2:.1f}) | "
          f"E_spl_eff: don={eff_don:6.1f} acc={eff_acc:6.1f} avg={0.5*(eff_don+eff_acc):6.1f} | "
          f"ratio coded/eff={E_spl_coded/(0.5*(eff_don+eff_acc)):.3f}")
    return mid_len, E_spl_coded, eff_don, eff_acc, flm_rna


if __name__ == "__main__":
    print("How E_spl_effective (= spliced_mass / true_rho_RNA) compares to coded half-triangle:")
    for ml in (300, 600, 1000, 2000, 5000):
        run_len(ml)
