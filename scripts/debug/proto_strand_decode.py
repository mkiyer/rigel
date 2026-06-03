#!/usr/bin/env python
"""Phase 2 validation — the strand clue (per-node gDNA-fraction posterior).

Demonstrates on the nrna_dc toy:
  - expressed t1 exons → π_g LOW (RNA), TIGHT (confident);
  - thin-count regions → π_g DIFFUSE (≈0.5, wide) — the Risk-A regularization;
  - transcribed t1 introns → π_g < 1 (nascent RNA detected, separated from gDNA) — previewing
    how the strand clue removes the count clue's nascent upper-bias in the Phase-3 joint.

Usage:  python scripts/debug/proto_strand_decode.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.strand_decode import strand_decode
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

_TS = {0: "NONE", 1: "POS", -1: "NEG", 2: "AMBIG"}
_RT = {0: "intergenic", 1: "intron", 2: "exon"}


def main():
    tmp = tempfile.mkdtemp(prefix="sd_")
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(2000, 2500), (3000, 3500), (4000, 4500),
        (5000, 5500), (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)], "abundance": 100}])
    sc.add_gene("gc", "-", [{"t_id": "tc", "exons": [(12000, 12500), (13000, 13500), (14000, 14500),
        (15000, 15500), (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)], "abundance": 0}])
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.65, seed=42)
    res = sc.build_oracle(n_fragments=2000, sim_config=cfg,
                          gdna_config=GDNAConfig(abundance=20, frag_mean=350, frag_std=100,
                                                 frag_min=100, frag_max=1000),
                          nrna_abundance=70)

    _, sm, _, _, payload = scan_and_buffer(str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto"))
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    kappa = fit_strand_balance(sm).kappa_rna

    sd_result = strand_decode(sub, kappa)
    rtype = coarse_type_array(ra.signature)
    c = sub.contained
    sense = np.where(ra.ts_class == -1, c.n_unspliced_neg, c.n_unspliced_pos)
    N = c.n_unspliced_pos + c.n_unspliced_neg

    print(f"\n=== Phase 2: strand clue (nrna_dc) — kappa_rna = {kappa:.4f} ===")
    print(f"  decodable (POS/NEG): {sd_result.n_decodable} of {ra.n_regions} regions")
    print(f"  {'idx':>3} {'ts':>5} {'rtype':>10} {'N':>4} {'sense':>5} {'pi_g':>6} {'sd':>5} {'gene':>7}")
    for i in range(ra.n_regions):
        if N[i] < 0.5 and not sd_result.decodable[i]:
            continue
        gene = "t1" if ra.start[i] < 11000 else ("t_ctrl" if ra.start[i] < 19500 else "")
        print(f"  {i:>3} {_TS.get(int(ra.ts_class[i]),'?'):>5} {_RT.get(int(rtype[i]),'?'):>10} "
              f"{int(N[i]):>4} {int(sense[i]):>5} {sd_result.pi_g[i]:>6.2f} "
              f"{np.sqrt(sd_result.pi_g_var[i]):>5.2f} {gene:>7}")

    dec = sd_result.decodable
    exon = (rtype == 2) & dec
    intron = (rtype == 1) & dec
    hi = exon & (N >= 40)  # well-covered exons (expressed t1)
    lo = dec & (N <= 10)   # thin-count regions
    print(f"\n  well-covered exons (N≥40): mean pi_g={sd_result.pi_g[hi].mean():.2f} "
          f"sd={np.sqrt(sd_result.pi_g_var[hi]).mean():.2f}  → RNA, confident")
    if lo.any():
        print(f"  thin regions (N≤10):       mean pi_g={sd_result.pi_g[lo].mean():.2f} "
              f"sd={np.sqrt(sd_result.pi_g_var[lo]).mean():.2f}  → diffuse (regularized)")
    if intron.any():
        print(f"  transcribed introns:       mean pi_g={sd_result.pi_g[intron].mean():.2f}  "
              f"→ <1 (nascent RNA separated from gDNA; corrects the count clue in Phase 3)")


if __name__ == "__main__":
    main()
