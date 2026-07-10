#!/usr/bin/env python3
"""Phase 1 — derive the FL-geometry constant for the count→density estimator.

Theory (both exact for uniform gDNA density ρ, starts/bp):
  * contained region:   E[contained count] = ρ · region_eff_len,   region_eff_len = E[max(0, L−ℓ)]
                        → ρ = contained_count / region_eff_len           (REFERENCE, exact)
  * boundary crossing:  E[#frags covering a point] = ρ · fl_mean
                        → ρ = one_side_crossing_flux / fl_mean           (A2 estimator)

So crossing-ρ / contained-ρ should be 1.0 on a uniform-density (NO-CAPTURE) genome. The dry-run
saw ~1.5× over. This script measures the ratio on count-observable nodes and DECOMPOSES it:
  if the bias is a biased fl_mean (gDNA pmf length-biased short), then
      ratio ≈ fl_mean_true / fl_mean_est       and       ratio·(fl_mean_est/fl_mean_true) ≈ 1.
Phase 0 (intergenic contained deposits, large regions → long fragments fit) should de-bias fl_mean.

Builds a NO-CAPTURE scenario: one multi-exon gene (→ intron/exon/intergenic regions + observable
boundaries) in a large genome, high stranded gDNA, zero nRNA. Run:

  python scripts/debug/phase1_fl_geometry_derivation.py
"""
from __future__ import annotations

import numpy as np
import pysam

from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import _check_region_payload_alignment, scan_and_buffer
from rigel.sim import Scenario
from rigel.sim.read_name import parse_origin
from rigel.sim.reads import GDNAConfig, ReadSimConfig
from rigel.splice import SpliceType

_TYPE = {0: "interg", 1: "intron", 2: "EXON"}
GDNA_FRAG_STD = 100.0


def build(gdna_frag_mean: float) -> object:
    tag = f"fl{int(gdna_frag_mean)}"
    sc = Scenario(
        f"phase1_flgeom_{tag}", genome_length=120000, seed=23,
        work_dir=f"/tmp/rigel_sim/phase1_flgeom_{tag}",
    )
    # multi-exon gene → intron + exon regions and intron↔intergenic / intron↔exon boundaries
    sc.add_gene("G", "+", [{"t_id": "G.1",
                            "exons": [(40000, 41200), (44000, 45200), (48000, 49200),
                                      (52000, 53200), (56000, 57200)],
                            "abundance": 60}])
    return sc.build_oracle(
        n_fragments=400000,
        sim_config=ReadSimConfig(strand_specificity=0.99, frag_mean=250, frag_std=50, seed=23),
        gdna_config=GDNAConfig(abundance=5000.0, frag_mean=gdna_frag_mean, frag_std=GDNA_FRAG_STD),
    )


def measure(res: object, label: str) -> None:
    print(f"\n############### CONDITION: {label} ###############")
    index = res.index
    _s, sm, fla, _b, pl = scan_and_buffer(
        str(res.bam_path), index, BamScanConfig(sj_strand_tag="auto")
    )
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    _check_region_payload_alignment(ra, pl)
    flm = build_fl_models(
        global_counts=fla.global_model.counts,
        rna_counts=fla.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl), max_size=fla.max_size,
    )
    gpmf = flm.gdna_pmf
    reg_eff = region_eff_length(ra.region_size_bp, gpmf)
    fl_mean_est = boundary_eff_length(gpmf)
    sub = CalibrationSubstrate.from_payload(pl, ra)

    ctype = coarse_type_array(np.asarray(ra.signature))
    rids = np.asarray(ra.ref_id)
    reg_obs, bnd_obs = count_observable_masks(np.asarray(ra.signature), rids)

    # ---- True gDNA fragment-length mean from the oracle BAM ----
    gdna_lengths: list[int] = []
    seen: set[str] = set()
    with pysam.AlignmentFile(str(res.bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.is_read1:
                continue
            if rec.query_name in seen:
                continue
            seen.add(rec.query_name)
            if parse_origin(rec.query_name).kind != "gdna":
                continue
            tl = abs(rec.template_length)
            if tl > 0:
                gdna_lengths.append(tl)
    fl_mean_true = float(np.mean(gdna_lengths)) if gdna_lengths else float("nan")

    # ---- Reference density: contained gDNA count / region_eff_len on observable regions ----
    # No-capture, observable (no exon bit) → contained unspliced ≈ pure gDNA (nrna=0).
    c = sub.contained
    cont_cnt = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.float64)
    obs_reg = np.where(reg_obs & (reg_eff > 1.0))[0]
    rho_ref = float(cont_cnt[obs_reg].sum() / reg_eff[obs_reg].sum())

    # per-type reference density (intergenic vs intron) — should match (uniform, no capture)
    print("\n=== Phase 1 FL-geometry derivation (NO CAPTURE, uniform ρ) ===")
    print(f"gDNA frag mean: oracle_true={fl_mean_true:.1f}  "
          f"pmf_estimate={fl_mean_est:.1f}  (est/true={fl_mean_est / fl_mean_true:.3f})")
    for t, name in ((0, "intergenic"), (1, "intron")):
        m = reg_obs & (ctype == t) & (reg_eff > 1.0)
        if m.any():
            rho_t = float(cont_cnt[m].sum() / reg_eff[m].sum())
            print(f"  contained ρ ({name:<10}) = {rho_t:.5f}  "
                  f"[{int(m.sum())} regions, {int(cont_cnt[m].sum())} gDNA frags]")
    print(f"  contained ρ (REFERENCE pooled) = {rho_ref:.5f} fragments/bp")

    # ---- Crossing density: one-side flux / fl_mean on observable boundaries ----
    # boundary r↔r+1 observable; left-region side = sub.right[r], right-region side = sub.left[r+1].
    obs_b = np.where(bnd_obs)[0]
    ends_a = np.asarray(ra.end)
    seam_pos = {}  # region r -> genomic position of boundary r↔r+1
    for r in obs_b:
        rs = r + 1
        if rs < ra.n_regions and rids[r] == rids[rs]:
            seam_pos[int(r)] = int(ends_a[r])  # boundary at region r's right edge

    # ---- Oracle: gDNA templates SPANNING each seam (template_start < b < template_end) ----
    span_cross = {r: 0 for r in seam_pos}
    seen2: set[str] = set()
    seam_items = sorted(seam_pos.items())
    seam_b = np.array([b for _, b in seam_items])
    seam_r = [r for r, _ in seam_items]
    with pysam.AlignmentFile(str(res.bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary or not rec.is_read1:
                continue
            if rec.query_name in seen2:
                continue
            seen2.add(rec.query_name)
            if parse_origin(rec.query_name).kind != "gdna":
                continue
            ts = min(rec.reference_start, rec.next_reference_start)
            te = ts + abs(rec.template_length)
            # template-span crossings
            lo = np.searchsorted(seam_b, ts, side="right")
            hi = np.searchsorted(seam_b, te, side="left")
            for k in range(lo, hi):
                span_cross[seam_r[k]] += 1

    rows = []
    for r in obs_b:
        rs = r + 1
        if rs >= ra.n_regions or rids[r] != rids[rs]:
            continue
        side_lo = float(sub.right.n_unspliced[r])    # left region's side of the seam
        side_hi = float(sub.left.n_unspliced[rs])    # right region's side of the seam
        rows.append((r, ctype[r], ctype[rs], side_lo, side_hi))

    flux_lo = np.array([x[3] for x in rows])
    flux_hi = np.array([x[4] for x in rows])
    # one-side estimator (A2): the pooled per-seam crossing count (avg of two sides) / fl_mean
    cross_rho_est_pooled = float(np.mean((flux_lo + flux_hi) / 2.0)) / fl_mean_est if rows else float("nan")
    cross_rho_true = float(np.mean((flux_lo + flux_hi) / 2.0)) / fl_mean_true if rows else float("nan")

    print(f"\n  observable boundaries: {len(rows)}")
    print(f"  crossing ρ (one-side flux / fl_mean_EST)  = {cross_rho_est_pooled:.5f}")
    print(f"  crossing ρ (one-side flux / fl_mean_TRUE) = {cross_rho_true:.5f}")

    print("\n--- RATIOS (target 1.0 on uniform density) ---")
    print(f"  crossing(est) / contained_ref  = {cross_rho_est_pooled / rho_ref:.4f}")
    print(f"  crossing(true)/ contained_ref  = {cross_rho_true / rho_ref:.4f}")
    print(f"  fl_mean_true / fl_mean_est     = {fl_mean_true / fl_mean_est:.4f}")
    print("\nINTERPRETATION: if crossing(est)/ref ≈ fl_mean_true/fl_mean_est and crossing(true)/ref ≈ 1,")
    print("the residual is a biased fl_mean estimate (FL-geometry), not a geometry constant.\n")

    # ---- Decisive: flux vs oracle template-span-crossing ----
    span_vals = np.array([span_cross.get(int(r), 0) for r, *_ in rows], dtype=float)
    flux_side_mean = (flux_lo + flux_hi) / 2.0
    pred_cover = rho_ref * fl_mean_true  # ρ·fl_mean = expected #templates covering a point
    print("\n--- DECISIVE: what does one-side flux count? ---")
    print(f"  ρ·fl_mean (expected templates covering a point) = {pred_cover:.0f}")
    print(f"  oracle templates spanning seam (mean)           = {span_vals.mean():.0f}")
    print(f"  accumulator one-side flux (mean)                = {flux_side_mean.mean():.0f}")
    print(f"  flux / oracle_span        = {flux_side_mean.mean() / max(span_vals.mean(), 1):.4f}")
    print(f"  oracle_span / ρ·fl_mean   = {span_vals.mean() / pred_cover:.4f}")

    # per-seam: flux/oracle_span (the deposit over-count constant c) by boundary type
    print(f"\n{'bnd':>5} {'seam':>13} {'flux_lo':>8} {'flux_hi':>8} "
          f"{'orcl_span':>9} {'flux/span':>9}")
    cvals = []
    for r, tlo, thi, slo, shi in rows[:25]:
        osp = span_cross.get(int(r), 0)
        cc = (slo + shi) / 2.0 / max(osp, 1)
        cvals.append(cc)
        print(f"{r:>5} {_TYPE[tlo] + '|' + _TYPE[thi]:>13} {slo:>8.0f} {shi:>8.0f} "
              f"{osp:>9} {cc:>9.4f}")
    cvals = np.array(cvals)
    print(f"\n  >>> deposit over-count c = flux/oracle_span: "
          f"mean={cvals.mean():.4f}  std={cvals.std():.4f}  "
          f"[exon-adjacent vs intron-adjacent should match]")


def main() -> None:
    # Two read/insert geometries: if the deposit over-count c varies with FL,
    # it must be SELF-CALIBRATED per library, not hardcoded.
    measure(build(gdna_frag_mean=350.0), "gDNA FL mean ≈ 350 (large mate gap)")
    measure(build(gdna_frag_mean=180.0), "gDNA FL mean ≈ 180 (small/no mate gap)")


if __name__ == "__main__":
    main()
