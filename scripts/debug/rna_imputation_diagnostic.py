"""RNA-imputation transcript-structure diagnostic (zero-gDNA).

Builds a MINIMAL zero-gDNA RNA-seq scenario that exposes the RNA-imputation error in calibrate(),
runs the REAL calibrate(), dumps every region+boundary node, identifies nodes with false gDNA
(f_g notably > 0 when the TRUE f_g is ~0 everywhere — zero gDNA), and traces each bad node back to
its RNA-prior imputation source (which flank, did the edge cross a TSS/TES or an exon<->intron, was
spliced included in the imputed RNA density).

Scenario (+ strand only for the targets, plus standard +/- multi-exon genes so the strand model trains):
  * a SINGLE-EXON transcript:  intergenic | exon | intergenic
  * a TWO-EXON transcript:     intergenic | exon1 | intron | exon2 | intergenic
  gdna_config=None (ZERO gDNA); nrna_abundance ~30; strand_specificity=0.99; ~8000 fragments.

The transcript-structure rules the imputation must obey (and the failures this exposes):
  (a) TSS/TES (intergenic<->exon) boundaries are a ZERO-RNA black hole — RNA must NOT impute across them.
  (b) exon<->intron (splice junction) boundaries carry one-sided FIXED spliced (mature) + two-sided
      unspliced (nascent/gDNA); imputation across them must carry ONLY the unspliced component.
  (c) an exon's contained mature is DENSE in the exon body but only a thin slice crosses a junction as
      spliced, so a junction crossing UNDER-REPRESENTS the exon's contained RNA -> the imputed RNA
      density is too low -> the exon reads false gDNA.

The script is a standalone repro: `python scripts/debug/rna_imputation_diagnostic.py`.
"""

from __future__ import annotations

import dataclasses
import tempfile
from pathlib import Path

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.rna_density_model import node_rna_density, rna_strand_densities
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NEG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.strand_deconv import cleaned_gdna_count, deconv_sides, strand_deconvolve
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType


def decode_signature(sig: int) -> str:
    """Decode a 4-bit region signature into a human-readable token list."""
    parts = []
    if sig & BIT_EXON_POS:
        parts.append("EXON+")
    if sig & BIT_EXON_NEG:
        parts.append("EXON-")
    if sig & BIT_INTRON_POS:
        parts.append("INTRON+")
    if sig & BIT_INTRON_NEG:
        parts.append("INTRON-")
    return "+".join(parts) if parts else "INTERGENIC"


_TS_NAME = {TS_NONE: "NONE", TS_POS: "POS", TS_NEG: "NEG", TS_AMBIG: "AMBIG"}


def build_and_calibrate(work_dir: Path):
    """Build the zero-gDNA scenario, scan, and return everything calibrate() needs + the result."""
    sc = Scenario("rna_imp_diag", genome_length=30000, seed=11, work_dir=work_dir)
    # TARGET 1: single-exon + transcript (intergenic | exon | intergenic).
    sc.add_gene("g_single", "+", [{"t_id": "T_single", "exons": [(3000, 4000)], "abundance": 120}])
    # TARGET 2: two-exon + transcript (intergenic | exon1 | intron | exon2 | intergenic).
    sc.add_gene(
        "g_two", "+", [{"t_id": "T_two", "exons": [(8000, 9000), (12000, 13000)], "abundance": 120}]
    )
    # Standard multi-exon +/- genes so the strand model trains (well away from the targets).
    sc.add_gene(
        "s_pos", "+",
        [{"t_id": "S_pos", "exons": [(16000, 16500), (17500, 18000), (19000, 19500)], "abundance": 120}],
    )
    sc.add_gene(
        "s_neg", "-",
        [{"t_id": "S_neg", "exons": [(22000, 22500), (23500, 24000), (25000, 25500)], "abundance": 120}],
    )
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(
            frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
            read_length=100, strand_specificity=0.99, seed=11,
        ),
        gdna_config=None,  # ZERO gDNA -> TRUE f_g == 0 everywhere
        nrna_abundance=30.0,
    )
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(pl),
        max_size=flm.max_size,
    )
    result = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
    return sc, idx, pl, ra, sm, fl, cfg, result


def main() -> dict:
    with tempfile.TemporaryDirectory() as td:
        work_dir = Path(td)
        sc, idx, pl, ra, sm, fl, cfg, result = build_and_calibrate(work_dir)

        substrate = CalibrationSubstrate.from_payload(pl, ra)
        c, left_v, right_v = substrate.contained, substrate.left, substrate.right
        R = ra.n_regions

        # ---- Solved per-node masses / f_g -----------------------------------
        g_mass = np.asarray(result.mass_gdna_contained, dtype=np.float64)
        r_mass = np.asarray(result.mass_rna_contained, dtype=np.float64)
        f_g = g_mass / np.maximum(g_mass + r_mass, 1e-12)

        sig = np.asarray(ra.signature).astype(np.int64)
        ts = np.asarray(ra.strand_class)
        start = np.asarray(ra.start)
        end = np.asarray(ra.end)

        cu = (c.n_unspliced_pos + c.n_unspliced_neg).astype(np.int64)
        cs = (c.n_spliced_sense + c.n_spliced_antisense).astype(np.int64)
        # crossing (boundary) counts attributed to this region's two flanking sides
        xu = (left_v.n_unspliced_pos + left_v.n_unspliced_neg
              + right_v.n_unspliced_pos + right_v.n_unspliced_neg).astype(np.int64)
        xs = (left_v.n_spliced_sense + left_v.n_spliced_antisense
              + right_v.n_spliced_sense + right_v.n_spliced_antisense).astype(np.int64)

        # ---- Replicate calibrate()'s RNA-prior inputs EXACTLY ----------------
        region_eff_len = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
        boundary_eff_len = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
        fl_mean = boundary_eff_length(fl.gdna_pmf)
        region_eff_len_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
        rna_boundary_side_eff_len = boundary_side_eff_length(fl.rna_pmf, ra.region_size_bp)
        rna_sense_frac = float(result.rna_sense_frac)
        gdna_od = float(result.gdna_strand_overdispersion)
        rna_od = float(result.rna_strand_overdispersion)

        node_density_raw = node_gdna_density(
            substrate, ra, region_eff_len, fl_mean, need_count_variance=False
        )
        _, left_split, right_split = strand_deconvolve(
            substrate, ra, rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_od, rna_strand_overdispersion=rna_od,
            deconv_quantile=cfg.calibration.gdna_deconv_quantile, n_grid=cfg.calibration.n_grid,
        )
        i0 = cfg.calibration.gdna_strand_info_scale

        def _raw(view):
            return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

        cleaned_left = cleaned_gdna_count(left_split, _raw(substrate.left), i0)
        cleaned_right = cleaned_gdna_count(right_split, _raw(substrate.right), i0)
        node_density = node_gdna_density(
            substrate, ra, region_eff_len, fl_mean, need_count_variance=False,
            gdna_counts=(_raw(substrate.contained), cleaned_left, cleaned_right),
        )
        left_dc, right_dc = deconv_sides(
            substrate, ra, node_density, boundary_eff_len,
            rna_sense_frac=rna_sense_frac, gdna_strand_overdispersion=gdna_od,
            rna_strand_overdispersion=rna_od,
            deconv_quantile=cfg.calibration.gdna_deconv_quantile,
            n_grid=cfg.calibration.n_grid, info_scale=i0,
        )

        # RNA prior uses the CONVERGED region f_g (calibrate passes prev_fg, == the converged f_g
        # at the last consumed pass) and the fixed deconv_sides anchors.
        rna_dens = rna_strand_densities(
            substrate, ra, region_eff_len_rna, rna_boundary_side_eff_len,
            gdna_frac=f_g,
            left_gdna_frac=left_dc.gdna_frac, right_gdna_frac=right_dc.gdna_frac,
            cleaned_left=cleaned_left, cleaned_right=cleaned_right,
        )
        rho_pos, rho_neg = node_rna_density(rna_dens, ra)

        # The region's OWN RNA density = mass_unspliced / region_eff_len_rna (what ρ̂ SHOULD match if
        # the flanks could see the contained mature). f_g, μ_s via calibrate's conversion.
        mass_u = np.maximum(np.asarray(c.mass_unspliced, dtype=np.float64), 1e-9)
        eff_rna = np.maximum(np.asarray(region_eff_len_rna, dtype=np.float64), 1e-9)
        own_rna_density = mass_u / eff_rna

        rho_hat = np.where(np.isfinite(rho_pos), rho_pos, np.where(np.isfinite(rho_neg), rho_neg, np.nan))
        # μ_s = clip(ρ̂·E_rna/M_u, 0, 1) — calibrate's _rna_prior_fraction
        mu_s = np.where(
            np.isfinite(rho_hat), np.clip(rho_hat * eff_rna / mass_u, 0.0, 1.0), np.nan
        )
        # τ_rna via the RNA var~mean curve (rebuild it the same way calibrate does, on this assembly)
        from rigel.calibration.rna_density_model import fit_rna_imputation_varmean

        rna_curve = fit_rna_imputation_varmean(rna_dens)
        geom_rna = (eff_rna / mass_u) ** 2
        rho_safe = np.where(np.isfinite(rho_hat), rho_hat, 0.0)
        sig2 = np.maximum(rna_curve.predict(rho_safe) * geom_rna, 1e-12)
        tau_rna = np.where(np.isfinite(rho_hat), np.minimum(1.0 / sig2, mass_u), 0.0)

        # ---- per-side RNA densities (the imputation predictors) --------------
        sd_pos = rna_dens[TS_POS]
        # left_density/right_density are shared across strands (single-strand region pools both genome
        # strands' unspliced + the motif-sense spliced); pull from the POS entry.
        left_rna_dens = np.asarray(sd_pos.left_density)
        right_rna_dens = np.asarray(sd_pos.right_density)
        left_ok = np.asarray(sd_pos.left_ok)
        right_ok = np.asarray(sd_pos.right_ok)
        # the left-side and right-side SPLICED-sense counts (the mature anchor that gates *_ok)
        l_spl = np.asarray(left_v.n_spliced_sense, dtype=np.float64)
        r_spl = np.asarray(right_v.n_spliced_sense, dtype=np.float64)

        # ---------------------------------------------------------------------
        # NODES DUMP (compact text table)
        # ---------------------------------------------------------------------
        lines = []
        hdr = (
            f"{'idx':>3} {'ref':>3} {'start':>6} {'end':>6} {'sig':>14} {'strand':>6} "
            f"{'cU':>5} {'cS':>5} {'xU':>5} {'xS':>5} {'f_g':>6} {'mGDNA':>8} {'mRNA':>9}"
        )
        lines.append(hdr)
        lines.append("-" * len(hdr))
        for i in range(R):
            lines.append(
                f"{i:>3} {int(ra.ref_id[i]):>3} {int(start[i]):>6} {int(end[i]):>6} "
                f"{decode_signature(int(sig[i])):>14} {_TS_NAME.get(int(ts[i]), '?'):>6} "
                f"{int(cu[i]):>5} {int(cs[i]):>5} {int(xu[i]):>5} {int(xs[i]):>5} "
                f"{f_g[i]:>6.3f} {g_mass[i]:>8.2f} {r_mass[i]:>9.2f}"
            )
        nodes_dump = "\n".join(lines)

        # ---------------------------------------------------------------------
        # Identify false-gDNA nodes (TRUE f_g == 0; flag f_g > 0.05 with real mass)
        # ---------------------------------------------------------------------
        FALSE_GDNA = 0.05
        bad = np.flatnonzero((f_g > FALSE_GDNA) & (g_mass + r_mass > 5.0))

        def _edge_kind(reg_sig: int, nbr_sig) -> str:
            """Classify a region<->neighbour edge by transcript structure (rules a/b/c)."""
            if nbr_sig is None:
                return "ref-terminal"
            reg_exon = bool(reg_sig & (BIT_EXON_POS | BIT_EXON_NEG))
            reg_intron = bool(reg_sig & (BIT_INTRON_POS | BIT_INTRON_NEG))
            nbr_exon = bool(nbr_sig & (BIT_EXON_POS | BIT_EXON_NEG))
            nbr_intron = bool(nbr_sig & (BIT_INTRON_POS | BIT_INTRON_NEG))
            nbr_intergenic = not (nbr_exon or nbr_intron)
            reg_intergenic = not (reg_exon or reg_intron)
            if (reg_exon and nbr_intergenic) or (reg_intergenic and nbr_exon):
                return "TSS/TES (intergenic<->exon: ZERO-RNA black hole)"
            if (reg_exon and nbr_intron) or (reg_intron and nbr_exon):
                return "exon<->intron (splice junction)"
            return "exon<->exon / other"

        mis = []
        ls_sig = np.concatenate([[None], sig[:-1].tolist()]) if R > 1 else [None]
        rs_sig = np.concatenate([sig[1:].tolist(), [None]]) if R > 1 else [None]
        # same-ref guards for neighbours
        same_left = np.zeros(R, bool)
        same_right = np.zeros(R, bool)
        if R > 1:
            same_left[1:] = ra.ref_id[1:] == ra.ref_id[:-1]
            same_right[:-1] = ra.ref_id[:-1] == ra.ref_id[1:]

        for i in bad:
            i = int(i)
            # which flank fed the imputation
            srcs = []
            if left_ok[i]:
                ln = int(ls_sig[i]) if (same_left[i] and ls_sig[i] is not None) else None
                srcs.append(
                    f"LEFT side (boundary {i-1}<->{i}): spliced_sense={l_spl[i]:.1f}, "
                    f"side_RNA_density={left_rna_dens[i]:.5f}, edge={_edge_kind(int(sig[i]), ln)}"
                )
            if right_ok[i]:
                rn = int(rs_sig[i]) if (same_right[i] and rs_sig[i] is not None) else None
                srcs.append(
                    f"RIGHT side (boundary {i}<->{i+1}): spliced_sense={r_spl[i]:.1f}, "
                    f"side_RNA_density={right_rna_dens[i]:.5f}, edge={_edge_kind(int(sig[i]), rn)}"
                )
            if not srcs:
                srcs.append("NO RNA-observable flank (ρ̂ = NaN -> τ_rna=0): prior OFF; node fell to "
                            "strand-likelihood + global foundation only")
            spliced_in = "yes (spliced_sense added to side numerator)" if (left_ok[i] or right_ok[i]) else "n/a"

            # structural violation classification
            reg_exon = bool(int(sig[i]) & (BIT_EXON_POS | BIT_EXON_NEG))
            if not (left_ok[i] or right_ok[i]) and reg_exon:
                struct = ("(c) contained-mature-undercount: exon node has NO same-strand "
                          "RNA-observable flank, so the RNA prior is OFF entirely — its dense "
                          "contained mature is invisible to the imputation -> reads false gDNA")
            elif reg_exon and np.isfinite(rho_hat[i]) and rho_hat[i] < own_rna_density[i] * 0.5:
                struct = ("(c) contained-mature-undercount: a junction crossing's thin spliced slice "
                          "under-represents the exon's dense contained mature (ρ̂ << own RNA density)")
            elif np.isfinite(rho_hat[i]):
                struct = ("(a/b) imputation source structurally wrong (check edge kind: TSS/TES "
                          "should be a ZERO-RNA black hole; exon<->intron must carry only unspliced)")
            else:
                struct = "RNA prior OFF (no anchor); node under-determined"

            mis.append({
                "node": (f"idx={i} sig={decode_signature(int(sig[i]))} ts={_TS_NAME.get(int(ts[i]))} "
                         f"[{int(start[i])},{int(end[i])}]"),
                "solved_fg": round(float(f_g[i]), 4),
                "rho_hat_s": (round(float(rho_hat[i]), 6) if np.isfinite(rho_hat[i]) else None),
                "region_own_rna_density": round(float(own_rna_density[i]), 6),
                "mu_s": (round(float(mu_s[i]), 4) if np.isfinite(mu_s[i]) else None),
                "tau_rna": round(float(tau_rna[i]), 4),
                "imputation_source": " | ".join(srcs),
                "structural_violation": struct,
                "spliced_included": spliced_in,
            })

        return {
            "nodes_dump": nodes_dump,
            "mis": mis,
            "rna_sense_frac": rna_sense_frac,
            "gdna_density_global": float(result.gdna_density_global),
            "total_gdna": float(g_mass.sum()),
            "total_rna": float(r_mass.sum()),
            "n_regions": R,
        }


if __name__ == "__main__":
    out = main()
    print(f"# rna_sense_frac={out['rna_sense_frac']:.3f}  "
          f"gdna_density_global={out['gdna_density_global']:.4g}  "
          f"total_gDNA={out['total_gdna']:.1f}  total_RNA={out['total_rna']:.1f}  "
          f"R={out['n_regions']}\n")
    print(out["nodes_dump"])
    print("\n# ---- MIS-ESTIMATED (false-gDNA) NODES ----")
    for m in out["mis"]:
        print(f"\n{m['node']}")
        print(f"  solved f_g            = {m['solved_fg']}")
        print(f"  rho_hat_s (imputed)   = {m['rho_hat_s']}")
        print(f"  region own RNA dens   = {m['region_own_rna_density']}")
        print(f"  mu_s / tau_rna        = {m['mu_s']} / {m['tau_rna']}")
        print(f"  spliced included      = {m['spliced_included']}")
        print(f"  imputation source     = {m['imputation_source']}")
        print(f"  structural violation  = {m['structural_violation']}")
