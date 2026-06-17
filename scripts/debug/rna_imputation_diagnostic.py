"""Zero-gDNA false-gDNA diagnostic (the Step-1 precision-rebuild gate).

Builds a MINIMAL zero-gDNA RNA-seq scenario (TRUE f_g == 0 everywhere), runs the REAL calibrate(), dumps
every region node with its solved gDNA fraction, and flags any node that reads FALSE gDNA (f_g notably > 0
with real mass). The Step-1 architecture (CALIBRATION_ARCHITECTURE.md) predicts: with the count prior and the
RNA prior both removed, every node is solved by its **strand likelihood** + the **global foundation** — and
since the strand alone calls the (stranded) nascent RNA, the introns AND exons must settle at f_g ≈ 0 (no
phantom gDNA), with ρ_global → 0 driving the unanchored nodes down too.

Scenario (+ strand only for the targets, plus standard +/- multi-exon genes so the strand model trains):
  * a SINGLE-EXON transcript:  intergenic | exon | intergenic
  * a TWO-EXON transcript:     intergenic | exon1 | intron | exon2 | intergenic
  gdna_config=None (ZERO gDNA); nrna_abundance ~30; strand_specificity=0.99; ~8000 fragments.

The introns are the historical failure (the removed count prior voted confident gDNA there); this script is
the regression repro. Standalone: `python scripts/debug/rna_imputation_diagnostic.py`.
"""

from __future__ import annotations

import dataclasses
import tempfile
from pathlib import Path

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
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

        def _region_class(s: int) -> str:
            exon = bool(s & (BIT_EXON_POS | BIT_EXON_NEG))
            intron = bool(s & (BIT_INTRON_POS | BIT_INTRON_NEG))
            if exon and intron:
                return "exon+intron"
            if exon:
                return "exon"
            if intron:
                return "intron"
            return "intergenic"

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
        mis = [
            {
                "node": (f"idx={int(i)} sig={decode_signature(int(sig[i]))} "
                         f"ts={_TS_NAME.get(int(ts[i]))} class={_region_class(int(sig[i]))} "
                         f"[{int(start[i])},{int(end[i])}]"),
                "solved_fg": round(float(f_g[i]), 4),
                "mass_gdna": round(float(g_mass[i]), 2),
                "mass_rna": round(float(r_mass[i]), 2),
            }
            for i in bad
        ]

        # Per-class worst-case false-gDNA (the architecture gate: introns + exons ≈ 0).
        by_class: dict[str, float] = {}
        for i in range(R):
            if g_mass[i] + r_mass[i] > 5.0:
                k = _region_class(int(sig[i]))
                by_class[k] = max(by_class.get(k, 0.0), float(f_g[i]))

        return {
            "nodes_dump": nodes_dump,
            "mis": mis,
            "max_fg_by_class": by_class,
            "rna_sense_frac": float(result.rna_sense_frac),
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
    print("\n# ---- max false-gDNA f_g by region class (TRUE f_g == 0; want all ≈ 0) ----")
    for k, v in sorted(out["max_fg_by_class"].items()):
        print(f"  {k:>14}: max f_g = {v:.4f}")
    print("\n# ---- MIS-ESTIMATED (false-gDNA, f_g > 0.05) NODES ----")
    if not out["mis"]:
        print("  (none — every node with real mass reads f_g ≈ 0; the Step-1 gate PASSES)")
    for m in out["mis"]:
        print(f"\n{m['node']}")
        print(f"  solved f_g  = {m['solved_fg']}   mGDNA={m['mass_gdna']}  mRNA={m['mass_rna']}")
