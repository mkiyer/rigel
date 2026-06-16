"""Phase-2b RNA-imputation ACTIVATION RATE — cheap structural measurement.

Decides whether the (planned) RNA-imputation var~mean component is worth prioritizing by
measuring, on representative scenarios, how many single-strand EXON nodes have an adjacent
count-observable boundary side carrying same-strand MATURE (spliced) crossing mass — i.e. the
nodes the RNA imputation would train on and fire at.

A structural pass over the accumulator payload (substrate + region_arrays + signatures + per-view
spliced/unspliced counts). NO sweep, NO quant. Mirrors how phase1_m2_calibrate_internals.py builds
the substrate.

Scenarios:
  * complex   — the complex-loci battery genome (complex_loci_benchmark.build), AMBIG-heavy.
  * capon     — capture-on, ss=0.99, with-nascent flagship-like (phase1_m2 capon builder).
  * capoff    — capture-off, ss=0.99, with-nascent (phase1_m2 capoff builder).

    python scripts/debug/phase2b_activation_rate.py [scenario=all|complex|capon|capoff]
"""
import dataclasses
import sys
import tempfile

import numpy as np

from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.run_fill import same_ref_left_right
from rigel.calibration.signature import (
    TS_POS,
    TS_NEG,
    TS_AMBIG,
    TS_NONE,
    BIT_EXON_POS,
    BIT_EXON_NEG,
)

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG

# the var~mean spline fit-vs-fallback threshold: variance_model.fit -> power-law below max(k,8)=18.
FIT_THRESHOLD = 18


# --------------------------------------------------------------------------------------------------
# Scenario builders — reuse the phase1_m2 / complex-loci genomes.
# --------------------------------------------------------------------------------------------------
def _scan(idx, bam):
    scan = dataclasses.replace(PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    return pl, ra


def build_phase1(kind, work):
    """capon / capoff — reuse phase1_m2_calibrate_internals.build_scenario."""
    import importlib.util
    import os

    here = os.path.dirname(__file__)
    spec = importlib.util.spec_from_file_location(
        "phase1_m2", os.path.join(here, "phase1_m2_calibrate_internals.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sc, pl, ra, sm, gpmf, rpmf, cfg, R, gpres = mod.build_scenario(kind, work)
    return sc, pl, ra


def build_complex(work):
    """complex battery genome — reuse complex_loci_benchmark.build, then scan."""
    import importlib.util
    import os

    here = os.path.dirname(__file__)
    spec = importlib.util.spec_from_file_location(
        "cxbench", os.path.join(here, "complex_loci_benchmark.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    from rigel.sim import GDNAConfig, ReadSimConfig

    sc, spans = mod.build(work, gdna=120.0, nrna=25.0)
    res = sc.build_oracle(
        n_fragments=max(80000, len(mod.LOCI) * 5000),
        sim_config=ReadSimConfig(
            frag_mean=250, frag_std=50, frag_min=80, frag_max=600, read_length=100,
            strand_specificity=0.99, seed=13,
        ),
        gdna_config=GDNAConfig(abundance=120.0, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000),
        nrna_abundance=25.0,
    )
    idx = res.index
    bam = str(res.bam_path)
    pl, ra = _scan(idx, bam)
    return sc, pl, ra


# --------------------------------------------------------------------------------------------------
# The measurement.
# --------------------------------------------------------------------------------------------------
def measure(kind, pl, ra):
    print(f"\n{'='*92}\n=== ACTIVATION SCENARIO: {kind} ===\n{'='*92}")
    substrate = CalibrationSubstrate.from_payload(pl, ra)
    R = ra.n_regions
    sig = np.asarray(ra.signature).astype(np.int64)
    ref_id = np.asarray(ra.ref_id)
    sc = np.asarray(ra.strand_class)  # TS_POS/TS_NEG/TS_AMBIG/TS_NONE

    rco, bco = count_observable_masks(sig, ref_id)  # region / right-boundary count-observable
    has_exon = (sig & _EXON_BITS) != 0

    # ---- (1) denominator: region-class counts ----
    n_total = R
    n_count_obs = int(rco.sum())  # intron / intergenic (no exon bit)
    is_pos_exon = (sc == TS_POS) & has_exon & ~rco
    is_neg_exon = (sc == TS_NEG) & has_exon & ~rco
    ss_exon = is_pos_exon | is_neg_exon  # single-strand exon nodes (NOT obs, NOT AMBIG)
    n_ss_exon = int(ss_exon.sum())
    n_pos_exon = int(is_pos_exon.sum())
    n_neg_exon = int(is_neg_exon.sum())
    is_ambig = sc == TS_AMBIG
    n_ambig = int(is_ambig.sum())
    n_none = int((sc == TS_NONE).sum())

    print("  -- (1) region-class counts (the denominator) --")
    print(f"    total regions                       : {n_total}")
    print(f"    count-observable (intron/intergenic): {n_count_obs}")
    print(f"      of which TS_NONE (intergenic)     : {n_none}")
    print(f"    single-strand EXON (TS_POS|TS_NEG)  : {n_ss_exon}  (POS={n_pos_exon}, NEG={n_neg_exon})")
    print(f"    AMBIG                               : {n_ambig}")

    # ---- per-side boundary observability for region r (mirror density_model) ----
    # boundary_count_observable[k] describes boundary (k, k+1).
    # region r LEFT side  uses boundary (r-1, r) -> bco[r-1]
    # region r RIGHT side uses boundary (r,   r+1) -> bco[r]
    left_same, right_same = same_ref_left_right(ref_id)
    left_anchor = np.zeros(R, dtype=bool)   # region r has an observable LEFT boundary side
    right_anchor = np.zeros(R, dtype=bool)  # region r has an observable RIGHT boundary side
    if R > 1:
        left_anchor[1:] = bco[:-1] & left_same[1:]
        right_anchor[:-1] = bco[:-1] & right_same[:-1]

    # ---- per-view spliced crossing counts, oriented to the region's strand ----
    # substrate.left  = right side of region r's LEFT boundary
    # substrate.right = left  side of region r's RIGHT boundary
    # n_spliced_sense/antisense are MOTIF-relative (transcript sense/antisense), valid even in AMBIG.
    # For a single-strand region the region's transcript strand IS the sense strand, so same-strand
    # MATURE crossing mass = n_spliced_sense on that boundary side. (For TS_NEG the region's '-' is its
    # sense; the scanner already oriented spliced to the transcript, so n_spliced_sense is correct for
    # either single strand.)
    left_spl_sense = np.asarray(substrate.left.n_spliced_sense, dtype=np.float64)
    right_spl_sense = np.asarray(substrate.right.n_spliced_sense, dtype=np.float64)
    left_spl_anti = np.asarray(substrate.left.n_spliced_antisense, dtype=np.float64)
    right_spl_anti = np.asarray(substrate.right.n_spliced_antisense, dtype=np.float64)

    # A side is "eligible" for region r if it is count-observable AND carries >0 same-strand sense
    # spliced crossing mass.
    left_eligible = left_anchor & (left_spl_sense > 0.0)
    right_eligible = right_anchor & (right_spl_sense > 0.0)

    # ---- (2) RNA-imputation activation numerator ----
    # eligible region = single-strand exon node with >=1 eligible (observable + same-strand sense
    # spliced>0) boundary side.
    region_eligible = ss_exon & (left_eligible | right_eligible)
    n_region_eligible = int(region_eligible.sum())
    # count of eligible (region, boundary-side, strand) PAIRS — each side counts once.
    n_pairs = int((ss_exon & left_eligible).sum() + (ss_exon & right_eligible).sum())
    activation_rate = (n_region_eligible / n_ss_exon) if n_ss_exon > 0 else float("nan")

    print("\n  -- (2) RNA-imputation activation (single-strand exon nodes w/ same-strand spliced anchor) --")
    print(f"    eligible single-strand exon regions : {n_region_eligible} / {n_ss_exon}"
          f"  (activation rate = {100*activation_rate:.1f}% of exon nodes)")
    print(f"    eligible (region, side, strand) PAIRS: {n_pairs}")
    # diagnostics: how many exon nodes have ANY observable boundary side at all (regardless of spliced)?
    n_any_obs_side = int((ss_exon & (left_anchor | right_anchor)).sum())
    n_pairs_anyspl = int(
        (ss_exon & left_anchor & ((left_spl_sense + left_spl_anti) > 0)).sum()
        + (ss_exon & right_anchor & ((right_spl_sense + right_spl_anti) > 0)).sum()
    )
    print(f"    [diag] exon nodes w/ ANY observable boundary side (ignoring spliced): {n_any_obs_side}")
    print(f"    [diag] eligible pairs counting ANY spliced (sense+anti) crossing    : {n_pairs_anyspl}")

    # ---- (3) vs the fit threshold ----
    above = n_pairs > FIT_THRESHOLD
    print("\n  -- (3) vs the var~mean fit threshold (power-law fallback below max(k,8)=18) --")
    print(f"    eligible-pair count = {n_pairs}  -> "
          + (f"ABOVE {FIT_THRESHOLD} (a real monotone-spline var~mean fit)"
             if above else f"AT/BELOW {FIT_THRESHOLD} (power-law / sparsity fallback)"))

    # ---- (4) gDNA-imputation pair count for comparison ----
    # gDNA imputation fires at any non-count-observable region (exon OR AMBIG) with >=1 adjacent
    # observable boundary side (regardless of spliced — it uses the unspliced crossing count).
    non_obs = ~rco
    gdna_region_eligible = non_obs & (left_anchor | right_anchor)
    n_gdna_region = int(gdna_region_eligible.sum())
    n_gdna_pairs = int((non_obs & left_anchor).sum() + (non_obs & right_anchor).sum())
    print("\n  -- (4) gDNA-imputation pairs for comparison --")
    print(f"    non-count-observable regions w/ >=1 observable boundary side: {n_gdna_region}")
    print(f"    gDNA-imputation (region, side) PAIRS                        : {n_gdna_pairs}")
    if n_pairs > 0:
        print(f"    RNA set is {n_gdna_pairs / max(n_pairs,1):.1f}x SMALLER as fraction "
              f"(RNA {n_pairs} vs gDNA {n_gdna_pairs} pairs)")

    # ---- one-line verdict ----
    rich = (n_pairs > FIT_THRESHOLD) and (activation_rate > 0.10)
    starved = n_pairs == 0 or activation_rate < 0.02
    verdict = "DATA-RICH" if rich else ("STARVED" if starved else "MARGINAL")
    print(f"\n  >>> VERDICT [{kind}]: RNA imputation is {verdict}  "
          f"(pairs={n_pairs} {'>' if above else '<='}18, "
          f"activation={100*activation_rate:.1f}%)")

    return {
        "kind": kind,
        "n_total": n_total,
        "n_count_obs": n_count_obs,
        "n_ss_exon": n_ss_exon,
        "n_ambig": n_ambig,
        "n_region_eligible": n_region_eligible,
        "n_pairs": n_pairs,
        "activation_rate": activation_rate,
        "above_18": above,
        "n_gdna_pairs": n_gdna_pairs,
        "verdict": verdict,
    }


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    kinds = ["complex", "capon", "capoff"] if which == "all" else [which]
    rows = []
    for kind in kinds:
        with tempfile.TemporaryDirectory() as work:
            sc = None
            try:
                if kind == "complex":
                    sc, pl, ra = build_complex(work)
                else:
                    sc, pl, ra = build_phase1(kind, work)
                rows.append(measure(kind, pl, ra))
            finally:
                try:
                    if sc is not None:
                        sc.cleanup()
                except Exception:
                    pass

    print(f"\n{'='*92}\n=== SUMMARY ===\n{'='*92}")
    hdr = f"  {'scenario':>9} {'#regions':>8} {'#obs':>6} {'#ssExon':>8} {'#AMBIG':>7} " \
          f"{'#RNApairs':>10} {'activation':>11} {'#gDNApairs':>11} {'verdict':>11}"
    print(hdr)
    for r in rows:
        print(f"  {r['kind']:>9} {r['n_total']:>8} {r['n_count_obs']:>6} {r['n_ss_exon']:>8} "
              f"{r['n_ambig']:>7} {r['n_pairs']:>10} {100*r['activation_rate']:>10.1f}% "
              f"{r['n_gdna_pairs']:>11} {r['verdict']:>11}")


if __name__ == "__main__":
    main()
