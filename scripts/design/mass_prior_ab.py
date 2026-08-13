#!/usr/bin/env python
"""⭐⭐⭐ **THE PRIOR AS A CONSERVED FRAGMENT COUNT** — plan 5.5 against the oracle, on real conditions.

**The question.** ``assemble_priors`` hands the EM a per-locus FRAGMENT COUNT, but ``boundary_unspliced_count``
is ``+1`` on every line a fragment crosses, so a fragment books ``max(K, 1)`` of them. The shipped
assembler repairs that with ``rho = Sum m / Sum S`` integrated over the locus span, and under capture the
repair over-calls by **+15.1 %** with a PERFECT deconvolution fed in (``prior_vs_oracle.py``, arm ``O``
against arm ``F``). The accumulator now carries the count directly — ``boundary_unspliced_mass``, which sums
to one per fragment — so the assembler can stop manufacturing one from a density.

This instrument scores that change end to end, against the same ``O`` and ``F`` arms, on a deterministic
subsample of real ladder conditions.

⛔ **IT MEASURES THE SHIPPED BANK. IT DOES NOT RE-DERIVE IT.** An earlier version of this file walked the
deposit rule itself in Python, because the bank did not exist yet; that walk is gone
(``TRAPS: converge-and-delete``) and with it the caveat that made it awkward — it enumerated no gap
hypotheses, so it deposited the ~1.5 % of fragments production HOLDS. The rule's own gates live in
``tests/native/test_conserved_mass.py`` (16 of them, four independent claims) and the C++ is gated on
byte-identity to the specification. Here there is one fragment stream: production's.

⭐ **WHAT IS STILL A PROTOTYPE** is the ASSEMBLER — :func:`assemble_priors_mass`. Plan 5.5 is not yet in
``src/``; this scores it before it goes there (``TRAPS: panel-before-src``).

⛔ **THE SUBSTRATE CHECK COMES FIRST.** A subsample that does not reproduce the defect cannot show it
gone (``TRAPS: can-the-benchmark-resolve-it``), so the shipped assembler is re-scored on the SAME
subsample as its own control rather than compared to the panel's recorded number
(``TRAPS: re-record-the-baseline``).

Usage::

    python scripts/design/mass_prior_ab.py                       # one capture-ON, one capture-OFF
    python scripts/design/mass_prior_ab.py --conditions A B --fragments 500000
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "tests" / "calibration"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# ⛔ THERE IS NOTHING TO DECODE. This module used to divide every mass bank by
# ``substrate.INV_LENGTH_SCALE`` (2^32) because the banks were FIXED POINT. The whole fixed-point layer
# was deleted at `94d283c0` under ONE NUMERIC CONVENTION — a count is an integer, a fraction is float64,
# and nothing decodes a bank — so the scale factor is gone and the banks are already in real units.
# ⚠ The recorded verdicts in this docstring are unchanged by that: they were computed from the same
# quantity, one representation earlier.

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: One of each capture arm, at a gDNA level where both components are abundant. ⛔ Not ``g00`` (no gDNA
#: to get wrong) and not ``g98`` (no RNA left). ⭐ STRANDED, because the ``O`` arm is fed oracle masses
#: and is therefore calibration-independent — nothing is bought by running it on the blind stratum.
DEFAULT_CONDITIONS = (
    "gdna_g50_ss_0.99_nrna_none_capture_on",
    "gdna_g50_ss_0.99_nrna_none_capture_off",
)


def subsample_bam(source: str, out: Path, target: int, total_fragments: int) -> tuple[int, float]:
    """Write a deterministic qname-hash subsample of a name-sorted BAM. Returns ``(reads, rate)``.

    ⭐ **Keyed on the QNAME, so it COMMUTES with the origin split.** Both mates share a qname and the
    origin is a function of the qname, so the same predicate applied to the whole and to each of the
    three partitions selects the same molecules — which is exactly what ``OracleTruth``'s sum-to-full
    identity needs. A positional or reservoir subsample would commute with nothing.

    ⚠ ``a_g / truth`` is a RATIO, so it is scale-invariant and a subsample preserves it. What a
    subsample does NOT preserve is anything with a fixed floor, which is why the shipped assembler is
    re-scored here rather than compared to the panel's number.
    """
    import pysam

    rate = min(1.0, target / max(int(total_fragments), 1))
    threshold = int(rate * (1 << 32))
    written = 0
    with pysam.AlignmentFile(source, "rb") as fin:
        with pysam.AlignmentFile(str(out), "wb", template=fin) as fout:
            for record in fin:
                digest = hashlib.blake2b(record.query_name.encode(), digest_size=4).digest()
                if int.from_bytes(digest, "big") < threshold:
                    fout.write(record)
                    written += 1
    return written, rate


def boundary_share(payload) -> np.ndarray:
    """``boundary_unspliced_mass / boundary_unspliced_count`` per line — the mean conserved share of a crossing.

    ⭐ This is the whole of plan 5.5's new input: one dimensionless scalar per line that converts an
    object-incidence total into a fragment count. It is 1.0 at a line whose flanking regions both exceed
    every fragment length (a crossing fragment can only cross that one line) and falls toward the region
    spacing where they do not — which is the K-inflation, per line.

    ⛔ Lines with no crossing get **1.0**, the identity. There is no mass there to rescale, and a 0 would
    delete whatever mass the calibration put on a line the accumulator never saw.
    """
    mass = np.asarray(payload.boundary_unspliced_mass, np.float64)
    count = np.asarray(payload.boundary_unspliced_count, np.float64).sum(axis=1)
    share = np.ones_like(mass)
    np.divide(mass, count, out=share, where=count > 0)
    return share


def assemble_priors_mass(calibration, region_arrays, multi_loci, share, eff_len_source):
    """⭐⭐⭐ **PLAN 5.5** — the prior as a CONSERVED FRAGMENT COUNT, with no density conversion::

        a_g(locus) = Sum_r share(r) * mass_gdna_region[r]
                   + Sum_l share(l) * mass_gdna_boundary[l] * (boundary_unspliced_mass[l] / boundary_count[l])

    and the same on the RNA masses, spliced withheld exactly as the shipped assembler withholds it.

    ⭐ ``mass_gdna_region[r]`` is already ``f_g(r) * contained_count[r]`` — one deposit per contained
    fragment — so the region term needs no arithmetic at all. Only the boundary term is rescaled, from the
    K-inflated incidence count onto the conserved mass. ``rho = Sum m / Sum S``, ``span_bp`` and the
    support-weighted pooling are all GONE: the count is in the bank, so nothing manufactures one.

    ⛔ **``gdna_eff_len`` is carried over unchanged** from ``eff_len_source``. Plan 5.5 changes the two
    pseudocounts and nothing else; re-deriving the third array would vary two things at once
    (``TRAPS: one-thing-varied``).

    ⚠ **Nothing is dropped for zero opportunity.** There is no denominator here — the mass IS the count
    — so dropping such mass would simply lose fragments the accumulator really deposited.

    ⭐ **Regions and boundaries are projected on their own axes**, matching the shipped assembler: a region owns
    the fragments contained in it, an boundary owns the fragments that cross it, and a locus collects the
    boundaries that touch its regions. No line is folded onto a flank region.
    """
    from rigel.calibration.priors import (
        LocusPriors,
        _boundary_locus_shares,
        _project_regions_to_loci,
        _sum_by_locus,
    )

    share = np.asarray(share, np.float64)
    n_loci = len(multi_loci)
    e_idx, e_lid, e_w = _boundary_locus_shares(region_arrays, multi_loci, n_loci)

    def component(mass_region, mass_boundary):
        region_part = _project_regions_to_loci(
            region_arrays, multi_loci, n_loci, {"m": np.asarray(mass_region, np.float64)}
        )["m"]
        boundary_part = _sum_by_locus(
            e_idx, e_lid, e_w, np.asarray(mass_boundary, np.float64) * share, n_loci
        )
        return np.maximum(region_part + boundary_part, 0.0)

    # the same spliced withholding the shipped assembler does: mass_rna_boundary is spliced-INCLUSIVE by an
    # existing per-boundary conservation convention, so the certified fraction is subtracted before the
    # unspliced competition sees it.
    rna_boundary_unspliced = np.maximum(
        np.asarray(calibration.mass_rna_boundary, np.float64)
        - np.asarray(calibration.mass_rna_spliced_boundary, np.float64),
        0.0,
    )
    return LocusPriors(
        gdna_prior_count=component(calibration.mass_gdna_region, calibration.mass_gdna_boundary),
        rna_prior_count=component(calibration.mass_rna_region, rna_boundary_unspliced),
        gdna_eff_len=np.asarray(eff_len_source.gdna_eff_len, np.float64).copy(),
    )


# ── one condition ─────────────────────────────────────────────────────────────────────────────────


def _condition_fragments(suite: Path, condition: str) -> int:
    """The condition's own true fragment total, from its truth summary. ⛔ Never counted from the BAM:
    that is a two-minute pass to learn a number the simulator already wrote down."""
    summary = json.loads((suite / condition / "truth_summary.json").read_text())
    return int(sum(summary["origin_counts"].values()))


def measure_one(args, index, pipeline_config, condition: str) -> dict:
    """Subsample, then run ``prior_vs_oracle``'s own arms on it. One scan, one fragment stream."""
    import prior_vs_oracle as PVO

    work = args.work_dir / condition
    work.mkdir(parents=True, exist_ok=True)
    tag = f"{condition}_sub"
    sub = work / f"{tag}.bam"

    t0 = time.perf_counter()
    if sub.is_file():
        print(f"    reusing {sub.name}", flush=True)
    else:
        reads, rate = subsample_bam(
            str(args.suite / condition / "sim_oracle.bam"), sub, args.fragments,
            _condition_fragments(args.suite, condition),
        )
        print(f"    subsample {reads:,} reads at {rate:.4%} ({time.perf_counter() - t0:.0f} s)",
              flush=True)

    result = PVO.measure_condition(
        str(sub), index, pipeline_config, work, tag, oracle_cache=None, drained_arm=False,
    )
    print(f"    P/O/F built: {result.n_loci:,} loci, {time.perf_counter() - t0:.0f} s", flush=True)
    return {"condition": condition, "result": result, "seconds": time.perf_counter() - t0}


# ── the report ────────────────────────────────────────────────────────────────────────────────────


def report(rows, args) -> bool:
    import prior_vs_oracle as PVO

    print()
    print("=" * 112)
    print("  ⭐⭐⭐ PLAN 5.5's PRIOR — the conserved mass, against the ORACLE arms")
    print(f"  target {args.fragments:,} fragments/condition   messages OFF   length_likelihood OFF"
          "   UNDRAINED   O and F from prior_vs_oracle.py")
    print("=" * 112)

    print()
    print("  ⛔ GATE A — the shipped noop gate (the override plumbing is inert)")
    for row in rows:
        bad = [f for f, ok in row["result"].noop_identical.items() if not ok]
        print(f"    {row['condition']:<44} " + ("✅ P reproduced byte-identically on all 3 arrays"
                                                if not bad else f"⛔ {bad}"))
        if bad:
            return False

    print()
    print("  ⭐⭐⭐ GATE B — THE SUBSTRATE: does the subsample reproduce the defect it must remove?")
    print("     ROADMAP 1.1 ④ on the FULL panel: gDNA arm rel 0.179 capture ON, 0.005 capture OFF, "
          "net +15.1 % under capture.")
    print(f"    {'condition':<44} {'F (truth)':>13} {'O shipped':>13} {'O/F':>7} {'rel':>7} "
          f"{'net':>13}")
    print("    " + "-" * 104)
    for row in rows:
        r = row["result"]
        row["O_old"] = PVO.score_arm(r.priors["O"].gdna_prior_count, r.f_gdna)
        s = row["O_old"]
        print(f"    {row['condition']:<44} {s.total_ref:>13,.0f} {s.total_arm:>13,.0f} "
              f"{s.total_arm / max(s.total_ref, 1e-9):>7.3f} {s.rel:>7.3f} {s.net_err:>+13,.0f}")

    print()
    print("  ⭐⭐⭐ THE ANSWER — same masses, same loci, same projection, conserved-mass rescale")
    print(f"    {'condition':<44} {'arm':<14} {'total':>13} {'O/F':>7} {'rel':>7} {'Σ|Δ|':>13} "
          f"{'net':>13}")
    print("    " + "-" * 104)
    for row in rows:
        r = row["result"]
        oracle_calibration = dataclasses.replace(
            r.calibration, **r.oracle.override_masses(r.region_arrays)
        )
        share = boundary_share(r.oracle.full)
        row["share"] = share

        def assemble(share_g, share_r, _c=oracle_calibration, _r=r):
            g = assemble_priors_mass(_c, _r.region_arrays, _r.multi_loci, share_g, _r.priors["O"])
            n = assemble_priors_mass(_c, _r.region_arrays, _r.multi_loci, share_r, _r.priors["O"])
            return g.gdna_prior_count, n.rna_prior_count

        pooled_g, pooled_r = assemble(share, share)
        # ⛔ THE ABLATION: share == 1 everywhere is the RAW per-object incidence sum — the K-inflated
        # count the whole change exists to remove. If it is not clearly worse under capture, the share
        # array is not what produced the row above it (``TRAPS: an-ablation-that-never-ran``).
        raw_g, _ = assemble(np.ones_like(share), np.ones_like(share))
        # ⭐ THE PLAN-8-Q4 CEILING: each component rescaled by its OWN origin's share. NOT implementable
        # — the origin split is the oracle — so it prices the one-pooled-share assumption, not a proposal.
        rna_parts = {k: v for k, v in r.oracle.parts.items() if k != "gdna"}
        split_g, split_r = assemble(
            boundary_share(r.oracle.parts["gdna"]), _pooled_share(rna_parts.values())
        )

        row["O_new"] = PVO.score_arm(pooled_g, r.f_gdna)
        row["O_raw"] = PVO.score_arm(raw_g, r.f_gdna)
        row["O_split"] = PVO.score_arm(split_g, r.f_gdna)
        row["O_new_rna"] = PVO.score_arm(pooled_r, r.f_rna_upper)
        row["O_old_rna"] = PVO.score_arm(r.priors["O"].rna_prior_count, r.f_rna_upper)
        arms = (
            ("O shipped", row["O_old"]),
            ("⛔ share≡1", row["O_raw"]),
            ("⭐ O mass", row["O_new"]),
            ("  ↳ own share", row["O_split"]),
        )
        for i, (label, s) in enumerate(arms):
            print(f"    {row['condition'] if i == 0 else '':<44} {label:<14} {s.total_arm:>13,.0f} "
                  f"{s.total_arm / max(s.total_ref, 1e-9):>7.3f} {s.rel:>7.3f} {s.abs_err:>13,.0f} "
                  f"{s.net_err:>+13,.0f}")
        print("    " + "-" * 104)
    print("    ⛔ `share≡1` is the ABLATION: no rescale at all. It must be clearly worse under capture.")
    print("    ⭐ `↳ own share` is the plan 8 Q4 CEILING, not a proposal — its distance from `O mass` is "
          "what one pooled share costs.")
    print()
    print("    ⚠ THE RNA ARM MUST BE READ WITH ITS REFERENCE. F_rna is spliced-INCLUSIVE and is an "
          "UPPER BOUND — `rna_prior_count`")
    print("      withholds the certified fraction — so `rel` there rewards a LARGER prior and a rule "
          "that removes K-inflation")
    print("      necessarily moves away from it. Read the net, not the rate:")
    for row in rows:
        print(f"      {row['condition']:<44} RNA  shipped rel {row['O_old_rna'].rel:.3f} "
              f"net {row['O_old_rna'].net_err:+,.0f}  →  mass rel {row['O_new_rna'].rel:.3f} "
              f"net {row['O_new_rna'].net_err:+,.0f}   (bound {row['O_old_rna'].total_ref:,.0f})")

    _report_q4(rows)
    return True


def _pooled_share(payloads) -> np.ndarray:
    """:func:`boundary_share` over several origin partitions summed — mass and count from one population."""
    payloads = list(payloads)
    mass = sum(np.asarray(p.boundary_unspliced_mass, np.float64) for p in payloads)
    count = sum(np.asarray(p.boundary_unspliced_count, np.float64).sum(axis=1) for p in payloads)
    share = np.ones_like(mass)
    np.divide(mass, count, out=share, where=count > 0)
    return share


def _report_q4(rows) -> None:
    """⭐ **PLAN 8 Q4** — is ONE share per line unbiased for the gDNA component?

    Plan 5.5 rescales both components at a line by the same ``mass/count``, which assumes the gDNA
    crossings and the RNA crossings there have the same mean conserved share. The origin-split oracle
    carries both separately, so the assumption is a measurement: the bias is
    ``Sum_l count_g[l]*share_pooled[l] - Sum_l mass_g[l]``, in fragments.
    """
    print()
    print("  ⭐ PLAN 8 Q4 — is ONE share per line unbiased for the gDNA arm? (from the origin split)")
    print(f"    {'condition':<44} {'Σ mass_g':>12} {'Σ count_g·share':>16} {'bias':>12} "
          f"{'as % of mass_g':>15} {'mwae |Δf_g|':>12}")
    print("    " + "-" * 116)
    for row in rows:
        parts = row["result"].oracle.parts
        share = row["share"]

        def bank(origin, _p=parts):
            payload = _p[origin]
            return (np.asarray(payload.boundary_unspliced_mass, np.float64),
                    np.asarray(payload.boundary_unspliced_count, np.float64).sum(axis=1))

        mass_g, count_g = bank("gdna")
        mass_m, count_m = bank("mrna")
        mass_n, count_n = bank("nrna")
        mass_r, count_r = mass_m + mass_n, count_m + count_n
        predicted, actual = float((count_g * share).sum()), float(mass_g.sum())

        total_count, total_mass = count_g + count_r, mass_g + mass_r
        live = (total_count > 0) & (total_mass > 0)
        f_count = np.divide(count_g, total_count, out=np.zeros_like(count_g), where=live)
        f_mass = np.divide(mass_g, total_mass, out=np.zeros_like(mass_g), where=live)
        weight = np.where(live, total_mass, 0.0)
        mwae = float((weight * np.abs(f_count - f_mass)).sum() / max(weight.sum(), 1e-12))
        print(f"    {row['condition']:<44} {actual:>12,.0f} {predicted:>16,.0f} "
              f"{predicted - actual:>+12,.0f} {(predicted - actual) / max(actual, 1e-9):>+15.2%} "
              f"{mwae:>12.4f}")
    print("    ⭐ `mwae |Δf_g|` compares the gDNA fraction of a line's INCIDENCES with the gDNA fraction "
          "of its conserved MASS,")
    print("       weighted by that line's mass. Small means one pooled share is unbiased for both.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=list(DEFAULT_CONDITIONS))
    ap.add_argument("--fragments", type=int, default=500_000)
    ap.add_argument("--work-dir", type=Path,
                    default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_mass_prior")
    args = ap.parse_args()

    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    index = TranscriptIndex.load(str(args.index))
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for condition in args.conditions:
        print(f"  … {condition}", flush=True)
        rows.append(measure_one(args, index, PipelineConfig(), condition))
    return 0 if report(rows, args) else 1


__all__ = ["assemble_priors_mass", "boundary_share", "subsample_bam"]

if __name__ == "__main__":
    raise SystemExit(main())
