"""THE calibration oracle — the single, validated ground-truth source for debugging/testing/benchmarking.

Principle: the oracle IS the production accumulator, partitioned by TRUE fragment origin. We split the sim
BAM into gdna / mrna / nrna by read-name origin (:func:`rigel.sim.read_name.parse_origin`), run the SAME
production scanner+accumulator on each partition, and assert the partitions sum to the full payload.
Because the accumulator deposits each fragment independently, this sum-to-full identity PROVES the
partition is the production payload split by origin — no reimplementation, nothing to get subtly wrong.

⭐ **THE TOLERANCE IS GONE (S5.f).** Sum-to-full used to be byte-exact on the integer channels and
approximate on ``boundary_mass_{left,right}`` — a float32 array, because the old accumulator split one
fragment's MASS across the objects it touched. The new one deposits ``+1`` on every object touched, so
every bank is an integer sum and the whole identity is exact. ``boundary_mass_tol`` had no successor
and is deleted rather than carried at zero: a tolerance parameter that must be zero is a claim that
some comparison is approximate, and none is.

This replaces the retired ``oracle_node_masses`` (in the deleted ``_metrics``/``oracle_*`` scripts), which
deposited WHOLE fragments by SPAN with no intron-cutting — an INCOMPATIBLE basis with the accumulator the
calibration actually consumes (per-base coverage, introns cut). That mismatch (e.g. it reported 0 RNA in
high-expression exons where the accumulator has the real unspliced exon-body mRNA) confounded earlier
"calibration error" conclusions.

⭐ **FIVE POPULATIONS ON THREE AXES**, each an integer count with two GENOME-strand columns::

    nodes             node_contained     the whole path lies inside the node
                      node_spanning      one segment covers the node whole
    contiguous edges  edge_unspliced     the mixture being deconvolved
                      edge_spliced       certified RNA -- gDNA cannot be spliced
    junction edges    sj_count           pure RNA by construction

⛔ **"Spliced" is a BANK now, not a channel.** The predecessor packed unspliced-± and spliced-sense/
antisense into one 4-column array, which put two strand conventions in one schema. gDNA fragments are
never spliced, and that is validated here as ``edge_spliced`` and ``sj_count`` being identically zero in
the gdna partition — a stronger statement than "columns 2 and 3 are zero", because it covers the
junction axis the old layout had no room for.
"""

from __future__ import annotations
import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.sim.read_name import parse_origin

ORIGINS = ("gdna", "mrna", "nrna")

#: Every per-object bank on the payload. Sum-to-full is asserted over ALL of them: a bank left out of
#: this tuple is a bank the oracle would silently stop validating.
_BANKS = (
    "node_contained_count",
    "node_spanning_count",
    "edge_unspliced_count",
    "edge_spliced_count",
    "sj_count",
    "node_start_count",
)

#: The banks a gDNA fragment can NEVER touch — it does not splice.
_RNA_ONLY_BANKS = ("edge_spliced_count", "sj_count")


def _split_bam(bam: str, out_dir: Path, tag: str) -> tuple[dict[str, str], dict[str, int]]:
    """Split a name-sorted BAM into per-origin BAMs. Both mates share a qname ⇒ same origin ⇒ same file;
    iteration order (hence name-sort) is preserved. EVERY read is written to exactly one partition — the
    'account for every fragment' guarantee — and the total is asserted to reconcile."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {k: str(out_dir / f"{tag}.{k}.bam") for k in ORIGINS}
    counts = {k: 0 for k in ORIGINS}
    n_in = 0
    with pysam.AlignmentFile(bam, "rb") as fin:
        w = {k: pysam.AlignmentFile(paths[k], "wb", template=fin) for k in ORIGINS}
        for r in fin:
            n_in += 1
            k = parse_origin(
                r.query_name
            ).kind  # raises on any unclassifiable read (no silent drop)
            w[k].write(r)
            counts[k] += 1
        for x in w.values():
            x.close()
    if sum(counts.values()) != n_in:
        raise AssertionError(f"oracle split dropped reads: in={n_in} out={sum(counts.values())}")
    return paths, counts


def _scan_payload(bam: str, index, cfg):
    from dataclasses import replace as dc

    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, _sm, _buf, payload = scan_and_buffer(bam, index, sc)
    return payload


@dataclass
class OracleTruth:
    """Validated per-origin accumulator payloads for one condition. Construct via :meth:`from_bam`, which
    runs the sum-to-full validation as a HARD gate (raises if the partition does not reconstruct the full
    payload — i.e. if the oracle is ever not trustworthy)."""

    full: object
    parts: dict  # origin -> payload
    read_counts: dict  # origin -> reads written (every input read accounted for)

    @classmethod
    def from_bam(
        cls,
        bam: str,
        index,
        cfg,
        work_dir: Path,
        tag: str,
        full_payload=None,
    ) -> "OracleTruth":
        """Split the BAM by origin, scan each partition, and validate sum-to-full.

        ``full_payload`` lets a caller that has ALREADY scanned the full BAM (e.g. to run
        ``calibrate`` on it) hand that payload in, skipping a redundant full re-scan. It must be
        the production scan of ``bam`` with the same ``cfg`` — sum-to-full then also PROVES the
        oracle partitions reconstruct the exact payload the calibration consumed."""
        paths, read_counts = _split_bam(bam, work_dir, tag)
        full = full_payload if full_payload is not None else _scan_payload(bam, index, cfg)
        parts = {k: _scan_payload(paths[k], index, cfg) for k in ORIGINS}
        return cls.from_parts(full, parts, read_counts)

    @classmethod
    def from_parts(cls, full, parts: dict, read_counts: dict | None = None) -> "OracleTruth":
        """Assemble from payloads that are ALREADY scanned, and validate — the entry point a cache
        needs.

        ⭐ **The validation is not optional here either, and that is the whole reason this exists
        rather than callers constructing the dataclass directly.** Splitting and re-scanning is
        minutes per condition, so any instrument that runs the panel more than once wants to persist
        the partitions; but a cached oracle that skipped sum-to-full would be a *silently* wrong
        truth source feeding everything downstream. Re-running the identity over loaded arrays costs
        milliseconds, so the cached path is exactly as gated as the scanned one.

        ⚠ ``read_counts`` is bookkeeping (every input read accounted for) and a cache need not carry
        it; it defaults to ``-1`` per origin, which is distinguishable from a real count of zero.
        """
        self = cls(
            full=full,
            parts=parts,
            read_counts=read_counts if read_counts is not None else {k: -1 for k in ORIGINS},
        )
        self._validate()
        return self

    def _validate(self) -> None:
        """Sum-to-full on EVERY bank, EXACTLY — no tolerance anywhere.

         ⭐ Every bank is an integer count now, so ``np.array_equal`` is the right comparison for all of
         them. The predecessor could only be exact on two of its four arrays because the other two were
         float32 fractional MASS; a tolerance is what hid this project's factor-of-2 bug for months
        and there is no longer any reason to carry one.
        """
        for bank in _BANKS:
            full = np.asarray(getattr(self.full, bank), np.int64)
            parts = sum(np.asarray(getattr(self.parts[k], bank), np.int64) for k in ORIGINS)
            if not np.array_equal(parts, full):
                raise AssertionError(
                    f"oracle INVALID: {bank} partitions do not sum to full "
                    f"(max|diff|={np.abs(parts - full).max()}). The partition is not the "
                    "production split, or the accumulator stopped depositing per fragment."
                )
        # gDNA is never spliced (physical), on EITHER spliced bank — including the junction axis, which
        # the old 4-channel layout had no room for.
        for bank in _RNA_ONLY_BANKS:
            g = int(np.asarray(getattr(self.parts["gdna"], bank), np.int64).sum())
            if g != 0:
                raise AssertionError(f"oracle INVALID: gdna partition has {g} deposits in {bank}.")

    # ---- per-NODE TRUE counts on the accumulator basis ----
    def node_unspliced(self):
        """``(G, R)`` per node: TRUE contained gDNA vs contained RNA count — the gDNA-vs-RNA competition
        basis the calibration deconvolves. ``R`` = mrna + nrna (exon-body + nascent).

        ⚠ There is no spliced term to exclude: ``node_contained`` is credited only when the fragment
        used no junction, so a node's contained population is unspliced by construction.
        """
        nc = lambda k: np.asarray(self.parts[k].node_contained_count, np.float64).sum(1)  # noqa: E731
        return nc("gdna"), nc("mrna") + nc("nrna")

    def node_true_fg(self):
        """Per-node TRUE gDNA fraction of the contained count (NaN where there is no contained mass)."""
        G, R = self.node_unspliced()
        tot = G + R
        return np.where(tot > 0, G / np.maximum(tot, 1e-12), np.nan), tot

    def node_pools(self) -> dict:
        """Per-node TRUE contained count by ORIGIN × GENOME strand.

        Calibration deconvolves the contained count into ``(RNA₊, RNA₋, gDNA)`` and cannot split mature
        from nascent — that is the downstream EM's job — so ``mat_*``/``nas_*`` here are the TRUE
        composition of the RNA calibration lumps together. All six components sum to the full per-node
        contained count (the validated sum-to-full identity).
        """
        nc = lambda k: np.asarray(self.parts[k].node_contained_count, np.float64)  # noqa: E731
        g, m, n = nc("gdna"), nc("mrna"), nc("nrna")
        return dict(
            gdna_pos=g[:, 0],
            gdna_neg=g[:, 1],
            mat_uns_pos=m[:, 0],
            mat_uns_neg=m[:, 1],
            nas_uns_pos=n[:, 0],
            nas_uns_neg=n[:, 1],
        )

    def edge_pools(self) -> dict:
        """Per-LINE TRUE crossing counts by ORIGIN × GENOME strand, plus the certified-RNA bank.

        ⭐ The exact mirror of :meth:`node_pools`, on the basis the solver's EDGE slots use — and it is
        ONE set of numbers per line, not a left/right pair. The predecessor summed ``left + right``
        because the old accumulator split one crossing across two faces; there is nothing to sum.

        ``*_spl`` is ``edge_spliced``: molecules that crossed this line CONTIGUOUSLY having spliced
        elsewhere. It is a different population from ``sj_count`` (:meth:`junction_flux`), which never
        crossed the line at all — it jumped.
        """
        eu = lambda k: np.asarray(self.parts[k].edge_unspliced_count, np.float64)  # noqa: E731
        es = lambda k: np.asarray(self.parts[k].edge_spliced_count, np.float64)  # noqa: E731
        g, m, n = eu("gdna"), eu("mrna"), eu("nrna")
        return dict(
            gdna_pos=g[:, 0],
            gdna_neg=g[:, 1],
            mat_uns_pos=m[:, 0],
            mat_uns_neg=m[:, 1],
            nas_uns_pos=n[:, 0],
            nas_uns_neg=n[:, 1],
            mat_spl=es("mrna").sum(1),
            nas_spl=es("nrna").sum(1),
        )

    def junction_flux(self) -> dict:
        """Per-JUNCTION TRUE flux by origin. ⚠ ``gdna`` is identically zero and is returned anyway —
        an all-zero row is the statement "gDNA does not splice", and omitting it would make the
        validator blind to a partition that suddenly produced one."""
        sj = lambda k: np.asarray(self.parts[k].sj_count, np.float64).sum(1)  # noqa: E731
        return {k: sj(k) for k in ORIGINS}

    def override_masses(self, region_arrays) -> dict:
        """The TRUE ``CalibrationResult`` mass arrays on all three axes, built DIRECTLY from the
        per-origin substrates — the exact schema ``calibrate`` assembles, so it can be fed via
        ``dataclasses.replace(cal, **override_masses(ra))`` as the perfect-calibration lever.

        gDNA = the gdna partition's count; RNA = the (mrna + nrna) count, spliced-INCLUSIVE on the edge
        axis to match ``chain_edge_deconv``'s ``rna = (1−f_g)·unspliced + spliced``. Conservation
        ``gdna + rna = the full object count`` holds on both axes because the partitions sum to the full
        payload (the validated identity).

        ⚠ ``mass_rna_junction`` is the FULL payload's junction flux, not the RNA partitions' — they are
        equal by the same identity, and taking it from ``full`` says so rather than re-deriving it.
        """
        from rigel.calibration.substrate import CalibrationSubstrate

        subs = {k: CalibrationSubstrate.from_payload(self.parts[k], region_arrays) for k in ORIGINS}
        full = CalibrationSubstrate.from_payload(self.full, region_arrays)

        def total(sub, population):
            return np.asarray(getattr(sub, population).count, np.float64).sum(1)

        rna_node = total(subs["mrna"], "node_contained") + total(subs["nrna"], "node_contained")
        rna_edge_unspliced = total(subs["mrna"], "edge_unspliced") + total(
            subs["nrna"], "edge_unspliced"
        )
        spliced_edge = total(full, "edge_spliced")
        return dict(
            mass_gdna_node=total(subs["gdna"], "node_contained"),
            mass_rna_node=rna_node,
            mass_gdna_edge=total(subs["gdna"], "edge_unspliced"),
            mass_rna_edge=rna_edge_unspliced + spliced_edge,
            mass_rna_spliced_edge=spliced_edge,
            mass_rna_junction=total(full, "junction"),
        )


def _main():
    ap = argparse.ArgumentParser()
    # ⚠ Defaults track the CURRENT panel. Both were stale for a deleted suite once and cost a session
    # its first hour; if they are wrong again the fix is here, not in the caller.
    ap.add_argument("condition", nargs="?", default="gdna_gdna100_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--suite", default=str(Path.home() / "Downloads/rigel_runs/suite/pilot"))
    ap.add_argument("--index", default=str(Path.home() / "Downloads/rigel_runs/suite/rigel_index"))
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
    index = TranscriptIndex.load(args.index)
    cfg = PipelineConfig()
    bam = f"{args.suite}/{args.condition}/sim_oracle.bam"
    print(f"=== ORACLE {args.condition} ===")
    orc = OracleTruth.from_bam(bam, index, cfg, wd, args.condition)
    print("VALIDATION PASSED: per-origin partitions sum to the full production payload.")

    G, R = orc.node_unspliced()
    print(
        f"\nTRUE contained count: gDNA={G.sum():,.0f}  RNA={R.sum():,.0f}  "
        f"(RNA = exon-body mRNA + nascent)"
    )

    # calibration accuracy on the CORRECT basis: compare per-region f_g
    from rigel.calibration import calibrate
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.splice_graph import (
        build_edge_flags_array,
        build_junction_geometry_arrays,
    )
    from dataclasses import replace as dc

    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, buffer, payload = scan_and_buffer(bam, index, sc)
    ra = RegionArrays.from_index(index)
    # ⭐ One object, one frame — the same call production makes. This harness used to build its
    # own mixture: the scanner's histogram as the anchor and the scanner's SPLICED_ANNOT category as
    # the RNA pool, neither of which production has used since S5.d/C2.1. A test harness that builds
    # a different model from the shipped one is calibrating something the tool does not ship.
    # ⭐ Both divisors, exactly as production passes them — a harness that builds a different model
    # from the shipped one is calibrating something the tool does not ship.
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index

    fl = build_fl_models(
        payload,
        junction_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    cal = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=cfg.calibration,
        junctions=build_junction_geometry_arrays(index),
        edge_flags=build_edge_flags_array(index),
    )
    cal_g = np.asarray(cal.mass_gdna_node, np.float64)
    cal_r = np.asarray(cal.mass_rna_node, np.float64)
    # cal contained total vs the payload's contained count (a node holds no spliced molecule)
    print(
        f"\ncal contained total (g+r)={(cal_g + cal_r).sum():,.0f}  vs TRUE contained={(G + R).sum():,.0f}"
    )
    true_fg, tot = orc.node_true_fg()
    cal_fg = np.where((cal_g + cal_r) > 0, cal_g / np.maximum(cal_g + cal_r, 1e-12), np.nan)
    ok = np.isfinite(true_fg) & np.isfinite(cal_fg)
    w = tot[ok]
    print("\n=== CALIBRATION ACCURACY on the correct (accumulator) basis ===")
    print(
        f"  contained gDNA mass: cal={cal_g.sum():,.0f}  true={G.sum():,.0f}  "
        f"net err={(cal_g - G).sum():+,.0f}  Σ|err|={np.abs(cal_g - G).sum():,.0f}"
    )
    mwae = float(np.sum(w * np.abs(cal_fg[ok] - true_fg[ok])) / max(w.sum(), 1))
    print(f"  mass-weighted |Δf_g| = {mwae:.4f}   (0 = perfect per-region gDNA fraction)")
    dirn = np.maximum(G - cal_g, 0.0)  # gDNA under-called (leaks to RNA)
    print(
        f"  directional gDNA under-call (RNA over-attribution) = {dirn.sum():,.0f}  "
        f"over-call = {np.maximum(cal_g - G, 0).sum():,.0f}"
    )


if __name__ == "__main__":
    _main()
