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

This replaces the retired ``oracle_region_masses`` (in the deleted ``_metrics``/``oracle_*`` scripts), which
deposited WHOLE fragments by SPAN with no intron-cutting — an INCOMPATIBLE basis with the accumulator the
calibration actually consumes (per-base coverage, introns cut). That mismatch (e.g. it reported 0 RNA in
high-expression exons where the accumulator has the real unspliced exon-body mRNA) confounded earlier
"calibration error" conclusions.

⭐ **FIVE POPULATIONS ON THREE AXES**, each an integer count with two GENOME-strand columns::

    regions             region_contained     the whole path lies inside the region
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
    "region_contained_count",
    "edge_unspliced_count",
    "edge_spliced_count",
    "sj_count",
    "region_start_count",
    # ⭐ The three CONSERVED MASS banks. Two of them were outside this tuple when they landed, which
    # meant the origin split was never validated on them — and they are exactly what `component_shares`
    # reads. Integer fixed point, so sum-to-full is byte-exact here like every other bank.
    "edge_unspliced_mass",
    "edge_spliced_mass",
    "sj_mass",
)

#: The banks a gDNA fragment can NEVER touch — it does not splice. ⭐ ``sj_mass`` joins them, and it is
#: a STRONGER statement than the two counts: the mass is where a spliced fragment's whole share goes
#: when its blocks cross no line, so a gDNA record leaking onto the junction axis shows up here as
#: fragment-scale mass rather than as a single incidence.
_RNA_ONLY_BANKS = ("edge_spliced_count", "sj_count", "sj_mass")

#: Origin -> code for the per-``frag_id`` truth array, so ``ORIGINS[code] == kind``. int8, because a
#: 10 M-fragment condition is then 10 MB and can be carried in a cache blob.
ORIGIN_CODE = {k: i for i, k in enumerate(ORIGINS)}


def frag_id_origins(bam: str, scan_config) -> tuple[np.ndarray, dict]:
    """⭐⭐ **EVERY FRAGMENT'S TRUE ORIGIN, KEYED BY THE SCANNER'S OWN ``frag_id``** — the join that
    lets an instrument ask "which of THIS multi-locus's EM candidates were really gDNA?".

    ``_split_bam`` answers origin questions on the *accumulator* axis: three BAMs, three payloads, one
    per-object bank each. It cannot answer a question about an EM UNIT, because a unit is a row in the
    scored CSR and its identity is a ``frag_id`` — and ``frag_id`` is assigned by the scanner, not by
    the read name. This function supplies the missing key.

    ⭐ **``frag_id`` is the index of the qname GROUP, over the records that survive the scanner's own
    record filter** (``bam_scanner.cpp`` pass 1: ``current_group.frag_id = frag_id++`` once per change
    of qname, after QC-fail / unmapped / duplicate / unpaired rejection). Both mates and every
    secondary alignment of one molecule share a qname, hence one ``frag_id``. This walk is the same
    rule, in the same order, over the same records — and ``rigel.annotate``'s BAM writer already
    re-derives it the same way, which is what makes it a contract rather than an implementation detail.

    ⛔ **``scan_config`` is REQUIRED and is not a convenience.** ``skip_duplicates`` decides whether a
    duplicate-flagged record is filtered, and a filtered record must NOT advance the counter. Reading
    it from the caller's own ``BamScanConfig`` is what stops this walk and the scan it is being joined
    to from disagreeing about where ``frag_id`` 5,000,000 is.

    ⛔⛔ **THE JOIN IS GATED BY A COUNT IDENTITY, NOT BY A STATISTICAL SMELL TEST** — see
    :func:`check_walk_alignment`. Both this walk and pass 1 stream the same file in the same order and
    stamp one id per qname change, so agreeing on ``n_records`` *and* ``n_groups`` is sufficient: two
    monotone counters over one sequence that end level cannot have differed in the middle without one
    of them double-stepping, which the record total would show.

    ⚠ **Do not reach instead for "gDNA never lands on a spliced unit" as the primary check.** It is
    carried as a secondary diagnostic and it is WEAK on this panel, measured: the simulator writes each
    population in a block, so a 10 M-fragment condition has only **15** origin transitions in BAM order
    and a shift of one fragment mislabels ~15 fragments in total. That detector sees a large slip and is
    blind to a small one.

    Returns ``(origins, diag)``: ``int8[n_groups]`` indexed by ``frag_id``, and the record accounting
    (``n_records`` / ``n_filtered`` / ``n_groups`` / ``n_transitions`` / per-origin totals) a caller
    must report — the per-origin totals are the denominators for "how many fragments of this origin
    never became an EM candidate at all", and ``n_transitions`` is what says how sensitive the
    secondary diagnostic is on this substrate.
    """
    from rigel.sim.read_name import parse_origin

    skip_duplicates = bool(scan_config.skip_duplicates)
    codes: list[int] = []
    n_records = n_filtered = 0
    current = None
    with pysam.AlignmentFile(bam, "rb") as fin:
        for rec in fin:
            n_records += 1
            flag = rec.flag
            if (flag & 0x200) or (flag & 0x4) or (skip_duplicates and (flag & 0x400)):
                n_filtered += 1
                continue
            if not (flag & 0x1):
                # The scanner throws on this; a walk that silently tolerated it would be counting
                # groups the scan never made.
                raise AssertionError(
                    f"{bam}: unpaired record {rec.query_name} — the production scanner refuses this "
                    "BAM, so a frag_id walk over it means nothing."
                )
            qname = rec.query_name
            if qname != current:
                current = qname
                codes.append(ORIGIN_CODE[parse_origin(qname).kind])
    origins = np.asarray(codes, dtype=np.int8)
    diag = {
        "n_records": n_records,
        "n_filtered": n_filtered,
        "n_groups": int(origins.shape[0]),
        "n_transitions": int((origins[1:] != origins[:-1]).sum()) if origins.size > 1 else 0,
        "totals": {k: int((origins == c).sum()) for k, c in ORIGIN_CODE.items()},
    }
    return origins, diag


def check_walk_alignment(walk: dict, stats) -> None:
    """⛔⛔ **THE ``frag_id`` JOIN'S ONE HARD GATE — raise, never warn.**

    ``stats.n_read_names`` is incremented once per qname group inside the scanner's own worker, so it
    IS the number of ``frag_id``\\ s pass 1 issued; ``stats.total`` is every record it read. A walk that
    matches both has visited the same records in the same order and cut them into the same number of
    groups, which is the whole of what "``frag_origin[frag_id]`` is that fragment's origin" requires.

    ⛔ A mismatch is not a small error: every unit's origin label shifts, ``Fo`` stays a plausible
    array of counts, and nothing downstream can tell. Hence an exception rather than a printed warning.
    """
    want_records, want_groups = int(stats.total), int(stats.n_read_names)
    if walk["n_records"] != want_records or walk["n_groups"] != want_groups:
        raise RuntimeError(
            f"the read-name walk does NOT reproduce the scanner's frag_id: records "
            f"{walk['n_records']:,} vs {want_records:,}, groups {walk['n_groups']:,} vs "
            f"{want_groups:,}. Every origin label would be shifted, and a shifted label is a "
            "plausible wrong number — refusing to score."
        )


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
    #: ⭐ held records whose ORIGIN attribution the lift could not determine — 0 unless ``drain_with``
    #: was used. It bounds the truth error exactly; a caller must REPORT it (`second_pass.lift_choices`).
    n_ambiguous: int = 0

    @classmethod
    def from_bam(
        cls,
        bam: str,
        index,
        cfg,
        work_dir: Path,
        tag: str,
        full_payload=None,
        drain_with=None,
    ) -> "OracleTruth":
        """Split the BAM by origin, scan each partition, and validate sum-to-full.

        ``full_payload`` lets a caller that has ALREADY scanned the full BAM (e.g. to run
        ``calibrate`` on it) hand that payload in, skipping a redundant full re-scan. It must be
        the production scan of ``bam`` with the same ``cfg`` — sum-to-full then also PROVES the
        oracle partitions reconstruct the exact payload the calibration consumed.

        ⭐⭐ ``drain_with`` makes the oracle valid for a DRAINED tally. Without it, every number an
        instrument reads off this class is an **undrained** one, because the second pass conditions on
        the whole tally and partitions drained independently do not sum to the whole drained
        (TRAPS: draining-breaks-the-oracle). Pass ``(undrained_whole, choices, region_types, junctions)`` — exactly what
        ``pipeline._drain_side_buffer(_lift=...)`` publishes — and each partition is drained by
        REPLAYING the whole's already-drawn choices inside it (`second_pass.lift_choices`).

        ⛔ **The two payloads are different objects and swapping them is the one way to get this
        wrong.** ``full_payload`` takes the **DRAINED** whole — it is what sum-to-full is asserted
        against and what calibration read. ``drain_with[0]`` takes the **UNDRAINED** whole — the drained
        bank is empty by design, so it cannot supply the key pool and every partition would raise. That
        failure is at least loud; ⭐ and in the other direction ``from_parts``' existing sum-to-full
        becomes the drain's own end-to-end identity gate, for free.
        """
        paths, read_counts = _split_bam(bam, work_dir, tag)
        full = full_payload if full_payload is not None else _scan_payload(bam, index, cfg)
        parts = {k: _scan_payload(paths[k], index, cfg) for k in ORIGINS}
        n_ambiguous = 0
        if drain_with is not None:
            from rigel.second_pass import drain, lift_choices

            undrained_whole, choices, region_types, junctions = drain_with
            # ⭐ ALL partitions in ONE call — the per-key queue's state is shared across them, so each
            # of the whole's choices is handed out exactly once (`lift_choices`' own docstring).
            lifted, n_ambiguous = lift_choices(
                undrained_whole, [parts[k] for k in ORIGINS], choices
            )
            parts = {
                k: drain(parts[k], ch, region_types=region_types, junctions=junctions)
                for k, ch in zip(ORIGINS, lifted)
            }
        return cls.from_parts(full, parts, read_counts, n_ambiguous=n_ambiguous)

    @classmethod
    def from_parts(
        cls, full, parts: dict, read_counts: dict | None = None, n_ambiguous: int = 0
    ) -> "OracleTruth":
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
            n_ambiguous=int(n_ambiguous),
        )
        self._validate()
        return self

    def _validate(self) -> None:
        """Sum-to-full on EVERY bank: integers EXACTLY, fractions to the representation and no further.

        ⭐⭐ **ONE CONVENTION, TWO STANDARDS, AND THE SPLIT IS THE POINT.** A COUNT is an integer and
        integer addition is associative, so its partitions must sum to the full payload with **no
        tolerance at all** — that is still the great majority of the banks. A FRACTION is float64, and
        summing three partitions re-associates the additions, so those agree to within the
        representation. ⛔ The budget is `n * EPS` scaled by the magnitude — derived from the deposit
        count, never fitted — and a COUNT bank that needed it fails here as it always did.

        ⛔⛔ **THE PREDECESSOR CAST EVERY BANK TO int64**, which was right when every bank was an
        integer and became **silent corruption** the moment three of them held fractions in (0, 1]:
        `int64` truncates them to zero, so the comparison was zeros against zeros — a vacuous pass that
        would have certified any partition at all. Caught 2026-08-11 while landing the convention.

        ⚠ The tolerance is REAL but it is not the old one: a float32 tolerance once hid this project's
        factor-of-2 bug for months. This one is 1e-16-scale and a factor of 2 is 1e16 times larger.
        """
        eps = float(np.finfo(np.float64).eps)
        for bank in _BANKS:
            full_raw = np.asarray(getattr(self.full, bank))
            if full_raw.dtype == np.float64:
                full = full_raw
                parts = sum(np.asarray(getattr(self.parts[k], bank), np.float64) for k in ORIGINS)
                budget = np.maximum(np.abs(full), 1.0) * max(full.size, 1) * eps
                bad = np.abs(parts - full) > budget
                if bad.any():
                    raise AssertionError(
                        f"oracle INVALID: {bank} partitions do not sum to full within the float64 "
                        f"representation (max|diff|={np.abs(parts - full).max():.3e} over "
                        f"{int(bad.sum())} cells). That is far larger than re-association, so the "
                        "partition is not the production split."
                    )
                continue
            full = full_raw.astype(np.int64)
            parts = sum(np.asarray(getattr(self.parts[k], bank), np.int64) for k in ORIGINS)
            if not np.array_equal(parts, full):
                raise AssertionError(
                    f"oracle INVALID: {bank} partitions do not sum to full "
                    f"(max|diff|={np.abs(parts - full).max()}). This bank is an integer COUNT and "
                    "integer addition is associative, so there is no tolerance to reach for: the "
                    "partition is not the production split, or the accumulator stopped depositing "
                    "per fragment."
                )
        # gDNA is never spliced (physical), on EITHER spliced bank — including the junction axis, which
        # the old 4-channel layout had no room for.
        for bank in _RNA_ONLY_BANKS:
            g = int(np.asarray(getattr(self.parts["gdna"], bank), np.int64).sum())
            if g != 0:
                raise AssertionError(f"oracle INVALID: gdna partition has {g} deposits in {bank}.")

    # ---- per-REGION TRUE counts on the accumulator basis ----
    def region_unspliced(self):
        """``(G, R)`` per region: TRUE contained gDNA vs contained RNA count — the gDNA-vs-RNA competition
        basis the calibration deconvolves. ``R`` = mrna + nrna (exon-body + nascent).

        ⚠ There is no spliced term to exclude: ``region_contained`` is credited only when the fragment
        used no junction, so a region's contained population is unspliced by construction.
        """
        nc = lambda k: np.asarray(self.parts[k].region_contained_count, np.float64).sum(1)  # noqa: E731
        return nc("gdna"), nc("mrna") + nc("nrna")

    def region_true_fg(self):
        """Per-region TRUE gDNA fraction of the contained count (NaN where there is no contained mass)."""
        G, R = self.region_unspliced()
        tot = G + R
        return np.where(tot > 0, G / np.maximum(tot, 1e-12), np.nan), tot

    def region_pools(self) -> dict:
        """Per-region TRUE contained count by ORIGIN × GENOME strand.

        Calibration deconvolves the contained count into ``(RNA₊, RNA₋, gDNA)`` and cannot split mature
        from nascent — that is the downstream EM's job — so ``mat_*``/``nas_*`` here are the TRUE
        composition of the RNA calibration lumps together. All six components sum to the full per-region
        contained count (the validated sum-to-full identity).
        """
        nc = lambda k: np.asarray(self.parts[k].region_contained_count, np.float64)  # noqa: E731
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

        ⭐ The exact mirror of :meth:`region_pools`, on the basis the solver's EDGE slots use — and it is
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

    def component_shares(self) -> dict:
        """⭐⭐ **THE TRUE PER-COMPONENT SHARE AT EVERY LINE, MEASURED — not modelled.**

        ``mass / count`` inside one origin partition IS that component's mean conserved fragment-mass
        per crossing, on the real partition, under the real placement. It needs **no fragment-length
        pmf, no uniform-placement assumption and no reach model** — which is the whole point, because
        the fitted pmfs carry a phantom length gap and an analytic share would import it.

        ⛔ **This is why the arm built on it is an ORACLE arm and not a MODEL arm.** Deriving the share
        from ``truth_fragment_lengths.tsv`` through ``crossing_eff_length`` would be a model evaluated
        on truth inputs, and it would defeat the purpose: the arm could then be wrong in the same
        direction as the thing it is supposed to price.

        ⚠ **1.0 where a component never crossed a line** — the identity, matching
        ``PopulationView.mass_per_crossing``. A zero there would delete mass rather than leave it
        unscaled, and the shipped accessor makes the same choice for the same reason.

        Returns ``{"gdna": float64[n_edges], "rna": float64[n_edges]}``.
        """
        out = {}
        for name, origins in (("gdna", ("gdna",)), ("rna", ("mrna", "nrna"))):
            mass = sum(
                np.asarray(self.parts[k].edge_unspliced_mass, np.float64) for k in origins
            )
            count = sum(
                np.asarray(self.parts[k].edge_unspliced_count, np.float64).sum(axis=1)
                for k in origins
            )
            share = np.ones_like(mass)
            np.divide(mass, count, out=share, where=count > 0)
            out[name] = share
        return out

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

        rna_region = total(subs["mrna"], "region_contained") + total(subs["nrna"], "region_contained")
        rna_edge_unspliced = total(subs["mrna"], "edge_unspliced") + total(
            subs["nrna"], "edge_unspliced"
        )
        spliced_edge = total(full, "edge_spliced")
        return dict(
            mass_gdna_region=total(subs["gdna"], "region_contained"),
            mass_rna_region=rna_region,
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

    G, R = orc.region_unspliced()
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
    # the RNA pool, neither of which production has used since S5.d/TRAPS: pure-and-length-censored.1. A test harness that builds
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
    cal_g = np.asarray(cal.mass_gdna_region, np.float64)
    cal_r = np.asarray(cal.mass_rna_region, np.float64)
    # cal contained total vs the payload's contained count (a region holds no spliced molecule)
    print(
        f"\ncal contained total (g+r)={(cal_g + cal_r).sum():,.0f}  vs TRUE contained={(G + R).sum():,.0f}"
    )
    true_fg, tot = orc.region_true_fg()
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
