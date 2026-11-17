"""The calibration oracle: the validated ground truth every truth-scoring instrument reads.

The oracle is the production accumulator, partitioned by true fragment origin: the sim BAM is split into gdna / mrna /
nrna by read-name origin (:func:`rigel.sim.read_name.parse_origin`), the production scanner and accumulator run on each
partition, and the partitions are asserted to sum to the full payload. Because the accumulator deposits each fragment
independently, that identity proves the partition is the production payload split by origin, with nothing reimplemented.
Four populations sit on three axes, each an integer count with two genome-strand columns: ``region_contained`` (the
whole path lies inside the region); ``boundary_unspliced`` (the mixture being deconvolved) and ``boundary_spliced``
(certified RNA) at a contiguous boundary; ``sj_count`` (pure RNA) at an sj. gDNA cannot splice, checked as
``boundary_spliced`` and the sj banks being identically zero in the gdna partition.

A helper module, so it has no row on the shelf; `_oracle_arms.py` caches and reads it.

Gate: ``tests/calibration/test_oracle.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from rigel.calibration.substrate import CalibrationSubstrate
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import DEFERRED_RECORD_FIELDS
from rigel.second_pass import drain
from rigel.sim.read_name import STRAND_FORWARD, parse_origin

ORIGINS = ("gdna", "mrna", "nrna")

#: The RNA reads again, split by transcript strand: an additional partition beside :data:`ORIGINS`, not a replacement.
#: The populations are ``{gDNA, RNA+, RNA−}``, so a per-population measurement needs per-slot RNA truth keyed by the
#: strand of the transcript the fragment came from. The payload's own columns are genome strand, which equals transcript
#: strand only at ``ss 0.99``, wrong on the unstranded half this partition exists for. The simulator writes the
#: transcript's ``f``/``r`` into the read name (`wgs_engine._t_strand_chars`), so no new simulation is needed.
#: ``ORIGINS`` stays the accumulator's three-way axis that instruments iterate, and since both partitions describe the
#: same reads, ``rna_pos + rna_neg == mrna + nrna`` is a real check.
RNA_STRAND_ORIGINS = ("rna_pos", "rna_neg")

#: Every per-object bank on the payload. Sum-to-full is asserted over all of them, so a bank missing here would go
#: unvalidated.
BANKS = (
    "region_contained_count",
    "boundary_unspliced_count",
    "boundary_spliced_count",
    "sj_count",
    "region_start_count",
    # The three conserved-mass banks, exactly what `component_shares` reads. They are float64, so sum-to-full holds
    # within the representation rather than byte-exact (see `_validate`).
    "boundary_unspliced_mass",
    "boundary_spliced_mass",
    "sj_mass",
)

#: The banks a gDNA fragment can never touch, since it does not splice. ``sj_mass`` is the strongest of the three: a
#: spliced fragment's whole share goes there when its blocks cross no boundary, so a gDNA record leaking onto the sj
#: axis shows as fragment-scale mass rather than a single incidence.
_RNA_ONLY_BANKS = ("boundary_spliced_count", "sj_count", "sj_mass")


def split_by_origin(bam: str, out_dir: Path, tag: str) -> tuple[dict[str, str], dict[str, int]]:
    """Split a name-sorted BAM into per-origin BAMs. Both mates share a qname and so an origin and a file; iteration
    order (hence name-sort) is preserved. Every read is written to exactly one partition, and the total is asserted to
    reconcile."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {k: str(out_dir / f"{tag}.{k}.bam") for k in ORIGINS}
    counts = {k: 0 for k in ORIGINS}
    n_in = 0
    with pysam.AlignmentFile(bam, "rb") as fin:
        w = {k: pysam.AlignmentFile(paths[k], "wb", template=fin) for k in ORIGINS}
        for r in fin:
            n_in += 1
            k = parse_origin(r.query_name).kind  # raises on any unclassifiable read (no silent drop)
            w[k].write(r)
            counts[k] += 1
        for x in w.values():
            x.close()
    if sum(counts.values()) != n_in:
        raise AssertionError(f"oracle split dropped reads: in={n_in} out={sum(counts.values())}")
    return paths, counts


def split_rna_by_transcript_strand(bam: str, out_dir, tag: str):
    """Split a name-sorted sim BAM's RNA reads into ``rna_pos`` / ``rna_neg`` by transcript strand. gDNA reads are
    excluded (they have no transcript), and every RNA read lands in exactly one partition; an unclassifiable strand
    raises, as in :func:`split_by_origin`. Returns ``(paths, counts)``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {k: str(out_dir / f"{tag}.{k}.bam") for k in RNA_STRAND_ORIGINS}
    counts = {k: 0 for k in RNA_STRAND_ORIGINS}
    n_rna = 0
    with pysam.AlignmentFile(bam, "rb") as fin:
        writers = {k: pysam.AlignmentFile(paths[k], "wb", template=fin) for k in RNA_STRAND_ORIGINS}
        for record in fin:
            origin = parse_origin(record.query_name)
            if origin.kind == "gdna":
                continue
            # an RNA name's strand letter is its transcript's; `parse_origin` refuses any other letter
            key = "rna_pos" if origin.strand == STRAND_FORWARD else "rna_neg"
            writers[key].write(record)
            counts[key] += 1
            n_rna += 1
        for writer in writers.values():
            writer.close()
    if sum(counts.values()) != n_rna:
        raise AssertionError("the RNA strand split dropped reads")
    return paths, counts


def _scan_payload(bam: str, index, cfg):
    _stats, _sm, _buf, payload = scan_and_buffer(bam, index, cfg.scan)
    return payload


def lift_choices(whole, parts, choices: np.ndarray):
    """Carry hypothesis choices made on the whole library onto its origin partitions.

    Splitting a BAM by origin and re-scanning reconstructs the pass-one payload exactly, because pass one deposits each
    fragment independently; that identity makes an origin-split oracle valid. The second pass breaks it: its multinomial
    is scored against the whole payload's densities, so draining three partitions separately would not sum to draining
    the whole. So the choices are scored and drawn once on the whole library, and each fragment's chosen hypothesis is
    replayed inside whichever partition holds it; every fragment then deposits independently given its choice, and each
    partition's drained tally is comparable with the whole's::

        scores  = score_held_fragments(whole, ...)
        choices = choose_hypotheses(scores, whole, seed=...)      # <- once, on the whole
        drained = drain(whole, choices, ...)                      # the tally the solver reads
        lifted, ambiguous = lift_choices(whole, partitions, choices)     # the truth, per origin
        for p, ch in zip(partitions, lifted):
            drain(p, ch, ...)

    Every partition must be lifted in one call, which is why the signature takes a sequence: a record's choice is
    consumed from a per-key queue whose state must be shared across partitions, or two could take the same entry and
    leave another unused. A one-partition caller passes ``[p]``.

    The key is the bank's own canonical sort key, a record's sequence of placements, each its
    ``DEFERRED_RECORD_FIELDS`` tuple (imported rather than restated), its observed introns, and each of its hypotheses'
    implied introns and strand, everything the C++ sorts on. :class:`DeferredFragments` guarantees that two records
    tying on that key are identical, so no tie-break is needed. Identical records have
    identical hypothesis sets (enumeration reads the span and the annotation, never the origin), so a local hypothesis
    index transfers between them.

    When several whole-library records share a key and are split across origins, no partition can know which it holds,
    and the assignment is greedy in canonical order. The deposits are interchangeable but the origin attribution is not,
    so the count of such records is returned; a caller must report it, and it bounds the truth error exactly. On
    distinct-span substrates it is 0 and the lift is exact.

    Returns ``(lifted, n_ambiguous)``: one ``int64`` array per partition, aligned with that partition's
    ``deferred.hypothesis_offsets``.

    :raises ValueError: if a partition holds a record ``whole`` does not, or more copies of one key than ``whole`` does;
    either means the payloads are not a partition of one scan.
    """
    w = whole.deferred
    parts = list(parts)
    choices = np.asarray(choices, np.int64)
    if choices.shape[0] != w.n_fragments:
        raise ValueError(
            f"lift_choices got {choices.shape[0]} choices for {w.n_fragments} whole-library held "
            "records. One choice per held record, in the bank's canonical order."
        )

    def _keys(d):
        scalars = list(zip(*(np.asarray(getattr(d, f), np.int64).tolist() for f in DEFERRED_RECORD_FIELDS)))

        def runs(offsets, values, stride):
            """Each item's run of ``values``, as a tuple of ``stride``-tuples, under CSR ``offsets``."""
            off = np.asarray(offsets, np.int64).tolist()
            flat = np.asarray(values, np.int64).tolist()
            return [
                tuple(tuple(flat[k : k + stride]) for k in range(stride * off[i], stride * off[i + 1], stride))
                for i in range(len(off) - 1)
            ]

        observed = runs(d.observed_intron_offsets, d.observed_introns, 2)
        implied = runs(d.hypothesis_intron_offsets, d.hypothesis_introns, 2)
        hyp_sj = np.asarray(d.hypothesis_sj_strand, np.int64).tolist()
        hyp_off = np.asarray(d.hypothesis_offsets, np.int64).tolist()
        # a placement as the native canonical key holds it: its scalar fields, its observed introns, and each of its
        # hypotheses' implied introns and strand, so records tie only when they are byte-identical
        placements = [
            (
                scalars[q],
                observed[q],
                tuple((implied[h], hyp_sj[h]) for h in range(hyp_off[q], hyp_off[q + 1])),
            )
            for q in range(len(scalars))
        ]
        offsets = np.asarray(d.placement_offsets, np.int64)
        return [
            tuple(placements[int(offsets[i]) : int(offsets[i + 1])])
            for i in range(d.n_fragments)
        ]

    # key -> the whole's choices for that key, in canonical order, consumed left to right across all partitions so each
    # entry is handed out exactly once.
    pool: dict[tuple, list[int]] = {}
    for i, k in enumerate(_keys(w)):
        pool.setdefault(k, []).append(int(choices[i]))
    taken: dict[tuple, int] = {}

    out, ambiguous = [], 0
    for p in parts:
        d = p.deferred
        lifted = np.zeros(d.n_fragments, np.int64)
        for j, k in enumerate(_keys(d)):
            run = pool.get(k)
            if run is None:
                raise ValueError(
                    f"lift_choices: a partition holds record {j} with key {k}, which is not in the whole "
                    "library's held bank. The payloads are not a partition of one scan — check they came "
                    "from the same BAM."
                )
            at = taken.get(k, 0)
            if at >= len(run):
                raise ValueError(
                    f"lift_choices: the partitions hold more records with key {k} ({at + 1}) than the "
                    f"whole library does ({len(run)}). They are not a partition of one scan."
                )
            if len(run) > 1:
                ambiguous += 1
            lifted[j] = run[at]
            taken[k] = at + 1
        out.append(lifted)
    return out, ambiguous


def lift_drain_parts(lift: dict, parts_list: list) -> tuple[list, int]:
    """Drain a list of partition payloads by replaying the whole's already-drawn choices.

    One call per exact partitioning of the whole (the choice queue is shared across the list, per :func:`lift_choices`).
    The list order matters for the greedy tie assignment, so a caller draining two partitionings (e.g. origins and RNA
    strands) should put any shared member (``gdna``) first in both, which hands it the identical choice slice. Returns
    ``(drained_parts, n_ambiguous)``; an empty ``lift`` (no held fragments) returns the parts unchanged with 0.
    """
    if not lift:
        return list(parts_list), 0
    lifted, n_amb = lift_choices(lift["undrained"], list(parts_list), lift["choices"])
    return [
        drain(p, ch, region_types=lift["region_types"], sj=lift["sj"])
        for p, ch in zip(parts_list, lifted)
    ], int(n_amb)


@dataclass
class OracleTruth:
    """Validated per-origin accumulator payloads for one condition. Every constructor runs the sum-to-full validation
    as a hard check and raises if the partition does not reconstruct the full payload."""

    full: object
    parts: dict  # origin -> payload
    read_counts: dict  # origin -> reads written (every input read accounted for)
    #: held records whose origin attribution the lift could not determine; 0 when no lift replayed held records. It
    #: bounds the truth error exactly, and a caller must report it (:func:`lift_choices`).
    n_ambiguous: int = 0
    #: Drained frame only (``full.drain is not None``): the gdna partition's deposits in each RNA-only bank. The
    #: whole-library drain draws a spliced hypothesis for some true-gDNA fragments, so production's tally violates "gDNA
    #: never splices" as a statement about deposits, and the truth records it rather than refusing to describe the frame
    #: the calibrator reads. ``None`` in the pass-one frame, where such deposits are impossible and remain a hard
    #: failure.
    gdna_spliced_leak: "dict | None" = None

    @classmethod
    def from_bam(cls, bam: str, index, cfg, work_dir: Path, tag: str) -> "OracleTruth":
        """Split the BAM by origin, scan the whole and each partition, and validate sum-to-full.

        Every number read off the result is in the pass-one frame: the second pass conditions on the whole tally, so
        independently drained partitions would not sum to the drained whole. The drained frame is
        :meth:`from_cached_parts`, which replays the whole's choices (:func:`lift_choices`).
        """
        paths, read_counts = split_by_origin(bam, work_dir, tag)
        full = _scan_payload(bam, index, cfg)
        parts = {k: _scan_payload(paths[k], index, cfg) for k in ORIGINS}
        return cls.from_parts(full, parts, read_counts)

    @classmethod
    def from_cached_parts(
        cls, drained_full, parts: dict, lift: dict, read_counts: dict | None = None
    ) -> "OracleTruth":
        """Assemble a drained-frame oracle from cached pass-one partitions.

        ``drained_full`` is the whole payload after ``pipeline.drain_side_buffer(..., lift=lift)``; ``parts`` are
        freshly loaded undrained partition payloads (drained here in place, so pass ones nothing else holds); ``lift``
        is the ``lift`` box the drain published. An empty ``lift`` means the side buffer was empty, so pass one is the
        drained frame and the parts stand. ``from_parts``' sum-to-full then runs on the drained frame, the lift's
        end-to-end identity check.
        """
        if not lift:
            return cls.from_parts(drained_full, parts, read_counts)
        drained, n_amb = lift_drain_parts(lift, [parts[k] for k in ORIGINS])
        return cls.from_parts(
            drained_full, dict(zip(ORIGINS, drained)), read_counts, n_ambiguous=n_amb
        )

    @classmethod
    def from_parts(
        cls, full, parts: dict, read_counts: dict | None = None, n_ambiguous: int = 0
    ) -> "OracleTruth":
        """Assemble from payloads that are already scanned, and validate: the entry point a cache needs.

        Splitting and re-scanning takes minutes per condition, so an instrument running the panel more than once
        persists the partitions; a cached oracle that skipped sum-to-full would be a silently wrong truth source.
        Re-running the identity over loaded arrays costs milliseconds, so the cached path is as checked as the scanned
        one. ``read_counts`` is bookkeeping a cache need not carry; it defaults to ``-1`` per origin, distinguishable
        from a real count of zero.
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
        """Sum-to-full on every bank: integers exactly, fractions to the representation and no further.

        A count is an integer and integer addition is associative, so its partitions must sum to the full payload with
        no tolerance; that is most banks. A fraction is float64, and summing three partitions re-associates the
        additions, so those agree within ``n * EPS`` scaled by the magnitude, derived from the deposit count. A count
        bank that needed the tolerance fails. Casting every bank to int64 would truncate the (0, 1] fraction banks to
        zero and certify any partition at all.
        """
        eps = float(np.finfo(np.float64).eps)
        for bank in BANKS:
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
        # gDNA is never spliced, but the frames differ in what that means for deposits. In the pass-one frame a spliced
        # deposit in the gdna partition is impossible and is a hard failure, on either spliced bank and the sj axis. In
        # the drained frame (``full.drain is not None``) the whole-library drain draws a spliced hypothesis for some
        # true-gDNA fragments, so such deposits are part of production's tally; the truth records them per bank
        # (``gdna_spliced_leak``) instead of refusing the frame the calibrator reads.
        leak = {
            bank: int(np.asarray(getattr(self.parts["gdna"], bank), np.int64).sum())
            for bank in _RNA_ONLY_BANKS
        }
        if getattr(self.full, "drain", None) is not None:
            self.gdna_spliced_leak = leak
        else:
            for bank, g in leak.items():
                if g != 0:
                    raise AssertionError(
                        f"oracle INVALID: gdna partition has {g} deposits in {bank}."
                    )

    def component_shares(self) -> dict:
        """The true per-component share at every boundary, measured rather than modelled.

        ``mass / count`` inside one origin partition is that component's mean conserved fragment mass per crossing, on
        the real partition under the real placement. It needs no fragment-length pmf, uniform-placement assumption or
        reach model, which is what makes an arm built on it an oracle arm rather than a model arm: a share derived from
        ``truth_fragment_lengths.tsv`` through ``crossing_eff_length`` would be a model on truth inputs and could err in
        the same direction as what it prices.

        1.0 where a component never crossed a boundary, the identity, matching ``PopulationView.mass_per_crossing``; a
        zero would delete mass rather than leave it unscaled.

        Returns ``{"gdna": float64[n_boundaries], "rna": float64[n_boundaries]}``.
        """
        out = {}
        for name, origins in (("gdna", ("gdna",)), ("rna", ("mrna", "nrna"))):
            mass = sum(
                np.asarray(self.parts[k].boundary_unspliced_mass, np.float64) for k in origins
            )
            count = sum(
                np.asarray(self.parts[k].boundary_unspliced_count, np.float64).sum(axis=1)
                for k in origins
            )
            share = np.ones_like(mass)
            np.divide(mass, count, out=share, where=count > 0)
            out[name] = share
        return out

    def override_masses(self, region_arrays) -> dict:
        """The true ``CalibrationResult`` mass arrays on all three axes, built directly from the per-origin substrates
        in the schema ``calibrate`` assembles, so ``dataclasses.replace(cal, **override_masses(ra))`` gives a perfect
        calibration.

        gDNA is the gdna partition's count; RNA is the (mrna + nrna) count, spliced-inclusive on the boundary axis to
        match ``chain_boundary_deconv``'s ``rna = (1−f_g)·unspliced + spliced``. ``gdna + rna`` equals the full object
        count on both axes because the partitions sum to the full payload. ``count_rna_sj`` is the full payload's sj
        flux, equal to the RNA partitions' by the same identity.
        """
        subs = {k: CalibrationSubstrate.from_payload(self.parts[k], region_arrays) for k in ORIGINS}
        full = CalibrationSubstrate.from_payload(self.full, region_arrays)

        def total(sub, population):
            return np.asarray(getattr(sub, population).count, np.float64).sum(1)

        rna_region = total(subs["mrna"], "region_contained") + total(
            subs["nrna"], "region_contained"
        )
        rna_boundary_unspliced = total(subs["mrna"], "boundary_unspliced") + total(
            subs["nrna"], "boundary_unspliced"
        )
        spliced_boundary = total(full, "boundary_spliced")
        return dict(
            count_gdna_region=total(subs["gdna"], "region_contained"),
            count_rna_region=rna_region,
            count_gdna_boundary=total(subs["gdna"], "boundary_unspliced"),
            count_rna_boundary=rna_boundary_unspliced + spliced_boundary,
            count_rna_spliced_boundary=spliced_boundary,
            count_rna_sj=total(full, "sj"),
        )
