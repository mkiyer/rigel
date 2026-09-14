"""rigel.scan_cache — scan once, calibrate many times.

    Gate: `tests/test_scan_cache.py`

Scanning is the expensive step and calibration is the one under development, so a calibration sweep
that re-scans every condition costs minutes where a cached one costs seconds. This module writes and
reads that cache.

WHAT IS STORED, AND WHAT DELIBERATELY IS NOT
--------------------------------------------
`calibrate()`'s inputs come from three places, and only one of them is expensive:

===========================  ==========================================  ==========
input                        origin                                      cached?
===========================  ==========================================  ==========
``payload``                  the scan                                    yes
``strand_model``             the scan                                    yes
``region_arrays``            ``RegionArrays.from_index``                 no
``boundary_flags``           ``build_boundary_flags_array``              no
``sj``                       ``build_sj_geometry_arrays``                no
``gdna_fl_pmf``/``rna_fl_pmf``  ``build_fl_models(payload)``             no — derived
``config``                   the thing you are varying                   no
``injected_priors``          fitted BY ``calibrate``                     no
===========================  ==========================================  ==========

Anything derivable from the index is rebuilt on load, never stored: it is a fraction of a second
against an index load that happens anyway, and a stored copy is how a cache goes stale against the
thing it describes.

There is no separate fragment-length row. Every fragment-length
histogram — the five length pools and the unconditional anchor they are EB-shrunk toward — is a field OF
the payload, so caching the payload caches them, in one frame, by construction. `build_fl_models`
remains the single source of truth for the derived pmfs, which are still not cached: freezing its
output would mean a change to the fl model silently does not reach a cached scan.

THE KEY NEEDS FOUR PARTS
------------------------
* ``graph_hash`` — the region partition plus the sj CSR. The payload already carries it.
* a PAYLOAD-SCHEMA digest — the accumulator's own field list. None of the other keys moves when the
  accumulator changes, so without it a cache written before an accumulator change is accepted and then
  fails deep inside the loader with a bare ``KeyError``.
* a REACH digest. ``reach`` is consumed by calibration and covered by neither ``partition_hash`` nor
  ``graph_hash`` — correctly, since neither the scan nor the accumulator reads it — so a reach-blind
  key would verify clean against an index rebuild that moved a large share of contiguous reaches while
  both existing hashes stayed byte-identical.
* the scan config, because two scans of one BAM under different settings are different tallies.

Not pickle. A pickle of numpy-holding dataclasses is fragile exactly across the schema changes this
cache has to survive; arrays go to ``.npz``, scalars and provenance to JSON.

Seeding a toy from a genome-scale scan needs ``InjectedCalibrationPriors``, which `calibrate` fits and
stashes in ``_debug["calibration_priors"]``, so that path requires `calibrate` to have run;
``test_population_priors_can_be_extracted_from_a_cached_scan`` covers it.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import typing
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from .scan_payload import AccumulatorPayload

if TYPE_CHECKING:  # pragma: no cover
    from .index import TranscriptIndex

__all__ = [
    "ScanCache",
    "ScanCacheKeyError",
    "calibration_inputs",
    "check_scan_config",
    "index_derived_inputs",
    "payload_schema_digest",
    "reach_digest",
    "read_scan_cache",
    "write_scan_cache",
]

MANIFEST_JSON = "manifest.json"
PAYLOAD_NPZ = "payload.npz"
STRAND_NPZ = "strand.npz"

#: The four reach columns the key must cover. Named here rather than globbed so that a NEW reach column
#: fails loudly instead of quietly escaping the digest.
REACH_COLUMNS = ("reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg")


class ScanCacheKeyError(RuntimeError):
    """The cache does not describe this index or this scan configuration."""


def _digest(*parts: bytes) -> str:
    h = hashlib.blake2b(digest_size=8)
    for part in parts:
        h.update(part)
    return h.hexdigest()


def reach_digest(index: "TranscriptIndex") -> str:
    """Hash the boundary reach columns — the part of the index calibration reads and no other key covers.

    Computed on demand, never stored beside the index. It is stored in the CACHE's manifest, which is a
    different thing: the cache is describing an external artifact it cannot recompute, exactly as
    `index.source_record` does for the FASTA and the GTF.
    """
    boundaries = index.edges_df
    if boundaries is None:
        raise ScanCacheKeyError("index has no edges_df loaded; cannot compute the reach digest")
    missing = [column for column in REACH_COLUMNS if column not in boundaries.columns]
    if missing:
        raise ScanCacheKeyError(f"edges.feather is missing reach columns {missing}")
    parts: list[bytes] = []
    for column in REACH_COLUMNS:
        values = np.ascontiguousarray(boundaries[column].to_numpy())
        parts.append(column.encode())
        parts.append(str(values.dtype).encode())
        parts.append(values.tobytes())
    return _digest(*parts)


def _payload_field_types() -> dict[str, type]:
    """``AccumulatorPayload``'s annotations resolved to real classes.

    ``dataclasses.fields()`` hands back annotation STRINGS under ``from __future__ import annotations``,
    so the nested banks cannot be reconstructed from them. Resolving the hints is what lets the read path
    stay generic instead of carrying a name→class table that a new bank could quietly fall out of.
    """
    return typing.get_type_hints(AccumulatorPayload)


def _nested_dataclass(annotation) -> type | None:
    """The dataclass an annotation names, looking through ``Optional`` / unions. ``None`` if none.

    ⛔ The union is not a corner case: ``AccumulatorPayload.drain`` is ``DrainQC | None``, and
    ``dataclasses.is_dataclass`` is False for that — so a plain check silently treats the whole bank as a
    scalar and every field inside it drops out of the schema key.
    """
    if dataclasses.is_dataclass(annotation):
        return annotation
    for argument in typing.get_args(annotation):
        if dataclasses.is_dataclass(argument):
            return argument
    return None


def _schema_names() -> list[str]:
    """Every name the cache is keyed by, NESTED BANKS INCLUDED, to any depth, in a stable order.

    This is what makes the digest below cover what it claims to, and it must stay FULLY recursive rather
    than one level deep: ``DrainQC`` nests a ``GapCensus`` inside itself, so a one-level walk leaves that
    census's fields out of the key entirely.
    """

    def walk(owner: type, prefix: str) -> list[str]:
        hints = typing.get_type_hints(owner)
        names: list[str] = []
        for field in dataclasses.fields(owner):
            names.append(prefix + field.name)
            nested = _nested_dataclass(hints.get(field.name))
            if nested is not None:
                names += walk(nested, f"{prefix}{field.name}__")
        return names

    return walk(AccumulatorPayload, "")


def payload_schema_digest() -> str:
    """Hash the schema the cached arrays were written under — ``AccumulatorPayload``'s field list and
    the fields of every bank nested inside it.

    No other key covers this. ``graph_hash`` describes the index, ``reach_digest`` the reaches,
    ``scan_config_digest`` the scan settings — none of them changes when the ACCUMULATOR changes. Adding
    a field to a population without this key means a cache written beforehand is accepted and then fails
    deep inside ``_payload_from_parts`` with a bare ``KeyError``, which reads as a bug in the cache
    rather than as a stale cache.

    It recurses because the banks nest: ``DeferredFragments`` puts thirteen array names inside one field
    and every one of them is an ``.npz`` key, and ``DrainQC`` nests a ``GapCensus`` inside itself. A
    top-level walk would leave a renamed field inside either one invisible to the key.
    :func:`_schema_names` is fully recursive and looks through ``Optional``.

    ⛔ It hashes the COLUMN COUNT, not only the name, and the obvious objection — "a dtype or shape
    change is already caught at load by ``_bank``'s assertions" — is false on the cache path. ``_bank``
    and ``_single_column_bank`` validate the C++ dict in
    :meth:`AccumulatorPayload.from_scan_result`; :func:`_payload_from_parts` puts the ``.npz`` arrays
    straight into the payload with no shape check at all. Collapsing a bank from ``[n, 2]`` to ``[n]``
    therefore leaves every field name identical, and a name-only digest would ACCEPT a stale cache and
    fail downstream with a shape error pointing nowhere near its cause.

    The column count is taken from the two axis TABLES rather than from the arrays, because the digest
    must be computable without a payload in hand. A bank moving between ``BANK_AXES`` and
    ``SINGLE_COLUMN_AXES`` IS the shape change, so the tables are the honest source.
    """
    return _digest(
        *(name.encode() for name in _schema_names()),
        *(shape.encode() for shape in _schema_shapes()),
        deposit_digest().encode(),
    )


def deposit_digest() -> str:
    """The deposit-BEHAVIOUR digest — a hash of what the accumulator DOES, not of what it is called.
    Scans a fixed tiny partition with a fixed fragment set and hashes the resulting banks.

    :func:`payload_schema_digest` hashes field NAMES and column counts, and a deposit-RULE change moves
    neither — so without this a cache written under the old rule is accepted by the key and silently
    serves OLD VALUES to NEW CODE. A rule change that alters
    no name at all, such as the sj-boundary rule, is caught by nothing else.

    It needs no version number (the project bans them) and no constant to maintain: it is a MEASUREMENT
    of the current code. Change any deposit rule and it moves; change none and it is stable across runs,
    processes and worker counts, because every channel is an integer and integer addition is
    associative.

    The fixture is deliberately awkward rather than minimal — two annotated sj, a short region whose far
    boundary a fragment may or may not reach, contained / crossing / spliced / sj-only fragments — so
    that a rule change confined to ONE of those cases still moves it.

    ⛔ It runs the NATIVE accumulator, not the specification. The cache holds what the production
    scanner deposited, so the key must certify THAT. Reading the reference here would also make `src`
    depend on `tests`, which is not installed with the package. The two are held byte-identical by
    ``tests/native/test_accumulator_native_parity.py``, and a test asserts this digest agrees across
    both — so a drift between them fails loudly rather than certifying the wrong artifact.
    """
    from rigel._bam_impl import Accumulator  # noqa: PLC0415

    #: Region-bound indices, not coordinates: the sj CSR is keyed by the LEFT BOUNDARY. 260 is bound 3
    #: and 1000 is bound 4; 1120 is bound 6 and 2000 is bound 7.
    region_bounds = np.array([0, 60, 200, 260, 1000, 1060, 1120, 2000, 2400], dtype=np.int64)
    accumulator = Accumulator(
        region_bounds=region_bounds,
        region_types=np.array([0, 2, 2, 1, 2, 2, 1, 2], dtype=np.uint8),
        max_length=1000,
        ref=0,
    )
    accumulator.set_sj(
        np.array(
            [0, 0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=np.int32
        ),  # per-donor-region_bound CSR offsets
        np.array([4, 7], dtype=np.int32),  # acceptor region_bound of each sj
        np.array([1, 1], dtype=np.int8),  # STRAND_POS
    )
    for start, end, introns in (
        (10, 50, ()),  # contained in one region, crosses nothing
        (30, 150, ()),  # crosses one boundary
        (30, 280, ()),  # crosses three boundaries
        (150, 1060, ((260, 1000),)),  # spliced; BOTH blocks cross a boundary
        (210, 1040, ((260, 1000),)),  # spliced; NEITHER block crosses a boundary
        (150, 2100, ((260, 1000), (1120, 2000))),  # two sj AND boundaries
        (1030, 2050, ((1120, 2000),)),  # one sj, one boundary
    ):
        # ``hypotheses=()`` is REQUIRED — the native binding has no default, unlike the specification,
        # whose default IS ``UNSPLICED_ONLY``. An empty set means "nothing to arbitrate", which the
        # accumulator treats as the unspliced-only set.
        accumulator.deposit(
            start=start, end=end, observed_introns=introns, sj_strand=1, hypotheses=()
        )

    parts: list[bytes] = []
    for name in sorted(_schema_names()):
        value = getattr(accumulator, name.split("__")[0], None)
        if isinstance(value, np.ndarray):
            parts.append(name.encode())
            parts.append(np.ascontiguousarray(value).tobytes())
    if not parts:
        raise RuntimeError(
            "deposit_digest hashed NO bank — the native accumulator exposed none of the payload's "
            "field names, so this key would be a constant and would certify nothing."
        )
    return _digest(*parts)


def _schema_shapes() -> list[str]:
    """``name:columns`` for every bank, in a stable order — the shape half of the key above."""
    from .scan_payload import BANK_AXES, N_STRAND_COLUMNS, SINGLE_COLUMN_AXES

    return [f"{name}:{N_STRAND_COLUMNS}" for name, _axis, _dtype in BANK_AXES] + [
        f"{name}:1" for name, _axis, _dtype in SINGLE_COLUMN_AXES
    ]


def _scan_config_digest(scan_config) -> str:
    fields = (
        dataclasses.asdict(scan_config)
        if dataclasses.is_dataclass(scan_config)
        else dict(scan_config)
    )
    return _digest(json.dumps(fields, sort_keys=True, default=str).encode())


@dataclasses.dataclass(frozen=True, slots=True)
class ScanCache:
    """Everything one BAM scan produced that calibration consumes."""

    payload: AccumulatorPayload  # the tally — the expensive artifact
    strand_model: object  # StrandModels, including its per-sj table
    provenance: dict  # the key, the BAM, the scan config, the counts

    # No fragment-length row is stored here. Every fragment-length histogram comes off `payload` — the
    # five length pools and the unconditional anchor they are EB-shrunk toward — so caching the payload
    # caches them in one frame. `fl.npz` is neither written nor read;
    # a cache that still has one on disk loads fine, since an extra file is not a key.


# ── strand model round-trip ──────────────────────────────────────────────────────────────────────
# The 2x2 is the MARGINAL of the per-sj table, and the strand OVERDISPERSION is fitted from the
# table, not the marginal. A cache that kept only the 2x2 would silently disable the dispersion estimate
# — which is one of the population priors the toy seed exists to carry.
_SJ_COLUMNS = ("ref_id", "start", "end", "motif_strand", "n_sense", "n_antisense")
_STRAND_SUBMODELS = ("exonic_spliced", "exonic")
_COUNT_FIELDS = ("pos_pos", "pos_neg", "neg_pos", "neg_neg")


def _strand_arrays(strand_model) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for name in _STRAND_SUBMODELS:
        sub = getattr(strand_model, name)
        out[f"{name}__counts"] = np.array(
            [getattr(sub, field) for field in _COUNT_FIELDS], dtype=np.int64
        )
        table = sub.sj_table
        out[f"{name}__has_table"] = np.array([table is not None], dtype=bool)
        if table is not None:
            for column in _SJ_COLUMNS:
                out[f"{name}__sj_{column}"] = np.ascontiguousarray(getattr(table, column))
    return out


def _strand_from_arrays(data) -> object:
    from .strand_model import SJStrandTable, StrandModel, StrandModels

    built = {}
    for name in _STRAND_SUBMODELS:
        counts = data[f"{name}__counts"]
        table = None
        if bool(data[f"{name}__has_table"][0]):
            table = SJStrandTable(**{c: data[f"{name}__sj_{c}"] for c in _SJ_COLUMNS})
        built[name] = StrandModel(
            **dict(zip(_COUNT_FIELDS, (int(v) for v in counts))), sj_table=table
        )
    return StrandModels(**built)


def write_scan_cache(
    cache_dir: str | Path,
    *,
    payload: AccumulatorPayload,
    strand_model,
    index: "TranscriptIndex",
    bam: str,
    scan_config,
) -> Path:
    """Persist one scan's calibration inputs under *cache_dir*, keyed to this index and config."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # A nested bank's arrays go to the .npz, NOT to the manifest. `dataclasses.asdict` on
    # `DeferredFragments` yields a dict of ndarrays, and the manifest is written with
    # `json.dumps(..., default=str)` — which would stringify each array to a TRUNCATED repr, silently, and
    # read back as text. Nested dataclasses are therefore split by the type of each sub-field: arrays are
    # prefixed `field__sub` and counters stay scalars.
    # The cache holds a SCAN, and a drained payload is not one. The second pass runs after the cache is
    # read, which is what lets one scan be drained repeatedly at different seeds without re-reading the
    # BAM. Writing a drained payload would bake one draw into the cache, and it would also serialise
    # `DrainQC.census_before` through `json.dumps(default=str)` as a stringified repr, silently — the
    # same truncation defect one level down.
    if payload.drain is not None:
        raise ValueError(
            "refusing to cache a DRAINED payload. The cache stores pass one so that the drain can be "
            "re-run at any seed; cache the payload `scan_and_buffer` returned and drain after loading."
        )

    arrays: dict[str, np.ndarray] = {}
    scalars: dict[str, object] = {}
    for field in dataclasses.fields(AccumulatorPayload):
        value = getattr(payload, field.name)
        if isinstance(value, np.ndarray):
            arrays[field.name] = np.ascontiguousarray(value)
        elif dataclasses.is_dataclass(value):
            nested: dict[str, object] = {}
            for sub in dataclasses.fields(value):
                sub_value = getattr(value, sub.name)
                if isinstance(sub_value, np.ndarray):
                    arrays[f"{field.name}__{sub.name}"] = np.ascontiguousarray(sub_value)
                else:
                    nested[sub.name] = sub_value
            scalars[field.name] = nested
        else:
            scalars[field.name] = value
    np.savez_compressed(cache_dir / PAYLOAD_NPZ, **arrays)
    np.savez_compressed(cache_dir / STRAND_NPZ, **_strand_arrays(strand_model))

    manifest = {
        # ── the key ──────────────────────────────────────────────────────────────────────────────
        "graph_hash": payload.graph_hash,
        "reach_digest": reach_digest(index),
        "payload_schema_digest": payload_schema_digest(),
        "scan_config_digest": _scan_config_digest(scan_config),
        # ── provenance: what this cache is OF ────────────────────────────────────────────────────
        "bam": str(Path(bam).resolve()),
        "scan_config": dataclasses.asdict(scan_config)
        if dataclasses.is_dataclass(scan_config)
        else dict(scan_config),
        "payload_scalars": scalars,
    }
    (cache_dir / MANIFEST_JSON).write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str)
    )
    return cache_dir


def read_scan_cache(cache_dir: str | Path, index: "TranscriptIndex", scan_config=None) -> ScanCache:
    """Load a cache and REFUSE it unless it describes this index — and, if ``scan_config`` is given,
    unless it was produced under that configuration.

    Pass ``scan_config`` whenever you have one. Without it this checks only that the manifest is
    consistent with ITSELF; a cache scanned under different settings is a different tally and will load
    silently. The check is :func:`check_scan_config`.
    """
    cache_dir = Path(cache_dir)
    manifest = json.loads((cache_dir / MANIFEST_JSON).read_text())

    expected_reach = reach_digest(index)
    if manifest["reach_digest"] != expected_reach:
        raise ScanCacheKeyError(
            f"cache reach digest {manifest['reach_digest']} != index reach digest {expected_reach}. "
            f"The boundary reaches moved. Neither partition_hash nor graph_hash covers reach, so this is "
            f"the only check that notices: a rebuild can move a large share of contiguous reaches "
            f"with both of those byte-identical. Re-scan against this index."
        )

    # Self-consistency: the manifest records the scan config AND its digest, so a tampered or
    # truncated manifest is caught here rather than surfacing as a mysteriously different tally.
    expected_schema = payload_schema_digest()
    if manifest.get("payload_schema_digest") != expected_schema:
        raise ScanCacheKeyError(
            f"cache payload_schema_digest {manifest.get('payload_schema_digest')!r} != "
            f"{expected_schema!r}. The ACCUMULATOR's schema moved, so these arrays are not the fields "
            f"this build reads — a missing one would otherwise surface as a bare KeyError far from here. "
            f"Re-scan; nothing derivable from the index needs rebuilding."
        )

    recorded_scan_digest = _scan_config_digest(manifest["scan_config"])
    if manifest["scan_config_digest"] != recorded_scan_digest:
        raise ScanCacheKeyError(
            f"cache scan_config_digest {manifest['scan_config_digest']} does not match the scan "
            f"config it records ({recorded_scan_digest}). The manifest is inconsistent with itself."
        )

    payload_scalars = manifest["payload_scalars"]
    if payload_scalars["graph_hash"] != manifest["graph_hash"]:
        raise ScanCacheKeyError(
            "cache manifest graph_hash disagrees with its own payload's graph_hash"
        )

    with np.load(cache_dir / PAYLOAD_NPZ) as data:
        arrays = {name: data[name] for name in data.files}
    with np.load(cache_dir / STRAND_NPZ) as data:
        strand = {name: data[name] for name in data.files}
    payload = _payload_from_parts(arrays, payload_scalars)
    if payload.graph_hash != manifest["graph_hash"]:
        raise ScanCacheKeyError(
            f"cache graph_hash {payload.graph_hash} != manifest graph_hash {manifest['graph_hash']}"
        )
    expected_graph = index.graph_hash
    if payload.graph_hash != expected_graph:
        raise ScanCacheKeyError(
            f"cache graph_hash {payload.graph_hash} != index graph_hash {expected_graph}. The region "
            f"partition or the sj CSR moved; this tally does not describe this index."
        )

    cache = ScanCache(
        payload=payload,
        strand_model=_strand_from_arrays(strand),
        provenance=manifest,
    )
    if scan_config is not None:
        check_scan_config(cache, scan_config)
    return cache


def _payload_from_parts(arrays: dict, scalars: dict) -> AccumulatorPayload:
    """Rebuild the payload from the ``.npz`` arrays and the manifest's scalars.

    Generic over the nested banks: each one's sub-fields are taken from the ``.npz`` when they are arrays
    and from the manifest when they are counters, so a bank that grows an array joins the round trip with
    no edit here. It grows the ``payload_schema_digest`` at the same time, which is what refuses a cache
    written before it existed instead of failing here with a bare ``KeyError``.
    """
    types = _payload_field_types()
    kwargs: dict[str, object] = {}
    for field in dataclasses.fields(AccumulatorPayload):
        nested = types.get(field.name)
        if field.name in arrays:
            kwargs[field.name] = arrays[field.name]
        elif dataclasses.is_dataclass(nested):
            recorded = scalars[field.name]
            kwargs[field.name] = nested(
                **{
                    sub.name: arrays[f"{field.name}__{sub.name}"]
                    if f"{field.name}__{sub.name}" in arrays
                    else recorded[sub.name]
                    for sub in dataclasses.fields(nested)
                }
            )
        else:
            kwargs[field.name] = scalars[field.name]
    return AccumulatorPayload(**kwargs)


def check_scan_config(cache: ScanCache, scan_config) -> None:
    """Refuse a cache produced under a different scan configuration.

    Called by :func:`read_scan_cache` when it is given a ``scan_config``. Two scans of one BAM under
    different settings are different tallies, and nothing else notices.
    """
    expected = _scan_config_digest(scan_config)
    if cache.provenance["scan_config_digest"] != expected:
        raise ScanCacheKeyError(
            f"cache scan_config_digest {cache.provenance['scan_config_digest']} != {expected}. Two "
            f"scans of one BAM under different settings are different tallies."
        )


def index_derived_inputs(index: "TranscriptIndex") -> dict:
    """The calibrate() arguments that come from the INDEX, rebuilt every time.

    Deliberately not cached: they rebuild in well under a second, against an index load that happens
    anyway, and a stored copy is how a cache goes stale against the thing it describes.
    """
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import (
        build_boundary_flags_array,
        build_contiguous_boundary_reach_arrays,
        build_mature_wall_distances,
        build_sj_geometry_arrays,
    )

    region_arrays = RegionArrays.from_index(index)
    return {
        "region_arrays": region_arrays,
        "boundary_flags": build_boundary_flags_array(index),
        # The two WALL inputs the measured-total exposure needs: how far a MATURE template continues
        # past each region bound (spliced bases, MAX over covering isoforms) and the NASCENT genomic
        # reach at each contiguous boundary. Both are annotation-only and sample-independent.
        # `build_mature_wall_distances` dominates this function's cost.
        "mature_walls": build_mature_wall_distances(index, region_arrays),
        "boundary_reach": build_contiguous_boundary_reach_arrays(index),
        # The SJ axis is index-derived too, and it is not optional: `calibrate` refuses an axis
        # whose length disagrees with the payload's `n_sj`, because one addressing a different graph
        # would place every splice on the wrong boundary.
        "sj": build_sj_geometry_arrays(index),
    }


def calibration_inputs(
    cache: ScanCache,
    index: "TranscriptIndex",
    *,
    drain_seed: int | None = None,
    lift_out: dict | None = None,
) -> dict:
    """Exactly the keyword arguments `calibrate` needs, in PRODUCTION'S FRAME.

    The payload is DRAINED here. The cache stores pass one so the drain can re-run at any seed;
    production drains the side buffer before ``build_fl_models``/``calibrate``, so an instrument that
    calibrated the undrained tally would be measuring a library the solver never sees — the certified
    spliced channel alone is materially larger once drained. ``drain_seed`` defaults to the production
    default (``PipelineConfig().second_pass_seed``); pass the run's own when scoring against a
    specific pipeline invocation. ``lift_out`` (a dict) receives the drain's ``_lift`` box —
    ``(undrained, choices, region_types, sj)`` — which is what `calibration._oracle.lift_drain_parts`
    needs to put cached ORACLE PARTITIONS into the same frame; it stays EMPTY when the side buffer
    held nothing, in which case pass one already IS the drained frame.

    The fl models are built exactly as production builds them — on the DRAINED payload, with
    ``region_lengths``/``region_types`` supplied so the two-pool contrast runs (omitting them is a
    supported fl fallback, but it is not what ships, and an instrument must not measure a different
    length model than production's).

    Every fragment-length histogram comes from the PAYLOAD — the five length pools and the unconditional
    anchor they are shrunk toward. One quantity, one source, one
    frame: the scanner's spliced histogram is transcript-space and requires a UNIQUE transcript, while
    the accumulator's `RNA_SPLICED` pool is a structural rule over a larger population; and the anchor
    is `deposited_lengths`, binned at the same `L`.

    The one thing that does NOT come from the payload is the RNA pool's de-tilt: "used an annotated
    sj" is a length-dependent selection, and how much so is a fact about the ANNOTATION.
    """
    from .calibration.fl import build_fl_models
    from .calibration.gdna_density import region_lengths_from_partition
    from .calibration.gdna_opportunity import gdna_opportunity_from_index
    from .calibration.sj_opportunity import crossing_probability_from_index
    from .calibration.splice_graph import build_region_partition_arrays
    from .config import PipelineConfig
    from .pipeline import _drain_side_buffer

    if drain_seed is None:
        drain_seed = PipelineConfig().second_pass_seed
    lift = {} if lift_out is None else lift_out
    payload = _drain_side_buffer(
        cache.payload, index, cache.strand_model, seed=int(drain_seed), _lift=lift
    )

    bounds, offsets, region_types = build_region_partition_arrays(index)
    fl_models = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
        region_lengths=region_lengths_from_partition(bounds, offsets, len(region_types)),
        region_types=region_types,
    )
    return {
        "payload": payload,
        "strand_model": cache.strand_model,
        "gdna_fl_pmf": fl_models.gdna_pmf,
        "rna_fl_pmf": fl_models.rna_pmf,
        **index_derived_inputs(index),
    }
