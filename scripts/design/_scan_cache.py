"""Scan once, calibrate many times: the scan cache every calibration instrument reads. A helper module, so it has
no row on the shelf; `build_scan_cache.py` writes the caches.

    Gate: `tests/test_scan_cache.py`

Scanning is the expensive step, so a calibration sweep reads a cached scan instead of re-scanning. The cache stores the
scan's two outputs, the accumulator payload and the strand model. What the index can derive (``RegionArrays``, the
boundary flags, the sj geometry) and the fragment-length pmfs are rebuilt on load, so no stored copy can go stale; the
length histograms are payload fields and are cached with it.

The key has five parts:

* ``graph_hash``: the region partition plus the sj CSR, carried by the payload;
* the payload-schema digest (:func:`payload_schema_digest`): the accumulator's and the strand model's fields, and what
  the deposit rule does;
* the reach digest: the splice-junction reaches calibration reads, which no other hash covers;
* the transcript-model digest (:func:`transcript_model_digest`): every transcript's reference, strand and exons, which
  resolution reads to decide which transcripts hold a fragment, and which two annotations can differ in with the
  graph and the reaches byte-identical;
* the scan settings that can change a tally.

A table derived from a cache (`calibration_oracle.py`'s ``slot_truth.npz``) records :func:`scan_cache_key`, the whole
key with the BAM's digest, and is refused against any other scan.

The BAM is not part of the key, since the caller names it, but the cache records its path, size, modification time and
a digest of its bytes, and refuses a BAM at the recorded path whose bytes differ. Arrays go to ``.npz`` and scalars and
provenance to JSON; a pickle would be fragile across the schema changes the cache must survive.

The bank axes and the strand-model classes are imported where they are read, not at module load, so the layout the key
hashes is always the one the classes have when it is computed.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import os
import typing
from pathlib import Path

import numpy as np

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import RegionType
from rigel.calibration.splice_graph import (
    EDGE_KIND_SJ,
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.native import Accumulator
from rigel.pipeline import drain_side_buffer, library_fl_models
from rigel.scan_payload import AccumulatorPayload

logger = logging.getLogger(__name__)

MANIFEST_JSON = "manifest.json"
PAYLOAD_NPZ = "payload.npz"
STRAND_NPZ = "strand.npz"

#: The four reach columns the key covers; an index missing one is refused.
REACH_COLUMNS = ("reach_lo_pos", "reach_hi_pos", "reach_lo_neg", "reach_hi_neg")


class ScanCacheKeyError(RuntimeError):
    """The cache does not describe this index or this scan configuration."""


def _digest(*parts: bytes) -> str:
    h = hashlib.blake2b(digest_size=8)
    for part in parts:
        h.update(part)
    return h.hexdigest()


def reach_digest(index: TranscriptIndex) -> str:
    """Hash the splice-junction rows' reach columns, the part of the index calibration reads that no other key covers.

    Only the sj rows: `build_sj_geometry_arrays`, the one reader of a reach, reads only those, and the contiguous rows
    hold 0. Computed on demand and stored only in the cache's manifest.
    """
    boundaries = index.edges_df
    if boundaries is None:
        raise ScanCacheKeyError("index has no edges_df loaded; cannot compute the reach digest")
    missing = [column for column in REACH_COLUMNS if column not in boundaries.columns]
    if missing:
        raise ScanCacheKeyError(f"edges.feather is missing reach columns {missing}")
    sj_rows = boundaries["kind"].to_numpy() == EDGE_KIND_SJ
    parts: list[bytes] = []
    for column in REACH_COLUMNS:
        values = np.ascontiguousarray(boundaries[column].to_numpy()[sj_rows])
        parts.append(column.encode())
        parts.append(str(values.dtype).encode())
        parts.append(values.tobytes())
    return _digest(*parts)


def transcript_model_digest(index: TranscriptIndex) -> str:
    """Hash every transcript's model, in ``t_index`` order: its reference, strand and synthetic flag, and its exons
    as the resolver holds them (the exon CSR the index built at load). Resolution reads which transcripts hold a
    fragment's blocks and junctions, and two annotations can share every region, boundary, junction and reach while
    differing in that."""
    t_df = index.t_df
    offsets, starts, ends, _cumsum = index.build_exon_csr()
    parts: list[bytes] = []
    for name, values in (
        ("ref", np.asarray(t_df["ref"].astype(str).to_numpy(), dtype=str)),
        ("strand", np.asarray(t_df["strand"].to_numpy(), dtype=np.int64)),
        ("is_synthetic", np.asarray(t_df["is_synthetic"].to_numpy(), dtype=bool)),
        ("exon_offsets", np.asarray(offsets, dtype=np.int64)),
        ("exon_starts", np.asarray(starts, dtype=np.int64)),
        ("exon_ends", np.asarray(ends, dtype=np.int64)),
    ):
        parts.append(name.encode())
        parts.append(np.ascontiguousarray(values).tobytes())
    return _digest(*parts)


def scan_cache_key(manifest: dict) -> str:
    """One digest of a cache's whole key and the BAM it was scanned from, from its manifest: what a table derived
    from the cache records, so the table is refused against any other scan."""
    key = {
        name: manifest.get(name)
        for name in (
            "graph_hash",
            "reach_digest",
            "payload_schema_digest",
            "transcript_model_digest",
        )
    }
    key["scan_settings"] = _scan_config_digest(manifest["scan_config"])
    key["bam"] = (manifest.get("bam") or {}).get("blake2b")
    return _digest(json.dumps(key, sort_keys=True).encode())


def bam_identity(bam: str | Path) -> dict:
    """What the cache records about the BAM it scanned: the resolved path, the size, the modification time and a
    digest of the bytes. The time lets an unchanged file skip the digest; the digest decides when the time moved."""
    path = Path(bam).resolve()
    stat = path.stat()
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, lambda: hashlib.blake2b(digest_size=8)).hexdigest()
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "blake2b": digest,
    }


def check_bam_identity(recorded, bam: str | Path | None) -> None:
    """Refuse *bam* unless its bytes are the ones the cache was scanned from.

    A manifest without an identity is refused before *bam* is looked at, so *bam* may be None then. A different
    size refuses at once. The same size and modification time accept without reading the file. Any other case is
    decided by the digest, so a copy under a new time still matches.
    """
    if not isinstance(recorded, dict) or not {"size", "mtime_ns", "blake2b"} <= recorded.keys():
        raise ScanCacheKeyError(
            "the cache records no identity for its BAM, only a path, and a path says nothing about the bytes. "
            "Re-scan."
        )
    path = Path(bam)
    stat = path.stat()
    if stat.st_size != recorded["size"]:
        raise ScanCacheKeyError(
            f"the BAM at {path} has {stat.st_size:,} bytes; the cache was scanned from {recorded['size']:,}. "
            f"Re-scan."
        )
    if stat.st_mtime_ns == recorded["mtime_ns"]:
        return
    if bam_identity(path)["blake2b"] != recorded["blake2b"]:
        raise ScanCacheKeyError(
            f"the BAM at {path} is not the one this cache was scanned from: its bytes differ. Re-scan."
        )


def _payload_field_types() -> dict[str, type]:
    """``AccumulatorPayload``'s annotations resolved to classes. ``dataclasses.fields()`` gives strings under
    ``from __future__ import annotations``, from which the nested banks cannot be rebuilt."""
    return typing.get_type_hints(AccumulatorPayload)


def _nested_dataclass(annotation) -> type | None:
    """The dataclass an annotation names, looking through unions; ``None`` if none.

    The union case matters: ``AccumulatorPayload.drain`` is ``DrainQC | None``, for which ``is_dataclass`` is False.
    """
    if dataclasses.is_dataclass(annotation):
        return annotation
    for argument in typing.get_args(annotation):
        if dataclasses.is_dataclass(argument):
            return argument
    return None


def _schema_names() -> list[str]:
    """Every name the cache is keyed by, nested banks included, in a stable order.

    The walk must stay fully recursive: ``DrainQC`` nests a ``GapCensus``, whose fields a one-level walk would miss.
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
    """Hash the schema the cached arrays were written under: the payload's fields and those of every bank nested in
    it (:func:`_schema_names`), each bank's column count, the strand model's layout (:func:`_strand_layout`) and the
    deposit rule's behaviour (:func:`deposit_digest`).

    No other key moves when the accumulator or the strand model changes. Without this, a cache written before a new
    payload field would fail in ``_payload_from_parts`` with a bare ``KeyError``, and one written before a new
    strand-model field would load it at its default.

    The column count is hashed as well as the name because :func:`_payload_from_parts` applies none of
    ``from_scan_result``'s shape checks, so a bank collapsed from ``[n, 2]`` to ``[n]`` would otherwise be accepted.
    It is read from the axis tables, so the digest needs no payload in hand.
    """
    return _digest(
        *(name.encode() for name in _schema_names()),
        *(shape.encode() for shape in _schema_shapes()),
        *(part.encode() for part in _strand_schema()),
        deposit_digest().encode(),
    )


def deposit_digest() -> str:
    """A hash of what the accumulator does: a fixed fragment set deposited into a fixed tiny partition, its banks
    hashed.

    A deposit-rule change moves no field name or column count, so without this a cache written under the old rule
    would serve old values to new code. One fragment set goes into one accumulator in one order, so the digest is
    stable across runs. The fixture is deliberately awkward (two annotated sj, a short region whose far boundary a
    fragment may or may not reach; contained, crossing, spliced and sj-only fragments) so that a change confined to
    one case still moves it.

    It runs the native accumulator, not the specification: the cache holds what the production scanner deposited.
    """
    # The sj CSR is keyed by region-bound index, not coordinate: 260 is bound 3, 1000 bound 4, 1120 bound 6 and
    # 2000 bound 7.
    region_bounds = np.array([0, 60, 200, 260, 1000, 1060, 1120, 2000, 2400], dtype=np.int64)
    accumulator = Accumulator(
        region_bounds=region_bounds,
        region_types=np.array(
            [
                RegionType.INTERGENIC,
                RegionType.EXON,
                RegionType.EXON,
                RegionType.INTRON,
                RegionType.EXON,
                RegionType.EXON,
                RegionType.INTRON,
                RegionType.EXON,
            ],
            dtype=np.uint8,
        ),
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
        (150, 1060, ((260, 1000),)),  # spliced; both blocks cross a boundary
        (210, 1040, ((260, 1000),)),  # spliced; neither block crosses a boundary
        (150, 2100, ((260, 1000), (1120, 2000))),  # two sj and boundaries
        (1030, 2050, ((1120, 2000),)),  # one sj, one boundary
        (100, 1060, ((100, 150), (260, 1000))),  # a leading intron, then an sj: each junction its own two ends
    ):
        # The native binding has no default for `hypotheses`; an empty set means the unspliced-only set.
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
    from rigel.scan_payload import BANK_AXES, N_STRAND_COLUMNS, SINGLE_COLUMN_AXES

    return [f"{name}:{N_STRAND_COLUMNS}" for name, _axis, _dtype in BANK_AXES] + [
        f"{name}:1" for name, _axis, _dtype in SINGLE_COLUMN_AXES
    ]


#: Settings that divide the scan's work, not its tally: the tally is identical at any worker count, so these stay
#: out of the key.
_WORK_ONLY_FIELDS = ("total_threads", "bgzf_threads")


def _scan_config_digest(scan_config) -> str:
    """The scan settings' part of the key: every field that can change a tally, and none that cannot."""
    fields = (
        dataclasses.asdict(scan_config)
        if dataclasses.is_dataclass(scan_config)
        else dict(scan_config)
    )
    tally = {k: v for k, v in fields.items() if k not in _WORK_ONLY_FIELDS}
    return _digest(json.dumps(tally, sort_keys=True, default=str).encode())


@dataclasses.dataclass(frozen=True, slots=True)
class ScanCache:
    """Everything one BAM scan produced that calibration consumes."""

    payload: AccumulatorPayload  # the tally — the expensive artifact
    strand_model: object  # StrandModels, including its per-sj table
    provenance: dict  # the key, the BAM, the scan config, the counts


# ── strand model round-trip ──────────────────────────────────────────────────────────────────────
# The per-sj table is stored, not only its 2x2 marginal, because the junction fit reads κ from the table.


def _strand_layout() -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """The strand half of the cache, read off the strand model's dataclasses: the sub-models, the int counters each
    stores as one int64 vector (in field order), and the per-sj table's columns.

    Derived rather than listed, so a new field is stored and moves :func:`payload_schema_digest` instead of loading
    at its default. A sub-model field that is neither an int counter nor ``sj_table`` has no place here and is refused.
    """
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels

    hints = typing.get_type_hints(StrandModel)
    counters = tuple(f.name for f in dataclasses.fields(StrandModel) if hints[f.name] is int)
    unstored = [
        f.name for f in dataclasses.fields(StrandModel) if f.name not in (*counters, "sj_table")
    ]
    if unstored:
        raise TypeError(
            f"StrandModel fields {unstored} are neither int counters nor sj_table, and the scan cache "
            f"stores a sub-model as exactly those, so it has nowhere to put them."
        )
    submodels = tuple(f.name for f in dataclasses.fields(StrandModels))
    columns = tuple(f.name for f in dataclasses.fields(SJStrandTable))
    return submodels, counters, columns


def _strand_schema() -> list[str]:
    """The strand half of the schema key: :func:`_strand_layout`, order included, since a sub-model's counters are
    read back by position."""
    submodels, counters, columns = _strand_layout()
    return [
        f"strand submodels:{','.join(submodels)}",
        f"strand counters:{','.join(counters)}",
        f"strand sj columns:{','.join(columns)}",
    ]


def _strand_arrays(strand_model) -> dict[str, np.ndarray]:
    submodels, counters, columns = _strand_layout()
    out: dict[str, np.ndarray] = {}
    for name in submodels:
        sub = getattr(strand_model, name)
        out[f"{name}__counts"] = np.array(
            [getattr(sub, field) for field in counters], dtype=np.int64
        )
        table = sub.sj_table
        out[f"{name}__has_table"] = np.array([table is not None], dtype=bool)
        if table is not None:
            for column in columns:
                out[f"{name}__sj_{column}"] = np.ascontiguousarray(getattr(table, column))
    return out


def _strand_from_arrays(data) -> object:
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels

    submodels, counters, columns = _strand_layout()
    built = {}
    for name in submodels:
        counts = data[f"{name}__counts"]
        table = None
        if bool(data[f"{name}__has_table"][0]):
            table = SJStrandTable(**{c: data[f"{name}__sj_{c}"] for c in columns})
        built[name] = StrandModel(**dict(zip(counters, (int(v) for v in counts))), sj_table=table)
    return StrandModels(**built)


def write_scan_cache(
    cache_dir: str | Path,
    *,
    payload: AccumulatorPayload,
    strand_model,
    index: TranscriptIndex,
    bam: str,
    scan_config,
) -> Path:
    """Persist one scan's calibration inputs under *cache_dir*, keyed to this index and config."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # The cache holds pass one, so the drain can re-run at any seed after loading.
    if payload.drain is not None:
        raise ValueError(
            "refusing to cache a DRAINED payload. The cache stores pass one so that the drain can be "
            "re-run at any seed; cache the payload `scan_and_buffer` returned and drain after loading."
        )

    # A nested bank's arrays go to the .npz as `field__sub` and its counters to the manifest: the manifest's
    # `json.dumps(default=str)` would write an array as a truncated repr.
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
        "transcript_model_digest": transcript_model_digest(index),
        # The scan settings; their part of the key is derived from this record when the cache is checked.
        "scan_config": dataclasses.asdict(scan_config)
        if dataclasses.is_dataclass(scan_config)
        else dict(scan_config),
        # ── provenance: what this cache is of ────────────────────────────────────────────────────
        "bam": bam_identity(bam),
        "payload_scalars": scalars,
    }
    (cache_dir / MANIFEST_JSON).write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str)
    )
    return cache_dir


def check_scan_cache(
    cache_dir: str | Path,
    index: TranscriptIndex,
    bam: str | Path | None = None,
    scan_config=None,
) -> None:
    """Refuse a cache, from its manifest alone, unless it describes this index and this BAM and, when *scan_config*
    is given, was scanned under those settings.

    These are the checks :func:`read_scan_cache` makes before it reads the arrays: the reach, the payload schema, the
    transcript models, the graph hash, the BAM and the scan settings. *bam* names the BAM the caller will use; without it
    the BAM at the recorded path is compared, if it is still there. A caller that holds settings should pass them,
    or a cache scanned under other settings reads as valid.
    """
    manifest = json.loads((Path(cache_dir) / MANIFEST_JSON).read_text())
    _check_manifest(manifest, index, bam, scan_config)


def _check_manifest(
    manifest: dict, index: TranscriptIndex, bam: str | Path | None, scan_config=None
) -> None:
    """Every refusal :func:`check_scan_cache` makes, in order: the index's key, the scan settings, the BAM."""
    _check_index_key(manifest, index)
    if scan_config is not None:
        _check_scan_settings(manifest, scan_config)
    _check_recorded_bam(manifest, bam)


def _check_index_key(manifest: dict, index: TranscriptIndex) -> None:
    """The reach, the payload schema, the transcript models and the graph hash the cache was keyed by."""
    expected_reach = reach_digest(index)
    if manifest["reach_digest"] != expected_reach:
        raise ScanCacheKeyError(
            f"cache reach digest {manifest['reach_digest']} != index reach digest {expected_reach}. "
            f"The boundary reaches moved. Neither partition_hash nor graph_hash covers reach, so this is "
            f"the only check that notices: a rebuild can move the splice-junction reaches "
            f"with both of those byte-identical. Re-scan against this index."
        )

    expected_schema = payload_schema_digest()
    if manifest.get("payload_schema_digest") != expected_schema:
        raise ScanCacheKeyError(
            f"cache payload_schema_digest {manifest.get('payload_schema_digest')!r} != "
            f"{expected_schema!r}. The ACCUMULATOR's or the STRAND MODEL's schema moved, so these arrays "
            f"are not the fields this build reads — a missing one would otherwise surface as a bare "
            f"KeyError far from here, or load at its default. Re-scan; nothing derivable from the index "
            f"needs rebuilding."
        )

    expected_models = transcript_model_digest(index)
    if manifest.get("transcript_model_digest") != expected_models:
        raise ScanCacheKeyError(
            f"cache transcript model digest {manifest.get('transcript_model_digest')!r} != index "
            f"{expected_models!r}. Resolution reads which transcripts hold a fragment's blocks and junctions, and the "
            f"graph and the reaches can stay byte-identical when the transcript models move. Re-scan against this "
            f"index."
        )

    if manifest["payload_scalars"]["graph_hash"] != manifest["graph_hash"]:
        raise ScanCacheKeyError(
            "cache manifest graph_hash disagrees with its own payload's graph_hash"
        )
    if manifest["graph_hash"] != index.graph_hash:
        raise ScanCacheKeyError(
            f"cache graph_hash {manifest['graph_hash']} != index graph_hash {index.graph_hash}. The region "
            f"partition or the sj CSR moved; this tally does not describe this index."
        )


def _check_scan_settings(manifest: dict, scan_config) -> None:
    """The tally settings the cache was scanned under. Both digests are computed here, from the recorded settings and
    the caller's, so a thread count on either side decides nothing (`_WORK_ONLY_FIELDS`)."""
    recorded_settings = _scan_config_digest(manifest["scan_config"])
    expected_settings = _scan_config_digest(scan_config)
    if recorded_settings != expected_settings:
        raise ScanCacheKeyError(
            f"cache scan settings digest {recorded_settings} != {expected_settings}. Two scans of one BAM "
            f"under different tally settings are different tallies."
        )


def _check_recorded_bam(manifest: dict, bam: str | Path | None) -> None:
    """The BAM's bytes against the identity the cache recorded. Without *bam*, the BAM at the recorded path is
    compared when it is still there."""
    recorded = manifest.get("bam")
    if bam is None and isinstance(recorded, dict):
        if not os.path.exists(recorded["path"]):
            logger.info(
                "the BAM this cache was scanned from is not at %s, so its identity is not checked",
                recorded["path"],
            )
            return
        bam = recorded["path"]
    check_bam_identity(recorded, bam)


def read_scan_cache(
    cache_dir: str | Path, index: TranscriptIndex, scan_config=None, bam: str | Path | None = None
) -> ScanCache:
    """Load a cache, refusing it as :func:`check_scan_cache` does.

    Pass *scan_config* whenever you have one: without it a cache scanned under other settings loads silently.
    """
    cache_dir = Path(cache_dir)
    manifest = json.loads((cache_dir / MANIFEST_JSON).read_text())
    _check_manifest(manifest, index, bam, scan_config)
    payload_scalars = manifest["payload_scalars"]

    with np.load(cache_dir / PAYLOAD_NPZ) as data:
        arrays = {name: data[name] for name in data.files}
    with np.load(cache_dir / STRAND_NPZ) as data:
        strand = {name: data[name] for name in data.files}
    payload = _payload_from_parts(arrays, payload_scalars)
    if payload.graph_hash != manifest["graph_hash"]:
        raise ScanCacheKeyError(
            f"cache graph_hash {payload.graph_hash} != manifest graph_hash {manifest['graph_hash']}"
        )

    return ScanCache(
        payload=payload,
        strand_model=_strand_from_arrays(strand),
        provenance=manifest,
    )


def _payload_from_parts(arrays: dict, scalars: dict) -> AccumulatorPayload:
    """Rebuild the payload from the ``.npz`` arrays and the manifest's scalars.

    Generic over the nested banks: a sub-field comes from the ``.npz`` when it is an array and from the manifest
    otherwise, so a bank that grows an array needs no edit here.
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


def index_derived_inputs(index: TranscriptIndex) -> dict:
    """The calibrate() arguments that come from the index, rebuilt every time rather than cached: they rebuild in
    under a second, and a stored copy could go stale against the index."""
    return {
        "region_arrays": RegionArrays.from_index(index),
        "boundary_flags": build_boundary_flags_array(index),
        # Required: `calibrate` refuses an sj axis whose length disagrees with the payload's `n_sj`.
        "sj": build_sj_geometry_arrays(index),
    }


def calibration_inputs(
    cache: ScanCache,
    index: TranscriptIndex,
    *,
    lift_out: dict | None = None,
) -> dict:
    """The keyword arguments `calibrate` needs, in production's frame.

    The payload is drained here, as production drains it before ``build_fl_models`` and ``calibrate``, and the fl
    models are production's own (`pipeline.library_fl_models`) on the drained payload, drained at production's
    ``second_pass_seed``. *lift_out*, a dict, receives the drain's ``lift`` box ``(undrained, choices, region_types,
    sj)``, which `_oracle.lift_drain_parts` needs to put cached oracle partitions in the same frame; it
    stays empty when the side buffer held nothing.
    """
    lift = {} if lift_out is None else lift_out
    payload = drain_side_buffer(
        cache.payload, index, cache.strand_model, seed=PipelineConfig().second_pass_seed, lift=lift
    )
    fl_models = library_fl_models(payload, index)
    return {
        "payload": payload,
        "strand_model": cache.strand_model,
        "gdna_fl_pmf": fl_models.gdna_pmf,
        "rna_fl_pmf": fl_models.rna_pmf,
        **index_derived_inputs(index),
    }
