"""rigel.scan_cache — scan once, calibrate many times.

    TODO item 2 (the cached substrate)   ·   `docs/testing/testing_plan.md`   ·   Gate: `tests/test_scan_cache.py`

⭐ **WHY.** Scanning is the expensive step and calibration is the one under development. On a real cfRNA
run: index load ~8 s, BAM scan ~2 s, **calibration ~66 s** (`CARRY_FORWARD.md` §1 fact 22) — and a 5 M
fragment simulated condition costs far more than that to scan. Caching the scan's output took a
24-condition sweep from ~13 min to ~9 s on the old path. This is that, rebuilt against the S4 payload.

WHAT IS STORED, AND WHAT DELIBERATELY IS NOT
--------------------------------------------
`calibrate()` takes eight things. They come from three different places, and only one of those is
expensive:

===========================  ==========================================  ==========
input                        origin                                      cached?
===========================  ==========================================  ==========
``payload``                  the scan                                    **yes**
``strand_model``             the scan                                    **yes**
``region_arrays``            ``RegionArrays.from_index``   (0.11 s)      no
``edge_flags``               ``build_edge_flags_array``    (0.04 s)      no
``junctions``                ``build_junction_geometry_arrays``          no
``gdna_fl_pmf``/``rna_fl_pmf``  ``build_fl_models(payload)``             no — derived
``config``                   the thing you are varying                   no
``injected_priors``          fitted BY ``calibrate``                     no — see below
===========================  ==========================================  ==========

⚠ **Anything derivable from the index is rebuilt on load, never stored.** 0.15 s against an 8.45 s index
load that happens anyway, and a stored copy is how a cache goes stale against the thing it describes
(`CARRY_FORWARD.md` §3 trap 25).

⭐ **There is no separate FL row any more, and that is C2.** Every fragment-length histogram — the two
pure pools and the unconditional anchor they are EB-shrunk toward — is a field OF the payload, so
caching the payload caches them, in one frame, by construction. The scanner's own histogram used to be
cached alongside it purely to serve as that anchor; `docs/FRAGMENT_LENGTH_AUDIT.md` D1–D3 is what that
cost. `build_fl_models` remains the single source of truth for the derived pmfs, which are still not
cached — freezing its output would mean a change to the FL model silently does not reach a cached scan.

THE KEY NEEDS THREE PARTS, AND THE MIDDLE ONE IS A GAP THAT WAS ALREADY LOGGED
-----------------------------------------------------------------------------
* ``graph_hash`` — the node partition plus the junction CSR. The payload already carries it.
* ⭐ **a PAYLOAD-SCHEMA digest.** The accumulator's own field list. None of the other keys moves when
  the accumulator changes, so without it a cache written before an accumulator change is accepted and
  then fails deep inside the loader with a bare ``KeyError``. S5.a made that concrete by adding
  ``length_sum`` to every population.
* ⭐ **a REACH digest.** `TODO.md`: ``reach`` is consumed by calibration and covered by **neither**
  ``partition_hash`` **nor** ``graph_hash`` — correctly, since neither the scan nor the accumulator reads
  it — and the gap "becomes live the moment something caches a calibration output". A cache loaded
  against an index is that moment. The 2026-07-30 rebuild moved ~38 % of human contiguous reaches while
  **both** existing hashes stayed byte-identical, so a reach-blind key would verify clean against a
  moved index.
* **the scan config**, because two scans of one BAM under different settings are different tallies.

⚠ **Not pickle.** The old `scripts/debug/calib_cache.py` pickled `calibrate()`'s kwargs; a pickle of
numpy-holding dataclasses is fragile exactly across schema changes, and the schema changes at S5 — a
pickle written today would not load after it. Arrays go to ``.npz``, scalars and provenance to JSON.

✅ **STEP 4 (the population-prior seed) IS LIVE as of S5.f/S6.** Seeding a toy from a genome-scale scan
needs ``InjectedCalibrationPriors``, which `calibrate` fits and stashes in
``_debug["calibration_priors"]`` — so it required `calibrate` to run, and it does.
``test_population_priors_can_be_extracted_from_a_cached_scan`` was a **strict** xfail naming S5 as the
blocker and is now a live test.

⛔ **It was still xfailing at S5.f for a DIFFERENT reason than the one recorded** — ``calibration_inputs``
had not been given ``junctions``, which `calibrate` gained at S5.f. A strict xfail proves only that
something fails, never that the recorded cause is still the cause; that is why the reason has to be
re-read when the blocker lifts rather than assumed.
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
    """Hash the edge reach columns — the part of the index calibration reads and no other key covers.

    ⚠ Computed on demand, never stored beside the index (`CARRY_FORWARD.md` §3 trap 25). It is stored in
    the CACHE's manifest, which is a different thing: the cache is describing an external artifact it
    cannot recompute, exactly as `index.source_record` does for the FASTA and the GTF.
    """
    edges = index.edges_df
    if edges is None:
        raise ScanCacheKeyError("index has no edges_df loaded; cannot compute the reach digest")
    missing = [column for column in REACH_COLUMNS if column not in edges.columns]
    if missing:
        raise ScanCacheKeyError(f"edges.feather is missing reach columns {missing}")
    parts: list[bytes] = []
    for column in REACH_COLUMNS:
        values = np.ascontiguousarray(edges[column].to_numpy())
        parts.append(column.encode())
        parts.append(str(values.dtype).encode())
        parts.append(values.tobytes())
    return _digest(*parts)


def _payload_field_types() -> dict[str, type]:
    """``AccumulatorPayload``'s annotations resolved to real classes.

    ⚠ ``dataclasses.fields()`` hands back annotation STRINGS under ``from __future__ import annotations``,
    so the nested banks cannot be reconstructed from them. Resolving the hints is what lets the read path
    stay generic instead of carrying a name→class table that a new bank could quietly fall out of.
    """
    return typing.get_type_hints(AccumulatorPayload)


def _nested_dataclass(annotation) -> type | None:
    """The dataclass an annotation names, **looking through ``Optional`` / unions**. ``None`` if none.

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

    ⭐ This is what makes the digest below cover what it claims to. ⚠ **Fully recursive, not one level.**
    It covered exactly one level until 2026-08-02, which was enough for ``ScanQC`` / ``GapCensus`` /
    ``DeferredFragments`` and stopped being enough the moment ``DrainQC`` nested a ``GapCensus`` inside
    itself — the same invisibility X8 was filed for, one level down.
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
    """Hash the schema the cached arrays were written under — ``AccumulatorPayload``'s field list **and
    the fields of every bank nested inside it**.

    ⛔ **No other key covers this, and the gap is not hypothetical.** ``graph_hash`` describes the index,
    ``reach_digest`` the reaches, ``scan_config_digest`` the scan settings — none of them changes when the
    ACCUMULATOR changes. S5.a added ``length_sum`` to every population, and without this key a cache
    written the day before would have been accepted and then failed deep inside ``_payload_from_parts``
    with a bare ``KeyError``, which reads as a bug in the cache rather than as a stale cache.

    ⭐ **THE NESTING IS WHY THIS RECURSES, and it was a real defect.** The digest used to hash the top-level
    field names ALONE, so a change *inside* ``ScanQC`` was invisible to it: a renamed qc field let a stale
    cache be accepted by the key and then fail deep in the loader with a bare ``TypeError`` — precisely the
    failure mode this digest exists to prevent. S1 made that worse rather than better, because
    ``DeferredFragments`` puts **thirteen** array names inside one field and every one of them is an
    ``.npz`` key. ⚠ It recursed exactly ONE level until 2026-08-02, when ``DrainQC`` nested a ``GapCensus``
    inside itself and put the same invisibility one level down; :func:`_schema_names` is fully recursive
    now, and looks through ``Optional``.

    Field names only: a dtype change is already caught at load by ``_bank``'s and
    ``DeferredFragments.from_dict``'s assertions, and names are what the ``.npz`` is keyed by.
    """
    return _digest(*(name.encode() for name in _schema_names()))


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
    strand_model: object  # StrandModels, including its per-junction table
    provenance: dict  # the key, the BAM, the scan config, the counts

    # ⛔ `fl_global_counts`, `fl_rna_counts` and `fl_max_size` were DELETED by C2
    # (docs/FRAGMENT_LENGTH_AUDIT.md). The first was the scanner's histogram, cached so it could be
    # the empirical-Bayes anchor — the frame mismatch itself, persisted. The second was **D5**: a
    # field named as if it were live, written and read back and consumed by nothing. The third
    # duplicated `payload.max_length`.
    #
    # ⭐ Every fragment-length histogram now comes off `payload`, so there is nothing left to cache
    # beside it. ⚠ `fl.npz` is no longer written or read; caches written before C2 still load, since
    # an extra file on disk is not a key.


# ── strand model round-trip ──────────────────────────────────────────────────────────────────────
# ⚠ The 2x2 is the MARGINAL of the per-junction table, and the strand OVERDISPERSION is fitted from the
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

    # ⛔ A NESTED BANK'S ARRAYS GO TO THE .npz, NOT TO THE MANIFEST. `dataclasses.asdict` on
    # `DeferredFragments` yields a dict of ndarrays, and the manifest is written with
    # `json.dumps(..., default=str)` — which would stringify each array to a TRUNCATED repr, silently, and
    # read back as text. Nested dataclasses are therefore split by the type of each sub-field: arrays are
    # prefixed `field__sub` and counters stay scalars.
    # ⛔ THE CACHE HOLDS A SCAN, AND A DRAINED PAYLOAD IS NOT ONE. The second pass runs *after* the cache
    # is read (`SPEC_SECOND_PASS.md` §2), which is what lets one scan be drained repeatedly at different
    # seeds without re-reading the BAM — the property P5 and P6 both need. Writing a drained payload would
    # bake one draw into the cache, and it would also serialise `DrainQC.census_before` through
    # `json.dumps(default=str)` as a stringified repr, silently, which is X9's defect one level down.
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

    ⚠ **Pass ``scan_config`` whenever you have one.** Without it this checks only that the manifest is
    consistent with ITSELF; a cache scanned under different settings is a different tally and will load
    silently. The check is :func:`check_scan_config`, which existed and had **no caller at all** until
    S6 — the docstring promised a guarantee the code did not provide (`CARRY_FORWARD.md` §3 trap 27).
    """
    cache_dir = Path(cache_dir)
    manifest = json.loads((cache_dir / MANIFEST_JSON).read_text())

    expected_reach = reach_digest(index)
    if manifest["reach_digest"] != expected_reach:
        raise ScanCacheKeyError(
            f"cache reach digest {manifest['reach_digest']} != index reach digest {expected_reach}. "
            f"The edge reaches moved. ⚠ Neither partition_hash nor graph_hash covers reach, so this is "
            f"the only check that notices — a 2026-07-30 rebuild moved ~38 % of human contiguous "
            f"reaches with both of those byte-identical. Re-scan against this index."
        )

    # ⚠ Self-consistency: the manifest records the scan config AND its digest, so a tampered or
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
            f"cache graph_hash {payload.graph_hash} != index graph_hash {expected_graph}. The node "
            f"partition or the junction CSR moved; this tally does not describe this index."
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

    ⭐ Generic over the nested banks: each one's sub-fields are taken from the ``.npz`` when they are arrays
    and from the manifest when they are counters, so a bank that grows an array joins the round trip with no
    edit here. ⚠ And it grows the ``payload_schema_digest`` at the same time, which is what refuses a cache
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

    ⚠ Called by :func:`read_scan_cache` when it is given a ``scan_config``. It had **no caller at all**
    until S6 while ``read_scan_cache``'s docstring claimed the guarantee it provides — two statements
    about one contract, disagreeing (`CARRY_FORWARD.md` §3 trap 27). Two scans of one BAM under
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

    ⚠ Deliberately not cached: 0.15 s together, against an 8.45 s index load that happens anyway. A
    stored copy is how a cache goes stale against the thing it describes (`CARRY_FORWARD.md` §3 trap 25).
    """
    from .calibration.region_arrays import RegionArrays
    from .calibration.splice_graph import (
        build_edge_flags_array,
        build_junction_geometry_arrays,
    )

    return {
        "region_arrays": RegionArrays.from_index(index),
        "edge_flags": build_edge_flags_array(index),
        # ⚠ The JUNCTION axis is index-derived too, and it is not optional: `calibrate` refuses an axis
        # whose length disagrees with the payload's `n_sj`, because one addressing a different graph
        # would place every splice on the wrong line. Omitting it here was an S5.f miss that only the
        # guard caught.
        "junctions": build_junction_geometry_arrays(index),
    }


def calibration_inputs(cache: ScanCache, index: "TranscriptIndex") -> dict:
    """Exactly the keyword arguments `calibrate` needs, with the index-derived ones REBUILT.

    ⭐ **Every fragment-length histogram comes from the PAYLOAD** — the two pure pools *and*, since
    C2.1, the unconditional anchor they are shrunk toward. One quantity, one source, one frame: the
    scanner's spliced histogram is transcript-space and requires a UNIQUE transcript, while the
    accumulator's `RNA_SPLICED` pool is a structural rule over a larger population and already
    excludes `sj_implicit` fragments; and the anchor is `deposited_lengths`, binned at the same `L`.
    """
    from .calibration.fl import build_fl_models

    fl_models = build_fl_models(cache.payload)
    return {
        "payload": cache.payload,
        "strand_model": cache.strand_model,
        "gdna_fl_pmf": fl_models.gdna_pmf,
        "rna_fl_pmf": fl_models.rna_pmf,
        **index_derived_inputs(index),
    }
