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
FL histograms                the scan                                    **yes, RAW**
``region_arrays``            ``RegionArrays.from_index``   (0.11 s)      no
``boundary_flags``           ``build_boundary_flags_array``(0.04 s)      no
``gdna_fl_pmf``/``rna_fl_pmf``  ``build_fl_models(...)``                 no — derived
``config``                   the thing you are varying                   no
``injected_priors``          fitted BY ``calibrate``                     no — see below
===========================  ==========================================  ==========

⚠ **Anything derivable from the index is rebuilt on load, never stored.** 0.15 s against an 8.45 s index
load that happens anyway, and a stored copy is how a cache goes stale against the thing it describes
(`CARRY_FORWARD.md` §3 trap 25).

⭐ **The FL histograms are cached RAW, not as the derived pmfs.** `build_fl_models` stays the single
source of truth; freezing its output would mean a change to the FL model silently does not reach a
cached scan.

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

⛔ **STEP 4 (the population-prior seed) IS NOT HERE, AND THAT IS DELIBERATE.** Seeding a toy from a
genome-scale scan needs ``InjectedCalibrationPriors``, which `calibrate` fits and stashes in
``_debug["calibration_priors"]`` — so extracting it requires `calibrate` to run, and it cannot until S5
rewires the consumers. `tests/test_scan_cache.py` carries a strict-xfail that documents the dependency
rather than pretending the seed path is verified.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
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
    "index_derived_inputs",
    "payload_schema_digest",
    "reach_digest",
    "read_scan_cache",
    "write_scan_cache",
]

MANIFEST_JSON = "manifest.json"
PAYLOAD_NPZ = "payload.npz"
STRAND_NPZ = "strand.npz"
FL_NPZ = "fl.npz"

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


def payload_schema_digest() -> str:
    """Hash ``AccumulatorPayload``'s own field list — the schema the cached arrays were written under.

    ⛔ **No other key covers this, and the gap is not hypothetical.** ``graph_hash`` describes the index,
    ``reach_digest`` the reaches, ``scan_config_digest`` the scan settings — none of them changes when the
    ACCUMULATOR changes. S5.a added ``length_sum`` to every population, and without this key a cache
    written the day before would have been accepted and then failed deep inside ``_payload_from_parts``
    with a bare ``KeyError``, which reads as a bug in the cache rather than as a stale cache.

    Field names only: a dtype change is already caught at load by ``_bank``'s assertion, and names are
    what the ``.npz`` is keyed by.
    """
    return _digest(*(field.name.encode() for field in dataclasses.fields(AccumulatorPayload)))


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
    fl_global_counts: np.ndarray  # int64[max_size + 1] — raw, not a pmf
    fl_rna_counts: np.ndarray  # int64[max_size + 1] — SPLICED_ANNOT, raw
    fl_max_size: int
    provenance: dict  # the key, the BAM, the scan config, the counts


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
    frag_length_models,
    index: "TranscriptIndex",
    bam: str,
    scan_config,
) -> Path:
    """Persist one scan's calibration inputs under *cache_dir*, keyed to this index and config."""
    from .splice import SpliceType

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    arrays: dict[str, np.ndarray] = {}
    scalars: dict[str, object] = {}
    for field in dataclasses.fields(AccumulatorPayload):
        value = getattr(payload, field.name)
        if isinstance(value, np.ndarray):
            arrays[field.name] = np.ascontiguousarray(value)
        elif dataclasses.is_dataclass(value):
            scalars[field.name] = dataclasses.asdict(value)
        else:
            scalars[field.name] = value
    np.savez_compressed(cache_dir / PAYLOAD_NPZ, **arrays)
    np.savez_compressed(cache_dir / STRAND_NPZ, **_strand_arrays(strand_model))

    global_counts = np.ascontiguousarray(frag_length_models.global_model.counts)
    rna_counts = np.ascontiguousarray(
        frag_length_models.category_models[SpliceType.SPLICED_ANNOT].counts
    )
    np.savez_compressed(cache_dir / FL_NPZ, global_counts=global_counts, rna_counts=rna_counts)

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
        "fl_max_size": int(frag_length_models.max_size),
        "payload_scalars": scalars,
    }
    (cache_dir / MANIFEST_JSON).write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str)
    )
    return cache_dir


def read_scan_cache(cache_dir: str | Path, index: "TranscriptIndex") -> ScanCache:
    """Load a cache and REFUSE it unless it describes this index and this scan configuration."""
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
    with np.load(cache_dir / FL_NPZ) as data:
        fl_global, fl_rna = data["global_counts"], data["rna_counts"]

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

    return ScanCache(
        payload=payload,
        strand_model=_strand_from_arrays(strand),
        fl_global_counts=fl_global,
        fl_rna_counts=fl_rna,
        fl_max_size=int(manifest["fl_max_size"]),
        provenance=manifest,
    )


def _payload_from_parts(arrays: dict, scalars: dict) -> AccumulatorPayload:
    from .scan_payload import ScanQC

    kwargs: dict[str, object] = {}
    for field in dataclasses.fields(AccumulatorPayload):
        if field.name in arrays:
            kwargs[field.name] = arrays[field.name]
        elif field.name == "qc":
            kwargs[field.name] = ScanQC(**scalars[field.name])
        else:
            kwargs[field.name] = scalars[field.name]
    return AccumulatorPayload(**kwargs)


def check_scan_config(cache: ScanCache, scan_config) -> None:
    """Refuse a cache produced under a different scan configuration."""
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
    from .calibration.splice_graph import build_boundary_flags_array

    return {
        "region_arrays": RegionArrays.from_index(index),
        "boundary_flags": build_boundary_flags_array(index),
    }


def calibration_inputs(cache: ScanCache, index: "TranscriptIndex") -> dict:
    """Exactly the keyword arguments `calibrate` needs, with the index-derived ones REBUILT.

    ⭐ **Both component histograms come from the PAYLOAD's five pure pools**, not from the scanner's
    category models. One quantity, one source: the scanner's spliced histogram is transcript-space and
    requires a UNIQUE transcript, while the accumulator's `RNA_SPLICED` pool is a structural rule over a
    larger population and already excludes `sj_implicit` fragments. Only the unconditional `global`
    histogram still comes from the scan, because no pool is unconditional.
    """
    from .calibration.fl import build_fl_models, gdna_fl_mass, rna_fl_mass

    fl_models = build_fl_models(
        global_counts=cache.fl_global_counts,
        rna_counts=rna_fl_mass(cache.payload),
        gdna_counts=gdna_fl_mass(cache.payload),
        max_size=cache.fl_max_size,
    )
    return {
        "payload": cache.payload,
        "strand_model": cache.strand_model,
        "gdna_fl_pmf": fl_models.gdna_pmf,
        "rna_fl_pmf": fl_models.rna_pmf,
        **index_derived_inputs(index),
    }
