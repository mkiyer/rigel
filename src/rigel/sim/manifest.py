"""Manifest helpers shared by simulation generation and analysis."""

from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "condition_dir_name",
    "condition_manifest_map",
    "gdna_label_for_rate",
    "load_manifest",
    "write_manifest",
]


def condition_dir_name(
    gdna_label: str,
    strand_specificity: float,
    nrna_label: str,
    capture_label: str | None = None,
) -> str:
    """Return the standard synthetic-suite condition directory name."""
    name = f"gdna_{gdna_label}_ss_{strand_specificity:.2f}_nrna_{nrna_label}"
    if capture_label is not None:
        name = f"{name}_capture_{capture_label}"
    return name


def gdna_label_for_rate(rate: float, labels: list[str] | None, index: int) -> str:
    """Return a stable gDNA label for a configured rate."""
    if labels is not None and index < len(labels):
        return labels[index]
    return f"r{rate:g}"


def _jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_manifest(outdir: Path, config: Any, conditions: list[dict[str, Any]]) -> Path:
    """Write a simulation manifest and return its path."""
    truth_abundances = "truth_abundances.tsv"
    for condition in conditions:
        if condition.get("truth_abundances"):
            truth_abundances = str(condition["truth_abundances"])
            break
    manifest = {
        "version": 1,
        "simulation": _jsonable(getattr(config, "simulation", {})),
        "gdna": _jsonable(getattr(config, "gdna", {})),
        "nrna": _jsonable(getattr(config, "nrna", {})),
        "capture": _jsonable(getattr(config, "capture", {})),
        "capture_configs": _jsonable(getattr(config, "capture_configs", [])),
        "abundance": _jsonable(getattr(config, "abundance", {})),
        "truth_abundances": truth_abundances,
        "conditions": _jsonable(conditions),
    }
    for key in ("genome", "gtf", "transcript_filter", "strand_specificities"):
        if hasattr(config, key):
            value = getattr(config, key)
            if key in {"genome", "gtf"} and value:
                value = str(Path(value).resolve())
            manifest[key] = _jsonable(value)
    path = outdir / "manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path


def load_manifest(path_or_dir: Path) -> dict[str, Any]:
    """Load a manifest from a file path or a simulation output directory."""
    path = path_or_dir / "manifest.json" if path_or_dir.is_dir() else path_or_dir
    if not path.exists():
        return {}
    with open(path) as handle:
        return json.load(handle)


def condition_manifest_map(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return condition metadata keyed by condition name."""
    conditions = manifest.get("conditions", []) if isinstance(manifest, dict) else []
    return {
        str(condition.get("name")): condition
        for condition in conditions
        if isinstance(condition, dict) and condition.get("name")
    }