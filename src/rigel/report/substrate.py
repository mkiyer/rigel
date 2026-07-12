"""Load the ``rigel quant`` report substrate into a structured bundle.

The substrate is everything ``rigel quant`` writes to its output directory:

* ``summary.json`` — the lean run manifest (schema v2), read eagerly.
* companion feather tables (``fragment_lengths``, ``gene_quant``,
  ``calibration_track``, …) — read **lazily** on first access, so a report only
  pays I/O + memory for the tables it actually uses. The big per-transcript /
  nascent tables (``quant``, ``nrna_quant``) are exposed for future features but
  never read unless something asks for them.

Only ``summary.json`` is required; every companion table is optional so the
report degrades gracefully on partial or older outputs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import pandas as pd

#: summary.json schema this builder targets. Reports still build on other
#: versions, but a mismatch is surfaced as a warning banner.
EXPECTED_SCHEMA_VERSION = 2


class SubstrateError(Exception):
    """Raised when the report substrate is missing or unreadable."""


@dataclass
class ReportSubstrate:
    """Everything the report builder needs from a quant output directory.

    ``summary`` is loaded up front; the feather tables are lazy ``cached_property``
    accessors that read from ``output_dir`` on first use and return ``None`` when
    the file is absent.
    """

    output_dir: Path
    summary: dict
    #: Non-fatal issues to surface in the report (schema drift, missing tables).
    warnings: list[str] = field(default_factory=list)

    def _feather(self, name: str) -> pd.DataFrame | None:
        path = self.output_dir / name
        if not path.exists():
            return None
        try:
            return pd.read_feather(path)
        except Exception as exc:  # pragma: no cover - corrupt file is unusual
            raise SubstrateError(f"Could not read {path.name}: {exc}") from exc

    # -- tables the report currently consumes --
    @cached_property
    def fragment_lengths(self) -> pd.DataFrame | None:
        return self._feather("fragment_lengths.feather")

    @cached_property
    def gene_quant(self) -> pd.DataFrame | None:
        return self._feather("gene_quant.feather")

    @cached_property
    def calibration_track(self) -> pd.DataFrame | None:
        return self._feather("calibration_track.feather")

    # -- available for future panels; not read unless accessed --
    @cached_property
    def quant(self) -> pd.DataFrame | None:
        return self._feather("quant.feather")

    @cached_property
    def nrna_quant(self) -> pd.DataFrame | None:
        return self._feather("nrna_quant.feather")

    @cached_property
    def loci(self) -> pd.DataFrame | None:
        return self._feather("loci.feather")

    @property
    def schema_version(self) -> int | None:
        return self.summary.get("schema_version")

    @property
    def sample_name(self) -> str:
        """A human label for the sample, derived from the input BAM stem."""
        bam = self.summary.get("input", {}).get("bam_file")
        return Path(bam).stem if bam else self.output_dir.name


def load_substrate(output_dir: str | Path) -> ReportSubstrate:
    """Load the report substrate from a ``rigel quant`` output directory.

    Reads ``summary.json`` and checks which companion tables are present (for
    warnings); the tables themselves load lazily on first access.

    Raises
    ------
    SubstrateError
        If the directory or ``summary.json`` is missing / unreadable.
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise SubstrateError(f"Not a directory: {output_dir}")

    summary_path = output_dir / "summary.json"
    if not summary_path.exists():
        raise SubstrateError(
            f"No summary.json in {output_dir}. Point `rigel report` at a "
            f"`rigel quant` output directory."
        )
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SubstrateError(f"Could not parse {summary_path}: {exc}") from exc

    warnings: list[str] = []
    sv = summary.get("schema_version")
    if sv is None:
        warnings.append(
            "summary.json has no schema_version (pre-v2 run); fragment-length "
            "distributions may be unavailable."
        )
    elif sv != EXPECTED_SCHEMA_VERSION:
        warnings.append(
            f"summary.json schema_version is {sv}; this builder targets "
            f"v{EXPECTED_SCHEMA_VERSION}. Some panels may be incomplete."
        )
    if not (output_dir / "fragment_lengths.feather").exists():
        warnings.append(
            "fragment_lengths.feather not found; the fragment-length "
            "distribution charts will be omitted."
        )

    return ReportSubstrate(output_dir=output_dir, summary=summary, warnings=warnings)
