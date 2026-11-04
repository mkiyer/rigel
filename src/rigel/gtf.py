"""rigel.gtf — a streaming GTF/GFF2 parser.

Reads GENCODE-style GTF files, plain or gzipped, into ``GTFRecord`` dataclasses for ``rigel.index``'s build.
"""

import gzip
import logging
import re
from collections.abc import Generator, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

logger = logging.getLogger(__name__)


AttrValue = str | list[str]


@dataclass(slots=True)
class GTFRecord:
    """One GTF line, its coordinates converted to 0-based half-open."""

    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: float | None
    strand: str
    phase: str
    attrs: dict[str, AttrValue] = field(default_factory=dict)
    tags: set[str] = field(default_factory=set)

    # `key "value";` or `key value;`, the last attribute's `;` optional, as many files omit it.
    # Group 1: key; group 2: quoted value; group 3: unquoted value.
    ATTR_TOKEN_PATTERN = re.compile(r'\s*([^";\s]+)\s+(?:"([^"]*)"|([^;]+))\s*(?:;|$)')

    @classmethod
    def from_str(cls, s: str) -> "GTFRecord":
        """Parse one GTF line."""
        fields = s.rstrip("\n").split("\t")
        if len(fields) < 9:
            raise ValueError(f"Malformed GTF line, expected 9 columns: {s}")

        # GTF is 1-based inclusive; convert to 0-based half-open
        try:
            start_pos = int(fields[3]) - 1
            end_pos = int(fields[4])
        except ValueError:
            raise ValueError(f"Non-integer coordinates in line: {fields}")
        if end_pos <= start_pos:
            raise ValueError(f"{fields[2]} ends ({fields[4]}) before it starts ({fields[3]})")

        # score: '.' -> None
        score_val = None
        if fields[5] != ".":
            try:
                score_val = float(fields[5])
            except ValueError as exc:
                raise ValueError(f"Non-numeric score in line: {fields[5]}") from exc

        # phase: 0, 1, 2 or '.'
        phase_val = fields[7]
        if phase_val != "." and phase_val not in {"0", "1", "2"}:
            raise ValueError(f"Invalid phase value: {phase_val}")

        parsed_attrs, parsed_tags = cls._parse_attrs(fields[8])

        return cls(
            seqname=fields[0],
            source=fields[1],
            feature=fields[2],
            start=start_pos,
            end=end_pos,
            score=score_val,
            strand=fields[6],
            phase=phase_val,
            attrs=parsed_attrs,
            tags=parsed_tags,
        )

    def __str__(self) -> str:
        """The record as a GTF line, 1-based."""
        score_str = "." if self.score is None else str(self.score)

        attr_parts = []
        for k, v in self.attrs.items():
            if isinstance(v, list):
                for vv in v:
                    attr_parts.append(f'{k} "{vv}";')
            else:
                attr_parts.append(f'{k} "{v}";')
        for k in self.tags:
            attr_parts.append(f'tag "{k}";')
        attr_field = " ".join(attr_parts)

        fields = [
            self.seqname,
            self.source,
            self.feature,
            str(self.start + 1),  # back to 1-based
            str(self.end),
            score_str,
            self.strand,
            self.phase,
            attr_field,
        ]
        return "\t".join(fields)

    @staticmethod
    def parse(
        line_iter: Iterable[str],
        *,
        parse_mode: Literal["strict", "warn-skip"] = "strict",
    ) -> Generator["GTFRecord", None, None]:
        """Parse GTF lines from an iterable, skipping blank and comment lines.

        Parameters
        ----------
        line_iter : Iterable[str]
            Input GTF lines.
        parse_mode : {"strict", "warn-skip"}
            - ``"strict"``: raise ``ValueError`` on malformed lines.
            - ``"warn-skip"``: log warning and skip malformed lines.
        """
        if parse_mode not in {"strict", "warn-skip"}:
            raise ValueError(f"Invalid parse_mode {parse_mode!r}; expected 'strict' or 'warn-skip'")

        for line_num, line in enumerate(line_iter, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                yield GTFRecord.from_str(line)
            except ValueError as e:
                if parse_mode == "strict":
                    raise ValueError(f"Malformed GTF at line {line_num}: {e}") from e
                logger.warning("Skipping malformed GTF line %d: %s", line_num, e)
                continue

    @classmethod
    def _parse_attrs(cls, raw_attr_field: str) -> tuple[dict[str, AttrValue], set[str]]:
        """Parse the attribute column into attrs and tags; a repeated key (such as GENCODE's ``ont``) becomes a list."""
        attrs: dict[str, AttrValue] = {}
        tags: set[str] = set()

        for match in cls.ATTR_TOKEN_PATTERN.finditer(raw_attr_field):
            key = match.group(1)
            value = match.group(2) if match.group(2) is not None else match.group(3).strip()

            if key == "tag":
                tags.add(value)
                continue

            if key in attrs:
                existing = attrs[key]
                if isinstance(existing, list):
                    existing.append(value)
                else:
                    attrs[key] = [existing, value]
            else:
                attrs[key] = value

        return attrs, tags

    @staticmethod
    def parse_file(
        filepath: str | Path,
        *,
        parse_mode: Literal["strict", "warn-skip"] = "strict",
    ) -> Generator["GTFRecord", None, None]:
        """Parse a GTF file, plain or gzipped (``.gz``)."""
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"GTF file not found: {filepath}")

        open_func = gzip.open if filepath.suffix == ".gz" else open
        with open_func(filepath, "rt") as f:
            yield from GTFRecord.parse(f, parse_mode=parse_mode)
