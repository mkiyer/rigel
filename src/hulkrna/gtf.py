import gzip
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set, Union

logger = logging.getLogger(__name__)


AttrValue = Union[str, List[str]]


@dataclass(slots=True)
class GTF:
    """
    Represents a single line in a GTF file
    """
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: Optional[float]
    strand: str
    phase: str
    attrs: Dict[str, AttrValue] = field(default_factory=dict)
    tags: Set[str] = field(default_factory=set)

    # Pre-compiled regex for performance
    # Matches: key "value";  or  key value;
    # Group 1: key, Group 2: quoted value, Group 3: unquoted value
    ATTR_TOKEN_PATTERN = re.compile(r'\s*([^";\s]+)\s+(?:"([^"]*)"|([^;]+))\s*;')

    @classmethod
    def from_str(cls, s: str) -> 'GTF':
        """
        Parses a raw GTF string line into a GTF object.
        """
        fields = s.rstrip('\n').split('\t')
        if len(fields) < 9:
            raise ValueError(f"Malformed GTF line, expected 9 columns: {s}")

        # Parse coordinates (Convert 1-based to 0-based)
        try:
            # GTF is 1-based inclusive. 
            # We convert to 0-based exclusive (Python/BED standard)
            start_pos = int(fields[3]) - 1 
            end_pos = int(fields[4])
        except ValueError:
            raise ValueError(f"Non-integer coordinates in line: {fields}")

        # Parse Score: '.' -> None
        score_val = None
        if fields[5] != '.':
            try:
                score_val = float(fields[5])
            except ValueError as exc:
                raise ValueError(f"Non-numeric score in line: {fields[5]}") from exc

        # Validate phase field (CDS: 0/1/2 or '.')
        phase_val = fields[7]
        if phase_val != '.' and phase_val not in {'0', '1', '2'}:
            raise ValueError(f"Invalid phase value: {phase_val}")

        # Parse Attributes
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
            tags=parsed_tags
        )
    
    def __str__(self) -> str:
        """
        Convert GTF object back to string representation.
        """
        # Convert score and phase back to string
        score_str = '.' if self.score is None else str(self.score)

        # Reconstruct attributes field
        attr_parts = []
        for k, v in self.attrs.items():
            if isinstance(v, list):
                for vv in v:
                    attr_parts.append(f'{k} "{vv}";')
            else:
                attr_parts.append(f'{k} "{v}";')
        for k in self.tags:
            attr_parts.append(f'tag "{k}";')
        attr_field = ' '.join(attr_parts)

        fields = [
            self.seqname,
            self.source,
            self.feature,
            str(self.start + 1),  # Convert back to 1-based
            str(self.end),
            score_str,
            self.strand,
            self.phase,
            attr_field
        ]
        return '\t'.join(fields)

    @staticmethod
    def parse(line_iter: Iterable[str]) -> Generator['GTF', None, None]:
        """
        Generator to parse lines from an iterator (file object or list).
        """
        for line_num, line in enumerate(line_iter, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                yield GTF.from_str(line)
            except ValueError as e:
                logger.warning("Skipping malformed GTF line %d: %s", line_num, e)
                continue


    @classmethod
    def _parse_attrs(cls, raw_attr_field: str) -> tuple[Dict[str, AttrValue], Set[str]]:
        """
        Parse the attribute column into attrs and tags.
        Preserves duplicate keys (e.g., multiple ont values in GENCODE) as lists.
        """
        attrs: Dict[str, AttrValue] = {}
        tags: Set[str] = set()

        for match in cls.ATTR_TOKEN_PATTERN.finditer(raw_attr_field):
            key = match.group(1)
            value = match.group(2) if match.group(2) is not None else match.group(3).strip()

            if key == 'tag':
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
    def parse_file(filepath: Union[str, Path]) -> Generator['GTF', None, None]:
        """
        Parse a GTF file and yield GTF objects.
        Handles both plain text and gzip files (.gz).
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"GTF file not found: {filepath}")
        
        open_func = gzip.open if filepath.suffix == '.gz' else open    
        with open_func(filepath, 'rt') as f:
            yield from GTF.parse(f)
