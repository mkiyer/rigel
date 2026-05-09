"""Count true-gDNA oracle boundary events by annotated BAM splice tag."""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.debug.oracle_boundary_density_split import (  # noqa: E402
    BASE_DEFAULT,
    SIDE_CLASSES,
    build_side_index,
    parse_gdna_qname,
)

CONDITION = "gdna_high_ss_0.99_nrna_none"


def _tag(read: pysam.AlignedSegment, name: str, default: str = "missing") -> str:
    try:
        return str(read.get_tag(name))
    except KeyError:
        return default


def main() -> None:
    base = BASE_DEFAULT
    region_df = pd.read_feather(base / "rigel_index" / "regions.feather")
    side_index = build_side_index(region_df)
    bam_path = base / CONDITION / "annotated.bam"
    q = 1
    ei_codes = {
        SIDE_CLASSES.index("exon_intron_internal"),
        SIDE_CLASSES.index("exon_intron_terminal"),
    }

    event_by_zs: Counter[str] = Counter()
    event_by_zc: Counter[str] = Counter()
    frag_by_zs: Counter[str] = Counter()
    missing_ei_events = 0
    total_ei_events = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            parsed = parse_gdna_qname(read.query_name)
            if parsed is None:
                continue
            ref_name, start, end = parsed
            positions = side_index.positions_by_ref.get(ref_name)
            if positions is None:
                continue
            class_codes = side_index.class_codes_by_ref[ref_name]
            left = int(np.searchsorted(positions, start + q, side="left"))
            right = int(np.searchsorted(positions, end - q, side="right"))
            if right <= left:
                continue
            n_ei = int(sum(1 for code in class_codes[left:right] if int(code) in ei_codes))
            if n_ei == 0:
                continue
            zs = _tag(read, "ZS", "missing")
            zc = _tag(read, "ZC", "missing")
            total_ei_events += n_ei
            event_by_zs[zs] += n_ei
            event_by_zc[zc] += n_ei
            frag_by_zs[zs] += 1
            if zs != "unspliced":
                missing_ei_events += n_ei

    print(f"Condition: {CONDITION}")
    print(f"Total true-gDNA physical exon-intron boundary events: {total_ei_events}")
    print(f"Events whose annotated splice tag is not unspliced: {missing_ei_events}")
    print("\nEvents by ZS splice tag:")
    for key, value in event_by_zs.most_common():
        print(f"  {key:<20s} {value:8d} ({value / total_ei_events:.3%})")
    print("\nFragments with EI events by ZS splice tag:")
    for key, value in frag_by_zs.most_common():
        print(f"  {key:<20s} {value:8d}")
    print("\nEvents by ZC ambiguity class:")
    for key, value in event_by_zc.most_common():
        print(f"  {key:<20s} {value:8d} ({value / total_ei_events:.3%})")


if __name__ == "__main__":
    main()
