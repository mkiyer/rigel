"""How many implicit-splice fragments does the S3 rule DEPOSIT, and how many does it defer?

    Rule: ``docs/ACCUMULATOR_DESIGN.md`` §9.1   ·   Ledger: ``docs/LEDGER.md`` S3

⭐ WHY THIS NUMBER IS NEEDED BEFORE THE SIDE BUFFER IS DESIGNED. A `SPLICE_IMPLICIT` fragment deposits iff
its candidate transcripts all imply the SAME intron set; otherwise `L` is undetermined and it is deferred
to the second pass as `path_ambiguous`. ⚠ The **empty** set counts as a hypothesis, so a single
retained-intron isoform among the candidates is enough to defer a fragment — and retained-intron isoforms
are common in GENCODE. If that defers most of the implicit population, the side buffer stops being a
tidy-up and becomes load-bearing. Nobody has measured it, so this measures it.

    OMP_NUM_THREADS=1 python scripts/design/implicit_splice_census.py INDEX BAM [BAM ...]

"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer


#: The denominators worth reading side by side, in the order that tells the story.
_KEYS = [
    "deposited",
    "sj_implicit_fragments",
    "dropped_ambiguous_path",
    "dropped_too_long",
    "dropped_strand_undefined",
    "dropped_empty",
    "unannotated_introns",
    "contradictory_sj_strand",
    "introns_absorbed",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index")
    ap.add_argument("bams", nargs="+")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    print(f"index  {args.index}   {len(index.nodes_df):,} nodes\n")

    header = f"{'library':<22}" + "".join(f"{key:>26}" for key in _KEYS)
    print(header)
    print("-" * len(header))

    for bam in args.bams:
        _, _, _, _, payload = scan_and_buffer(
            bam, index, BamScanConfig(sj_strand_tag="auto", total_threads=args.threads)
        )
        qc = {f.name: getattr(payload.qc, f.name) for f in dataclasses.fields(payload.qc)}
        row = f"{Path(bam).parts[-3][:22]:<22}"
        for key in _KEYS:
            row += f"{qc.get(key, 0):>26,}"
        print(row)

        implicit = qc["sj_implicit_fragments"]
        ambiguous = qc["dropped_ambiguous_path"]
        total = implicit + ambiguous
        print(
            f"  implicit splices seen: {total:,} — "
            f"{implicit:,} deposited ({100 * implicit / total if total else 0:.1f} %), "
            f"{ambiguous:,} deferred to the second pass "
            f"({100 * ambiguous / total if total else 0:.1f} %)"
        )
        print(
            f"  deferred as a share of everything accepted: "
            f"{100 * ambiguous / max(qc['deposited'] + ambiguous, 1):.2f} %\n"
        )


if __name__ == "__main__":
    main()
