"""How much library mass does the arbitration RESOLVE, and how much goes to the side buffer?

    Rule: ``docs/SPEC_GAP_PATHS.md`` §0–§2   ·   Plan: ``docs/PLAN_TWO_PASS.md``   ·   Ledger: S1

⭐ WHY THIS NUMBER IS NEEDED BEFORE THE SECOND PASS IS BUILT. A fragment arrives at the accumulator with
its hypothesis SET — every explanation of what its unsequenced mate gaps contain, the empty set being the
genomic one. Exactly one survivor deposits; two or more and ``L`` is undetermined, so the fragment is held
WHOLE for the second pass. ⚠ A single retained-intron isoform among the candidates is enough to hold a
fragment, and retained-intron isoforms are common in GENCODE. If that holds most of the population the
side buffer stops being a tidy-up and becomes load-bearing, and the second pass has to be built for scale.

⚠ **The held population is the LONG one** — a longer gap admits more hypotheses — so between S1 and S3 the
tally is deliberately thinner and biased short. ⛔ Do not read any accuracy number off this state; the gate
here is conservation, not accuracy (``PLAN_TWO_PASS.md`` §4).

    OMP_NUM_THREADS=1 python scripts/design/implicit_splice_census.py INDEX BAM [BAM ...]

"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer


#: The denominators worth reading side by side, in the order that tells the story. ⚠ `sj_implicit_fragments`
#: used to lead this list and is GONE: D1 is deleted, so "this L was partly inferred" no longer selects
#: anything — a fragment deposits when exactly ONE hypothesis survives, however it got there.
_KEYS = [
    "deposited",
    "deferred_undetermined_gap",
    "dropped_too_long",
    "dropped_strand_undefined",
    "dropped_empty",
    "unannotated_introns",
    "contradictory_sj_strand",
    "introns_absorbed",
]

#: ⭐ The umbrella census, which is the interesting half: the deferred total split by WHICH QUESTION is
#: open. `rna_or_gdna` is one bit and it is the composition question calibration exists to answer;
#: `which_introns` is certified RNA with only the structure open; `both` is both at once.
_CENSUS_KEYS = [
    "gap_resolved_spliced",
    "gap_deferred_rna_or_gdna",
    "gap_deferred_which_introns",
    "gap_deferred_both",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index")
    ap.add_argument("bams", nargs="+")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    print(f"index  {args.index}   {len(index.nodes_df):,} nodes\n")

    header = f"{'library':<22}" + "".join(f"{key:>26}" for key in _KEYS + _CENSUS_KEYS)
    print(header)
    print("-" * len(header))

    for bam in args.bams:
        _, _, _, payload = scan_and_buffer(
            bam, index, BamScanConfig(sj_strand_tag="auto", total_threads=args.threads)
        )
        qc = {f.name: getattr(payload.qc, f.name) for f in dataclasses.fields(payload.qc)}
        census = {
            f.name: getattr(payload.gap_resolution, f.name)
            for f in dataclasses.fields(payload.gap_resolution)
        }
        row = f"{Path(bam).parts[-3][:22]:<22}"
        for key in _KEYS:
            row += f"{qc.get(key, 0):>26,}"
        for key in _CENSUS_KEYS:
            row += f"{census.get(key, 0):>26,}"
        print(row)

        umbrella = sum(census.values())
        held = payload.gap_resolution.deferred
        print(
            f"  gaps needing resolution: {umbrella:,} — "
            f"{census['gap_resolved_spliced']:,} resolved "
            f"({100 * census['gap_resolved_spliced'] / umbrella if umbrella else 0:.1f} %), "
            f"{held:,} held for the second pass "
            f"({100 * held / umbrella if umbrella else 0:.1f} %)"
        )
        print(
            f"  held as a share of everything offered: "
            f"{100 * held / max(qc['deposited'] + held, 1):.2f} %"
        )
        # ⭐ Conservation, on the two identities the payload already refuses at its door — printed because a
        # census is only readable if the reader can see that nothing fell out of it.
        print(
            f"  bank holds {payload.deferred.n_fragments:,} fragments "
            f"({payload.deferred.n_hypotheses:,} hypotheses); counter says "
            f"{qc['deferred_undetermined_gap']:,}; census subclasses sum to {held:,}\n"
        )


if __name__ == "__main__":
    main()
