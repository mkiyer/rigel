"""Design a hybrid-capture probe panel over a suite reference. ⭐ This is also the DENSITY STEP.

    TODO item 2, requirement (a)

⭐ **The panel is not scaffolding — it is one of the seven requirements.** `TODO.md` item 2 asks for *a
density step, not just a uniform background*, because over a run of flat nodes a relayed message decays
geometrically per hop, so a uniform scenario cannot distinguish "the relay works" from "the global prior
reached it". A capture panel that covers some gene groups and not others IS that step, and a sharp one:
`CARRY_FORWARD.md` §1 fact 15 measures hybrid capture at **~1000x** on exons and 1.0 off target, so the
captured/uncaptured boundary is a 3-decade cliff in gDNA density with real transcripts on both sides.

⚠ The whole-genome (existing-reference) simulator **requires** a panel when capture is enabled — unlike
the generated-mini-genome path, it does not design one (`whole_genome.py:170`,
`require_probes_when_enabled=True`). So this is a prerequisite for the capture axis, not an extra.

`capture_fraction` is a scenario parameter, declared in the config exactly like `gdna_rate` — it is
*which half of the annotation is on-panel*, not a tuned constant.

    python scripts/sim/design_suite_probes.py --gtf genes.gtf -o panel.tsv --capture-fraction 0.5
"""

from __future__ import annotations

import argparse
from pathlib import Path

from rigel.sim.capture.design import write_random_capture_probes
from rigel.sim.whole_genome import load_transcripts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--capture-fraction", type=float, required=True, help="share of GENE GROUPS on panel")
    ap.add_argument("--probe-length", type=int, default=120)
    ap.add_argument("--probe-density", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    transcripts = load_transcripts(args.gtf, transcript_filter="all")
    result = write_random_capture_probes(
        transcripts,
        args.out,
        capture_fraction=args.capture_fraction,
        probe_length=args.probe_length,
        probe_density=args.probe_density,
        seed=args.seed,
    )
    print(f"transcripts       {len(transcripts):,}")
    print(f"eligible          {result.n_eligible:,} over {result.n_eligible_genes:,} gene groups")
    print(f"gene groups ON    {result.n_captured_genes:,}  "
          f"({100 * result.n_captured_genes / max(result.n_eligible_genes, 1):.1f} %)")
    print(f"probes            {result.n_probes:,}")
    print(f"written           {result.path}")


if __name__ == "__main__":
    main()
