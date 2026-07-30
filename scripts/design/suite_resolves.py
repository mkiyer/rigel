"""Can this benchmark suite resolve the axis it is about to be used to judge? Answer BEFORE running it.

    TODO item 2   ·   `CARRY_FORWARD.md` §3 traps 15 and 16

⛔ **WHY THIS EXISTS, AND WHY IT IS WRITTEN BEFORE THE SUITE.** The 32-condition `ambig_dense_10mb` suite
was used for months to judge a partition change it was **structurally incapable of seeing**: its fine node
set was row-for-row identical to its merged region set (1,698 == 1,698). It also had `frag_std = 0`, so
nothing fragment-length-dependent was exercised at all, and it was Poisson by construction. The standing
rule that came out of that: *before running a benchmark, prove it can resolve the axis you are changing.*
This script is that proof, made executable, and it is written first so it cannot be retrofitted to whatever
the suite turned out to be.

⭐ **THERE ARE NO TUNED THRESHOLDS HERE, DELIBERATELY** — a pass mark would be a magic number, and the
standing rule forbids one. Every requirement is scored against its **degenerate value**: the number a suite
scores when it is structurally blind to that requirement. Those are boundaries, not choices —

* a partition that cannot be resolved has `nodes / merged regions` **exactly 1.000` and **0** hidden termini;
* a suite with no fragment-length variation has variance **exactly 0**;
* a Poisson simulator has overdispersion **exactly 0** (measured `< 5e-5`);
* a reference with no single-stranded node to train the population prior on has **0** of them.

A requirement passes iff its measurement is strictly on the non-degenerate side. The **human index** is
printed beside it as the calibration reference — not as a threshold, but so a number that is technically
non-zero and practically useless is visible as such.

    python scripts/design/suite_resolves.py SUITE_INDEX [--compare-index HUMAN_INDEX] [--suite SUITE_DIR]

Without `--suite` only the structural requirements run and the rest say so. That is the useful mode while
the reference exists and the reads do not.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)
from rigel.calibration.splice_graph import (
    EDGE_KIND_CONTIGUOUS,
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
)

POS_BITS = BIT_EXON_POS | BIT_INTRON_POS
NEG_BITS = BIT_EXON_NEG | BIT_INTRON_NEG
TERMINUS_POS = FLAG_TSS_POS | FLAG_TES_POS
TERMINUS_NEG = FLAG_TSS_NEG | FLAG_TES_NEG


@dataclass(frozen=True)
class Verdict:
    """One requirement, its measurement, and the value a blind suite would score."""

    key: str  # the letter from `TODO.md` item 2
    requirement: str
    measured: float
    degenerate: float
    unit: str
    note: str = ""

    @property
    def resolves(self) -> bool:
        return self.measured > self.degenerate

    def render(self, reference: float | None) -> str:
        mark = "OK  " if self.resolves else "⛔ NO"
        ref = "" if reference is None else f"  human {reference:,.4g}"
        return (
            f"  {mark}  ({self.key}) {self.requirement:<44} "
            f"{self.measured:>12,.4g} {self.unit:<10} blind {self.degenerate:g}{ref}\n"
            f"        {self.note}"
        )


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# Structural requirements — properties of the REFERENCE, answerable before a single read exists
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def merged_region_boundaries(nodes: pd.DataFrame) -> np.ndarray:
    """Rebuild the v7 merged partition: a region ends where the reference or the signature changes."""
    signature = nodes["signature"].to_numpy(np.uint8)
    ref = nodes["ref_name"].astype(str).to_numpy()
    boundary = np.ones(len(nodes), dtype=bool)
    boundary[1:] = (signature[1:] != signature[:-1]) | (ref[1:] != ref[:-1])
    return boundary


def requirement_g_partition(nodes: pd.DataFrame, edges: pd.DataFrame) -> list[Verdict]:
    """(g) Can the suite see a partition change at all?

    ⛔ The deleted suite scored **exactly** the degenerate value on both of these: 1,698 nodes against
    1,698 merged regions, and therefore zero termini hidden by the merge.
    """
    boundary = merged_region_boundaries(nodes)
    n_merged = int(boundary.sum())
    ratio = len(nodes) / max(n_merged, 1)

    contiguous = edges["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
    dst = edges.loc[contiguous, "dst"].to_numpy(np.int64)
    flags = edges.loc[contiguous, "flags"].to_numpy(np.uint16)
    is_terminus = (flags & (TERMINUS_POS | TERMINUS_NEG)) != 0
    hidden = int((is_terminus & ~boundary[dst]).sum())

    return [
        Verdict(
            "g",
            "partition resolution (nodes / merged regions)",
            ratio,
            1.0,
            "x",
            f"{len(nodes):,} nodes against {n_merged:,} merged regions. "
            f"At 1.000x the two partitions are the same object and no partition change is visible.",
        ),
        Verdict(
            "g",
            "transcript termini the old merge hid",
            hidden,
            0.0,
            "seams",
            f"{100 * hidden / max(int(is_terminus.sum()), 1):.2f} % of terminus seams. "
            f"This is the visibility the v8 partition exists to buy.",
        ),
    ]


def requirement_d_interior_termini(nodes: pd.DataFrame, edges: pd.DataFrame) -> list[Verdict]:
    """(d) Alternative TSS/TES that fall strictly INSIDE an exon.

    A terminus seam whose flanking nodes are BOTH exonic on the terminus's own strand sits strictly
    inside an exon of some other isoform — which is the case a merged partition cannot represent and
    the case a real annotation is full of. A generated mini-genome with one isoform per gene has none.
    """
    signature = nodes["signature"].to_numpy(np.uint8)
    contiguous = edges["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
    dst = edges.loc[contiguous, "dst"].to_numpy(np.int64)
    flags = edges.loc[contiguous, "flags"].to_numpy(np.uint16)

    interior = np.zeros(dst.size, dtype=bool)
    for terminus_bits, exon_bit in ((TERMINUS_POS, BIT_EXON_POS), (TERMINUS_NEG, BIT_EXON_NEG)):
        exonic_both_sides = ((signature[dst - 1] & exon_bit) != 0) & ((signature[dst] & exon_bit) != 0)
        interior |= ((flags & terminus_bits) != 0) & exonic_both_sides

    total_terminus = max(int(((flags & (TERMINUS_POS | TERMINUS_NEG)) != 0).sum()), 1)
    return [
        Verdict(
            "d",
            "termini strictly inside an exon",
            int(interior.sum()),
            0.0,
            "seams",
            f"{100 * interior.sum() / total_terminus:.2f} % of terminus seams have EXON on both sides "
            f"on the terminus's own strand — an alternative TSS/TES within another isoform's exon.",
        )
    ]


def requirement_e_single_stranded(nodes: pd.DataFrame) -> list[Verdict]:
    """(e) Ample single-stranded nodes — the population prior trains on them.

    ⚠ Scored PER REFERENCE and by DISTANCE, not pooled. The stated failure is not "the suite has no
    single-stranded nodes anywhere", it is *an isolated both-strand node with no single-stranded
    neighbours* — a starved toy. So the measurement is: from each both-stranded node, how many nodes
    away is the nearest single-stranded one? A reference with both-stranded nodes and none at all
    scores an infinite distance, which is the degenerate case.
    """
    signature = nodes["signature"].to_numpy(np.uint8)
    ref = nodes["ref_name"].astype(str).to_numpy()
    has_pos = (signature & POS_BITS) != 0
    has_neg = (signature & NEG_BITS) != 0
    single = has_pos ^ has_neg
    both = has_pos & has_neg

    worst_distance = 0.0
    starved_refs = 0
    for name in np.unique(ref[both]) if both.any() else []:
        on_ref = ref == name
        single_idx = np.flatnonzero(on_ref & single)
        both_idx = np.flatnonzero(on_ref & both)
        if single_idx.size == 0:
            starved_refs += 1
            worst_distance = float("inf")
            continue
        # Distance in NODES to the nearest single-stranded node on the same reference.
        slot = np.searchsorted(single_idx, both_idx)
        sentinel = np.iinfo(np.int64).max
        left = np.where(slot > 0, both_idx - single_idx[np.maximum(slot - 1, 0)], sentinel)
        right = np.where(
            slot < single_idx.size,
            single_idx[np.minimum(slot, single_idx.size - 1)] - both_idx,
            sentinel,
        )
        worst_distance = max(worst_distance, float(np.minimum(left, right).max()))

    fraction = float(single.sum()) / max(len(nodes), 1)
    return [
        Verdict(
            "e",
            "single-stranded nodes",
            int(single.sum()),
            0.0,
            "nodes",
            f"{100 * fraction:.1f} % of nodes; {int(both.sum()):,} both-stranded, "
            f"{int((~has_pos & ~has_neg).sum()):,} intergenic.",
        ),
        Verdict(
            "e",
            "both-stranded nodes (a both-strand test needs some)",
            int(both.sum()),
            0.0,
            "nodes",
            "Without these there is no both-strand stress case for the single-stranded nodes to support.",
        ),
        Verdict(
            "e",
            "single-stranded nodes PER both-stranded node",
            float(single.sum()) / max(int(both.sum()), 1) if both.any() else 0.0,
            0.0,
            "ratio",
            f"This is what 'ample' means: {starved_refs} reference(s) carry a both-strand node with NO "
            f"single-stranded node at all, and the worst node-distance from a both-strand node to the "
            f"nearest single-stranded one is {worst_distance:g}.",
        ),
    ]


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# Empirical requirements — properties of the SIMULATED DATA
#
# ⭐ All four read the suite's own truth files, not a BAM scan: `truth_summary.json` carries the
# realised per-kind fragment-length moments, `truth_abundances.tsv` carries pre- and post-capture
# fragments per transcript, and `manifest.json` carries the condition grid. Nothing here needs rigel to
# run, which matters because calibration is red until S5.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def load_conditions(suite: Path) -> list[dict]:
    manifest = json.loads((suite / "manifest.json").read_text())
    return list(manifest.get("conditions", []))


def requirement_b_fragment_length_variance(suite: Path, conditions: list[dict]) -> list[Verdict]:
    """(b) Fragment-length variance. ⛔ The deleted suite had EXACTLY zero — every fragment 200 bp."""
    worst = None
    detail = ""
    for condition in conditions:
        summary_path = suite / condition["name"] / "truth_summary.json"
        if not summary_path.exists():
            continue
        lengths = json.loads(summary_path.read_text()).get("fragment_lengths", {})
        for kind in ("rna", "gdna"):
            stats = lengths.get(kind) or {}
            if not stats.get("n"):
                continue
            std = float(stats["std"])
            if worst is None or std < worst:
                worst = std
                detail = (
                    f"least-variable pool: {kind} in {condition['name']} — "
                    f"mean {stats['mean']:.1f}, sd {std:.1f}, range {stats['min']}-{stats['max']}"
                )
    return [
        Verdict(
            "b",
            "fragment-length sd, WORST pool over all conditions",
            0.0 if worst is None else worst,
            0.0,
            "bp",
            detail or "no truth_summary.json found",
        )
    ]


def requirement_a_density_step(suite: Path, conditions: list[dict]) -> list[Verdict]:
    """(a) A density STEP, not a uniform background.

    Over a run of flat nodes a relayed message decays geometrically per hop, so a uniform scenario
    cannot distinguish "the relay works" from "the global prior reached it". Capture supplies the step:
    on-panel transcripts are enriched, off-panel ones are not, and the ratio between them IS the step.
    Measured from the truth as post-capture fragments per unit pre-capture abundance.
    """
    best_ratio = 0.0
    detail = "no capture-on condition with a truth table"
    for condition in conditions:
        if not condition.get("capture_enabled"):
            continue
        table = suite / condition["name"] / "truth_abundances.tsv"
        if not table.exists():
            continue
        frame = pd.read_csv(table, sep="\t")
        expressed = frame[frame["pre_capture_total_rna"] > 0]
        if expressed.empty:
            continue
        enrichment = (
            expressed["post_capture_total_rna_fragments"] / expressed["pre_capture_total_rna"]
        )
        enriched = enrichment[enrichment > 0]
        if enriched.empty:
            continue
        high, low = enriched.quantile(0.9), enriched.quantile(0.1)
        ratio = float(high / low) if low > 0 else float("inf")
        if ratio > best_ratio:
            best_ratio = ratio
            detail = (
                f"{condition['name']}: p90/p10 of post-capture fragments per unit abundance = "
                f"{ratio:,.1f}x over {len(enriched):,} expressed transcripts. "
                f"At 1.0x the field is uniform and the relay cannot be told from the prior."
            )
    return [Verdict("a", "capture density step (p90/p10 enrichment)", best_ratio, 1.0, "x", detail)]


def requirement_f_low_gdna_corner(conditions: list[dict]) -> list[Verdict]:
    """(f) A low-gDNA x strong-capture corner. Real libraries live at 1-10 % gDNA."""
    corner = [
        c
        for c in conditions
        if c.get("capture_enabled") and 0.0 < float(c.get("gdna_rate", 0.0)) <= 0.10
    ]
    rates = sorted({float(c.get("gdna_rate", 0.0)) for c in conditions})
    return [
        Verdict(
            "f",
            "conditions at 0 < gDNA <= 10 % WITH capture on",
            len(corner),
            0.0,
            "conds",
            f"gDNA rates present: {rates}. 0/100/300 % conditions cannot see the corner real "
            f"libraries actually live in.",
        )
    ]


def requirement_c_overdispersion(conditions: list[dict]) -> list[Verdict]:
    """(c) Non-Poisson counts.

    ⚠ Reported as ESTIMABILITY, not as a number. With one draw per parameter set there is nothing to
    estimate overdispersion against: the simulator draws multinomial at fixed abundance, so counts are
    Poisson given truth and any `omega` computed from a single condition is measuring the estimator, not
    the data. It becomes measurable when the suite carries **replicate conditions** — identical
    parameters, different `sim_seed` — because then per-object variance across replicates is directly
    comparable to the mean.
    """
    signature: dict[tuple, int] = {}
    for condition in conditions:
        key = (
            condition.get("gdna_label"),
            condition.get("strand_specificity"),
            condition.get("nrna_label"),
            condition.get("capture_label"),
        )
        signature[key] = signature.get(key, 0) + 1
    replicated = sum(1 for count in signature.values() if count > 1)
    return [
        Verdict(
            "c",
            "replicate condition pairs (needed to estimate omega)",
            replicated,
            0.0,
            "pairs",
            "The simulator is Poisson by construction (`wgs_engine._accumulate_pool` draws multinomial "
            "at fixed abundance; measured omega < 5e-5), so overdispersion must be built IN and then "
            "measured across replicates. Deferred by owner ruling 2026-07-30.",
        )
    ]


def empirical_verdicts(suite: Path) -> list[Verdict]:
    conditions = load_conditions(suite)
    return [
        *requirement_a_density_step(suite, conditions),
        *requirement_b_fragment_length_variance(suite, conditions),
        *requirement_c_overdispersion(conditions),
        *requirement_f_low_gdna_corner(conditions),
    ]


def structural_verdicts(index_dir: Path) -> list[Verdict]:
    nodes = pd.read_feather(index_dir / "nodes.feather")
    edges = pd.read_feather(index_dir / "edges.feather")
    return [
        *requirement_g_partition(nodes, edges),
        *requirement_d_interior_termini(nodes, edges),
        *requirement_e_single_stranded(nodes),
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("suite_index", type=Path)
    ap.add_argument("--compare-index", type=Path, default=None, help="e.g. the human index, for scale")
    ap.add_argument("--suite", type=Path, default=None, help="Suite dir; enables the empirical checks")
    args = ap.parse_args()

    verdicts = structural_verdicts(args.suite_index)
    reference = None
    if args.compare_index is not None:
        reference = {v.requirement: v.measured for v in structural_verdicts(args.compare_index)}

    print(f"suite index  {args.suite_index}")
    if args.compare_index:
        print(f"compared to  {args.compare_index}")
    print("\nSTRUCTURAL REQUIREMENTS — properties of the reference, answerable with no reads at all")
    for verdict in verdicts:
        print(verdict.render(None if reference is None else reference.get(verdict.requirement)))

    print("\nEMPIRICAL REQUIREMENTS — properties of the simulated data, read from the suite's own truth")
    empirical: list[Verdict] = []
    if args.suite is None:
        print("  SKIP  no --suite given. These need simulated reads and are not inferable from the index.")
    else:
        empirical = empirical_verdicts(args.suite)
        for verdict in empirical:
            print(verdict.render(None))

    failed = [v for v in verdicts + empirical if not v.resolves]
    print()
    if failed:
        print("⛔ THIS SUITE CANNOT RESOLVE:")
        for verdict in failed:
            print(f"  ({verdict.key}) {verdict.requirement}")
        raise SystemExit(1)
    print("✅ every structural requirement resolves — no number from it is disqualified on these grounds")


if __name__ == "__main__":
    main()
