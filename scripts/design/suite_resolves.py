"""Can this benchmark suite resolve the axis it is about to be used to judge? Answer BEFORE running it.

    TODO item 2 · and 16

⛔ **WHY THIS EXISTS, AND WHY IT IS WRITTEN BEFORE THE SUITE.** The 32-condition `ambig_dense_10mb` suite
was used for months to judge a partition change it was **structurally incapable of seeing**: its fine region
set was row-for-row identical to its merged region set (1,698 == 1,698). It also had `frag_std = 0`, so
nothing fragment-length-dependent was exercised at all, and it was Poisson by construction. The standing
rule that came out of that: *before running a benchmark, prove it can resolve the axis you are changing.*
This script is that proof, made executable, and it is written first so it cannot be retrofitted to whatever
the suite turned out to be.

⭐ **THERE ARE NO TUNED THRESHOLDS HERE, DELIBERATELY** — a pass mark would be a magic number, and the
standing rule forbids one. Every requirement is scored against its **degenerate value**: the number a suite
scores when it is structurally blind to that requirement. Those are boundaries, not choices —

* a partition that cannot be resolved has `regions / merged regions` **exactly 1.000` and **0** hidden termini;
* a suite with no fragment-length variation has variance **exactly 0**;
* a Poisson simulator has overdispersion **exactly 0** (measured `< 5e-5`);
* a reference with no single-stranded region to train the population prior on has **0** of them.

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

    key: str # the letter from
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


def merged_region_boundaries(regions: pd.DataFrame) -> np.ndarray:
    """Rebuild the v7 merged partition: a region ends where the reference or the signature changes."""
    signature = regions["signature"].to_numpy(np.uint8)
    ref = regions["ref_name"].astype(str).to_numpy()
    boundary = np.ones(len(regions), dtype=bool)
    boundary[1:] = (signature[1:] != signature[:-1]) | (ref[1:] != ref[:-1])
    return boundary


def requirement_g_partition(regions: pd.DataFrame, boundaries: pd.DataFrame) -> list[Verdict]:
    """(g) Can the suite see a partition change at all?

    ⛔ The deleted suite scored **exactly** the degenerate value on both of these: 1,698 regions against
    1,698 merged regions, and therefore zero termini hidden by the merge.
    """
    boundary = merged_region_boundaries(regions)
    n_merged = int(boundary.sum())
    ratio = len(regions) / max(n_merged, 1)

    contiguous = boundaries["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
    dst = boundaries.loc[contiguous, "dst"].to_numpy(np.int64)
    flags = boundaries.loc[contiguous, "flags"].to_numpy(np.uint16)
    is_terminus = (flags & (TERMINUS_POS | TERMINUS_NEG)) != 0
    hidden = int((is_terminus & ~boundary[dst]).sum())

    return [
        Verdict(
            "g",
            "partition resolution (regions / merged regions)",
            ratio,
            1.0,
            "x",
            f"{len(regions):,} regions against {n_merged:,} merged regions. "
            f"At 1.000x the two partitions are the same object and no partition change is visible.",
        ),
        Verdict(
            "g",
            "transcript termini the old merge hid",
            hidden,
            0.0,
            "boundaries",
            f"{100 * hidden / max(int(is_terminus.sum()), 1):.2f} % of terminus boundaries. "
            f"This is the visibility the v8 partition exists to buy.",
        ),
    ]


def requirement_d_interior_termini(regions: pd.DataFrame, boundaries: pd.DataFrame) -> list[Verdict]:
    """(d) Alternative TSS/TES that fall strictly INSIDE an exon.

    A terminus boundary whose flanking regions are BOTH exonic on the terminus's own strand sits strictly
    inside an exon of some other isoform — which is the case a merged partition cannot represent and
    the case a real annotation is full of. A generated mini-genome with one isoform per gene has none.
    """
    signature = regions["signature"].to_numpy(np.uint8)
    contiguous = boundaries["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
    dst = boundaries.loc[contiguous, "dst"].to_numpy(np.int64)
    flags = boundaries.loc[contiguous, "flags"].to_numpy(np.uint16)

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
            "boundaries",
            f"{100 * interior.sum() / total_terminus:.2f} % of terminus boundaries have EXON on both sides "
            f"on the terminus's own strand — an alternative TSS/TES within another isoform's exon.",
        )
    ]


def requirement_e_single_stranded(regions: pd.DataFrame) -> list[Verdict]:
    """(e) Ample single-stranded regions — the population prior trains on them.

    ⚠ Scored PER REFERENCE and by DISTANCE, not pooled. The stated failure is not "the suite has no
    single-stranded regions anywhere", it is *an isolated both-strand region with no single-stranded
    neighbours* — a starved toy. So the measurement is: from each both-stranded region, how many regions
    away is the nearest single-stranded one? A reference with both-stranded regions and none at all
    scores an infinite distance, which is the degenerate case.
    """
    signature = regions["signature"].to_numpy(np.uint8)
    ref = regions["ref_name"].astype(str).to_numpy()
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
        # Distance in REGIONS to the nearest single-stranded region on the same reference.
        slot = np.searchsorted(single_idx, both_idx)
        sentinel = np.iinfo(np.int64).max
        left = np.where(slot > 0, both_idx - single_idx[np.maximum(slot - 1, 0)], sentinel)
        right = np.where(
            slot < single_idx.size,
            single_idx[np.minimum(slot, single_idx.size - 1)] - both_idx,
            sentinel,
        )
        worst_distance = max(worst_distance, float(np.minimum(left, right).max()))

    fraction = float(single.sum()) / max(len(regions), 1)
    return [
        Verdict(
            "e",
            "single-stranded regions",
            int(single.sum()),
            0.0,
            "regions",
            f"{100 * fraction:.1f} % of regions; {int(both.sum()):,} both-stranded, "
            f"{int((~has_pos & ~has_neg).sum()):,} intergenic.",
        ),
        Verdict(
            "e",
            "both-stranded regions (a both-strand test needs some)",
            int(both.sum()),
            0.0,
            "regions",
            "Without these there is no both-strand stress case for the single-stranded regions to support.",
        ),
        Verdict(
            "e",
            "single-stranded regions PER both-stranded region",
            float(single.sum()) / max(int(both.sum()), 1) if both.any() else 0.0,
            0.0,
            "ratio",
            f"This is what 'ample' means: {starved_refs} reference(s) carry a both-strand region with NO "
            f"single-stranded region at all, and the worst region-distance from a both-strand region to the "
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

    Over a run of flat regions a relayed message decays geometrically per hop, so a uniform scenario
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


def requirement_h_length_gap_regime(suite: Path, conditions: list[dict]) -> list[Verdict]:
    """(h) ⭐ **G-S7 — does the suite exercise a NARROWED gDNA<->RNA length gap?**

    ``mu_g - mu_r`` is the **only** thing that identifies the fragment-length channel: at equal
    component means ``(count, sum 1/L)`` carries exactly zero information about composition at any
    depth. So a suite that presents one single value of that gap cannot resolve anything the length
    channel does, and a suite whose capture arm leaves the gap untouched has never tested the tool
    against the regime real capture produces.

    ⛔ **The panel this replaces scored the degenerate value exactly.** Its post-capture truth summaries
    were byte-identical to the pre-capture ones — gDNA mean 195.57 and RNA 217.24 on BOTH capture arms —
    so the gap was −21.66 bp twice and the narrowing was 0.00. Hybrid capture hybridises probes to
    sequence, so it selects for length, and it narrows the gap whenever gDNA is the shorter component
    because the short tail it removes is disproportionately gDNA.

    ⚠ **Directional, no threshold.** A magnitude would be inventing the capture efficiency curve. The
    degenerate value is a *measured* 0.00, not a chosen one.
    """
    def gap(name: str) -> float | None:
        path = suite / name / "truth_summary.json"
        if not path.exists():
            return None
        lengths = json.loads(path.read_text()).get("fragment_lengths", {})
        gdna, rna = lengths.get("gdna") or {}, lengths.get("rna") or {}
        if not gdna.get("n") or not rna.get("n"):
            return None
        return abs(float(gdna["mean"]) - float(rna["mean"]))

    arms: dict[tuple, dict[str, str]] = {}
    for condition in conditions:
        key = (
            condition.get("gdna_label"),
            condition.get("strand_specificity"),
            condition.get("nrna_label"),
            condition.get("gdna_strand_overdispersion"),
        )
        arms.setdefault(key, {})[str(condition.get("capture_label"))] = str(condition["name"])

    narrowings: list[tuple[str, float, float]] = []
    for key, by_label in sorted(arms.items(), key=lambda kv: str(kv[0])):
        if "on" not in by_label or "off" not in by_label:
            continue
        gap_off, gap_on = gap(by_label["off"]), gap(by_label["on"])
        if gap_off is None or gap_on is None:
            continue
        narrowings.append((f"{key[0]} ss{key[1]}", gap_off, gap_on))

    observed = sorted({round(value, 2) for _n, off, on in narrowings for value in (off, on)})
    worst = min((off - on for _n, off, on in narrowings), default=0.0)
    detail = (
        "; ".join(f"{name}: {off:.2f} -> {on:.2f} bp" for name, off, on in narrowings)
        or "no gDNA-bearing condition with both capture arms"
    )
    return [
        Verdict(
            "h",
            "capture narrows |mu_g - mu_r|, WORST condition",
            worst,
            0.0,
            "bp",
            f"{detail}. At exactly 0.00 the capture arm is length-neutral and the length channel has "
            f"never been tested against the regime real capture produces.",
        ),
        Verdict(
            "h",
            "distinct |mu_g - mu_r| regimes presented",
            len(observed),
            1.0,
            "values",
            f"observed gaps (bp): {observed}. One value is one regime, and the channel the gap "
            f"identifies cannot be resolved against a constant.",
        ),
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
        *requirement_h_length_gap_regime(suite, conditions),
    ]


def structural_verdicts(index_dir: Path) -> list[Verdict]:
    regions = pd.read_feather(index_dir / "regions.feather")
    boundaries = pd.read_feather(index_dir / "edges.feather")
    return [
        *requirement_g_partition(regions, boundaries),
        *requirement_d_interior_termini(regions, boundaries),
        *requirement_e_single_stranded(regions),
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
