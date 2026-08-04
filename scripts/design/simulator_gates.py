r"""The simulator's acceptance gates — G-S1..G-S6, scored on the panel's own per-fragment truth.

    `docs/TESTING.md` §3

⭐ **Every gate here is DIRECTIONAL or an absolute count.** Not one of them carries a threshold: a pass
mark on "how much longer is a captured fragment" would be inventing the capture efficiency curve, which
is `docs/TRAPS.md` G1. The *sign* is physics — hybrid capture hybridises probes to
sequence, so a short fragment presents less sequence, binds worse, and is captured less efficiently —
and the sign is all that is asserted. The magnitude is whatever `binding_per_base` and the probe length
imply.

| | gate | form |
|---|---|---|
| G-S1 | gDNA fragments on an RNA-only reference | absolute count, must be 0 |
| G-S2 | genomic references carrying gDNA | count, must be >= 2, each non-zero |
| G-S3 | gDNA mean length, capture on vs off | strictly greater under capture |
| G-S4 | gDNA mean length, on-target vs off-target | on-target strictly longer |
| G-S5 | \|mu_gdna - mu_rna\|, capture on vs off | strictly narrower under capture |
| G-S6 | gDNA fragments longer than their own reference | absolute count, must be 0 |

⚠ **G-S4 is scored on PROBE OVERLAP, not on the start's territory, and the difference is not cosmetic**
(`docs/TRAPS.md` F6). Conditioned on being captured, a fragment whose start is in the intron is one that
was **long enough to reach the probe**, so an intronic start selects long fragments by construction
(weight ~ w^2/2) while an exonic start does not (weight ~ p^2/2, flat in w). So an exonic-start population
reads SHORTER than an intronic-start one under any capture model of this form, and gating on that pair is
gating on geometry. The population that physically binds is the one that **overlaps a probe**, so that is
what is scored; the start-territory table is printed underneath as the diagnostic it is.

⚠ **G-S4 is a REGRESSION GUARD, not a falsification** (`docs/TRAPS.md` A3): it passed with the capture
defect present, because the engine's per-fragment conditional was right and only the marginal was being
discarded. G-S3 is the gate that falsifies.

⚠ **G-S3 and G-S5 read `truth_fragment_lengths.tsv`**, the panel's own post-capture empirical truth,
and never the configured `frag_mean`. The configured value describes a library that was never
sequenced; post-capture is the baseline (`docs/TRAPS.md` F5).

    python scripts/design/simulator_gates.py --suite ~/Downloads/rigel_runs/suite/pilot \\
        --reference ~/Downloads/rigel_runs/suite/reference [--genomic-refs chr21 chr22]

Exits non-zero if any gate fails.
"""

from __future__ import annotations

import argparse
import array
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


# ── interval machinery ──────────────────────────────────────────────────────────────────────────


def merge(intervals: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    """Merge and return ``(starts, ends)`` as sorted, disjoint, half-open arrays."""
    if not intervals:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    ordered = sorted(intervals)
    out: list[list[int]] = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    arr = np.asarray(out, dtype=np.int64)
    return arr[:, 0], arr[:, 1]


def overlaps_any(
    frag_start: np.ndarray,
    frag_end: np.ndarray,
    iv_start: np.ndarray,
    iv_end: np.ndarray,
) -> np.ndarray:
    """Boolean: does ``[frag_start, frag_end)`` intersect any merged interval?

    The intervals are disjoint and sorted, so the only candidate is the last one starting at or
    before ``frag_end - 1`` — it is the only interval that can reach back into the fragment.
    """
    if iv_start.size == 0:
        return np.zeros(len(frag_start), dtype=bool)
    idx = np.searchsorted(iv_start, frag_end, side="left") - 1
    hit = idx >= 0
    result = np.zeros(len(frag_start), dtype=bool)
    result[hit] = iv_end[idx[hit]] > frag_start[hit]
    return result


def contains(point: np.ndarray, iv_start: np.ndarray, iv_end: np.ndarray) -> np.ndarray:
    """Boolean: is ``point`` inside any merged interval?"""
    if iv_start.size == 0:
        return np.zeros(len(point), dtype=bool)
    idx = np.searchsorted(iv_start, point, side="right") - 1
    hit = idx >= 0
    result = np.zeros(len(point), dtype=bool)
    result[hit] = iv_end[idx[hit]] > point[hit]
    return result


# ── inputs ──────────────────────────────────────────────────────────────────────────────────────


def read_fai(path: Path) -> dict[str, int]:
    return {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in path.read_text().splitlines()
        if line
    }


def read_annotation(gtf: Path) -> tuple[dict[str, list], dict[str, list]]:
    """Return per-reference ``(exon intervals, transcript-span intervals)`` from a GTF."""
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 5)
            if len(fields) < 5:
                continue
            ref, feature = fields[0], fields[2]
            if feature == "exon":
                exons[ref].append((int(fields[3]) - 1, int(fields[4])))
            elif feature == "transcript":
                spans[ref].append((int(fields[3]) - 1, int(fields[4])))
    return exons, spans


def read_probe_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Genomic probe blocks per reference, from the panel's BED12."""
    probes: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track")):
            continue
        fields = line.split("\t")
        ref, chrom_start = fields[0], int(fields[1])
        sizes = [int(v) for v in fields[10].rstrip(",").split(",") if v]
        offsets = [int(v) for v in fields[11].rstrip(",").split(",") if v]
        for size, offset in zip(sizes, offsets):
            if size > 0:
                probes[ref].append((chrom_start + offset, chrom_start + offset + size))
    return probes


def read_truth_length_stats(path: Path) -> dict[str, dict[str, float]]:
    """``kind -> {n, mean, sd}`` from a condition's ``truth_fragment_lengths.tsv``."""
    totals: dict[str, tuple[float, float, float]] = {}
    acc: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    with open(path) as handle:
        next(handle)
        for line in handle:
            kind, length_text, count_text, _fraction = line.rstrip("\n").split("\t")
            length, count = float(length_text), float(count_text)
            entry = acc[kind]
            entry[0] += count
            entry[1] += count * length
            entry[2] += count * length * length
    for kind, (n, s1, s2) in acc.items():
        totals[kind] = (n, s1, s2)
    stats: dict[str, dict[str, float]] = {}
    for kind, (n, s1, s2) in totals.items():
        mean = s1 / n if n else float("nan")
        var = s2 / n - mean * mean if n else float("nan")
        stats[kind] = {"n": n, "mean": mean, "sd": var**0.5 if var > 0 else 0.0}
    return stats


def read_gdna_origins(bam: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per-reference ``(starts, ends)`` of every gDNA fragment, from the oracle BAM's read names.

    ⭐ **The read name IS the per-fragment truth** — ``gdna:chr22:21469960-21470127:f:0`` names the
    reference and the molecule's own coordinates. One record per fragment is taken by filtering to
    R1 (flag 0x40), which is cheaper and exact where a de-duplicating set over ten million query
    names is neither.
    """
    starts: dict[str, array.array] = defaultdict(lambda: array.array("q"))
    ends: dict[str, array.array] = defaultdict(lambda: array.array("q"))
    view = subprocess.Popen(
        ["samtools", "view", "-f", "64", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    assert view.stdout is not None
    for line in view.stdout:
        if not line.startswith("gdna:"):
            continue
        qname = line.split("\t", 1)[0]
        _tag, ref, interval, _strand, _index = qname.split(":")
        start_text, end_text = interval.split("-")
        starts[ref].append(int(start_text))
        ends[ref].append(int(end_text))
    view.stdout.close()
    if view.wait() != 0:
        raise SystemExit(f"samtools view failed on {bam}")
    return {
        ref: (np.frombuffer(starts[ref], dtype=np.int64), np.frombuffer(ends[ref], dtype=np.int64))
        for ref in starts
    }


# ── per-condition census ────────────────────────────────────────────────────────────────────────


def gdna_census(
    origins: dict[str, tuple[np.ndarray, np.ndarray]],
    ref_lengths: dict[str, int],
    exon_iv: dict[str, tuple[np.ndarray, np.ndarray]],
    span_iv: dict[str, tuple[np.ndarray, np.ndarray]],
    probe_iv: dict[str, tuple[np.ndarray, np.ndarray]],
) -> dict:
    """Decompose a condition's gDNA fragments by reference, probe overlap and start territory."""
    empty = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))
    per_ref: dict[str, int] = {}
    sums: dict[str, list[float]] = defaultdict(lambda: [0.0, 0.0])  # [n, sum(length)]
    over_reference = 0

    for ref, (start, end) in origins.items():
        length = end - start
        per_ref[ref] = int(len(start))
        over_reference += int(np.count_nonzero(end > ref_lengths.get(ref, 0)))

        on_target = overlaps_any(start, end, *probe_iv.get(ref, empty))
        for label, mask in (("on_target", on_target), ("off_target", ~on_target)):
            sums[label][0] += int(mask.sum())
            sums[label][1] += float(length[mask].sum())

        in_exon = contains(start, *exon_iv.get(ref, empty))
        in_span = contains(start, *span_iv.get(ref, empty))
        territory = {
            "exonic": in_exon,
            "intronic": in_span & ~in_exon,
            "intergenic": ~in_span,
        }
        for label, mask in territory.items():
            sums[label][0] += int(mask.sum())
            sums[label][1] += float(length[mask].sum())

    return {
        "per_ref": per_ref,
        "n_total": int(sum(per_ref.values())),
        "over_reference": over_reference,
        "means": {
            label: {"n": int(n), "mean": (total / n if n else float("nan"))}
            for label, (n, total) in sums.items()
        },
    }


def base_key(condition: dict) -> tuple:
    """The condition identity with the capture arm removed — G-S3 and G-S5 pair on this."""
    return (
        condition["gdna_label"],
        condition["strand_specificity"],
        condition["nrna_label"],
        condition.get("gdna_strand_overdispersion", 0.0),
    )


# ── the gates ───────────────────────────────────────────────────────────────────────────────────


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[str, bool, str]] = []

    def gate(self, name: str, passed: bool, detail: str) -> None:
        self.rows.append((name, passed, detail))

    def emit(self) -> bool:
        print("\n══ GATES ══")
        width = max(len(name) for name, _, _ in self.rows)
        for name, passed, detail in self.rows:
            print(f"  {'PASS' if passed else 'FAIL'}  {name:<{width}}  {detail}")
        n_failed = sum(1 for _, passed, _ in self.rows if not passed)
        print(f"\n  {len(self.rows) - n_failed}/{len(self.rows)} gates pass")
        return n_failed == 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suite", type=Path, required=True, help="pilot directory (holds manifest.json)")
    ap.add_argument("--reference", type=Path, required=True, help="reference directory")
    ap.add_argument(
        "--genomic-refs",
        nargs="*",
        default=None,
        help="references that carry gDNA. Default: the manifest's own record of it.",
    )
    args = ap.parse_args()

    manifest = json.loads((args.suite / "manifest.json").read_text())
    conditions = manifest["conditions"]
    ref_lengths = read_fai(args.reference / "genome.fa.fai")

    genomic_refs = args.genomic_refs
    if genomic_refs is None:
        genomic_refs = (manifest.get("gdna") or {}).get("genomic_refs")
    if not genomic_refs:
        raise SystemExit(
            "no genomic reference set: pass --genomic-refs, or regenerate the panel with a "
            "simulator whose scenario declares `gdna.genomic_refs` (recorded in the manifest)"
        )
    genomic = set(genomic_refs)
    rna_only = sorted(set(ref_lengths) - genomic)
    print(f"genomic references  {len(genomic)}  {sorted(genomic)}")
    print(f"RNA-only references {len(rna_only)}  (first: {rna_only[:3]}{'...' if rna_only else ''})")

    exons, spans = read_annotation(args.reference / "genes.gtf")
    exon_iv = {ref: merge(iv) for ref, iv in exons.items()}
    span_iv = {ref: merge(iv) for ref, iv in spans.items()}
    probe_bed = args.reference / "capture_panel.bed"
    probe_iv = {ref: merge(iv) for ref, iv in read_probe_bed(probe_bed).items()} if probe_bed.exists() else {}
    probe_bases = sum(int((e - s).sum()) for s, e in probe_iv.values())
    print(f"probe bases         {probe_bases:,}  over {len(probe_iv)} references")

    census: dict[str, dict] = {}
    lengths: dict[str, dict[str, dict[str, float]]] = {}
    print("\n══ PER CONDITION ══")
    for condition in conditions:
        name = condition["name"]
        lengths[name] = read_truth_length_stats(args.suite / name / "truth_fragment_lengths.tsv")
        bam = args.suite / name / "sim_oracle.bam"
        if not bam.exists():
            raise SystemExit(f"missing oracle BAM: {bam}")
        origins = read_gdna_origins(bam)
        census[name] = gdna_census(origins, ref_lengths, exon_iv, span_iv, probe_iv)
        stats = lengths[name]
        gdna_mean = stats.get("gdna", {}).get("mean", float("nan"))
        rna_mean = stats.get("rna", {}).get("mean", float("nan"))
        print(
            f"  {name:<48} gdna n={census[name]['n_total']:>9,} mean={gdna_mean:7.2f}"
            f"  rna mean={rna_mean:7.2f}  refs={len(census[name]['per_ref'])}"
        )

    report = Report()

    # ── G-S1 ────────────────────────────────────────────────────────────────────────────────────
    leaked = {
        name: {ref: n for ref, n in entry["per_ref"].items() if ref not in genomic and n}
        for name, entry in census.items()
    }
    n_leaked = sum(sum(refs.values()) for refs in leaked.values())
    n_leaked_refs = len({ref for refs in leaked.values() for ref in refs})
    report.gate(
        "G-S1 gDNA on an RNA-only reference",
        n_leaked == 0,
        f"{n_leaked:,} fragments on {n_leaked_refs} RNA-only references (must be 0)",
    )

    # ── G-S2 ────────────────────────────────────────────────────────────────────────────────────
    genomic_with_gdna = {
        ref
        for entry in census.values()
        for ref, n in entry["per_ref"].items()
        if ref in genomic and n
    }
    contaminated = [name for name, entry in census.items() if entry["n_total"]]
    per_condition_ok = all(
        len({ref for ref, n in census[name]["per_ref"].items() if ref in genomic and n}) >= 2
        for name in contaminated
    )
    report.gate(
        "G-S2 genomic references carrying gDNA",
        len(genomic_with_gdna) >= 2 and per_condition_ok and bool(contaminated),
        f"{len(genomic_with_gdna)} of {len(genomic)} genomic references, "
        f"every gDNA condition >= 2: {per_condition_ok} (need >= 2)",
    )

    # ── G-S3 and G-S5: pair each base condition's capture arms ──────────────────────────────────
    pairs: dict[tuple, dict[str, str]] = defaultdict(dict)
    for condition in conditions:
        pairs[base_key(condition)][condition["capture_label"]] = condition["name"]

    g3_rows: list[tuple[str, float, float]] = []
    g5_rows: list[tuple[str, float, float]] = []
    for key, arms in sorted(pairs.items(), key=lambda kv: str(kv[0])):
        if "on" not in arms or "off" not in arms:
            continue
        off, on = lengths[arms["off"]], lengths[arms["on"]]
        label = f"{key[0]} ss{key[1]:.2f}"
        if census[arms["on"]]["n_total"] and census[arms["off"]]["n_total"]:
            g3_rows.append((label, off["gdna"]["mean"], on["gdna"]["mean"]))
            gap_off = abs(off["gdna"]["mean"] - off["rna"]["mean"])
            gap_on = abs(on["gdna"]["mean"] - on["rna"]["mean"])
            g5_rows.append((label, gap_off, gap_on))

    print("\n══ G-S3  gDNA mean length, capture off -> on ══")
    for label, off_mean, on_mean in g3_rows:
        print(f"  {label:<24} {off_mean:7.2f} -> {on_mean:7.2f}   {on_mean - off_mean:+7.2f} bp")
    report.gate(
        "G-S3 gDNA mean length rises under capture",
        bool(g3_rows) and all(on > off for _, off, on in g3_rows),
        f"{sum(1 for _, off, on in g3_rows if on > off)}/{len(g3_rows)} conditions strictly greater",
    )

    print("\n══ G-S5  |mu_gdna - mu_rna|, capture off -> on ══")
    for label, gap_off, gap_on in g5_rows:
        print(f"  {label:<24} {gap_off:7.2f} -> {gap_on:7.2f}   {gap_on - gap_off:+7.2f} bp")
    report.gate(
        "G-S5 the gDNA<->RNA length gap narrows under capture",
        bool(g5_rows) and all(gap_on < gap_off for _, gap_off, gap_on in g5_rows),
        f"{sum(1 for _, off, on in g5_rows if on < off)}/{len(g5_rows)} conditions strictly narrower",
    )

    # ── G-S4 ────────────────────────────────────────────────────────────────────────────────────
    print("\n══ G-S4  gDNA mean length by probe overlap (capture arms only) ══")
    g4_rows: list[tuple[str, float, float]] = []
    for condition in conditions:
        name = condition["name"]
        if condition["capture_label"] != "on" or not census[name]["n_total"]:
            continue
        means = census[name]["means"]
        on_target, off_target = means["on_target"], means["off_target"]
        print(
            f"  {name:<48} on-target n={on_target['n']:>8,} mean={on_target['mean']:7.2f}"
            f"   off-target n={off_target['n']:>8,} mean={off_target['mean']:7.2f}"
        )
        if on_target["n"] and off_target["n"]:
            g4_rows.append((name, off_target["mean"], on_target["mean"]))
    report.gate(
        "G-S4 on-target gDNA is longer than off-target",
        bool(g4_rows) and all(on > off for _, off, on in g4_rows),
        f"{sum(1 for _, off, on in g4_rows if on > off)}/{len(g4_rows)} conditions strictly longer",
    )

    print("\n══ diagnostic: gDNA mean length by the territory the START lands in ══")
    print("  ⚠ NOT a gate. An intronic start conditioned on capture selects LONG fragments by")
    print("    geometry (it had to reach the probe); an exonic start does not. See the docstring.")
    for condition in conditions:
        name = condition["name"]
        if not census[name]["n_total"]:
            continue
        means = census[name]["means"]
        cells = "  ".join(
            f"{label} n={means[label]['n']:>8,} mean={means[label]['mean']:7.2f}"
            for label in ("intergenic", "intronic", "exonic")
        )
        print(f"  {name:<48} {cells}")

    # ── G-S6 ────────────────────────────────────────────────────────────────────────────────────
    over = sum(entry["over_reference"] for entry in census.values())
    report.gate(
        "G-S6 gDNA fragments longer than their reference",
        over == 0,
        f"{over:,} fragments run past their reference end (must be 0)",
    )

    return 0 if report.emit() else 1


if __name__ == "__main__":
    sys.exit(main())
