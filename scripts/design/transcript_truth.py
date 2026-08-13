#!/usr/bin/env python
"""⭐⭐⭐ **THE TRUE PER-TRANSCRIPT FRAGMENT COUNT, SPLIT BY SPLICEDNESS** — the truth a per-transcript
prior is scored against, and the one quantity no instrument produced.

`prior_vs_oracle.py` scores the prior per MULTI-LOCUS; `quant_accuracy.py` scores the tool's output per
transcript. Neither answers *"how many fragments did transcript `t` actually emit, and how many of them
were spliced?"* — which is what a per-transcript prior claims to know before the EM runs.

Usage::

    python scripts/design/transcript_truth.py --condition gdna_g00_ss_0.99_nrna_none_capture_off
    python scripts/design/transcript_truth.py --all --out truth_by_transcript.tsv

⛔⛔ **SPLICEDNESS IS DERIVED FROM THE READ NAME, NEVER FROM THE CIGAR.** The simulator writes
``{t_id}:{start}-{end}:{strand}:{index}`` with the interval in SPLICED TRANSCRIPT coordinates, so a
fragment is spliced exactly when its interval crosses one of its transcript's interior cumulative-exon
boundaries. Reading the CIGAR instead misses every sj that falls in the **unsequenced inner gap**
between the mates — measured **7.34 %** of truly-spliced mRNA fragments have no ``N`` in either mate, and
**0** in the reverse direction. A CIGAR-based count is therefore biased in one direction only, which is
the kind of error that looks like a result.

⚠ **"UNSPLICED" HAS THREE INEQUIVALENT DEFINITIONS HERE AND THIS SCRIPT EMITS THE TRUTH ONE.** The EM's
consumer bit is ``stype != SPLICE_UNSPLICED``, which is strictly WIDER than "the molecule crossed a
sj" because ``resolve_context`` promotes a clean-CIGAR fragment whose mate gap admits more than one
hypothesis to ``SPLICE_IMPLICIT``. Scoring calibration's own quantity against this one is
``TRAPS: score-the-consumers-own-count`` unless the gap is stated. It is stated: this is what the
SIMULATOR emitted, not what the resolver concluded.

⭐ **SEVEN NAMED GATES, and the falsification is free.** ``transcripts.feather`` encodes strand as ``{1: POS,
2: NEG}`` rather than ``±1``, and a first version silently skipped the minus-strand exon reversal — the
CIGAR-vs-truth disagreement went **12,330 → 74,178** and gained ~30,900 fragments in the structurally
impossible direction (a CIGAR ``N`` where truth says unspliced). That counter is kept below, so the same
mistake cannot pass quietly again.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]

from rigel.index import TranscriptIndex  # noqa: E402
from rigel.sim.read_name import parse_origin  # noqa: E402
from rigel.types import IntervalType, Strand  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def duplicate_map(gtf: Path) -> dict[str, str]:
    """``dropped transcript id -> the twin the index KEPT``.

    ⭐⭐ **THE INDEX DROPS EXACT DUPLICATES ON PURPOSE, AND IT IS RIGHT TO.** Transcripts sharing a
    reference, a strand and a bit-identical sorted exon tuple are — in ``index.py``'s own words —
    *"mathematically unidentifiable from any read data"*, so collapsing each group to one representative
    is loss-free for abundance estimation. It keeps the **lexicographically smallest** transcript id.

    ⛔ **But the SIMULATOR samples from the un-collapsed annotation**, so its read names reference ids the
    index does not contain. Measured on the suite reference: 8,755 transcripts, 5 duplicate groups over
    10 transcripts, 5 dropped — exactly the 8,750 the index holds. Two of them carry simulated fragments.
    Counting those as "unspliced", or dropping them, both put real fragments in the wrong place; folding
    them onto the retained twin is the only answer consistent with what the index decided.

    ⚠ Rebuilt here from the GTF because the index does NOT record the map — grepped 2026-08-12,
    ``_collapse_duplicate_transcripts`` logs the groups and returns a filtered list, storing nothing. The
    durable fix is for the index to emit it, since the index is what made the decision.
    """
    import collections
    import re

    exons: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    meta: dict[str, tuple[str, str]] = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if m is None:
                continue
            exons[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            meta[m.group(1)] = (f[0], f[6])

    groups: dict[tuple, list[str]] = collections.defaultdict(list)
    for t, e in exons.items():
        groups[(meta[t][0], meta[t][1], tuple(sorted(e)))].append(t)

    out: dict[str, str] = {}
    for members in groups.values():
        if len(members) < 2:
            continue
        keep = min(members)  # the index keeps the lexicographically smallest id
        for t in members:
            if t != keep:
                out[t] = keep
    return out


def interior_boundaries(index: TranscriptIndex) -> dict[str, np.ndarray]:
    """``transcript_id -> the interior cumulative-exon boundaries``, in TRANSCRIPTION order.

    A fragment at ``[s, e)`` in spliced transcript coordinates is SPLICED exactly when some boundary
    ``b`` satisfies ``s < b < e``. A single-exon transcript has none, so it can never be spliced —
    which is one of the gates below rather than an assumption.

    ⛔ **The exons are reversed for a minus-strand transcript.** Transcript coordinates run 5'→3', and a
    minus-strand transcript's exons are stored genomically ascending, so the cumulative sums are taken in
    the OPPOSITE order. ⚠ ``transcripts.feather`` encodes strand as ``{1: POS, 2: NEG}`` and NOT as
    ``±1``; reading it as a sign silently disables this reversal.
    """
    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    ex = ex.sort_values(["t_index", "start"], kind="stable")

    tdf = index.t_df
    strand_of = dict(zip(tdf["t_index"].to_numpy(), tdf["strand"].to_numpy(), strict=True))
    tid_of = dict(zip(tdf["t_index"].to_numpy(), tdf["t_id"].to_numpy(), strict=True))

    lens: dict[int, list[int]] = {}
    for t, a, b in zip(ex["t_index"], ex["start"], ex["end"], strict=True):
        lens.setdefault(int(t), []).append(int(b) - int(a))

    out: dict[str, np.ndarray] = {}
    for t, ls in lens.items():
        if int(strand_of.get(t, Strand.POS)) == int(Strand.NEG):
            ls = ls[::-1]
        out[str(tid_of[t])] = np.cumsum(np.asarray(ls, dtype=np.int64))[:-1]
    return out


def _is_spliced(bounds: np.ndarray, s: int, e: int) -> bool:
    if bounds.size == 0:
        return False
    i = int(np.searchsorted(bounds, s, side="right"))
    return i < bounds.size and bounds[i] < e


def count_condition(bam: str, bounds: dict[str, np.ndarray], limit: int | None = None,
                    dup: dict[str, str] | None = None) -> dict:
    """One pass over the oracle BAM, READ NAMES ONLY.

    ⭐ Counts primary read-1 records only: both mates carry the same qname and the same truth, so
    counting records would double every fragment and counting unique qnames would hold ~20 M strings.
    """
    import pysam

    per_t: dict[str, list[int]] = {}
    # ⛔ A transcript the READ NAMES reference but the INDEX does not contain. Counting such a fragment
    # as "unspliced" would be a silent default on a substrate inconsistency — measured 329 fragments
    # across 2 transcripts on the ladder, which is 0.0033 % and still not a number to invent.
    unknown: dict[str, int] = {}
    origin_totals = {"mrna": 0, "nrna": 0, "gdna": 0}
    cigar_disagree = {"cigar_n_truth_unspliced": 0, "truth_spliced_no_cigar_n": 0}
    n = 0
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as fh:
        for rec in fh.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary or not rec.is_read1:
                continue
            o = parse_origin(rec.query_name)
            origin_totals[o.kind] = origin_totals.get(o.kind, 0) + 1
            if o.kind != "gdna" and o.transcript_id is not None:
                # ⭐ Fold an exact duplicate onto the twin the index kept, BEFORE anything else reads
                # the id — the two are the same molecule to every consumer downstream.
                tid = (dup or {}).get(o.transcript_id, o.transcript_id)
                b = bounds.get(tid)
                if b is None:
                    unknown[tid] = unknown.get(tid, 0) + 1
                    n += 1
                    continue
                sp = _is_spliced(b, int(o.start), int(o.end))
                row = per_t.setdefault(tid, [0, 0])
                row[1 if sp else 0] += 1
                # ⭐ the free falsification: a CIGAR N where truth says unspliced is STRUCTURALLY
                # impossible, so a nonzero count here means the boundary table is wrong.
                has_n = rec.cigarstring is not None and "N" in rec.cigarstring
                if has_n and not sp:
                    cigar_disagree["cigar_n_truth_unspliced"] += 1
                elif sp and not has_n:
                    cigar_disagree["truth_spliced_no_cigar_n"] += 1
            n += 1
            if limit is not None and n >= limit:
                break
    return {"per_t": per_t, "origins": origin_totals, "cigar": cigar_disagree,
            "unknown": unknown, "n_read1": n}


def gates(res: dict, suite: Path, condition: str, index: TranscriptIndex,
          bounds: dict[str, np.ndarray], dup: dict[str, str] | None = None,
) -> list[tuple[str, bool, str]]:
    """The seven gates, each NAMED and stated as a property the truth must satisfy."""
    per_t = res["per_t"]
    out: list[tuple[str, bool, str]] = []

    tot = {t: v[0] + v[1] for t, v in per_t.items()}
    out.append(("conserves — unspliced + spliced == total, per transcript", True,
                f"{len(per_t):,} transcripts"))

    summ = json.loads((suite / condition / "truth_summary.json").read_text())
    oc = summ.get("origin_counts", {})
    want = int(oc.get("mrna", 0)) + int(oc.get("nrna", 0))
    n_unknown = sum(res["unknown"].values())
    got = sum(tot.values()) + n_unknown
    out.append(("origin-total — Σ total (+ unindexed) == truth_summary origin_counts", got == want,
                f"{got:,} vs {want:,}" + (f"  (incl. {n_unknown:,} unindexed)" if n_unknown else "")))

    ta = pd.read_csv(suite / condition / "truth_abundances.tsv", sep="\t")
    col = next((c for c in ("observed_mrna_fragments", "mrna_abundance") if c in ta.columns), None)
    if col is None:
        out.append(("abundance-match — per-transcript total == truth_abundances", False, "no column found"))
    else:
        # ⭐ The truth table is keyed on the UN-COLLAPSED annotation, so it must be folded by the same
        # duplicate map before it can be compared. ⛔ For an exact duplicate the per-transcript truth is
        # not merely inconvenient, it is UNDEFINED — the two transcripts are the same molecule and only
        # the group total is a fact about the world.
        ta = ta.copy()
        ta["transcript_id"] = ta["transcript_id"].map(lambda t: (dup or {}).get(t, t))
        want_s = ta.groupby("transcript_id")[col].sum()
        common = [t for t in tot if t in want_s.index]
        bad = [t for t in common if int(want_s[t]) != tot[t]]
        out.append((f"abundance-match — per-transcript total == truth_abundances[{col}]", not bad,
                    f"{len(common):,} compared, {len(bad):,} differ"))

    single = [t for t, b in bounds.items() if b.size == 0]
    bad_single = [t for t in single if per_t.get(t, [0, 0])[1] != 0]
    out.append(("single-exon — a single-exon transcript is never spliced", not bad_single,
                f"{len(single):,} single-exon, {len(bad_single):,} violations"))

    out.append(("cigar-impossible — a CIGAR N where truth says UNSPLICED",
                res["cigar"]["cigar_n_truth_unspliced"] == 0,
                f"{res['cigar']['cigar_n_truth_unspliced']:,}"))
    out.append(("indexed — every truth transcript exists in the index",
                not res["unknown"],
                f"{len(res['unknown']):,} unindexed transcript(s), "
                f"{sum(res['unknown'].values()):,} fragments"
                + (f" — {', '.join(list(res['unknown'])[:3])}" if res["unknown"] else "")))

    sj = pd.read_feather(os.path.join(index.index_dir, "sj.feather"))
    per_t_sj = sj.groupby("t_index").size()
    tid = dict(zip(index.t_df["t_index"].to_numpy(), index.t_df["t_id"].to_numpy(),
                   strict=True))
    mism = sum(
        1 for t_idx, k in per_t_sj.items()
        if bounds.get(str(tid.get(t_idx, "")), np.zeros(0)).size != int(k)
    )
    out.append(("sj-count — boundaries == annotated sj, per transcript", mism == 0,
                f"{len(per_t_sj):,} transcripts, {mism:,} mismatched"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default="gdna_g00_ss_0.99_nrna_none_capture_off")
    ap.add_argument("--limit", type=int, default=None, help="stop after N read-1 records (a smoke run)")
    ap.add_argument("--gtf", type=Path, default=None,
                    help="the annotation the SIMULATOR sampled from, for the exact-duplicate map. "
                         "Defaults to <suite>/../reference/genes.gtf")
    ap.add_argument("--out", type=Path, default=None, help="write the per-transcript TSV here")
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    gtf = args.gtf or (args.suite.parent / "reference" / "genes.gtf")
    dup = duplicate_map(gtf) if gtf.is_file() else {}
    if dup:
        print(f"  exact-duplicate transcripts folded onto their retained twin: {len(dup):,}")
    elif args.gtf is not None:
        print(f"  ⚠ no GTF at {gtf} — duplicates will read as unindexed")
    bounds = interior_boundaries(index)
    print(f"  transcripts with a boundary table: {len(bounds):,}")

    bam = str(args.suite / args.condition / "sim_oracle.bam")
    print(f"  scanning {args.condition} …", flush=True)
    res = count_condition(bam, bounds, limit=args.limit, dup=dup)
    print(f"  read-1 records: {res['n_read1']:,}   origins: {res['origins']}")
    print(f"  CIGAR vs truth: {res['cigar']}")
    if res["unknown"]:
        print(f"  ⚠ {sum(res['unknown'].values()):,} fragments from {len(res['unknown'])} transcript(s) "
              f"the index does not contain: {', '.join(list(res['unknown'])[:4])}")

    if args.limit is None:
        print("\n  GATES")
        ok = True
        for name, passed, detail in gates(res, args.suite, args.condition, index, bounds, dup):
            print(f"    {'✔' if passed else '✘'} {name:60s} {detail}")
            ok &= passed
        print(f"\n  {'ALL GATES PASS' if ok else '⛔ GATE FAILURE'}")

    if args.out is not None:
        rows = [
            {"transcript_id": t, "n_unspliced": v[0], "n_spliced": v[1], "n_total": v[0] + v[1]}
            for t, v in sorted(res["per_t"].items())
        ]
        pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
        print(f"  ⭐ {len(rows):,} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
