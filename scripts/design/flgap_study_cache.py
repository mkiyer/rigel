#!/usr/bin/env python
"""⭐⭐⭐ **SCAN ONCE, INTERROGATE MANY TIMES — the flgap pair, ready to question in a second.**

Everything downstream of the BAM that a prior arm needs, cached: the scored ``multi_loci``, the
calibration, the DRAINED payloads and their lifted origin partitions, the oracle's truth masses and
per-component shares, and the F arm. Building costs ~5 min per condition (two scans + a calibration +
a drain); loading costs ~1 s.

⭐ **This is what makes the post-calibration work cheap.** The interrogation that found the
edge-attribution defect ran a dozen questions against this cache; each would have been a five-minute
scan otherwise. ⛔ ``priors.py`` is deliberately **excluded from the cache key**, so editing the
assembler does NOT invalidate a scan — testing a change to ``assemble_priors`` is a one-second loop.

⛔ **KEYED, because a silently stale cache here would poison every number read off it.** The key is
the scan cache's own manifest PLUS a content hash of every source file that produces these artifacts
(``src/rigel/calibration/**`` except ``priors.py``, plus ``scan_payload`` / ``second_pass`` /
``pipeline`` / ``locus``). A mismatch REBUILDS; it never warns and proceeds. Same contract as
``read_scan_cache``, whose refusal is what taught this repo that ``reach`` is covered by no other hash.

⚠ **The drained arm is admissible on flgap_short and marginal on flgap_long.** Sum-to-full is exact on
every drained partition, but the gDNA partition's spliced leak — gDNA cannot splice — is **1** record on
flgap_short against **1,010 / 5,827** on flgap_long. Both numbers are stored in the blob
(``gdna_spliced_leak``, ``lift_ambiguous``) and a caller must report them
(`TRAPS: an-equal-length-panel-defeats-the-lift`).

Usage::

    python scripts/design/flgap_study_cache.py                    # build/refresh all four conditions
    python scripts/design/flgap_study_cache.py --rebuild          # force
    python scripts/design/flgap_study_cache.py --list             # what is cached, and its key state

    # in an analysis script:
    import flgap_study_cache as FC
    st = FC.load("flgap_long", "gdna_g50_ss_0.50_nrna_none_capture_on")
    cal, ml, f = st["cal"], st["multi_loci"], st["f_gdna"]
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import pickle
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
for _p in (_REPO / "scripts" / "design", _REPO / "tests" / "calibration"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

RUNS = Path.home() / "Downloads" / "rigel_runs"
INDEX = RUNS / "suite" / "rigel_index"
CACHE = RUNS / "arms" / "flgap_study"
ORIGINS = ("gdna", "mrna", "nrna")

#: The four conditions the campaign gates on — the flgap PAIR, both capture arms.
#: ⛔ Never the ladder alone: its realised gDNA/RNA length gap is +1.5–2.1 % and it is structurally
#: blind to this whole axis.
CONDITIONS = [
    (panel, f"gdna_g50_ss_0.50_nrna_none_capture_{arm}")
    for panel in ("flgap_long", "flgap_short")
    for arm in ("on", "off")
]

#: Source files whose output is stored here. ⛔ ``priors.py`` is deliberately ABSENT — the arms are
#: rebuilt on load so that editing the assembler does not invalidate a 5-minute scan.
_KEYED_SOURCES = (
    "src/rigel/scan_payload.py",
    "src/rigel/second_pass.py",
    "src/rigel/pipeline.py",
    "src/rigel/locus.py",
)


def _source_key() -> str:
    h = hashlib.blake2b(digest_size=8)
    for rel in _KEYED_SOURCES:
        h.update((_REPO / rel).read_bytes())
    for f in sorted((_REPO / "src" / "rigel" / "calibration").rglob("*.py")):
        if f.name != "priors.py":
            h.update(f.read_bytes())
    return h.hexdigest()


def _key(panel: str, cond: str) -> dict:
    manifest = RUNS / "suite" / panel / "oracle_cache" / cond / "_main" / "manifest.json"
    return {"source": _source_key(), "scan": json.loads(manifest.read_text())}


def build(panel: str, cond: str) -> dict:
    """Scan, calibrate, score, drain and lift — everything an arm needs, undrained AND drained."""
    from _oracle import OracleTruth
    from prior_vs_oracle import _calibrate_and_prior, _oracle_parts, fragment_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag, scan_and_buffer
    from rigel.second_pass import drain as sp_drain, lift_choices

    index = TranscriptIndex.load(str(INDEX))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    bam = str(RUNS / "suite" / panel / cond / "sim_oracle.bam")
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))

    with tempfile.TemporaryDirectory() as work:
        stats, sm, buf, payload = scan_and_buffer(bam, index, scan)
        parts = _oracle_parts(
            bam, index, scan, cfg, work, cond, RUNS / "suite" / panel / "oracle_cache"
        )
        oracle = OracleTruth.from_parts(payload, parts)  # sum-to-full, HARD gate
        cal, _fl, ml, p_arm = _calibrate_and_prior(payload, sm, buf, stats, index, ra, cfg)
        out = {
            "panel": panel,
            "cond": cond,
            "cal": cal,
            "multi_loci": ml,
            "p_arm": p_arm,
            "override": oracle.override_masses(ra),
            "shares": oracle.component_shares(),
            "n_held": int(payload.deferred.n_fragments),
        }
        out["f_gdna"], out["f_rna_upper"], out["f_dropped"] = fragment_truth(oracle, ra, ml)

        if out["n_held"]:
            lift: dict = {}
            payload_d = _drain_side_buffer(
                payload, index, sm, seed=cfg.second_pass_seed, _lift=lift
            )
            lifted, n_amb = lift_choices(
                lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"]
            )
            parts_d = {
                k: sp_drain(
                    parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"]
                )
                for k, ch in zip(ORIGINS, lifted)
            }
            # ⛔ gDNA cannot splice. Reported, never swallowed — it is what makes flgap_short the
            # admissible panel for a drained measurement and flgap_long a cross-check.
            leak = sum(
                int(np.asarray(getattr(parts_d["gdna"], b), np.int64).sum())
                for b in ("edge_spliced_count", "sj_count")
                if hasattr(parts_d["gdna"], b)
            )
            oracle_d = OracleTruth(
                full=payload_d,
                parts=parts_d,
                read_counts={k: -1 for k in ORIGINS},
                n_ambiguous=n_amb,
            )
            # ⚠ A SECOND SCAN, and it buys one thing: a fresh fragment BUFFER. The first was consumed
            # by scoring, and a buffer is not re-scannable; the payload is cacheable and the buffer is not.
            s2, sm2, buf2, _p2 = scan_and_buffer(bam, index, scan)
            cal_d, _fld, ml_d, p_d = _calibrate_and_prior(payload_d, sm2, buf2, s2, index, ra, cfg)
            out.update(
                {
                    "cal_d": cal_d,
                    "multi_loci_d": ml_d,
                    "p_arm_d": p_d,
                    "override_d": oracle_d.override_masses(ra),
                    "shares_d": oracle_d.component_shares(),
                    "lift_ambiguous": int(n_amb),
                    "gdna_spliced_leak": int(leak),
                }
            )
            out["f_gdna_d"], out["f_rna_upper_d"], _ = fragment_truth(oracle_d, ra, ml_d)
    return out


def load(panel: str, cond: str, *, rebuild: bool = False) -> dict:
    """Cached :func:`build`. ⛔ Refuses and rebuilds on any key mismatch — never warns and proceeds."""
    d = CACHE / panel / cond
    blob, man = d / "study.pkl", d / "manifest.json"
    key = _key(panel, cond)
    if not rebuild and blob.is_file() and man.is_file():
        try:
            if json.loads(man.read_text()) == key:
                with blob.open("rb") as fh:
                    return pickle.load(fh)
            print(f"  [cache] KEY MISMATCH {panel}/{cond} — rebuilding", file=sys.stderr)
        except Exception as exc:  # noqa: BLE001 — a bad cache is a rebuild, not a crash
            print(f"  [cache] unreadable ({exc}) — rebuilding", file=sys.stderr)
    print(f"  [cache] building {panel}/{cond} (~5 min)...", file=sys.stderr)
    out = build(panel, cond)
    d.mkdir(parents=True, exist_ok=True)
    with blob.open("wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    man.write_text(json.dumps(key, indent=2))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rebuild", action="store_true", help="force a rebuild of every condition")
    ap.add_argument("--list", action="store_true", help="report cache state without building")
    args = ap.parse_args()

    if args.list:
        print(f"{'condition':<50} {'state':<12} {'size':>8}")
        print("-" * 72)
        for panel, cond in CONDITIONS:
            d = CACHE / panel / cond
            blob, man = d / "study.pkl", d / "manifest.json"
            if not blob.is_file():
                state, size = "MISSING", "-"
            else:
                size = f"{blob.stat().st_size / 1e6:.0f} MB"
                try:
                    state = "ok" if json.loads(man.read_text()) == _key(panel, cond) else "STALE"
                except Exception:  # noqa: BLE001
                    state = "UNREADABLE"
            print(f"{panel + '/' + cond:<50} {state:<12} {size:>8}")
        return 0

    for panel, cond in CONDITIONS:
        st = load(panel, cond, rebuild=args.rebuild)
        print(
            f"{panel}/{cond}: {len(st['multi_loci'])} loci  held={st['n_held']:,}  "
            f"lift_ambiguous={st.get('lift_ambiguous', 0):,}  "
            f"gdna_spliced_leak={st.get('gdna_spliced_leak', 0):,}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
