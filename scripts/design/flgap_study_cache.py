#!/usr/bin/env python
"""⭐⭐⭐ **SCAN ONCE, INTERROGATE MANY TIMES — the flgap pair, ready to question in a second.**

Everything downstream of the BAM that a prior arm needs, cached: the scored ``multi_loci``, the
calibration, the DRAINED payloads and their lifted origin partitions, the oracle's truth masses and
per-component shares, the F arm, and ⭐ the per-unit ``(true origin, is_spliced)`` pair the **Fo** arm
joins on — for BOTH the undrained and the drained scoring stage, because a drain changes the
calibration and therefore the candidate sets. Building costs ~5 min per condition (two scans + a
calibration + a drain + a ~30 s read-name walk); loading costs ~1 s.

⭐ **This is what makes the post-calibration work cheap.** The interrogation that found the
boundary-attribution defect ran a dozen questions against this cache; each would have been a five-minute
scan otherwise. ⛔ ``priors.py`` is deliberately **excluded from the cache key**, so editing the
assembler does NOT invalidate a scan — testing a change to ``assemble_priors`` is a one-second loop.

⛔ **KEYED, because a silently stale cache here would poison every number read off it.** The key is
the scan cache's own manifest PLUS a content hash of every source file that produces these artifacts —
see ``_KEYED_SOURCES``. A mismatch REBUILDS; it never warns and proceeds. Same contract as
``read_scan_cache``, whose refusal is what taught this repo that ``reach`` is covered by no other hash.

⭐⭐ **``priors.py`` STAYS OUT OF THE KEY, AND THAT ONLY WORKS IF NOTHING PRODUCED BY IT IS STORED.**
The blob therefore holds ``(cal, multi_loci, override, shares)`` and the raw per-region ``start_*``
counts, and every ARM — P, O, S, F — is assembled at read time from those. ⛔ It used to store the
assembled ``p_arm`` and the projected ``f_gdna``: both come out of ``priors.py``, so an assembler edit
served a fresh O beside a stale P and F and the comparison between them was meaningless. That is
``TRAPS: a-guard-outlives-its-divisor``'s shape — the exclusion was justified by "the arms are rebuilt
on load" and then artifacts stopped being rebuilt on load.

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
#: rebuilt on load so that editing the assembler does not invalidate a 5-minute scan, which is only
#: sound because nothing ``priors.py`` produces is stored (see the module docstring).
#: ⚠ ``_oracle.py`` and ``prior_vs_oracle.py`` were missing until 2026-08-08, so the truth masses, the
#: shares and the per-unit origin — all of them produced there — could go stale with the key still
#: reading "ok". They are hashed WHOLE rather than per-function: over-triggering costs a rebuild,
#: under-triggering costs a wrong number. ⭐ If that loop ever becomes painful, move the four producing
#: functions into their own small module and key that, rather than narrowing this list.
_KEYED_SOURCES = (
    "src/rigel/scan_payload.py",
    "src/rigel/second_pass.py",
    "src/rigel/pipeline.py",
    "src/rigel/locus.py",
    "src/rigel/scan.py",
    "src/rigel/sim/read_name.py",
    "tests/calibration/_oracle.py",
    "scripts/design/prior_vs_oracle.py",
)


def _source_key() -> str:
    h = hashlib.blake2b(digest_size=8)
    # ⛔⛔ THIS FILE IS PART OF ITS OWN KEY. What the blob CONTAINS is decided here, so a build that
    # starts storing a new artifact must invalidate the blobs that lack it — otherwise `load` returns
    # an old dict and the caller gets a KeyError, or worse, silently skips an arm. Landing
    # ``unit_origin`` (2026-08-08) is exactly that case: the calibration sources were untouched and
    # every stale blob still verified "ok".
    # ⚠ The price is that editing a comment here costs a ~20-minute rebuild. Accepted: this file is
    # small and rarely edited, and the failure it prevents is silent.
    h.update(Path(__file__).read_bytes())
    for rel in _KEYED_SOURCES:
        h.update((_REPO / rel).read_bytes())
    for f in sorted((_REPO / "src" / "rigel" / "calibration").rglob("*.py")):
        if f.name != "priors.py":
            h.update(f.read_bytes())
    return h.hexdigest()


def _key(panel: str, cond: str) -> dict:
    manifest = RUNS / "suite" / panel / "oracle_cache" / cond / "_main" / "manifest.json"
    return {"source": _source_key(), "scan": json.loads(manifest.read_text())}


def _start_counts(oracle, origins) -> dict:
    """``{"gdna", "rna"}`` per-region first-base counts — the RAW input to ``F``'s projection."""
    starts = {k: np.asarray(oracle.parts[k].region_start_count, np.float64) for k in origins}
    return {"gdna": starts["gdna"], "rna": starts["mrna"] + starts["nrna"]}


def _gate_p_is_recomputable(priors_mod, fields, p_arm, cal, ra, ml) -> None:
    """⛔ The blob stores no ``P``, so a reader rebuilds it as ``assemble_priors(cal, ra, ml)``. That is
    only a legitimate substitution if it reproduces the call PRODUCTION made, byte for byte — which is
    what this asserts at build time, once, on the real objects.

    ⭐ Without it, dropping ``p_arm`` from the blob would be an unproven claim that the shipped prior is
    a pure function of the three things that ARE stored (TRAPS: byte-identity-gate)."""
    rebuilt = priors_mod.assemble_priors(cal, ra, ml)
    bad = [f for f in fields if not np.array_equal(getattr(rebuilt, f), getattr(p_arm, f))]
    if bad:
        raise AssertionError(
            f"assemble_priors(cal, ra, ml) does not reproduce the prior quant_from_buffer built: "
            f"{bad}. P cannot be rebuilt from the cache and must be stored again."
        )


def build(panel: str, cond: str) -> dict:
    """Scan, calibrate, score, drain and lift — everything an arm needs, undrained AND drained."""
    from _oracle import ORIGINS as _OR, OracleTruth, check_walk_alignment, frag_id_origins
    from prior_vs_oracle import (
        PRIOR_FIELDS,
        _calibrate_and_prior,
        _oracle_parts,
        unit_origins,
    )

    import rigel.calibration.priors as PRIORS

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
        # ⭐ The frag_id -> true-origin key, walked ONCE: it is a property of the BAM alone, so the
        # undrained and drained arms share it (the drain moves deposits, never a read name).
        frag_origin, walk = frag_id_origins(bam, scan)
        check_walk_alignment(walk, stats)
        cal, _fl, ml, p_arm, units = _calibrate_and_prior(payload, sm, buf, stats, index, ra, cfg)
        _gate_p_is_recomputable(PRIORS, PRIOR_FIELDS, p_arm, cal, ra, ml)
        out = {
            "panel": panel,
            "cond": cond,
            "cal": cal,
            "multi_loci": ml,
            "override": oracle.override_masses(ra),
            "shares": oracle.component_shares(),
            # ⭐ RAW per-region start counts, not the projected F — the projection lives in priors.py,
            # which is outside the key, so it must run at read time. See the module docstring.
            "starts": _start_counts(oracle, _OR),
            "n_held": int(payload.deferred.n_fragments),
            # ⭐ The Fo arm's inputs, PRE-JOINED to keep the blob small: int8 + bool per unit rather
            # than an int64 frag_id per unit plus the whole 10 M-entry walk. ⛔ ``walk`` itself is
            # carried because its per-origin totals are the denominators for "how many fragments of
            # this origin never became an EM candidate at all", which ``overlap_truth`` reports.
            "walk": walk,
            "unit_origin": unit_origins(units["frag_ids"], frag_origin),
            "unit_is_spliced": np.asarray(units["is_spliced"], bool),
            "n_units": int(units["n_units"]),
        }

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
                    parts[k], ch, region_types=lift["region_types"], sj=lift["sj"]
                )
                for k, ch in zip(ORIGINS, lifted)
            }
            # ⛔ gDNA cannot splice. Reported, never swallowed — it is what makes flgap_short the
            # admissible panel for a drained measurement and flgap_long a cross-check.
            leak = sum(
                int(np.asarray(getattr(parts_d["gdna"], b), np.int64).sum())
                for b in ("boundary_spliced_count", "sj_count")
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
            cal_d, _fld, ml_d, p_d, units_d = _calibrate_and_prior(
                payload_d, sm2, buf2, s2, index, ra, cfg
            )
            out.update(
                {
                    "cal_d": cal_d,
                    "multi_loci_d": ml_d,
                    "override_d": oracle_d.override_masses(ra),
                    "shares_d": oracle_d.component_shares(),
                    "lift_ambiguous": int(n_amb),
                    "gdna_spliced_leak": int(leak),
                    # ⚠ The DRAINED run scores its own fragments — the drain changes the calibration,
                    # hence the FL models, hence the candidate sets — so it gets its own unit arrays.
                    # The origin key is shared: a drain moves deposits, never a read name.
                    "unit_origin_d": unit_origins(units_d["frag_ids"], frag_origin),
                    "unit_is_spliced_d": np.asarray(units_d["is_spliced"], bool),
                    "n_units_d": int(units_d["n_units"]),
                    "starts_d": _start_counts(oracle_d, _OR),
                }
            )
            _gate_p_is_recomputable(PRIORS, PRIOR_FIELDS, p_d, cal_d, ra, ml_d)
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
