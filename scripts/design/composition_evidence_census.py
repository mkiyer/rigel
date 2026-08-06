"""P0 — how much of the library's unspliced mass reaches the solver with NO composition evidence?

       Hook: `bp_solver.node_sweep(_capture=...)`.

⭐ **THE QUESTION.** `node_init` gives every slot its own composition precision `tau_lam` from three live
sources — the structural lock, the intron-factory density deconvolution, and the strand Beta-Binomial.
A slot with `tau_lam == 0` and no structural lock has **no own evidence at all**: its gDNA/RNA split is
decided entirely by neighbour messages and the population prior. This script measures how much of the
library sits there, and splits it by the axes that say WHY.

⚠ **The strand source is identically zero at κ = ½** (the gDNA fraction cancels
out of the Beta-Binomial mean, verified to 5.6e-17) and is gated OFF on AMBIG slots by the Schur
argument (`node_init.py`, approach E). So the prediction under test is that the `ss0.50` conditions
carry materially more no-evidence mass than the `ss0.99` ones, concentrated on AMBIG slots.

⛔ **MASS-WEIGHTED, NEVER NODE-COUNT-WEIGHTED.** records a claim that
survived for months because it was bp-weighted (0.8738) when the estimator is mass-weighted (0.9596).
80.5 % of partition nodes carry zero fragments; counting them would drown the answer.

⭐ **No production code is changed.** `calibrate(_debug={})` already fills `_debug["capture"]` with
`_tau0_lam`, the per-slot counts and the strand-freedom masks — an inert diagnostic hook that exists.

    python scripts/design/composition_evidence_census.py --index IDX --cache-root CACHE_DIR
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import EDGE, NODE  # noqa: E402
from rigel.calibration.node_geometry import g1_locked  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

#: ⭐ The census asks the SOLVER's question by IMPORTING the solver's predicate, not by restating its
#: constant (`TRAPS.md` A11 — the previous arrangement had this number written out here, in
#: ``pass0_vs_oracle.py`` and in the solver, each claiming to match the others).
from rigel.calibration.node_init import has_own_composition_evidence  # noqa: E402


def census_one(index: TranscriptIndex, cache_dir: Path, inject_kappa: float | None = None) -> dict:
    """Calibrate one cached condition and return its composition-evidence census.

    ``inject_kappa`` is the FALSIFICATION handle: re-run with κ replaced and **nothing else**, by
    harvesting the condition's own fitted priors and substituting one field. A census that does not move
    when the strand channel is switched off is measuring something other than what it claims
    (`falsification_needs_perturbation`). ⭐ Verified 2026-07-31: injecting κ=0.5 into
    ``gdna100 ss0.99 capture_off`` reproduces ``gdna100 ss0.50 capture_off``'s numbers **exactly**
    (49.4 % / 63.6 %) — same payload, one variable, landing on the natural experiment's value.
    """
    cache = read_scan_cache(cache_dir, index)
    inputs = calibration_inputs(cache, index)
    debug: dict = {}
    injected = None
    if inject_kappa is not None:
        seed: dict = {}
        calibrate(**inputs, config=CalibrationConfig(), _debug=seed)
        injected = dataclasses.replace(
            seed["calibration_priors"], rna_sense_frac=float(inject_kappa)
        )
    result = calibrate(
        **inputs,
        config=CalibrationConfig(),
        _debug=debug,
        injected_priors=injected,
    )
    cap = debug["capture"]
    chain = debug["chain"]

    tau = np.asarray(cap["_tau0_lam"], np.float64)
    count = np.asarray(cap["count"], np.float64).sum(axis=1)  # unspliced mass per slot, both strands
    free_pos = np.asarray(cap["free_pos"], bool)
    free_neg = np.asarray(cap["free_neg"], bool)
    kind = np.asarray(chain.kind)

    is_node = kind == NODE
    is_edge = kind == EDGE
    #   struct_lock   = the G1 class (composition CERTAIN, pinned {0,0,1} at var 0)
    #   single_strand = free_pos XOR free_neg         (the strand λ-term is gated to these)
    # ⛔ struct_lock is BOTH AXES, and it comes from the ONE definition
    # (`node_geometry.g1_locked`) rather than being re-derived here. It was `(~solvable) & is_node`,
    # which filed every structurally-locked EDGE — an intergenic<->exon seam, where RNA cannot cross a
    # gene boundary — as a slot with NO EVIDENCE rather than as one that is certain.
    # ⚠⚠ NOT the same mask as `node_init.strand_evidence`'s own `struct_lock`, which is node-only ON
    # PURPOSE (it governs whether a slot may EMIT certainty into its messages). See `g1_locked`.
    struct_lock = g1_locked(free_pos, free_neg)
    single_strand = free_pos ^ free_neg
    ambig = free_pos & free_neg

    # ⭐ THE QUANTITY. Not "tau == 0" -- a structurally locked node is CERTAIN, not uninformed, and
    # lumping the two would report a pure-gDNA intergenic node as a failure of the solver.
    no_evidence = ~has_own_composition_evidence(tau) & (~struct_lock)

    total = float(count.sum())

    def share(mask: np.ndarray) -> float:
        return float(count[mask].sum()) / total if total > 0 else 0.0

    return {
        "condition": cache_dir.name,
        "kappa": float(result.rna_sense_frac),
        "f_gdna": _f_gdna(result),
        "mass": total,
        "no_evidence": share(no_evidence),
        "no_evidence_node": share(no_evidence & is_node),
        "no_evidence_edge": share(no_evidence & is_edge),
        "no_evidence_ambig": share(no_evidence & ambig),
        "no_evidence_ss": share(no_evidence & single_strand),
        "struct_lock": share(struct_lock),
        "mass_ambig": share(ambig),
        "mass_ss": share(single_strand),
        "mass_edge": share(is_edge),
        # ⚠ The conditional shares: of the mass on AMBIG slots, how much has no evidence? This is the
        # number the Schur gate predicts should be ~1 wherever the intron factory is silent.
        "ambig_no_ev_frac": _cond(count, no_evidence, ambig),
        "ss_no_ev_frac": _cond(count, no_evidence, single_strand),
    }


def _cond(count: np.ndarray, numer: np.ndarray, denom: np.ndarray) -> float:
    """Mass share of `numer` WITHIN `denom`; 0 where the denominator carries no mass (never floored)."""
    d = float(count[denom].sum())
    return float(count[numer & denom].sum()) / d if d > 0 else 0.0


def _f_gdna(result) -> float:
    """The library gDNA fraction as the LEDGER reports it. ⚠ This is an incidence-weighted sum, not a
    fragment count — quoted here only to key the row to the
    baseline table."""
    g = float(np.asarray(result.mass_gdna_node).sum() + np.asarray(result.mass_gdna_edge).sum())
    r = float(np.asarray(result.mass_rna_node).sum() + np.asarray(result.mass_rna_edge).sum())
    return g / (g + r) if (g + r) > 0 else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--cache-root", type=Path, required=True, help="dir of per-condition scan caches")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument(
        "--inject-kappa",
        type=float,
        default=None,
        help="FALSIFICATION: re-run with this kappa and nothing else changed. 0.5 switches the "
             "strand channel off; the no-evidence share must jump.",
    )
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    dirs = sorted(d for d in args.cache_root.iterdir() if (d / "payload.npz").exists())
    if args.conditions:
        dirs = [d for d in dirs if d.name in set(args.conditions)]

    rows = [census_one(index, d, args.inject_kappa) for d in dirs]

    print(
        f"\n{'condition':46s} {'kappa':>6s} {'f_gdna':>7s} "
        f"{'NO-EV':>7s} {'node':>7s} {'edge':>7s} {'AMBIG':>7s} {'1-str':>7s} {'lock':>6s}"
    )
    print("-" * 46 + " " + "-" * 63)
    for r in rows:
        print(
            f"{r['condition']:46s} {r['kappa']:6.4f} {r['f_gdna']:7.4f} "
            f"{r['no_evidence']:7.1%} {r['no_evidence_node']:7.1%} {r['no_evidence_edge']:7.1%} "
            f"{r['no_evidence_ambig']:7.1%} {r['no_evidence_ss']:7.1%} {r['struct_lock']:6.1%}"
        )

    print(f"\n{'condition':46s} {'mass AMBIG':>11s} {'mass 1-str':>11s} "
          f"{'no-ev|AMBIG':>12s} {'no-ev|1-str':>12s}")
    print("-" * 46 + " " + "-" * 50)
    for r in rows:
        print(
            f"{r['condition']:46s} {r['mass_ambig']:11.1%} {r['mass_ss']:11.1%} "
            f"{r['ambig_no_ev_frac']:12.1%} {r['ss_no_ev_frac']:12.1%}"
        )

    # ⭐ The prediction under test, evaluated rather than eyeballed.
    unst = [r for r in rows if abs(r["kappa"] - 0.5) < 0.05]
    strd = [r for r in rows if abs(r["kappa"] - 0.5) >= 0.05]
    if unst and strd:
        mu = np.mean([r["no_evidence"] for r in unst])
        ms = np.mean([r["no_evidence"] for r in strd])
        print(
            f"\nno-evidence mass share:  unstranded (kappa~0.5) {mu:.1%}   "
            f"stranded {ms:.1%}   ratio {mu / ms if ms > 0 else float('inf'):.2f}x"
        )


if __name__ == "__main__":
    main()
