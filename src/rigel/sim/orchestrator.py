"""The one condition-grid loop, shared by both simulator frontends.

``suite.main`` (synthetic mini-genome suites) and ``whole_genome.run_simulation`` (simulate from an
existing reference) sweep the same grid — nascent x gDNA rate x gDNA strand overdispersion x strand
specificity x capture — so it is implemented once here and a new axis is wired once.

Per condition, :func:`run_condition_grid` deep-copies the transcripts back to the base abundances,
applies the nascent mode, resolves the RNA/gDNA fragment split (:func:`resolve_depths`), runs a
:class:`~rigel.sim.wgs_engine.WholeGenomeSimulator` on a capture sampler shared by every condition
using that panel, writes the post-capture truth from the oracle BAM, and appends one manifest entry
for the caller to write. Conditions are resumable: a condition whose existence key is already on
disk is skipped, the key being the oracle BAM whenever one is requested, since that is the artifact
the instruments read.

Every condition draws a distinct seed through :func:`capture_paired_condition_seed`, which keys on
``(gdna, strand specificity)`` only, so the capture, overdispersion and nascent variants of one base
condition are paired for controlled comparison while the main axes are decorrelated.
"""

from __future__ import annotations

import copy
import hashlib
import time
from dataclasses import dataclass, replace
from itertools import product
from pathlib import Path

from .manifest import condition_dir_name
from .truth import write_post_capture_truth
from .capture import CaptureSampler
from .whole_genome import (
    WholeGenomeSimulator,
    apply_nrna_fragment_share,
    apply_nrna_ratio,
    apply_sparse_nrna,
    write_truth_abundances,
)

__all__ = [
    "stable_seed",
    "capture_paired_condition_seed",
    "run_condition_grid",
    "ConditionDepths",
    "resolve_depths",
]


@dataclass(frozen=True, slots=True)
class ConditionDepths:
    """How one condition's fragment budget splits between RNA and gDNA. The mature/nascent split
    inside ``n_rna`` is never imposed: nascent entities are transcripts in the same multinomial, so
    their share follows from molecules and lengths and is read off the realised origin counts
    afterwards."""

    n_rna: int
    n_gdna: int

    @property
    def total(self) -> int:
        return self.n_rna + self.n_gdna


def resolve_depths(sim, *, gdna_rate: float) -> ConditionDepths:
    """Split one condition's fragment budget between RNA and gDNA. Two modes, by config:

    * ``n_total_fragments is None`` (the default) — the RNA depth is fixed and gDNA is added on top
      at ``rate x n_rna``. A ``rate`` of 1.0 is therefore a 50 % gDNA library at *twice* the depth
      of a ``rate`` of 0.
    * ``n_total_fragments`` set — the total is fixed and ``rate`` decides only the split:
      ``n_rna = total/(1+rate)``, ``n_gdna = total − n_rna``.

    The second mode is what reaches the high-gDNA end at all: under the first, a 98 % gDNA library
    is ``rate = 49``, i.e. 490 M fragments against a 10 M RNA depth. It is also the more faithful
    model, since a sequencing run has a fixed budget and the contamination fraction decides how that
    budget is split rather than how much extra is generated. The accepted trade-off is that the RNA
    side thins as gDNA rises, so per-transcript accuracy degrades at high contamination — a property
    of such libraries rather than an artifact.

    Nascent comes out of the RNA share and never on top of it, structurally: the nascent entities
    are rows of the one RNA multinomial, so ``n_rna`` is the whole RNA budget and the mature/nascent
    split inside it is realised, not imposed. Raises when the requested split leaves no RNA at all.

    Gates: ``tests/test_sim_fixed_total_depth.py``.
    """
    total = getattr(sim, "n_total_fragments", None)

    if total is None:
        # ── additive: RNA depth fixed, gDNA added on top ──────────────────────────────────────
        n_rna = int(sim.n_rna_fragments)
        return ConditionDepths(n_rna, round(gdna_rate * n_rna))

    # ── fixed total: `rate` decides the split only ────────────────────────────────────────────
    total = int(total)
    n_rna = round(total / (1.0 + float(gdna_rate)))
    # gDNA is the remainder, so the total is conserved by construction rather than by luck. What
    # breaks conservation is computing the two shares independently with truncation (`int()` on
    # both), which drifts by a fragment at some rates; rounding both and taking the remainder does
    # not. That truncation is the perturbation the gate uses.
    n_gdna = total - n_rna
    if n_rna <= 0:
        raise ValueError(
            f"n_total_fragments={total} at gdna_rate={gdna_rate} leaves {n_rna} RNA fragments. "
            "A condition with no RNA is not an RNA-seq library; raise the total or lower the rate."
        )
    return ConditionDepths(n_rna, n_gdna)


def stable_seed(base_seed: int, *parts: object) -> int:
    """Derive a reproducible 32-bit seed from a base seed and string parts."""
    text = "\0".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") & 0xFFFF_FFFF


def capture_paired_condition_seed(
    base_seed: int,
    gdna_label: str,
    strand_specificity: float,
) -> int:
    """Seed shared by the capture and nascent variants of one ``(gdna, ss)`` base condition.

    The nascent label is deliberately not an argument, so the variants of one base condition start
    from the same stream. Nascent rows are drawn in the same multinomial as the mature rows,
    so a nascent-on cell and its nascent-off twin share the seed but not a bit-identical mature
    stream: turning nascent on re-allocates the RNA budget, as it physically must. The gDNA stream
    is unaffected.
    """
    seed_name = condition_dir_name(gdna_label, strand_specificity, "_paired")
    return stable_seed(base_seed, seed_name)


def run_condition_grid(
    *,
    outdir: Path,
    genome_path: Path,
    transcripts: list,
    base_abundances: list[tuple[float, float]],
    sim,
    gdna,
    nrna,
    genomic_refs: list[str],
    gdna_pairs: list[tuple[str, float]],
    gdna_od_pairs: list[tuple[str, float]],
    strand_specificities: list[float],
    nrna_pairs: list[tuple[str, str, object, int]],
    capture_scenarios: list,
    include_capture_in_names: bool,
    base_seed: int,
    oracle_bam: bool = True,
    skip_existing: bool = True,
    emit_fastq: bool = True,
    selected_conditions: set[str] | None = None,
    capture_meta_by_label: dict[str, dict] | None = None,
) -> list[dict]:
    """Run the full condition grid and return the per-condition manifest entries.

    ``nrna_pairs`` entries are ``(label, mode, value, index)`` (see ``whole_genome._build_nrna_pairs``).
    The modes reach the entities by two different routes and the difference is the point:
    ``additive_ratio`` and ``fragment_share`` pool each entity's molecules from its contributors
    (`whole_genome.assign_nrna_to_entities`, so nascent tracks mature and cannot exceed it), while
    ``sparse`` writes each entity's absolute abundance directly (`whole_genome.apply_sparse_nrna`, so
    most entities get exactly zero and the rest are independent of the mature level). ``file`` leaves
    the loaded abundances alone. Every mode then shares one thing: the RNA budget is one multinomial
    over mature and nascent rows, so the fragment split is realised rather than allocated.
    ``capture_meta_by_label`` supplies the suite's probe-provenance fields per
    capture label (empty for the reference-driven path). The caller writes the manifest.
    """
    capture_meta_by_label = capture_meta_by_label or {}
    # One CaptureSampler per capture scenario, built once and reused by every condition that shares
    # it. The probe layout and the per-width partition depend only on the panel and the templates —
    # not on abundance, gDNA rate, strand or nascent — so rebuilding per condition would recompute
    # minutes of identical numbers each time.
    from rigel.index import load_reference_lengths

    _ref_lengths = load_reference_lengths(genome_path)
    samplers = {
        scenario.label: CaptureSampler.from_config(scenario.config, transcripts, _ref_lengths)
        for scenario in capture_scenarios
    }
    conditions: list[dict] = []
    cond_num = 0
    total = (
        len(nrna_pairs)
        * len(gdna_pairs)
        * len(gdna_od_pairs)
        * len(strand_specificities)
        * len(capture_scenarios)
    )

    for nrna_label, nrna_mode, nrna_value, nrna_index in nrna_pairs:
        cond_transcripts = copy.deepcopy(transcripts)
        for t, (base_mrna, base_nrna) in zip(cond_transcripts, base_abundances):
            t.abundance = base_mrna
            t.nrna_abundance = base_nrna

        nrna_ratio: float | None = None
        nrna_abundance_range: tuple[float, float] | None = None
        nrna_share: float | None = None
        if nrna_mode == "additive_ratio":
            nrna_ratio = float(nrna_value or 0.0)
            apply_nrna_ratio(cond_transcripts, nrna_ratio)
        elif nrna_mode == "fragment_share":
            # the config states the nascent share of RNA fragments; the molecular ratio is solved
            # from the annotation (`whole_genome.apply_nrna_fragment_share`) and recorded per condition
            nrna_share = float(nrna_value or 0.0)
            nrna_ratio = apply_nrna_fragment_share(cond_transcripts, nrna_share, sim)
        elif nrna_mode == "sparse":
            # nascent is absent from most gene spans and independent of the mature level where it
            # is present; the fragment share is emergent and recorded below
            nrna_abundance_range = tuple(nrna_value)  # type: ignore[arg-type]
            nrna_ratio = apply_sparse_nrna(
                cond_transcripts,
                nrna_abundance_range,
                on_fraction=nrna.on_fraction,
                seed=nrna.seed + nrna_index,
            )

        molecular_truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        write_truth_abundances(cond_transcripts, outdir / molecular_truth_name)

        for gdna_label, gdna_rate in gdna_pairs:
            for (gdna_od_label, gdna_od), strand_spec in product(
                gdna_od_pairs, strand_specificities
            ):
                for capture_scenario in capture_scenarios:
                    capture_label = capture_scenario.label if include_capture_in_names else None
                    cond_name = condition_dir_name(
                        gdna_label,
                        strand_spec,
                        nrna_label,
                        capture_label,
                        gdna_strand_overdispersion=gdna_od,
                    )
                    if selected_conditions and cond_name not in selected_conditions:
                        continue
                    cond_num += 1

                    depths = resolve_depths(sim, gdna_rate=gdna_rate)
                    n_rna, n_gdna = depths.n_rna, depths.n_gdna

                    condition_seed = capture_paired_condition_seed(
                        base_seed, gdna_label, strand_spec
                    )
                    cond_dir = outdir / cond_name
                    truth_abundances_name = f"{cond_name}/truth_abundances.tsv"
                    truth_fl_name = f"{cond_name}/truth_fragment_lengths.tsv"
                    truth_summary_name = f"{cond_name}/truth_summary.json"
                    probe_meta = capture_meta_by_label.get(capture_scenario.label, {})

                    print(
                        f"\n  [{cond_num}/{total}] {cond_name}: RNA={n_rna:,} gDNA={n_gdna:,} "
                        f"SS={strand_spec:.2f} nRNA={nrna_label} capture={capture_scenario.label}",
                        flush=True,
                    )

                    cond_entry: dict = {
                        "name": cond_name,
                        "gdna_label": gdna_label,
                        "gdna_rate": gdna_rate,
                        "gdna_strand_overdispersion": gdna_od,
                        "gdna_strand_overdispersion_label": gdna_od_label,
                        "strand_specificity": strand_spec,
                        "nrna_label": nrna_label,
                        "nrna_mode": nrna_mode,
                        "nrna_ratio": nrna_ratio,
                        "nrna_fragment_share": nrna_share,
                        "nrna_abundance_range": nrna_abundance_range,
                        "nrna_on_fraction": (nrna.on_fraction if nrna_mode == "sparse" else None),
                        "capture_label": capture_scenario.label,
                        "capture_enabled": bool(capture_scenario.config.probes),
                        "capture_config": capture_scenario.config,
                        "capture_probe_source": probe_meta.get("source"),
                        "capture_probe_panel": probe_meta.get("panel"),
                        "capture_probe_tsv": probe_meta.get("tsv"),
                        "capture_probe_bed": probe_meta.get("bed"),
                        "n_rna": n_rna,
                        "n_gdna": n_gdna,
                        "n_total": n_rna + n_gdna,
                        "seed": condition_seed,
                        "truth_kind": "post_capture_empirical",
                        "pre_capture_abundances": molecular_truth_name,
                        "post_capture_abundances": truth_abundances_name,
                        "post_capture_fragment_lengths": truth_fl_name,
                        "molecular_truth_abundances": molecular_truth_name,
                        "truth_abundances": truth_abundances_name,
                        "truth_fragment_lengths": truth_fl_name,
                        "truth_summary": truth_summary_name,
                        "fastq_r1": f"{cond_name}/sim_R1.fq.gz",
                        "fastq_r2": f"{cond_name}/sim_R2.fq.gz",
                    }

                    # The existence key is the oracle BAM whenever one is requested, because the BAM
                    # is the artifact the instruments read. Keying on ``sim_R1.fq.gz`` instead makes
                    # a panel whose FASTQs were dropped — by ``--no-fastq``, or by hand to reclaim
                    # disk — silently re-simulate every condition and rewrite the oracle.
                    # The trade, stated: on a panel whose BAM exists but whose FASTQs are gone, a run
                    # that wants FASTQs will skip instead of regenerating them. Delete the BAM, or
                    # pass ``--no-skip-existing``, to force it.
                    _exists_key = cond_dir / ("sim_oracle.bam" if oracle_bam else "sim_R1.fq.gz")
                    if skip_existing and _exists_key.exists():
                        print("    Output exists, skipping", flush=True)
                        cond_entry["oracle_bam"] = (
                            f"{cond_name}/sim_oracle.bam" if oracle_bam else None
                        )
                    else:
                        print("    Simulating...", end="", flush=True)
                        t0 = time.monotonic()
                        cond_sim = replace(sim, sim_seed=condition_seed)
                        simulator = WholeGenomeSimulator(
                            genome_path,
                            cond_transcripts,
                            cond_sim,
                            replace(gdna, strand_overdispersion=gdna_od),
                            genomic_refs=genomic_refs,
                            strand_specificity=strand_spec,
                            capture_config=capture_scenario.config,
                            capture_sampler=samplers[capture_scenario.label],
                        )
                        _, _, bam_path = simulator.simulate_and_write(
                            cond_dir,
                            n_rna,
                            n_gdna,
                            oracle_bam=oracle_bam,
                            prefix="sim",
                            n_workers=sim.n_workers,
                        )
                        simulator.close()
                        cond_entry["oracle_bam"] = (
                            f"{cond_name}/sim_oracle.bam" if bam_path else None
                        )
                        print(f" done ({time.monotonic() - t0:.1f}s)", flush=True)

                    bam_source = cond_dir / "sim_oracle.bam"
                    truth_summary = write_post_capture_truth(
                        cond_transcripts,
                        outdir / truth_abundances_name,
                        outdir / truth_fl_name,
                        outdir / truth_summary_name,
                        bam_path=bam_source if bam_source.exists() else None,
                        fastq_path=cond_dir / "sim_R1.fq.gz",
                        condition=cond_name,
                        molecular_truth=molecular_truth_name,
                        gdna_strand_overdispersion=gdna_od,
                    )
                    # The FASTQs are dropped only after the truth is written, never before, so an
                    # interrupted run cannot lose the origin counts. No calibration instrument reads
                    # one — nothing under ``scripts/design/`` or ``src/rigel/calibration/`` opens a
                    # FASTQ — and ``write_post_capture_truth`` prefers the oracle BAM and returns
                    # before touching the FASTQ path (``sim.truth._iter_origins_from_source``), so
                    # dropping them cannot change a truth file. They are roughly half a suite's
                    if not emit_fastq:
                        for _fq in (cond_dir / "sim_R1.fq.gz", cond_dir / "sim_R2.fq.gz"):
                            if _fq.exists():
                                _fq.unlink()

                    origin_counts = truth_summary["origin_counts"]
                    cond_entry["n_mrna_observed"] = int(origin_counts.get("mrna", 0))
                    cond_entry["n_nrna_observed"] = int(origin_counts.get("nrna", 0))
                    cond_entry["n_gdna_observed"] = int(origin_counts.get("gdna", 0))
                    conditions.append(cond_entry)

    return conditions
