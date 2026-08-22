"""The shared condition-grid orchestrator for the simulator suites.

Both ``suite.main`` (synthetic mini-genome suites) and ``whole_genome.run_simulation`` (simulate
from an existing reference) sweep the same condition grid — nRNA × gDNA-rate × gDNA-strand-
overdispersion × strand-specificity × capture — and for each condition deep-copy the transcripts,
apply the nRNA mode, allocate fragment counts, run the :class:`whole_genome.WholeGenomeSimulator`,
write the post-capture truth, and record a manifest entry. That loop lived in two near-identical
copies (so every new axis had to be wired twice). This is the single implementation both call.

Seeding (intentional): every condition gets a *distinct*
per-condition seed via :func:`capture_paired_condition_seed`, keyed by ``(gdna, ss, nrna)`` so that
capture- and overdispersion-variants of one base condition are paired (share a seed) for controlled
comparison while the main axes are decorrelated. This is the scheme ``suite.main`` already used; it
also fixes ``run_simulation``'s former same-seed-for-all-conditions behaviour.
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
    """How one condition's fragment budget splits between RNA and gDNA. ⭐ The mature/nascent split
    INSIDE ``n_rna`` is never imposed (owner, 2026-08-19): nascent entities are transcripts in the same
    multinomial, so their share follows from molecules and lengths and is read off the realised
    origin counts afterwards."""

    n_rna: int
    n_gdna: int

    @property
    def total(self) -> int:
        return self.n_rna + self.n_gdna


def resolve_depths(sim, *, gdna_rate: float) -> ConditionDepths:
    """Split one condition's fragment budget between RNA and gDNA.

    ⭐ **TWO MODES, AND THE SECOND ONE EXISTS BECAUSE THE FIRST CANNOT REACH THE HIGH-gDNA END.**

    * ``n_total_fragments is None`` (the default, and every pre-existing config) — the legacy
      behaviour, unchanged: the RNA depth is fixed and gDNA is added ON TOP at ``rate x n_rna``. A
      ``rate`` of 1.0 is therefore a 50 % gDNA library at *twice* the depth of a ``rate`` of 0.
    * ``n_total_fragments`` set — the TOTAL is fixed and ``rate`` decides only the SPLIT:
      ``n_rna = total/(1+rate)``, ``n_gdna = total − n_rna``.

    ⛔ The second mode is not a convenience. Under the first, a 98 % gDNA library is ``rate = 49``,
    i.e. 490 M fragments against a 10 M RNA depth — unsimulatable. It is also the more faithful model:
    a sequencing run has a fixed budget and the contamination fraction decides how it is *split*, not
    how much extra is generated. ⚠ The accepted trade-off (owner, 2026-08-03) is that the RNA side
    thins as gDNA rises, so per-transcript abundance accuracy degrades at the top of the ladder — that
    is a property of such libraries, not an artifact.

    ⚠ **Nascent comes OUT of the RNA share, never on top** — and since 2026-08-19 that is structural
    rather than arithmetic: the nascent entities are rows of the ONE RNA multinomial, so ``n_rna`` is
    the whole RNA budget and the mature/nascent split inside it is realised, not imposed.

    Gates: ``tests/test_sim_fixed_total_depth.py``.
    """
    total = getattr(sim, "n_total_fragments", None)

    if total is None:
        # ── legacy: RNA depth fixed, gDNA added on top ────────────────────────────────────────
        n_rna = int(sim.n_rna_fragments)
        return ConditionDepths(n_rna, round(gdna_rate * n_rna))

    # ── fixed total: `rate` decides the split only ────────────────────────────────────────────
    total = int(total)
    n_rna = round(total / (1.0 + float(gdna_rate)))
    # ⭐ gDNA is the REMAINDER, so the total is conserved by construction rather than by luck.
    # ⚠ Two plausible-sounding claims about WHY were both written here and both measured FALSE:
    # ``round(T−x)`` and ``round(rate·n_rna)`` each conserve the total at every rung of the panel
    # ladder, so neither is what this guards against. What actually breaks conservation is computing
    # the two shares INDEPENDENTLY WITH TRUNCATION (``int()`` on both), which drifts −1 at four of the
    # nine rungs. That is the perturbation the gate uses, and it is the only one of the three that
    # fires — a reminder that "obviously safer" arithmetic deserves a measurement like anything else.
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
    nrna_label: str,
) -> int:
    """Seed shared by the capture **and nascent** variants of one ``(gdna, ss)`` base condition.

    ``nrna_label`` is intentionally **excluded** from the seed so the variants of one base condition
    start from the same stream. ⚠ Since 2026-08-19 nascent rows are drawn in the SAME multinomial as
    the mature rows, so a nascent-on cell and its nascent-off twin share the seed but NOT a
    bit-identical mature stream — turning nascent on re-allocates the RNA budget, as it physically
    must; the gDNA stream is unaffected.
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
    ⭐ The three live modes reach the entities by two different routes and the difference is the point:
    ``additive_ratio`` and ``fragment_share`` POOL each entity's molecules from its contributors
    (`whole_genome.assign_nrna_to_entities`, so nascent tracks mature and cannot exceed it), while
    ``sparse`` writes each ENTITY's absolute abundance directly (`whole_genome.apply_sparse_nrna`, so
    most entities get exactly zero and the rest are independent of the mature level). ``file`` leaves
    the loaded abundances alone. Every mode then shares one thing: the RNA budget is ONE multinomial
    over mature and nascent rows, so the fragment split is realised rather than allocated.
    ``capture_meta_by_label`` supplies the suite's probe-provenance fields per
    capture label (empty for the reference-driven path). The caller writes the manifest.
    """
    capture_meta_by_label = capture_meta_by_label or {}
    # ⭐⭐ ONE CaptureSampler PER CAPTURE SCENARIO, built once and reused by every condition that shares
    # it. The probe layout and the per-(width) partition depend only on the panel and the templates —
    # not on abundance, gDNA rate, strand or nascent — so rebuilding per condition recomputed 260 s of
    # identical numbers each time (measured 2026-08-19).
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
            # ⭐ the config states the nascent share of RNA FRAGMENTS; the MOLECULAR ratio is solved
            # from the annotation (`whole_genome.apply_nrna_fragment_share`) and recorded per condition
            nrna_share = float(nrna_value or 0.0)
            nrna_ratio = apply_nrna_fragment_share(cond_transcripts, nrna_share, sim)
        elif nrna_mode == "sparse":
            # ⭐ nascent is ABSENT from most gene spans and INDEPENDENT of the mature level where it
            # is present (owner, 2026-08-22); the fragment share is emergent and recorded below
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
                        base_seed, gdna_label, strand_spec, nrna_label
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
                        "nrna_on_fraction": (
                            nrna.on_fraction if nrna_mode == "sparse" else None
                        ),
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

                    # ⭐⭐ THE EXISTENCE KEY IS THE ORACLE BAM WHENEVER ONE IS REQUESTED, and that is a
                    # correctness fix rather than a convenience. It used to be ``sim_R1.fq.gz``, so a
                    # panel whose FASTQs had been dropped — deliberately by ``--no-fastq``, or by hand
                    # to reclaim disk — would silently RE-SIMULATE every condition on the next run,
                    # which for the 36-condition ladder is hours and a rewritten oracle. The BAM is the
                    # artifact every instrument actually reads, so it is the honest key.
                    # ⚠ The trade, stated: on a panel whose BAM exists but whose FASTQs are gone, a run
                    # that WANTS FASTQs will skip instead of regenerating them. Delete the BAM, or pass
                    # ``--no-skip-existing``, to force it.
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
                    # ⭐⭐ **THE FASTQs ARE DROPPED AFTER THE TRUTH IS WRITTEN** (owner, 2026-08-07).
                    # No calibration instrument reads one — nothing under ``scripts/design/`` or
                    # ``src/rigel/calibration/`` opens a FASTQ — and ``write_post_capture_truth``
                    # prefers the oracle BAM and returns before touching the FASTQ path
                    # (``sim.truth._iter_origins_from_source``), so this cannot change a truth file.
                    # ⚠ They are 30 G of a 67 G suite and ARE read by ``scripts/benchmarking/`` for the
                    # third-party tool comparison; a panel built with this off cannot be benchmarked
                    # against another tool without re-simulating. ⛔ Deleted only AFTER the truth is
                    # written, never before, so an interrupted run cannot lose the origin counts.
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
