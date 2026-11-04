"""The condition-grid loop of ``whole_genome.run_simulation``: nascent x gDNA rate x gDNA strand overdispersion x
strand specificity x capture.

Per condition, :func:`run_condition_grid` resets the transcripts to the base abundances, applies the nascent mode,
splits the fragment budget (:func:`resolve_depths`), simulates on the capture sampler its panel shares, writes the
truth and returns a manifest entry. A condition whose output already exists is skipped. The seed keys on
``(gdna, strand specificity)`` only (:func:`capture_paired_condition_seed`), so the capture, overdispersion and
nascent variants of one base condition are paired.
"""

from __future__ import annotations

import copy
import hashlib
import time
from dataclasses import dataclass, replace
from itertools import product
from pathlib import Path

from .manifest import condition_dir_name
from .truth import write_molar_truth, write_post_capture_truth
from .capture import CaptureSampler
from .wgs_engine import WholeGenomeSimulator
from .whole_genome import apply_nrna_ratio, apply_sparse_nrna

__all__ = [
    "stable_seed",
    "capture_paired_condition_seed",
    "run_condition_grid",
    "ConditionDepths",
    "resolve_depths",
]


@dataclass(frozen=True, slots=True)
class ConditionDepths:
    """How one condition's fragment budget splits between RNA and gDNA. The nascent share of ``n_rna`` is not
    imposed: nascent entities are rows of the one RNA multinomial, so it follows from molecules and lengths."""

    n_rna: int
    n_gdna: int

    @property
    def total(self) -> int:
        return self.n_rna + self.n_gdna


def resolve_depths(sim, *, gdna_rate: float) -> ConditionDepths:
    """Split one condition's fragment budget between RNA and gDNA:

    * ``n_total_fragments is None`` (the default): the RNA depth is fixed and gDNA is added at ``rate x n_rna``,
      so a ``rate`` of 1.0 is a 50 % gDNA library at twice the depth of ``rate`` 0.
    * ``n_total_fragments`` set: the total is fixed and ``rate`` decides only the split,
      ``n_rna = total/(1+rate)``, ``n_gdna = total − n_rna``. This mirrors a sequencing run's fixed budget and
      reaches high gDNA fractions at a bounded depth; the RNA side thins as gDNA rises.

    Nascent comes out of ``n_rna``, never on top of it. Raises when the split leaves no RNA.

    Gates: ``tests/test_sim_fixed_total_depth.py``.
    """
    total = getattr(sim, "n_total_fragments", None)

    if total is None:
        # additive: RNA depth fixed, gDNA added on top
        n_rna = int(sim.n_rna_fragments)
        return ConditionDepths(n_rna, round(gdna_rate * n_rna))

    # fixed total: `rate` decides the split only
    total = int(total)
    n_rna = round(total / (1.0 + float(gdna_rate)))
    # gDNA must stay the remainder: two shares computed independently drift from the total by a fragment.
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

    A nascent-on cell and its nascent-off twin share the seed but not an identical mature stream, since nascent
    rows are drawn in the same multinomial and re-allocate the RNA budget.
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
    emit_fastq: bool = True,
) -> list[dict]:
    """Run the full condition grid and return the per-condition manifest entries for the caller to write.

    ``nrna_pairs`` entries are ``(label, mode, value, index)`` (``whole_genome._build_nrna_pairs``).
    ``additive_ratio`` pools each entity's molecules from its contributors (`whole_genome.apply_nrna_ratio`),
    ``sparse`` sets each entity's absolute abundance (`whole_genome.apply_sparse_nrna`), and ``file`` keeps the
    loaded abundances.
    """
    # One CaptureSampler per capture scenario, shared by its conditions: the probe layout and the per-width
    # partition depend only on the panel and the templates.
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
        if nrna_mode == "additive_ratio":
            nrna_ratio = float(nrna_value or 0.0)
            apply_nrna_ratio(cond_transcripts, nrna_ratio)
        elif nrna_mode == "sparse":
            nrna_abundance_range = tuple(nrna_value)  # type: ignore[arg-type]
            nrna_ratio = apply_sparse_nrna(
                cond_transcripts,
                nrna_abundance_range,
                on_fraction=nrna.on_fraction,
                seed=nrna.seed + nrna_index,
            )

        molecular_truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        write_molar_truth(cond_transcripts, outdir / molecular_truth_name)

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
                        "nrna_abundance_range": nrna_abundance_range,
                        "nrna_on_fraction": (nrna.on_fraction if nrna_mode == "sparse" else None),
                        "capture_label": capture_scenario.label,
                        "capture_enabled": bool(capture_scenario.config.probes),
                        "capture_config": capture_scenario.config,
                        "n_rna": n_rna,
                        "n_gdna": n_gdna,
                        "n_total": n_rna + n_gdna,
                        "seed": condition_seed,
                        # each table's key is the one its truth_summary.json uses
                        "truth_kind": "post_capture_empirical",
                        "pre_capture_abundances": molecular_truth_name,
                        "post_capture_abundances": truth_abundances_name,
                        "post_capture_fragment_lengths": truth_fl_name,
                        "truth_summary": truth_summary_name,
                        "fastq_r1": f"{cond_name}/sim_R1.fq.gz",
                        "fastq_r2": f"{cond_name}/sim_R2.fq.gz",
                    }

                    # The resume key is the oracle BAM whenever one is requested, so a condition whose BAM
                    # exists is skipped even if its FASTQs are gone; delete the BAM to regenerate it.
                    _exists_key = cond_dir / ("sim_oracle.bam" if oracle_bam else "sim_R1.fq.gz")
                    if _exists_key.exists():
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
                        pre_capture_abundances=molecular_truth_name,
                        gdna_strand_overdispersion=gdna_od,
                    )
                    # The FASTQs are dropped only after the truth is written, since it falls back to them.
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
