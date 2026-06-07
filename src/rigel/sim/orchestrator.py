"""The shared condition-grid orchestrator for the simulator suites.

Both ``suite.main`` (synthetic mini-genome suites) and ``whole_genome.run_simulation`` (simulate
from an existing reference) sweep the same condition grid — nRNA × gDNA-rate × gDNA-strand-
overdispersion × strand-specificity × capture — and for each condition deep-copy the transcripts,
apply the nRNA spike-in, allocate fragment counts, run the :class:`whole_genome.WholeGenomeSimulator`,
write the post-capture truth, and record a manifest entry. That loop lived in two near-identical
copies (so every new axis had to be wired twice). This is the single implementation both call.

Seeding (intentional, see docs/sim/architecture_redesign.md P4): every condition gets a *distinct*
per-condition seed via :func:`capture_paired_condition_seed`, keyed by ``(gdna, ss, nrna)`` so that
capture- and overdispersion-variants of one base condition are paired (share a seed) for controlled
comparison while the main axes are decorrelated. This is the scheme ``suite.main`` already used; it
also fixes ``run_simulation``'s former same-seed-for-all-conditions behaviour.
"""

from __future__ import annotations

import copy
import hashlib
import time
from dataclasses import replace
from itertools import product
from pathlib import Path

from .manifest import condition_dir_name
from .truth import write_post_capture_truth
from .whole_genome import (
    WholeGenomeSimulator,
    apply_nrna_ratio,
    apply_random_nrna_fraction,
    write_truth_abundances,
)

__all__ = ["stable_seed", "capture_paired_condition_seed", "run_condition_grid"]


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
    """Seed shared by the capture/overdispersion variants of one base ``(gdna, ss, nrna)`` condition."""
    seed_name = condition_dir_name(gdna_label, strand_specificity, nrna_label)
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
    gdna_pairs: list[tuple[str, float]],
    gdna_od_pairs: list[tuple[str, float]],
    strand_specificities: list[float],
    nrna_pairs: list[tuple[str, str, object, int]],
    capture_scenarios: list,
    include_capture_in_names: bool,
    base_seed: int,
    oracle_bam: bool = True,
    skip_existing: bool = True,
    selected_conditions: set[str] | None = None,
    capture_meta_by_label: dict[str, dict] | None = None,
) -> list[dict]:
    """Run the full condition grid and return the per-condition manifest entries.

    ``nrna_pairs`` entries are ``(label, mode, value, index)`` (see ``whole_genome._build_nrna_pairs``):
    mode ``additive_ratio`` / ``random_fraction`` spike nRNA onto the base abundances and allocate
    explicit mRNA/nRNA counts; any other mode (``file``) leaves abundances as loaded and allocates a
    single RNA pool. ``capture_meta_by_label`` supplies the suite's probe-provenance fields per
    capture label (empty for the reference-driven path). The caller writes the manifest.
    """
    capture_meta_by_label = capture_meta_by_label or {}
    conditions: list[dict] = []
    cond_num = 0
    total = (
        len(nrna_pairs) * len(gdna_pairs) * len(gdna_od_pairs)
        * len(strand_specificities) * len(capture_scenarios)
    )

    for nrna_label, nrna_mode, nrna_value, nrna_index in nrna_pairs:
        cond_transcripts = copy.deepcopy(transcripts)
        for t, (base_mrna, base_nrna) in zip(cond_transcripts, base_abundances):
            t.abundance = base_mrna
            t.nrna_abundance = base_nrna

        nrna_ratio: float | None = None
        nrna_ratio_range: tuple[float, float] | None = None
        if nrna_mode == "additive_ratio":
            nrna_ratio = float(nrna_value or 0.0)
            apply_nrna_ratio(cond_transcripts, nrna_ratio)
        elif nrna_mode == "random_fraction":
            nrna_ratio_range = tuple(nrna_value)  # type: ignore[arg-type]
            nrna_ratio = apply_random_nrna_fraction(
                cond_transcripts, nrna_ratio_range,
                eligible_fraction=nrna.eligible_fraction, seed=nrna.seed + nrna_index,
            )

        molecular_truth_name = f"truth_abundances_nrna_{nrna_label}.tsv"
        write_truth_abundances(cond_transcripts, outdir / molecular_truth_name)

        for gdna_label, gdna_rate in gdna_pairs:
            for (gdna_od_label, gdna_od), strand_spec in product(
                gdna_od_pairs, strand_specificities
            ):
                for capture_scenario in capture_scenarios:
                    capture_label = (
                        capture_scenario.label if include_capture_in_names else None
                    )
                    cond_name = condition_dir_name(
                        gdna_label, strand_spec, nrna_label, capture_label,
                        gdna_strand_overdispersion=gdna_od,
                    )
                    if selected_conditions and cond_name not in selected_conditions:
                        continue
                    cond_num += 1

                    if nrna_mode in {"additive_ratio", "random_fraction"}:
                        n_mrna: int | None = sim.n_rna_fragments
                        n_nrna: int | None = round(n_mrna * float(nrna_ratio or 0.0))
                        n_rna = n_mrna + n_nrna
                        n_gdna = round(gdna_rate * n_mrna)
                    else:  # file mode: a single RNA pool, no explicit mRNA/nRNA split
                        n_mrna = n_nrna = None
                        n_rna = sim.n_rna_fragments
                        n_gdna = round(gdna_rate * n_rna)

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
                        "nrna_ratio_range": nrna_ratio_range,
                        "nrna_eligible_fraction": (
                            nrna.eligible_fraction if nrna_mode == "random_fraction" else None
                        ),
                        "capture_label": capture_scenario.label,
                        "capture_enabled": bool(capture_scenario.config.probes),
                        "capture_config": capture_scenario.config,
                        "capture_probe_source": probe_meta.get("source"),
                        "capture_probe_panel": probe_meta.get("panel"),
                        "capture_probe_tsv": probe_meta.get("tsv"),
                        "capture_probe_bed": probe_meta.get("bed"),
                        "n_mrna": n_mrna,
                        "n_nrna": n_nrna,
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

                    if skip_existing and (cond_dir / "sim_R1.fq.gz").exists():
                        print("    Output exists, skipping", flush=True)
                        cond_entry["oracle_bam"] = (
                            f"{cond_name}/sim_oracle.bam" if oracle_bam else None
                        )
                    else:
                        print("    Simulating...", end="", flush=True)
                        t0 = time.monotonic()
                        cond_sim = replace(sim, sim_seed=condition_seed)
                        simulator = WholeGenomeSimulator(
                            genome_path, cond_transcripts, cond_sim,
                            replace(gdna, strand_overdispersion=gdna_od),
                            strand_specificity=strand_spec,
                            capture_config=capture_scenario.config,
                        )
                        _, _, bam_path = simulator.simulate_and_write(
                            cond_dir, n_rna, n_gdna,
                            oracle_bam=oracle_bam, prefix="sim",
                            n_mrna=n_mrna, n_nrna=n_nrna, n_workers=sim.n_workers,
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
                    origin_counts = truth_summary["origin_counts"]
                    cond_entry["n_mrna_observed"] = int(origin_counts.get("mrna", 0))
                    cond_entry["n_nrna_observed"] = int(origin_counts.get("nrna", 0))
                    cond_entry["n_gdna_observed"] = int(origin_counts.get("gdna", 0))
                    conditions.append(cond_entry)

    return conditions
