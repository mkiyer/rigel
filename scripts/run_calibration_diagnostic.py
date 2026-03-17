#!/usr/bin/env python3
"""Run a single simulation scenario and produce calibration diagnostics.

Uses the Rigel Scenario API to create a synthetic BAM + index with
known ground truth, then runs the calibration diagnostics pipeline
to produce a self-contained HTML report.

Usage::

    # Stranded sample with 30% gDNA contamination
    python scripts/run_calibration_diagnostic.py \\
        -c scripts/calibration_validation_small.yaml \\
        -o diag_output/ --ss 1.0 --gdna-fraction 0.3 --kappa 50

    # Use config defaults (first combination from the sweep grid)
    python scripts/run_calibration_diagnostic.py \\
        -c scripts/calibration_validation_small.yaml -o diag_output/

    # Multiple scenarios in one go
    python scripts/run_calibration_diagnostic.py \\
        -c scripts/calibration_validation_small.yaml -o diag_output/ \\
        --scenarios stranded_clean stranded_dirty unstranded_dirty
"""
from __future__ import annotations

import argparse
import logging
import sys
import textwrap
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Ensure src/ is importable when run from scripts/
# ---------------------------------------------------------------------------
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))
# Also make scripts/ importable (for calibration_diagnostics)
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from rigel.sim import GDNAConfig, Scenario, SimConfig

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Re-use nRNA group parsing from the sweep infrastructure
# ---------------------------------------------------------------------------
from synthetic_sim_sweep import (
    NrnaGroup,
    compute_nrna_groups,
    parse_nrna_coords,
)


# ---------------------------------------------------------------------------
# Predefined diagnostic scenarios
# ---------------------------------------------------------------------------

SCENARIOS = {
    "stranded_clean": {
        "label": "Stranded, no gDNA",
        "strand_specificity": 1.0,
        "gdna_fraction": 0.0,
        "gdna_strand_kappa": 50,
        "n_rna_fragments": 20000,
    },
    "stranded_low_gdna": {
        "label": "Stranded, 10% gDNA",
        "strand_specificity": 1.0,
        "gdna_fraction": 0.1,
        "gdna_strand_kappa": 50,
        "n_rna_fragments": 20000,
    },
    "stranded_dirty": {
        "label": "Stranded, 30% gDNA",
        "strand_specificity": 1.0,
        "gdna_fraction": 0.3,
        "gdna_strand_kappa": 50,
        "n_rna_fragments": 20000,
    },
    "stranded_very_dirty": {
        "label": "Stranded, 50% gDNA",
        "strand_specificity": 1.0,
        "gdna_fraction": 0.5,
        "gdna_strand_kappa": 50,
        "n_rna_fragments": 20000,
    },
    "partial_strand_dirty": {
        "label": "Partial strand (SS=0.75), 30% gDNA",
        "strand_specificity": 0.75,
        "gdna_fraction": 0.3,
        "gdna_strand_kappa": 50,
        "n_rna_fragments": 20000,
    },
    "unstranded_dirty": {
        "label": "Unstranded (SS=0.5), 30% gDNA",
        "strand_specificity": 0.5,
        "gdna_fraction": 0.3,
        "gdna_strand_kappa": 1000,
        "n_rna_fragments": 20000,
    },
    "high_overdispersion": {
        "label": "Stranded, 30% gDNA, high strand overdispersion (κ=5)",
        "strand_specificity": 1.0,
        "gdna_fraction": 0.3,
        "gdna_strand_kappa": 5,
        "n_rna_fragments": 20000,
    },
}


def _build_scenario_from_config(
    config: dict,
    *,
    strand_specificity: float = 1.0,
    gdna_fraction: float = 0.3,
    gdna_strand_kappa: float = 50.0,
    n_rna_fragments: int = 20000,
    work_dir: Path,
) -> tuple:
    """Build a Scenario + configs from YAML config dict.

    Returns (scenario, sim_config, gdna_config, n_rna_fragments,
             gdna_fraction, run_params_dict).
    """
    genome_length = config.get("genome_length", 500000)
    seed = config.get("seed", 42)
    read_length = config.get("read_length", 150)
    rna_params = config.get("rna", {})
    gdna_params = config.get("gdna", {})
    transcripts_config = config["transcripts"]

    # Parse nRNA groups
    nrnas_config = config.get("nrnas")
    if nrnas_config:
        nrna_groups = parse_nrna_coords(nrnas_config, transcripts_config)
    else:
        user_groups = config.get("nrna_groups")
        nrna_groups = compute_nrna_groups(transcripts_config, user_groups)

    # Get pattern (first pattern entry, or default abundances)
    patterns = config.get("patterns", [])
    sweep = config.get("sweep", {})
    pattern = patterns[0] if patterns else {}

    # Build Scenario
    sc = Scenario(
        "diag",
        genome_length=genome_length,
        seed=seed,
        work_dir=work_dir,
    )

    # Add genes grouped by nRNA group
    for grp_label, grp in nrna_groups.items():
        nrna_ab = float(pattern.get(grp_label, sweep.get(grp_label, 0)))
        t_list = []
        for t_id in grp.t_ids:
            t_def = transcripts_config[t_id]
            ab = float(pattern.get(t_id, sweep.get(t_id, 0)))
            nrna = nrna_ab if t_id == grp.carrier else 0.0
            t_list.append({
                "t_id": t_id,
                "exons": [tuple(e) for e in t_def["exons"]],
                "abundance": ab,
                "nrna_abundance": nrna,
            })
        sc.add_gene(grp_label, grp.strand, t_list)

    # Configs
    sim_cfg = SimConfig(
        frag_mean=rna_params.get("frag_mean", 250),
        frag_std=rna_params.get("frag_std", 50),
        frag_min=rna_params.get("frag_min", 80),
        frag_max=rna_params.get("frag_max", 600),
        read_length=read_length,
        strand_specificity=strand_specificity,
        seed=seed,
    )

    gdna_cfg = None
    if gdna_fraction > 0:
        gdna_cfg = GDNAConfig(
            abundance=1.0,
            frag_mean=gdna_params.get("frag_mean", 350),
            frag_std=gdna_params.get("frag_std", 80),
            frag_min=gdna_params.get("frag_min", 100),
            frag_max=gdna_params.get("frag_max", 1000),
            strand_kappa=gdna_strand_kappa if gdna_strand_kappa > 0 else None,
        )

    run_params = {
        "genome_length": genome_length,
        "seed": seed,
        "read_length": read_length,
        "strand_specificity": f"{strand_specificity:.4f}",
        "n_rna_fragments": n_rna_fragments,
        "gdna_fraction": f"{gdna_fraction:.2f}",
        "gdna_strand_kappa": gdna_strand_kappa,
        "rna_frag_mean": rna_params.get("frag_mean", 250),
        "rna_frag_std": rna_params.get("frag_std", 50),
        "n_transcripts": len(transcripts_config),
        "n_expressed": sum(
            1 for t_id in transcripts_config
            if float(pattern.get(t_id, sweep.get(t_id, 0))) > 0
        ),
        "n_nrna_groups": len(nrna_groups),
    }
    if gdna_cfg:
        run_params["gdna_frag_mean"] = gdna_cfg.frag_mean
        run_params["gdna_frag_std"] = gdna_cfg.frag_std

    return sc, sim_cfg, gdna_cfg, n_rna_fragments, gdna_fraction, run_params


def run_single_diagnostic(
    config: dict,
    output_dir: Path,
    *,
    scenario_params: dict | None = None,
    scenario_name: str = "default",
) -> Path:
    """Simulate + run calibration diagnostics for one scenario.

    Parameters
    ----------
    config : dict
        Parsed YAML config.
    output_dir : Path
        Output directory for this scenario.
    scenario_params : dict or None
        Override params (strand_specificity, gdna_fraction, etc).
    scenario_name : str
        Name for labels/logging.

    Returns
    -------
    Path to generated HTML report.
    """
    from calibration_diagnostics import run_diagnostics

    if scenario_params is None:
        scenario_params = {}

    work_dir = output_dir / "sim_artifacts"

    ss = float(scenario_params.get("strand_specificity", 1.0))
    gdna_frac = float(scenario_params.get("gdna_fraction", 0.3))
    kappa = float(scenario_params.get("gdna_strand_kappa", 50))
    n_rna = int(scenario_params.get("n_rna_fragments", 20000))

    logger.info("=" * 60)
    logger.info(f"Scenario: {scenario_name}")
    logger.info(f"  SS={ss}, gDNA frac={gdna_frac}, κ={kappa}, n_rna={n_rna}")
    logger.info("=" * 60)

    # Build scenario from config
    sc, sim_cfg, gdna_cfg, n_rna_frags, gfrac, run_params = \
        _build_scenario_from_config(
            config,
            strand_specificity=ss,
            gdna_fraction=gdna_frac,
            gdna_strand_kappa=kappa,
            n_rna_fragments=n_rna,
            work_dir=work_dir,
        )

    # Add scenario label to run params
    run_params["scenario"] = scenario_name
    label = scenario_params.get("label", scenario_name)
    run_params["description"] = label

    # Simulate
    logger.info("Running simulation...")
    result = sc.build_oracle(
        n_fragments=n_rna_frags,
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        nrna_abundance=0.0,
        n_rna_fragments=n_rna_frags,
        gdna_fraction=gfrac if gfrac > 0 else None,
    )
    logger.info(f"BAM: {result.bam_path}")
    logger.info(f"Index: {result.index_dir}")

    # Run diagnostics
    report_path = run_diagnostics(
        bam_path=str(result.bam_path),
        index_dir=str(result.index_dir),
        output_dir=str(output_dir),
        strand_specificity=ss,
        max_iterations=50,
        convergence_tol=1e-4,
        parse_truth=True,
        emit_html=True,
        run_params=run_params,
    )

    html_path = output_dir / "calibration_report.html"
    logger.info(f"Report: {html_path}")
    return html_path


def main():
    parser = argparse.ArgumentParser(
        description="Run synthetic simulation + calibration diagnostics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Predefined scenarios:
              stranded_clean        - SS=1.0, 0%% gDNA
              stranded_low_gdna     - SS=1.0, 10%% gDNA
              stranded_dirty        - SS=1.0, 30%% gDNA
              stranded_very_dirty   - SS=1.0, 50%% gDNA
              partial_strand_dirty  - SS=0.75, 30%% gDNA
              unstranded_dirty      - SS=0.5, 30%% gDNA
              high_overdispersion   - SS=1.0, 30%% gDNA, κ=5

            Examples:
              # Run a single predefined scenario
              python scripts/run_calibration_diagnostic.py \\
                  -c scripts/calibration_validation_small.yaml \\
                  -o diag_out/ --scenarios stranded_dirty

              # Run all predefined scenarios
              python scripts/run_calibration_diagnostic.py \\
                  -c scripts/calibration_validation_small.yaml \\
                  -o diag_out/ --all

              # Custom parameters
              python scripts/run_calibration_diagnostic.py \\
                  -c scripts/calibration_validation_small.yaml \\
                  -o diag_out/ --ss 0.9 --gdna-fraction 0.2 --kappa 20
        """),
    )
    parser.add_argument("-c", "--config", required=True,
                        help="YAML config file (same format as sweep configs)")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory")
    parser.add_argument("--scenarios", nargs="+", choices=list(SCENARIOS.keys()),
                        help="Run predefined scenario(s)")
    parser.add_argument("--all", action="store_true",
                        help="Run all predefined scenarios")
    parser.add_argument("--ss", type=float, default=None,
                        help="Strand specificity (for custom scenario)")
    parser.add_argument("--gdna-fraction", type=float, default=None,
                        help="gDNA fraction (for custom scenario)")
    parser.add_argument("--kappa", type=float, default=None,
                        help="gDNA strand overdispersion κ (for custom scenario)")
    parser.add_argument("--n-rna", type=int, default=None,
                        help="Number of RNA fragments (for custom scenario)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    out = Path(args.output_dir)

    # Determine which scenarios to run
    if args.all:
        scenario_list = list(SCENARIOS.keys())
    elif args.scenarios:
        scenario_list = args.scenarios
    elif any(v is not None for v in [args.ss, args.gdna_fraction, args.kappa, args.n_rna]):
        # Custom scenario from CLI params
        custom = {
            "label": "Custom scenario",
            "strand_specificity": args.ss or 1.0,
            "gdna_fraction": args.gdna_fraction if args.gdna_fraction is not None else 0.3,
            "gdna_strand_kappa": args.kappa or 50,
            "n_rna_fragments": args.n_rna or 20000,
        }
        scenario_list = None  # signal custom mode
    else:
        # Default: run a representative subset
        scenario_list = [
            "stranded_clean", "stranded_dirty", "partial_strand_dirty",
        ]

    html_paths = []

    if scenario_list is not None:
        for name in scenario_list:
            params = SCENARIOS[name]
            scenario_out = out / name
            scenario_out.mkdir(parents=True, exist_ok=True)
            html = run_single_diagnostic(
                config, scenario_out,
                scenario_params=params,
                scenario_name=name,
            )
            html_paths.append(html)
    else:
        # Custom scenario
        custom_out = out / "custom"
        custom_out.mkdir(parents=True, exist_ok=True)
        html = run_single_diagnostic(
            config, custom_out,
            scenario_params=custom,
            scenario_name="custom",
        )
        html_paths.append(html)

    # Print summary
    print("\n" + "=" * 60)
    print("CALIBRATION DIAGNOSTIC REPORTS")
    print("=" * 60)
    for p in html_paths:
        print(f"  {p}")
    print()
    if html_paths:
        print(f"Open in browser:  open {html_paths[0]}")


if __name__ == "__main__":
    main()
