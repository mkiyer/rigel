#!/usr/bin/env python3
"""Ablation study: fragment length distribution recovery across
strand specificity × RNA FL dist × gDNA FL dist.

Generates YAML configs for each (RNA, gDNA) distribution pair,
runs the sweep, and combines all results into a single TSV.
"""
import csv
import itertools
import subprocess
import sys
import textwrap
from pathlib import Path


# ── Parameter space ──────────────────────────────────────────────────
RNA_DISTS = [
    {"label": "rna_short", "frag_mean": 150, "frag_std": 20, "frag_min": 80, "frag_max": 400},
    {"label": "rna_typical", "frag_mean": 200, "frag_std": 30, "frag_min": 80, "frag_max": 600},
    {"label": "rna_long", "frag_mean": 300, "frag_std": 50, "frag_min": 100, "frag_max": 800},
]

GDNA_DISTS = [
    {"label": "gdna_short", "frag_mean": 250, "frag_std": 50, "frag_min": 100, "frag_max": 600},
    {"label": "gdna_typical", "frag_mean": 350, "frag_std": 80, "frag_min": 100, "frag_max": 1000},
    {"label": "gdna_long", "frag_mean": 500, "frag_std": 120, "frag_min": 100, "frag_max": 1200},
]

STRAND_SPECIFICITIES = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
GDNA_FRACTIONS = [0.1, 0.3, 0.5]

# Per-config: 7 SS × 3 gdna_frac = 21 runs
# 9 FL combos × 21 = 189 total runs

YAML_TEMPLATE = textwrap.dedent("""\
# Auto-generated ablation config: {rna_label} × {gdna_label}
genome_length: 50000
seed: 42
read_length: 150

rna:
  frag_mean: {rna_mean}
  frag_std: {rna_std}
  frag_min: {rna_min}
  frag_max: {rna_max}

gdna:
  frag_mean: {gdna_mean}
  frag_std: {gdna_std}
  frag_min: {gdna_min}
  frag_max: {gdna_max}

transcripts:
  TA1:
    strand: "+"
    exons: [[1000, 1020], [5000, 5500], [12000, 13000]]
  TA2:
    strand: "+"
    exons: [[2000, 2500], [5000, 5500], [12000, 13000]]
  TA4:
    strand: "+"
    exons: [[4500, 5500], [12000, 13000]]
  TB1:
    strand: "-"
    exons: [[9000, 10000], [14000, 16000]]
  TCF:
    strand: "+"
    exons: [[20000, 24000]]
  TCR:
    strand: "-"
    exons: [[23000, 26000], [30000, 31000], [35000, 36000]]

nrnas:
  NTA1: ["+", 1000, 13000]

sweep:
  strand_specificity: [{ss_list}]
  n_rna_fragments: 10000
  gdna_fraction: [{gf_list}]
  NTA1: 0
  TD: 0

patterns:
  - {{TA1: 100, TA2: 50, TA4: 50, TB1: 50, TCF: 0, TCR: 0}}
""")


def generate_yaml(rna, gdna, outdir):
    """Generate a sweep YAML for one (RNA, gDNA) distribution pair."""
    name = f"{rna['label']}__{gdna['label']}"
    content = YAML_TEMPLATE.format(
        rna_label=rna["label"],
        gdna_label=gdna["label"],
        rna_mean=rna["frag_mean"],
        rna_std=rna["frag_std"],
        rna_min=rna["frag_min"],
        rna_max=rna["frag_max"],
        gdna_mean=gdna["frag_mean"],
        gdna_std=gdna["frag_std"],
        gdna_min=gdna["frag_min"],
        gdna_max=gdna["frag_max"],
        ss_list=", ".join(str(s) for s in STRAND_SPECIFICITIES),
        gf_list=", ".join(str(g) for g in GDNA_FRACTIONS),
    )
    yaml_path = outdir / f"{name}.yaml"
    yaml_path.write_text(content)
    return name, yaml_path


def run_sweep(yaml_path, sweep_outdir):
    """Run one sweep config."""
    cmd = [
        sys.executable, "scripts/synthetic_sim_sweep.py",
        "-c", str(yaml_path),
        "-o", str(sweep_outdir),
        "-v",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}", file=sys.stderr)
        return False
    return True


def combine_results(outdir, pairs):
    """Combine all sweep_results.tsv into one master TSV."""
    combined = outdir / "ablation_results.tsv"
    all_rows = []
    fieldnames = None

    for name, _ in pairs:
        tsv = outdir / name / "sweep_results.tsv"
        if not tsv.exists():
            print(f"  MISSING: {tsv}")
            continue
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if fieldnames is None:
                fieldnames = ["rna_dist", "gdna_dist"] + reader.fieldnames
            for row in reader:
                # Parse dist labels from name
                rna_label, gdna_label = name.split("__")
                row["rna_dist"] = rna_label
                row["gdna_dist"] = gdna_label
                all_rows.append(row)

    with open(combined, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nCombined {len(all_rows)} rows → {combined}")
    return combined


def analyze(combined_tsv):
    """Print ablation analysis summary."""
    with open(combined_tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    print("\n" + "=" * 90)
    print("ABLATION STUDY: Fragment Length Distribution Recovery")
    print(f"Total runs: {len(rows)}")
    print("=" * 90)

    # ── gDNA model recovery by (rna_dist, gdna_dist, strand_specificity) ──
    print("\n" + "─" * 90)
    print("gDNA MODEL MEAN RECOVERY")
    print("─" * 90)
    print(f"{'RNA dist':<14} {'gDNA dist':<14} {'SS':>5} {'gdna_f':>7} "
          f"{'true':>6} {'est':>7} {'bias':>7} {'est_std':>8} {'true_std':>9} {'n_obs':>6}")
    print("-" * 90)

    # Group by (rna_dist, gdna_dist, ss)
    groups = {}
    for r in rows:
        if not r.get("gdna_fl_est_mean"):
            continue
        key = (r["rna_dist"], r["gdna_dist"], r["strand_specificity"])
        groups.setdefault(key, []).append(r)

    for key in sorted(groups):
        rna_d, gdna_d, ss = key
        vals = groups[key]
        true_mean = float(vals[0]["gdna_fl_true_mean"])
        true_std = float(vals[0]["gdna_fl_true_std"])
        ests = [float(v["gdna_fl_est_mean"]) for v in vals]
        est_stds = [float(v["gdna_fl_est_std"]) for v in vals]
        n_obs = [int(v["gdna_fl_n_obs"]) for v in vals]
        avg_est = sum(ests) / len(ests)
        avg_std = sum(est_stds) / len(est_stds)
        avg_n = sum(n_obs) / len(n_obs)
        bias = avg_est - true_mean
        fracs = ", ".join(sorted(set(v["gdna_fraction"] for v in vals)))
        print(f"{rna_d:<14} {gdna_d:<14} {ss:>5} {fracs:>7} "
              f"{true_mean:>6.0f} {avg_est:>7.1f} {bias:>+7.1f} "
              f"{avg_std:>8.1f} {true_std:>9.1f} {avg_n:>6.0f}")

    # ── Summary table: gDNA bias by (FL pair, SS) ──
    print("\n" + "─" * 90)
    print("SUMMARY: Average gDNA Mean Bias by Distribution Pair × Strand Specificity")
    print("─" * 90)

    # Collect unique values
    all_rna = sorted(set(r["rna_dist"] for r in rows))
    all_gdna = sorted(set(r["gdna_dist"] for r in rows))
    all_ss = sorted(set(r["strand_specificity"] for r in rows), key=float)

    # Header
    header = f"{'RNA':<14} {'gDNA':<14} {'sep':>5}"
    for ss in all_ss:
        header += f" {'ss='+ss:>8}"
    header += f" {'n_est':>6}"
    print(header)
    print("-" * len(header))

    for rna_d, gdna_d in itertools.product(all_rna, all_gdna):
        # Compute separation in std units
        rna_mean = float(next(
            (r["rna_fl_true_mean"] for r in rows
             if r["rna_dist"] == rna_d), 0))
        gdna_mean = float(next(
            (r["gdna_fl_true_mean"] for r in rows
             if r["gdna_dist"] == gdna_d), 0))
        gdna_std = float(next(
            (r["gdna_fl_true_std"] for r in rows
             if r["gdna_dist"] == gdna_d), 1))
        rna_std = float(next(
            (r["rna_fl_true_std"] for r in rows
             if r["rna_dist"] == rna_d), 1))
        sep = abs(gdna_mean - rna_mean) / max(gdna_std, rna_std, 1)

        line = f"{rna_d:<14} {gdna_d:<14} {sep:>5.1f}"
        n_est_total = 0
        for ss in all_ss:
            subset = [
                r for r in rows
                if r["rna_dist"] == rna_d
                and r["gdna_dist"] == gdna_d
                and r["strand_specificity"] == ss
                and r.get("gdna_fl_est_mean")
            ]
            if subset:
                bias = sum(
                    float(r["gdna_fl_est_mean"]) - float(r["gdna_fl_true_mean"])
                    for r in subset
                ) / len(subset)
                n_est_total += len(subset)
                line += f" {bias:>+8.1f}"
            else:
                line += f" {'n/a':>8}"
        line += f" {n_est_total:>6}"
        print(line)

    # ── RNA model bias summary ──
    print("\n" + "─" * 90)
    print("SUMMARY: Average RNA Mean Bias by Distribution Pair × Strand Specificity")
    print("─" * 90)
    header = f"{'RNA':<14} {'gDNA':<14} {'sep':>5}"
    for ss in all_ss:
        header += f" {'ss='+ss:>8}"
    print(header)
    print("-" * len(header))

    for rna_d, gdna_d in itertools.product(all_rna, all_gdna):
        rna_mean = float(next(
            (r["rna_fl_true_mean"] for r in rows
             if r["rna_dist"] == rna_d), 0))
        gdna_mean = float(next(
            (r["gdna_fl_true_mean"] for r in rows
             if r["gdna_dist"] == gdna_d), 0))
        gdna_std = float(next(
            (r["gdna_fl_true_std"] for r in rows
             if r["gdna_dist"] == gdna_d), 1))
        rna_std = float(next(
            (r["rna_fl_true_std"] for r in rows
             if r["rna_dist"] == rna_d), 1))
        sep = abs(gdna_mean - rna_mean) / max(gdna_std, rna_std, 1)

        line = f"{rna_d:<14} {gdna_d:<14} {sep:>5.1f}"
        for ss in all_ss:
            subset = [
                r for r in rows
                if r["rna_dist"] == rna_d
                and r["gdna_dist"] == gdna_d
                and r["strand_specificity"] == ss
                and r.get("rna_fl_est_mean")
            ]
            if subset:
                bias = sum(
                    float(r["rna_fl_est_mean"]) - float(r["rna_fl_true_mean"])
                    for r in subset
                ) / len(subset)
                line += f" {bias:>+8.1f}"
            else:
                line += f" {'n/a':>8}"
        print(line)

    # ── gDNA model population rate ──
    print("\n" + "─" * 90)
    print("gDNA MODEL POPULATION RATE (runs with non-empty gDNA model / total)")
    print("─" * 90)
    header = f"{'RNA':<14} {'gDNA':<14}"
    for ss in all_ss:
        header += f" {'ss='+ss:>8}"
    print(header)
    print("-" * len(header))

    for rna_d, gdna_d in itertools.product(all_rna, all_gdna):
        line = f"{rna_d:<14} {gdna_d:<14}"
        for ss in all_ss:
            subset = [
                r for r in rows
                if r["rna_dist"] == rna_d
                and r["gdna_dist"] == gdna_d
                and r["strand_specificity"] == ss
            ]
            populated = sum(1 for r in subset if r.get("gdna_fl_est_mean"))
            total = len(subset)
            line += f" {populated}/{total:>5}" if total else f" {'n/a':>8}"
        print(line)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="FL distribution ablation study")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Skip sweep, just analyze existing results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pairs = []
    for rna, gdna in itertools.product(RNA_DISTS, GDNA_DISTS):
        name, yaml_path = generate_yaml(rna, gdna, outdir)
        pairs.append((name, yaml_path))

    if not args.analyze_only:
        total = len(pairs)
        for i, (name, yaml_path) in enumerate(pairs, 1):
            sweep_out = outdir / name
            print(f"\n[{i}/{total}] {name}")
            if (sweep_out / "sweep_results.tsv").exists():
                print("  Already completed, skipping.")
                continue
            ok = run_sweep(yaml_path, sweep_out)
            if ok:
                print(f"  Done.")
            else:
                print(f"  FAILED")

    combined = combine_results(outdir, pairs)
    analyze(combined)


if __name__ == "__main__":
    main()
