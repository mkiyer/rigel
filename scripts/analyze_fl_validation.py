#!/usr/bin/env python3
"""Analyze fragment length distribution recovery from sweep results."""
import csv
import sys
from pathlib import Path

def main():
    tsv = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        '/Users/mkiyer/Downloads/rigel_runs/synthetic_sim/fl_validation_v2/sweep_results.tsv')
    with open(tsv) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    print("=" * 80)
    print("FRAGMENT LENGTH DISTRIBUTION RECOVERY ANALYSIS")
    print("=" * 80)
    print(f"\nTrue RNA FL: mean=200, std=30")
    print(f"True gDNA FL: mean=350, std=80")
    print(f"Total runs: {len(rows)}")

    # ---- RNA model ----
    print("\n" + "=" * 80)
    print("RNA FRAGMENT LENGTH MODEL")
    print("=" * 80)

    rna_nonempty = sum(1 for r in rows if r['rna_fl_est_mean'])
    print(f"Runs with RNA FL estimate: {rna_nonempty}/{len(rows)}")
    print(f"Runs with EMPTY rna_model: {len(rows) - rna_nonempty}/{len(rows)}")

    # RNA by gdna_fraction
    print("\n--- RNA Mean Recovery by gDNA Fraction ---")
    print(f"{'gdna_frac':>10} {'n_runs':>7} {'avg_est':>8} {'bias':>7} {'avg_std':>8} {'true_std':>8}")
    print("-" * 55)
    rna_by_frac = {}
    for r in rows:
        if not r['rna_fl_est_mean']:
            continue
        gf = float(r['gdna_fraction'])
        rna_by_frac.setdefault(gf, []).append({
            'est_mean': float(r['rna_fl_est_mean']),
            'est_std': float(r['rna_fl_est_std']),
            'n_obs': int(r['rna_fl_n_obs']),
        })
    for gf in sorted(rna_by_frac):
        vals = rna_by_frac[gf]
        means = [v['est_mean'] for v in vals]
        stds = [v['est_std'] for v in vals]
        avg_m = sum(means) / len(means)
        avg_s = sum(stds) / len(stds)
        print(f"{gf:>10.2f} {len(vals):>7} {avg_m:>8.1f} {avg_m-200:>+7.1f} {avg_s:>8.1f} {'30':>8}")

    # RNA by NTA1
    print("\n--- RNA Mean Recovery by NTA1 Level ---")
    print(f"{'NTA1':>5} {'n_runs':>7} {'avg_est':>8} {'bias':>7} {'avg_std':>8}")
    print("-" * 40)
    rna_by_nta1 = {}
    for r in rows:
        if not r['rna_fl_est_mean']:
            continue
        nta1 = r['NTA1']
        rna_by_nta1.setdefault(nta1, []).append(float(r['rna_fl_est_mean']))
    for nta1 in sorted(rna_by_nta1, key=lambda x: float(x)):
        vals = rna_by_nta1[nta1]
        avg_m = sum(vals) / len(vals)
        print(f"{nta1:>5} {len(vals):>7} {avg_m:>8.1f} {avg_m-200:>+7.1f}")

    # RNA detailed per-run
    print("\n--- RNA Detailed per-run ---")
    print(f"{'run':>4} {'gdna_f':>7} {'NTA1':>5} {'n_obs':>6} {'true':>6} {'est':>7} {'err':>7} {'t_std':>6} {'e_std':>7}")
    print("-" * 65)
    for i, r in enumerate(rows):
        n_obs = r.get('rna_fl_n_obs', '0')
        est = r.get('rna_fl_est_mean', '') or '-'
        err = r.get('rna_fl_mean_err', '') or '-'
        es = r.get('rna_fl_est_std', '') or '-'
        gf = float(r['gdna_fraction'])
        nta1 = r['NTA1']
        print(f"{i+1:>4} {gf:>7.2f} {nta1:>5} {n_obs:>6} {'200':>6} {est:>7} {err:>7} {'30':>6} {es:>7}")

    # ---- gDNA model ----
    print("\n" + "=" * 80)
    print("gDNA FRAGMENT LENGTH MODEL")
    print("=" * 80)

    gdna_nonempty = sum(1 for r in rows if r['gdna_fl_est_mean'])
    print(f"Runs with gDNA FL estimate: {gdna_nonempty}/{len(rows)}")
    print(f"Runs with EMPTY gdna_model: {len(rows) - gdna_nonempty}/{len(rows)}")

    # ---- gDNA recovery by contamination level ----
    print("\n--- gDNA Mean Recovery by Contamination Level ---")
    print(f"{'gdna_frac':>10} {'n_runs':>7} {'avg_est':>8} {'bias':>7} {'avg_std':>8} {'min_est':>8} {'max_est':>8}")
    print("-" * 65)

    by_frac = {}
    for r in rows:
        if not r['gdna_fl_est_mean']:
            continue
        gf = float(r['gdna_fraction'])
        by_frac.setdefault(gf, []).append({
            'est_mean': float(r['gdna_fl_est_mean']),
            'est_std': float(r['gdna_fl_est_std']),
            'n_obs': int(r['gdna_fl_n_obs']),
        })

    for gf in sorted(by_frac):
        vals = by_frac[gf]
        means = [v['est_mean'] for v in vals]
        stds = [v['est_std'] for v in vals]
        avg_m = sum(means) / len(means)
        avg_s = sum(stds) / len(stds)
        print(f"{gf:>10.2f} {len(vals):>7} {avg_m:>8.1f} {avg_m-350:>+7.1f} {avg_s:>8.1f} {min(means):>8.1f} {max(means):>8.1f}")

    # ---- Per-run detail for gDNA ----
    print("\n--- gDNA Detailed per-run ---")
    print(f"{'run':>4} {'gdna_f':>7} {'NTA1':>5} {'n_obs':>6} {'true':>6} {'est':>7} {'err':>7} {'t_std':>6} {'e_std':>7}")
    print("-" * 65)

    for i, r in enumerate(rows):
        gf = float(r['gdna_fraction'])
        if gf == 0:
            continue
        nta1 = r['NTA1']
        n_obs = r['gdna_fl_n_obs']
        est = r['gdna_fl_est_mean'] or '-'
        err = r['gdna_fl_mean_err'] or '-'
        es = r['gdna_fl_est_std'] or '-'
        print(f"{i+1:>4} {gf:>7.2f} {nta1:>5} {n_obs:>6} {'350':>6} {est:>7} {err:>7} {'80':>6} {es:>7}")

    # ---- Impact of nRNA on gDNA estimation ----
    print("\n--- Impact of nRNA Contamination on gDNA FL Estimation ---")
    print(f"{'gdna_f':>7} {'NTA1':>5} {'gdna_est_mean':>14} {'bias':>7} {'gdna_count_err':>15}")
    print("-" * 55)

    for r in sorted(rows, key=lambda x: (float(x['gdna_fraction']), float(x['NTA1']))):
        gf = float(r['gdna_fraction'])
        if gf == 0:
            continue
        nta1 = r['NTA1']
        est = r['gdna_fl_est_mean']
        if est:
            bias = f"{float(est)-350:+.1f}"
            est_s = f"{float(est):.1f}"
        else:
            bias = '-'
            est_s = '-'
        gdna_exp = float(r['gdna_expected'])
        gdna_obs = float(r['gdna_observed'])
        gdna_cerr = f"{gdna_obs - gdna_exp:+.1f}" if gdna_exp > 0 else '-'
        print(f"{gf:>7.2f} {nta1:>5} {est_s:>14} {bias:>7} {gdna_cerr:>15}")

    # ---- Summary of issues ----
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Compute summary metrics
    rna_ests = [float(r['rna_fl_est_mean']) for r in rows if r['rna_fl_est_mean']]
    gdna_ests = [float(r['gdna_fl_est_mean']) for r in rows if r['gdna_fl_est_mean']]
    
    if rna_ests:
        avg_rna = sum(rna_ests) / len(rna_ests)
        print(f"RNA model: {len(rna_ests)}/{len(rows)} populated, "
              f"avg estimated mean={avg_rna:.1f} (true=200, bias={avg_rna-200:+.1f})")
    else:
        print(f"RNA model: 0/{len(rows)} populated")
    
    if gdna_ests:
        avg_gdna = sum(gdna_ests) / len(gdna_ests)
        print(f"gDNA model: {len(gdna_ests)}/{len(rows)} populated, "
              f"avg estimated mean={avg_gdna:.1f} (true=350, bias={avg_gdna-350:+.1f})")
    else:
        print(f"gDNA model: 0/{len(rows)} populated")

if __name__ == "__main__":
    main()
