#!/usr/bin/env python3
"""Diagnostic: trace mix_models internals from the sweep infrastructure."""
import sys
import yaml
import tempfile
from pathlib import Path
import numpy as np
import pysam

from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from rigel.pipeline import run_pipeline

sys.path.insert(0, str(Path(__file__).parent))
from synthetic_sim_sweep import parse_nrna_coords

with open("scripts/frag_length_validation.yaml") as f:
    cfg = yaml.safe_load(f)

read_length = cfg["read_length"]
rna_fp = cfg["rna"]
gdna_fp = cfg["gdna"]
transcripts_config = cfg["transcripts"]
nrna_config = cfg.get("nrna_coords", {})
nrna_groups = parse_nrna_coords(nrna_config, transcripts_config)
patterns_config = cfg.get("patterns", [])

params = dict(patterns_config[0])
params["gdna_fraction"] = 0.20
params["n_rna_fragments"] = 10000
params["strand_specificity"] = 1.0

n_rna = int(params["n_rna_fragments"])
gdna_frac = float(params["gdna_fraction"])
n_gdna = int(n_rna * gdna_frac)
ss = float(params["strand_specificity"])
seed = 42

with tempfile.TemporaryDirectory(prefix="rigel_diag_") as tmpdir:
    sc = Scenario(
        "diag",
        genome_length=cfg.get("genome_length", 100000),
        seed=seed,
        work_dir=Path(tmpdir) / "work",
    )

    for grp_label, grp in nrna_groups.items():
        nrna_ab = float(params.get(grp_label, 0))
        t_list = []
        for t_id in grp.t_ids:
            t_def = transcripts_config[t_id]
            ab = float(params.get(t_id, 0))
            nrna = nrna_ab if t_id == grp.carrier else 0.0
            t_list.append({
                "t_id": t_id,
                "exons": [tuple(e) for e in t_def["exons"]],
                "abundance": ab,
                "nrna_abundance": nrna,
            })
        sc.add_gene(grp_label, grp.strand, t_list)

    gdna_config = GDNAConfig(
        abundance=1.0,
        frag_mean=gdna_fp.get("frag_mean", 350),
        frag_std=gdna_fp.get("frag_std", 100),
        frag_min=gdna_fp.get("frag_min", 100),
        frag_max=gdna_fp.get("frag_max", 1000),
    )
    sim_cfg = SimConfig(
        frag_mean=rna_fp.get("frag_mean", 200),
        frag_std=rna_fp.get("frag_std", 30),
        frag_min=rna_fp.get("frag_min", 80),
        frag_max=rna_fp.get("frag_max", 450),
        read_length=read_length,
        strand_specificity=ss,
        seed=seed,
    )

    result = sc.build_oracle(
        n_fragments=0,
        sim_config=sim_cfg,
        gdna_config=gdna_config,
        nrna_abundance=0.0,
        n_rna_fragments=n_rna,
        gdna_fraction=gdna_frac,
    )

    print(f"Simulation: n_rna={n_rna}, n_gdna={n_gdna}, total={n_rna+n_gdna}")
    print(f"  RNA frag: mean={rna_fp['frag_mean']}, std={rna_fp['frag_std']}")
    print(f"  gDNA frag: mean={gdna_fp['frag_mean']}, std={gdna_fp['frag_std']}")
    print(f"  strand_specificity={ss}")

    # Count gDNA strands from oracle BAM read names
    n_gdna_fwd = 0
    n_gdna_rev = 0
    seen = set()
    with pysam.AlignmentFile(str(result.bam_path), "rb") as bam:
        for read in bam:
            qname = read.query_name
            if not qname.startswith("gdna:"):
                continue
            if qname in seen:
                continue
            seen.add(qname)
            strand_char = qname.split(":")[2]
            if strand_char == "f":
                n_gdna_fwd += 1
            else:
                n_gdna_rev += 1

    n_gdna_bam = n_gdna_fwd + n_gdna_rev
    print(f"\n--- Oracle BAM gDNA strand counts ---")
    print(f"  Forward (+): {n_gdna_fwd}")
    print(f"  Reverse (-): {n_gdna_rev}")
    print(f"  Total: {n_gdna_bam}")
    print(f"  Fwd ratio: {n_gdna_fwd/n_gdna_bam:.3f}")

    # Run pipeline
    pr = run_pipeline(result.bam_path, result.index)

    fl = pr.frag_length_models
    from rigel.splice import SpliceType

    print(f"\n--- Fragment Length Model Observations ---")
    print(f"  Global: {fl.global_model.total_weight:.0f}")
    spliced = fl.category_models[SpliceType.SPLICED_ANNOT].total_weight
    unspliced = fl.category_models[SpliceType.UNSPLICED].total_weight
    print(f"  SPLICED_ANNOT: {spliced:.0f}")
    print(f"  UNSPLICED: {unspliced:.0f}")
    print(f"  Unspliced same-strand: {fl.unspliced_same_strand.total_weight:.0f} "
          f"mean={fl.unspliced_same_strand.mean:.1f}")
    print(f"  Unspliced opp-strand: {fl.unspliced_opp_strand.total_weight:.0f} "
          f"mean={fl.unspliced_opp_strand.mean:.1f}")
    print(f"  Intergenic: {fl.intergenic.total_weight:.0f} "
          f"mean={fl.intergenic.mean:.1f}")

    same_opp = fl.unspliced_same_strand.total_weight + fl.unspliced_opp_strand.total_weight
    print(f"\n  same+opp={same_opp:.0f}, UNSPLICED={unspliced:.0f}, "
          f"diff={unspliced - same_opp:.0f} (unspliced NOT classified into same/opp)")

    S = fl.unspliced_same_strand.total_weight
    A = fl.unspliced_opp_strand.total_weight
    sspec = pr.strand_models.strand_specificity
    p1s = pr.strand_models.exonic_spliced.p_r1_sense

    print(f"\n--- Mixing Inputs ---")
    print(f"  S (sense)={S:.0f}, A (antisense)={A:.0f}")
    print(f"  strand_specificity={sspec:.4f}, p_r1_sense={p1s:.4f}")

    denom = 2.0 * sspec - 1.0
    R_est = max(0.0, (S - A) / denom) if denom > 1e-6 else 0
    G_est = (S + A) - R_est
    W_sense = min(max((G_est * 0.5 / S) if S > 0 else 0, 0), 1)
    W_anti = min(max((G_est * 0.5 / A) if A > 0 else 0, 0), 1)

    print(f"  R_est={R_est:.1f}, G_est={G_est:.1f}")
    print(f"  W_sense={W_sense:.4f}, W_anti={W_anti:.4f}")

    n_inter = fl.intergenic.total_weight
    genic_gdna = n_gdna - n_inter
    print(f"\n--- Expected vs Estimated ---")
    print(f"  True gDNA total={n_gdna}, intergenic={n_inter:.0f}, genic={genic_gdna:.0f}")
    print(f"  Expected gDNA antisense ~= {0.5 * genic_gdna:.0f}")
    print(f"  Actual A (antisense pool)={A:.0f}")

    print(f"\n--- Mixed Model Results ---")
    print(f"  gDNA model: weight={fl.gdna_model.total_weight:.0f}, "
          f"mean={fl.gdna_model.mean:.1f}, std={fl.gdna_model.std:.1f}")
    print(f"  RNA model: weight={fl.rna_model.total_weight:.0f}, "
          f"mean={fl.rna_model.mean:.1f}, std={fl.rna_model.std:.1f}")

    H_anti = fl.unspliced_opp_strand.counts.copy()
    H_sense = fl.unspliced_same_strand.counts.copy()
    H_inter = fl.intergenic.counts.copy()

    x = np.arange(len(H_anti))
    anti_mean = np.average(x, weights=H_anti) if H_anti.sum() > 0 else 0
    inter_mean = np.average(x, weights=H_inter) if H_inter.sum() > 0 else 0
    sg = W_sense * H_sense
    sg_mean = np.average(x, weights=sg) if sg.sum() > 0 else 0

    print(f"\n--- Component Analysis ---")
    print(f"  Antisense pool mean (should be ~350): {anti_mean:.1f}")
    print(f"  Intergenic pool mean: {inter_mean:.1f}")
    print(f"  W_sense * H_sense mean: {sg_mean:.1f}")

    at = (W_anti * H_anti).sum()
    st = sg.sum()
    it = H_inter.sum()
    tot = at + st + it
    cg = H_inter + W_sense * H_sense + W_anti * H_anti
    cm = np.average(x, weights=cg) if tot > 0 else 0

    print(f"\n--- gDNA estimate decomposition ---")
    print(f"  Intergenic: {it:.0f} ({100*it/tot:.1f}%)")
    print(f"  Antisense:  {at:.0f} ({100*at/tot:.1f}%)")
    print(f"  Sense:      {st:.0f} ({100*st/tot:.1f}%)")
    print(f"  combined gDNA mean: {cm:.1f} (true: 350)")
