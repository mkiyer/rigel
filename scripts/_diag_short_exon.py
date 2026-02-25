"""Diagnostic: T1 tiny unique exon + large shared exon, 10:1 ratio."""
from hulkrna.sim import Scenario, SimConfig, run_benchmark
from hulkrna.pipeline import run_pipeline
from collections import Counter
import tempfile
import pathlib
import pysam

tmp = tempfile.mkdtemp()
sc = Scenario("diag", genome_length=25000, seed=42,
              work_dir=pathlib.Path(tmp) / "diag")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(1000, 1050), (2000, 22000)], "abundance": 1000},
    {"t_id": "t2", "exons": [(1200, 1500), (2000, 22000)], "abundance": 100},
])
cfg = SimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                read_length=100, strand_specificity=1.0, seed=42)
result = sc.build_oracle(n_fragments=1000, sim_config=cfg)
gt = result.ground_truth_auto()

print("=== GEOMETRY ===")
for t in result.transcripts:
    el = sum(e - s for s, e in t.exons)
    print(f"  {t.t_id}: exons={t.exons}, exonic_len={el}bp")

print("\n=== GROUND TRUTH ===")
for tid, c in sorted(gt.items()):
    print(f"  {tid}: {c}")

# Spliced vs shared analysis
pair_data = {}
with pysam.AlignmentFile(str(result.bam_path), "rb") as bam:
    for read in bam:
        qn = read.query_name
        cigar = read.cigartuples or []
        has_splice = any(op == 3 for op, _ in cigar)
        if qn not in pair_data:
            pair_data[qn] = {
                "has_splice": has_splice,
                "min_start": read.reference_start,
                "max_end": read.reference_end or read.reference_start,
            }
        else:
            pair_data[qn]["has_splice"] |= has_splice
            pair_data[qn]["min_start"] = min(
                pair_data[qn]["min_start"], read.reference_start
            )
            pair_data[qn]["max_end"] = max(
                pair_data[qn]["max_end"],
                read.reference_end or read.reference_start,
            )

spliced = Counter()
shared_only = Counter()
total = Counter()
for qn, info in pair_data.items():
    tid = qn.split(":")[0]
    total[tid] += 1
    if info["has_splice"]:
        spliced[tid] += 1
    if info["min_start"] >= 2000 and info["max_end"] <= 22000:
        shared_only[tid] += 1

print("\n=== FRAGMENT PLACEMENT ===")
for tid in sorted(total):
    print(
        f"  {tid}: total={total[tid]}, "
        f"spliced={spliced.get(tid, 0)}, "
        f"shared_only={shared_only.get(tid, 0)}"
    )

pr = run_pipeline(result.bam_path, result.index, sj_strand_tag="ts", seed=42)
bench = run_benchmark(result, pr, scenario_name="big_shared_exon")

print(f"\n{bench.summary()}")
print("\n=== PER-TRANSCRIPT ===")
for ta in bench.transcripts:
    eff = float(pr.estimator.effective_lengths[ta.t_index])
    print(
        f"  {ta.t_id}: expected={ta.expected}, observed={ta.observed:.1f}, "
        f"abs_diff={ta.abs_diff:.1f}, rel_error={ta.rel_error:.2%}, "
        f"eff_len={eff:.0f}"
    )

sc.cleanup()
