"""Check what SJ strand tags minimap2 writes on spliced reads."""
from hulkrna.sim import Scenario, SimConfig
import tempfile
import pysam
from pathlib import Path

tmp = Path(tempfile.mkdtemp())
sc = Scenario("tag_check", genome_length=5000, seed=42, work_dir=tmp / "check")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (2000, 2300)],
     "abundance": 100},
])
r = sc.build(n_fragments=200, sim_config=SimConfig(seed=42))

# Check tags on spliced reads
n_spliced = 0
n_xs = 0
n_ts = 0
n_both = 0
xs_correct = 0
ts_correct = 0
ts_needs_flip = 0

with pysam.AlignmentFile(str(r.bam_path), "rb") as bam:
    for read in bam:
        if read.is_unmapped or not read.cigartuples:
            continue
        has_sj = any(op == 3 for op, _ in read.cigartuples)
        if not has_sj:
            continue
        n_spliced += 1
        tags = dict(read.get_tags())
        has_XS = "XS" in tags
        has_ts = "ts" in tags

        if has_XS:
            n_xs += 1
        if has_ts:
            n_ts += 1
        if has_XS and has_ts:
            n_both += 1

        if n_spliced <= 15:
            xs = tags.get("XS", "N/A")
            ts = tags.get("ts", "N/A")
            rev = read.is_reverse
            # For a + strand gene, XS should be "+"
            # If ts is alignment-relative:
            #   forward read on + gene: ts should be "+"
            #   reverse read on + gene: ts should be "-" (opposite)
            xs_ok = xs == "+"
            if xs_ok:
                xs_correct += 1

            # Check if ts gives correct ref strand directly
            ts_gives_ref = ts == "+"
            # Check if ts needs flip (ts is alignment-relative)
            ts_flipped = (ts == "+" and not rev) or (ts == "-" and rev)

            print(
                f"  read={read.query_name[:25]:25s} rev={rev!s:5s} "
                f"XS={xs:3s} ts={ts:3s} "
                f"XS_correct={xs_ok!s:5s} "
                f"ts_direct={ts_gives_ref!s:5s} "
                f"ts_flipped={ts_flipped!s:5s}"
            )

print(f"\nTotal spliced reads: {n_spliced}")
print(f"  with XS: {n_xs} ({n_xs/max(n_spliced,1)*100:.0f}%)")
print(f"  with ts: {n_ts} ({n_ts/max(n_spliced,1)*100:.0f}%)")
print(f"  with both: {n_both}")

# Detailed check on ts semantics
print("\n--- ts tag semantics check ---")
n_ts_direct_correct = 0
n_ts_flip_correct = 0
with pysam.AlignmentFile(str(r.bam_path), "rb") as bam:
    for read in bam:
        if read.is_unmapped or not read.cigartuples:
            continue
        has_sj = any(op == 3 for op, _ in read.cigartuples)
        if not has_sj:
            continue
        tags = dict(read.get_tags())
        if "ts" not in tags:
            continue
        ts = tags["ts"]
        rev = read.is_reverse
        # Gene is on + strand, so reference SJ strand should be +
        # Test: does ts directly give "+"?
        if ts == "+":
            n_ts_direct_correct += 1
        # Test: does flipping for reverse reads give "+"?
        if (ts == "+" and not rev) or (ts == "-" and rev):
            n_ts_flip_correct += 1

print(f"  ts directly gives correct ref strand: {n_ts_direct_correct}/{n_ts}")
print(f"  ts after alignment-flip gives correct: {n_ts_flip_correct}/{n_ts}")
if n_ts_direct_correct > n_ts_flip_correct:
    print("  ==> ts appears REFERENCE-RELATIVE (like XS)")
elif n_ts_flip_correct > n_ts_direct_correct:
    print("  ==> ts appears ALIGNMENT-RELATIVE (needs flip)")
else:
    print("  ==> INCONCLUSIVE (both equal)")

sc.cleanup()
