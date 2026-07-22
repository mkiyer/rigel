"""Reproduce the antisense-intronic ss_0.65 leak and check the MECHANISM: is Fix #1's t2 leak just the phantom
gDNA (a sink) becoming correctly-RNA (which the EM then leaks to the antisense single-exon at weak strand), or a
genuine calibration over-shift? Prints t2 leak + the gDNA/nRNA pool totals. Run with fix on (current) and off
(git stash src/rigel/calibration/bp_solver.py) to A/B."""
import sys
from pathlib import Path
sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests")
from rigel.sim import Scenario
from scenarios.conftest import build_and_run

work = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/scratch_anti")
work.mkdir(parents=True, exist_ok=True)
sc = Scenario("anti_intron", genome_length=8000, seed=1234, work_dir=work / "anti_intron")
sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(1000, 1200), (2000, 2200), (5000, 5600)], "abundance": 100}])
sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0}])
sc.add_gene("g_ctrl", "+", [{"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0}])
bench = build_and_run(sc, nrna_abundance=50, strand_specificity=0.65, n_fragments=2000,
                      scenario_name="anti_intron_repro")
t2 = next(t for t in bench.transcripts if t.t_id == "t2")
t1 = next(t for t in bench.transcripts if t.t_id == "t1")
print("\n===== anti_intron ss=0.65 nrna=50 gdna=0 =====")
print(f"  t2 (NEG single-exon, expected 0) observed = {t2.observed:.0f}")
print(f"  t1 (POS multi-exon, expected ~mRNA)         = {t1.observed:.0f}")
print(f"  n_gdna_pipeline (PHANTOM — should be ~0)     = {bench.n_gdna_pipeline:.0f}")
print(f"  n_nrna_pipeline                              = {bench.n_nrna_pipeline:.0f}")
print(f"  n_gdna_expected                              = {bench.n_gdna_expected:.0f}")
sc.cleanup()
