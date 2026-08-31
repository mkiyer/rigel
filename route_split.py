"""THE TWO-ESTIMAND ROUTING: geometry eats the uniform-frame law (what ships), the SCORER eats the
integrated realized law. Prediction from the channel split: the capture-ON wins survive and the +150k
geometry regression never happens, because geometry is left in its own frame."""
import importlib.util, sys, dataclasses
from pathlib import Path
import numpy as np
D = Path("/Users/mkiyer/proj/rigel/scripts/design"); sys.path.insert(0, str(D))
sys.path.insert(0, "/private/tmp/claude-503/-Users-mkiyer-proj-rigel/51e9e3c4-4f62-4adf-bf01-2a5434d5341a/scratchpad")
def sib(n):
    k = n[:-3]
    if k not in sys.modules:
        sp = importlib.util.spec_from_file_location(k, D/n); m = importlib.util.module_from_spec(sp)
        sys.modules[k] = m; sp.loader.exec_module(m)
    return sys.modules[k]
QA = sib("quant_accuracy.py")
import integrated_fl as IF
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
import rigel.calibration.fl as FL
import rigel.pipeline as PL
from rigel.frag_length_model import FragmentLengthModel
from rigel.calibration.splice_graph import build_region_partition_arrays
from rigel.calibration.gdna_density import region_lengths_from_partition

panel, idx = Path(sys.argv[1]), sys.argv[2]
index = TranscriptIndex.load(idx); pc = PipelineConfig()
b, off, rt = build_region_partition_arrays(index)
rl = region_lengths_from_partition(b, off, len(rt))
real_build, real_score = FL.build_fl_models, PL._score_fragments

def run_routed(cond):
    box = {"g": None}
    def build(payload, **kw):
        m = real_build(payload, **kw)
        if kw.get("gdna_opportunity") is not None and box["g"] is None:
            gi, _ = IF.integrated(payload, index, kw["gdna_opportunity"], rl, rt,
                                  payload.ref_region_offsets, payload.ref_boundary_offsets, m.rna_pmf)
            if gi is not None:
                g = np.asarray(gi, float)[:1001]; box["g"] = g / g.sum()
        return m                      # geometry keeps the fitted (uniform-frame) pair
    def score(buffer, index_, sm, rna_fl, gdna_fl, *a, **kw):
        if box["g"] is not None:
            gdna_fl = FragmentLengthModel.from_pmf(box["g"], 1000)   # scorer gets the realized law
        return real_score(buffer, index_, sm, rna_fl, gdna_fl, *a, **kw)
    FL.build_fl_models, PL._score_fragments = build, score
    try:
        return {x["axis"]: x for x in QA.run_condition("base", panel, index, cond, pc, None)}, box
    finally:
        FL.build_fl_models, PL._score_fragments = real_build, real_score

print(f"{'condition':44s} {'base tx':>12s} {'ROUTED tx':>12s} {'delta':>9s} {'base fp':>10s} {'ROUTED fp':>10s}")
BASE = {"gdna_g05_ss_0.99_nrna_mid_capture_on": (2570436, 482611),
        "gdna_g50_ss_0.99_nrna_mid_capture_on": (1667577, 393451),
        "gdna_g50_ss_0.99_nrna_mid_capture_off": (542188, 97376)}
for cond, (bt, bf) in BASE.items():
    out, box = run_routed(cond)
    assert box["g"] is not None, "integrated never produced a pmf"
    t = float(out["transcript"]["count_abs_err"]); f = float(out["transcript"]["fp_mass"])
    print(f"{cond:44s} {bt:12,.0f} {t:12,.0f} {t-bt:+9,.0f} {bf:10,.0f} {f:10,.0f}")
