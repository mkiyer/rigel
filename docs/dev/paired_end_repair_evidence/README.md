# Evidence for `PAIRED_END_REPAIR.md` (2026-10-15/16)

Diagnosis only: every run reads an existing wheel and writes outside the repository's source.

| file | what it is |
|---|---|
| `arms.json` | one row per run: build, condition, which law or ruler was injected, transcript and gene error, pools, the fitted gDNA law, the ruler's spread |
| `laws/<run>.npz` | each run's fitted length laws (gDNA uniform-frame and realized, RNA, global) and the five raw pool histograms |
| `rulers/<run>.npz` | the EM's capture lengths a run produced (`own`), its `fl` lengths, and calibration's region and boundary efficiencies |
| `contrast/<input>.json` | the gDNA length contrast's internals: purities, separation, internal density, pools before and after the drain |
| `stage_run.py` | runs one build on one condition, pinned and fractional; `RULER_FROM=<run dir>`, `RULER_FROM=fl` or `GDNA_PMF_FROM=<run dir>` injects one input |
| `compare_arms.py` | splits a transcript-error change into within-gene and between-gene parts |
| `ruler_shape.py` | each ruler's within-gene shape against the simulator's capture factor |
| `audit_contrast.py` | the reviewer's contrast audit, parameterised by BAM |
| `counterexample_baseline.py` | the design's counterexamples through the baseline resolver |
| `mate_disagreement_census.py` | how often a pair's primary mates disagree about a junction, by kind and true origin |
| `census_*.json` | its results: minimap2 `g50 ss.99` OFF and ON (whole BAMs), VCaP (first 3 M templates) |

Builds: the supplied wheels in `$RIGEL_RUNS/prototypes/2026-10-15_pairing/site_{base,proto,s2,s3b}`
(stage 1 is `site_proto`) and the continuation's `<scratch>/rigel-pairing-finish/site_{endpoints,final}`. A run:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel && export OMP_NUM_THREADS=1
P=$RIGEL_RUNS/prototypes/2026-10-15_pairing
PYTHONPATH=$P/harness RIGEL_PROTO_SITE=$P/site_proto GDNA_PMF_FROM=<oracle run dir> \
  python stage_run.py <out dir> gdna_g50_ss_0.99_nrna_mid_capture_on
```

`SUITE=$RIGEL_RUNS/suite/ladder` points a run at the oracle BAMs instead of the aligned panel. Each run
takes 4–6 minutes on the panel and about 1.5 minutes on an oracle BAM, at about 1.3 GB.
