# Wave-1 independent review evidence

Review-only artifacts for `../PAIRED_END_WAVE1_INDEPENDENT_REVIEW.md`, at source commit
`41899ae4b8209893f7b5ae0cf893304e6796dabf`. No production implementation or installation.

- `audit.py`: finite query-map, geometry, support and probability counterexamples; native scorer probes;
  nine independently enumerated emission fixtures; a native resolver probe for individually annotated
  junctions unsupported by any common transcript.
- `audit.json`: the executed audit's output. Assertions check the stated counterexamples and baseline
  behavior. They are review probes, not a replacement for future failing scanner regression fixtures.
- `validation.txt`: targeted existing pytest results after the review artifacts were written.
- `runtime.json`: interpreter and loaded native-extension paths/hashes for reproduction.

Run in the existing `rigel` conda environment from the repository root:

```sh
OMP_NUM_THREADS=1 python -I docs/dev/paired_end_wave1_review_evidence/audit.py > docs/dev/paired_end_wave1_review_evidence/audit.json
OMP_NUM_THREADS=1 python -m pytest tests/test_pipeline_routing.py tests/test_em_pseudocounts.py tests/test_second_pass_scoring.py tests/test_second_pass_placements.py tests/test_transcript_space_fl.py tests/test_docs_boundary.py -q
```

These tests exercise the installed native extensions; the independent maps do not call the production
path helpers. Resolver inputs are already assembled blocks, and scorer inputs are already buffered
readings. Neither proves that the proposed scanner can build those inputs from BAM records.

The targeted suite has **185 passed, 1 failed**. The failure is
`tests/test_docs_boundary.py::test_no_numbered_rule_labels`, naming only pre-existing untracked files:

- `docs/dev/RNA_RELEASE_0_8_0_CAPTURE_REFERENCE_INDEPENDENT_REVIEW.md`: two numbered experiment labels;
- `docs/dev/RNA_RELEASE_0_8_0_MODEFREE_INDEPENDENT_REVIEW.md`: a numbered blocker label.

Those files were present in the initial workspace status and were not edited by this review. The new
report and audit are not listed by the failing check. No full-suite or BAM benchmark result is claimed.
The original baseline counterexample script was also executed successfully. Recorded VCaP/panel
measurements were read from the supplied documents and evidence; no library-level A/B was rerun.
