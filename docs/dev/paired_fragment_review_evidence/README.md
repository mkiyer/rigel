# Paired-fragment design review: reproducible evidence

2026-10-15 · Source checkout `41899ae4b8209893f7b5ae0cf893304e6796dabf`.
The installed `rigel` Python module resolves to this checkout. No source or compiled extension was changed.

- `paired_fragment_audit.py` creates six synthetic BAM fixtures in `/private/tmp`, scans them with the current
  native scanner, drains calibration, and reports lengths. Its result is `paired_fragment_results.json`.
- `length_contrast_audit.py` scans the existing minimap2 gDNA-half, stranded capture-OFF BAM, discards the EM buffer,
  and repeats calibration's drain and length fitting. It does not run EM. Paths are explicit in the script.
  `length_contrast.json` and `length_contrast_arrays.npz` preserve its output and intermediate arrays.
- `test_paired_fragment_contract.py` reruns the native synthetic experiment and asserts the required paired contract.
  Its **four expected baseline failures and two passing controls** are captured in `contract_tests.txt`. These are
  deliberately failing falsification tests, not a claim that a repair has passed. They are outside `tests/`.
- `docs_tests.txt` records the repository documentation check: ten passed, one failed because of numbered labels
  in three other review files. None of the files changed by this review is an offender. The additional plan-review
  file was present by the time the check ran; it was left untouched.

Both runs used the `rigel` environment, `OMP_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1`, and
`BamScanConfig(total_threads=1)`. Full-library audit elapsed time was 326.7 seconds. No new alignment, simulation
panel, production fix or end-to-end abundance comparison was performed.

The synthetic fixture establishes that incomplete reported alternatives become 100-bp fragment deposits, even
beside complete 300-bp primary pairs. A real pair at a remote location restores the correct length. The fixture's
reference and reads are deliberately simple; it tests representation and deposit behavior, not alignment search.

The full-library output establishes the drain's differential pool spike and the contrast's amplification. It does
not identify the provenance of every 100-bp library observation or measure how much a complete-pair repair improves
abundance. Those require a separate repair A/B after the design gate.

Pool order is intergenic contained, intronic contained, intron/exon crossing, intergenic/exon crossing, and spliced
RNA. The raw signed contrast is not a probability distribution: its negative mass is retained in this evidence.
`contrast_counts` precedes the later coupled/census and empirical-Bayes steps; `gdna_pmf` is the final uniform-frame
PMF. Small differences from the packet's empirical `gdna` histogram reflect these distinct stages. All raw pool
counts reproduced the supplied aligned length table.

The review also re-read the supplied oracle and aligned `fragment_lengths.tsv` tables to compute category totals,
means and fractions at 100 bp. Existing oracle mean 216.69 and aligned empirical gDNA mean 99.19 are comparisons of
those reported histograms, not a claim that a valid paired-acceptance correction must force the two means to match.

## Subsequent pairing-plan review

`mate_pairing_proposal_reviewed.txt` preserves the proposal before the revised execution plan was written.
`pair_geometry_audit.py` builds six additional synthetic paired-BAM cases in a temporary directory and runs the
current native scanner, with one scan thread. `pair_geometry_results.json` contains its output:

- a 4,100-bp pair in one shared exon reaches the EM buffer, although it deposits no length;
- an outward primary pair deposits 400 bp;
- a protruding overlap passes the proposal's facing predicate and deposits 150 bp;
- a spliced mate and a continuous mate that reads through its skipped intron deposit 125 bp, although the
  continuous mate alone aligns 150 bases;
- genuine paired 100-bp and clipped 50-bp controls retain those lengths.

These are geometry/representation probes, not biologically simulated sequence evidence or end-to-end EM runs.
No source repair was built. They establish that completeness alone and the proposal's shared-transcript/facing
shortcuts cannot certify compatible paired paths. The revised plan specifies independent path-enumeration and
mutation gates for the future implementation.

## 5′-clip policy review

`mate_pairing_before_clip_review.txt` preserves the plan before the blanket-withholding policy was reconsidered.
`five_prime_clip_audit.py` is an independent, pure-Python check of a conditional endpoint model, not a native fix.
`five_prime_clip_results.json` records its output. Under one templated path base per clipped query base, infer the
forward endpoint from its aligned 5′ anchor minus its clip, and the reverse endpoint from its aligned 5′ boundary
plus its clip. All 98,692 enumerated molecule/read/clip combinations recover their true endpoints and lengths;
adding clips to the aligned union span instead fails in 9,075 combinations, including a 100-bp overlapping pair
incorrectly turned into 110 bp.

Additional examples show an RNA extension crossing an intron, a candidate path with insufficient upstream sequence,
and failures of unconditional exactness with non-template prefixes or sequencing indels. These examples establish
the geometry under stated assumptions; they do not measure real-data clip causes or prove a production parser.
The owner specified the intended input scope: remove non-template 5′ sequence and preserve the insert's 5′ bases.

The review read the supplied `policy_census.py`, `vcap_split.py` and three `census_*.txt` outputs in the mate-pairing
prototype directory. These inspect primary pairs in the first 3M name groups and do not count final accepted-fragment
loss under the revised assembler. The review did not rerun that census or establish clip-to-junction enrichment.

## Review of the staged implementation

`pairing_stage3_audit.py` runs nine small synthetic probes against the supplied `site_s3b` wheel. Its results are in
`pairing_stage3_audit_results.json`, including the imported module path and SHA-256. It confirms a 300-bp ordinary
pair, and reproduces HI-internal duplicate NM inflation, geometry-only clip collapse, a mate with no M/=/X support,
a tolerated three-base splice read-through, conflicting-strand representative selection, and the observed-path
versus EM-length mismatch on a novel intergene splice. These test representation, not alignment accuracy.

`pairing_stage3_fate_audit.py` replays the two supplied disappeared-from-buffer BAM subsets through the native
assembly seam. `pairing_stage3_fates_mm2.json` finds 12,055 unpaired and 16 still-placed templates in the 12,071
minimap2 subset; `pairing_stage3_fates_vcap.json` finds 46,865 unpaired and 10 still-placed templates in the 46,875
VCaP subset. All those still-placed fragments have multiple placements and no transcript candidates. These subset
counts do not replace a full per-fragment outcome census.

`pairing_stage3_review_tests.txt` records a fresh focused run: 36 passed and the five documented stage-4 gates
failed. The full-suite and mutation results were read from the supplied logs, not rerun by this review. Neither
worktree source nor main-tree production source was edited or rebuilt. The complete review and proposed execution
order are in `../MATE_PAIRING_IMPLEMENTATION_REVIEW.md`.

To run the new probes, use the `rigel` Python executable with `PYTHONDONTWRITEBYTECODE=1`, `OMP_NUM_THREADS=1`,
`PYTHONPATH=$RIGEL_RUNS/prototypes/2026-10-15_pairing/harness`, and
`RIGEL_PROTO_SITE=$RIGEL_RUNS/prototypes/2026-10-15_pairing/site_s3b`. The harness removes the editable
installation's finder so the supplied wheel is actually imported. The fate audit takes `mm2` or `vcap` as its
argument. The scripts create only temporary synthetic inputs; capture stdout to retain their JSON results.

`pairing_stage3_review_docs_tests.txt` records ten documentation checks passing and the existing numbered-label
failure in two older review documents. None of this review's files is listed as an offender. All local links in
the implementation review were checked, and the main tree's production/source/test diff remains empty.
