# NEXT SESSION — develop the message-propagation policy on the test chromosome

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0/§1, the
    rulings `DESIGN.md`, the lessons `TRAPS.md`, the substrate `TESTING.md` §0a/§0b.
    ⛔ Delete this file when the plan it points at is executed.

    Written 2026-08-19 at the end of a long session, for the migration the owner announced.

## ⭐⭐⭐ FIRST COMMAND OF THE SESSION — prove you can run and regenerate everything

    python scripts/design/preflight.py          # ~2 min, or --fast to skip the instrument sweep
    python -m pytest tests/ -q                  # 0 failed / 3,554 passed / 9 xfailed, 3,563 collected

⛔ **Do not start work on a ✘.** Preflight checks the toolchain, both references, both panels (five
oracle partitions and a certified `slot_truth` each) and every instrument, and for anything missing it
prints the command that regenerates it — a missing DERIVED artifact is a command you have not run, not
damage. ⭐ Its first run caught the 8 test-chromosome scenarios sitting uncertified.

## WHERE THINGS STAND

✅ **Simulator** — takes its transcriptome from a rigel INDEX, so nascent RNA is the index's own
single-exon ENTITY and capture binds by genomic overlap, strand-agnostically (`TESTING.md` §0).
✅ **Panel** — 16 conditions, rebuilt in 24.5 min, **nascent ON in every one**, so no row is
`nrna = 0` restated. 6/6 simulator gates; 16/16 COMPOSITION + FIELD certified (12 real per-reference
field verifications — the capture-ON and zero-gDNA rows are VACUOUS by design).
✅ **Per-arm truth** — the oracle cache carries five partitions: the three `ORIGINS` plus
`rna_pos`/`rna_neg` (RNA by TRANSCRIPT strand), gated `n_rna_pos + n_rna_neg == n_mrna + n_nrna` at
every slot, max gap **0**.
✅ **The hop-currency MAP, per arm** (`hop_currency.py`, 16 conditions in 23 s) — Stage 2 of the
rebuild is DONE. The arms genuinely disagree; the map is `ROADMAP.md` §0.
✅ **The standing BENCHMARK** — `relay_pool_ab.py --out` + `benchmark_report.py` (HTML). Artifacts:
`<ladder>/benchmark/baseline_2026-08-19_shipped.{log,tsv}` and `benchmark_ladder.html`.
✅ **The method-development TEST CHROMOSOME** — 8 transcripts over two stages, `TESTING.md` §0a.

## WHAT THE NEXT SESSION DOES

⭐⭐⭐ **STAGE 3 — the THIRD POLICY, developed on the test chromosome, one transcript at a time.**
The arithmetic it must implement is `EQUATIONS.md` §3.5e (the owner's worked SPLICE OUT / SPLICE IN and
the terminus case); the rulings are `DESIGN.md` §0c.0; the hop table is keyed on
`object class × {sj, term, sj+term}` (`TRAPS: an-object-class-does-not-see-a-terminus`).

**The loop, and it is seconds long:**

    python scripts/sim/build_test_reference.py                     # after editing the GTF/abundances
    rigel index --fasta <T>/test_chr.fa --gtf <T>/test_chr.gtf --no-mappability --no-tsv -o <T>/idx
    python scripts/sim/simulate_reads.py --config scripts/sim/configs/test_reference.yaml -j 8
    python scripts/design/build_scan_cache.py --index <T>/idx --suite <T>/scenarios
    # the 4 g00 rows need the per-condition path — `pass0_vs_oracle.py` holds zero-gDNA rows out:
    python scripts/design/pass0_vs_oracle.py --suite <T>/scenarios --index <T>/idx \
        --oracle-cache <T>/scenarios/oracle_cache --_prewarm <condition>
    python scripts/design/relay_pool_ab.py --suite <T>/scenarios --index <T>/idx \
        --oracle-cache <T>/scenarios/oracle_cache --out rows.tsv
    python scripts/design/benchmark_report.py rows.tsv --transcripts <T>/scenarios -o page.html

with `<T> = ~/Downloads/rigel_runs/test_reference`. Simulation of all 8 scenarios is **~18 s**.

⛔ **REPORT EVERY SCENARIO. DO NOT POOL** (owner, 2026-08-19: *"DON'T POOL SCENARIOS. WE NEED TO SEE
EVERY SCENARIO. ONLY POOL AT THE END."*). `benchmark_report.py` enforces the shape.

## THE MEASUREMENT THAT SHOULD DRIVE THE DESIGN

On the test chromosome the relay HELPS on all 8 scenarios; on the 16-condition ladder it COSTS the
three in-scope strata **1.4–1.7×** while winning the control (0.124×) and the deferred stratum (0.324×).
⭐ **That gap is the most useful thing on the table**: whatever `RelayPolicy` gets wrong on the real
panel is NOT yet represented on the test chromosome. The next structures to add are the ones that would
close it — overlapping genes on opposite strands, alternate TSS/TES, a long multi-exon gene whose exons
are many hops apart (the ladder puts 23–29 % of imputed mass ≥ 9 hops from any measured gDNA).

## THE TRAPS THIS SESSION PAID FOR

* ⛔ **A donor's LIBRARY-level priors must match the condition's strand/capture axes.** A `ss_0.99`
  donor injected into `ss_0.50` reported **82,581** false positives where a matching donor reports 0.
  `relay_pool_ab.py --donor` now refuses a mismatch. ⚠ Stage 2 removed the need entirely — with real
  splice junctions κ fits from the data.
* ⛔ **A test fixture that builds its subject with `__new__` tests the rebuild, not the subject.** One
  bypassed `__init__` and 242 tests failed on a new attribute while production was fine.
* ⛔ **A duplicate key in a dict literal is a LOST EDIT, not an error.** A jargon-gate allowance silently
  did nothing because the key already existed later in the same dict.
* ⛔ **`--collect-only` every count, never adjust one.** Every suite move this session was accounted.

## STILL OPEN

* The `ROADMAP.md` §1 rank 3 spike-and-slab: at a TERMINUS into an EXON under capture NEITHER currency
  works (~10 % excess on every arm) — the residual the map names and no policy choice fixes.
* `relay_pool_ab.py`'s docstring promises a `--table pipeline` (the EM-level nascent/annotated split)
  that does not exist. Build it or remove the promise.
* Prose in ~18 instruments still describes the retired 36-condition ladder and the deleted
  `pilot`/`flgap` panels. ⚠ Historical measurement stamps are PROVENANCE and must stay; what is worth
  fixing is any prose that presents a retired panel as current. All LIVE CODE and usage examples were
  repointed on 2026-08-19 (24 substitutions, verified by an AST scan that now returns empty).
