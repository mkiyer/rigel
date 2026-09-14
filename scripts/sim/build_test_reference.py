#!/usr/bin/env python3
"""BUILD THE METHOD-DEVELOPMENT TEST REFERENCE — one chromosome described in ONE YAML file, everything
else derived from it.

⭐⭐⭐ **THIS IS THE BENCHMARK THE MESSAGE-PROPAGATION POLICY IS DEVELOPED AGAINST** (owner, 2026-08-19):
*"over the course of this method development, we will add transcripts to a 'test chromosome' … the test
reference chromosome, test transcript GTF, and test abundances comprise a critical method development
benchmark."* It grows one structure at a time, and the policy must solve the new structure AND
everything already there.

**THE ONE HAND-EDITED FILE (owner ruling, 2026-09-02): `test_chr.yaml`.** It is the `rigel sim`
scenario schema — ``genes → {gene_id, strand, transcripts: [{t_id, exons, abundance, nrna_abundance}]}``,
the same ``genes`` shape the toy harness's ``ToySpec`` uses, exons 0-based half-open — plus two keys the
scenario schema lacks: ``probed`` per gene (the capture panels are DESIGNED from it) and ``shadow_genes``
(unannotated transcription on the BLANK contig, simulator-only). ⛔ **Every other file in
`test_reference/` is RENDERED from it by this script and must never be hand-edited**:

===========================  ========================================================================
`test_chr.gtf`               one `exon` line per exon — what `rigel index` reads
`test_shadow.gtf`            the shadow transcripts, reference `test_blank` — what the SIMULATOR reads
`test_abundances.tsv`        `transcript_id / mrna_abundance / nrna_abundance` (no comment lines — both
                             readers key on the first line)
`test_probes.bed`            the BENIGN capture panel: each probed gene's exon UNION tiled 8 × 125 bp
                             seamless from each union piece's start, never spanning an sj
`test_probes_sparse.bed`     the SPARSE adversarial panel: one 125 bp probe CENTRED on every annotated
                             exon of a probed gene (faces depleted, mid-exon enriched)
`test_probes_junction.bed`   the JUNCTION adversarial panel: one BED12 two-block probe (62 + 63 bp,
                             contiguous only in cDNA) per sj of a probed gene, and nothing else
`test_chr.fa` + the index    written to the runs dir only (`--out`)
===========================  ========================================================================

The renders are written into the repo directory beside the YAML (versioned, so a reader sees the
annotation without running anything) AND copied into the runs dir. ``--check`` re-renders in memory and
refuses if any checked-in render has drifted from the YAML — the suite gate
(`tests/test_test_reference_renders.py`) runs exactly that.

⛔ **WHY THE FASTA IS DERIVED RATHER THAN VERSIONED.** A spliced transcript needs a GT..AG at every
intron or the aligner and the simulator disagree with the annotation. Deriving the sequence FROM the
annotation makes that impossible to get wrong: the chromosome is regenerated from a fixed seed and
every declared intron gets its motif injected.

⭐ **The nascent RNA needs no declaration.** `rigel index` creates a single-exon nascent ENTITY spanning
each multi-exon transcript; give it abundance through its contributor's ``nrna_abundance``.

⭐ **Strands (owner, 2026-09-02): every gene carries an explicit strand and the chromosome keeps EQUAL
representation of + and −** — a sign error in any strand-dependent rule is invisible on a one-strand
chromosome. This is not the both-stranded (overlapping, opposite-strand) locus, which is a later step;
the builder refuses a chromosome whose strand counts differ by more than one.

Usage::

    python scripts/sim/build_test_reference.py                 # render + build into the default runs dir
    python scripts/sim/build_test_reference.py --out DIR       # elsewhere
    python scripts/sim/build_test_reference.py --check         # are the checked-in renders in sync?
    python scripts/sim/build_test_reference.py --self-test     # no I/O
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rigel.sim.genome import MutableGenome  # noqa: E402
from rigel.transcript import Transcript  # noqa: E402
from rigel.types import Interval, Strand  # noqa: E402

REF_NAME = "test_chr"
#: the chromosome's random sequence is fixed by this seed, so the same annotation always gives the same
#: FASTA (the YAML may override both the seed and the length)
GENOME_SEED = 20260819
GENOME_LENGTH = 1_500_000
#: ⭐⭐ THE BLANK CHROMOSOME (owner design, 2026-08-29): a second contig with NO annotation, so `rigel index`
#: knows nothing on it — while the SIMULATOR draws unannotated "shadow" transcripts there. Its own seed;
#: its own length. gDNA is drawn on it too (`gdna.genomic_refs`).
BLANK_REF_NAME = "test_blank"
BLANK_GENOME_SEED = 20260829
BLANK_GENOME_LENGTH = 1_500_000

#: the probe geometry — the benign tiling and the two adversarial designs (owner directive 2026-09-01)
PROBE_LEN = 125
PROBES_PER_PIECE = 8
JUNCTION_BLOCKS = (62, 63)

HERE = Path(__file__).resolve().parent / "test_reference"
DEFAULT_SPEC = HERE / "test_chr.yaml"
DEFAULT_OUT = Path.home() / "Downloads" / "rigel_runs" / "test_reference"


# ── the spec ─────────────────────────────────────────────────────────────────────────────────────


@dataclass
class Gene:
    gene_id: str
    strand: str
    probed: bool
    transcripts: list[Transcript]


@dataclass
class Spec:
    genes: list[Gene]
    shadow_genes: list[Gene]
    genome_length: int = GENOME_LENGTH
    seed: int = GENOME_SEED
    blank_length: int = BLANK_GENOME_LENGTH
    blank_seed: int = BLANK_GENOME_SEED
    ref_name: str = REF_NAME
    blank_ref_name: str = BLANK_REF_NAME
    problems: list[str] = field(default_factory=list)

    @property
    def transcripts(self) -> list[Transcript]:
        return [t for g in self.genes for t in g.transcripts]

    @property
    def shadows(self) -> list[Transcript]:
        return [t for g in self.shadow_genes for t in g.transcripts]

    @property
    def abundances(self) -> dict[str, tuple[float, float]]:
        return {
            t.t_id: (float(t.abundance or 0.0), float(t.nrna_abundance or 0.0))
            for t in self.transcripts + self.shadows
        }


def _tx(t_id: str, g_id: str, ref: str, strand: str, exons, abundance: float, nrna: float) -> Transcript:
    t = Transcript(
        ref=ref,
        strand=Strand.from_str(strand),
        exons=[Interval(int(a), int(b)) for a, b in exons],
        t_id=t_id,
        g_id=g_id,
        abundance=float(abundance),
        nrna_abundance=float(nrna),
    )
    t.length = t.compute_length()
    return t


def spec_from_dict(cfg: dict) -> Spec:
    """The YAML as a :class:`Spec`. Shape problems (a gene without a strand, an exon that is not a
    pair) are collected into ``spec.problems`` rather than raised, so the builder reports every problem
    at once like it always has."""
    problems: list[str] = []

    def genes_of(key: str, ref: str, with_probed: bool) -> list[Gene]:
        out: list[Gene] = []
        for g in cfg.get(key) or []:
            gid = str(g.get("gene_id", "?"))
            strand = g.get("strand")
            if strand not in ("+", "-"):
                problems.append(f"{gid}: strand must be '+' or '-', got {strand!r}")
                strand = "+"
            probed = bool(g.get("probed", False)) if with_probed else False
            if not with_probed and "probed" in g:
                problems.append(f"{gid}: a shadow gene cannot be probed (a panel cannot target what the annotation does not know)")
            txs: list[Transcript] = []
            for td in g.get("transcripts") or []:
                tid = str(td.get("t_id", "?"))
                exons = td.get("exons") or []
                if not exons or any(not (isinstance(e, (list, tuple)) and len(e) == 2) for e in exons):
                    problems.append(f"{tid}: exons must be a non-empty list of [start, end) pairs")
                    continue
                txs.append(_tx(tid, gid, ref, strand, exons, td.get("abundance", 0.0),
                               td.get("nrna_abundance", 0.0)))
            if not txs and not problems:
                problems.append(f"{gid}: a gene needs at least one transcript")
            out.append(Gene(gid, strand, probed, txs))
        return out

    ref = str(cfg.get("ref_name", REF_NAME))
    blank = cfg.get("blank") or {}
    blank_ref = str(blank.get("ref_name", BLANK_REF_NAME))
    spec = Spec(
        genes=genes_of("genes", ref, True),
        shadow_genes=genes_of("shadow_genes", blank_ref, False),
        genome_length=int(cfg.get("genome_length", GENOME_LENGTH)),
        seed=int(cfg.get("seed", GENOME_SEED)),
        blank_length=int(blank.get("genome_length", BLANK_GENOME_LENGTH)),
        blank_seed=int(blank.get("seed", BLANK_GENOME_SEED)),
        ref_name=ref,
        blank_ref_name=blank_ref,
        problems=problems,
    )
    return spec


def load_spec(path: Path) -> Spec:
    import yaml

    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")
    cfg = yaml.safe_load(path.read_text()) or {}
    return spec_from_dict(cfg)


# ── the checks ───────────────────────────────────────────────────────────────────────────────────


def intron_spans(transcripts: list[Transcript]) -> list[tuple[int, int, Strand]]:
    """Every intron the annotation declares, as ``(start, end, strand)`` — what needs a splice motif."""
    out: list[tuple[int, int, Strand]] = []
    for t in transcripts:
        exons = sorted(t.exons, key=lambda e: e.start)
        for left, right in zip(exons, exons[1:]):
            if right.start > left.end:
                out.append((left.end, right.start, t.strand))
    return out


def inject_motifs(genome: MutableGenome, introns: list[tuple[int, int, Strand]]) -> int:
    """GT..AG at every intron, reverse-complemented on the minus strand (CT..AC in genomic terms).

    ⛔ The motif is written in GENOMIC orientation, which is what an aligner reads. Getting this
    backwards on ``−`` is invisible in the FASTA and shows up only as unaligned spliced reads.
    """
    for start, end, strand in introns:
        donor, acceptor = ("GT", "AG") if strand != Strand.NEG else ("CT", "AC")
        genome.edit(start, donor)
        genome.edit(end - 2, acceptor)
    return len(introns)


def check_transcripts(transcripts: list[Transcript], abundances: dict, shadows: list[Transcript] | None = None,
                      *, genome_length: int = GENOME_LENGTH, blank_length: int = BLANK_GENOME_LENGTH) -> list[str]:
    """Everything that would make a transcript unusable, reported together rather than one per run.

    ``transcripts`` must live on ``test_chr`` (the annotated chromosome, the one the index is built from);
    ``shadows`` must live on ``test_blank`` (the BLANK chromosome — unannotated by construction); ids are
    unique across BOTH sets (a shadow the index would know is not a shadow); every transcript of either
    set has an abundance row and every row names a transcript of one of them; a shadow carries no
    nascent abundance (a shadow IS the unannotated RNA)."""
    problems: list[str] = []
    seen: dict[str, Transcript] = {}
    shadows = shadows or []
    for t, expected_ref, length in [(t, REF_NAME, genome_length) for t in transcripts] + [
        (t, BLANK_REF_NAME, blank_length) for t in shadows
    ]:
        if t.ref != expected_ref:
            problems.append(f"{t.t_id}: reference {t.ref!r}, expected {expected_ref!r}")
        if t.t_id in seen:
            problems.append(f"{t.t_id}: declared twice")
        seen[t.t_id] = t
        exons = sorted(t.exons, key=lambda e: e.start)
        for e in exons:
            if e.start < 0 or e.end > length:
                problems.append(f"{t.t_id}: exon [{e.start}, {e.end}) outside [0, {length})")
            if e.end <= e.start:
                problems.append(f"{t.t_id}: empty exon [{e.start}, {e.end})")
        for left, right in zip(exons, exons[1:]):
            if right.start < left.end:
                problems.append(f"{t.t_id}: exons overlap at {left.end}/{right.start}")
            elif right.start - left.end < 4:
                # ⛔ a GT..AG needs 4 bases; a shorter gap cannot carry a motif
                problems.append(f"{t.t_id}: intron [{left.end}, {right.start}) is under 4 bp")
        if t.t_id not in abundances:
            problems.append(f"{t.t_id}: no row in the abundances")
    for t in shadows:
        if t.t_id in abundances and abundances[t.t_id][1] != 0:
            problems.append(f"{t.t_id}: a shadow carries nascent abundance (a shadow IS the unannotated RNA)")
    for t_id in abundances:
        if t_id not in seen:
            problems.append(f"{t_id}: has an abundance but is not in either gene list")
    return problems


def check_strand_balance(genes: list[Gene]) -> list[str]:
    """⭐ Equal representation of the two strands (owner, 2026-09-02): the counts may differ by at most
    one. A one-strand chromosome cannot see a sign error."""
    pos = sum(1 for g in genes if g.strand == "+")
    neg = len(genes) - pos
    if genes and abs(pos - neg) > 1:
        return [f"strand balance: {pos} genes on + and {neg} on − — the chromosome must carry both equally"]
    return []


def check_spec(spec: Spec) -> list[str]:
    return list(spec.problems) + check_transcripts(
        spec.transcripts, spec.abundances, spec.shadows,
        genome_length=spec.genome_length, blank_length=spec.blank_length,
    ) + check_strand_balance(spec.genes)


# ── the renders ──────────────────────────────────────────────────────────────────────────────────

_DERIVED = "DERIVED from test_chr.yaml by build_test_reference.py — do not edit; edit the YAML and re-run it."


def render_gtf(transcripts: list[Transcript], what: str) -> str:
    """One `exon` line per exon, GTF 1-based inclusive coordinates, the attribute form every consumer reads."""
    lines = [f"## {what}", f"## {_DERIVED}"]
    for t in transcripts:
        s = Strand.to_str(t.strand)
        for e in sorted(t.exons, key=lambda e: e.start):
            lines.append(f"{t.ref}\trigel_test\texon\t{e.start + 1}\t{e.end}\t.\t{s}\t.\t"
                         f'gene_id "{t.g_id}"; transcript_id "{t.t_id}";')
    return "\n".join(lines) + "\n"


def _num(x: float) -> str:
    return str(int(x)) if float(x) == int(x) else repr(float(x))


def render_abundances(spec: Spec) -> str:
    """⛔ No comment lines — both readers key on the first line."""
    lines = ["transcript_id\tmrna_abundance\tnrna_abundance"]
    for t in spec.transcripts + spec.shadows:
        lines.append(f"{t.t_id}\t{_num(t.abundance or 0.0)}\t{_num(t.nrna_abundance or 0.0)}")
    return "\n".join(lines) + "\n"


def _bed12(chrom: str, s: int, e: int, name: str, sizes, starts) -> str:
    return (f"{chrom}\t{s}\t{e}\t{name}\t0\t.\t{s}\t{e}\t0\t{len(sizes)}\t"
            f"{','.join(map(str, sizes))}\t{','.join(map(str, starts))}")


def exon_union(gene: Gene) -> list[tuple[int, int]]:
    """The gene's exonic sequence as disjoint pieces, merged across its isoforms."""
    ivs = sorted((e.start, e.end) for t in gene.transcripts for e in t.exons)
    out: list[list[int]] = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def probe_lines(genes: list[Gene], mode: str) -> list[str]:
    """The three panels, each from the probed genes alone. Coincident probes (two isoforms sharing an
    exon or an sj) are emitted once. ``mode`` is ``benign`` / ``sparse`` / ``junction``."""
    out: list[str] = []
    seen: set[tuple] = set()

    def emit(line: str) -> None:
        f = line.split("\t")
        key = (f[1], f[2], f[10], f[11])
        if key not in seen:
            seen.add(key)
            out.append(line)

    for g in genes:
        if not g.probed:
            continue
        ref = g.transcripts[0].ref
        if mode == "benign":
            for k, (s, e) in enumerate(exon_union(g), 1):
                for j in range(PROBES_PER_PIECE):
                    a = s + PROBE_LEN * j
                    if a + PROBE_LEN > e:
                        break  # a piece shorter than the full tiling gets what fits
                    emit(_bed12(ref, a, a + PROBE_LEN, f"probe_{g.gene_id}_e{k}_{j}", [PROBE_LEN], [0]))
        elif mode == "sparse":
            for t in g.transcripts:
                for k, e in enumerate(sorted(t.exons, key=lambda e: e.start)):
                    c = (e.start + e.end) // 2
                    lo, hi = c - PROBE_LEN // 2, c - PROBE_LEN // 2 + PROBE_LEN
                    emit(_bed12(ref, lo, hi, f"sparse_{t.t_id}_e{k}", [PROBE_LEN], [0]))
        elif mode == "junction":
            b1, b2 = JUNCTION_BLOCKS
            for t in g.transcripts:
                ex = sorted(t.exons, key=lambda e: e.start)
                for k, (lo, hi) in enumerate(zip(ex, ex[1:])):
                    a, b = lo.end - b1, hi.start + b2
                    emit(_bed12(ref, a, b, f"junction_{t.t_id}_sj{k}", [b1, b2], [0, hi.start - a]))
        else:
            raise ValueError(mode)
    return out


_PROBE_HEADERS = {
    "benign": ("# THE BENIGN CAPTURE PANEL — each probed gene's exon UNION tiled 8 x 125 bp seamless from each\n"
               "# union piece's start, so no probe spans an sj (a split probe suppresses exactly the\n"
               "# boundary-crossing gDNA). BED12, 0-based half-open."),
    "sparse": ("# THE SPARSE ADVERSARIAL PANEL (owner directive 2026-09-01) — one 125 bp probe CENTRED on\n"
               "# every annotated exon of a probed gene: faces depleted, mid-exon enriched. BED12."),
    "junction": ("# THE JUNCTION ADVERSARIAL PANEL (owner directive 2026-09-01) — one BED12 two-block probe\n"
                 "# (62 + 63 bp, contiguous only in cDNA) per sj of a probed gene and nothing else: spliced\n"
                 "# fragments enriched over unspliced/gDNA at the same face."),
}


def render_probes(spec: Spec, mode: str) -> str:
    return _PROBE_HEADERS[mode] + f"\n# {_DERIVED}\n" + "\n".join(probe_lines(spec.genes, mode)) + "\n"


def render_all(spec: Spec) -> dict[str, str]:
    """Every rendered file, by name."""
    return {
        "test_chr.gtf": render_gtf(spec.transcripts, "THE TEST CHROMOSOME'S ANNOTATION — what `rigel index` reads"),
        "test_shadow.gtf": render_gtf(spec.shadows, "SHADOW TRANSCRIPTS on the blank contig — what the SIMULATOR reads and the index never sees"),
        "test_abundances.tsv": render_abundances(spec),
        "test_probes.bed": render_probes(spec, "benign"),
        "test_probes_sparse.bed": render_probes(spec, "sparse"),
        "test_probes_junction.bed": render_probes(spec, "junction"),
    }


def check_renders(spec: Spec, render_dir: Path = HERE) -> list[str]:
    """Which checked-in renders differ from what the YAML renders now."""
    drift = []
    for name, text in render_all(spec).items():
        p = render_dir / name
        if not p.is_file() or p.read_text() != text:
            drift.append(name)
    return drift


# ── the build ────────────────────────────────────────────────────────────────────────────────────


def _write_two_record_fasta(genomes: list[MutableGenome], path: Path) -> Path:
    """One FASTA, two records, 80-column wrap, + ``samtools faidx``. ⚠ The file keeps the name
    ``test_chr.fa`` so every panel config's ``genome:`` path stays valid; it carries BOTH contigs."""
    import pysam

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for g in genomes:
            f.write(f">{g.name}\n")
            seq = g.seq
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    for stale in (path.with_suffix(path.suffix + ".fai"),):
        if stale.exists():
            stale.unlink()
    pysam.faidx(str(path))
    return path


def build(spec_path: Path, out: Path, render_dir: Path | None = HERE) -> dict:
    """Render every derived file (into the repo directory when ``render_dir`` is given, and into
    ``out``), write the derived FASTA (both contigs), and report what the benchmark holds."""
    spec = load_spec(spec_path)
    problems = check_spec(spec)
    if problems:
        raise ValueError(
            f"{len(problems)} problem(s) in the test reference:\n  " + "\n  ".join(problems)
        )
    renders = render_all(spec)
    out.mkdir(parents=True, exist_ok=True)
    for d in ([render_dir] if render_dir is not None else []) + [out]:
        for name, text in renders.items():
            (Path(d) / name).write_text(text)
    genome = MutableGenome(spec.genome_length, seed=spec.seed, name=spec.ref_name)
    introns = intron_spans(spec.transcripts)
    inject_motifs(genome, introns)
    blank = MutableGenome(spec.blank_length, seed=spec.blank_seed, name=spec.blank_ref_name)
    shadow_introns = intron_spans(spec.shadows)
    inject_motifs(blank, shadow_introns)
    fasta = _write_two_record_fasta([genome, blank], out / f"{spec.ref_name}.fa")
    multi = sum(1 for t in spec.transcripts if len(t.exons) > 1)
    pos = sum(1 for g in spec.genes if g.strand == "+")
    return {"fasta": fasta, "n_genes": len(spec.genes), "n_pos": pos, "n_neg": len(spec.genes) - pos,
            "n_probed": sum(1 for g in spec.genes if g.probed),
            "n_transcripts": len(spec.transcripts), "n_multi_exon": multi,
            "n_introns": len(introns), "length": spec.genome_length,
            "n_nascent_entities_expected": multi,
            "n_shadow_transcripts": len(spec.shadows), "n_shadow_introns": len(shadow_introns),
            "blank_ref": spec.blank_ref_name, "blank_length": spec.blank_length,
            "renders": list(renders)}


# ── the self-test ────────────────────────────────────────────────────────────────────────────────


def self_test() -> int:
    """⛔ Every check perturbed, no I/O against the real reference."""
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    def tx(t_id, exons, strand=Strand.POS, ref=REF_NAME, ab=1.0, nr=0.0):
        return _tx(t_id, "g_" + t_id, ref, Strand.to_str(strand), exons, ab, nr)

    # ── the genome is reproducible, and the seed is what makes it so
    a = MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq
    b = MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq
    check("the chromosome is reproducible from its seed", a == b)
    check("a different seed gives a different chromosome",
          MutableGenome(2000, seed=GENOME_SEED + 1, name=REF_NAME).seq != a)

    # ── motifs: written, in genomic orientation, and only where the annotation says
    g = MutableGenome(2000, seed=1, name=REF_NAME)
    n = inject_motifs(g, intron_spans([tx("T", [(100, 200), (400, 500)])]))
    check("one intron gives one motif pair", n == 1)
    check("a + intron gets GT..AG", g[200:202] == "GT" and g[398:400] == "AG")
    gm = MutableGenome(2000, seed=1, name=REF_NAME)
    inject_motifs(gm, intron_spans([tx("T", [(100, 200), (400, 500)], strand=Strand.NEG)]))
    check("a − intron gets CT..AC in genomic orientation",
          gm[200:202] == "CT" and gm[398:400] == "AC")
    check("a single-exon transcript declares no intron", intron_spans([tx("S", [(10, 900)])]) == [])

    # ── the checks fire on what would break a build, and pass on what would not
    good = [tx("T1", [(1000, 2000), (5000, 6000)])]
    ab = {"T1": (100.0, 25.0)}
    check("a well-formed transcript with an abundance passes", check_transcripts(good, ab) == [])
    check("an exon past the end of the chromosome is caught",
          any("outside" in p for p in check_transcripts([tx("X", [(0, GENOME_LENGTH + 1)])], {"X": (1.0, 0.0)})))
    check("the YAML's genome length is what bounds an exon",
          any("outside" in p for p in check_transcripts([tx("X", [(0, 3000)])], {"X": (1.0, 0.0)}, genome_length=2000)))
    check("overlapping exons are caught",
          any("overlap" in p for p in check_transcripts([tx("X", [(100, 500), (400, 900)])], {"X": (1.0, 0.0)})))
    check("an intron too short to carry a motif is caught",
          any("under 4 bp" in p for p in check_transcripts([tx("X", [(100, 200), (202, 400)])], {"X": (1.0, 0.0)})))
    check("a transcript with no abundance row is caught",
          any("no row" in p for p in check_transcripts(good, {})))
    check("an abundance for a transcript that is in neither gene list is caught",
          any("not in either" in p for p in check_transcripts(good, {**ab, "GHOST": (1.0, 0.0)})))
    check("a duplicate transcript id is caught",
          any("twice" in p for p in check_transcripts(good + good, ab)))
    check("the wrong reference name is caught",
          any("expected" in p for p in check_transcripts([tx("X", [(10, 90)], ref="chr1")], {"X": (1.0, 0.0)})))

    # ── the BLANK chromosome and its shadows (2026-08-29)
    sh = [tx("shadow_A", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)]
    both = {"T1": (100.0, 25.0), "shadow_A": (5.0, 0.0)}
    check("a shadow on the blank chromosome with an abundance row passes", check_transcripts(good, both, sh) == [])
    check("a shadow placed on the ANNOTATED chromosome is caught",
          any("expected 'test_blank'" in p for p in check_transcripts(good, both, [tx("shadow_A", [(1000, 2000), (5000, 6000)])])))
    check("an annotated transcript placed on the blank chromosome is caught",
          any("expected 'test_chr'" in p for p in check_transcripts([tx("T1", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)], both, sh)))
    check("a shadow id the annotation already declares is caught (it would not be a shadow)",
          any("twice" in p for p in check_transcripts(good, both, [tx("T1", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)])))
    check("a shadow with no abundance row is caught", any("no row" in p for p in check_transcripts(good, ab, sh)))
    check("a shadow with nascent abundance is caught",
          any("shadow IS" in p for p in check_transcripts(good, {**both, "shadow_A": (5.0, 1.0)}, sh)))
    gb = MutableGenome(8000, seed=BLANK_GENOME_SEED, name=BLANK_REF_NAME)
    inject_motifs(gb, intron_spans(sh))
    check("a shadow intron gets its motif on the BLANK genome", gb[2000:2002] == "GT" and gb[4998:5000] == "AG")
    check("the two contigs are distinct sequences",
          MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq != MutableGenome(2000, seed=BLANK_GENOME_SEED, name=BLANK_REF_NAME).seq)

    # ── zero transcripts is the STARTING state and must be legal
    check("zero transcripts is legal", check_transcripts([], {}) == [])
    check("zero transcripts declares no intron", intron_spans([]) == [])

    # ── the YAML spec (2026-09-02): shape problems are collected, not raised
    cfg = {"genome_length": 50_000, "genes": [
        {"gene_id": "gA", "strand": "+", "probed": True,
         "transcripts": [{"t_id": "A", "exons": [[1000, 2000], [9000, 10000]], "abundance": 10, "nrna_abundance": 3}]},
        {"gene_id": "gB", "strand": "-", "transcripts": [{"t_id": "B", "exons": [[20000, 21000]], "abundance": 5}]},
    ], "shadow_genes": [
        {"gene_id": "gS", "strand": "-", "transcripts": [{"t_id": "S", "exons": [[100, 600]], "abundance": 2}]},
    ]}
    spec = spec_from_dict(cfg)
    check("a well-formed spec has no problems", check_spec(spec) == [])
    check("the spec's abundances carry mature AND nascent", spec.abundances["A"] == (10.0, 3.0))
    check("a shadow lands on the blank contig", spec.shadows[0].ref == BLANK_REF_NAME)
    bad = spec_from_dict({**cfg, "genes": [{**cfg["genes"][0], "strand": "*"}] + cfg["genes"][1:]})
    check("a gene without a legal strand is caught", any("strand must be" in p for p in bad.problems))
    bad = spec_from_dict({**cfg, "genes": [{**cfg["genes"][0], "transcripts": [{"t_id": "A", "exons": [[1, 2, 3]]}]}]})
    check("a malformed exon pair is caught", any("pairs" in p for p in bad.problems))
    bad = spec_from_dict({**cfg, "shadow_genes": [{**cfg["shadow_genes"][0], "probed": True}]})
    check("a probed shadow is caught", any("shadow gene cannot be probed" in p for p in bad.problems))

    # ── strand balance (owner, 2026-09-02): equal representation, the counts differ by at most one
    gpos = Gene("g", "+", False, [tx("P", [(10, 90)])])
    gneg = Gene("g", "-", False, [tx("N", [(10, 90)])])
    check("balanced strands pass", check_strand_balance([gpos, gneg, gpos]) == [])
    check("a one-strand chromosome is refused", check_strand_balance([gpos, gpos, gpos]) != [])
    check("an empty chromosome is balanced", check_strand_balance([]) == [])

    # ── the renders
    gtf = render_gtf(spec.transcripts, "x")
    check("the GTF is one exon line per exon, 1-based inclusive",
          "test_chr\trigel_test\texon\t1001\t2000\t.\t+\t.\tgene_id \"gA\"; transcript_id \"A\";" in gtf
          and gtf.count("\texon\t") == 3)
    check("a − exon renders with its strand", "\t20001\t21000\t.\t-\t" in gtf)
    tsv = render_abundances(spec)
    check("the abundances render with no comment lines and the header first",
          tsv.splitlines()[0] == "transcript_id\tmrna_abundance\tnrna_abundance" and "A\t10\t3" in tsv and "S\t2\t0" in tsv)
    ben = probe_lines(spec.genes, "benign")
    check("the benign panel tiles only the probed gene, 8 per exon piece", len(ben) == 16 and all("gA" in line for line in ben))
    check("the benign tiling starts at the exon start and never spans the sj",
          ben[0].split("\t")[1:3] == ["1000", "1125"] and ben[7].split("\t")[1:3] == ["1875", "2000"])
    spa = probe_lines(spec.genes, "sparse")
    check("the sparse panel centres one probe per exon", [line.split("\t")[1:3] for line in spa] == [["1438", "1563"], ["9438", "9563"]])
    jun = probe_lines(spec.genes, "junction")
    check("the junction panel is one two-block probe per sj",
          len(jun) == 1 and jun[0].split("\t")[9:12] == ["2", "62,63", "0,7062"] and jun[0].split("\t")[1:3] == ["1938", "9063"])
    # two isoforms sharing exons: the union tiles once, coincident probes appear once
    iso = spec_from_dict({"genome_length": 50_000, "genes": [
        {"gene_id": "gI", "strand": "+", "probed": True, "transcripts": [
            {"t_id": "U", "exons": [[1000, 2000], [9000, 10000]], "abundance": 1},
            {"t_id": "T", "exons": [[500, 2000], [9000, 10000]], "abundance": 1}]},
        {"gene_id": "gJ", "strand": "-", "transcripts": [{"t_id": "J", "exons": [[30000, 31000]], "abundance": 1}]}]})
    check("the exon union merges overlapping isoform exons", exon_union(iso.genes[0]) == [(500, 2000), (9000, 10000)])
    ben2 = probe_lines(iso.genes, "benign")
    check("the benign tiling follows the union, 8 per piece, from the union start",
          len(ben2) == 16 and ben2[0].split("\t")[1] == "500")
    check("a shared sj yields ONE junction probe", len(probe_lines(iso.genes, "junction")) == 1)
    check("the sparse panel keeps each isoform's own exon centre", len(probe_lines(iso.genes, "sparse")) == 3)
    check("an unprobed gene contributes to no panel",
          not any("gJ" in line or "_J_" in line for m in ("benign", "sparse", "junction") for line in probe_lines(iso.genes, m)))
    # the check refuses drift
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        for name, text in render_all(spec).items():
            (d / name).write_text(text)
        check("in-sync renders pass the check", check_renders(spec, d) == [])
        (d / "test_chr.gtf").write_text("edited by hand\n")
        check("a hand-edited render is caught by the check", check_renders(spec, d) == ["test_chr.gtf"])

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", type=Path, default=DEFAULT_SPEC, help="the ONE hand-edited file")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--no-repo-render", action="store_true",
                    help="do not refresh the versioned renders beside the YAML (only the runs dir)")
    ap.add_argument("--check", action="store_true",
                    help="render in memory and refuse if a checked-in render has drifted from the YAML")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if args.check:
        spec = load_spec(args.spec)
        problems = check_spec(spec)
        drift = check_renders(spec, args.spec.parent)
        for p in problems:
            print(f"   ⛔ {p}")
        for name in drift:
            print(f"   ⛔ {name} differs from what {args.spec.name} renders — re-run this script")
        if not problems and not drift:
            print(f"   ✔ every render in {args.spec.parent} matches {args.spec.name}")
        return 1 if (problems or drift) else 0

    info = build(args.spec, args.out, None if args.no_repo_render else args.spec.parent)
    print(f"\n⭐ THE TEST REFERENCE — {REF_NAME}, {info['length']:,} bp, from {args.spec}")
    print(f"   fasta        {info['fasta']}")
    print(f"   renders      {', '.join(info['renders'])}  (repo dir + {args.out})")
    print(f"   genes        {info['n_genes']}  (+ {info['n_pos']} / − {info['n_neg']}; {info['n_probed']} probed)")
    print(f"   transcripts  {info['n_transcripts']}  ({info['n_multi_exon']} multi-exon)")
    print(f"   introns      {info['n_introns']} (each with a GT..AG injected)")
    print(f"   ⭐ nascent entities `rigel index` will create: {info['n_nascent_entities_expected']}"
          " — one single-exon transcript spanning each multi-exon one")
    print(f"   ⭐⭐ BLANK chromosome {info['blank_ref']}, {info['blank_length']:,} bp — NO annotation; "
          f"{info['n_shadow_transcripts']} SHADOW transcripts ({info['n_shadow_introns']} shadow introns with motifs), "
          "never given to the index")
    if info["n_transcripts"] == 0:
        print(f"\n   ⭐ ZERO TRANSCRIPTS — the benchmark's starting state. Add a gene to {args.spec.name} and re-run this.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
