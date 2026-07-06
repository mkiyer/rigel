#!/usr/bin/env python
"""Generate an AMBIG-dense synthetic reference (genome.fa + genes.gtf) to stress the
2-D calibration solver. Drops into <suite>/reference/ for `simulate_suite.py --skip-existing`.

Design (per maintainer spec 2026-07-06):
  - ~100 single-stranded ANCHOR genes (prior anchor); mix of single-exon (1 isoform) + multi-exon.
  - >=200 AMBIG genes across many opposite-strand overlap TOPOLOGIES, tagged for later dissection.
  - single-exon transcripts: exactly ONE isoform.  NO exact-duplicate transcripts (global check).
  - a splice-junction POSITION supports only ONE strand (no +/- junction shares a coordinate).
  - long/short/MICRO exon mix, incl. exons < RNA fragment length (200) -> no contained frags -> relay.
  - genome ~10 Mb, scaled so exonic < 5% (rest intergenic+intronic).
Emits class_tags.tsv (transcript_id, gene_id, topology, strand, role) for dissection grouping.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

from rigel.sim.synthetic_genome import (
    GeneDef, TranscriptDef, ExonDef, REF_NAME,
    generate_random_genome, inject_splice_sites, write_fasta, write_gtf,
)

SEED = 20260706
FL_RNA = 200          # RNA fragment length (defines the "short exon" threshold)
N_ANCHOR = 100
N_AMBIG_LOCI = 115    # ~2 genes/locus -> >=200 AMBIG genes (multi_partner adds a 3rd)
MIN_INTRON = 800
MAX_INTRON = 14000
INTERGENIC_LO = 6000
INTERGENIC_HI = 28000

# ── exon-length mixture: microexon / short(<FL) / medium / long ──────────────
def sample_exon_len(rng) -> int:
    r = rng.random()
    if r < 0.12:  return int(rng.integers(20, 60))      # microexon: pure relay (< FL, tiny)
    if r < 0.38:  return int(rng.integers(60, FL_RNA))   # short: no contained fragment at FL=200
    if r < 0.78:  return int(rng.integers(FL_RNA, 700))  # medium
    return int(rng.integers(700, 2500))                  # long


def build_exons(rng, start: int, n_exons: int) -> tuple[list[ExonDef], int]:
    """Place n_exons left->right from `start`; return (exons, end_exclusive)."""
    exons = []
    pos = start
    for i in range(n_exons):
        el = sample_exon_len(rng)
        exons.append(ExonDef(pos, pos + el))
        pos += el
        if i < n_exons - 1:
            pos += int(rng.integers(MIN_INTRON, MAX_INTRON))
    return exons, pos


def _sig(exons: list[ExonDef]) -> tuple:
    return tuple((e.start, e.end) for e in exons)


class Builder:
    def __init__(self, seed=SEED):
        self.rng = np.random.default_rng(seed)
        self.genes: list[GeneDef] = []
        self.tags: list[tuple] = []          # (t_id, gene_id, topology, strand, role)
        self._sigs: set = set()              # global (strand, exon-tuple) dedup
        self._junc: dict = {}                # junction coord -> strand (one-strand-per-position)
        self._gid = 0

    # ---- registration with global constraints ----
    def _register_tx(self, tx: TranscriptDef, topology: str, role: str) -> bool:
        key = (tx.strand, _sig(tx.exons))
        if key in self._sigs:
            return False  # exact duplicate -> reject
        # one-strand-per-junction: every intron boundary coord must not be claimed by other strand
        jcoords = []
        for i in range(len(tx.exons) - 1):
            for c in (tx.exons[i].end, tx.exons[i + 1].start):
                if self._junc.get(c, tx.strand) != tx.strand:
                    return False  # collision with opposite strand -> reject
                jcoords.append(c)
        self._sigs.add(key)
        for c in jcoords:
            self._junc[c] = tx.strand
        self.tags.append((tx.t_id, tx.gene_id, topology, tx.strand, role))
        return True

    def _new_gene(self, strand: str, name_suffix: str) -> tuple[str, str]:
        self._gid += 1
        gid = f"G{self._gid:04d}"
        return gid, f"{gid}_{name_suffix}"

    def _isoforms(self, rng, base: list[ExonDef], gid, gname, strand, topology, role,
                  max_iso=3) -> list[TranscriptDef]:
        """Single-exon -> exactly 1 isoform. Multi-exon -> 1..max_iso distinct isoforms
        via exon-skip / alt-first / alt-last, all globally deduped."""
        txs = []
        if len(base) == 1:
            tx = TranscriptDef(f"{gid}.1", gid, gname, strand, [ExonDef(base[0].start, base[0].end)])
            if self._register_tx(tx, topology, role):
                txs.append(tx)
            return txs
        n_iso = int(rng.integers(1, max_iso + 1))
        variants = [base]
        for _ in range(n_iso * 3):
            if len(variants) >= n_iso:
                break
            v = [ExonDef(e.start, e.end) for e in base]
            mode = rng.integers(0, 3)
            if mode == 0 and len(v) > 2:          # exon skip (interior)
                j = int(rng.integers(1, len(v) - 1)); del v[j]
            elif mode == 1 and len(v) >= 2:        # alt first exon (truncate 5')
                v = v[1:]
            elif len(v) >= 2:                      # alt last exon (truncate 3')
                v = v[:-1]
            if _sig(v) not in {_sig(x) for x in variants}:
                variants.append(v)
        for k, ex in enumerate(variants[:n_iso], 1):
            tx = TranscriptDef(f"{gid}.{k}", gid, gname, strand, ex)
            if self._register_tx(tx, topology, role):
                txs.append(tx)
        return txs

    def add_gene(self, strand, exons_list_or_base, topology, role, name_suffix, multi_iso=True):
        gid, gname = self._new_gene(strand, name_suffix)
        base = exons_list_or_base
        txs = self._isoforms(self.rng, base, gid, gname, strand, topology, role,
                             max_iso=(3 if multi_iso else 1))
        if not txs:
            return None
        g = GeneDef(gid, gname, strand,
                    min(t.start for t in txs), max(t.end for t in txs), txs)
        self.genes.append(g)
        return g

    # ─────────────────────────────────────────────────────────────────────
    # ANCHOR single-stranded gene
    def anchor(self, start: int) -> int:
        rng = self.rng
        strand = "+" if rng.random() < 0.5 else "-"
        single = rng.random() < 0.35
        n_exons = 1 if single else int(rng.integers(2, 7))
        exons, end = build_exons(rng, start, n_exons)
        self.add_gene(strand, exons, "anchor", "anchor", "anchor", multi_iso=not single)
        return end

    # ─────────────────────────────────────────────────────────────────────
    # AMBIG topologies. Each returns the locus end coordinate.
    def _ambig_exon_exon(self, start):
        """EX+ INTERSECT EX-: an exon of A(+) overlaps an exon of B(-)."""
        rng = self.rng
        aex, aend = build_exons(rng, start, int(rng.integers(2, 5)))
        # pick an interior/last exon of A to be the overlap host
        hi = rng.integers(0, len(aex))
        host = aex[hi]
        ov_len = int(rng.integers(120, max(140, min(600, host.end - host.start))))
        ov_start = host.start + int(rng.integers(0, max(1, (host.end - host.start) - ov_len)))
        # B(-) exon overlaps [ov_start, ov_start+ov_len], offset so junctions differ
        bex = [ExonDef(ov_start + 30, ov_start + ov_len + 30)]
        if rng.random() < 0.6:  # give B a second exon (a junction) placed clear of A
            bex = [ExonDef(ov_start + 30, ov_start + ov_len + 30),
                   ExonDef(aend + MIN_INTRON + 500, aend + MIN_INTRON + 500 + sample_exon_len(rng))]
        self.add_gene("+", aex, "exon_exon", "sense", "eeA")
        end2 = max(e.end for e in bex)
        self.add_gene("-", bex, "exon_exon", "antisense", "eeB")
        return max(aend, end2)

    def _ambig_exon_over_intron(self, start):
        """EX(-) inside IN(+): the region-236 class. A(+) has a long intron holding B(-)'s exon(s)."""
        rng = self.rng
        e1 = ExonDef(start, start + sample_exon_len(rng))
        intron_len = int(rng.integers(4000, MAX_INTRON))
        i0 = e1.end
        i1 = i0 + intron_len
        e2 = ExonDef(i1, i1 + sample_exon_len(rng))
        aex = [e1, e2]
        if rng.random() < 0.5:  # 3-exon A (two introns; B sits in the first)
            e3s = e2.end + int(rng.integers(MIN_INTRON, MAX_INTRON))
            aex.append(ExonDef(e3s, e3s + sample_exon_len(rng)))
        # B(-) exon(s) fully inside A's first intron (i0, i1), offset from i0/i1
        bstart = i0 + int(rng.integers(300, 800))
        nB = int(rng.integers(1, 3))
        bex, bend = build_exons(rng, bstart, nB)
        if bend >= i1 - 200:  # keep B inside the intron
            bex = [ExonDef(bstart, min(bstart + sample_exon_len(rng), i1 - 300))]
        self.add_gene("+", aex, "exon_over_intron", "sense_intron_host", "eiA")
        self.add_gene("-", bex, "exon_over_intron", "antisense_exon", "eiB")
        return max(e.end for e in aex)

    def _ambig_intron_intron(self, start):
        """IN(+) INTERSECT IN(-): overlapping introns of both strands (pure nascent-vs-nascent 2-D)."""
        rng = self.rng
        # A(+): exon - long intron - exon
        e1 = ExonDef(start, start + sample_exon_len(rng))
        iA = int(rng.integers(6000, MAX_INTRON))
        e2 = ExonDef(e1.end + iA, e1.end + iA + sample_exon_len(rng))
        aex = [e1, e2]
        # B(-): exon - long intron - exon, its intron overlapping A's intron, exons offset & clear
        b1 = ExonDef(e1.end + 400, e1.end + 400 + int(rng.integers(60, 200)))
        iB = int(rng.integers(6000, MAX_INTRON))
        b2s = min(b1.end + iB, e2.start - 400)
        b2 = ExonDef(b2s, b2s + int(rng.integers(60, 300)))
        bex = [b1, b2]
        self.add_gene("+", aex, "intron_intron", "sense", "iiA")
        self.add_gene("-", bex, "intron_intron", "antisense", "iiB")
        return max(e2.end, b2.end)

    def _ambig_nested(self, start):
        """B(-) multi-exon gene nested entirely within an intron of A(+)."""
        rng = self.rng
        e1 = ExonDef(start, start + sample_exon_len(rng))
        iA = int(rng.integers(9000, MAX_INTRON))
        e2 = ExonDef(e1.end + iA, e1.end + iA + sample_exon_len(rng))
        aex = [e1, e2]
        # nested B(-): 2-3 exons inside (e1.end, e2.start)
        bstart = e1.end + int(rng.integers(600, 1500))
        bex, bend = build_exons(rng, bstart, int(rng.integers(2, 4)))
        if bend >= e2.start - 400:
            scale = (e2.start - 400 - bstart) / max(1, bend - bstart)
            bex = [ExonDef(bstart, bstart + max(40, int((e.end - e.start) * scale))) for e in bex[:1]]
        self.add_gene("+", aex, "nested", "sense_host", "nsA")
        self.add_gene("-", bex, "nested", "antisense_nested", "nsB")
        return max(e.end for e in aex)

    def _ambig_staggered(self, start):
        """Partial/staggered overlap: A(+) and B(-) offset, overlapping in the middle."""
        rng = self.rng
        aex, aend = build_exons(rng, start, int(rng.integers(2, 5)))
        mid = (start + aend) // 2
        bex, bend = build_exons(rng, mid + int(rng.integers(200, 1200)), int(rng.integers(2, 5)))
        self.add_gene("+", aex, "staggered", "sense", "stA")
        self.add_gene("-", bex, "staggered", "antisense", "stB")
        return max(aend, bend)

    def _ambig_multi_partner(self, start):
        """A(+) long gene overlapped by TWO -strand genes at different sub-regions."""
        rng = self.rng
        aex, aend = build_exons(rng, start, int(rng.integers(4, 7)))
        self.add_gene("+", aex, "multi_partner", "sense_host", "mpA")
        # two - partners over different A exons
        for k, hi in enumerate(sorted(rng.choice(len(aex), size=min(2, len(aex)), replace=False))):
            host = aex[hi]
            bs = host.start + 25
            be = min(host.end + int(rng.integers(0, 400)), aend)
            if be - bs < 80:
                be = bs + 100
            self.add_gene("-", [ExonDef(bs, be)], "multi_partner", f"antisense_{k+1}", f"mpB{k+1}")
        return aend

    def ambig_locus(self, start: int, which: int) -> int:
        fns = [self._ambig_exon_exon, self._ambig_exon_over_intron, self._ambig_intron_intron,
               self._ambig_nested, self._ambig_staggered, self._ambig_multi_partner]
        return fns[which % len(fns)](start)

    # ─────────────────────────────────────────────────────────────────────
    def layout(self):
        rng = self.rng
        pos = int(rng.integers(INTERGENIC_LO, INTERGENIC_HI))
        # interleave anchors and ambig loci in a shuffled order
        slots = ["anchor"] * N_ANCHOR + ["ambig"] * N_AMBIG_LOCI
        rng.shuffle(slots)
        ai = 0
        for s in slots:
            if s == "anchor":
                pos = self.anchor(pos)
            else:
                pos = self.ambig_locus(pos, ai); ai += 1
            pos += int(rng.integers(INTERGENIC_LO, INTERGENIC_HI))
        return pos + int(rng.integers(INTERGENIC_LO, INTERGENIC_HI))


def main(outdir: str):
    ref_dir = Path(outdir) / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    b = Builder()
    genome_end = b.layout()

    # exonic fraction -> scale genome length so exonic < 5%
    exon_bp = sum(e.end - e.start for g in b.genes for tx in g.transcripts for e in tx.exons)
    # (double-counts shared exons across isoforms; use unique intervals for a fair fraction)
    uniq = set()
    for g in b.genes:
        for tx in g.transcripts:
            for e in tx.exons:
                uniq.add((e.start, e.end))
    uniq_exon_bp = sum(e - s for s, e in uniq)
    genome_length = max(genome_end, int(uniq_exon_bp / 0.045))  # pad so exonic <= 4.5%

    all_tx = [tx for g in b.genes for tx in g.transcripts]
    n_ambig_genes = sum(1 for g in b.genes if any(t[1] == g.gene_id and t[2] != "anchor" for t in b.tags))
    n_anchor = sum(1 for g in b.genes if g.transcripts and b._tag_of(g) == "anchor") if hasattr(b, "_tag_of") else \
               sum(1 for g in b.genes if any(t[1] == g.gene_id and t[2] == "anchor" for t in b.tags))

    # genome + motifs
    seq = generate_random_genome(genome_length, SEED)
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
    inject_splice_sites(arr, all_tx)
    seq = arr.tobytes().decode("ascii")
    write_fasta(seq, REF_NAME, ref_dir)
    write_gtf(b.genes, REF_NAME, ref_dir)

    # class tags sidecar
    with open(ref_dir / "class_tags.tsv", "w") as f:
        f.write("transcript_id\tgene_id\ttopology\tstrand\trole\n")
        for t in b.tags:
            f.write("\t".join(map(str, t)) + "\n")

    # exon-length histogram for the report
    elens = np.array([e.end - e.start for _, e in [((), ExonDef(s, e)) for s, e in uniq]])
    micro = int((elens < 60).sum()); short = int(((elens >= 60) & (elens < FL_RNA)).sum())
    med = int(((elens >= FL_RNA) & (elens < 700)).sum()); lng = int((elens >= 700).sum())

    from collections import Counter
    topo = Counter(t[2] for t in b.tags)
    print(f"genome_length      : {genome_length:,} bp")
    print(f"unique exonic bp   : {uniq_exon_bp:,} ({100*uniq_exon_bp/genome_length:.2f}% of genome)")
    print(f"genes              : {len(b.genes)}  ({n_anchor} anchor, {len(b.genes)-n_anchor} AMBIG)")
    print(f"transcripts        : {len(all_tx)}  (multi-exon {sum(1 for t in all_tx if len(t.exons)>1)})")
    print(f"exon len (unique)  : micro<60={micro}  short<{FL_RNA}={short}  med<700={med}  long>=700={lng}")
    print(f"topology counts    : {dict(topo)}")
    print(f"single-exon genes w/ >1 iso (must be 0): "
          f"{sum(1 for g in b.genes if len(g.transcripts)>1 and all(len(t.exons)==1 for t in g.transcripts))}")
    print(f"junction coords claimed: {len(b._junc)}  (one-strand-per-position enforced)")
    print(f"reference written  : {ref_dir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "Downloads/rigel_runs/ambig_dense_10mb"))
