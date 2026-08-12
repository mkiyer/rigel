"""Unit tests for the net fragment-flow deconvolution math (rigel.sim.net_flow).

Uses the canonical worked example: one locus with isoforms T1/T2/T3 + a gDNA component,
truth = 10 fragments each, and the tool observes T1=12, T2=2, T3=8, gDNA=18. The net-flow
reduction must (a) satisfy observed-expected = Σ net inflow, (b) sum to zero over the locus,
and (c) decompose each transcript's surplus/deficit into gDNA-source vs RNA-isoform-source.
"""

from rigel.sim.net_flow import FlowData, _flow_marginals, _net_flow_rows


def _example() -> FlowData:
    # components: 0=T1, 1=T2, 2=T3 (rna), 3=gdna@0 — all in locus 0.
    # rows=true, cols=assigned; row sums (expected) = 10 each (truth),
    # col sums (observed) = T1:12, T2:2, T3:8, gDNA:18.
    flow = {
        (0, 0): 8,
        (0, 3): 2,  # T1 true: 8 stay, 2 -> gDNA
        (1, 0): 4,
        (1, 1): 2,
        (1, 2): 2,
        (1, 3): 2,  # T2 true: leaks out widely
        (2, 2): 6,
        (2, 3): 4,  # T3 true: 6 stay, 4 -> gDNA
        (3, 3): 10,  # gDNA true: all correct
    }
    return FlowData(
        condition="unit",
        flow=flow,
        comp_name={0: "T1", 1: "T2", 2: "T3", 3: "gdna@0"},
        comp_kind={0: "rna", 1: "rna", 2: "rna", 3: "gdna"},
        comp_locus={0: 0, 1: 0, 2: 0, 3: 0},
        total_gdna_true=10,
    )


def test_marginals_match_truth_and_observed():
    fd = _example()
    expected, observed = _flow_marginals(fd.flow)
    assert [expected[c] for c in (0, 1, 2, 3)] == [10, 10, 10, 10]
    assert [observed[c] for c in (0, 1, 2, 3)] == [12, 2, 8, 18]


def test_net_inflow_identity():
    """observed[c] - expected[c] == Σ_a net(a→c) for every component."""
    fd = _example()
    expected, observed = _flow_marginals(fd.flow)

    def net(a, b):
        return fd.flow.get((a, b), 0) - fd.flow.get((b, a), 0)

    for c in (0, 1, 2, 3):
        net_inflow = sum(net(a, c) for a in (0, 1, 2, 3) if a != c)
        assert observed[c] - expected[c] == net_inflow


def test_per_transcript_decomposition():
    fd = _example()
    tx_rows, locus_rows = _net_flow_rows(fd)
    by_tx = {r["transcript_id"]: r for r in tx_rows}

    # Δ = net_from_gdna + net_from_rna_isoforms + cross_locus, cross_locus == 0 here.
    for r in tx_rows:
        assert r["delta"] == r["net_from_gdna"] + r["net_from_rna_isoforms"] + r["cross_locus"]
        assert r["cross_locus"] == 0

    assert (
        by_tx["T1"]["delta"],
        by_tx["T1"]["net_from_gdna"],
        by_tx["T1"]["net_from_rna_isoforms"],
    ) == (2, -2, 4)
    assert (
        by_tx["T2"]["delta"],
        by_tx["T2"]["net_from_gdna"],
        by_tx["T2"]["net_from_rna_isoforms"],
    ) == (-8, -2, -6)
    assert (
        by_tx["T3"]["delta"],
        by_tx["T3"]["net_from_gdna"],
        by_tx["T3"]["net_from_rna_isoforms"],
    ) == (-2, -4, 2)

    # Conservation: Σ transcript Δ + gDNA Δ == 0 over the closed locus.
    locus = locus_rows[0]
    gdna_delta = locus["gdna_observed"] - locus["gdna_expected"]
    assert sum(r["delta"] for r in tx_rows) + gdna_delta == 0


def test_locus_row_net_gdna_to_rna():
    fd = _example()
    _, locus_rows = _net_flow_rows(fd)
    locus = locus_rows[0]
    assert locus["gdna_expected"] == 10
    assert locus["gdna_observed"] == 18
    # gDNA gained 8 net from RNA ⇒ RNA siphoned into gDNA ⇒ net_gdna_to_rna negative.
    assert locus["net_gdna_to_rna"] == -8
    assert locus["rna_delta"] == -8
    assert locus["n_transcripts"] == 3
