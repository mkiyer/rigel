"""Tests for ``scripts/design/calibration_oracle.py``'s gate on gDNA outside the genomic references.

Which references are genomic is an input to the simulator (the manifest's ``gdna.genomic_refs``), never an
inference from where gDNA landed: an annotated reference can be an RNA-only spike-in. A gate that derives the set
from the realized gDNA counts certifies whatever it is shown, so these tests plant gDNA where the manifest says
none can be.
"""

from __future__ import annotations

import numpy as np
import pytest

import calibration_oracle as co

#: two references: a genome and an RNA-only spike-in backbone
REF_IDS = {"chr1": 0, "ERCC-00002": 1}


def _manifest(genomic_refs, n_gdna=1000, name="c"):
    return {
        "gdna": {"genomic_refs": genomic_refs},
        "conditions": [{"name": name, "n_gdna": n_gdna}],
    }


def _outside_gate(verdicts):
    (gate,) = [v for v in verdicts if v["gate"] == "gate exact-zeros:gdna outside genomic refs"]
    return gate


def _slots():
    """Four slots: two on the genome, two on the spike-in."""
    ref = np.array([0, 0, 1, 1])
    eff_g = np.full(4, 100.0)
    label = np.array(["R intron", "R exon", "R exon", "R exon"])
    return ref, eff_g, label


def test_gdna_planted_on_a_spike_in_fails_the_gate():
    ref, eff_g, label = _slots()
    n_g = np.array([7.0, 9.0, 0.0, 3.0])  # three gDNA fragments on the RNA-only spike-in
    gdna_refs = co.genomic_ref_ids(_manifest(["chr1"]), "c", REF_IDS)
    gate = _outside_gate(co.gdna_scope_verdicts(False, n_g, eff_g, label, ref, gdna_refs))
    assert not gate["vacuous"]
    assert not gate["ok"] and gate["mass"] == 3.0


def test_gdna_confined_to_the_genomic_references_passes_the_gate():
    ref, eff_g, label = _slots()
    n_g = np.array([7.0, 9.0, 0.0, 0.0])
    gdna_refs = co.genomic_ref_ids(_manifest(["chr1"]), "c", REF_IDS)
    gate = _outside_gate(co.gdna_scope_verdicts(False, n_g, eff_g, label, ref, gdna_refs))
    assert gate["ok"] and not gate["vacuous"]


def test_a_zero_gdna_condition_holds_every_reference_to_zero():
    ref, eff_g, label = _slots()
    n_g = np.array([1.0, 0.0, 0.0, 0.0])  # one gDNA fragment in a condition that simulated none
    gdna_refs = co.genomic_ref_ids(_manifest(["chr1"], n_gdna=0), "c", REF_IDS)
    assert gdna_refs.size == 0
    gate = _outside_gate(co.gdna_scope_verdicts(False, n_g, eff_g, label, ref, gdna_refs))
    assert not gate["ok"]


def test_the_genomic_set_does_not_follow_the_realized_counts():
    """The set is the manifest's, whatever the counts: no argument carries them."""
    gdna_refs = co.genomic_ref_ids(_manifest(["chr1"]), "c", REF_IDS)
    assert gdna_refs.tolist() == [0]


def test_a_genomic_reference_the_index_does_not_hold_is_refused():
    with pytest.raises(ValueError, match="chr9"):
        co.genomic_ref_ids(_manifest(["chr1", "chr9"]), "c", REF_IDS)


def test_a_condition_the_manifest_does_not_list_is_refused():
    with pytest.raises(KeyError, match="other"):
        co.genomic_ref_ids(_manifest(["chr1"]), "other", REF_IDS)
