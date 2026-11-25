"""`policy_benchmark.py` reads its two bars per stratum, strand half x capture, and never pools capture ON with OFF.

A strand half's "beats silence on k/n" and "worst row" over capture ON and OFF together is read over two strata: a
policy that wins every OFF row and loses every ON row reads as a draw.
"""

from __future__ import annotations

import io
from contextlib import redirect_stdout

import policy_benchmark as pb

STRANDED_OFF = "gdna_g50_ss_0.99_nascent_mid_capture_off"
STRANDED_ON = "gdna_g50_ss_0.99_nascent_mid_capture_on"
UNSTRANDED_ON = "gdna_g50_ss_0.50_nascent_mid_capture_on"


def _row(silent, transfer):
    return {"silent": {"total": float(silent)}, "transfer": {"total": float(transfer)}}


def _bars(table) -> str:
    buf = io.StringIO()
    with redirect_stdout(buf):
        pb._print_bars(table, ["silent", "transfer"])
    return buf.getvalue()


def _lines(out: str) -> list[str]:
    return [line for line in out.splitlines() if "beats silence" in line]


def test_one_on_row_and_one_off_row_are_two_strata():
    out = _bars({STRANDED_OFF: _row(100, 50), STRANDED_ON: _row(100, 150)})
    lines = _lines(out)
    assert len(lines) == 2, out
    (off,) = [line for line in lines if "capture OFF" in line]
    (on,) = [line for line in lines if "capture ON" in line]
    assert "1/1" in off and "0.50x" in off
    assert "0/1" in on and "1.50x" in on


def test_no_line_counts_rows_of_two_capture_states():
    out = _bars({STRANDED_OFF: _row(100, 50), STRANDED_ON: _row(100, 150)})
    assert all("/2" not in line for line in _lines(out)), out


def test_the_deferred_stratum_is_marked():
    out = _bars({UNSTRANDED_ON: _row(100, 50)})
    (line,) = _lines(out)
    assert "DEFERRED" in line, out


def test_the_bars_are_returned_per_stratum():
    rows = pb.stratum_bars(
        {STRANDED_OFF: _row(100, 50), STRANDED_ON: _row(100, 150)}, ["silent", "transfer"]
    )
    assert {(r["stratum"], r["better"], r["n"]) for r in rows} == {
        (("stranded", "capture OFF"), 1, 1),
        (("stranded", "capture ON"), 0, 1),
    }
