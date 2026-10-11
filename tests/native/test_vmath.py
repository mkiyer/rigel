"""``native/vmath.h`` — the series `log` and `exp` over an array — against numpy, in ulps, on a million points
spanning the reader's arguments and the floating-point edges the kernel can reach."""

import numpy as np

from rigel import _honest_reader_impl as HR


def _ulps(got, want):
    return np.abs(got - want) / np.spacing(np.abs(want))


def test_vlog_is_within_two_ulps_of_numpy_on_the_readers_arguments():
    rng = np.random.default_rng(1)
    x = np.exp(rng.uniform(-40.0, 40.0, 1_000_000))
    x = np.concatenate(
        [
            x,
            [
                1.0,
                2.0,
                0.5,
                np.sqrt(2.0),
                np.sqrt(0.5),
                np.nextafter(1.0, 0.0),
                np.nextafter(1.0, 2.0),
                np.finfo(float).tiny,
                np.finfo(float).max,
                1e-300,
                1e300,
            ],
        ]
    )
    out = np.empty_like(x)
    HR.vlog(x, out)
    want = np.log(x)
    u = _ulps(out, want)
    u[want == 0.0] = np.abs(out[want == 0.0]) / np.spacing(1.0)
    assert out[x == 1.0].tolist() == [0.0], "log 1 must be exactly 0"
    assert u.max() <= 2.0, f"max error {u.max():.2f} ulp at x={x[np.argmax(u)]!r}"


def _ln2_distance_of_the_double():
    """|M_LN2 − ln 2|, with ln 2 = 2·atanh(1/3) summed to forty digits."""
    from decimal import Decimal, getcontext

    getcontext().prec = 45
    s = Decimal(1) / Decimal(3)
    ln2, term, k = Decimal(0), s, 0
    while term > Decimal(10) ** -44:
        ln2 += term / (2 * k + 1)
        term *= s * s
        k += 1
    ln2 *= 2
    return float(abs(Decimal(np.log(2.0)) - ln2))


def test_vexp_is_within_its_derived_ulp_bound_and_exact_at_the_edges():
    """3 ulp plus the double ln 2's own distance from ln 2 multiplied by the reduction's k = round(x / ln 2)."""
    rng = np.random.default_rng(2)
    x = np.concatenate(
        [
            rng.uniform(-745.0, 0.0, 1_000_000),
            rng.uniform(-40.0, 0.0, 1_000_000),
            [0.0, -1e-300, np.log(2.0), -np.log(2.0), -np.inf, -1e308, -745.0, -708.0, -709.0],
        ]
    )
    out = np.empty_like(x)
    HR.vexp(x, out)
    want = np.exp(x)
    assert out[x == 0.0].tolist() == [1.0], "exp 0 must be exactly 1"
    assert np.all(out[x == -np.inf] == 0.0)
    normal = want >= np.finfo(float).tiny
    rel = np.abs(out[normal] - want[normal]) / want[normal]
    eps = np.finfo(
        float
    ).eps  # 2^-52: three ulps at mantissa 1 is 3·eps relative, the largest an ulp can be
    allowed = 3.0 * eps + (np.abs(x[normal]) / np.log(2.0) + 0.5) * _ln2_distance_of_the_double()
    worst = np.argmax(rel - allowed)
    assert np.all(rel <= allowed), (
        f"relative error {rel[worst]:.3e} at x={x[normal][worst]!r} exceeds the bound {allowed[worst]:.3e}"
    )
    live = normal & (x >= -36.0)
    assert rel[live[normal]].max() <= 2e-15, (
        "within T nats of a sum's maximum the relative error must stay under 2e-15"
    )
    assert np.all(out[~normal] == 0.0) or np.all(
        np.abs(out[~normal] - want[~normal]) <= np.finfo(float).tiny
    ), "below the smallest normal the result must be 0 or a faithful denormal"
