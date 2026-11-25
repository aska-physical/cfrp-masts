#!/usr/bin/env python3
"""Must-fail controls for the m=1 fcc tube-clearance checker.

The distance formula here is not the one in the checker. It uses the scalar
triple product for skew lines and clamps to the four endpoints otherwise.
"""

import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fcc_end_tube_clearance import (
    CheckFailure,
    bars_exact,
    certificate,
    decide,
    fcc_m1,
    need,
    seg_seg_d2,
    share_vertex,
    tube_clears,
)


def cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def dot(p, q):
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1], p[2] - q[2])


def point_seg_d2(p, a, b) -> Fraction:
    ab = sub(b, a)
    ab2 = dot(ab, ab)
    t = Fraction(dot(sub(p, a), ab), ab2)
    if t < 0:
        t = Fraction(0, 1)
    elif t > 1:
        t = Fraction(1, 1)
    pa = sub(p, a)
    return Fraction(dot(pa, pa), 1) - 2 * t * Fraction(dot(pa, ab), 1) + t * t * Fraction(ab2, 1)


def independent_seg_d2(a, b, c, d) -> Fraction:
    d1 = sub(b, a)
    d2 = sub(d, c)
    cr = cross(d1, d2)
    denom = dot(cr, cr)
    best = min(
        point_seg_d2(a, c, d),
        point_seg_d2(b, c, d),
        point_seg_d2(c, a, b),
        point_seg_d2(d, a, b),
    )
    if denom == 0:
        return best
    # Feet on the infinite lines: solve the same 2x2, then accept only interior feet.
    r = sub(a, c)
    aa = dot(d1, d1)
    ee = dot(d2, d2)
    bb = dot(d1, d2)
    cr_ = dot(d1, r)
    dr = dot(d2, r)
    u = Fraction((-cr_) * ee + bb * dr, denom)
    v = Fraction(aa * dr - bb * cr_, denom)
    if Fraction(0, 1) <= u <= 1 and Fraction(0, 1) <= v <= 1:
        num = dot(sub(c, a), cr)
        line_d2 = Fraction(num * num, denom)
        if line_d2 < best:
            best = line_d2
    return best


def expect_fail(fn, code: str) -> None:
    try:
        fn()
    except CheckFailure as exc:
        need(exc.code == code, "WRONG_CODE_" + exc.code)
        return
    raise CheckFailure("EXPECTED_" + code)


def test_dual_scan_matches() -> None:
    """Every non-incident pair, both distance formulas, same minimum."""
    from fcc_end_tube_clearance import SHIFTS, endpoint

    dots, rows, bars = fcc_m1()
    best = None
    n = len(bars)
    for i in range(n):
        a, b = endpoint(dots, rows, bars[i])
        for j in range(i, n):
            for shift in SHIFTS:
                if i == j and shift <= (0, 0, 0):
                    continue
                if share_vertex(dots, rows, bars[i], (0, 0, 0), bars[j], shift):
                    continue
                c, d = endpoint(dots, rows, bars[j], shift)
                left = seg_seg_d2(a, b, c, d)
                right = independent_seg_d2(a, b, c, d)
                need(left == right, "DUAL_PAIR")
                if best is None or left < best:
                    best = left
    need(best == 1, "DUAL_MIN")
    # Opposite octahedron edges stay a full bar apart. They are not the witness.
    opp = independent_seg_d2((1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0))
    need(opp == 2, "OPPOSITE_EDGES")


def test_closed_form_and_dual_distance() -> None:
    dots, rows, bars = fcc_m1()
    need(len(bars) == 15, "BARS")
    need(bars_exact(dots, rows, bars), "LENGTHS")
    cert = certificate(Fraction(0, 1))
    need(cert["min_nonincident_d2"] == "1", "MIN")
    need(cert["r_out_max_squared"] == "1/4", "RMAX")
    need(cert["ball_radius_squared"] == "1/2", "BALL")
    need(cert["pairs_at_minimum"] >= 1, "HITS")
    need(cert["clears"] is True, "ZERO")
    # Hand witness: edge (+e_x, -e_y) and the z linking bar.
    a, b = (1, 0, 0), (0, -1, 0)
    c, d = (0, 0, 1), (1, -1, 1)
    need(independent_seg_d2(a, b, c, d) == 1, "HAND")
    need(seg_seg_d2(a, b, c, d) == independent_seg_d2(a, b, c, d), "DUAL_WITNESS")
    ends = cert["witness_endpoints"]
    need(ends["mid1"] == ["1/2", "-1/2", "0"], "MID1")
    need(ends["mid2"] == ["1/2", "-1/2", "1"], "MID2")


def test_pass_fail_split() -> None:
    need(decide(Fraction(0, 1)) == "PASS", "PASS0")
    need(decide(Fraction(1, 5)) == "PASS", "PASS_THIN")
    # Strict: radius squared 1/4 meets the witness. The node ball is wider and also meets it.
    expect_fail(lambda: decide(Fraction(1, 4)), "TUBE_COLLISION")
    expect_fail(lambda: decide(Fraction(1, 2)), "TUBE_COLLISION")
    need(tube_clears(Fraction(1, 4), Fraction(1, 1)) is False, "LIMIT_CLOSED")
    need(tube_clears(Fraction(1, 5), Fraction(1, 1)) is True, "OPEN")


def test_incidence_is_load_bearing() -> None:
    dots, rows, bars = fcc_m1()
    # Two octahedron edges that share +e_x. Their centerlines meet, so an
    # unfiltered minimum would be 0 and every tube would "collide".
    a, b = endpoint_of(dots, rows, bars[0])
    # Find a second bar that shares an endpoint with bar 0.
    found = False
    for bar in bars[1:]:
        if share_vertex(dots, rows, bars[0], (0, 0, 0), bar, (0, 0, 0)):
            c, d = endpoint_of(dots, rows, bar)
            need(seg_seg_d2(a, b, c, d) == 0, "JOINT_DISTANCE")
            found = True
            break
    need(found, "NO_JOINT")


def endpoint_of(dots, rows, bar):
    i, j, n = bar
    start = dots[i]
    far = (
        dots[j][0] + n[0] * rows[0][0] + n[1] * rows[1][0] + n[2] * rows[2][0],
        dots[j][1] + n[0] * rows[0][1] + n[1] * rows[1][1] + n[2] * rows[2][1],
        dots[j][2] + n[0] * rows[0][2] + n[1] * rows[1][2] + n[2] * rows[2][2],
    )
    return start, far


def test_broken_end_state_is_rejected() -> None:
    dots, rows, bars = fcc_m1()
    moved = (dots[0], (dots[1][0] + 1, dots[1][1], dots[1][2])) + dots[2:]
    need(bars_exact(moved, rows, bars) is False, "BROKEN")
    bad = bars[:14] + ((0, 1, (0, 0, 0)),)
    expect_fail(lambda: certificate(Fraction(0, 1), dots, rows, bad), "BARS_NOT_EXACT")


# Lattices with exact bar lengths that fail exactly one later guard.
# Found by enumerating integer solutions of the three linking-bar equations.
DET_ROWS = ((1, -1, 0), (-1, 1, 0), (-1, -1, 2))
BALL_ROWS = ((1, -1, 0), (1, 3, 0), (-1, 0, 3))
DEG_ROWS = ((2, -1, -1), (-1, 2, -1), (0, 1, 3))
BOX_ROWS = ((2, -1, -1), (-1, 3, 0), (1, -1, 2))
# shortest outside-box vector squared is 9, so left = 9 - 8 - 1 = 0
GAP_EQ_ROWS = ((-2, -2, -2), (-2, -2, 0), (-2, 1, -2))
# shortest outside-box vector squared is 18, so left^2 = 64 = 32 * q at q = 2
SQUARE_EQ_ROWS = ((-2, -2, -2), (-2, -2, 1), (-2, 1, -2))


def test_each_guard_has_a_witness() -> None:
    from fcc_end_tube_clearance import (
        assert_search_box,
        min_nonincident,
        shortest_lattice_outside_box,
    )

    dots, rows, bars = fcc_m1()
    expect_fail(lambda: need(False, "NEED_BODY"), "NEED_BODY")
    expect_fail(lambda: certificate(Fraction(0, 1), dots, rows, bars[:14]), "BAR_COUNT")
    expect_fail(lambda: certificate(Fraction(0, 1), dots, DET_ROWS, bars), "DENSITY_DET")
    expect_fail(lambda: certificate(Fraction(0, 1), dots, BALL_ROWS, bars), "BALL_OVERLAP")
    expect_fail(lambda: certificate(Fraction(0, 1), dots, DEG_ROWS, bars), "CONTACT_COUNT")
    expect_fail(lambda: certificate(Fraction(0, 1), dots, BOX_ROWS, bars), "SEARCH_BOX_SQUARE")
    expect_fail(lambda: min_nonincident(dots, rows, ()), "NO_PAIR")
    # Overlapping non-incident segments. Distance is zero, so a non-strict
    # positivity test would accept them.
    cross_dots = ((0, 0, 0), (2, 0, 0), (1, 0, 0), (3, 0, 0))
    cross_rows = ((10, 0, 0), (0, 10, 0), (0, 0, 10))
    cross_bars = ((0, 1, (0, 0, 0)), (2, 3, (0, 0, 0)))
    expect_fail(lambda: min_nonincident(cross_dots, cross_rows, cross_bars), "ZERO_CLEARANCE")
    expect_fail(lambda: shortest_lattice_outside_box(rows, 1), "NO_OUTSIDE_VECTOR")
    expect_fail(lambda: assert_search_box(GAP_EQ_ROWS, Fraction(1, 1)), "SEARCH_BOX_GAP")
    expect_fail(lambda: assert_search_box(SQUARE_EQ_ROWS, Fraction(2, 1)), "SEARCH_BOX_SQUARE")
    # Empty dot set: the contact scan never sees a pair.
    from fcc_end_tube_clearance import ball_contacts

    degree, ball_min = ball_contacts((), rows)
    need(ball_min is None and degree == [], "EMPTY_BALLS")
    expect_fail(
        lambda: certificate(Fraction(0, 1), (), ((1, 0, 0), (0, 1, 0), (0, 0, 1)), ()),
        "BAR_COUNT",
    )


def test_receipt_matches() -> None:
    raw = json.loads(Path(__file__).with_name("fcc_end_tube_clearance_certificate.json").read_text())
    cert = certificate(Fraction(0, 1))
    for key in (
        "min_nonincident_d2",
        "r_out_max_squared",
        "ball_radius_squared",
        "pairs_at_minimum",
        "det_L_over_K3",
        "bars",
    ):
        need(cert[key] == raw[key], "RECEIPT")
    # Do not pin checker bytes here. aletheia mutate re-serialises the checker
    # for its null mutant, and a hash pin would fail that control.
    need(raw["novelty"] == "N0", "NOVELTY")
    need(raw["scope"]["path"] == "NOT_SHOWN", "PATH")
    need(raw["scope"]["bench_measurement"] is False, "BENCH")


def main() -> None:
    test_closed_form_and_dual_distance()
    test_dual_scan_matches()
    test_pass_fail_split()
    test_incidence_is_load_bearing()
    test_broken_end_state_is_rejected()
    test_each_guard_has_a_witness()
    test_receipt_matches()
    print("r5 ok")


if __name__ == "__main__":
    main()
