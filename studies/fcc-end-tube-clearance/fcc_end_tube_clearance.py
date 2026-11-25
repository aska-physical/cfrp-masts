#!/usr/bin/env python3
"""Exact non-incident tube clearance at the proved m=1 fcc end state.

Source, on Foundation-0/incentre-lattices main: results/fold-m1/PROOF.md.
Dots sit at ±e_i. Lattice rows are (2,-1,-1), (1,2,-1), (1,-1,2), in units
of K = 2 - sqrt(2). Each cab bar has squared length 2, which is (2R)^2 with
R = sqrt(2) - 1. This file does not read that repository. It rebuilds the
closed form and certifies segment distances in fractions.Fraction.

Scope: the end state only. The path is not certified. Shared vertices are
joints, not strut collisions. No bench measurement. Novelty default N0.
"""

from __future__ import annotations

import json
from fractions import Fraction

CHECK_NAME = "fcc-end-tube-clearance"
SOURCE = "Foundation-0/incentre-lattices results/fold-m1/PROOF.md (fcc row)"


class CheckFailure(Exception):
    def __init__(self, code: str):
        super().__init__(code)
        self.code = code


def need(ok: bool, code: str) -> None:
    if not ok:
        raise CheckFailure(code)


def add(p, q):
    return (p[0] + q[0], p[1] + q[1], p[2] + q[2])


def sub(p, q):
    return (p[0] - q[0], p[1] - q[1], p[2] - q[2])


def dot(p, q) -> int:
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]


def smul(n: int, p):
    return (n * p[0], n * p[1], n * p[2])


def lattice_point(n, rows):
    return add(add(smul(n[0], rows[0]), smul(n[1], rows[1])), smul(n[2], rows[2]))


def fcc_m1():
    """Proved m=1 fcc end state, integer coordinates in units of K."""
    dots = (
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    )
    rows = (
        (2, -1, -1),
        (1, 2, -1),
        (1, -1, 2),
    )
    # 12 octahedron edges: every non-antipodal pair. Antipodes are (0,1), (2,3), (4,5).
    bars = []
    antipode = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
    for i in range(6):
        for j in range(i + 1, 6):
            if antipode[i] != j:
                bars.append((i, j, (0, 0, 0)))
    # Linking bar on axis a: +e_a (even index 2a) to -e_a of the next cell.
    for a in range(3):
        bars.append((2 * a, 2 * a + 1, (1 if a == 0 else 0, 1 if a == 1 else 0, 1 if a == 2 else 0)))
    return dots, rows, tuple(bars)


def endpoint(dots, rows, bar, shift=(0, 0, 0)):
    i, j, n = bar
    origin = lattice_point(shift, rows)
    a = add(dots[i], origin)
    b = add(add(dots[j], lattice_point(n, rows)), origin)
    return a, b


def bars_exact(dots, rows, bars, length2: int = 2) -> bool:
    for bar in bars:
        a, b = endpoint(dots, rows, bar)
        if dot(sub(b, a), sub(b, a)) != length2:
            return False
    return True


def det3(rows) -> int:
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def share_vertex(dots, rows, bar1, shift1, bar2, shift2) -> bool:
    a1, b1 = endpoint(dots, rows, bar1, shift1)
    a2, b2 = endpoint(dots, rows, bar2, shift2)
    ends1 = (a1, b1)
    ends2 = (a2, b2)
    return any(p == q for p in ends1 for q in ends2)


def point_seg_d2(p, a, b) -> Fraction:
    ab = sub(b, a)
    ab2 = dot(ab, ab)
    if ab2 == 0:
        d = sub(p, a)
        return Fraction(dot(d, d), 1)
    t = Fraction(dot(sub(p, a), ab), ab2)
    if t < 0:
        t = Fraction(0, 1)
    elif t > 1:
        t = Fraction(1, 1)
    # q = a + t (b-a); |p-q|^2
    # (p-a - t ab) · same
    pa = sub(p, a)
    # |pa|^2 - 2 t pa·ab + t^2 |ab|^2
    return Fraction(dot(pa, pa), 1) - 2 * t * Fraction(dot(pa, ab), 1) + t * t * Fraction(ab2, 1)


def seg_seg_d2(a, b, c, d) -> Fraction:
    """Exact squared distance between segments ab and cd."""
    d1 = sub(b, a)
    d2 = sub(d, c)
    r = sub(a, c)
    aa = dot(d1, d1)
    ee = dot(d2, d2)
    bb = dot(d1, d2)
    best = min(
        point_seg_d2(a, c, d),
        point_seg_d2(b, c, d),
        point_seg_d2(c, a, b),
        point_seg_d2(d, a, b),
    )
    denom = aa * ee - bb * bb
    if denom != 0 and aa != 0 and ee != 0:
        # d1·r + aa u - bb v = 0
        # -d2·r - bb u + ee v = 0
        cr = dot(d1, r)
        dr = dot(d2, r)
        u = Fraction(ee * cr - bb * (-dr), denom)  # wait, second eq is -dr - bb u + ee v = 0
        # From eq1: aa u - bb v = -cr
        # From eq2: -bb u + ee v = dr     because -dr - bb u + ee v = 0 => -bb u + ee v = dr
        # u = ( -cr * ee - (-bb) * dr ) / (aa*ee - bb*bb)? Use Cramer.
        # | aa  -bb | |u|   |-cr|
        # |-bb   ee | |v| = | dr|
        u = Fraction((-cr) * ee - (-bb) * dr, denom)
        v = Fraction(aa * dr - (-bb) * (-cr), denom)
        if Fraction(0, 1) <= u <= 1 and Fraction(0, 1) <= v <= 1:
            # w = r + u d1 - v d2
            # compute |w|^2 exactly
            # w = r + u d1 - v d2
            # Do it componentwise with Fractions
            w = []
            for k in range(3):
                w.append(Fraction(r[k], 1) + u * d1[k] - v * d2[k])
            interior = w[0] * w[0] + w[1] * w[1] + w[2] * w[2]
            if interior < best:
                best = interior
    return best


SHIFTS = tuple(
    (i, j, k)
    for i in (-1, 0, 1)
    for j in (-1, 0, 1)
    for k in (-1, 0, 1)
)


def min_nonincident(dots, rows, bars):
    """Minimum squared centerline distance between bars that do not share a vertex.

    Periodic images use lattice shifts in {-1,0,1}^3. The caller must check that
    every nonzero lattice vector outside that coefficient box is too long to hide
    a closer pair (see lattice_box_is_exhaustive).
    """
    best = None
    witness = None
    n = len(bars)
    for i in range(n):
        a, b = endpoint(dots, rows, bars[i])
        for j in range(n):
            for shift in SHIFTS:
                if i == j and shift == (0, 0, 0):
                    continue
                # unordered pair: keep i < j, or i == j and shift in the positive half
                if j < i:
                    continue
                if i == j and shift < (0, 0, 0):
                    continue
                if share_vertex(dots, rows, bars[i], (0, 0, 0), bars[j], shift):
                    continue
                c, d = endpoint(dots, rows, bars[j], shift)
                d2 = seg_seg_d2(a, b, c, d)
                if best is None or d2 < best:
                    best = d2
                    witness = {
                        "bar_i": i,
                        "bar_j": j,
                        "shift": list(shift),
                        "d2": str(d2),
                    }
    need(best is not None, "NO_PAIR")
    need(best > 0, "ZERO_CLEARANCE")
    return best, witness


def shortest_lattice_outside_box(rows, coeff_limit: int = 4) -> int:
    """Minimum |s|^2 over max|n_i| <= coeff_limit, excluding the open search box max|n_i| <= 1.

    Vectors outside this coefficient box are not measured. The mast cell uses
    coeff_limit 4. Its shortest supercell vectors sit inside max|n_i| <= 1.
    """
    best = None
    for i in range(-coeff_limit, coeff_limit + 1):
        for j in range(-coeff_limit, coeff_limit + 1):
            for k in range(-coeff_limit, coeff_limit + 1):
                if max(abs(i), abs(j), abs(k)) <= 1:
                    continue
                s = lattice_point((i, j, k), rows)
                d2 = dot(s, s)
                if best is None or d2 < best:
                    best = d2
    need(best is not None, "NO_OUTSIDE_VECTOR")
    return best


def assert_search_box(rows, found_d2: Fraction, coeff_limit: int = 4) -> None:
    """Reject a lattice whose unsearched images could hide a closer strut pair.

    A bar image shifted by s is at least |s| - 2*sqrt(2) from a fundamental bar,
    because each bar has squared length 2. That lower bound beats found_d2 when
    left = |s|^2_min_outside - 8 - found_d2 is positive and left^2 > 32 found_d2.
    """
    outside = shortest_lattice_outside_box(rows, coeff_limit)
    left = Fraction(outside, 1) - 8 - found_d2
    need(left > 0, "SEARCH_BOX_GAP")
    need(left * left > 32 * found_d2, "SEARCH_BOX_SQUARE")


def tube_clears(r_out_squared: Fraction, min_d2: Fraction) -> bool:
    """Open tubes of radius r clear when (2r)^2 < min centerline distance squared."""
    return (4 * r_out_squared) < min_d2


def witness_geometry(dots, rows, bars, witness) -> dict:
    """Endpoints and mid-span feet for the closest pair. Coordinates are fractions in units of K."""
    i, j, shift = witness["bar_i"], witness["bar_j"], tuple(witness["shift"])
    a, b = endpoint(dots, rows, bars[i])
    c, d = endpoint(dots, rows, bars[j], shift)
    return {
        "seg1": [list(a), list(b)],
        "seg2": [list(c), list(d)],
        "mid1": [str(Fraction(a[k] + b[k], 2)) for k in range(3)],
        "mid2": [str(Fraction(c[k] + d[k], 2)) for k in range(3)],
    }


def count_at_minimum(dots, rows, bars, min_d2: Fraction) -> int:
    n = len(bars)
    hits = 0
    for i in range(n):
        a, b = endpoint(dots, rows, bars[i])
        for j in range(i, n):
            for shift in SHIFTS:
                if i == j and shift <= (0, 0, 0):
                    continue
                if share_vertex(dots, rows, bars[i], (0, 0, 0), bars[j], shift):
                    continue
                c, d = endpoint(dots, rows, bars[j], shift)
                if seg_seg_d2(a, b, c, d) == min_d2:
                    hits += 1
    return hits


def ball_contacts(dots, rows):
    """Exact contact count for balls of squared diameter 2. Returns degrees, min pair d^2."""
    n = len(dots)
    degree = [0] * n
    min_d2 = None
    for i in range(n):
        for j in range(n):
            for shift in SHIFTS:
                if i == j and shift == (0, 0, 0):
                    continue
                # count each unordered pair once
                if j < i:
                    continue
                if i == j and shift < (0, 0, 0):
                    continue
                delta = sub(add(dots[j], lattice_point(shift, rows)), dots[i])
                d2 = dot(delta, delta)
                if min_d2 is None or d2 < min_d2:
                    min_d2 = d2
                if d2 == 2:
                    degree[i] += 1
                    if i != j:
                        degree[j] += 1
                    else:
                        # a ball touching its own image: one contact, already counted
                        pass
    return degree, min_d2


def certificate(r_out_squared: Fraction = Fraction(0, 1), dots=None, rows=None, bars=None) -> dict:
    if dots is None:
        dots, rows, bars = fcc_m1()
    need(len(bars) == 15, "BAR_COUNT")
    need(bars_exact(dots, rows, bars), "BARS_NOT_EXACT")
    need(det3(rows) == 12, "DENSITY_DET")
    degree, ball_min = ball_contacts(dots, rows)
    need(ball_min >= 2, "BALL_OVERLAP")
    need(degree == [12, 12, 12, 12, 12, 12], "CONTACT_COUNT")
    # (2R)^2 = 2 K^2 and det L = 2 * 6 * K^3 = 12 K^3, the proved density identity.
    min_d2, witness = min_nonincident(dots, rows, bars)
    assert_search_box(rows, min_d2)
    geom = witness_geometry(dots, rows, bars, witness)
    hits = count_at_minimum(dots, rows, bars, min_d2)
    clears = tube_clears(r_out_squared, min_d2)
    # r_max^2 = min_d2 / 4. Report it as a reduced fraction. It is not a measured radius.
    return {
        "check": CHECK_NAME,
        "oracle": "exact",
        "novelty": "N0",
        "provenance": "source_claim",
        "source": SOURCE,
        "units": "K = 2 - sqrt(2); coordinates are integers in that unit",
        "scope": {
            "end_state": "m=1 fcc closed form",
            "path": "NOT_SHOWN",
            "shared_vertices": "joints, excluded",
            "bench_measurement": False,
            "periodic_shifts": "[-1,1]^3",
        },
        "bars": 15,
        "det_L_over_K3": 12,
        "bar_length_squared": "2",
        "min_nonincident_d2": str(min_d2),
        "pairs_at_minimum": hits,
        "min_d2_over_bar_length_squared": str(min_d2 / 2),
        "r_out_max_squared": str(min_d2 / 4),
        "ball_radius_squared": "1/2",
        "witness": witness,
        "witness_endpoints": geom,
        "r_out_squared_under_test": str(r_out_squared),
        "clears": clears,
    }


def decide(r_out_squared: Fraction) -> str:
    """PASS a tube that clears. FAIL one that meets a non-incident strut."""
    cert = certificate(r_out_squared)
    if cert["clears"]:
        return "PASS"
    raise CheckFailure("TUBE_COLLISION")


def canonical(doc: dict) -> str:
    return json.dumps(doc, sort_keys=True, separators=(",", ":"))


def main() -> None:
    cert = certificate(Fraction(0, 1))
    print(
        json.dumps(
            {
                "min_nonincident_d2": cert["min_nonincident_d2"],
                "r_out_max_squared": cert["r_out_max_squared"],
                "pairs_at_minimum": cert["pairs_at_minimum"],
                "clears": cert["clears"],
            },
            sort_keys=True,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
