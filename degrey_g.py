"""De Grey's 1581-vertex non-4-colourable graph G, built from the explicit
recipe in section 5.1 of arXiv:1804.02385.

    1. S   39 seed points, listed exactly below
    2. S_a all images of S under the 6 rotations by multiples of 60 degrees
           and/or negation of the y-coordinate            -> 397 vertices
    3. S_b S_a rotated anticlockwise about the origin by 2·arcsin(1/4)
    4. Y   S_a union S_b, with (1/3, 0) and (-1/3, 0) deleted
    5. Y_a Y rotated anticlockwise about (-2, 0) by  pi/2 + arcsin(1/8)
    6. Y_b Y rotated anticlockwise about (-2, 0) by  pi/2 - arcsin(1/8)
    7. G   Y_a union Y_b                                  -> 1581 vertices

The coordinates need √3 and √11 (so √33 comes along); the first rotation needs
√15 and the last two need √7.  Hence Q(√3, √5, √7, √11), of degree 16.

Exact cos/sin used:
    60 degrees            1/2,   √3/2
    2·arcsin(1/4)         7/8,   √15/8
    pi/2 + arcsin(1/8)   -1/8,   3√7/8
    pi/2 - arcsin(1/8)    1/8,   3√7/8
"""

from field import get_field
from udg import UDG

PRIMES = (3, 5, 7, 11)


def C(fld, *terms):
    """A coordinate from (radicand, numerator, denominator) terms.

    C(f, (33, 1, 6), (1, -3, 6)) is (√33 - 3)/6.
    """
    out = fld.zero
    for rad, num, den in terms:
        piece = fld.rat(num, den) if rad == 1 else fld.surd(rad, num, den)
        out = fld.add(out, piece)
    return out


def seed_points(fld):
    """The 39 points of S, transcribed from section 5.1.

    Note √12 = 2√3 and 1/√12 = √3/6 and 1/√3 = √3/3.
    """
    z = ()          # empty term list -> 0
    P = []
    add = P.append

    add((C(fld), C(fld)))                                          # (0, 0)
    add((C(fld, (1, 1, 3)), C(fld)))                               # (1/3, 0)
    add((C(fld, (1, 1, 1)), C(fld)))                               # (1, 0)
    add((C(fld, (1, 2, 1)), C(fld)))                               # (2, 0)
    add((C(fld, (33, 1, 6), (1, -3, 6)), C(fld)))                  # ((√33-3)/6, 0)
    add((C(fld, (1, 1, 2)), C(fld, (3, 1, 6))))                    # (1/2, 1/√12)
    add((C(fld, (1, 1, 1)), C(fld, (3, 1, 3))))                    # (1, 1/√3)
    add((C(fld, (1, 3, 2)), C(fld, (3, 1, 2))))                    # (3/2, √3/2)

    add((C(fld, (1, 7, 6)), C(fld, (11, 1, 6))))                   # (7/6, √11/6)
    add((C(fld, (1, 1, 6)), C(fld, (3, 2, 6), (11, -1, 6))))       # (1/6, (√12-√11)/6)
    add((C(fld, (1, 5, 6)), C(fld, (3, 2, 6), (11, -1, 6))))       # (5/6, (√12-√11)/6)

    add((C(fld, (1, 2, 3)), C(fld, (11, 1, 6), (3, -1, 6))))       # (2/3, (√11-√3)/6)
    add((C(fld, (1, 2, 3)), C(fld, (3, 3, 6), (11, -1, 6))))       # (2/3, (3√3-√11)/6)
    add((C(fld, (33, 1, 6)), C(fld, (3, 1, 6))))                   # (√33/6, 1/√12)

    add((C(fld, (33, 1, 6), (1, 3, 6)), C(fld, (3, 1, 3))))
    add((C(fld, (33, 1, 6), (1, 1, 6)), C(fld, (3, 3, 6), (11, -1, 6))))
    add((C(fld, (33, 1, 6), (1, -1, 6)), C(fld, (3, 3, 6), (11, -1, 6))))

    add((C(fld, (33, 1, 6), (1, 1, 6)), C(fld, (11, 1, 6), (3, -1, 6))))
    add((C(fld, (33, 1, 6), (1, -1, 6)), C(fld, (11, 1, 6), (3, -1, 6))))

    add((C(fld, (33, 1, 6), (1, -2, 6)), C(fld, (3, 2, 6), (11, -1, 6))))
    add((C(fld, (33, 1, 6), (1, -4, 6)), C(fld, (3, 2, 6), (11, -1, 6))))

    add((C(fld, (33, 1, 12), (1, 13, 12)), C(fld, (11, 1, 12), (3, -1, 12))))
    add((C(fld, (33, 1, 12), (1, 11, 12)), C(fld, (3, 1, 12), (11, 1, 12))))

    add((C(fld, (33, 1, 12), (1, 9, 12)), C(fld, (11, 1, 4), (3, -1, 4))))
    add((C(fld, (33, 1, 12), (1, 9, 12)), C(fld, (3, 3, 12), (11, 1, 12))))

    add((C(fld, (33, 1, 12), (1, 7, 12)), C(fld, (3, 1, 12), (11, 1, 12))))
    add((C(fld, (33, 1, 12), (1, 7, 12)), C(fld, (3, 3, 12), (11, -1, 12))))

    add((C(fld, (33, 1, 12), (1, 5, 12)), C(fld, (3, 5, 12), (11, -1, 12))))
    add((C(fld, (33, 1, 12), (1, 5, 12)), C(fld, (11, 1, 12), (3, -1, 12))))

    add((C(fld, (33, 1, 12), (1, 3, 12)), C(fld, (11, 3, 12), (3, -5, 12))))
    add((C(fld, (33, 1, 12), (1, 3, 12)), C(fld, (3, 1, 12), (11, 1, 12))))

    add((C(fld, (33, 1, 12), (1, 3, 12)), C(fld, (3, 3, 12), (11, -1, 12))))
    add((C(fld, (33, 1, 12), (1, 1, 12)), C(fld, (11, 1, 12), (3, -1, 12))))

    add((C(fld, (33, 1, 12), (1, -1, 12)), C(fld, (3, 3, 12), (11, -1, 12))))
    add((C(fld, (33, 1, 12), (1, -3, 12)), C(fld, (11, 1, 12), (3, -1, 12))))

    add((C(fld, (1, 15, 12), (33, -1, 12)), C(fld, (11, 1, 4), (3, -1, 4))))
    add((C(fld, (1, 15, 12), (33, -1, 12)), C(fld, (3, 7, 12), (11, -3, 12))))

    add((C(fld, (1, 13, 12), (33, -1, 12)), C(fld, (3, 3, 12), (11, -1, 12))))
    add((C(fld, (1, 11, 12), (33, -1, 12)), C(fld, (11, 1, 12), (3, -1, 12))))

    return P


# ------------------------------------------------------------- assembly ----

def _dedup(points):
    seen, out = set(), []
    for p in points:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def build_Sa(fld):
    """S closed under the 6 multiples of 60 degrees and y-negation."""
    cos60, sin60 = fld.rat(1, 2), fld.surd(3, 1, 2)
    pts = []
    for p in seed_points(fld):
        q = p
        for _ in range(6):
            pts.append(q)
            pts.append((q[0], fld.neg(q[1])))          # reflect in the x-axis
            q = fld.rotate(q, cos60, sin60)
    return _dedup(pts)


def build_Y(fld, Sa=None):
    """S_a together with S_b, minus the two points (±1/3, 0)."""
    Sa = Sa or build_Sa(fld)
    cos_t, sin_t = fld.rat(7, 8), fld.surd(15, 1, 8)   # 2·arcsin(1/4)
    Sb = [fld.rotate(p, cos_t, sin_t) for p in Sa]
    drop = {(C(fld, (1, 1, 3)), fld.zero), (C(fld, (1, -1, 3)), fld.zero)}
    return [p for p in _dedup(Sa + Sb) if p not in drop]


def build_G(fld=None):
    """The 1581-vertex graph G."""
    fld = fld or get_field(PRIMES)
    Y = build_Y(fld)
    centre = (fld.rat(-2), fld.zero)
    s = fld.surd(7, 3, 8)                              # sin = 3√7/8 for both
    Ya = [fld.rotate(p, fld.rat(-1, 8), s, centre) for p in Y]   # pi/2 + arcsin(1/8)
    Yb = [fld.rotate(p, fld.rat(1, 8), s, centre) for p in Y]    # pi/2 - arcsin(1/8)
    return UDG(fld, _dedup(Ya + Yb))


if __name__ == "__main__":
    import time
    fld = get_field(PRIMES)
    t0 = time.time()
    Sa = build_Sa(fld)
    print(f"S    : {len(seed_points(fld))} points   (expected 39)")
    print(f"S_a  : {len(Sa)} points   (expected 397)")
    Y = build_Y(fld, Sa)
    print(f"Y    : {len(Y)} points   (expected 791)")
    G = build_G(fld)
    print(f"G    : {len(G)} points   (expected 1581)   [{time.time()-t0:.1f}s]")
    t0 = time.time()
    print(f"       {len(G.edges)} edges   [{time.time()-t0:.1f}s]")
