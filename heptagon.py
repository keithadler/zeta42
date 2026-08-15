"""Haugland's 21-vertex 7-fold symmetric unit-distance graph, exactly.

From arXiv:2608.04542 (J. K. Haugland, August 2026): three concentric unit-edge
star polygons {7/1}, {7/2}, {7/3}, phased so that P_j, Q_j, R_j form a unit
equilateral triangle for each j.  It seeds a 2131-vertex 5-chromatic graph that
contains no Moser spindle -- the first 5-chromatic construction that does not.

Placing them exactly is the whole trick.  A unit-edge {7/k} star polygon has
circumradius 1/(2*sin(k*pi/7)), which is not a quadratic surd -- but

    V_j = zeta^j / (zeta^k - 1)        with zeta = zeta_7

has consecutive differences V_(j+k) - V_j = zeta^j exactly, so every edge is a
root of unity and every vertex is in Q(zeta_7).  No square roots needed, and no
floating point anywhere.

The phase offsets below were found by exhaustive search over zeta_7^a * zeta_6^b
and then confirmed exactly; there are two solutions, mirror images of each
other.  Both need a sixth root of unity for the equilateral triangles, so the
ambient field is Q(zeta_42), of degree 12.
"""

from cyclo import get_cyclofield
from udg import UDG

N = 42          # lcm(7, 6): sevenfold symmetry and equilateral triangles


def field():
    return get_cyclofield(N)


def star(fld, k, phase=0):
    """Unit-edge star polygon {7/k}, optionally rotated by zeta_42^phase."""
    z7 = fld.zeta(6)                       # zeta_42^6 = zeta_7
    denom = fld.sub(fld.zeta(6 * k), fld.one)
    base = fld.inv(denom)
    if phase:
        base = fld.mul(base, fld.zeta(phase))
    return [fld.mul(fld.zeta(6 * j), base) for j in range(7)]


def build_H(fld=None, mirror=False):
    """The 21-vertex graph: heptagon P, heptagrams Q and R, phased to triangle.

    Phases in units of zeta_42: zeta_7^a = zeta_42^(6a), zeta_6^b = zeta_42^(7b).
        Q *= zeta_7^4 * zeta_6^2        R *= zeta_7^1 * zeta_6^1
    (mirror=True takes the other solution, zeta_6^4 and zeta_6^5.)
    """
    fld = fld or field()
    qb, rb = (4, 2), (1, 1)
    if mirror:
        qb, rb = (4, 4), (1, 5)
    P = star(fld, 1)
    Q = star(fld, 2, phase=6 * qb[0] + 7 * qb[1])
    R = star(fld, 3, phase=6 * rb[0] + 7 * rb[1])
    return P, Q, R


def graph_H(fld=None, mirror=False):
    fld = fld or field()
    P, Q, R = build_H(fld, mirror)
    return UDG(fld, P + Q + R)


def unit_vectors(fld, g):
    """The distinct edge directions of g, as field elements of modulus 1.

    Haugland's construction uses these as the generating "arcs": every edge of
    the big graph is a translation by one of them.
    """
    seen = {}
    for i, j in g.edges:
        for a, b in ((i, j), (j, i)):
            d = fld.sub(g.points[b], g.points[a])
            if d not in seen:
                seen[d] = True
    return list(seen)


if __name__ == "__main__":
    fld = field()
    P, Q, R = build_H(fld)

    print("unit-edge star polygons, exact:")
    for nm, S, k in (("P {7/1}", P, 1), ("Q {7/2}", Q, 2), ("R {7/3}", R, 3)):
        edge_ok = all(fld.is_unit(S[j], S[(j + k) % 7]) for j in range(7))
        rad = abs(complex(*fld.pfloat(S[0])))
        print(f"  {nm}: 7 vertices, radius {rad:.6f}, all step-{k} edges unit: {edge_ok}")

    print("\nequilateral triangles P_j Q_j R_j:")
    tri = all(fld.is_unit(P[j], Q[j]) and fld.is_unit(P[j], R[j])
              and fld.is_unit(Q[j], R[j]) for j in range(7))
    print(f"  all 7 are unit equilateral: {tri}")

    g = graph_H(fld)
    print(f"\nH: {g.summary()}   (expected 21 vertices)")

    uv = unit_vectors(fld, g)
    print(f"distinct unit edge-directions: {len(uv)}")

    from colour import chromatic_number
    print(f"chromatic number of H: {chromatic_number(g, 2, 6)}")
