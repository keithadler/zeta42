"""De Grey's J -> K -> L assembly, with the H-copies tracked through it.

Following "The chromatic number of the plane is at least 5" (arXiv:1804.02385):

    H  7 vertices   regular hexagon of side 1 plus centre
    J  31 vertices  13 copies of H: centre, 6 at distance 1, 6 at distance √3
    K  61 vertices  J together with J rotated about the origin by 2·arcsin(1/4)
    L  121 vertices K together with K rotated about A by 2·arcsin(1/8)

The rotation angles are exactly the ones that move a vertex at radius 2 (resp.
radius 4) a distance of 1: cos = 7/8, sin = √15/8 and cos = 31/32,
sin = 3√7/32.  So everything lives in Q(√3, √5, √7).

The property that matters, and that verify_L() checks with a SAT call, is:
in every 4-colouring of L at least one of the 52 copies of H contains a
monochromatic triple (an equilateral triangle of side √3 among its ring).
"""

from field import get_field
from udg import UDG, HEX_RING, hexagon, lattice_pt
from colour import ColourInstance

# lattice coordinates of the 6 points at distance √3 from the origin
SQRT3_RING = ((1, 1), (2, -1), (1, -2), (-1, -1), (-2, 1), (-1, 2))


class Assembly:
    """A point set together with the copies of H it is built from."""

    def __init__(self, fld, points=(), hcopies=()):
        self.fld = fld
        self.points = list(points)
        self.hcopies = list(hcopies)     # each a 7-tuple: centre, then ring

    def rotated(self, cos_t, sin_t, centre=None):
        f = self.fld
        pts = [f.rotate(p, cos_t, sin_t, centre) for p in self.points]
        hcs = [tuple(f.rotate(p, cos_t, sin_t, centre) for p in h)
               for h in self.hcopies]
        return Assembly(f, pts, hcs)

    def union(self, other):
        seen = set()
        pts = []
        for p in list(self.points) + list(other.points):
            if p not in seen:
                seen.add(p)
                pts.append(p)
        hcs = []
        seen_h = set()
        for h in list(self.hcopies) + list(other.hcopies):
            key = frozenset(h)
            if key not in seen_h:
                seen_h.add(key)
                hcs.append(h)
        return Assembly(self.fld, pts, hcs)

    def graph(self):
        return UDG(self.fld, self.points)

    def triples(self, graph):
        """The two √3-triples of each H copy, as vertex-id triples."""
        out = []
        for h in self.hcopies:
            ring = [graph.index_of(p) for p in h[1:]]
            if any(i is None for i in ring):
                raise ValueError("H copy not contained in graph")
            out.append(tuple(ring[i] for i in (0, 2, 4)))
            out.append(tuple(ring[i] for i in (1, 3, 5)))
        return out

    def summary(self):
        return f"{len(self.points)} vertices, {len(self.hcopies)} copies of H"


def h_copy(fld, centre):
    """A copy of H as (centre, r0..r5) in lattice orientation."""
    return tuple(hexagon(fld, centre))


def build_J(fld):
    """13 copies of H: one central, 6 at distance 1, 6 at distance √3."""
    centres = [fld.origin]
    centres += [lattice_pt(fld, a, b) for a, b in HEX_RING]
    centres += [lattice_pt(fld, a, b) for a, b in SQRT3_RING]
    hcs = [h_copy(fld, c) for c in centres]
    pts = []
    seen = set()
    for h in hcs:
        for p in h:
            if p not in seen:
                seen.add(p)
                pts.append(p)
    return Assembly(fld, pts, hcs)


def build_K(fld, J=None):
    """J union J rotated about the origin by 2·arcsin(1/4)."""
    J = J or build_J(fld)
    cos_t, sin_t = fld.rat(7, 8), fld.surd(15, 1, 8)
    return J.union(J.rotated(cos_t, sin_t))


def build_L(fld, K=None):
    """K union K rotated about A = (-2, 0) by 2·arcsin(1/8).

    A is a linking vertex at distance 2 from the centre; its opposite B is at
    distance 4 from A, so this rotation moves B exactly one unit, to B'.
    """
    K = K or build_K(fld)
    A = lattice_pt(fld, -2, 0)
    cos_t, sin_t = fld.rat(31, 32), fld.surd(7, 3, 32)
    return K.union(K.rotated(cos_t, sin_t, centre=A))


def verify_L(fld=None, verbose=True):
    """Check de Grey's key lemma about L with a SAT call.

    Asserts there is NO 4-colouring of L in which every copy of H avoids a
    monochromatic √3-triple.  UNSAT is the result the construction needs.
    """
    fld = fld or get_field((3, 5, 7))
    L = build_L(fld)
    g = L.graph()
    inst = ColourInstance(g, 4)
    triples = L.triples(g)
    for t in triples:
        inst.forbid_monochromatic(t)
    sat, _ = inst.solve()
    if verbose:
        print(f"L: {L.summary()}, {len(g.edges)} edges")
        print(f"   forbidding {len(triples)} monochromatic triples")
        print(f"   4-colourable with no mono triple: {sat}  (expected False)")
    return not sat


if __name__ == "__main__":
    fld = get_field((3, 5, 7))
    J = build_J(fld)
    K = build_K(fld, J)
    L = build_L(fld, K)
    for name, a in (("J", J), ("K", K), ("L", L)):
        g = a.graph()
        print(f"{name}: {a.summary()}, {len(g.edges)} edges")
    print()
    print("lemma holds:", verify_L(fld))
