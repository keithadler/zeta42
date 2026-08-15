"""Unit-distance graphs over exact points.

Edge detection is a hybrid: a float grid narrows the search to candidate pairs
in near-linear time, then every candidate is confirmed with exact field
arithmetic.  Floats only ever *narrow*, never decide, so no spurious edge can
survive; the grid cell size (1) with a full 3x3 halo cannot miss a real one.
"""

from collections import defaultdict

from field import get_field

EPS = 1e-6


class UDG:
    """A set of exact points plus the induced unit-distance edges."""

    def __init__(self, fld, points=()):
        self.fld = fld
        self.points = []
        self._index = {}
        self._edges = None
        for p in points:
            self.add(p)

    # ------------------------------------------------------------ building --
    def add(self, p):
        i = self._index.get(p)
        if i is None:
            i = len(self.points)
            self._index[p] = i
            self.points.append(p)
            self._edges = None
        return i

    def add_all(self, pts):
        return [self.add(p) for p in pts]

    def index_of(self, p):
        return self._index.get(p)

    def __len__(self):
        return len(self.points)

    def __contains__(self, p):
        return p in self._index

    # ------------------------------------------------------------- edges ----
    @property
    def edges(self):
        if self._edges is None:
            self._edges = self._compute_edges()
        return self._edges

    def _compute_edges(self):
        fld = self.fld
        fl = [fld.pfloat(p) for p in self.points]
        cells = defaultdict(list)
        for i, (x, y) in enumerate(fl):
            cells[(int(x // 1.0), int(y // 1.0))].append(i)

        edges = set()
        for (cx, cy), members in cells.items():
            neigh = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    neigh.extend(cells.get((cx + dx, cy + dy), ()))
            for i in members:
                xi, yi = fl[i]
                for j in neigh:
                    if j <= i:
                        continue
                    xj, yj = fl[j]
                    d2 = (xi - xj) ** 2 + (yi - yj) ** 2
                    if abs(d2 - 1.0) > EPS:
                        continue
                    if fld.is_unit(self.points[i], self.points[j]):
                        edges.add((i, j))
        return sorted(edges)

    def adjacency(self):
        adj = [[] for _ in self.points]
        for i, j in self.edges:
            adj[i].append(j)
            adj[j].append(i)
        return adj

    def subgraph(self, ids):
        return UDG(self.fld, [self.points[i] for i in ids])

    def summary(self):
        return f"{len(self)} vertices, {len(self.edges)} edges"


# --------------------------------------------------------------- pieces -----

def tri_lattice(fld, radius2_max):
    """Unit triangular lattice points p with |p|^2 <= radius2_max.

    Basis (1,0), (1/2, √3/2); the point a*u + b*v has |p|^2 = a^2 + ab + b^2.
    """
    s3 = fld.surd(3)
    pts = []
    lim = int(radius2_max ** 0.5) + 2
    for a in range(-2 * lim, 2 * lim + 1):
        for b in range(-2 * lim, 2 * lim + 1):
            if a * a + a * b + b * b > radius2_max:
                continue
            pts.append((fld.add(fld.rat(a), fld.rat(b, 2)),
                        fld.scale(s3, b, 2)))
    return pts


def lattice_pt(fld, a, b):
    """The triangular-lattice point a*(1,0) + b*(1/2, √3/2)."""
    return (fld.add(fld.rat(a), fld.rat(b, 2)), fld.scale(fld.surd(3), b, 2))


HEX_RING = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))


def hexagon(fld, centre=None):
    """H: regular hexagon of side 1 plus its centre -- 7 points, 12 edges."""
    if centre is None:
        centre = fld.origin
    return [centre] + [fld.padd(centre, lattice_pt(fld, a, b)) for a, b in HEX_RING]


def moser_spindle(fld=None):
    """The 7-vertex, 11-edge Moser spindle.  chi = 4.

    Two unit rhombi sharing the apex at the origin; the second is the first
    rotated about the apex so the far tips (at radius √3) land unit distance
    apart -- cos t = 5/6, sin t = √11/6.
    """
    if fld is None:
        fld = get_field((3, 11))
    s3 = fld.surd(3)
    a = fld.origin
    b = (fld.rat(1, 2), fld.scale(s3, 1, 2))
    c = (fld.rat(-1, 2), fld.scale(s3, 1, 2))
    d = (fld.zero, s3)
    cos_t, sin_t = fld.rat(5, 6), fld.surd(11, 1, 6)
    rot = [fld.rotate(p, cos_t, sin_t) for p in (b, c, d)]
    return UDG(fld, [a, b, c, d] + rot)
