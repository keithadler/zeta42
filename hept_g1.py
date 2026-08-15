"""Haugland's G1: the 740-vertex core of the spindle-free construction.

From arXiv:2608.04542, section on the graph G1:

    T5   all vertices lying on a path of <= 5 unit steps from A=0 to B=i*sqrt3
    T6   the same with <= 6 steps
    G0   the vertices of T6 that are in T5, or adjacent to >= 7 vertices of T5
    G1   the 7-core of G0  (repeatedly delete vertices of degree < 7)

and the property that makes it useful: forcing A and B to the *same* colour
makes G1 non-4-colourable, so {A, B} is a forced non-monochromatic pair.

Representation.  The 84 arcs all lie in (1/7)Z[zeta_42], so multiplying by 7
turns every reachable point into a 12-tuple of small integers.  Points are
therefore plain integer tuples here -- fast to hash and light on memory --
and only converted back to field elements when the graph is handed to the
exact unit-distance machinery.

Since the arc set is closed under negation, step-distance is a genuine metric,
so a vertex on a path of length <= n from A to B is exactly a vertex with
d(A,v) + d(v,B) <= n.  And every vertex of such a path has min(d_A, d_B) <= 3
for n <= 6, which is what keeps the enumeration finite: we only ever need
balls of radius 3.
"""

import json
import os
import time

from cyclo import get_cyclofield
from heptagon import field, graph_H, unit_vectors
from udg import UDG

SCALE = 7          # arcs live in (1/7)Z[zeta_42]


def arc_lattice(fld=None):
    """The 84 arcs as integer vectors, plus the scale used."""
    fld = fld or field()
    g = graph_H(fld)
    arcs = unit_vectors(fld, g)
    ints = [tuple(int(c * SCALE) for c in a) for a in arcs]
    assert all(all((c * SCALE).denominator == 1 for c in a) for a in arcs)
    return fld, arcs, ints


def _add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def _sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def balls(steps, radius):
    """{point: min number of steps} for every point within `radius` of 0."""
    zero = (0,) * len(steps[0])
    dist = {zero: 0}
    frontier = [zero]
    for r in range(1, radius + 1):
        nxt = []
        for p in frontier:
            for s in steps:
                q = _add(p, s)
                if q not in dist:
                    dist[q] = r
                    nxt.append(q)
        frontier = nxt
    return dist


def _reachable_within(v, k, ball, steps):
    """Is d(0, v) <= k, given `ball` holds exact distances up to radius R?"""
    d = ball.get(v)
    if d is not None and d <= k:
        return True
    R = max(ball.values())
    if k <= R:
        return False
    # meet in the middle: v = w + (v - w) with d(w) <= k-R and d(v-w) <= R
    need = k - R
    for w, dw in ball.items():
        if dw <= need and _sub(v, w) in ball:
            return True
    return False


def build_T(steps, B, n, radius=3, log=print):
    """Vertices on some path of <= n steps from 0 to B."""
    ball = balls(steps, radius)
    log(f"  ball(<= {radius}) around each endpoint: {len(ball)} points")

    # candidates: everything within `radius` of A=0 or of B
    cands = {}
    for p, d in ball.items():
        cands[p] = min(cands.get(p, 99), d)          # d_A known
    ballB = {_add(B, p): d for p, d in ball.items()}

    out = set()
    for v, dA in ball.items():
        if _reachable_within(_sub(v, B), n - dA, ball, steps):
            out.add(v)
    for v, dB in ballB.items():
        if _reachable_within(v, n - dB, ball, steps):
            out.add(v)
    return out, ball


CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "g1-cache.json")


def load_G1(log=print):
    """G1 from cache if present, else build it and cache the result.

    The enumeration takes ~220s (83581-point balls, twice), which is far too
    slow to sit in front of every experiment.
    """
    fld = field()
    if os.path.exists(CACHE):
        with open(CACHE) as fh:
            d = json.load(fh)
        core = {tuple(v): set() for v in d["vertices"]}
        steps = [tuple(s) for s in d["steps"]]
        for v in core:
            for s in steps:
                w = _add(v, s)
                if w in core:
                    core[v].add(w)
        log(f"G1 from cache: {len(core)} vertices")
        return fld, core, tuple(d["B"]), None
    fld, core, B, Bf = build_G1(log=log)
    _, _, steps = arc_lattice(fld)
    with open(CACHE, "w") as fh:
        json.dump({"vertices": [list(v) for v in sorted(core)],
                   "steps": [list(s) for s in steps],
                   "B": list(B)}, fh)
    log(f"cached to {CACHE}")
    return fld, core, B, Bf


def build_G1(log=print):
    t0 = time.time()
    fld, arcs, steps = arc_lattice()
    log(f"arcs: {len(steps)} unit vectors, integer lattice Z^{len(steps[0])}")

    # B = i*sqrt3 = zeta_6 - zeta_6^5, with zeta_6 = zeta_42^7
    z6, z6_5 = fld.zeta(7), fld.zeta(35)
    Bf = fld.sub(z6, z6_5)
    assert fld.norm2(Bf) == fld.rat(3), "B should satisfy |B|^2 = 3"
    B = tuple(int(c * SCALE) for c in Bf)

    log("T5:")
    T5, ball = build_T(steps, B, 5, log=log)
    log(f"  |T5| = {len(T5)}   (paper: 1042)")

    log("T6:")
    T6, _ = build_T(steps, B, 6, log=log)
    log(f"  |T6| = {len(T6)}   (paper: 12856)")

    # G0: vertices of T6 in T5, or adjacent to >= 7 vertices of T5
    arcset = set(steps)
    T5set = T5
    G0 = set()
    for v in T6:
        if v in T5set:
            G0.add(v)
            continue
        deg = 0
        for s in steps:
            if _add(v, s) in T5set:
                deg += 1
                if deg >= 7:
                    break
        if deg >= 7:
            G0.add(v)
    log(f"G0: {len(G0)} vertices")

    # G1: the 7-core
    adj = {v: set() for v in G0}
    for v in G0:
        for s in steps:
            w = _add(v, s)
            if w in adj:
                adj[v].add(w)
    core = dict(adj)
    changed = True
    while changed:
        changed = False
        for v in [v for v, ns in core.items() if len(ns) < 7]:
            for w in core[v]:
                core[w].discard(v)
            del core[v]
            changed = True
    log(f"G1: {len(core)} vertices, {sum(len(n) for n in core.values())//2} edges"
        f"   (paper: 740 vertices, 3985 edges)")
    log(f"[{time.time()-t0:.1f}s]")
    return fld, core, B, Bf


if __name__ == "__main__":
    build_G1()
