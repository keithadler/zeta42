"""Rung one of the B-hunt: Haugland's G1 template in the pentagonal world.

Target object B: a unit-distance graph with a designated pair at distance
phi = (sqrt5+1)/2 that is non-monochromatic in every 5-colouring.  (Substitute
B into the 28 phi-edges of Parts' G31 and the result is a 6-chromatic
unit-distance graph.  B is the only missing piece of that program.)

Template (arXiv:2608.04542, one level up): vertices on short arc-paths
between the designated pair, then a density filter and a k-core.  Fuel is the
height-1 unit-arc alphabet of Q(zeta_5) -- every w/conj(w) with w of
coefficient height <= 1 -- symmetrised under negation so path distance is a
metric.  The designated pair is (0, 1 + zeta), which sits at distance phi
exactly (|1 + zeta|^2 = phi + 1 = phi^2).

Gate: the built graph's designated pair gets the 1000-sample dial.  Raw
Minkowski powers floor at ~77/1000 mono; the template earns escalation to the
height-2 alphabet (154 arcs) only if it lands materially below that.
"""

import time
from fractions import Fraction
from itertools import product
from math import lcm

import numpy as np

from pentagon import field, phi2
from hept_g1 import balls, build_T, _add
from udg import UDG
from rigidity import sample_colourings, directed_query


def height1_arcs(fld):
    """All w/conj(w), w of height <= 1, symmetrised under negation."""
    arcs = set()
    for coeffs in product((-1, 0, 1), repeat=4):
        if not any(coeffs):
            continue
        w = tuple(map(Fraction, coeffs))
        u = fld.div(w, fld.conj(w))
        arcs.add(u)
        arcs.add(fld.neg(u))
    assert all(fld.norm2(u) == fld.one for u in arcs)
    return sorted(arcs)


def to_lattice(fld, arcs, extra=()):
    """Scale arcs (and extra points) to a common integer lattice."""
    D = 1
    for el in list(arcs) + list(extra):
        for c in el:
            D = lcm(D, c.denominator)
    vec = lambda el: tuple(int(c * D) for c in el)
    return D, [vec(a) for a in arcs], [vec(e) for e in extra]


def build(fld, steps, B_int, n, deg_min, core_k, log=print):
    """T_n between the pair, degree filter against T_(n-1), then a core."""
    Tsmall, _ = build_T(steps, B_int, n - 1, radius=3, log=lambda s: None)
    Tbig, _ = build_T(steps, B_int, n, radius=3, log=lambda s: None)
    log(f"  T{n-1}: {len(Tsmall)}   T{n}: {len(Tbig)}")

    G0 = set()
    for v in Tbig:
        if v in Tsmall:
            G0.add(v)
            continue
        deg = sum(1 for s in steps if _add(v, s) in Tsmall)
        if deg >= deg_min:
            G0.add(v)

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
        for v in [v for v, ns in core.items() if len(ns) < core_k]:
            for w in core[v]:
                core[w].discard(v)
            del core[v]
            changed = True
    log(f"  G0: {len(G0)} -> {core_k}-core: {len(core)}")
    return core


def dial_pair(fld, core, D, A_int, B_int, n_samples=1000, log=print):
    if A_int not in core or B_int not in core:
        log("  designated pair fell out of the core -- dead parameter set")
        return None
    pts = sorted(core)
    tofield = lambda v: tuple(Fraction(x, D) for x in v)
    g = UDG(fld, [tofield(v) for v in pts])
    ai, bi = g.index_of(tofield(A_int)), g.index_of(tofield(B_int))
    assert fld.dist2(g.points[ai], g.points[bi]) == phi2(fld)
    log(f"  graph: {g.summary()}")
    cols = sample_colourings(g, 5, n_samples, log=lambda s: None)
    a = np.array(cols)
    f = int((a[:, ai] == a[:, bi]).sum())
    log(f"  designated phi-pair mono: {f}/{len(cols)}   "
        f"(raw-power floor ~77/1000)")
    if f == 0:
        verdict, dt = directed_query(g, 5, ai, bi, budget_s=120)
        log(f"  never mono in samples -> directed query: {verdict.upper()} [{dt:.1f}s]")
    return f, g


def main():
    fld = field()
    arcs = height1_arcs(fld)
    print(f"height-1 alphabet, symmetrised: {len(arcs)} arcs", flush=True)

    B_el = fld.add(fld.one, fld.zeta(1))          # |B| = phi exactly
    D, steps, (B_int,) = to_lattice(fld, arcs, extra=[B_el])
    A_int = (0,) * 4
    print(f"lattice denominator: {D}", flush=True)

    t0 = time.time()
    ball = balls(steps, 3)
    print(f"ball(<=3): {len(ball)} points  [{time.time()-t0:.0f}s]", flush=True)
    reach = B_int in ball
    print(f"B reachable within 3 steps: {reach} "
          f"(d = {ball.get(B_int, '>3')})", flush=True)

    for n, deg_min, core_k in ((5, 6, 6), (5, 7, 7), (6, 7, 7), (6, 8, 8)):
        print(f"\nparams: paths <= {n}, degree filter >= {deg_min}, "
              f"{core_k}-core", flush=True)
        core = build(fld, steps, B_int, n, deg_min, core_k,
                     log=lambda m: print(m, flush=True))
        if core:
            dial_pair(fld, core, D, A_int, B_int,
                      log=lambda m: print(m, flush=True))


if __name__ == "__main__":
    main()
