"""Spindled zeta_15: the first construction with chromatic twists AND native
golden-ratio pairs in one arithmetic.

Measured background (this session): raw lattice balls are 3-chromatic in every
field; chromatic strength has always come from Moser-style twists (rotation by
w, |w| = 1, cos = 5/6, sin = sqrt(11)/6), which move every point at radius
sqrt(3) by exactly 1; and Q(zeta_15) is the first searched world whose balls
natively carry thousands of exact phi-pairs AND unit triangles.  Here we
apply the twists inside Q(zeta_15, sqrt(-11)) and watch three dials:

    chi of the union            (does the twist push 3 -> 4 here as it did
                                 for the plain triangular lattice?)
    exact phi-pair count        (does the forcing distance survive the union?)
    aggregate phi tests         (k = 4 first; k = 5 is the prize)

Usage:  ./.venv/bin/python spindle15.py [radius] [max_pts]
"""

import sys
import time
from collections import deque

from ext_imag import get_imagext
from udg import UDG
from colour import chromatic_number, ColourInstance, verify_colouring


def build_ball(E, radius, max_pts):
    """BFS ball of the 30 unit zeta_15-steps, in the extension field."""
    steps = [E.zeta(k) for k in range(15)] + [E.neg(E.zeta(k)) for k in range(15)]
    seen = {E.zero}
    q = deque([E.zero])
    while q and len(seen) < max_pts:
        p = q.popleft()
        for st in steps:
            np_ = E.add(p, st)
            if np_ in seen:
                continue
            if abs(E.to_complex(np_)) > radius:
                continue
            seen.add(np_)
            q.append(np_)
    return sorted(seen)


def phi_pairs_exact(E, g):
    """Vertex pairs at exactly distance phi: d^2 = phi + 1."""
    b = E.base
    phi_el = b.add(b.add(b.zeta(3), b.zeta(12)), b.one)      # 2cos(pi/5)... = phi
    phi2_el = E.lift(b.add(phi_el, b.one))                   # phi^2 = phi + 1
    PHI2 = 2.618033988749895
    fl = [E.pfloat(p) for p in g.points]
    out = []
    for i in range(len(g)):
        xi, yi = fl[i]
        for j in range(i + 1, len(g)):
            dx, dy = xi - fl[j][0], yi - fl[j][1]
            if abs(dx * dx + dy * dy - PHI2) > 1e-6:
                continue
            if E.dist2(g.points[i], g.points[j]) == phi2_el:
                out.append((i, j))
    return out


def aggregate_phi(E, g, k, pairs, log=print):
    inst = ColourInstance(g, k)
    for u, v in pairs:
        for c in range(k):
            inst.add_clause([-inst.var(u, c), -inst.var(v, c)])
    t0 = time.time()
    sat, col = inst.solve()
    dt = time.time() - t0
    if sat:
        ok, _ = verify_colouring(g, col)
        assert ok
        log(f"  k={k}: SAT in {dt:.1f}s -- colouring dodges all phi-pairs")
    else:
        log(f"  k={k}: UNSAT in {dt:.1f}s -- AGGREGATE PHI-FORCER at k={k}")
    return not sat


def main():
    radius = float(sys.argv[1]) if len(sys.argv) > 1 else 3.2
    max_pts = int(sys.argv[2]) if len(sys.argv) > 2 else 800

    E = get_imagext(15, 11)
    print(f"field: {E}", flush=True)

    # Moser twist: w = 5/6 + sqrt(-11)/6, |w| = 1 exactly
    w = E.add(E.rat(5, 6), E.sqrtnegD(1, 6))
    assert E.norm2(w) == E.one, "twist is not unimodular"
    winv = E.conj(w)

    t0 = time.time()
    base = build_ball(E, radius, max_pts)
    print(f"base ball: {len(base)} points  [{time.time()-t0:.0f}s]", flush=True)

    # union: base + twist about origin + inverse twist about a lattice point
    c2 = E.zeta(0)                                   # centre (1,0)
    pts = list(base)
    seen = set(base)
    for p in base:
        for q in (E.rotate_mul(p, w),
                  E.rotate_mul(p, winv, centre=c2)):
            if q not in seen:
                seen.add(q)
                pts.append(q)
    t0 = time.time()
    g = UDG(E, pts)
    print(f"spindled union: {g.summary()}  "
          f"mean degree {2*len(g.edges)/len(g):.2f}  [{time.time()-t0:.0f}s]",
          flush=True)

    t0 = time.time()
    chi = chromatic_number(g, 3, 7)
    print(f"chromatic number: {chi}  [{time.time()-t0:.1f}s]", flush=True)

    t0 = time.time()
    pairs = phi_pairs_exact(E, g)
    print(f"exact phi-pairs: {len(pairs)}  [{time.time()-t0:.0f}s]", flush=True)

    if pairs:
        print("aggregate phi tests:", flush=True)
        for k in range(max(chi, 4), 6):
            forced = aggregate_phi(E, g, k, pairs,
                                   log=lambda m: print(m, flush=True))
            if not forced:
                break


if __name__ == "__main__":
    main()
