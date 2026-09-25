"""The ninefold analogue of Haugland's heptagonal seed, in Q(zeta_9) = Q(zeta_18).

Question: `heptagon.py` places three unit-edge star polygons {7/1}, {7/2},
{7/3} in Q(zeta_42), phased (by exhaustive search over roots of unity) so that
P_j, Q_j, R_j are unit equilateral triangles; the 42 non-root edge directions
over denominator 7 are the fuel of the whole G1/G3 construction.  Does the
same approach work with nine-fold symmetry?

The field is smaller than for seven: zeta_6 = -zeta_9^6 already lives in
Q(zeta_9), so equilateral triangles cost nothing and the ambient field is
Q(zeta_18) = Q(zeta_9), of degree 6 (vs 12).  The unit-edge star {9/k} is,
exactly as for 7,

    V_j = zeta_9^j / (zeta_9^k - 1),     V_(j+k) - V_j = zeta_9^j.

Answer (all exact except (1), which is a continuous check):

1. There is no ninefold Haugland seed.  Forget arithmetic and rotate the
   stars freely: a unit equilateral triangle with one vertex on each of three
   {9/k} circumcircles does not exist for any of the four star triples; the
   closest miss is 0.068 in radius ({9/2},{9/3},{9/4}).  The same test
   returns an exact solution for n = 7 (the control), and none for n = 11, 13.
   Haugland's phasing is a heptagonal coincidence, not a family.

2. The exact phase search over zeta_18 finds one best configuration (and its
   mirror image): {9/1}, {9/2}, {9/4} with the middle star rotated by pi/9.
   It has 27 vertices, 81 edges, is 6-regular and 3-chromatic, and all 18 of
   its edge directions are roots of unity: it is a patch of the Z[zeta_18]
   lattice.  The seed hands the template no new arcs.

3. Unit-arc census (v / conj v, v in Z[zeta_9] of height <= 3, closed under
   the 18 roots of unity):  denominator 1: 18, 3: 0, 7: 36, 13: 36,
   19: 468, 37: 324, ...  Control: the same method recovers Q(zeta_5)'s
   80 + 10 = 90 arcs over 11.  Arithmetic reading: an arc over p is
   pi / conj(pi) for a prime pi | p that conjugation moves.  In Q(zeta_42) the
   prime 7 is totally ramified in Q(zeta_7) *and* split in Q(sqrt -3), which
   is what makes the denominator-7 arcs cheap and plentiful.  In Q(zeta_9) the
   only ramified prime is 3, which also ramifies in Q(sqrt -3): conjugation
   fixes its prime, so 3 contributes no arcs at all.  Primes 7 and 13 are
   inert-then-split (two primes of degree 3) -> 2 x 18 = 36 each, and 19
   (= 1 mod 9) splits completely -> (3^3 - 1) x 18 = 468.

4. The G1 template (`hept_g1.py`: T_n = vertices on arc-paths of <= n steps
   between a designated pair, G0 = T_(n-1) plus the T_n vertices with >= k
   neighbours in it, then the k-core) on every small-denominator alphabet,
   for the pair 0 -- i*sqrt3 (Haugland's) and 0 -- 2 (the distance-2
   campaign's).  Directed queries: can the pair share a colour at k = 4, 5?

     alphabet (arcs)     pair     n  core   T_n     core size     chi  forced?
     roots only (18)     i*sqrt3  6  7      438     empty
     1+7 (54)            i*sqrt3  6  6,7    5478    empty
     1+7 (54)            2        6  6      4425    369v/1378e    3    no/no
     1+7 (54)            2        6  7      4425    empty
     1+7+13 (90)         i*sqrt3  6  5      15642   776v/2945e    3    no/no
     1+7+13 (90)         i*sqrt3  6  6,7    15642   empty
     1+7+13 (90)         2        6  6      12289   621v/2314e    3    no/no
     1+7+13 (90)         2        6  7      12289   empty
     1+19 (486)          i*sqrt3  5  7..14  8726    empty (G0 = T4)
     1+19 (486)          2        5  7..14  7559    empty (G0 = T4)

   Haugland's 84-arc heptagonal alphabet at the same depth: T5 = 1042,
   T6 = 12856, 7-core 740v/3985e, 5-chromatic, pair forced at k = 4.
   Every surviving Q(zeta_9) core is 3-chromatic with its pair free at
   k = 4 and 5.  The 486-arc system is the sharpest case: nominally six times
   the heptagonal alphabet, yet no T5 vertex has even 7 neighbours in T4.

Verdict: the approach that produces the heptagonal world has no ninefold
solution.  The seed does not exist; the best nonagonal substitute adds no
arcs; and the arcs Q(zeta_9) does have (over 7, 13, 19) build path lattices
that stay 3-chromatic under the template.  This is the same fate as
Q(zeta_5) (`pent_b.py`, `pent_b2.py`), and point 3 suggests the common
cause for the open question in note.md section 4: the denominator-7 system
of Q(zeta_42) is fed by a prime that is totally ramified in the n-part *and*
split in the CM direction.  For n = 9 no prime is both.

Usage:
    ./.venv/bin/python nonagon.py              # (1)-(3), ~2 min at height 3
    ./.venv/bin/python nonagon.py --template   # (4), ~15 min
"""

import cmath
import itertools
import math
import sys
import time
from collections import Counter
from fractions import Fraction

from cyclo import get_cyclofield
from udg import UDG

N = 18          # lcm(9, 6): ninefold symmetry and equilateral triangles


def field():
    return get_cyclofield(N)


def star(fld, k, phase=0):
    """Unit-edge star polygon {9/k}, rotated by zeta_18^phase."""
    base = fld.inv(fld.sub(fld.zeta(2 * k), fld.one))       # zeta_18^2 = zeta_9
    if phase:
        base = fld.mul(base, fld.zeta(phase))
    return [fld.mul(fld.zeta(2 * j), base) for j in range(9)]


def denominator(u):
    return math.lcm(*(c.denominator for c in u))


def roots_of_unity(fld):
    return {fld.zeta(k) for k in range(N)}


# ------------------------------------------------------- (1) continuous ----

def triangle_miss(n, ks):
    """Smallest | |c| - r3 | over unit equilateral triangles a, b, c with
    |a| = r1, |b| = r2 -- zero iff Haugland's phasing exists for {n/ks}."""
    best = math.inf
    for r1, r2, r3 in itertools.permutations(
            [1 / (2 * math.sin(k * math.pi / n)) for k in ks]):
        ct = (r2 * r2 - r1 * r1 - 1) / (2 * r1)       # b = r1 + e^{it}
        if abs(ct) > 1:
            continue
        for t in (math.acos(ct), -math.acos(ct)):
            for s in (1, -1):
                c = r1 + cmath.exp(1j * (t + s * math.pi / 3))
                best = min(best, abs(abs(c) - r3))
    return best


# ---------------------------------------------------- (2) exact search -----

def phase_search(fld):
    """All three-star unions over zeta_18 phases, ranked by edge count.

    star(k, p + 2) is star(k, p) relabelled (zeta_18^2 = zeta_9 permutes the
    vertices), so only phase parity matters: 2^3 phasings per star triple.
    """
    out = []
    for ks in itertools.combinations((1, 2, 3, 4), 3):
        for ph in itertools.product((0, 1), repeat=3):
            pts = [p for k, q in zip(ks, ph) for p in star(fld, k, q)]
            if len(set(pts)) == len(pts):
                out.append((len(UDG(fld, pts).edges), ks, ph))
    return sorted(out, reverse=True)


def graph_H9(fld=None, mirror=False):
    """The best ninefold seed: {9/1}, {9/2} (turned by pi/9), {9/4}."""
    fld = fld or field()
    ph = (1, 0, 1) if mirror else (0, 1, 0)
    return UDG(fld, [p for k, q in zip((1, 2, 4), ph) for p in star(fld, k, q)])


def edge_directions(fld, g):
    out = set()
    for i, j in g.edges:
        d = fld.sub(g.points[j], g.points[i])
        out.add(d)
        out.add(fld.neg(d))
    return out


# ------------------------------------------------------- (3) arc census ----

def arc_census(fld, height=3):
    """{denominator: arcs}, from v / conj(v) over v of coefficient height
    <= `height`, closed under the roots of unity."""
    roots = roots_of_unity(fld)
    found = set()
    for c in itertools.product(range(-height, height + 1), repeat=fld.dim):
        if any(c):
            v = tuple(map(Fraction, c))
            found.add(fld.div(v, fld.conj(v)))
    closed = {fld.mul(u, r) for u in found for r in roots}
    by = {}
    for u in closed:
        by.setdefault(denominator(u), []).append(u)
    return {D: sorted(v) for D, v in sorted(by.items())}


def main():
    print("(1) Haugland's triangle phasing, continuous rotations allowed:")
    for n in (7, 9):
        for ks in itertools.combinations(range(1, n // 2 + 1), 3):
            m = triangle_miss(n, ks)
            print(f"  n={n} {ks}: miss {m:.3e}"
                  + ("   <-- exists" if m < 1e-9 else ""))

    fld = field()
    print(f"\n(2) exact phase search in {fld}:")
    for e, ks, ph in phase_search(fld)[:4]:
        print(f"  stars {ks} phase parities {ph}: {e} edges")
    g = graph_H9(fld)
    from colour import chromatic_number
    dirs = edge_directions(fld, g)
    degs = Counter(len(a) for a in g.adjacency())
    print(f"  H9: {g.summary()}, degrees {dict(degs)}, "
          f"chi {chromatic_number(g, 2, 6)}")
    print(f"  edge directions: {len(dirs)}, of which roots of unity: "
          f"{len(dirs & roots_of_unity(fld))}")

    print("\n(3) unit-arc census, height <= 3:")
    t0 = time.time()
    census = arc_census(fld, 3)
    for D in list(census)[:7]:
        print(f"  denominator {D}: {len(census[D])}")
    print(f"  denominator 3: {len(census.get(3, []))}   [{time.time()-t0:.0f}s]")


# ---------------------------------------------------- (4) G1 template ------

def _add(u, v):
    return tuple(x + y for x, y in zip(u, v))


def _sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def path_sets(steps, B, ns, radius):
    """{n: T_n} for the pair 0 -- B: points v with d(0,v) + d(v,B) <= n.

    Exact distances are tabulated up to `radius`; beyond it d(x) <= k is
    decided by peeling steps (x - s, k - 1), which is cheap because every
    point on a short path lies within `radius` of one endpoint.  (The
    meet-in-the-middle of `hept_g1.build_T` is quadratic in the ball and too
    slow for the 486-arc alphabet.)  Needs max(ns) <= 2 * radius + 1.
    """
    from hept_g1 import balls
    assert max(ns) <= 2 * radius + 1
    ball = balls(steps, radius)

    def le(x, k):
        d = ball.get(x)
        if d is not None:
            return d <= k
        if k <= radius:
            return False
        return any(le(_sub(x, s), k - 1) for s in steps)

    out = {}
    for n in ns:
        T = set()
        for v, d in ball.items():
            if d <= n:
                if le(_sub(v, B), n - d):
                    T.add(v)
                w = _add(v, B)
                if le(w, n - d):
                    T.add(w)
        out[n] = T
    return out


def kcore(steps, verts, k):
    adj = {v: {w for w in (_add(v, s) for s in steps) if w in verts}
           for v in verts}
    changed = True
    while changed:
        changed = False
        for v in [v for v, ns in adj.items() if len(ns) < k]:
            for w in adj[v]:
                adj[w].discard(v)
            del adj[v]
            changed = True
    return adj


def template(fld, arcs, B_el, n, ks, radius, log=print):
    """Haugland's G1 recipe on `arcs` between 0 and B_el, for each core k."""
    from colour import chromatic_number
    from rigidity import directed_query

    S = math.lcm(*(denominator(a) for a in arcs + [B_el]))
    steps = [tuple(int(c * S) for c in a) for a in arcs]
    B, A = tuple(int(c * S) for c in B_el), (0,) * fld.dim
    T = path_sets(steps, B, (n - 1, n), radius)
    small, big = T[n - 1], T[n]
    for k in ks:
        G0 = {v for v in big if v in small
              or sum(_add(v, s) in small for s in steps) >= k}
        core = kcore(steps, G0, k)
        line = (f"  n={n} {k}-core: T{n-1}={len(small)} T{n}={len(big)} "
                f"G0={len(G0)} core={len(core)}")
        if A in core and B in core:
            g = UDG(fld, [tuple(Fraction(c, S) for c in v) for v in sorted(core)])
            ai = g.index_of(tuple(Fraction(c, S) for c in A))
            bi = g.index_of(tuple(Fraction(c, S) for c in B))
            line += f" ({len(g.edges)} edges) chi={chromatic_number(g, 2, 6)}"
            for kk in (4, 5):
                verdict, _ = directed_query(g, kk, ai, bi, budget_s=300)
                line += f"  k={kk} pair-mono {verdict}"
        elif core:
            line += "  (pair fell out)"
        log(line)


def template_main():
    fld = field()
    t0 = time.time()
    census = arc_census(fld, 3)
    roots = census[1]
    log = lambda m: print(m, flush=True)
    log(f"census [{time.time()-t0:.0f}s]")
    pairs = (("i*sqrt3", fld.sub(fld.zeta(3), fld.zeta(15))),
             ("2", fld.rat(2)))
    runs = (("roots only", roots, 6, (7,), 3),
            ("1+7", roots + census[7], 6, (6, 7), 3),
            ("1+7+13", roots + census[7] + census[13], 6, (5, 6, 7), 3),
            ("1+19", roots + census[19], 5, (7, 10, 14), 2))
    for name, arcs, n, ks, radius in runs:
        for pname, B_el in pairs:
            log(f"\n{name} ({len(arcs)} arcs), pair 0 -- {pname}")
            template(fld, arcs, B_el, n, ks, radius, log=log)


if __name__ == "__main__":
    if "--template" in sys.argv:
        template_main()
    else:
        main()
