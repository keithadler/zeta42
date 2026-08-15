"""Overnight validation: chi(G3) = 5, checked directly.

Fast half: G3 is 5-colourable (seconds, verified independently).
Slow half: G3 is not 4-colourable -- the ~8h-class UNSAT.  Symmetry-broken
encoding (sound: fixing a clique's colours, no other constraints present).

Also exports the 4-colourability CNF to g3_4col.cnf so kissat can be pointed
at it as an independent engine when it frees up.
"""

import time
from fractions import Fraction

from hept_g1 import load_G1, SCALE
from cycloext import get_cycloext
from udg import UDG
from colour import ColourInstance, k_colourable, verify_colouring


def build_g3(log=print):
    fld, core, B, _ = load_G1(log=lambda s: None)
    E = get_cycloext(42, 5)
    lift = lambda v: (tuple(Fraction(x, SCALE) for x in v), E.base.zero)
    g1_pts = [lift(v) for v in sorted(core)]

    c1, c2 = E.zeta(14), E.zeta(7)
    V1 = [E.rotate_mul(p, E.zeta(-7), centre=c1) for p in g1_pts]
    V2 = [E.rotate_mul(p, E.zeta(7), centre=c2) for p in g1_pts]
    seen, g2_pts = set(), []
    for p in V1 + V2:
        if p not in seen:
            seen.add(p)
            g2_pts.append(p)

    b = E.base
    isqrt3 = b.sub(b.zeta(7), b.zeta(35))
    w = (b.rat(7, 8), tuple(x * Fraction(1, 8) for x in isqrt3))
    assert E.norm2(w) == E.one
    centre = E.rat(-1)
    rot = [E.rotate_mul(p, w, centre=centre) for p in g2_pts]
    seen, g3_pts = set(), []
    for p in g2_pts + rot:
        if p not in seen:
            seen.add(p)
            g3_pts.append(p)
    g3 = UDG(E, g3_pts)
    log(f"G3: {g3.summary()}  (expected 2131 / 12530)")
    assert (len(g3), len(g3.edges)) == (2131, 12530)
    return g3


def main():
    t0 = time.time()
    g3 = build_g3(log=lambda m: print(m, flush=True))

    ok5, col5 = k_colourable(g3, 5)
    good = verify_colouring(g3, col5)[0] if col5 else False
    print(f"5-colourable: {ok5} (verified {good})  [{time.time()-t0:.0f}s]",
          flush=True)

    inst = ColourInstance(g3, 4, break_symmetry=True)
    nv = inst.pool.top
    with open("g3_4col.cnf", "w") as fh:
        fh.write(f"p cnf {nv} {len(inst.clauses)}\n")
        for cl in inst.clauses:
            fh.write(" ".join(map(str, cl)) + " 0\n")
    print(f"exported g3_4col.cnf ({nv} vars, {len(inst.clauses)} clauses)",
          flush=True)

    print("solving 4-colourability (expect UNSAT, ~hours)...", flush=True)
    t0 = time.time()
    sat, _ = inst.solve()
    print(f"4-colourable: {sat}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"chi(G3) = 5: {'CONFIRMED' if not sat else 'FAILED -- investigate'}",
          flush=True)


if __name__ == "__main__":
    main()
