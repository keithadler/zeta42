"""Aggregate mono-distance forcing at k=5: the two-distance lift.

The question, per graph G: can G be properly 5-coloured while NO pair of
vertices at distance sqrt(3) shares a colour?

  UNSAT -> every 5-colouring of G has some mono sqrt3-pair: G is the level-5
           aggregate forcing gadget, the exact analogue of de Grey's L one
           level up, and the hard half of a 6-chromatic construction
           (sqrt3 spindles at cos 5/6, sin sqrt11/6 finish it, in-field).
  SAT   -> G is floppy in the aggregate sense too; the solve time and the
           survivor colouring are the data.

This is NOT the specific-pair sweep: that asked whether one named pair is
forced (none is, measured).  Here the constraint is collective -- some pair,
anywhere.  Validation of the instrument: verify_L() in degrey.py is this test
at level 4 and confirms de Grey's published property in 0.09s.
"""

import sys
import time

from colour import ColourInstance, verify_colouring


def dist_pairs(g, d2_int):
    """All vertex pairs at exact squared distance d2_int (field-exact)."""
    fld = g.fld
    target = fld.rat(d2_int) if hasattr(fld, "rat") else None
    pts = g.points
    n = len(pts)
    # float prescreen, exact confirm -- same discipline as edge detection
    import math
    fl = [fld.pfloat(p) for p in pts]
    out = []
    d = math.sqrt(d2_int)
    for i in range(n):
        xi, yi = fl[i]
        for j in range(i + 1, n):
            dx, dy = xi - fl[j][0], yi - fl[j][1]
            if abs(dx * dx + dy * dy - d2_int) > 1e-6:
                continue
            if fld.dist2(pts[i], pts[j]) == target:
                out.append((i, j))
    return out


def aggregate_test(g, k=5, d2=3, log=print):
    pairs = dist_pairs(g, d2)
    log(f"  sqrt({d2})-pairs: {len(pairs)}")
    if not pairs:
        log("  vacuous (no such pairs) -- skip")
        return None
    inst = ColourInstance(g, k)
    for u, v in pairs:
        for c in range(k):
            inst.add_clause([-inst.var(u, c), -inst.var(v, c)])
    t0 = time.time()
    sat, col = inst.solve()
    dt = time.time() - t0
    if sat:
        ok, _ = verify_colouring(g, col)
        mono = sum(1 for u, v in pairs if col[u] == col[v])
        assert ok and mono == 0
        log(f"  SAT in {dt:.1f}s -- 5-colourable avoiding all mono-sqrt{d2} "
            f"pairs (verified)")
    else:
        log(f"  UNSAT in {dt:.1f}s -- AGGREGATE FORCING GADGET FOUND")
    return sat


def main():
    from parts510 import load

    shelf = []
    for base in ("510", "517", "529", "553"):
        _, g, _ = load(base)
        shelf.append((base, g))
    from degrey_g import build_G
    shelf.append(("deGreyG-1581", build_G()))

    for name, g in shelf:
        print(f"{name}: {g.summary()}", flush=True)
        aggregate_test(g, log=lambda m: print(m, flush=True))
        print(flush=True)


if __name__ == "__main__":
    main()
