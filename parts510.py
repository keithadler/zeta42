"""Load Parts' 510-vertex 5-chromatic graph from Heule's CNP-SAT repository,
exactly, and cross-validate it against our own unit-distance machinery.

Data: external/CNP-SAT/vtx/510.vtx (Mathematica-syntax exact coordinates) and
edge/510.edge (DIMACS edge list, 1-indexed).  The radicals occurring are
Sqrt[3,5,11,15,33,55,165] -- all products of {3,5,11} -- so every coordinate
lies in Q(sqrt3, sqrt5, sqrt11) and fits `field.Field((3,5,11))` directly.

The validation that matters: we parse only the *vertices*, recompute the edge
set with our exact is_unit test, and demand it equal Heule's published edge
list exactly.  Agreement validates our parser, our edge detection, and their
data in one shot; any mismatch names the offending pair.

Same trick works for 517/529/553 by passing a different basename.
"""

import os
import re

import sympy
from sympy.parsing.mathematica import parse_mathematica

from field import get_field
from udg import UDG

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "external", "CNP-SAT")
PRIMES = (3, 5, 11)


def _to_field(fld, expr):
    """A sympy expression -> an exact element of Q(sqrt3, sqrt5, sqrt11)."""
    expr = sympy.expand(sympy.radsimp(sympy.sqrtdenest(expr)))
    if any(f.is_Pow and not f.base.is_Rational
           for t in (expr.as_ordered_terms() if expr.is_Add else [expr])
           for f in t.as_ordered_factors()):
        # a nested radical survived one pass; denest again after expansion
        expr = sympy.expand(sympy.sqrtdenest(sympy.radsimp(expr)))
    out = fld.zero
    terms = expr.as_ordered_terms() if expr.is_Add else [expr]
    for t in terms:
        coeff, rad = sympy.Rational(1), 1
        for f in t.as_ordered_factors():
            if f.is_Rational:
                coeff *= f
            elif f.is_Pow and abs(f.exp) == sympy.Rational(1, 2):
                base = f.base if f.exp > 0 else 1 / f.base
                if not base.is_Rational or base <= 0:
                    raise ValueError(f"unexpected radical {f!r} in {expr!r}")
                p, q = base.p, base.q             # sqrt(p/q) = sqrt(p*q)/q
                rad *= p * q
                coeff /= q
            else:
                raise ValueError(f"unexpected factor {f!r} in {expr!r}")
        # fold square factors of rad into coeff so rad is square-free
        d = 2
        while d * d <= rad:
            while rad % (d * d) == 0:
                rad //= d * d
                coeff *= d
            d += 1
        piece = (fld.rat(coeff.p, coeff.q) if rad == 1
                 else fld.surd(rad, coeff.p, coeff.q))
        out = fld.add(out, piece)
    return out


def load(basename="510"):
    """Returns (field, UDG, published_edges)."""
    fld = get_field(PRIMES)
    pts = []
    with open(os.path.join(REPO, "vtx", f"{basename}.vtx")) as fh:
        for line in fh:
            line = line.strip().rstrip(",")
            if not line:
                continue
            m = re.match(r"^\{(.*)\}$", line)
            body = m.group(1)
            # split on the top-level comma
            depth, cut = 0, None
            for i, ch in enumerate(body):
                if ch in "([{":
                    depth += 1
                elif ch in ")]}":
                    depth -= 1
                elif ch == "," and depth == 0:
                    cut = i
                    break
            xs, ys = body[:cut], body[cut + 1:]
            pts.append((_to_field(fld, parse_mathematica(xs)),
                        _to_field(fld, parse_mathematica(ys))))

    edges = set()
    with open(os.path.join(REPO, "edge", f"{basename}.edge")) as fh:
        for line in fh:
            if line.startswith("e "):
                _, a, b = line.split()
                i, j = int(a) - 1, int(b) - 1
                edges.add((min(i, j), max(i, j)))

    return fld, UDG(fld, pts), edges


def validate(basename="510", log=print):
    fld, g, published = load(basename)
    log(f"{basename}: parsed {len(g)} exact vertices")
    ours = set(g.edges)
    log(f"  our exact edge detection: {len(ours)} edges")
    log(f"  published edge list:      {len(published)} edges")
    missing = published - ours
    extra = ours - published
    if missing:
        log(f"  MISSING from ours: {sorted(missing)[:5]} ...")
    if extra:
        log(f"  EXTRA in ours:     {sorted(extra)[:5]} ...")
    ok = not missing and not extra
    log(f"  edge sets identical: {ok}")
    return fld, g, ok


if __name__ == "__main__":
    import time
    t0 = time.time()
    fld, g, ok = validate("510")
    print(f"[{time.time()-t0:.1f}s]")
    if ok:
        from colour import k_colourable, verify_colouring
        t0 = time.time()
        ok5, col5 = k_colourable(g, 5)
        good = verify_colouring(g, col5)[0] if col5 else False
        print(f"5-colourable: {ok5} (verified {good})  [{time.time()-t0:.1f}s]")
        print("(non-4-colourability is Heule's published, DRAT-checked result;"
              " rerun with --slow logic later if we want our own proof)")
