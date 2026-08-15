"""Baseline checks.  Run this before trusting any search result.

    ./.venv/bin/python test_baseline.py

Everything here has a known right answer from the literature, so a failure
means the machinery is broken rather than the mathematics being interesting.
"""

import sys

from field import get_field, chord_rotation
from udg import UDG, moser_spindle, hexagon, tri_lattice
from colour import chromatic_number, k_colourable, verify_colouring
from degrey import build_J, build_K, build_L, verify_L

FAILURES = []


def check(label, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}: {got}" + ("" if ok else f"  (expected {want})"))
    if not ok:
        FAILURES.append(label)


print("exact arithmetic")
f357 = get_field((3, 5, 7))
check("√3·√5 == √15", f357.mul(f357.surd(3), f357.surd(5)) == f357.surd(15), True)
check("(√7)² == 7", f357.mul(f357.surd(7), f357.surd(7)) == f357.rat(7), True)
f = get_field((3, 5, 7, 11))
check("chord rotation r²=3  cos", f.fmt(chord_rotation(f, 3)[0]), "5/6")
check("chord rotation r²=3  sin", f.fmt(chord_rotation(f, 3)[1]), "1/6·√11")
check("chord rotation r²=4  cos", f.fmt(chord_rotation(f, 4)[0]), "7/8")
check("chord rotation r²=4  sin", f.fmt(chord_rotation(f, 4)[1]), "1/8·√15")
check("chord rotation r²=16 cos", f.fmt(chord_rotation(f, 16)[0]), "31/32")
check("chord rotation r²=16 sin", f.fmt(chord_rotation(f, 16)[1]), "3/32·√7")

print("\nknown small graphs")
H = UDG(f357, hexagon(f357))
check("H vertices", len(H), 7)
check("H edges", len(H.edges), 12)
check("H chromatic number", chromatic_number(H, 2, 6), 3)

ms = moser_spindle()
check("Moser spindle vertices", len(ms), 7)
check("Moser spindle edges", len(ms.edges), 11)
check("Moser spindle chromatic number", chromatic_number(ms, 2, 6), 4)
check("Moser spindle 3-colourable", k_colourable(ms, 3)[0], False)
_, col = k_colourable(ms, 4)
check("returned 4-colouring is proper", verify_colouring(ms, col)[0], True)

print("\nde Grey assembly (arXiv:1804.02385)")
J = build_J(f357)
K = build_K(f357, J)
L = build_L(f357, K)
check("J vertices", len(J.points), 31)
check("J copies of H", len(J.hcopies), 13)
check("K vertices", len(K.points), 61)
check("K copies of H", len(K.hcopies), 26)
check("L vertices", len(L.points), 121)
check("L copies of H", len(L.hcopies), 52)
check("L lemma: no 4-colouring avoids a mono triple", verify_L(f357, verbose=False), True)

print("\nlattice")
check("|p|²≤4 triangular lattice size", len(tri_lattice(f357, 4)), 19)

print("\nde Grey's G (structure only — colouring check needs --slow)")
from degrey_g import build_G, build_Sa, build_Y, seed_points, PRIMES
fg = get_field(PRIMES)
check("S seed points", len(seed_points(fg)), 39)
check("S_a vertices", len(build_Sa(fg)), 397)
check("Y vertices", len(build_Y(fg)), 791)
G = build_G(fg)
check("G vertices", len(G), 1581)
check("G edges", len(G.edges), 7877)
check("G is 5-colourable", k_colourable(G, 5)[0], True)

if "--slow" in sys.argv:
    # ~7 minutes on an M3: proving UNSAT is the expensive direction
    print("\n  running the slow check (G is not 4-colourable)...")
    import time as _t
    _t0 = _t.time()
    check("G is NOT 4-colourable", k_colourable(G, 4)[0], False)
    print(f"        took {_t.time()-_t0:.0f}s")
else:
    print("  skip  G is NOT 4-colourable  (pass --slow; takes ~7 min)")

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILURE(S): {', '.join(FAILURES)}")
    sys.exit(1)
print("all baseline checks passed")
