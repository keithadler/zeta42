"""Rung two of the B-hunt: the G1 template on the full 90-arc D=11 alphabet.

The census showed Q(zeta_5) has exactly 90 unit vectors with denominator 11
(vs 34 at height 1, and vs the 84-over-7 that fueled the level-4 build).
Extraction is vectorised: 4|sum a_i zeta^i|^2 = R(a) + S(a) sqrt5 with R, S
integer quadratic forms, so the norm equation |w|^2 = 121 is R = 484, S = 0.

Then the faithful template: T-sets between the designated phi-pair (0, 1+zeta),
degree filter, 7-core, and the 1000-sample dial on the pair.
"""

import time
from fractions import Fraction

import numpy as np

from pentagon import field
from pent_b import to_lattice, dial_pair, build
from hept_g1 import balls


def d11_arcs(fld):
    b = 22
    rng = np.arange(-b, b + 1)
    A0, A1, A2, A3 = np.meshgrid(rng, rng, rng, rng, indexing="ij")
    coeff = [A0.astype(np.int64), A1.astype(np.int64),
             A2.astype(np.int64), A3.astype(np.int64)]
    R4 = {0: 4, 1: -1, 2: -1, 3: -1, 4: -1}
    S4 = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
    Rv = np.zeros_like(coeff[0])
    Sv = np.zeros_like(coeff[0])
    for i in range(4):
        Rv += 4 * coeff[i] * coeff[i]
        for j in range(i + 1, 4):
            k = (j - i) % 5
            Rv += 2 * R4[k] * coeff[i] * coeff[j]
            Sv += 2 * S4[k] * coeff[i] * coeff[j]
    m = (Rv == 4 * 121) & (Sv == 0)
    vecs = np.stack([c[m] for c in coeff], axis=1)
    arcs = [tuple(Fraction(int(x), 11) for x in v) for v in vecs]
    assert all(fld.norm2(u) == fld.one for u in arcs)
    return arcs


def main():
    fld = field()
    t0 = time.time()
    arcs = d11_arcs(fld)
    print(f"D=11 alphabet: {len(arcs)} arcs  [{time.time()-t0:.0f}s]", flush=True)

    B_el = fld.add(fld.one, fld.zeta(1))
    D, steps, (B_int,) = to_lattice(fld, arcs, extra=[B_el])
    A_int = (0,) * 4
    print(f"lattice denominator: {D}", flush=True)

    t0 = time.time()
    ball = balls(steps, 3)
    print(f"ball(<=3): {len(ball)}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"d(0,B) = {ball.get(B_int, '>3')}", flush=True)

    for n, deg_min, core_k in ((5, 7, 7), (6, 7, 7), (6, 8, 8)):
        print(f"\nparams: n={n}, deg>={deg_min}, {core_k}-core", flush=True)
        t0 = time.time()
        core = build(fld, steps, B_int, n, deg_min, core_k,
                     log=lambda m: print(m, flush=True))
        print(f"  [{time.time()-t0:.0f}s]", flush=True)
        if core and A_int in core and B_int in core:
            dial_pair(fld, core, D, A_int, B_int,
                      log=lambda m: print(m, flush=True))
        elif core:
            print(f"  core nonempty ({len(core)}) but pair fell out", flush=True)


if __name__ == "__main__":
    main()
