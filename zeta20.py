"""Q(zeta_20): a test of the ramified-and-split hypothesis, with golden-ratio pairs.

Hypothesis (from `nonagon.py`, point 3): the denominator-7 arc system that
makes Q(zeta_42) dense is fed by a prime that is totally ramified in the
n-part (7 in Q(zeta_7)) and split in the CM direction (7 = (2+sqrt-3)(2-sqrt-3)).
Q(zeta_9) and Q(zeta_5) have no such prime and their template cores are empty
or 3-chromatic.  Q(zeta_20) does: 5 is totally ramified in Q(zeta_5) (e = 4)
and split in Q(i) (5 = (2+i)(2-i)).  It also contains phi = |1 + zeta_5|, the
distance the B-hunt needs, but NOT zeta_3 -- it has no unit equilateral
triangles at all, which the heptagonal world has in abundance.

Pass/fail was fixed before running: PASS = a non-empty template 7-core that
keeps the designated pair and needs >= 4 colours; FAIL = empty or 3-chromatic
cores (the Q(zeta_9) / Q(zeta_5) outcome).

Results.

1. Arcs over 5, exhaustive (|w|^2 = 25 in every embedding, w in Z[zeta_20]
   with coefficients in [-5,5], stable at [-7,7]): 160 with exact
   denominator 5, i.e. 20 roots x 8 = 20 x {(pi/conj pi)^a : a = +-1..+-4},
   exactly the ramification count.  By depth a (least t with
   (1 - zeta_5)^t u integral): 40 each at a = 1, 2, 3, 4.  The doubly
   ramified prime 2 contributes 0, as 3 does in Q(zeta_9).

2. The G1 template (T-sets <= 6 steps, degree filter, k-core; `nonagon.template`):

     alphabet (arcs)   pair  k-core   T5     T6      core            chi  forced
     roots (20)        phi   5,6,7    166    710     empty
     roots (20)        2     5        71     571     31v/80e         3    no
     depth<=1 (60)     phi   5        1558   9968    3978v/29486e    4    no
                             6                       2944v/21456e    4    no
                             7                       2132v/14986e    4    no
     depth<=1 (60)     2     5        475    6019    833v/4419e      4    no
                             7                       213v/892e       4    no
     depth<=2 (100)    phi   5        3144   31016   8662v/67916e    4    no
                             7                       4436v/32506e    4    no
     depth<=2 (100)    2     7        779    17399   213v/892e       4    no
     all (180)         phi   7        5348   91564   4580v/33678e    4    no
                             8                       2980v/20966e    4    no
                             10                      empty
     all (180)         2     5        1131   55099   1429v/7299e     4    no
                             7                       213v/892e       4    no
                             8                       empty
     ("forced" = pair cannot share a colour at k = 4 and at k = 5; every
      entry answered SAT, i.e. not forced, at both.)

   Q(zeta_42) at the same depth, for scale: 84 arcs, T5 = 1042, 7-core
   740v/3985e, chi 5, pair forced at k = 4.

Verdict: the hypothesis PASSES as stated -- the first non-heptagonal world in
this repo whose template 7-cores survive with the pair inside and are
4-chromatic, and the depth-1 arcs alone already give it.  It does not reach
forcing: every core stays at chi = 4 and every pair is free at k = 4, 5,
where the heptagonal G1 is 5-chromatic and forced.  The missing triangles are
the obvious suspect for that ceiling; Q(zeta_60) (adds zeta_3) is the test.

Usage:
    ./.venv/bin/python zeta20.py            # (1), ~1 min
    ./.venv/bin/python zeta20.py --template # (2), a few hours (chi proofs)
"""

import itertools
import sys
import time
from collections import Counter
from fractions import Fraction

import numpy as np

from cyclo import get_cyclofield
from nonagon import denominator, template

N = 20


def field():
    return get_cyclofield(N)


def arcs_over(fld, D, bound):
    """Every unit vector w/D, w in Z[zeta_n] with coefficients in [-bound, bound].

    |w|^2 = D^2 is rational, so it holds in every complex embedding: a float
    screen over one embedding per conjugate pair, then an exact check.
    """
    n, d = fld.n, fld.dim
    js = [j for j in range(1, n // 2) if np.gcd(j, n) == 1]
    Z = np.array([[np.exp(2j * np.pi * j * k / n) for k in range(d)] for j in js])
    rng = np.arange(-bound, bound + 1)
    tail = np.array(list(itertools.product(rng, repeat=d - 2)), dtype=np.int64)
    Zt = tail @ Z[:, 2:].T
    out = set()
    for a0 in rng:
        for a1 in rng:
            n2 = np.abs(Zt + (a0 * Z[:, 0] + a1 * Z[:, 1])[None, :]) ** 2
            for row in tail[np.all(np.abs(n2 - D * D) < 1e-6, axis=1)]:
                w = tuple(Fraction(int(x)) for x in (a0, a1, *row))
                if fld.norm2(w) == fld.rat(D * D):
                    out.add(tuple(x / D for x in w))
    return sorted(out)


def depth(fld, u):
    """Least t with (1 - zeta_5)^t * u integral: the ramification exponent."""
    pi = fld.sub(fld.one, fld.zeta(N // 5))
    t = 0
    while denominator(u) > 1:
        u = fld.mul(u, pi)
        t += 1
    return t


def main():
    fld = field()
    t0 = time.time()
    arcs = arcs_over(fld, 5, 5)
    dep = {u: depth(fld, u) for u in arcs}
    print(f"{fld}: {len(arcs)} unit vectors w/5  [{time.time()-t0:.0f}s]")
    print(f"  exact denominator 5: {sum(denominator(u) == 5 for u in arcs)}")
    print(f"  by ramification depth: {sorted(Counter(dep.values()).items())}")
    print(f"  w/2 with exact denominator 2: "
          f"{sum(denominator(u) == 2 for u in arcs_over(fld, 2, 3))}")
    if "--template" not in sys.argv:
        return

    log = lambda m: print(m, flush=True)
    roots = [u for u in arcs if dep[u] == 0]
    pairs = (("phi", fld.add(fld.one, fld.zeta(N // 5))), ("2", fld.rat(2)))
    runs = (("roots", roots, (5, 6, 7)),
            ("depth<=1", [u for u in arcs if dep[u] <= 1], (5, 6, 7)),
            ("depth<=2", [u for u in arcs if dep[u] <= 2], (5, 6, 7)),
            ("all", arcs, (5, 6, 7, 8, 10)))
    for name, al, ks in runs:
        for pname, B_el in pairs:
            log(f"\n{name} ({len(al)} arcs), pair 0 -- {pname}")
            template(fld, al, B_el, 6, ks, 3, log=log)


if __name__ == "__main__":
    main()
