"""Idea A: does the primitive Hankel polynomial P(X) factor over Z? Any factor Q | P is nonzero at zeta(s)
and may be much smaller there, and the criterion applies to it."""
import math, sys
from math import gcd, lcm
from flint import fmpz_poly, arb, ctx
from profile import build, staircase
from flint import fmpq_mat
def prim(s, e, h):
    a, b = build(s, e, h)
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    D = (-(B.inv()*A)).charpoly()*B.det()
    co = [D[i] for i in range(D.degree()+1)]; g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    return fmpz_poly([int(c.p)*(l//int(c.q))//g for c in co])
def logval(P, s):
    bits = max(abs(int(P[i])).bit_length() for i in range(P.degree()+1))
    ctx.prec = 2*bits + 3000
    return float(abs(P(arb(s).zeta())).log())
for s, K, N, r in [(5, 20, 0, 0), (5, 40, 3, 5), (7, 20, 0, 0), (7, 40, 5, 4), (3, 20, 0, 0)]:
    e = staircase(K, [(N, r)] if N else [], (N, K)); h = K - N
    P = prim(s, e, h)
    c, facs = P.factor()
    total = logval(P, s)
    desc = ", ".join(f"deg {f.degree()}^{m}" for f, m in facs)
    print(f"s={s} K={K} N={N}: deg P={P.degree()}, log P(zeta)={total:+.1f}; factors: {desc}", flush=True)
    if len(facs) > 1 or facs[0][1] > 1:
        for f, m in facs:
            print(f"    factor deg {f.degree()}: log|f(zeta{s})| = {logval(f, s):+.2f}", flush=True)
