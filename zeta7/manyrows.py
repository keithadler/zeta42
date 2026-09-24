"""Lead 1: more rows than poles. deg Delta <= m = #poles, so for fixed m irrationality only needs
log P(zeta s) -> -infinity (any rate). Delta computed by exact interpolation of det(A + xB)."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpq_mat, fmpq_poly, arb, ctx
from profile import build, staircase
def delta_poly(s, e, h):
    a, b = build(s, e, h)
    m = sum(1 for v in e.values() if v == -1)
    d = min(m, h)
    xs = list(range(d+1))
    ys = [fmpq_mat(h, h, [a[i+j] + x*b[i+j] for i in range(h) for j in range(h)]).det() for x in xs]
    D = fmpq_poly([0])
    for i, xi in enumerate(xs):        # Lagrange interpolation
        L = fmpq_poly([1])
        for k, xk in enumerate(xs):
            if k != i: L = L * fmpq_poly([-xk, 1]) / (xi - xk)
        D += L * ys[i]
    return D
def score(s, e, h):
    D = delta_poly(s, e, h)
    co = [D[i] for i in range(D.degree()+1)]
    g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    bits = max(max(abs(int(c.p)).bit_length(), int(c.q).bit_length()) for c in co)
    ctx.prec = 2*bits + 3000
    z = arb(s).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= z
    assert v > 0
    return D.degree(), float(v.log()) - (math.log(g) - math.log(l))
if __name__ == "__main__":
    s = int(sys.argv[1]); m = int(sys.argv[2]); hs = [int(x) for x in sys.argv[3].split(',')]
    e = {j: -1 for j in range(1, m+1)}
    for h in hs:
        deg, lP = score(s, e, h)
        print(f"s={s} poles 1..{m} h={h:3d}: deg {deg}  log P(zeta{s}) = {lP:+10.2f}", flush=True)
