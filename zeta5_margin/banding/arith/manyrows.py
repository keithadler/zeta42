"""Lead 1: more rows than poles. deg Delta <= m = #poles, so for fixed m irrationality only needs
log P(zeta s) -> -infinity (any rate). Delta computed by exact interpolation of det(A + xB)."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpq_mat, fmpq_poly, fmpz_mat, arb, ctx
from profile import build, staircase
def delta_poly(s, e, h):
    """Delta(X) = det(A + X B), computed memory-light: clear denominators (L), take exact integer
    determinants det(L A + x L B) at x = 0..d (flint multimodular), interpolate, divide by L^h.
    deg Delta <= d = min(h, #poles), so h > #poles is allowed."""
    a, b = build(s, e, h)
    L = 1
    for c in a + b: L = math.lcm(L, int(c.q))
    ai = [int(c.p)*(L//int(c.q)) for c in a]; bi = [int(c.p)*(L//int(c.q)) for c in b]
    m = sum(1 for v in e.values() if v == -1)
    d = min(m, h)
    xs = list(range(d+1))
    ys = [fmpz_mat(h, h, [ai[i+j] + x*bi[i+j] for i in range(h) for j in range(h)]).det() for x in xs]
    # Newton divided differences over Q, then expand
    coef = [fmpq(int(y)) for y in ys]
    for k in range(1, d+1):
        for i in range(d, k-1, -1):
            coef[i] = (coef[i] - coef[i-1]) / (xs[i] - xs[i-k])
    D = fmpq_poly([coef[d]])
    for i in range(d-1, -1, -1):
        D = D * fmpq_poly([-xs[i], 1]) + coef[i]
    return D / fmpq(L)**h
def delta_poly_interp(s, e, h, a, b):
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
