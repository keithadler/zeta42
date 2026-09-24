"""Lead 2: the |Vandermonde|^4 (beta = 4) analogue of the Hankel determinant.
de Bruijn: int_{(0,oo)^h} prod_{i<k} (t_i - t_k)^4 prod dmu_X(t_i) = h! Pf[(b-a) mu_X(t^{a+b-1} W/D)]_{a,b<2h}.
Positive at X = zeta(s) (positive integrand), polynomial in X with deg <= min(h, #poles).
Pf^2 = det and primitive(Pf)^2 = primitive(Pf^2) (Gauss), so log P_Pf(zeta) = log P_det(zeta) / 2."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpq_mat, fmpq_poly, arb, ctx
from profile import build
def pf_score(s, e, h):
    a, b = build(s, e, 2*h)                  # moments of index 0..4h-2 (need up to 4h-3)
    n = 2*h; m = sum(1 for v in e.values() if v == -1)
    def M(x):
        ent = []
        for i in range(n):
            for j in range(n):
                k = i + j - 1
                ent.append(fmpq(0) if k < 0 else (j - i) * (a[k] + x*b[k]))
        return fmpq_mat(n, n, ent)
    d = 2*min(h, m)
    xs = list(range(d+1)); ys = [M(fmpq(x)).det() for x in xs]
    D = fmpq_poly([0])
    for i, xi in enumerate(xs):
        L = fmpq_poly([1])
        for k, xk in enumerate(xs):
            if k != i: L = L * fmpq_poly([-xk, 1]) / (xi - xk)
        D += L * ys[i]
    co = [D[i] for i in range(D.degree()+1)]
    g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    bits = max(max(abs(int(c.p)).bit_length(), int(c.q).bit_length()) for c in co)
    ctx.prec = 2*bits + 3000
    z = arb(s).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= z
    assert v > 0
    return D.degree()//2, (float(v.log()) - (math.log(g) - math.log(l))) / 2
if __name__ == "__main__":
    s = int(sys.argv[1]); m = int(sys.argv[2])
    e = {j: -1 for j in range(1, m+1)}
    for h in map(int, sys.argv[3].split(',')):
        deg, lP = pf_score(s, e, h)
        print(f"beta=4 s={s} poles 1..{m} h={h:2d}: deg {deg:2d}  log P(zeta{s}) = {lP:+9.2f}  /deg^2 {lP/max(deg,1)**2:+.4f}", flush=True)
