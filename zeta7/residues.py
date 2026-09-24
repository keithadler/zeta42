"""Designed residues: measure w(y) R(y^2) with R(t) = sum_j c_j/(t + j^2), any c_j with R >= 0 on [0, oo).
mu(t^e/(t+j^2)) = sum_{i<e} (-j^2)^{e-1-i} mu(t^i) + (-j^2)^e M_j,  M_j = j^{s-1}(X - H_j) - 1/(s-1) + 1/(2j).
The X-part of each entry is (1/720) int_0^1 (-log z)^6 P_e(z)/(1-z) dz with P_e(z) = sum_j c_j (-j^2)^e j^6 z^j."""
import math, sys
from math import gcd, lcm, comb
from flint import fmpq, fmpq_mat, fmpq_poly, arb, ctx
import mpmath
def build_res(s, c, h):
    js = sorted(c); E = 2*h - 1; sf = math.factorial(s-1)
    mus = [fmpq(-1)**k * fmpq.bernoulli(2*k+2) * fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf) for k in range(E)]
    Hs = {}; acc = fmpq(0)
    for j in range(1, max(js)+1):
        acc += fmpq(1, j**s); Hs[j] = acc
    a, b = [], []
    for e in range(E):
        ae = fmpq(0); be = fmpq(0)
        for j in js:
            cj = fmpq(c[j]); q = fmpq(-j*j)
            poly = sum((q**(e-1-i) * mus[i] for i in range(e)), fmpq(0))
            ae += cj * (poly + q**e * (-fmpq(j)**(s-1)*Hs[j] - fmpq(1, s-1) + fmpq(1, 2*j)))
            be += cj * q**e * fmpq(j)**(s-1)
        a.append(ae); b.append(be)
    return a, b
def R_nonneg(c):
    """R(t) = sum c_j/(t+j^2) >= 0 on [0, oo), checked exactly: R = N/D with D > 0, so check N >= 0.
    N is a polynomial; count its real roots on [0, oo) (Sturm via flint) and check N(0), leading sign."""
    from flint import fmpz_poly
    D = fmpz_poly([1])
    for j in c: D *= fmpz_poly([j*j, 1])
    N = fmpz_poly([0])
    for j, cj in c.items():
        N += cj * (D // fmpz_poly([j*j, 1]))
    if N == 0: return False
    if N(0) < 0: return False
    lead = N[N.degree()]
    if lead < 0: return False
    # no sign change on (0, oo): no root of odd multiplicity there (exact: square-free factors + Sturm)
    import sympy
    x = sympy.symbols('x')
    P = sympy.Poly([int(N[i]) for i in range(N.degree(), -1, -1)], x)
    for fac, mult in P.sqf_list()[1]:
        if mult % 2 == 1 and fac.count_roots(0, None) - (1 if fac.eval(0) == 0 else 0) > 0:
            return False
    return True

def score_res(s, c, h):
    a, b = build_res(s, c, h)
    m = len(c); d = min(m, h)
    xs = list(range(d+1))
    ys = [fmpq_mat(h, h, [a[i+j] + x*b[i+j] for i in range(h) for j in range(h)]).det() for x in xs]
    D = fmpq_poly([0])
    for i, xi in enumerate(xs):
        L = fmpq_poly([1])
        for k, xk in enumerate(xs):
            if k != i: L = L * fmpq_poly([-xk, 1]) / (xi - xk)
        D += L * ys[i]
    co = [D[i] for i in range(D.degree()+1)]
    g, l = 0, 1
    for cc in co: g = gcd(g, int(cc.p)); l = lcm(l, int(cc.q))
    bits = max(max(abs(int(cc.p)).bit_length(), int(cc.q).bit_length()) for cc in co)
    ctx.prec = 2*bits + 3000
    z = arb(s).zeta(); v = arb(0); zp = arb(1)
    for cc in co: v += arb(int(cc.p))/arb(int(cc.q))*zp; zp *= z
    assert v > 0
    return D.degree(), float(v.log()) - (math.log(g) - math.log(l))
if __name__ == "__main__":
    # validation: c_j = (-1)^j C(2K,K-j) j^2 is the product family 1/D_K (poles 1..K, W = 1)
    K = 20
    c = {j: (-1)**(j-1) * comb(2*K, K-j) * j*j for j in range(1, K+1)}
    print("nonneg:", R_nonneg(c), " score:", score_res(7, c, K)[1]/K**2, "(expect +0.8200 = plain poles 1..20)")
