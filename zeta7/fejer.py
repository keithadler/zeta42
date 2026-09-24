"""Complete positive family: numerator N(t) >= 0 on [0,oo) iff N(y^2) = |Pi(iy)|^2 (Fejer-Riesz / Markov-Lukacs).
Residues c_j ~ Pi(j) Pi(-j) / D'(-j^2). Integer-root Pi = the product family; here: non-product Pi."""
import math
from math import comb
import sympy
from flint import fmpz_poly
from hankel_s import det_poly
import hankel_s
from sweep import show
from math import gcd, lcm
from flint import arb, ctx
z, t = sympy.symbols('z t')
def N_from_Pi(Pi):
    """N(t) with N(y^2) = Pi(iy) Pi(-iy)."""
    y = sympy.symbols('y', real=True)
    e = sympy.expand(Pi.subs(z, sympy.I*y) * Pi.subs(z, -sympy.I*y))
    e = sympy.expand(e.subs(y, sympy.sqrt(t)))
    P = sympy.Poly(e, t)
    return fmpz_poly([int(P.coeff_monomial(t**k)) for k in range(P.degree()+1)])
def score_W(s, K, W):
    h, a, b = hankel_s.entries(s, K, 0, 1, W)
    from flint import fmpq, fmpq_mat, fmpq_poly
    npol = sum(1 for j in range(1, K+1) if W(-j*j) != 0); d = min(h, npol)
    xs = list(range(d+1))
    ys = [fmpq_mat(h, h, [a[i+j] + x*b[i+j] for i in range(h) for j in range(h)]).det() for x in xs]
    D = fmpq_poly([0])
    for i, xi in enumerate(xs):
        L = fmpq_poly([1])
        for k, xk in enumerate(xs):
            if k != i: L = L * fmpq_poly([-xk, 1]) / (xi - xk)
        D += L * ys[i]
    co = [D[i] for i in range(D.degree()+1)]; g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    ctx.prec = 64*h + 6*K*K + 3000
    zz = arb(s).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= zz
    assert v > 0
    return (float(v.log()) - (math.log(g) - math.log(l))) / K**2
if __name__ == "__main__":
  K = 20
  cands = {}
  for L in (3, 5, 7):
      cands[f"product prod_(i<=L)(z+i), L={L} (control)"] = sympy.prod([z + i for i in range(1, L+1)])
      cands[f"Legendre-type sum_k C(L,k)C(L+k,k)C(z,k), L={L}"] = sum(comb(L, k)*comb(L+k, k)*sympy.binomial(z, k) for k in range(L+1))
      cands[f"Apery-type sum_k C(L,k)^2 C(L+k,k) C(z+k,k), L={L}"] = sum(comb(L, k)**2*comb(L+k, k)*sympy.expand_func(sympy.binomial(z+k, k)) for k in range(L+1))
      cands[f"central sum_k C(L,k)^2 C(z,k), L={L}"] = sum(comb(L, k)**2*sympy.binomial(z, k) for k in range(L+1))
  for name, Pi in cands.items():
      Pi = sympy.expand(sympy.expand_func(Pi))
      den = sympy.Poly(Pi, z).clear_denoms()[0]
      Pi = sympy.expand(Pi * den)
      W = N_from_Pi(Pi)
      print(f"{name:52s} s=5 {score_W(5, K, W):+.4f}   s=7 {score_W(7, K, W):+.4f}", flush=True)
