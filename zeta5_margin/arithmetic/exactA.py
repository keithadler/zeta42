"""Exact dA for the split x_c = 5/2: sympy integration of exact piecewise polynomials.
inner:  int_{5/2}^{3} R(x) x^-3 dx      outer removed: int_{1/3}^{2/5} (R0 - d - 2 lam floor(1/y) + sum (2lam - jy)_+) dy
R and the outer integrand are evaluated exactly with Fractions; on each piece between
breakpoints they are polynomials (degree <= 2 in x; <= 1 in y), recovered by interpolation
and cross-checked at an extra point."""
from fractions import Fraction as F
from math import floor
import sympy as sp
al, lam = F(3, 40), F(37, 40); H = lam + 3*al
def fl(v): return floor(v)
def ell(x, z): return fl(x - z) + fl(x + z) + 1
def Gamma(x):
    T = fl(2*H*x); s = H*x - F(T, 2); q = fl(2*x); npl = (2*x - q)/2
    fr = lambda v: v - fl(v)
    br = sorted(set([F(0), F(1, 2)] + [c for c in (fr(al*x), 1 - fr(al*x), fr(x), 1 - fr(x)) if 0 < c < F(1, 2)]))
    tot = F(0)
    for u, v in zip(br, br[1:]):
        z = (u + v)/2; l = ell(x, z); b = 3*ell(al*x, z)
        tot += (v - u)*(T - b)*(T + b - l - 5)
    return tot + s*(2*T - q - 5) + max(s - npl, F(0))
def J(u):
    m = fl(2*u); return m*u - F(m*(m+1), 4)
def R(x): return -Gamma(x) - (2*lam*x*fl(x) - 12*lam*x*fl(al*x) - 2*J(lam*x))
pp = lambda u: max(u, F(0))
def outer(y):
    R0 = 8 - 9*y - 8*al - 5*min(al, 1 - 2*y) - 5*pp(1 + al - 3*y)
    d = pp(1 + 4*al - 3*y - pp(1 + al - 3*y))
    return R0 - d - 2*lam*fl(1/y) + sum(pp(2*lam - j*y) for j in range(1, 6))
X = sp.Symbol('x', positive=True)
def piece_poly(f, a, b, deg):
    pts = [a + (b - a)*F(k + 1, deg + 3) for k in range(deg + 2)]
    poly = sp.interpolate([(sp.Rational(p.numerator, p.denominator), sp.Rational(f(p).numerator, f(p).denominator)) for p in pts[:-1]], X)
    chk = pts[-1]; assert poly.subs(X, sp.Rational(chk.numerator, chk.denominator)) == sp.Rational(f(chk).numerator, f(chk).denominator), "not polynomial on piece"
    return sp.expand(poly)
def breakpoints(a, b, fns):
    """candidate breakpoints: fine rational grid refined where the piecewise form changes."""
    cand = set([a, b])
    # jumps of floor(c x) and z-break crossings, all rational
    for c in (2*H, F(2), 2*lam, F(1), al):
        for k in range(0, 40):
            v = F(k)/c
            if a < v < b: cand.add(v)
    for k1 in range(0, 4):
        for k2 in range(-1, 5):
            for s1 in (1, -1):
                for s2 in (1, -1):
                    den = s1*al - s2
                    for h in (F(0), F(1, 2)):
                        for v in ((k1 - k2 + (s2 - s1)*0)/den, (h + k2)/1, (h + k1)/al):
                            if a < v < b: cand.add(v)
    return sorted(cand)
xs = breakpoints(F(5, 2), F(3), None)
inner = 0
for u, v in zip(xs, xs[1:]):
    pol = piece_poly(R, u, v, 2)
    inner += sp.integrate(pol/X**3, (X, sp.Rational(u.numerator, u.denominator), sp.Rational(v.numerator, v.denominator)))
ys = sorted(set([F(1, 3), F(2, 5), (1 + al)/3, (1 + 4*al)/3, (1 - al)/2] + [F(1, k) for k in (2, 3)] + [2*lam/j for j in range(1, 6)]))
ys = [y for y in ys if F(1, 3) <= y <= F(2, 5)]
outer_int = 0
for u, v in zip(ys, ys[1:]):
    pol = piece_poly(outer, u, v, 1)
    outer_int += sp.integrate(pol, (X, sp.Rational(u.numerator, u.denominator), sp.Rational(v.numerator, v.denominator)))
dA = sp.nsimplify(sp.simplify(inner - outer_int))
print("inner int_{5/2}^3 R x^-3 dx =", sp.simplify(inner))
print("outer int_{1/3}^{2/5}        =", sp.simplify(outer_int))
print("dA =", sp.simplify(inner - outer_int), "=", sp.N(inner - outer_int, 30))
