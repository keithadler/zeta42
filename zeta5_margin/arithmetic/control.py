import sympy as sp
from fractions import Fraction as F
import exactA as E
al, lam, H = E.al, E.lam, E.H
def all_breaks(a, b):
    c = set([a, b])
    for k in range(0, 100):
        for cc in (2*H, F(2), 2*lam, F(1), al):          # floor jumps
            c.add(F(k)/cc)
        c.add(F(k) + F(1, 2)); c.add((F(k) + F(1, 2))/al)  # frac = 1/2
    for A in range(0, 3):
        for B in range(0, 25):
            c.add(F(B - A)/(1 - al))                     # frac(al x) = frac(x)
            c.add(F(1 + A + B)/(1 + al))                 # frac(al x) = 1 - frac(x)
    for T in range(0, 60):
        for q in range(0, 60):
            c.add(F(T - q)/(2*(H - 1)))                  # s = n_+
    return sorted(v for v in c if a <= v <= b)
X = E.X
tot = 0
xs = all_breaks(F(3), F(20))
for u, v in zip(xs, xs[1:]):
    pol = E.piece_poly(E.R, u, v, 2)
    tot += sp.integrate(pol/X**3, (X, sp.Rational(u.numerator, u.denominator), sp.Rational(v.numerator, v.denominator)))
paper = sp.Rational(322437603634266857629, 7535670527041937280000)
print(f"inner [3,20]: {len(xs)-1} pieces, matches paper exactly: {sp.simplify(tot - paper) == 0}")
tot2 = 0
xs2 = all_breaks(F(5, 2), F(3))
for u, v in zip(xs2, xs2[1:]):
    pol = E.piece_poly(E.R, u, v, 2)
    tot2 += sp.integrate(pol/X**3, (X, sp.Rational(u.numerator, u.denominator), sp.Rational(v.numerator, v.denominator)))
print(f"inner [5/2,3] with full breakpoint list: {sp.simplify(tot2)}  (earlier 10433/48000)")
