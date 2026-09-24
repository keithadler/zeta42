"""Idea B: per-prime decomposition of the arithmetic side. -log content = sum_p -v_p(content) log p.
Compare zeta(7) with zeta(5) on the same construction, grouped by x = K/p."""
import math
from sympy import primerange
from flint import fmpz, fmpq_mat
from profile import build
def delta_coeffs(s, e, h):
    a, b = build(s, e, h)
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    D = (-(B.inv()*A)).charpoly()*B.det()
    return [D[i] for i in range(D.degree()+1)]
def vp(n, p):
    n = abs(int(n))
    if n == 0: return 10**9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k
def profile_of(s, e, h, pmax):
    co = [c for c in delta_coeffs(s, e, h) if c != 0]
    out = {}
    for p in primerange(2, pmax):
        v = min(vp(c.p, p) - vp(c.q, p) for c in co)
        out[p] = -v*math.log(p)          # contribution to -log content
    return out
K = 40
e = {j: -1 for j in range(1, K+1)}
pr5 = profile_of(5, e, K, 6*K); pr7 = profile_of(7, e, K, 6*K)
bands = [(0, 1/8), (1/8, 1/4), (1/4, 1/2), (1/2, 1), (1, 2), (2, 4), (4, 10), (10, 100)]
print(f"K={K}, plain construction; contributions to -log content / K^2 by band of p/K")
print(f"{'p/K band':>14s} {'zeta5':>9s} {'zeta7':>9s} {'7 - 5':>9s}")
for lo, hi in bands:
    a5 = sum(v for p, v in pr5.items() if lo < p/K <= hi)/K**2
    a7 = sum(v for p, v in pr7.items() if lo < p/K <= hi)/K**2
    print(f"{lo:6.3f}-{hi:<6.3f} {a5:+9.4f} {a7:+9.4f} {a7-a5:+9.4f}")
print(f"{'total':>14s} {sum(pr5.values())/K**2:+9.4f} {sum(pr7.values())/K**2:+9.4f} {(sum(pr7.values())-sum(pr5.values()))/K**2:+9.4f}")
