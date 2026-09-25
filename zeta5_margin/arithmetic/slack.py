"""Per-prime arithmetic slack of Fauzan's construction at finite K = 40n.
true_p  = v_p(S_K) + v_p^G(Delta_K)      (content of F_K = S_K Delta_K: best possible normaliser)
paper_p = the valuation the proof credits, from its limiting functions:
          inner (K/p >= 3):  L_p ~ -p R(K/p)
          outer (K/p < 3) :  L_p ~ -K (R0(y) - d(y) - 2 lam floor(1/y) + sum_j (2lam - j y)_+),  y = p/K
Slack in log units: (true_p - paper_p) log p.  Reported per range of x = K/p, divided by K^2."""
import sys, math, json
import numpy as np
from math import floor, log
from sympy import primerange
from flint import fmpq, fmpz
from hankel_s import entries, det_poly
al, lam = 3/40, 37/40; Hh = lam + 3*al
def ell(x, z): return floor(x - z) + floor(x + z) + 1
ZG = (np.arange(4000) + 0.5) / 8000          # midpoint grid on (0, 1/2)
def Gamma(x):
    T = floor(2*Hh*x); s = Hh*x - T/2; q = floor(2*x); npl = (2*x - q)/2
    tot = 0.0
    for z in ZG:
        l = ell(x, z); b = 3*ell(al*x, z)
        tot += (T - b)*(T + b - l - 5)
    return tot/8000 + s*(2*T - q - 5) + max(s - npl, 0)
def J(u):
    m = floor(2*u); return m*u - m*(m+1)/4
def Nl(x): return 2*lam*x*floor(x) - 12*lam*x*floor(al*x) - 2*J(lam*x)
def R(x): return -Gamma(x) - Nl(x)
def R0(y):
    pp = lambda u: max(u, 0)
    if 1/3 < y < 1/2: return 8 - 9*y - 8*al - 5*min(al, 1 - 2*y) - 5*pp(1 + al - 3*y)
    if 1/2 < y < 1: return 7*(1 - y) - 6*min(al, 1 - y) - 6*pp(1 + al - 2*y) + pp(1 + 4*al - 2*y)
    return 0.0
def dcor(y):
    pp = lambda u: max(u, 0)
    return pp(1 + 4*al - 3*y - pp(1 + al - 3*y)) if 1/3 < y < 1/2 else 0.0
def outer_minusL(K, p):
    y = p/K
    return K*(R0(y) - dcor(y) - 2*lam*floor(1/y) + sum(max(2*lam - j*y, 0) for j in range(1, 6)))
def vp(x, p):
    x = abs(int(x)); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v
def legendre(m, p):
    v, q = 0, p
    while q <= m: v += m//q; q *= p
    return v
def vS(K, N, h, p):
    return (2*h*legendre(K, p) - 12*h*legendre(N, p)
            - 2*sum(legendre(2*i, p) for i in range(1, h)) + (h - 1)*vp(4, p))
def run(n):
    K, N = 40*n, 3*n; h = K - N
    hh, a, b = entries(5, K, N, 6)
    D = det_poly(hh, a, b)
    co = [D[i] for i in range(h + 1)]
    rows = []
    for p in primerange(2, 2*h + 1):
        vG = min(vp(c.p, p) - vp(c.q, p) for c in co if c != 0)
        true = vS(K, N, h, p) + vG
        x = K/p
        if p < 7: paper = None
        elif x >= 3: paper = -p*R(x)
        else: paper = -outer_minusL(K, p)
        rows.append((p, true, paper))
    return K, rows
if __name__ == "__main__":
    out = {}
    for n in map(int, sys.argv[1:]):
        K, rows = run(n)
        out[K] = rows
        bands = {"p<7": [0, 0], "x>=10": [0, 0], "5<=x<10": [0, 0], "3<=x<5": [0, 0], "x<3 (outer)": [0, 0]}
        for p, t, pap in rows:
            x = K/p
            key = "p<7" if p < 7 else "x>=10" if x >= 10 else "5<=x<10" if x >= 5 else "3<=x<5" if x >= 3 else "x<3 (outer)"
            bands[key][0] += -t*log(p)
            if pap is not None: bands[key][1] += -pap*log(p)
        tot_t = sum(v[0] for v in bands.values()); tot_p = sum(v[1] for k, v in bands.items() if k != "p<7")
        print(f"K={K}: -log cont(F)/K^2 = {tot_t/K**2:.4f}   (paper's A_* = 1.3175)")
        for k, (t, pap) in bands.items():
            print(f"   {k:12s} true {t/K**2:+.4f}   paper-credited {pap/K**2:+.4f}   slack {(pap - t)/K**2 if k != 'p<7' else float('nan'):+.4f}")
        sys.stdout.flush()
    json.dump({str(k): v for k, v in out.items()}, open("slack.json", "w"))
