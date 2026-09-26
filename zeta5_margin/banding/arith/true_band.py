"""True per-prime valuations of F = S * Delta for a banded construction.
S = prod_j j^(-2 h eps_j) * 4^(h-1) / prod_{i<h} ((2i)!)^2   (the paper's S_K for its own profile)."""
import sys, json, math
from math import log
from sympy import primerange
from manyrows import delta_poly
def build(ex, m0, n, alpha):
    K = 40*n; e = {}
    for i, v in enumerate(ex, start=1):
        for j in range((i-1)*n+1, i*n+1):
            if v: e[j] = v
    for j in range(len(ex)*n+1, K+1): e[j] = -1
    if m0: e[0] = m0
    m = sum(1 for v in e.values() if v == -1)
    return e, K, round(alpha*m)
def vp(x, p):
    x = abs(int(x)); v = 0
    if x == 0: return 10**9
    while x % p == 0: x //= p; v += 1
    return v
def leg(m, p):
    v, q = 0, p
    while q <= m: v += m//q; q *= p
    return v
def vS(e, h, p):
    v = 0
    for j, ej in e.items():
        if j > 0: v += -2*h*ej*vp(j, p)
    v += (h-1)*vp(4, p) - 2*sum(leg(2*i, p) for i in range(1, h))
    return v
def logS(e, h):
    s = sum(-2*h*ej*math.log(j) for j, ej in e.items() if j > 0)
    return s + (h-1)*math.log(4) - 2*sum(math.lgamma(2*i+1) for i in range(1, h))
def run(ex, m0, n, alpha):
    e, K, h = build(ex, m0, n, alpha)
    D = delta_poly(5, e, h)
    co = [D[i] for i in range(D.degree()+1)]
    rows = []
    for p in primerange(2, 2*h + 2):
        vG = min(vp(c.p, p) - vp(c.q, p) for c in co if c != 0)
        rows.append((p, vS(e, h, p) + vG, vG))
    return K, h, rows, logS(e, h)
if __name__ == "__main__":
    ex = tuple(int(x) for x in sys.argv[1].split(',')); m0 = int(sys.argv[2]); alpha = float(sys.argv[3])
    for n in map(int, sys.argv[4:]):
        K, h, rows, lS = run(ex, m0, n, alpha)
        tag = "_".join(map(str, ex)) + f"_t{m0}_a{alpha}"
        json.dump({"K": K, "h": h, "rows": rows, "logS": lS}, open(f"true_{tag}_K{K}.json", "w"))
        minus_log_cont = -sum(t*log(p) for p, t, _ in rows)
        print(f"{tag} K={K} h={h}: -log cont(F)/K^2 = {minus_log_cont/K**2:.4f}", flush=True)
