"""Limiting per-prime credits of the generalised bounds, sampled with a large prime p and K ~ x p.
inner (p <= K/2):  R_gen(x) = -(v_p(S) + gamma_in) / p       (paper: R(x))
outer (K/2 < p <= K): g(y) = -(v_p(S) + gamma_out) / K        (paper: R0 - d - 2 lam floor(1/y) + sum(...))
Float arithmetic: this estimates constants; exact integrals come after."""
import numpy as np
from sympy import nextprime
from gout_general import gout

def profile(ex, alpha, n):
    K = 40*n; e = {}
    for i, v in enumerate(ex, start=1):
        for j in range((i-1)*n+1, i*n+1):
            if v: e[j] = v
    for j in range(len(ex)*n+1, K+1): e[j] = -1
    m = sum(1 for v in e.values() if v == -1)
    return e, K, round(alpha*m)

def vp_int(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def vS(e, K, h, p):
    v = 0
    for j in range(p, K+1, p):
        ej = e.get(j, 0)
        if ej: v += -2*h*ej*vp_int(j, p)
    s, q = 0, p
    while q <= 2*h:
        s += sum((2*i)//q for i in range(1, h)); q *= p
    return v - 2*s          # p >= 7: no 4-term

def inner_gamma(e, K, h, p):
    m = (p-1)//2
    ell = np.zeros(m+1); b = np.zeros(m+1); nz = 0.0; mK = 0
    for j in range(1, K+1):
        r = j % p
        pole = e.get(j, 0) == -1
        if r == 0:
            mK += 1
            if not pole: nz += e.get(j, 0) + 1
            continue
        a = min(r, p-r); ell[a] += 1
        if not pole: b[a] += (e.get(j, 0) + 1)/2
    best = -np.inf
    for L0 in range(0, int(4*K/p) + 14):
        # paper's rule generalised: equalise Z_c = L_c + b_c (ties: larger ell first), L_c >= 0
        tot_rows = h - L0
        lo, hi = -1e6, 1e6
        for _ in range(200):
            th = (lo+hi)/2
            Lc = np.maximum(0, np.floor(th - b[1:]))
            if Lc.sum() > tot_rows: hi = th
            else: lo = th
        Lc = np.maximum(0, np.floor(lo - b[1:])); rem = int(tot_rows - Lc.sum())
        order = np.lexsort((np.arange(m), -ell[1:], Lc + b[1:]))
        Lc[order[:rem]] += 1
        L = np.concatenate([[L0], Lc])
        base = b[1:] - (ell[1:] + 4)/2
        Phi = L[1:] + base
        srt = np.sort(Phi); lo1 = srt[0]; lo2 = srt[1]
        zsrc = 2*L0 + nz - mK + 0.5
        tot = 0.0
        amin = np.argmin(Phi)
        for c in range(m):
            cap = min(lo2 if c == amin else lo1, zsrc)
            n_ = int(L[c+1])
            if n_ == 0: continue
            i = np.arange(n_)
            tot += np.minimum(i + base[c], cap).sum()
        i = np.arange(L0)
        tot += np.minimum(2*i + nz - mK + 0.5, lo1).sum()
        best = max(best, 2*tot)
    return best

def R_gen(ex, alpha, x, p0=2003):
    p = int(nextprime(p0)); n = max(1, round(x*p/40))
    e, K, h = profile(ex, alpha, n)
    return -(vS(e, K, h, p) + inner_gamma(e, K, h, p)) / p, K/p

def g_gen(ex, alpha, y, K0=8000):
    n = round(K0/40); e, K, h = profile(ex, alpha, n)
    p = int(nextprime(int(y*K)))
    return -(vS(e, K, h, p) + gout(e, K, h, p)) / K, p/K

if __name__ == "__main__":
    import sys; sys.path.insert(0, "../../arithmetic")
    from slack import R, outer_minusL
    for x in (2.6, 3.0, 3.7, 5.3, 8.1, 12.6):
        v, xx = R_gen((5, 5, 5), 1.0, x, p0=1009)
        print(f"inner x={xx:.4f}: generalised {v:+.4f}   paper R(x) {R(xx):+.4f}", flush=True)
    for y in (0.55, 0.7, 0.9):
        v, yy = g_gen((5, 5, 5), 1.0, y, K0=4000); v = float(v)
        print(f"outer y={yy:.4f}: generalised {v:+.4f}   paper {float(outer_minusL(1, yy)):+.4f}", flush=True)
