"""Vectorised limiting functions (same bounds as limits.py) and the asymptotic constant A.
A = int_{xc}^{inf} R_gen(x) x^-3 dx  +  int_{y in outer range} g(y) dy  +  int_1^{2 lam} sum_j (2lam - j y)_+ dy
(the last term: p > K, scalar only; v_p(Delta) = 0 there)."""
import numpy as np, sys
from sympy import nextprime, prevprime, primerange
import os
from gout_general import gout
from limits import profile, vS

def exps(ex, n, K):
    e = np.full(K+1, -1.0); e[0] = 0
    for i, v in enumerate(ex, start=1): e[(i-1)*n+1:i*n+1] = v
    return e

def inner_gamma_fast(ex, n, K, h, p):
    e = exps(ex, n, K); j = np.arange(1, K+1); ej = e[1:]
    r = j % p; a = np.minimum(r, p - r); m = (p-1)//2
    pole = ej == -1
    zc = r == 0
    mK = int(zc.sum()); nz = float((ej[zc & ~pole] + 1).sum())
    ell = np.bincount(a[~zc], minlength=m+1)[1:].astype(float)
    b = np.bincount(a[~zc & ~pole], weights=(ej[~zc & ~pole] + 1)/2, minlength=m+1)[1:]
    base = b - (ell + 4)/2
    best = -np.inf
    for L0 in range(0, int(4*K/p) + 14):
        R = h - L0
        lo, hi = -1e6, 1e6
        for _ in range(100):
            th = (lo+hi)/2
            if np.maximum(0, np.floor(th - b)).sum() > R: hi = th
            else: lo = th
        Lc = np.maximum(0, np.floor(lo - b)); rem = int(R - Lc.sum())
        order = np.lexsort((np.arange(m), -ell, Lc + b)); Lc[order[:rem]] += 1
        Phi = Lc + base
        srt = np.sort(Phi); zsrc = 2*L0 + nz - mK + 0.5
        cap = np.minimum(np.where(np.arange(m) == np.argmin(Phi), srt[1], srt[0]), zsrc)
        k = np.clip(np.floor(cap - base) + 1, 0, Lc)
        tot = (k*base + k*(k-1)/2 + (Lc - k)*cap).sum()
        i = np.arange(L0); tot += np.minimum(2*i + nz - mK + 0.5, srt[0]).sum()
        best = max(best, 2*tot)
    return best

def R_gen(ex, alpha, x, p):
    n = max(1, round(x*p/40)); e, K, h = profile(ex, alpha, n)
    return -(vS(e, K, h, p) + inner_gamma_fast(ex, n, K, h, p)) / p, K/p

def g_gen(ex, alpha, y, K0):
    n = round(K0/40); e, K, h = profile(ex, alpha, n)
    p = int(nextprime(int(y*K)))
    if p > K: p = int(prevprime(K + 1))
    return -float(vS(e, K, h, p) + gout(e, K, h, p)) / K, p/K

def A_const(ex, alpha, lam, xc, outer_fn=None, p=1009, K0=4000, xmax=80, dx=0.01, log=print):
    cache = f"inner_{'_'.join(map(str, ex))}_a{alpha}_xc{xc}_p{p}.npz"
    xx, vals = np.array([]), np.array([])
    if os.path.exists(cache):
        z = np.load(cache); xx, vals = z["xx"], z["vals"]
    start = xc if len(xx) == 0 else xx.max() + dx
    if start < xmax:
        xs = np.arange(start, xmax, dx)
        vals = []; xx = []
        for x in xs:
            v, x_ = R_gen(ex, alpha, x, p); vals.append(v); xx.append(x_)
        xx = np.concatenate([z["xx"], xx]) if os.path.exists(cache) else np.array(xx)
        vals = np.concatenate([z["vals"], vals]) if os.path.exists(cache) else np.array(vals)
        np.savez(cache, xx=xx, vals=vals); o = np.argsort(xx); xx, vals = xx[o], vals[o]
    inner = np.trapezoid(vals/xx**3, xx)
    per = (xx >= xx.max() - 40)          # one full period of the band/floor pattern
    tailF = np.trapezoid(vals[per]/xx[per], xx[per]) / (xx[per].max() - xx[per].min())
    inner += tailF/xx[-1]
    ys = np.arange(max(0.5, 1/xc), 1.0, 0.0025)
    if outer_fn is None:
        gv = np.array([g_gen(ex, alpha, y, K0) for y in ys]); yy = gv[:, 1]; gg = gv[:, 0]
    else:
        yy = ys; gg = np.array([outer_fn(y) for y in ys])
    outer = np.trapezoid(gg, yy)
    big_y = np.linspace(1, 2*lam, 4001)
    big = np.trapezoid([sum(max(2*lam - j*y, 0) for j in range(1, 6)) for y in big_y], big_y)
    log(f"   inner [{xc}, inf) {inner:+.5f} (tail mean F {tailF:+.3f})   outer {outer:+.5f}   p>K {big:+.5f}")
    return inner + outer + big

if __name__ == "__main__":
    sys.path.insert(0, "../../arithmetic")
    from slack import outer_minusL
    import time; t = time.time()
    # control: paper, inner from 5/2, outer (0.4, 0.5] from the paper's formula, (0.5, 1] generalised
    lam = 37/40
    A_ctrl = A_const((5, 5, 5), 1.0, lam, 2.5)
    # add paper outer on (0.4, 0.5]
    ys = np.linspace(0.4, 0.5, 2001); add = np.trapezoid([float(outer_minusL(1, y)) for y in ys], ys)
    print(f"paper control: A = {A_ctrl + add:.5f}   exact A_* - 187/9600 = {9928298118277006344769/7535670527041937280000 - 187/9600:.5f}   [{time.time()-t:.0f}s]", flush=True)
