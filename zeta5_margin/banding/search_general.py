"""Near-equilibrium nested arcsine measure for a band profile, then rationalise for certify_general.py.
Intervals: innermost around the minimiser q of V, outer ones spreading geometrically (left, towards
the equilibrium support's lower end) and by a power law (right, towards its upper end); weights by
the convex program as in ../search.py; shape parameters tuned on a coarse grid at 64 components."""
import json, sys
from fractions import Fraction as Fr
import numpy as np, cvxpy as cp
from scipy.optimize import minimize_scalar
from profiles import PROFILES
from equilibrium_general import make_V, solve

def U_arc(t, a, b):
    m = (a + b)/2; inside = (t >= a) & (t <= b)
    return np.where(inside, np.log((b - a)/4), np.log((np.abs(t - m) + np.sqrt(np.abs((t - a)*(t - b))))/2))

def grid(a, b, base=6000):
    pts = [np.geomspace(1e-9, 2, base)]
    for x in np.concatenate([a, b]):
        d = np.geomspace(1e-7, 0.05, 50); pts += [x*(1 + d), x*(1 - d)]
    t = np.unique(np.concatenate(pts)); return t[(t > 0) & (t <= 2)]

def weights(V, lam, a, b):
    t = grid(a, b); Um = np.stack([U_arc(t, a[j], b[j]) for j in range(len(a))], 1)
    L = np.log((b - a)/4); c, z = cp.Variable(len(a), nonneg=True), cp.Variable(); S = cp.cumsum(c)
    negI = cp.sum(cp.multiply(L[1:] - L[:-1], cp.square(S[:-1]))) - lam**2*L[-1]
    pr = cp.Problem(cp.Minimize(lam*z + negI), [2*Um@c - V(t) <= z, cp.sum(c) == lam]); pr.solve(solver="CLARABEL")
    cv = np.maximum(c.value, 0); return cv*lam/cv.sum(), pr.value

def intervals(n, q, lo, hi, w1, w2, half=0.0026):
    s = np.linspace(0, 1, n)
    a1, b1 = q - half, q + half
    a = np.exp(np.log(a1) + (np.log(lo) - np.log(a1))*s**w1)
    b = b1 + (hi - b1)*s**w2
    return a, b

def best_measure(prof, lam, n):
    V = make_V(prof); lamf = float(lam)
    q = minimize_scalar(V, bounds=(1e-6, 0.5), method="bounded").x
    R, (lo, hi) = solve(V, lamf, n=800)
    best = None
    for w1 in (0.8, 1.0, 1.3):
        for w2 in (1.0, 1.4, 1.8, 2.3):
            for f in (1.0, 0.6, 1.05):
                a, b = intervals(64, q, lo*f if f != 1.05 else lo, hi*(1.02 if f == 1.05 else 1.0), w1, w2)
                if not (np.all(np.diff(a) < 0) and np.all(np.diff(b) > 0) and a[0] > 0): continue
                try: c, val = weights(V, lamf, a, b)
                except Exception: continue
                if best is None or val < best[0]: best = (val, w1, w2, f)
    _, w1, w2, f = best
    a, b = intervals(n, q, lo*f if f != 1.05 else lo, hi*(1.02 if f == 1.05 else 1.0), w1, w2)
    c, val = weights(V, lamf, a, b)
    return a, b, c, val, R, V

def rationalise(prof, lam, a, b, c, V, slack=1e-5, DEN=10**12):
    A_ = [Fr(round(x*DEN), DEN) for x in a]; B_ = [Fr(round(x*DEN), DEN) for x in b]
    C_ = [Fr(max(1, round(x*DEN)), DEN) for x in c]; C_[int(np.argmax(c))] += lam - sum(C_)
    af, bf, cf = (np.array([float(x) for x in X]) for X in (A_, B_, C_))
    t = grid(af, bf, 40000); M = np.max(2*sum(cj*U_arc(t, aj, bj) for aj, bj, cj in zip(af, bf, cf)) - V(t))
    M0 = Fr(int(np.ceil((M + slack)*10**6)), 10**6)
    return {"denominator": DEN, "rows": [[int(x*DEN) for x in r] for r in zip(A_, B_, C_)],
            "M0": [M0.numerator, M0.denominator], "lam": [lam.numerator, lam.denominator],
            "profile": [[[u0.numerator, u0.denominator], [u1.numerator, u1.denominator], [e.numerator, e.denominator]] for u0, u1, e in prof]}

if __name__ == "__main__":
    name, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    prof, lam = PROFILES[name]
    a, b, c, val, R, V = best_measure(prof, lam, n)
    d = rationalise(prof, lam, a, b, c, V)
    json.dump(d, open(out, "w"), indent=0)
    print(f"{name}: equilibrium R = {R:.5f};  {n}-component measure estimate {val:.5f}  -> {out}", flush=True)
