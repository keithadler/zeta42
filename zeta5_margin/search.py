"""Build a nested arcsine comparison measure close to the equilibrium, then rationalise it.

Intervals: the paper's 16 nested supports, interpolated (monotone, in log scale) to n
components.  Weights: for fixed intervals the certified quantity
    lam * max_t (2 U^rho(t) - V(t)) - I(rho)
is convex in c (-I is a positive combination of squared partial sums), so it is solved
exactly by a conic program on a fine t-grid.  The result is rounded to denominator 10^12,
the mass is restored exactly, and M0 is set a little above the grid maximum; `certify.py`
then decides rigorously.

    python search.py 256 measure256.json     (~10 min)
"""
import json
import sys
from fractions import Fraction as Fr

import cvxpy as cp
import numpy as np
from scipy.interpolate import PchipInterpolator

from field import V, U_arcsine, LAM, CSTAR, real_constant

DEN = 10**12
PAPER = json.load(open("paper16.json"))["rows"]
PA = np.array([r[0] for r in PAPER]) / DEN
PB = np.array([r[1] for r in PAPER]) / DEN


def intervals(n, warp=1.0):
    x = np.linspace(0, 1, 16)
    s = np.linspace(0, 1, n) ** warp
    return (np.exp(PchipInterpolator(x, np.log(PA))(s)),
            np.exp(PchipInterpolator(x, np.log(PB))(s)))


def grid(a, b, base=8000):
    pts = [np.geomspace(1e-8, 2, base)]
    for x in np.concatenate([a, b]):
        d = np.geomspace(1e-7, 0.05, 60)
        pts += [x * (1 + d), x * (1 - d)]
    t = np.unique(np.concatenate(pts))
    return t[(t > 0) & (t <= 2)]


def weights(a, b):
    t = grid(a, b)
    Um = np.stack([U_arcsine(t, a[j], b[j]) for j in range(len(a))], 1)
    L = np.log((b - a) / 4)
    c, z = cp.Variable(len(a), nonneg=True), cp.Variable()
    S = cp.cumsum(c)
    negI = cp.sum(cp.multiply(L[1:] - L[:-1], cp.square(S[:-1]))) - LAM**2 * L[-1]
    cp.Problem(cp.Minimize(LAM * z + negI),
               [2 * Um @ c - V(t) <= z, cp.sum(c) == LAM]).solve(solver="CLARABEL")
    cv = np.maximum(c.value, 0)
    return cv * LAM / cv.sum()


def rationalise(a, b, c, slack):
    A = [Fr(round(x * DEN), DEN) for x in a]
    B = [Fr(round(x * DEN), DEN) for x in b]
    C = [Fr(max(1, round(x * DEN)), DEN) for x in c]
    C[int(np.argmax(c))] += Fr(37, 40) - sum(C)
    af, bf, cf = (np.array([float(x) for x in X]) for X in (A, B, C))
    U, M = real_constant(af, bf, cf, grid(af, bf, 40000))
    M0 = Fr(int(np.ceil((M + slack) * 10**6)), 10**6)
    rows = [[int(x * DEN) for x in r] for r in zip(A, B, C)]
    return {"denominator": DEN, "rows": rows, "M0": [M0.numerator, M0.denominator]}, U


if __name__ == "__main__":
    n, out = int(sys.argv[1]), sys.argv[2]
    a, b = intervals(n, 1.25 if n <= 64 else 1.0)
    c = weights(a, b)
    d, U = rationalise(a, b, c, 1e-5 if n > 128 else 2e-5)
    json.dump(d, open(out, "w"), indent=0)
    print(f"{out}: {n} components, float estimate U = {U:.6f}; now run certify.py {out}")
