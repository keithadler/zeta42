"""Equilibrium value R = sup_sigma [ I(sigma) - int V dsigma ] over positive measures of mass lam,
for the field of a general profile:  V(t) = 2 pi sqrt t - int_0^1 eps(u) log(t + u^2) du.
Then  log Delta_K(zeta5) = kappa K^2 log K + R K^2 + O(K log K),  kappa = 2 lam^2 + 2 lam int eps."""
import numpy as np
from profiles import PROFILES

def Fanti(t, u):
    """int_0^u log(t + v^2) dv."""
    t = np.asarray(t, float); s = np.sqrt(t)
    with np.errstate(divide="ignore", invalid="ignore"):
        at = np.where(s > 0, 2 * s * np.arctan(u / np.where(s > 0, s, 1)), 0.0)
        return u * np.log(t + u * u) - 2 * u + at if u > 0 else 0 * t

def make_V(prof):
    def V(t):
        t = np.asarray(t, float)
        out = 2 * np.pi * np.sqrt(t)
        for u0, u1, e in prof:
            out = out - float(e) * (Fanti(t, float(u1)) - Fanti(t, float(u0)))
        return out
    return V

def _G(x):
    ax = np.abs(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ax > 0, x * x * np.log(ax) / 2 - 3 * x * x / 4, 0.0)

def solve(V, lam, n=1600, tmax=2.0):
    e = np.concatenate([[0], np.geomspace(2e-4, np.sqrt(tmax), n)]) ** 2
    p, q = e[:-1], e[1:]; w = q - p
    P, Q, R_, S = p[:, None], q[:, None], p[None, :], q[None, :]
    E = (_G(Q - R_) - _G(Q - S) - _G(P - R_) + _G(P - S)) / (w[:, None] * w[None, :])
    g, wq = np.polynomial.legendre.leggauss(16)
    v = (V((p + q)[:, None] / 2 + w[:, None] / 2 * g[None, :]) * wq).sum(1) / 2
    act = np.ones(len(v), bool)
    for _ in range(2000):
        idx = np.where(act)[0]
        M = np.zeros((len(idx) + 1,) * 2)
        M[:-1, :-1] = 2 * E[np.ix_(idx, idx)]; M[:-1, -1] = -1; M[-1, :-1] = 1
        sol = np.linalg.solve(M, np.concatenate([v[idx], [lam]]))
        x = np.zeros(len(v)); x[idx] = sol[:-1]; Fc = sol[-1]
        if (x < 0).any():
            act[np.argmin(x)] = False; continue
        gg = 2 * E @ x - v
        viol = (~act) & (gg > Fc + 1e-12)
        if viol.any():
            act[np.argmax(np.where(viol, gg, -np.inf))] = True; continue
        break
    supp = x > 1e-14
    return x @ E @ x - v @ x, (e[:-1][supp].min(), e[1:][supp].max())

if __name__ == "__main__":
    base = None
    for name, (prof, lam) in PROFILES.items():
        V = make_V(prof); lamf = float(lam)
        kappa = 2 * lamf**2 + 2 * lamf * float(sum((u1 - u0) * e for u0, u1, e in prof))
        R, (lo, hi) = solve(V, lamf)
        base = R if base is None else base
        print(f"{name:28s} lam={lamf:.4f}  kappa={kappa:+.4f}  R = {R:.5f}  (vs paper {R-base:+.5f})  supp [{lo:.1e}, {hi:.3f}]", flush=True)
