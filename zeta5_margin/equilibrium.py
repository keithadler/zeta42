"""The best real constant any comparison measure can give: the weighted equilibrium.

For any positive rho of mass lam, the paper's bound is U(rho) = lam*max(2U^rho - V) - I(rho) + C_*,
and U(rho) >= U* := C_* + sup_sigma [I(sigma) - int V dsigma], with equality at the
equilibrium measure.  We compute U* with piecewise-constant densities (exact cell-cell
log energies) and an active-set solve of the KKT system  2(E x)_k - v_k = F on the support.

    python equilibrium.py      ->  U* ~ -1.38544  (paper's measure: -1.366996)
"""
import numpy as np
from field import V, LAM, CSTAR


def _G(x):
    ax = np.abs(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ax > 0, x * x * np.log(ax) / 2 - 3 * x * x / 4, 0.0)


def solve(n):
    e = np.concatenate([[0], np.geomspace(2e-4, np.sqrt(1.2), n)]) ** 2
    p, q = e[:-1], e[1:]
    w = q - p
    P, Q, R, S = p[:, None], q[:, None], p[None, :], q[None, :]
    E = (_G(Q - R) - _G(Q - S) - _G(P - R) + _G(P - S)) / (w[:, None] * w[None, :])
    g, wq = np.polynomial.legendre.leggauss(16)
    v = (V((p + q)[:, None] / 2 + w[:, None] / 2 * g[None, :]) * wq).sum(1) / 2
    act = np.ones(len(v), bool)
    for _ in range(500):
        idx = np.where(act)[0]
        M = np.zeros((len(idx) + 1,) * 2)
        M[:-1, :-1] = 2 * E[np.ix_(idx, idx)]
        M[:-1, -1] = -1
        M[-1, :-1] = 1
        sol = np.linalg.solve(M, np.concatenate([v[idx], [LAM]]))
        x = np.zeros(len(v))
        x[idx] = sol[:-1]
        F = sol[-1]
        if (x < 0).any():
            act[np.argmin(x)] = False
            continue
        g_ = 2 * E @ x - v
        viol = (~act) & (g_ > F + 1e-12)
        if viol.any():
            act[np.argmax(np.where(viol, g_, -np.inf))] = True
            continue
        break
    return x @ E @ x - v @ x + CSTAR


if __name__ == "__main__":
    for n in (200, 400, 800, 1600):
        print(f"{n:5d} cells: U* ~ {solve(n):.6f}", flush=True)
