"""The external field V of Fauzan's real estimate and nested-arcsine potentials (floats).

Used only for searching; every accepted number is re-checked in `certify.py`.
"""
import numpy as np

ALPHA, LAM = 3 / 40, 37 / 40
CSTAR = 2.653035990340          # C_* of the working edition, eq. (Cstar)


def V(t):
    t = np.asarray(t, float)
    s = np.sqrt(t)
    return (np.log1p(t) - 6 * ALPHA * np.log(t + ALPHA**2) - 2 + 12 * ALPHA
            + 2 * s * (np.pi + np.arctan(1 / s) - 6 * np.arctan(ALPHA / s)))


def U_arcsine(t, a, b):
    """Log potential of the unit arcsine measure on [a, b] (paper's sign)."""
    m = (a + b) / 2
    inside = (t >= a) & (t <= b)
    out = np.log((np.abs(t - m) + np.sqrt(np.abs((t - a) * (t - b)))) / 2)
    return np.where(inside, np.log((b - a) / 4), out)


def energy_nested(a, b, c):
    """I(rho) for nested supports: sum (S_j^2 - S_{j-1}^2) log((b_j - a_j)/4)."""
    S = np.cumsum(c)
    return np.sum((S**2 - np.concatenate([[0], S[:-1]])**2) * np.log((b - a) / 4))


def real_constant(a, b, c, t):
    """lam * max_t (2U - V) - I(rho) + C_*  on the grid t (float estimate of U)."""
    Ur = sum(cj * U_arcsine(t, aj, bj) for aj, bj, cj in zip(a, b, c))
    M = np.max(2 * Ur - V(t))
    return LAM * M - energy_nested(a, b, c) + CSTAR, M
