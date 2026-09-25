"""Exact inner valuation bound gamma_p^in (working edition, eqs. allocation, innerweights,
zeroweights, gammain), with the reserved zero block L0 as a free parameter.
Returns None if the allocation makes some L_a < 0."""
from fractions import Fraction as F
def ellA(A, a, p): return sum(1 for j in range(1, A+1) if j % p in (a, p-a))
def gamma_in(K, N, p, L0):
    h = K - N; m = (p-1)//2; mN, mK = N//p, K//p
    lK = {a: ellA(K, a, p) for a in range(1, m+1)}
    b = {a: 3*ellA(N, a, p) for a in range(1, m+1)}
    tot = h - L0 + 3*(N - mN)
    T, E = divmod(tot, m)
    order = sorted(range(1, m+1), key=lambda a: (-lK[a], a))
    eps = {a: 0 for a in order}
    for a in order[:E]: eps[a] = 1
    L = {a: T - b[a] + eps[a] for a in order}
    if min(L.values()) < 0: return None
    # hypotheses of the inner proposition (working edition), made explicit:
    if not (p >= 7 and p*p > 5*K): return None
    if not (2*(T + 1) < p + 1 and 5 + 4*L0 + 12*mN < p + 1): return None          # degree checks
    wmax = max((F(i + b[a]) - F(lK[a] + 4, 2) for a in order for i in range(L[a])), default=None)
    if wmax is not None and not (2*L0 + 6*mN - mK + F(1, 2) >= wmax): return None  # zero source dominates
    Z = {a: T + eps[a] for a in order}
    assert L0 + sum(L.values()) == h
    g = F(0)
    for a in order:
        for i in range(L[a]): g += i + b[a] - F(lK[a] + 4, 2)
    zc = min(Z[c] - F(lK[c] + 4, 2) for c in order)
    for i in range(L0): g += min(2*i + 6*mN - mK + F(1, 2), zc)
    return 2*g
