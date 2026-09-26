"""Generalised inner bound, any profile e and rows h, prime p with p^2 > 5K, p > max small j.
Classes a = 1..m (m=(p-1)/2): ell_a = #{j <= K : j = +-a mod p},
   b_a = sum over non-pole members j of (e_j + 1)/2      (paper: 3*ell_N(a))
Zero class: m_K = floor(K/p) poles (no small members since p > max small j).
Half-bounds of a row (the paper's source table, eqs. ordinarysource/zerosource):
   row (a,i) at its own source a :  i + b_a - (ell_a+4)/2
   any row at source c (not own) :  Phi_c = L_c + b_c - (ell_c+4)/2      (zero of order L_c there)
   ordinary row at zero source   :  2 L_0 - m_K + 1/2
   zero row (0,i) at zero source :  2 i - m_K + 1/2
Weight of a row = min of its half-bounds over all sources (valid by the ultrametric inequality);
gamma = 2 * sum of weights.  L = (L_0, L_1..L_m), sum = h, is free: the paper's allocation is one
choice; `optimise` improves it by single-row moves.  Degree conditions (p large) are checked
by the caller."""
from fractions import Fraction as F

def classes(e, K, p):
    m = (p - 1)//2
    poles = {j for j, v in e.items() if v == -1 and j > 0}
    ell, b = [0]*(m+1), [F(0)]*(m+1)
    nz = F(0)                      # zero-class numerator: sum of (e_j+1) over non-pole j = 0 mod p
    for j in range(1, K+1):
        r = j % p
        if r == 0:
            if j not in poles: nz += e.get(j, 0) + 1
            continue
        a = min(r, p - r)
        ell[a] += 1
        if j not in poles: b[a] += F(e.get(j, 0) + 1, 2)
    mK = sum(1 for j in range(p, K+1, p))
    return m, ell, b, (mK, nz)

def gamma(L, ell, b, mz):
    mK, nz = mz
    m = len(L) - 1
    Phi = [None] + [L[c] + b[c] - F(ell[c] + 4, 2) for c in range(1, m+1)]
    order = sorted(range(1, m+1), key=lambda c: Phi[c])
    lo1, lo2 = order[0], order[1] if m > 1 else None
    zsrc = 2*L[0] + nz - mK + F(1, 2)
    tot = F(0)
    for a in range(1, m+1):
        other = Phi[lo2] if a == lo1 else Phi[lo1]
        cap = min(other, zsrc)
        base = b[a] - F(ell[a] + 4, 2)
        for i in range(L[a]): tot += min(i + base, cap)
    capz = Phi[lo1]
    for i in range(L[0]): tot += min(2*i + nz - mK + F(1, 2), capz)
    return 2*tot

def paper_alloc(ell, b, h, L0):
    m = len(ell) - 1
    tot = h - L0 + sum(b[1:])
    T, E = divmod(tot, m)
    order = sorted(range(1, m+1), key=lambda a: (-ell[a], a))
    L = [L0] + [0]*m
    for k, a in enumerate(order): L[a] = int(T - b[a] + (1 if k < E else 0))
    return L if min(L) >= 0 and sum(L) == h else None

def start_alloc(ell, b, h, L0):
    """Integer allocation equalising Phi_c = L_c + b_c - (ell_c+4)/2 as far as possible."""
    m = len(ell) - 1
    L = [L0] + [0]*m
    for _ in range(h - L0):
        c = min(range(1, m+1), key=lambda c: (L[c] + b[c] - F(ell[c] + 4, 2), c))
        L[c] += 1
    return L

def optimise(L, ell, b, mK, iters=100000):
    L = list(L); g = gamma(L, ell, b, mK); m = len(L) - 1
    improved = True
    while improved:
        improved = False
        for s in range(m+1):
            if L[s] == 0: continue
            for d in range(m+1):
                if d == s: continue
                L[s] -= 1; L[d] += 1
                g2 = gamma(L, ell, b, mK)
                if g2 > g: g = g2; improved = True; break
                L[s] += 1; L[d] -= 1
            if improved: break
    return L, g
