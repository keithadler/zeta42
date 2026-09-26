"""Generalised outer bound for K/2 < p <= K (K < 2p), any profile e (dict j -> exponent, -1 = pole)
and any number of rows h >= #poles.  Per class a = 1..(p-1)/2 (members j <= K with j = +-a mod p):
   ell   = #members,  delta = [a is a member and not a pole],  e_a = exponent of the small member a
   rows  = #poles in the class = ell - delta
   early rows i < ell-2 : w = min(0, i + (e_a+1)*delta/2 - (ell+4)/2)       (paper: 3*delta = 6*delta/2)
   late rows            : w = 0
 zero class (j = p, a pole): one row, w = -1/2.   Extra rows (h - #poles, basis D_tail t^k): w = 0.
 gamma = 2*sum(w) - min(r, z),  z = #zero weights,  r = max(0, 2h + 2 - 2p + deg W - #poles)."""
from fractions import Fraction as F
def gout(e, K, h, p):
    assert K < 2*p <= 2*K
    poles = {j for j, v in e.items() if v == -1 and j > 0}
    degW = sum(v for j, v in e.items() if v > 0 and j > 0)
    tot = F(0); z = 0; nrows = 0
    for a in range(1, (p-1)//2 + 1):
        mem = [j for j in range(1, K+1) if j % p in (a, p-a)]
        ell = len(mem)
        delta = 1 if (a in mem and a not in poles) else 0
        ea = e.get(a, 0) if delta else 0
        npo = sum(1 for j in mem if j in poles)
        assert npo == ell - delta, (p, a, mem)
        for i in range(npo):
            if i < ell - 2:
                w = min(F(0), i + F((ea+1)*delta, 2) - F(ell+4, 2))
            else:
                w = F(0)
            tot += w; z += (w == 0); nrows += 1
    if p in poles:
        tot += F(-1, 2); nrows += 1
    extra = h - len(poles); assert extra >= 0 and nrows == len(poles)
    z += extra
    r = max(0, 2*h + 2 - 2*p + degW - len(poles))
    return 2*tot - min(r, z)
if __name__ == "__main__":
    import sys; sys.path.insert(0, "../../arithmetic")
    from gamma_out import gamma_parts
    for n in (1, 2, 3):
        K, N = 40*n, 3*n; h = K - N
        e = {j: 5 for j in range(1, N+1)}; e.update({j: -1 for j in range(N+1, K+1)})
        for p in (q for q in range(K//2+1, K+1) if all(q % d for d in range(2, int(q**.5)+1))):
            base, r, z = gamma_parts(K, N, p); paper = base - min(r, z)
            g = gout(e, K, h, p)
            assert g == paper, (K, p, g, paper)
    print("generalised outer bound reproduces the paper's gamma_out for K/2 < p <= K at K = 40, 80, 120")
