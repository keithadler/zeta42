"""Numerator zeros at half-integers: W = prod_{i odd <= M} (4t + i^2)^r, poles at t = -j^2, 1 <= j <= K."""
import profile as P
from flint import fmpz_poly
def score_W(s, K, h, Wextra):
    orig = P.fmpz_poly
    e = {j: -1 for j in range(1, K+1)}
    # inject extra numerator by wrapping build: multiply W inside via a sentinel profile entry
    a_b = None
    old_build = P.build
    def build(s, e, h):
        poles = sorted(j for j, v in e.items() if v == -1)
        W = Wextra; Den = fmpz_poly([1])
        for j in poles: Den *= fmpz_poly([j*j, 1])
        import math
        from flint import fmpq, fmpq_poly, fmpz
        E = 2*h-1; maxk = max(0, W.degree()+E-1-Den.degree()); sf = math.factorial(s-1)
        mus = [fmpq(-1)**k*fmpq.bernoulli(2*k+2)*fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf) for k in range(maxk+1)]
        Hs = {}; acc = fmpq(0)
        for j in poles: acc += fmpq(1, j**s); Hs[j] = acc
        base = {}
        for j in poles:
            den = fmpq(1)
            for k in poles:
                if k != j: den *= (k*k-j*j)
            base[j] = fmpq(W(fmpz(-j*j)))/den
        Wq, Dq = fmpq_poly(W), fmpq_poly(Den); a, b = [], []
        for ee in range(E):
            Pe = (Wq*fmpq_poly([0]*ee+[1]))//Dq
            ae = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0)); be = fmpq(0)
            for j in poles:
                c = base[j]*fmpq(-j*j)**ee; js = fmpq(j)**(s-1)
                be += c*js; ae += c*(-js*Hs[j]-fmpq(1, s-1)+fmpq(1, 2*j))
            a.append(ae); b.append(be)
        return a, b
    P.build = build
    try: return P.score(s, e, h)
    finally: P.build = old_build
K = 40
for M, r in [(0, 0), (3, 2), (5, 2), (5, 4), (9, 2), (9, 3), (15, 1), (15, 2), (21, 1)]:
    W = fmpz_poly([1])
    for i in range(1, M+1, 2): W *= fmpz_poly([i*i, 4])**r
    lD, lc, lP = score_W(7, K, K, W)
    print(f"s=7 poles 1..{K} h={K}, W=prod_(i odd<={M}) (4t+i^2)^{r}: logP/h^2 {lP/K**2:+.4f}  (logD {lD/K**2:+.4f}, -logc {-lc/K**2:+.4f})", flush=True)
