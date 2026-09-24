"""Diagnostic only (not valid constructions): which ingredient carries the arithmetic cost?"""
import math
from math import gcd, lcm
from flint import fmpq, fmpq_mat
import profile as P
orig_bern = P.fmpq.bernoulli
def content_only(s, e, h, killH=False, killB=False, killC=False):
    src = P.build
    a, b = [], []
    # re-run build with patched pieces
    import types
    code = P.build.__code__
    g = dict(P.build.__globals__)
    return None
def run(s, K, killH, killB, killC):
    e = {j: -1 for j in range(1, K+1)}; h = K
    poles = list(range(1, K+1))
    import math
    from flint import fmpz_poly, fmpq_poly, fmpz
    Den = fmpz_poly([1])
    for j in poles: Den *= fmpz_poly([j*j, 1])
    E = 2*h-1; maxk = max(0, E-1-Den.degree()); sf = math.factorial(s-1)
    mus = [fmpq(0) if killB else fmpq(-1)**k*fmpq.bernoulli(2*k+2)*fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf) for k in range(maxk+1)]
    Hs = {}; acc = fmpq(0)
    for j in poles: acc += fmpq(1, j**s); Hs[j] = fmpq(0) if killH else acc
    base = {}
    for j in poles:
        den = fmpq(1)
        for k in poles:
            if k != j: den *= (k*k-j*j)
        base[j] = 1/den
    a, b = [], []
    for ee in range(E):
        Pe = fmpq_poly([0]*ee+[1]) // fmpq_poly(Den)
        ae = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0)); be = fmpq(0)
        for j in poles:
            c = base[j]*fmpq(-j*j)**ee; js = fmpq(j)**(s-1)
            be += c*js
            ae += c*(-js*Hs[j] + (0 if killC else (-fmpq(1, s-1)+fmpq(1, 2*j))))
        a.append(ae); b.append(be)
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    D = (-(B.inv()*A)).charpoly()*B.det()
    g_, l_ = 0, 1
    for i in range(D.degree()+1):
        c = D[i]; g_ = gcd(g_, int(c.p)); l_ = lcm(l_, int(c.q))
    return -(math.log(g_)-math.log(l_))/K**2
for s in (5, 7):
    K = 40
    print(f"s={s}: -log content/K^2  full {run(s,K,0,0,0):+.4f} | H:=0 {run(s,K,1,0,0):+.4f} | Bernoulli:=0 {run(s,K,0,1,0):+.4f} | const:=0 {run(s,K,0,0,1):+.4f} | only residues (all 0) {run(s,K,1,1,1):+.4f}", flush=True)
