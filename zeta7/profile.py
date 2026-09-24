"""Hankel construction for zeta(s) with an arbitrary exponent profile.
Measure: w(y) * prod_j (t + j^2)^{e_j}, t = y^2, with e_j in {-1, 0, 1, 2, ...}.
Poles (e_j = -1) carry X = zeta(s); numerator factors (e_j > 0) may sit below or above the poles.
h <= #poles rows. The primitive P of Delta(X) is basis-independent, so this family is the whole knob set."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpz, fmpz_poly, fmpq_poly, fmpq_mat, arb, ctx

def build(s, e, h):
    """e: dict j -> exponent. Returns Hankel data a (const part), b (X part) for indices 0..2h-2."""
    poles = sorted(j for j, v in e.items() if v == -1)
    assert 0 < h <= len(poles)
    W = fmpz_poly([1]); Den = fmpz_poly([1])
    for j, v in e.items():
        if v > 0: W *= fmpz_poly([j*j, 1])**v
        elif v == -1: Den *= fmpz_poly([j*j, 1])
    E = 2*h - 1
    maxk = max(0, W.degree() + E - 1 - Den.degree())
    sf = math.factorial(s-1)
    mus = [fmpq(-1)**k * fmpq.bernoulli(2*k+2) * fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf)
           for k in range(maxk+1)]
    Hs = {}; acc = fmpq(0)
    for j in range(1, max(poles)+1):
        acc += fmpq(1, j**s); Hs[j] = acc
    base = {}
    for j in poles:
        den = fmpq(1)
        for k in poles:
            if k != j: den *= (k*k - j*j)
        base[j] = fmpq(W(fmpz(-j*j))) / den
    Wq, Dq = fmpq_poly(W), fmpq_poly(Den)
    a, b = [], []
    for ee in range(E):
        Pe = (Wq * fmpq_poly([0]*ee+[1])) // Dq
        ae = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0))
        be = fmpq(0)
        for j in poles:
            c = base[j] * fmpq(-j*j)**ee
            js = fmpq(j)**(s-1)
            be += c*js
            ae += c*(-js*Hs[j] - fmpq(1, s-1) + fmpq(1, 2*j))
        a.append(ae); b.append(be)
    return a, b

def score(s, e, h):
    """Returns (log Delta(zeta s), log content, log P(zeta s))."""
    a, b = build(s, e, h)
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    D = (-(B.inv()*A)).charpoly() * B.det()
    co = [D[i] for i in range(D.degree()+1)]
    g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    bits = max(max(abs(int(c.p)).bit_length(), int(c.q).bit_length()) for c in co)
    ctx.prec = 2*bits + 64*h + 2000
    z = arb(s).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= z
    assert v > 0, "Gram positivity failed"
    lD = float(v.log()); lc = math.log(g) - math.log(l)
    return lD, lc, lD - lc

def staircase(K, steps, poles, above=()):
    """steps: [(N, r)] adds r to e_j for j<=N; poles: (lo, hi] get -1; above: [(M, r)] adds r for K<j<=M."""
    e = {}
    for N, r in steps:
        for j in range(1, N+1): e[j] = e.get(j, 0) + r
    lo, hi = poles
    for j in range(lo+1, hi+1): e[j] = -1
    for M, r in above:
        for j in range(hi+1, M+1): e[j] = e.get(j, 0) + r
    return e

if __name__ == "__main__":
    # validation: paper's zeta(5) n=1 (D_3^6/D_40 = net exponent 5 on j<=3)
    print("zeta5 paper:", score(5, staircase(40, [(3, 5)], (3, 40)), 37)[2], "(expect -265.13)")
