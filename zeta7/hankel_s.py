"""Generalised Fauzan-type Hankel construction for zeta(s), s odd.
Weight w(y) = y^s F^{(s-1)}(y) / ((s-1)!/2),  F = 1/(e^{2 pi y} - 1)  (positive on (0,oo)).
  mu(t^k)         = (-1)^k B_{2k+2} (2k+s)! / ((2k+2)! (s-1)!)
  mu(1/(t+j^2))   = j^{s-1} (X - H_j^{(s)}) - 1/(s-1) + 1/(2j),   X = zeta(s)
G(X)_{a,b} = mu_X( D_N(t)^r t^{a+b} / D_K(t) ),  a,b < h = K - N,  Delta(X) = det G(X), deg h.
Delta(zeta(s)) > 0 (Gram determinant of a positive weight). P = primitive part of Delta in Z[X].
Irrationality of zeta(s) follows if log P(zeta(s)) <= -c n^2 (degree h is linear in n)."""
import sys, time, math
from math import gcd, lcm
from flint import fmpq, fmpz, fmpz_poly, fmpq_poly, fmpq_mat, arb, arb_poly, ctx

def Dpoly(m):
    p = fmpz_poly([1])
    for j in range(1, m+1): p *= fmpz_poly([j*j, 1])
    return p

def entries(s, K, N, r, W=None):
    """W: optional numerator polynomial (must vanish on -j^2 for j<=N); default D_N^r."""
    h = K - N
    DN, DK = Dpoly(N), Dpoly(K)
    W = DN**r if W is None else W
    E = 2*h - 1
    maxk = max(0, W.degree() + E - 1 - K)
    sf = math.factorial(s-1)
    mus = [fmpq(-1)**k * fmpq.bernoulli(2*k+2) * fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf)
           for k in range(maxk+1)]
    Hs = [fmpq(0)]
    for j in range(1, K+1): Hs.append(Hs[-1] + fmpq(1, j**s))
    base = {}
    for j in range(N+1, K+1):
        num = fmpq(W(fmpz(-j*j)))
        den = fmpq(1)
        for k in range(1, K+1):
            if k != j: den *= (k*k - j*j)
        base[j] = num/den
    Wq, DKq = fmpq_poly(W), fmpq_poly(DK)
    a, b = [], []
    for e in range(E):
        Pe = (Wq * fmpq_poly([0]*e+[1])) // DKq
        ae = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0))
        be = fmpq(0)
        for j in range(N+1, K+1):
            c = base[j] * fmpq(-j*j)**e
            js = fmpq(j)**(s-1)
            be += c * js
            ae += c * (-js*Hs[j] - fmpq(1, s-1) + fmpq(1, 2*j))
        a.append(ae); b.append(be)
    return h, a, b

def det_poly(h, a, b):
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    return (-(B.inv() * A)).charpoly() * B.det()

def run(s, K, N, r, n=None):
    t0 = time.time()
    h, a, b = entries(s, K, N, r)
    D = det_poly(h, a, b)
    assert D.degree() == h
    co = [D[i] for i in range(h+1)]
    g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    P = [int(c.p) * (l // int(c.q)) // g for c in co]
    bits = max(abs(x).bit_length() for x in P)
    ctx.prec = bits + 64*h + 2000
    z = arb(s).zeta()
    Pz = arb_poly(P)(z)
    assert Pz > 0, Pz            # Gram positivity (sign is fixed by the construction)
    logP = float(Pz.log()); logH = bits*math.log(2)
    out = dict(s=s, K=K, N=N, r=r, h=h, logP=logP, logH=logH, secs=time.time()-t0)
    msg = f"s={s} K={K} N={N} r={r} h={h}: log P(zeta{s}) = {logP:10.2f}  log H(P) = {logH:10.2f}  /K^2: {logP/K**2:+.4f}"
    if n: msg += f"  /n^2: {logP/n**2:+.2f}"
    print(msg + f"  ({out['secs']:.1f}s)", flush=True)
    return out

if __name__ == "__main__":
    s, kK, kN, r = map(int, sys.argv[1:5])
    for n in map(int, sys.argv[5:]):
        run(s, kK*n, kN*n, r, n)
