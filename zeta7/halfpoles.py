"""zeta(7) with poles at half-integers: points a = j/2 (j in half-units), t = y^2, measure w(y) R(t),
R = prod (t + j^2/4)^{e_j}. Pole moments (Hermite for real a, checked to 40 digits):
  j = 2i even: i^6 (X - H_i) - 1/6 + 1/j
  j odd      : j^6 ((127/64) X - 2 O_j) - 1/6 + 1/j,   O_j = sum_{m odd <= j} 1/m^7
Equivalent to the weight y^7 F^(6)(y/2) with integer poles."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpq_poly, fmpq_mat, arb, ctx
S = 7
def build(e, h):
    poles = sorted(j for j, v in e.items() if v == -1)
    W = fmpq_poly([1]); Den = fmpq_poly([1])
    for j, v in e.items():
        f = fmpq_poly([fmpq(j*j, 4), 1])
        if v > 0: W *= f**v
        elif v == -1: Den *= f
    E = 2*h - 1; sf = math.factorial(S-1)
    maxk = max(0, W.degree() + E - 1 - Den.degree())
    mus = [fmpq(-1)**k*fmpq.bernoulli(2*k+2)*fmpq(math.factorial(2*k+S), math.factorial(2*k+2)*sf) for k in range(maxk+1)]
    H = [fmpq(0)]; O = [fmpq(0)]
    for k in range(1, max(poles)+1):
        H.append(H[-1] + fmpq(1, k**S)); O.append(O[-1] + (fmpq(1, k**S) if k % 2 else 0))
    kap, nu = {}, {}
    for j in poles:
        if j % 2 == 0:
            i = j//2; kap[j] = fmpq(i**6); nu[j] = -fmpq(i**6)*H[i]
        else:
            kap[j] = fmpq(j**6*127, 64); nu[j] = -2*fmpq(j**6)*O[j]
        nu[j] += -fmpq(1, 6) + fmpq(1, j)
    base = {}
    for j in poles:
        q = fmpq(-j*j, 4); den = fmpq(1)
        for k in poles:
            if k != j: den *= (fmpq(k*k, 4) + q)
        base[j] = W(q)/den
    a, b = [], []
    for ee in range(E):
        Pe = (W*fmpq_poly([0]*ee+[1])) // Den
        ae = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0)); be = fmpq(0)
        for j in poles:
            c = base[j]*fmpq(-j*j, 4)**ee
            be += c*kap[j]; ae += c*nu[j]
        a.append(ae); b.append(be)
    return a, b
def score(e, h):
    a, b = build(e, h)
    m = sum(1 for v in e.values() if v == -1); d = min(m, h)
    A = fmpq_mat(h, h, [a[i+j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i+j] for i in range(h) for j in range(h)])
    if d == h and B.rank() == h:
        D = (-(B.inv()*A)).charpoly()*B.det()
    else:
        xs = list(range(d+1)); ys = [(A + x*B).det() for x in xs]; D = fmpq_poly([0])
        for i, xi in enumerate(xs):
            L = fmpq_poly([1])
            for k, xk in enumerate(xs):
                if k != i: L = L*fmpq_poly([-xk, 1])/(xi - xk)
            D += L*ys[i]
    co = [D[i] for i in range(D.degree()+1)]; g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    bits = max(max(abs(int(c.p)).bit_length(), int(c.q).bit_length()) for c in co)
    ctx.prec = 2*bits + 3000
    z = arb(S).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= z
    assert v > 0
    return D.degree(), float(v.log()) - (math.log(g) - math.log(l))
if __name__ == "__main__":
    for m in (20, 40):
        sets = {f"integer poles 1..{m} (control)": [2*i for i in range(1, m+1)],
                f"half-units 1..{m} (a=0.5..{m/2})": list(range(1, m+1)),
                f"odd half-units only (a=0.5,1.5,..)": list(range(1, 2*m, 2))}
        for name, J in sets.items():
            e = {j: -1 for j in J}
            deg, lP = score(e, m)
            print(f"{m} poles, {name:36s} h={m}: log P/m^2 {lP/m**2:+.4f}", flush=True)
