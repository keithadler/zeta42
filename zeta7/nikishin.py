"""'Cube' idea: mixed (multiple-orthogonality) determinants over a Nikishin-type system.
Row types r: functions t^i * g_r(t), g_0 = R (base rational weight), g_1 = f*R, g_2 = f2*f*R, ...
with f(t) = sum_{j in S} c_j/(t+j^2), c_j > 0 (Markov function of a measure on the negative axis, poles
disjoint from R's). Columns t^k, k < N = sum n_r. M_{(r,i),k} = mu_X(t^{i+k} g_r), linear in X.
By Andreief, det M = (1/N!) int det[phi_r(t_a)] det[t_a^k] prod dmu, and for an AT (Chebyshev) system both
factors have fixed sign on ordered points, so det M(zeta s) != 0. We also verify nonvanishing numerically."""
import math, sys
from math import gcd, lcm
from flint import fmpq, fmpz, fmpz_poly, fmpq_poly, fmpz_mat, arb, ctx

def mu_moments(s, W, poles, count):
    """mu_X(t^e W/prod_{j in poles}(t+j^2)) for e < count, as (const, Xcoef) rationals."""
    poles = sorted(poles)
    Den = fmpz_poly([1])
    for j in poles: Den *= fmpz_poly([j*j, 1])
    sf = math.factorial(s-1)
    maxk = max(0, W.degree() + count - 1 - Den.degree())
    mus = [fmpq(-1)**k*fmpq.bernoulli(2*k+2)*fmpq(math.factorial(2*k+s), math.factorial(2*k+2)*sf) for k in range(maxk+1)]
    H = [fmpq(0)]
    for k in range(1, max(poles)+1): H.append(H[-1] + fmpq(1, k**s))
    base = {}
    for j in poles:
        den = fmpq(1)
        for k in poles:
            if k != j: den *= (k*k - j*j)
        base[j] = fmpq(W(fmpz(-j*j))) / den
    Wq, Dq = fmpq_poly(W), fmpq_poly(Den)
    out = []
    for e in range(count):
        Pe = (Wq*fmpq_poly([0]*e+[1])) // Dq
        a = sum((Pe[k]*mus[k] for k in range(Pe.degree()+1)), fmpq(0)); b = fmpq(0)
        for j in poles:
            c = base[j]*fmpq(-j*j)**e; js = fmpq(j)**(s-1)
            b += c*js; a += c*(-js*H[j] - fmpq(1, s-1) + fmpq(1, 2*j))
        out.append((a, b))
    return out

def markov_times(W, poles, fpoles, fc):
    """(W/prod poles) * sum_j fc_j/(t+j^2) as a single numerator over poles U fpoles."""
    num = fmpz_poly([0])
    Fden = fmpz_poly([1])
    for j in fpoles: Fden *= fmpz_poly([j*j, 1])
    for j, c in zip(fpoles, fc):
        num += c * (Fden // fmpz_poly([j*j, 1]))
    return W*num, sorted(set(poles) | set(fpoles))

def score(s, types, N):
    """types: list of (W, poles, n_rows). Returns (degree, log|P(zeta s)|)."""
    rows_a, rows_b = [], []
    for W, poles, nr in types:
        m = mu_moments(s, W, poles, nr + N - 1)
        for i in range(nr):
            rows_a.append([m[i+k][0] for k in range(N)]); rows_b.append([m[i+k][1] for k in range(N)])
    assert len(rows_a) == N
    L = 1
    for r in rows_a + rows_b:
        for c in r: L = lcm(L, int(c.q))
    Ai = [[int(c.p)*(L//int(c.q)) for c in r] for r in rows_a]
    Bi = [[int(c.p)*(L//int(c.q)) for c in r] for r in rows_b]
    npoles = len(set().union(*[set(p) for _, p, _ in types]))
    d = min(N, npoles)
    xs = list(range(d+1))
    ys = [fmpz_mat(N, N, [Ai[i][k] + x*Bi[i][k] for i in range(N) for k in range(N)]).det() for x in xs]
    coef = [fmpq(int(y)) for y in ys]
    for k in range(1, d+1):
        for i in range(d, k-1, -1): coef[i] = (coef[i]-coef[i-1])/(xs[i]-xs[i-k])
    D = fmpq_poly([coef[d]])
    for i in range(d-1, -1, -1): D = D*fmpq_poly([-xs[i], 1]) + coef[i]
    co = [D[i] for i in range(D.degree()+1)]
    g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    P = fmpz_poly([int(c.p)*(l//int(c.q))//g for c in co])
    bits = max(abs(int(P[i])).bit_length() for i in range(P.degree()+1))
    ctx.prec = 2*bits + 4000
    v = P(arb(s).zeta())
    assert not v.contains(0), "determinant vanishes at zeta(s) (or not certified nonzero)"
    return P.degree(), float(abs(v).log())

if __name__ == "__main__":
    one = fmpz_poly([1])
    for s in (5, 7):
        K = 20
        base = (one, list(range(1, K+1)))
        # control: single measure, N rows = Hankel Gram (known: s=5 -0.0828, s=7 +0.8200 at N=20)
        d, lp = score(s, [(one, base[1], 20)], 20)
        print(f"s={s} control Gram, 20 poles, N=20: deg {d} logP/N^2 {lp/400:+.4f}", flush=True)
        for L2, split in [(10, (14, 6)), (10, (10, 10)), (20, (14, 6)), (20, (10, 10))]:
            fp = list(range(K+1, K+L2+1)); fc = [1]*L2
            W1, p1 = markov_times(one, base[1], fp, fc)
            d, lp = score(s, [(one, base[1], split[0]), (W1, p1, split[1])], 20)
            print(f"s={s} Nikishin: R=1/D_{K}, f=sum_(j={K+1}..{K+L2}) 1/(t+j^2), rows {split}: deg {d} logP/N^2 {lp/400:+.4f}", flush=True)
