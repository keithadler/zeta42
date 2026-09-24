import sys, math, itertools
from math import gcd, lcm
from flint import arb, ctx
from hankel_s import entries, det_poly, Dpoly
def parts(s, K, N, r, W=None):
    h, a, b = entries(s, K, N, r, W); D = det_poly(h, a, b)
    co = [D[i] for i in range(h+1)]; g, l = 0, 1
    for c in co: g = gcd(g, int(c.p)); l = lcm(l, int(c.q))
    ctx.prec = 64*h + 4*K*K + 3000
    z = arb(s).zeta(); v = arb(0); zp = arb(1)
    for c in co: v += arb(int(c.p))/arb(int(c.q))*zp; zp *= z
    logD = float(v.log()); logc = math.log(g) - math.log(l)
    return logD, logc, logD - logc
def show(tag, K, res):
    lD, lc, lP = res
    print(f"{tag} | logDelta/K^2 {lD/K**2:+8.4f}  -logcontent/K^2 {-lc/K**2:+8.4f}  logP/K^2 {lP/K**2:+8.4f}", flush=True)
if __name__ == "__main__":
    s = int(sys.argv[1]); K = int(sys.argv[2])
    Ns = [int(x) for x in sys.argv[3].split(',')]; rs = [int(x) for x in sys.argv[4].split(',')]
    for N, r in itertools.product(Ns, rs):
        if N < K: show(f"s={s} K={K:3d} N={N:3d} r={r:2d}", K, parts(s, K, N, r))
