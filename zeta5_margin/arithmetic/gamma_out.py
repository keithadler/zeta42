import json, sys
from slack import vS
def gamma_parts(K, N, p):
    v = K - p*(K//p); u = max(0, N + v - p + 1); t = min(N, v) + u
    r = max(0, K + 4*N - 2*p + 2)
    if K < 2*p: base = -7*(K-p) + 6*t - 1; z = p - 1 - N + u
    else:       base = -7*(K-p) + 3 + 12*N + 5*t; z = p + u
    return base, r, z
for K in map(int, sys.argv[1:]):
    N = 3*K//40; h = K - N
    rows = json.load(open(f"rows{K}.json"))
    print(f"K={K}  (y=p/K; vG = true v_p^G(Delta); 2Sw = paper's weight sum; loss_true = 2Sw - vG; paper loses min(r,z))")
    for p, true, _ in rows:
        if not (K/3 < p <= K): continue
        vG = true - vS(K, N, h, p)
        base, r, z = gamma_parts(K, N, p)
        print(f"  p={p:4d} y={p/K:.3f}  vG={vG:5d}  2Sw={base:5d}  r={r:4d} z={z:4d}  paper loss={min(r,z):4d}  true loss={base-vG:4d}")
