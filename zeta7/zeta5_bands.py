"""zeta(5) banding: per-j numerator exponents e_1..e_B on the first B values of j (each e_j in
-1..8; -1 = pole), poles for B < j <= K, rows h = alpha * #poles. Hill-climb from the paper's
choice (e = 5 on j <= 3). Score: log P(zeta5) / K^2 at fixed K (paper at n=1: -265.13/1600)."""
import sys, multiprocessing as mp
from manyrows import score
K, B = 40, 8
def profile(ex, K):
    e = {j: v for j, v in enumerate(ex, start=1) if v != 0}
    for j in range(len(ex)+1, K+1): e[j] = -1
    return e
def ev(arg):
    ex, a = arg
    e = profile(ex, K); m = sum(1 for v in e.values() if v == -1)
    if m < 30: return None
    try: return score(5, e, round(a*m))[1] / K**2
    except Exception: return None
def nbrs(ex, a):
    out = []
    for i in range(B):
        for d in (-1, 1):
            v = ex[i] + d
            if -1 <= v <= 8:
                x = list(ex); x[i] = v; out.append((tuple(x), a))
    for a2 in (round(a - 0.05, 2), round(a + 0.05, 2)):
        if 1.0 <= a2 <= 1.2: out.append((ex, a2))
    return out
if __name__ == "__main__":
    pool = mp.Pool(4); seen = {}
    starts = [((5, 5, 5, -1, -1, -1, -1, -1), 1.0), ((5, 5, 5, -1, -1, -1, -1, -1), 1.1),
              ((6, 5, 4, 3, -1, -1, -1, -1), 1.1), ((4, 4, 4, 4, -1, -1, -1, -1), 1.1)]
    for st in starts:
        cur = st; cv = ev(cur); seen[cur] = cv
        while True:
            nb = [n for n in nbrs(*cur) if n not in seen]
            vals = pool.map(ev, nb)
            for n, v in zip(nb, vals): seen[n] = v
            c = [(v, n) for n, v in zip(nb, vals) if v is not None]
            if not c or min(c)[0] >= cv: break
            cv, cur = min(c)
            print(f"  step -> {cur}  {cv:+.4f}  (log P = {cv*K*K:+.2f})", flush=True)
        print(f"start {st} -> {cur}  score {cv:+.4f}  (log P = {cv*K*K:+.2f})", flush=True)
    top = sorted((v, k) for k, v in seen.items() if v is not None)[:6]
    print(f"evaluated {len(seen)}; best:")
    for v, k in top: print(f"  {v:+.4f}  (log P = {v*K*K:+.2f})  e={k[0]} alpha={k[1]}")
