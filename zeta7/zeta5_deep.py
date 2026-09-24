"""Deeper zeta(5) search. Stage 1 (n=1, fast): hill-climb over 12 band exponents (-1..8), an optional
factor t^m (m = 0..2), and alpha in [1.0, 1.35]. Stage 2: re-score the 8 best distinct profiles at n=2
(bands scaled: band i covers j in ((i-1)n, in]), since n=1 can flatter a profile."""
import sys, time, multiprocessing as mp
from manyrows import score
B = 12
def build(ex, m0, n):
    K = 40*n; e = {}
    for i, v in enumerate(ex, start=1):
        for j in range((i-1)*n+1, i*n+1):
            if v: e[j] = v
    for j in range(B*n+1, K+1): e[j] = -1
    if m0: e[0] = m0
    return e, K
def ev(arg):
    ex, m0, a, n = arg
    e, K = build(ex, m0, n); m = sum(1 for v in e.values() if v == -1)
    if m < 26*n: return None
    try: return score(5, e, round(a*m))[1] / K**2
    except Exception: return None
def nbrs(ex, m0, a):
    out = []
    for i in range(B):
        for d in (-1, 1):
            v = ex[i] + d
            if -1 <= v <= 8:
                x = list(ex); x[i] = v; out.append((tuple(x), m0, a))
    for m2 in (m0-1, m0+1):
        if 0 <= m2 <= 2: out.append((ex, m2, a))
    for a2 in (round(a-0.05, 2), round(a+0.05, 2)):
        if 1.0 <= a2 <= 1.35: out.append((ex, m0, a2))
    return out
if __name__ == "__main__":
    pool = mp.Pool(4); seen = {}
    pad = lambda t: tuple(t) + (-1,)*(B-len(t))
    starts = [(pad((4, 4, 3, 2, 0, 0)), 0, 1.2), (pad((5, 5, 3, 1)), 0, 1.1),
              (pad((5, 5, 5)), 1, 1.1), (pad((6, 5, 4, 3, 2, 1, 0, 0)), 0, 1.2)]
    t0 = time.time()
    for st in starts:
        cur = st; cv = ev((*cur, 1)); seen[cur] = cv
        while True:
            nb = [x for x in nbrs(*cur) if x not in seen]
            vals = pool.map(ev, [(*x, 1) for x in nb])
            for x, v in zip(nb, vals): seen[x] = v
            c = [(v, x) for x, v in zip(nb, vals) if v is not None]
            if not c or min(c)[0] >= cv: break
            cv, cur = min(c)
        print(f"stage1 start {st} -> {cur}  n=1 log P = {cv*1600:+.2f}  ({time.time()-t0:.0f}s)", flush=True)
    top = sorted((v, k) for k, v in seen.items() if v is not None)
    uniq = []
    for v, k in top:
        if k not in [u[1] for u in uniq]: uniq.append((v, k))
        if len(uniq) == 8: break
    print(f"stage1 evaluated {len(seen)}; top 8 at n=1:", flush=True)
    for v, k in uniq: print(f"  n=1 log P = {v*1600:+.2f}  e={k[0]} t^{k[1]} alpha={k[2]}", flush=True)
    res = pool.map(ev, [(*k, 2) for v, k in uniq])
    print("stage2 at n=2 (paper -833.33, extra rows -923.81, best banded so far -1064.40):", flush=True)
    for (v1, k), v2 in sorted(zip(uniq, res), key=lambda z: z[1] if z[1] is not None else 9):
        if v2 is not None:
            print(f"  n=2 log P = {v2*6400:+.2f}  (n=1 {v1*1600:+.2f})  e={k[0]} t^{k[1]} alpha={k[2]}", flush=True)
