"""Score = log P(zeta s) / h^2 (negative is good). Hill-climb over block exponent profiles for zeta(s). A shape is a tuple of block exponents in
{-1..6} over 1..K split into B equal blocks (-1 = poles), plus a row fraction hf = h/#poles."""
import random, sys, json, multiprocessing as mp
from profile import score
S = int(sys.argv[1]) if len(sys.argv) > 1 else 7
K, B = 40, 8
def profile(shape, K):
    bs = K // len(shape); e = {}
    for i, v in enumerate(shape):
        for j in range(i*bs+1, (i+1)*bs+1):
            if v: e[j] = v
    return e
def evaluate(arg):
    shape, hf, K = arg
    e = profile(shape, K); npoles = sum(1 for v in e.values() if v == -1)
    h = max(1, round(hf*npoles))
    if npoles < 15: return None
    try: lD, lc, lP = score(S, e, h)
    except Exception as ex: return None
    return lP / h**2   # size-invariant: log P scales like h^2
def neighbours(shape, hf):
    out = []
    for i in range(len(shape)):
        for d in (-1, 1):
            v = shape[i] + d
            if -1 <= v <= 6:
                s2 = list(shape); s2[i] = v; out.append((tuple(s2), hf))
    for hf2 in (hf - 0.1, hf + 0.1):
        if 0.5 <= hf2 <= 1.0: out.append((shape, round(hf2, 2)))
    return out
if __name__ == "__main__":
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    pool = mp.Pool(4); seen = {}
    starts = [((4,) + (-1,)*7, 1.0), ((3, 1) + (-1,)*6, 0.9), ((-1,)*8, 1.0),
              ((2, -1, 2, -1, 2, -1, 2, -1), 1.0), ((5, 5, -1, -1, -1, -1, 0, 0), 1.0)]
    starts += [(tuple(rng.choice([-1, -1, -1, 0, 1, 2, 3, 4]) for _ in range(B)), 1.0) for _ in range(5)]
    best_all = None
    for st in starts:
        cur = st; cv = pool.map(evaluate, [(cur[0], cur[1], K)])[0]
        if cv is None: continue
        seen[cur] = cv
        while True:
            nb = [n for n in neighbours(*cur) if n not in seen]
            vals = pool.map(evaluate, [(n[0], n[1], K) for n in nb])
            for n, v in zip(nb, vals): seen[n] = v
            cands = [(v, n) for n, v in zip(nb, vals) if v is not None]
            if not cands or min(cands)[0] >= cv: break
            cv, cur = min(cands)
        print(f"start {st} -> local min {cur} score {cv:+.4f}", flush=True)
        if best_all is None or cv < best_all[0]: best_all = (cv, cur)
    top = sorted((v, k) for k, v in seen.items() if v is not None)[:8]
    print("evaluated", len(seen), "shapes; best:")
    for v, k in top: print(f"  {v:+.4f}  {k}")
    json.dump([[v, list(k[0]), k[1]] for v, k in top], open(f"search_s{S}_top.json", "w"))
