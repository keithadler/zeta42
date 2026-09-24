"""Scale a banded zeta(5) profile from n=1 to size n: band i (j = i at n=1) becomes j in ((i-1)n, i n].
Poles for j > B n up to K = 40n; rows h = alpha * #poles. Prints log P(zeta5)."""
import sys, time, multiprocessing as mp
from manyrows import score
def run(arg):
    n, ex, a = arg
    K = 40*n; e = {}
    for i, v in enumerate(ex, start=1):
        for j in range((i-1)*n+1, i*n+1):
            if v: e[j] = v
    for j in range(len(ex)*n+1, K+1): e[j] = -1
    m = sum(1 for v in e.values() if v == -1)
    t = time.time(); d, lP = score(5, e, round(a*m))
    return n, ex, a, round(a*m), m, lP, time.time()-t
if __name__ == "__main__":
    n = int(sys.argv[1])
    cfgs = [((4, 4, 3, 2, 0, 0), 1.2), ((5, 4, 2, 2, 0, 0), 1.2), ((5, 5, 3, 1), 1.1)]
    with mp.Pool(3) as pool:
        for n_, ex, a, h, m, lP, dt in pool.imap(run, [(n, ex, a) for ex, a in cfgs]):
            print(f"n={n_} bands={ex} alpha={a} poles={m} h={h}: log P(zeta5) = {lP:+.2f}  /n^2 = {lP/n_**2:+.2f}  ({dt:.0f}s)", flush=True)
