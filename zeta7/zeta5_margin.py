"""zeta(5) margin: the paper's construction (K=40n, N=3n, D_N^6/D_K, i.e. net exponent 5 on j<=N)
with h = alpha*(K-N) rows instead of K-N, and an optional numerator factor t^m.
Reports log P_K(zeta5) (primitive integer polynomial) and log P/n^2; the paper's own values are
-265 (n=1), -833 (n=2), -1804 (n=3)."""
import sys, time, multiprocessing as mp
from manyrows import score
def run(args):
    n, alpha, m = args
    K, N = 40*n, 3*n
    e = {j: 5 for j in range(1, N+1)}
    for j in range(N+1, K+1): e[j] = -1
    if m: e[0] = m
    h = round(alpha*(K-N))
    t0 = time.time(); deg, lP = score(5, e, h)
    return n, alpha, m, h, deg, lP, time.time()-t0
if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1].split(',')]
    alphas = [float(x) for x in sys.argv[2].split(',')]
    ms = [int(x) for x in sys.argv[3].split(',')]
    jobs = [(n, a, m) for n in ns for a in alphas for m in ms]
    with mp.Pool(int(sys.argv[4]) if len(sys.argv) > 4 else 4) as pool:
        for n, a, m, h, deg, lP, dt in pool.imap(run, jobs):
            print(f"n={n} K={40*n} alpha={a:.2f} h={h:3d} t^{m}: deg {deg:3d}  log P(zeta5) = {lP:+10.2f}  /n^2 = {lP/n**2:+9.2f}  ({dt:.0f}s)", flush=True)
