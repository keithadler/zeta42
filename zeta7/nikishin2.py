"""Nikishin rows with f's poles inside the numerator-zero region (j <= N): f*R has no new poles.
Base: R = D_N^r / D_(N,K], i.e. W0 = D_N^r, poles N+1..K. f = sum_{j<=N} c_j/(t+j^2), c_j > 0.
Row types: t^i W0/D and t^i W0 f/D (and optionally W0 f^2 for a 3-level 'cube')."""
from flint import fmpz_poly
from hankel_s import Dpoly
from nikishin import score
def W_times_markov(W, js, cs):
    num = fmpz_poly([0]); F = fmpz_poly([1])
    for j in js: F *= fmpz_poly([j*j, 1])
    for j, c in zip(js, cs): num += c*(F // fmpz_poly([j*j, 1]))
    q, r = divmod(W*num, F)
    assert r == 0
    return q
for s, K, N, r in [(5, 20, 2, 4), (7, 20, 2, 4), (7, 20, 3, 4)]:
    poles = list(range(N+1, K+1)); W0 = Dpoly(N)**r; m = len(poles)
    for n in (m, m+2):
        d, lp = score(s, [(W0, poles, n)], n)
        print(f"s={s} K={K} W=D_{N}^{r}: control Gram N={n}: deg {d} logP/N^2 {lp/n**2:+.4f}", flush=True)
    js = list(range(1, N+1))
    for cs in ([1]*N, [j*j for j in js], [1]+[0]*(N-1) if N > 1 else [1]):
        W1 = W_times_markov(W0, js, cs)
        for split in ((m-4, 4), (m-2, 4), (m//2, m//2)):
            n = sum(split)
            d, lp = score(s, [(W0, poles, split[0]), (W1, poles, split[1])], n)
            print(f"s={s} K={K} W=D_{N}^{r}: f c={cs}, rows {split}: deg {d} logP/N^2 {lp/n**2:+.4f}", flush=True)
    # 3-level: rows t^i W0, t^i W0 f, t^i W0 f^2
    W1 = W_times_markov(W0, js, [1]*N); W2 = W_times_markov(W1, js, [1]*N)
    for split in ((m-6, 3, 3), (m//3 + 2, m//3, m - 2*(m//3) - 2)):
        n = sum(split)
        d, lp = score(s, [(W0, poles, split[0]), (W1, poles, split[1]), (W2, poles, split[2])], n)
        print(f"s={s} K={K} W=D_{N}^{r}: 3-level rows {split}: deg {d} logP/N^2 {lp/n**2:+.4f}", flush=True)
