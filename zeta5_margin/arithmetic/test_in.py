"""Check the extended inner bound gamma_p^in (all hypotheses of the inner proposition enforced)
against the true valuation v_p^G(Delta_K) at K = 40, 80, 120, including the new window
2K/5 <= p < K/3 ... i.e. 5/2 <= K/p < 3.  Expect: violations: 0."""
import json
from slack import vS
from gamma_in import gamma_in
bad = 0
import os
from slack import run
for K in (40, 80, 120):
    N = 3*K//40; h = K - N
    if not os.path.exists(f"rows{K}.json"):
        _, rows = run(K//40); json.dump(rows, open(f"rows{K}.json", "w"))
    rows = json.load(open(f"rows{K}.json"))
    for p, true, _ in rows:
        if p < 7 or p > K/2.3: continue
        vG = true - vS(K, N, h, p)
        best = None
        for L0 in range(0, 12):
            g = gamma_in(K, N, p, L0)
            if g is None: continue
            if g > vG: bad += 1; print(f"  VIOLATION K={K} p={p} L0={L0}: gamma_in {g} > true {vG}")
            best = g if best is None or g > best else best
        tag = "paper range x>=3" if K/p >= 3 else "NEW window" if K/p > 2.465 else "beyond window"
        print(f"K={K:3d} p={p:3d} x={K/p:5.2f} [{tag:16s}]  true vG={vG:5d}  best gamma_in={str(best):>6s}")
print("violations:", bad)
