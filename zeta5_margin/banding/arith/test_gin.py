"""Validity and tightness of the generalised inner bound against true valuations."""
import json, sys
from fractions import Fraction as F
sys.path.insert(0, "../../arithmetic")
from gin_general import classes, gamma, paper_alloc, optimise, start_alloc
from true_band import build
def cases():
    for tag, ex, alpha in (("paper", (5,5,5), 1.0), ("bands", (4,4,3,2,0,0), 1.2)):
        for n in (1, 2, 3):
            try: d = json.load(open(f"true_{'_'.join(map(str,ex))}_t0_a{alpha}_K{40*n}.json"))
            except FileNotFoundError: continue
            e, K, h = build(ex, 0, n, alpha)
            yield tag, e, K, h, d
viol = 0
for tag, e, K, h, d in cases():
    small = max(j for j, v in e.items() if v != -1)
    for p, t, vG in d["rows"]:
        if not (p*p > 5*K and p > small and p <= K/2): continue
        m, ell, b, mK = classes(e, K, p)
        best_paper = None
        for L0 in range(0, 12):
            if 5 + 4*L0 >= p + 1: continue
            L = paper_alloc(ell, b, h, L0) or start_alloc(ell, b, h, L0)
            g = gamma(L, ell, b, mK)
            if best_paper is None or g > best_paper[0]: best_paper = (g, L)
        if best_paper is None: continue
        Lo, go = optimise(best_paper[1], ell, b, mK)
        if 5 + 4*Lo[0] >= p + 1: go = None
        v = (best_paper[0] > vG) or (go is not None and go > vG); viol += v
        print(f"{tag:5s} K={K:3d} p={p:3d} x={K/p:5.2f}: true {vG:5d}  paper-alloc {str(best_paper[0]):>7s}  optimised {str(go):>7s}  {'VIOLATION' if v else ''}", flush=True)
print("violations:", viol)
