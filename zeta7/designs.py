from math import comb
from residues import score_res, R_nonneg
K = 20
base = lambda j: (-1)**(j-1) * comb(2*K, K-j) * j*j
designs = {
  "product 1/D_K (f=1)":            lambda j: 1,
  "f = C(2K,K-j)":                  lambda j: comb(2*K, K-j),
  "f = C(K+j,2j)":                  lambda j: comb(K+j, 2*j),
  "f = C(K+j,j)^2 (Apery-like)":    lambda j: comb(K+j, j)**2,
  "f = C(K+j,2j)^2":                lambda j: comb(K+j, 2*j)**2,
  "f = C(K,j)^2":                   lambda j: comb(K, j)**2,
  "f = j^2":                        lambda j: j*j,
  "f = j^4":                        lambda j: j**4,
  "f = C(K+j,K-j)*C(K,j)":          lambda j: comb(K+j, K-j)*comb(K, j),
  "f = reflect C(K-1,j-1)^2":       lambda j: comb(K-1, j-1)**2,
  "f = C(2j,j)":                    lambda j: comb(2*j, j),
}
for name, f in designs.items():
    c = {j: base(j) * f(j) for j in range(1, K+1)}
    c = {j: v for j, v in c.items() if v != 0}
    ok = R_nonneg(c)
    if not ok:
        print(f"{name:32s} R changes sign -> not a positive measure, skipped", flush=True); continue
    out = []
    for s in (5, 7):
        best = min((score_res(s, c, h)[1]/K**2, h) for h in (K, K+2, K+4))
        out.append(f"s={s}: {best[0]:+.4f} (h={best[1]})")
    print(f"{name:32s} " + "   ".join(out), flush=True)
