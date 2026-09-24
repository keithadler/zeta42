from sweep import parts, show
from hankel_s import Dpoly
K = 40
for N1, r1, N2, r2 in [(3,4,8,2),(3,6,8,1),(4,4,10,1),(2,6,6,2),(4,3,8,2),(5,4,10,1),(3,5,12,1),(5,3,7,2),(2,4,5,2),(5,5,6,1)]:
    W = Dpoly(N1)**r1 * Dpoly(N2)**r2
    show(f"s=7 K={K} D_{N1}^{r1}*D_{N2}^{r2}", K, parts(7, K, N2, 0, W))
