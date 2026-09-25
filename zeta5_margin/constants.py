"""What a certified real constant U does to the theorem, in exact arithmetic.

The arithmetic side is untouched: A_* and A_M are the working edition's (Appendix B).
The final bound is  limsup K^-2 log Q_{K,M}(zeta5) <= A_M + U   (eq. finalasymptotic),
and the approximation corollary needs  eps < -(A_1e5 + U),  eps > lam/c,  rho*c < mu.

    python constants.py
"""
from fractions import Fraction as F

LAM = F(37, 40)
ASTAR = F(9928298118277006344769, 7535670527041937280000)
RHO = F(6933, 500)                       # the paper's relative-norm constant, unchanged


def A(M):
    return ASTAR + 7 * LAM / M - (F(2923, 240) - F(1, 4)) / M**2 + F(32, M**3)


CASES = {
    "paper (16 components)": (F(-2733991, 2000000), F(139, 5), F(7907, 100),
                              F(79, 1600), F(75, 4), 260),
    "measure64.json":        (F(-1383657, 1000000), F(54, 1), F(105, 1),
                              F(105, 1600), F(57, 4), 198),
    "measure256.json":       (F(-692469, 500000), F(113, 2), F(1077, 10),
                              F(1077, 16000), F(55, 4), 191),
}

if __name__ == "__main__":
    for name, (U, c200, c1e5, eps, c, mu) in CASES.items():
        m200, m1e5 = -1600 * (A(200) + U), -1600 * (A(100000) + U)
        assert c200 < m200 and c1e5 < m1e5, name
        assert eps < -(A(100000) + U) and eps - LAM / c > 0 and RHO * c - mu < 0, name
        print(f"{name:22s} U <= {float(U):.6f}")
        print(f"   Q_(40n,200)(zeta5)    < exp(-{float(c200):g} n^2)   (margin {float(m200):.3f})")
        print(f"   Q_(40n,100000)(zeta5) < exp(-{float(c1e5):g} n^2)   (margin {float(m1e5):.3f})")
        print(f"   |zeta(5) - a/b| > b^-{mu}   (eps = {eps}, c = {c})")
