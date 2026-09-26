"""Constructions compared on the real side.  A profile is a list of (u0, u1, eps): the measure
carries prod_j (t + j^2)^eps_j, and with u = j/K the exponent density eps(u) is piecewise constant.
eps = -1 marks poles.  lam = h/K (rows per K).  A factor t^m in the numerator only contributes
O(K log K) and is ignored at order K^2."""
from fractions import Fraction as F

def bands(ex, rows_factor):
    """Banded profile as in zeta7/zeta5_deep.py: band i = ((i-1)/40, i/40] has exponent ex[i-1];
    poles on (len(ex)/40, 1]; h = rows_factor * #poles."""
    prof = [(F(i, 40), F(i + 1, 40), F(e)) for i, e in enumerate(ex) if e != 0]
    npole_frac = 1 - F(len(ex), 40)
    prof.append((F(len(ex), 40), F(1), F(-1)))
    return prof, F(rows_factor) * npole_frac

PROFILES = {
    "paper (5 on j<=3n)":             ([(F(0), F(3, 40), F(5)), (F(3, 40), F(1), F(-1))], F(37, 40)),
    "paper + 10% rows":               ([(F(0), F(3, 40), F(5)), (F(3, 40), F(1), F(-1))], F(11, 10) * F(37, 40)),
    "bands (5,5,3,1) x1.10":          bands((5, 5, 3, 1), F(11, 10)),
    "bands (4,4,3,2,0,0) x1.20":      bands((4, 4, 3, 2, 0, 0), F(6, 5)),
    "bands (5,4,2,1,0) t x1.15":      bands((5, 4, 2, 1, 0), F(23, 20)),
}
