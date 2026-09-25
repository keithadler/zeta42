"""Rigorous check of a nested arcsine comparison measure for Fauzan's zeta(5) proof.

Same argument as certificates/verify.py (long-mathematics/zeta5-irrationality, verify_A):
  * nested supports  =>  U^rho is monotone between consecutive endpoints a_j, b_j,
    so on a cell [l,r] inside one such interval, sup U^rho <= max(U(l), U(r));
  * V is decreasing then increasing (one minimum q in [q_-,q_+]), so inf V on a cell
    is V(r), V(l) or V(q) as in the paper;
  * cells are bisected adaptively until 2*sup U - inf V < M0 holds rigorously.
Arithmetic: python-flint arb balls (rigorous, outward-rounded).  All data rational.
Also: t >= 2 is covered by the paper's tail argument, which needs supp in (0,2) and
M0 above the field-tail bound at t=2 (checked here too).
"""
import sys, json
from fractions import Fraction as F
from flint import arb, ctx
ctx.prec = 200
alpha = F(3, 40); lam = F(37, 40)
A = lambda x: arb(x.numerator) / arb(x.denominator)
PI = arb.pi()

def V(t):
    if t == 0: return -12*A(alpha)*A(alpha).log() - 2 + 12*A(alpha)
    r = A(t).sqrt()
    return ((1 + A(t)).log() - 6*A(alpha)*(A(t) + A(alpha)**2).log() - 2 + 12*A(alpha)
            + 2*r*(PI + (1/r).atan() - 6*(A(alpha)/r).atan()))

def check(rows, M0, max_cells=10**7, log=print):
    n = len(rows)
    assert sum(c for a, b, c in rows) == lam, "mass"
    assert all(0 < a < b < 2 and c > 0 and b - a > F(1, 225) for a, b, c in rows), "lengths"
    assert all(rows[j+1][0] < rows[j][0] < rows[j][1] < rows[j+1][1] for j in range(n-1)), "nesting"
    caps = [(A(b - a)/4).log() for a, b, c in rows]
    cache = {}
    def U(t):
        if t in cache: return cache[t]
        s = arb(0)
        for (a, b, c), cap in zip(rows, caps):
            if a <= t <= b: u = cap
            else: u = ((A(abs(t - (a+b)/2)) + A((t-a)*(t-b)).sqrt())/2).log()
            s += A(c)*u
        cache[t] = s; return s
    # minimiser bracket of V (paper's values), and V* lower bound
    qm, qp = F(59205077, 10**10), F(59205079, 10**10)
    dphi = lambda t: 2*(PI + (1/A(t).sqrt()).atan() - 6*(A(alpha)/A(t).sqrt()).atan())
    assert dphi(qm) < 0 and dphi(qp) > 0, "minimiser bracket"
    par = PI + (1/A(qp).sqrt()).atan() - 6*(A(alpha)/A(qm).sqrt()).atan()
    assert par < 0
    Vstar = ((1+A(qm)).log() - 6*A(alpha)*(A(qp)+A(alpha)**2).log() - 2 + 12*A(alpha)
             + 2*A(qp).sqrt()*par)
    ends = sorted(set([F(0), F(2), qm, qp] + [a for a, b, c in rows] + [b for a, b, c in rows]))
    thr = A(M0)
    stack = [(l, r) for l, r in zip(ends, ends[1:])]
    ncell = 0; worst = None
    while stack:
        l, r = stack.pop()
        # sup U on the cell <= the larger upper endpoint of the two end values
        up = max(U(l).upper(), U(r).upper())
        lowV = V(r) if r <= qm else V(l) if l >= qp else Vstar
        bound = 2*arb(up) - lowV
        if bound < thr:
            ncell += 1
            u = bound.upper()
            if worst is None or u > worst: worst = u
            continue
        if r - l < F(1, 10**15) or ncell > max_cells:
            raise AssertionError(f"cell [{float(l)},{float(r)}] fails: {bound} vs {M0}")
        m = (l + r)/2
        stack += [(l, m), (m, r)]
    # t >= 2: paper Lemma (confinement) needs 13/10 log t + 2 alpha^3/t - 2 pi sqrt t + sqrt2/K < M0 at t = 2
    tail2 = arb(13)/10*arb(2).log() + 2*A(alpha)**3/2 - 2*PI*arb(2).sqrt()
    assert tail2 < thr, "tail at t=2"
    # energy (nested closed form) and the real constant
    S = F(0); I = arb(0)
    for (a, b, c), cap in zip(rows, caps):
        old = S; S += c; I += A(S*S - old*old)*cap
    Cstar = (-2*A(lam) + 12*A(alpha)*A(lam)*(1 - A(alpha).log()) + 3*A(lam)**2
             - 2*A(lam)**2*(2*A(lam)).log())
    real = A(lam)*thr - I + Cstar
    gap = thr - arb(worst)
    log(f"  {n} components, {ncell} cells PASS; worst 2U-V upper = {float(worst):.12f}, "
        f"below M0 = {float(M0)} by {float(gap.lower()):.2e}")
    log(f"  U = lam*M0 - I(rho) + C*  in  {real.str(12)}   (upper {float(real.upper()):.9f})")
    return real.upper()

if __name__ == "__main__":
    raw = json.load(open(sys.argv[1])); den = raw["denominator"]
    rows = [tuple(F(v, den) for v in r) for r in raw["rows"]]
    M0 = F(*raw["M0"]) if "M0" in raw else F(-1329, 200)
    check(rows, M0)
