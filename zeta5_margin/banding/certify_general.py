"""Rigorous real-side certificate for a general band profile (generalises ../certify.py).

Proves   2 U^rho(t) - V(t) < M0   for all t >= 0, where
    V(t) = 2 pi sqrt t - sum_b eps_b int_{u0_b}^{u1_b} log(t + u^2) du,
for a nested rational arcsine measure rho of mass lam supported in (0, 2), and returns
    R_cert = lam*M0 - I(rho)   >=   sup_sigma [ I(sigma) - int V dsigma ],
the real constant in  log Delta_K(zeta5) <= kappa K^2 log K + R K^2 + o(K^2)  (same argument as
the working edition's Section 6 with its field V; the scalar S_K and C_* are not used here).

Cells:  sup U^rho on [l, r] <= max(U(l), U(r)) (nesting => monotone between endpoints);
        inf V on [l, r] >= V_inc(l) + V_dec(r), splitting V into the increasing part
        2 pi sqrt t - sum_{eps<0} eps int log  and the decreasing part  - sum_{eps>0} eps int log.
Tail t >= 2:  U^rho(t) <= lam log t and  sum eps int log(t+u^2) <= E+ log(t+1) - E- log t,
        so 2U - V <= (2 lam - E-) log t + E+ log(t+1) - 2 pi sqrt t, whose derivative is < 0 when
        (2 lam + E+) <= pi sqrt t; both that and the value at t = 2 are checked.
Arithmetic: python-flint arb, 200 bits.  All data rational."""
import json, sys
from fractions import Fraction as F
from flint import arb, ctx
ctx.prec = 200
A = lambda x: arb(x.numerator) / arb(x.denominator)
PI = arb.pi()

def Fanti(t, u):          # int_0^u log(t + v^2) dv, exact for rational t >= 0, u > 0
    if u == 0: return arb(0)
    if t == 0: return 2*A(u)*A(u).log() - 2*A(u)
    s = A(t).sqrt()
    return A(u)*(A(t) + A(u)**2).log() - 2*A(u) + 2*s*(A(u)/s).atan()

def check(prof, lam, rows, M0, log=print):
    assert sum(c for a, b, c in rows) == lam, "mass"
    n = len(rows)
    assert all(0 < a < b < 2 and c > 0 and b - a > F(1, 225) for a, b, c in rows), "lengths"
    assert all(rows[j+1][0] < rows[j][0] < rows[j][1] < rows[j+1][1] for j in range(n-1)), "nesting"
    caps = [(A(b - a)/4).log() for a, b, c in rows]
    cache = {}
    def U(t):
        if t not in cache:
            s = arb(0)
            for (a, b, c), cap in zip(rows, caps):
                s += A(c)*(cap if a <= t <= b else ((A(abs(t - (a+b)/2)) + A((t-a)*(t-b)).sqrt())/2).log())
            cache[t] = s
        return cache[t]
    def part(t, sign):
        v = 2*PI*A(t).sqrt() if sign > 0 else arb(0)
        for u0, u1, e in prof:
            if (e < 0) == (sign > 0) and e != 0:
                v -= A(e)*(Fanti(t, u1) - Fanti(t, u0))
        return v
    Vc = {}
    def Vp(t, s):
        k = (t, s)
        if k not in Vc: Vc[k] = part(t, s)
        return Vc[k]
    thr = A(M0)
    ends = sorted(set([F(0), F(2)] + [a for a, b, c in rows] + [b for a, b, c in rows]))
    stack = list(zip(ends, ends[1:])); ncell = 0; worst = None
    while stack:
        l, r = stack.pop()
        up = max(U(l).upper(), U(r).upper())
        bound = 2*arb(up) - (Vp(l, +1) + Vp(r, -1))
        if bound < thr:
            ncell += 1; u = bound.upper(); worst = u if worst is None or u > worst else worst
            continue
        if r - l < F(1, 10**15): raise AssertionError(f"cell [{float(l)}, {float(r)}] fails")
        m = (l + r)/2; stack += [(l, m), (m, r)]
    Ep = sum((u1 - u0)*e for u0, u1, e in prof if e > 0)
    Em = -sum((u1 - u0)*e for u0, u1, e in prof if e < 0)
    assert A(2*lam + Ep) <= PI*arb(2).sqrt(), "tail monotonicity"
    tail2 = A(2*lam - Em)*arb(2).log() + A(Ep)*arb(3).log() - 2*PI*arb(2).sqrt()
    assert tail2 < thr, "tail at t = 2"
    S = F(0); I = arb(0)
    for (a, b, c), cap in zip(rows, caps):
        old = S; S += c; I += A(S*S - old*old)*cap
    Rc = A(lam)*thr - I
    log(f"  {n} components, {ncell} cells PASS (worst 2U-V below M0 by {float((thr - arb(worst)).lower()):.1e}); "
        f"R_cert = lam*M0 - I(rho) in {Rc.str(10)}  (upper {float(Rc.upper()):.7f})")
    return Rc.upper()

if __name__ == "__main__":
    d = json.load(open(sys.argv[1])); den = d["denominator"]
    rows = [tuple(F(v, den) for v in r) for r in d["rows"]]
    prof = [(F(*u0), F(*u1), F(*e)) for u0, u1, e in d["profile"]]
    check(prof, F(*d["lam"]), rows, F(*d["M0"]))
