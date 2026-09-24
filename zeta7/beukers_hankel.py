"""Idea C: Hankel determinants of the Beukers-type measure on (0,1):
  int_0^1 z^k (-log z)^{s-1}/(1-z) dz = (s-1)! (zeta(s) - H_k^{(s)}).
Positive measure, so Delta(zeta s) > 0; X enters every moment with the same coefficient, so B has
rank 1 and Delta(X) = q X - p is a linear form. Criterion: log |Delta| - log content -> -oo."""
import math, sys
from math import gcd
from flint import fmpq, fmpq_mat, arb, ctx
def run(s, h):
    H = [fmpq(0)]
    for k in range(1, 2*h): H.append(H[-1] + fmpq(1, k**s))
    A = fmpq_mat(h, h, [-H[i+j] for i in range(h) for j in range(h)])
    ones = fmpq_mat(h, h, [1]*(h*h))
    d0 = A.det(); d1 = (A + ones).det()
    q = d1 - d0; p = -d0                      # Delta(X) = q X - p  (factor (s-1)!^h dropped)
    g = gcd(gcd(int(q.p), int(p.p)), 0); l = math.lcm(int(q.q), int(p.q))
    Q = int(q.p)*(l//int(q.q))//g; P = int(p.p)*(l//int(p.q))//g
    ctx.prec = 2*max(abs(Q).bit_length(), abs(P).bit_length()) + 2000
    v = arb(Q)*arb(s).zeta() - arb(P)
    assert v > 0
    return math.log(abs(Q)), float(v.log())
for s in (3, 5, 7):
    for h in (4, 8, 12, 16, 24, 32):
        lq, lv = run(s, h)
        print(f"s={s} h={h:2d}: log q = {lq:9.2f}   log(qX - p) = {lv:+10.2f}   /h^2 {lv/h**2:+.4f}", flush=True)
