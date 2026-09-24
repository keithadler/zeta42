# Weight w_s(y) = y^s F^{(s-1)}(y) / c_s,  F = 1/(e^{2pi y}-1),  c_s = (s-1)!/2
# Find r_s(j) := int_0^oo w_s(y)/(y^2+j^2) dy - j^{s-1}(zeta(s) - H_j^{(s)})
from mpmath import mp, mpf, quad, polylog, exp, pi, zeta, inf, factorial, identify, nsum
mp.dps = 40
def Fd(y, m):  # m-th derivative of F, m even: (2pi)^m Li_{-m}(e^{-2pi y})
    return (2*pi)**m * polylog(-m, exp(-2*pi*y))
for s in (3, 5, 7):
    c = factorial(s-1)/2
    for j in range(1, 7):
        I = quad(lambda y: y**s*Fd(y, s-1)/(y*y+j*j)/c, [0, 0.05, 0.5, 2, 10, inf])
        H = sum(mpf(1)/k**s for k in range(1, j+1))
        r = I - j**(s-1)*(zeta(s)-H)
        print(s, j, mp.nstr(r, 25), identify(r))
