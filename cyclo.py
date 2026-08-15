"""Exact arithmetic in a cyclotomic field Q(zeta_n).

Why this exists: `field.py` handles multiquadratic fields Q(√p1,...,√pk), and
every known 5-chromatic unit-distance construction -- de Grey, Heule, Parts --
lives in one.  But 7-fold symmetry does not:

    degree of cos(2*pi/n) over Q = phi(n)/2
    n = 3,4,6 -> 1     n = 5,8,10,12 -> 2     n = 7 -> 3

cos(2*pi/7) is a root of 8x^3 + 4x^2 - 4x - 1, so a regular heptagon cannot be
written with square roots at all.  Haugland's 2026 spindle-free construction
needs one.  This module opens that arithmetic.

A point here is a single field element interpreted as a complex number, not a
coordinate pair.  That makes rotation a multiplication and squared distance a
multiplication by a conjugate:

    |z|^2 = z * conj(z),   conj(zeta) = zeta^-1

Elements are tuples of phi(n) Fractions against the power basis
zeta^0 ... zeta^(phi(n)-1), reduced modulo the nth cyclotomic polynomial.
They are canonical, so equality is tuple equality and they hash.

The interface deliberately mirrors `field.Field` (zero, one, origin, padd,
psub, is_unit, pfloat), so `udg.UDG` works over either without changes -- it
only ever asks whether two points are exactly one apart.
"""

import cmath
from fractions import Fraction
from functools import lru_cache

_Q0 = Fraction(0)


# ------------------------------------------------- cyclotomic polynomial ----

def _poly_divmod(num, den):
    """Divide integer polynomials (low-order-first lists).  Exact division."""
    num = list(num)
    out = [0] * (len(num) - len(den) + 1)
    for i in range(len(out) - 1, -1, -1):
        c = num[i + len(den) - 1] // den[-1]
        out[i] = c
        if c:
            for j, d in enumerate(den):
                num[i + j] -= c * d
    return out


@lru_cache(maxsize=None)
def cyclotomic_poly(n):
    """Phi_n(x) as a low-order-first list of ints.  Phi_n = (x^n-1)/prod Phi_d."""
    num = [-1] + [0] * (n - 1) + [1]              # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            num = _poly_divmod(num, cyclotomic_poly(d))
    return num


class CycloField:
    """Q(zeta_n), with elements as complex numbers in the power basis."""

    __slots__ = ("n", "dim", "_phi", "_pow", "zero", "one", "_zc")

    def __init__(self, n):
        self.n = n
        phi = cyclotomic_poly(n)
        self.dim = len(phi) - 1
        self._phi = phi
        d = self.dim

        # zeta^k in the power basis, for k up to what multiplication needs
        top = max(2 * d, n + 1)
        pw = [[_Q0] * d for _ in range(top)]
        for k in range(min(d, top)):
            pw[k][k] = Fraction(1)
        for k in range(d, top):
            # zeta^k = zeta * zeta^(k-1), then fold zeta^d back in via Phi_n
            prev = pw[k - 1]
            shifted = [_Q0] + prev[:-1]
            carry = prev[-1]
            if carry:
                # zeta^d = -(c_0 + c_1 zeta + ... + c_{d-1} zeta^{d-1})
                for i in range(d):
                    shifted[i] -= carry * phi[i]
            pw[k] = shifted
        self._pow = [tuple(r) for r in pw]

        self.zero = (_Q0,) * d
        self.one = self.rat(1)
        self._zc = cmath.exp(2j * cmath.pi / n)

    # ------------------------------------------------------- constructors --
    def rat(self, num, den=1):
        v = [_Q0] * self.dim
        v[0] = Fraction(num, den)
        return tuple(v)

    def zeta(self, k=1):
        """zeta_n^k as a field element."""
        return self._pow[k % self.n]

    # ---------------------------------------------------------- arithmetic --
    @staticmethod
    def add(a, b):
        return tuple(x + y for x, y in zip(a, b))

    @staticmethod
    def sub(a, b):
        return tuple(x - y for x, y in zip(a, b))

    @staticmethod
    def neg(a):
        return tuple(-x for x in a)

    def mul(self, a, b):
        d = self.dim
        raw = [_Q0] * (2 * d - 1)
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if y:
                    raw[i + j] += x * y
        out = [_Q0] * d
        pw = self._pow
        for k, c in enumerate(raw):
            if not c:
                continue
            if k < d:
                out[k] += c
            else:
                for i, e in enumerate(pw[k]):
                    if e:
                        out[i] += c * e
        return tuple(out)

    def conj(self, a):
        """Complex conjugate: zeta -> zeta^-1."""
        out = [_Q0] * self.dim
        pw = self._pow
        for i, c in enumerate(a):
            if not c:
                continue
            for j, e in enumerate(pw[(self.n - i) % self.n]):
                if e:
                    out[j] += c * e
        return tuple(out)

    def inv(self, a):
        """Multiplicative inverse, via the rational norm a * conj-product."""
        # brute force over the basis: solve the dim x dim rational system
        d = self.dim
        cols = [self.mul(a, tuple(Fraction(1) if k == i else _Q0 for k in range(d)))
                for i in range(d)]
        # gaussian elimination on the matrix whose columns are cols
        M = [[cols[c][r] for c in range(d)] + [Fraction(1) if r == 0 else _Q0]
             for r in range(d)]
        for col in range(d):
            piv = next((r for r in range(col, d) if M[r][col]), None)
            if piv is None:
                raise ZeroDivisionError("not invertible")
            M[col], M[piv] = M[piv], M[col]
            f = M[col][col]
            M[col] = [x / f for x in M[col]]
            for r in range(d):
                if r != col and M[r][col]:
                    g = M[r][col]
                    M[r] = [x - g * y for x, y in zip(M[r], M[col])]
        return tuple(M[r][d] for r in range(d))

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    # -------------------------------------------------------------- points --
    @property
    def origin(self):
        return self.zero

    def padd(self, p, q):
        return self.add(p, q)

    def psub(self, p, q):
        return self.sub(p, q)

    def norm2(self, p):
        """|p|^2 as a field element (in fact rational)."""
        return self.mul(p, self.conj(p))

    def dist2(self, p, q):
        return self.norm2(self.sub(p, q))

    def is_unit(self, p, q):
        """Exact |p - q| == 1."""
        return self.dist2(p, q) == self.one

    def rotate(self, p, k, centre=None):
        """Rotate about `centre` by 2*pi*k/n -- multiplication by zeta^k."""
        if centre is None:
            return self.mul(p, self.zeta(k))
        return self.add(self.mul(self.sub(p, centre), self.zeta(k)), centre)

    def to_complex(self, a):
        return sum(float(c) * self._zc ** i for i, c in enumerate(a) if c)

    def pfloat(self, p):
        z = self.to_complex(p)
        return (z.real, z.imag)

    def __repr__(self):
        return f"CycloField(zeta_{self.n}, degree {self.dim})"


@lru_cache(maxsize=None)
def get_cyclofield(n):
    return CycloField(n)
