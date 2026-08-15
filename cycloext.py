"""Q(zeta_n, sqrt(m)): a quadratic extension of a cyclotomic field.

Needed because Haugland's G3 spindles a forced pair at distance sqrt(3):
the chord-1 rotation at radius sqrt(3) has cos = 5/6, sin = sqrt(11)/6, and
sqrt(11) is not in Q(zeta_42).  So the final graph of the reproduction lives
in Q(zeta_42, sqrt(11)), degree 24.

An element is a pair (a, b) of CycloField elements meaning a + b*sqrt(m),
with m a positive square-free integer (sqrt(m) real, so complex conjugation
acts componentwise):

    (a + b√m)(c + d√m) = (ac + m·bd) + (ad + bc)√m
    conj(a + b√m)      = conj(a) + conj(b)√m
    |z|^2              = z · conj(z)

The interface mirrors CycloField/Field (zero, one, origin, padd, psub,
is_unit, rotate-by-multiplication, pfloat) so udg.UDG works unchanged.
"""

from functools import lru_cache
from math import sqrt as _fsqrt

from cyclo import get_cyclofield


class CycloExtField:
    """Q(zeta_n, sqrt(m)) as pairs over CycloField(n)."""

    __slots__ = ("base", "m", "n", "zero", "one", "_sqm")

    def __init__(self, n, m):
        self.base = get_cyclofield(n)
        self.n = n
        self.m = m
        b = self.base
        self.zero = (b.zero, b.zero)
        self.one = (b.one, b.zero)
        self._sqm = _fsqrt(m)

    # ------------------------------------------------------- constructors --
    def lift(self, a):
        """A CycloField element -> extension element a + 0*sqrt(m)."""
        return (a, self.base.zero)

    def rat(self, num, den=1):
        return (self.base.rat(num, den), self.base.zero)

    def zeta(self, k=1):
        return (self.base.zeta(k), self.base.zero)

    def sqrtm(self, num=1, den=1):
        """(num/den) * sqrt(m)."""
        return (self.base.zero, self.base.rat(num, den))

    # ---------------------------------------------------------- arithmetic --
    def add(self, x, y):
        b = self.base
        return (b.add(x[0], y[0]), b.add(x[1], y[1]))

    def sub(self, x, y):
        b = self.base
        return (b.sub(x[0], y[0]), b.sub(x[1], y[1]))

    def neg(self, x):
        b = self.base
        return (b.neg(x[0]), b.neg(x[1]))

    def mul(self, x, y):
        b = self.base
        a, bb = x
        c, d = y
        re = b.add(b.mul(a, c), tuple(x * self.m for x in b.mul(bb, d)))
        im = b.add(b.mul(a, d), b.mul(bb, c))
        return (re, im)

    def conj(self, x):
        b = self.base
        return (b.conj(x[0]), b.conj(x[1]))

    # -------------------------------------------------------------- points --
    @property
    def origin(self):
        return self.zero

    def padd(self, p, q):
        return self.add(p, q)

    def psub(self, p, q):
        return self.sub(p, q)

    def norm2(self, p):
        return self.mul(p, self.conj(p))

    def dist2(self, p, q):
        return self.norm2(self.sub(p, q))

    def is_unit(self, p, q):
        return self.dist2(p, q) == self.one

    def rotate_mul(self, p, w, centre=None):
        """Rotate by the unit complex number w (an extension element)."""
        if centre is None:
            return self.mul(p, w)
        return self.add(self.mul(self.sub(p, centre), w), centre)

    def to_complex(self, x):
        b = self.base
        return b.to_complex(x[0]) + self._sqm * b.to_complex(x[1])

    def pfloat(self, p):
        z = self.to_complex(p)
        return (z.real, z.imag)

    def __repr__(self):
        return f"CycloExtField(zeta_{self.n}, sqrt{self.m}, degree {2 * self.base.dim})"


@lru_cache(maxsize=None)
def get_cycloext(n, m):
    return CycloExtField(n, m)
