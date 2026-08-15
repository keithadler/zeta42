"""Q(zeta_n, sqrt(-D)): an IMAGINARY quadratic extension of a cyclotomic field.

Why this exists: the Moser spindle's chromatic twist is rotation by
cos t = 5/6, sin t = sqrt(11)/6, i.e. multiplication by

    w = 5/6 + i*sqrt(11)/6 = 5/6 + sqrt(-11)/6.

In fields without i (any Q(zeta_odd)), i*sqrt(11) is not a real-extension
element -- but sqrt(-11) is a clean imaginary quadratic generator.  Elements
here are pairs (a, b) over CycloField(n) meaning a + b*sqrt(-D):

    (a + b√-D)(c + d√-D) = (ac - D·bd) + (ad + bc)√-D
    conj(a + b√-D)       = conj(a) - conj(b)√-D        <- the sign flip
    |z|^2 = z · conj(z)  (real, as always)

Interface mirrors CycloExtField so udg.UDG works unchanged.
"""

from functools import lru_cache
from math import sqrt as _fsqrt

from cyclo import get_cyclofield


class ImagExtField:
    __slots__ = ("base", "D", "n", "zero", "one", "_sqD")

    def __init__(self, n, D):
        assert D > 0, "pass D for sqrt(-D)"
        self.base = get_cyclofield(n)
        self.n = n
        self.D = D
        b = self.base
        self.zero = (b.zero, b.zero)
        self.one = (b.one, b.zero)
        self._sqD = _fsqrt(D)

    # ------------------------------------------------------- constructors --
    def lift(self, a):
        return (a, self.base.zero)

    def rat(self, num, den=1):
        return (self.base.rat(num, den), self.base.zero)

    def zeta(self, k=1):
        return (self.base.zeta(k), self.base.zero)

    def sqrtnegD(self, num=1, den=1):
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
        ac = b.mul(x[0], y[0])
        bd = b.mul(x[1], y[1])
        return (b.sub(ac, b.scale_rat(bd, self.D) if hasattr(b, "scale_rat")
                else b.mul(bd, b.rat(self.D))),
                b.add(b.mul(x[0], y[1]), b.mul(x[1], y[0])))

    def conj(self, x):
        b = self.base
        return (b.conj(x[0]), b.neg(b.conj(x[1])))

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
        if centre is None:
            return self.mul(p, w)
        return self.add(self.mul(self.sub(p, centre), w), centre)

    def to_complex(self, x):
        b = self.base
        return b.to_complex(x[0]) + b.to_complex(x[1]) * (1j * self._sqD)

    def pfloat(self, p):
        z = self.to_complex(p)
        return (z.real, z.imag)

    def __repr__(self):
        return f"ImagExtField(zeta_{self.n}, sqrt(-{self.D}))"


@lru_cache(maxsize=None)
def get_imagext(n, D):
    return ImagExtField(n, D)
