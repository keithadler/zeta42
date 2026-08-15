"""Exact arithmetic in a real multiquadratic field Q(sqrt(p1), ..., sqrt(pk)).

Unit-distance graph work is unforgiving of floating point: two points that are
0.9999999997 apart are either an edge or not, and guessing wrong invalidates
everything downstream.  So every coordinate here is exact.

A field is fixed by a tuple of pairwise-coprime square-free integers.  An
element is a tuple of Fractions against the basis indexed by subsets:

    Field((3, 5, 7)) -> basis 1, √3, √5, √15, √7, √21, √35, √105

with basis[i] * basis[j] = (product of shared primes) * basis[i ^ j].
Elements are canonical tuples, so they compare and hash directly.

Fields you need for the classical constructions:
    (3, 11)          Moser spindle          cos t = 5/6,   sin t = √11/6
    (3, 5, 7)        de Grey's J -> K -> L  cos = 7/8, √15/8; 31/32, 3√7/32
    (3, 5, 7, 11)    both at once
"""

from fractions import Fraction
from functools import lru_cache
from math import sqrt

_Q0 = Fraction(0)


class Field:
    """Q(sqrt(p) for p in primes), as a vector space over Q."""

    __slots__ = ("primes", "dim", "_shared", "_sqrt", "zero", "one", "_names")

    def __init__(self, primes):
        primes = tuple(sorted(primes))
        self.primes = primes
        self.dim = 1 << len(primes)
        # basis[i] * basis[j] = _shared[i & j] * basis[i ^ j]
        self._shared = []
        self._names = []
        for m in range(self.dim):
            prod = 1
            for b, p in enumerate(primes):
                if m >> b & 1:
                    prod *= p
            self._shared.append(prod)
            self._names.append("1" if prod == 1 else f"√{prod}")
        self._sqrt = [sqrt(s) for s in self._shared]
        self.zero = (_Q0,) * self.dim
        self.one = self.rat(1)

    # ------------------------------------------------------- constructors --
    def rat(self, n, d=1):
        v = [_Q0] * self.dim
        v[0] = Fraction(n, d)
        return tuple(v)

    def surd(self, radicand, num=1, den=1):
        """(num/den) * sqrt(radicand); radicand must be a basis product."""
        try:
            m = self._shared.index(radicand)
        except ValueError:
            raise ValueError(
                f"√{radicand} is not in Q({', '.join('√%d' % p for p in self.primes)})"
            ) from None
        v = [_Q0] * self.dim
        v[m] = Fraction(num, den)
        return tuple(v)

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
        out = [_Q0] * self.dim
        shared = self._shared
        for i, x in enumerate(a):
            if not x:
                continue
            for j, y in enumerate(b):
                if not y:
                    continue
                out[i ^ j] += x * y * shared[i & j]
        return tuple(out)

    @staticmethod
    def scale(a, num, den=1):
        f = Fraction(num, den)
        return tuple(x * f for x in a)

    def to_float(self, a):
        return sum(float(c) * s for c, s in zip(a, self._sqrt) if c)

    def fmt(self, a):
        parts = []
        for c, n in zip(a, self._names):
            if not c:
                continue
            parts.append(str(c) if n == "1" else f"{c}·{n}")
        return " + ".join(parts) if parts else "0"

    # -------------------------------------------------------------- points --
    @property
    def origin(self):
        return (self.zero, self.zero)

    def padd(self, p, q):
        return (self.add(p[0], q[0]), self.add(p[1], q[1]))

    def psub(self, p, q):
        return (self.sub(p[0], q[0]), self.sub(p[1], q[1]))

    def dist2(self, p, q):
        d = self.psub(p, q)
        return self.add(self.mul(d[0], d[0]), self.mul(d[1], d[1]))

    def is_unit(self, p, q):
        """Exact |p - q| == 1.  The only edge test that is allowed to matter."""
        return self.dist2(p, q) == self.one

    def rotate(self, p, cos_t, sin_t, centre=None):
        if centre is None:
            centre = self.origin
        dx, dy = self.psub(p, centre)
        rx = self.sub(self.mul(cos_t, dx), self.mul(sin_t, dy))
        ry = self.add(self.mul(sin_t, dx), self.mul(cos_t, dy))
        return self.padd((rx, ry), centre)

    def pfloat(self, p):
        return (self.to_float(p[0]), self.to_float(p[1]))

    def __repr__(self):
        return f"Field({self.primes})"


@lru_cache(maxsize=None)
def get_field(primes):
    return Field(tuple(sorted(primes)))


# ------------------------------------------------------------- rotations ----

def chord_rotation(fld, radius2, chord2=1):
    """cos/sin of the rotation moving a point at |p|^2 = radius2 by |chord|.

    A rotation by t moves a point at radius r a distance 2*r*sin(t/2), so
    sin^2(t/2) = chord2 / (4 * radius2), giving
        cos t = 1 - chord2 / (2 * radius2)
        sin t = sqrt(1 - cos^2 t)
    Both must land in `fld` or this raises -- which is exactly the check you
    want before trusting a candidate rotation.
    """
    c = Fraction(1) - Fraction(chord2, 2 * radius2)
    s2 = 1 - c * c                      # a rational
    num, den = s2.numerator, s2.denominator
    # sqrt(num/den) = sqrt(num*den)/den; split off the square part
    square_part, rem = 1, num * den
    d = 2
    while d * d <= rem:
        while rem % (d * d) == 0:
            square_part *= d
            rem //= d * d
        d += 1
    cos_t = fld.rat(c.numerator, c.denominator)
    if rem == 1:
        sin_t = fld.rat(square_part, den)
    else:
        sin_t = fld.surd(rem, square_part, den)
    return cos_t, sin_t
