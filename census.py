"""Cross-world coincidence census.

Question: glue the heptagonal arc lattice (Haugland's world, Q(zeta_42)) and
the classical multiquadratic world (Parts' 510, Q(sqrt3,sqrt5,sqrt11)) at the
natural anchor -- both contain 0 and (1,0) -- and count unit distances BETWEEN
the worlds.  Rich interface -> a genuinely new place to hunt primitives.
Near-empty interface -> the mixed-field direction dies this afternoon.

Method: float64 prescreen of all cross pairs, then 80-digit mpmath
verification of every candidate.  80 digits is not a formal proof (the two
fields share sqrt3, so the honest exact test needs the degree-96 compositum,
built only if the census warrants it), but at heights this small a false
80-digit coincidence is not a realistic concern for a strategic count.

A cross pair is TRIVIAL if both endpoints are shared points (present in both
worlds); those just re-detect the common triangular skeleton.
"""

import time

import numpy as np
from mpmath import mp, mpc, mpf, sqrt as msqrt, exp as mexp, pi as mpi

mp.dps = 80

from hept_g1 import arc_lattice, balls, SCALE
from parts510 import load


def world_A(radius):
    """Arc-lattice points as 80-digit complex numbers (and the raw vectors)."""
    fld, arcs, steps = arc_lattice()
    ball = balls(steps, radius)
    zeta = mexp(2j * mpi / 42)
    zpow = [zeta ** k for k in range(12)]
    pts = []
    for vec in ball:
        z = mpc(0)
        for k, c in enumerate(vec):
            if c:
                z += mpf(c) / SCALE * zpow[k]
        pts.append(z)
    return list(ball.keys()), pts


def world_B():
    """The 510's exact points at 80 digits."""
    fld, g, _ = load("510")
    s3, s5, s11 = msqrt(3), msqrt(5), msqrt(11)
    basis = [mpf(1), s3, s5, s3 * s5, s11, s3 * s11, s5 * s11, s3 * s5 * s11]

    def ev(coord):
        return sum(mpf(c.numerator) / mpf(c.denominator) * b
                   for c, b in zip(coord, basis) if c)

    return [(ev(p[0]), ev(p[1])) for p in g.points]


def census(radius=2):
    t0 = time.time()
    avecs, A = world_A(radius)
    B = world_B()
    print(f"world A (arc ball r<={radius}): {len(A)} points")
    print(f"world B (Parts 510):           {len(B)} points")

    Ax = np.array([float(z.real) for z in A])
    Ay = np.array([float(z.imag) for z in A])
    Bx = np.array([float(x) for x, y in B])
    By = np.array([float(y) for x, y in B])

    # shared points: identical in both worlds (float then 80-digit)
    shared_A = set()
    shared_B = set()
    for j in range(len(B)):
        d2 = (Ax - Bx[j]) ** 2 + (Ay - By[j]) ** 2
        for i in np.nonzero(d2 < 1e-12)[0]:
            if abs(A[i] - mpc(B[j][0], B[j][1])) < mpf(10) ** -60:
                shared_A.add(int(i))
                shared_B.add(j)
    print(f"shared points (in both worlds): {len(shared_B)}")

    # cross unit pairs
    cand = []
    chunk = 20000
    for lo in range(0, len(A), chunk):
        hi = min(lo + chunk, len(A))
        dx = Ax[lo:hi, None] - Bx[None, :]
        dy = Ay[lo:hi, None] - By[None, :]
        d2 = dx * dx + dy * dy
        ii, jj = np.nonzero(np.abs(d2 - 1.0) < 1e-9)
        cand.extend(zip((ii + lo).tolist(), jj.tolist()))
    print(f"float candidates: {len(cand)}")

    hits, both_shared, one_shared = [], 0, 0
    one = mpf(1)
    tol = mpf(10) ** -60
    for i, j in cand:
        d = A[i] - mpc(B[j][0], B[j][1])
        if abs(abs(d) - one) < tol:
            ia, jb = i in shared_A, j in shared_B
            if ia and jb:
                both_shared += 1        # the common triangular skeleton
            elif ia or jb:
                one_shared += 1         # an existing edge of one world
            else:
                hits.append((i, j))     # A-only point <-> B-only point: NEW
    print(f"verified cross unit pairs: {both_shared + one_shared + len(hits)}")
    print(f"  both endpoints shared (skeleton):     {both_shared}")
    print(f"  one endpoint shared (existing edges): {one_shared}")
    print(f"  NEITHER shared — genuinely new:       {len(hits)}")

    if hits:
        print("\nsample non-trivial coincidences (A-vector -> B-index):")
        for i, j in hits[:10]:
            print(f"  A{avecs[i]} ~ B[{j}]  at ({float(A[i].real):.6f},"
                  f"{float(A[i].imag):.6f})-({float(B[j][0]):.6f},"
                  f"{float(B[j][1]):.6f})")
        # how many distinct A and B points participate?
        pa = len({i for i, _ in hits})
        pb = len({j for _, j in hits})
        print(f"distinct participants: {pa} A-points, {pb} B-points")
    print(f"[{time.time()-t0:.0f}s]")
    return hits


if __name__ == "__main__":
    import sys
    r = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    census(r)
