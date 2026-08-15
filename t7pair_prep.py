"""Export the T7 10-core DIRECTED-PAIR instance for kissat.

Question encoded: can the designated pair (0,0)-(2,0) share a colour in some
proper 5-colouring of the 18,959-vertex 10-core?

    SAT   -> pair still not forced (the mono-colouring is the witness;
             its solve time is the strain measurement)
    UNSAT -> the pair is FORCED non-mono at k=5: that is the object B,
             modulo the usual core-reduction triviality check.

Base colouring clauses plus colour-equality on the pair; no aggregate
constraints (that is t7core_dial.cnf's job).
"""

import time
from fractions import Fraction

import cmath
import numpy as np

from campaign_a import numpy_ball
from hept_g1 import arc_lattice, SCALE
from udg import UDG
from colour import ColourInstance
from rung1b import carve_T, kcore_ids


def log(*a):
    print(*a, flush=True)


def main():
    t0 = time.time()
    fld, arcs, steps = arc_lattice()
    zeta = cmath.exp(2j * cmath.pi / 42)
    basis = np.array([zeta ** k for k in range(12)])
    dist_map, dim = numpy_ball(steps, 4.0, basis, SCALE, 2_000_000,
                               log=lambda *a: None)
    log(f"ball: {len(dist_map)}  [{time.time()-t0:.0f}s]")

    B = np.zeros(dim, dtype=np.int16)
    B[0] = 2 * SCALE
    tofield = lambda vb: tuple(Fraction(int(c), SCALE)
                               for c in np.frombuffer(vb, dtype=np.int16))
    T7 = carve_T(dist_map, dim, B, 7)
    g7 = UDG(fld, [tofield(vb) for vb in T7])
    core = g7.subgraph(kcore_ids(g7, 10))
    log(f"10-core: {core.summary()}  [{time.time()-t0:.0f}s]")

    A_id = core.index_of(tofield(bytes(np.zeros(dim, dtype=np.int16))))
    B_id = core.index_of(tofield(bytes(B)))
    assert A_id is not None and B_id is not None
    log(f"designated pair: A={A_id} B={B_id}")

    inst = ColourInstance(core, 5, break_symmetry=False)
    for c in range(5):
        inst.add_clause([-inst.var(A_id, c), inst.var(B_id, c)])
        inst.add_clause([inst.var(A_id, c), -inst.var(B_id, c)])
    with open("t7pair_dial.cnf", "w") as fh:
        fh.write(f"p cnf {inst.pool.top} {len(inst.clauses)}\n")
        for cl in inst.clauses:
            fh.write(" ".join(map(str, cl)) + " 0\n")
    log(f"exported t7pair_dial.cnf ({inst.pool.top} vars, "
        f"{len(inst.clauses)} clauses)")


if __name__ == "__main__":
    main()
