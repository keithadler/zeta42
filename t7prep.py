"""Export the T7 10-core aggregate instance (rung 2's arena)."""
import time
from fractions import Fraction
import numpy as np, cmath
from campaign_a import numpy_ball
from hept_g1 import arc_lattice, SCALE
from udg import UDG
from colour import ColourInstance
from twodist import dist_pairs
from rung1b import carve_T, kcore_ids

t0 = time.time()
fld, arcs, steps = arc_lattice()
basis = np.array([cmath.exp(2j*cmath.pi/42)**k for k in range(12)])
dist_map, dim = numpy_ball(steps, 4.0, basis, SCALE, 2_000_000, log=lambda *a: None)
print(f"ball: {len(dist_map)}  [{time.time()-t0:.0f}s]", flush=True)
B = np.zeros(dim, dtype=np.int16); B[0] = 2*SCALE
T7 = carve_T(dist_map, dim, B, 7)
tofield = lambda vb: tuple(Fraction(int(c), SCALE) for c in np.frombuffer(vb, dtype=np.int16))
g = UDG(fld, [tofield(vb) for vb in T7])
print(f"T7: {g.summary()}  [{time.time()-t0:.0f}s]", flush=True)
core = g.subgraph(kcore_ids(g, 10))
pairs = dist_pairs(core, 4)
print(f"10-core: {core.summary()}, {len(pairs)} d2-pairs  [{time.time()-t0:.0f}s]", flush=True)
inst = ColourInstance(core, 5)
for u, v in pairs:
    for c in range(5):
        inst.add_clause([-inst.var(u,c), -inst.var(v,c)])
with open("t7core_dial.cnf", "w") as fh:
    fh.write(f"p cnf {inst.pool.top} {len(inst.clauses)}\n")
    for cl in inst.clauses:
        fh.write(" ".join(map(str, cl)) + " 0\n")
print(f"exported t7core_dial.cnf ({inst.pool.top} vars, {len(inst.clauses)} clauses)", flush=True)
A_id = core.index_of(tofield(bytes(np.zeros(dim, dtype=np.int16))))
B_id = core.index_of(tofield(bytes(B)))
print(f"designated pair in 10-core: A={A_id} B={B_id}", flush=True)
