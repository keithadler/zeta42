"""Rung 1b: finish the T6 rung without the aggregate marathon.

The T6 8-core aggregate dial is delegated to kissat via exported CNF (launched
separately from this script's output).  Here we run everything that was queued
behind it:

  1. rebuild the T6 8-core (from scratch; ~1 min with the numpy ball)
  2. export the aggregate instance for kissat        -> t6core_dial.cnf
  3. directed query on the designated pair (0,0)-(2,0) at k=5
  4. 9-core and 10-core sizes (does the world hold even deeper density?)
  5. |T7| and its 8/9/10-core sizes (the next rung's material, no dials)
"""

import time
from fractions import Fraction

import numpy as np
import cmath

from campaign_a import numpy_ball
from hept_g1 import arc_lattice, SCALE
from udg import UDG
from colour import ColourInstance
from twodist import dist_pairs
from rigidity import directed_query


def log(*a):
    print(*a, flush=True)


def carve_T(dist_map, dim, B, n):
    def dfb(vb):
        v = np.frombuffer(vb, dtype=np.int16)
        return dist_map.get(bytes((v - B).astype(np.int16)))
    return [vb for vb, d0 in dist_map.items()
            if (d := dfb(vb)) is not None and d0 + d <= n]


def kcore_ids(g, k):
    adj = g.adjacency()
    keep = set(range(len(g)))
    changed = True
    while changed:
        changed = False
        for v in list(keep):
            if sum(1 for u in adj[v] if u in keep) < k:
                keep.discard(v)
                changed = True
    return sorted(keep)


def main():
    t0 = time.time()
    fld, arcs, steps = arc_lattice()
    zeta = cmath.exp(2j * cmath.pi / 42)
    basis = np.array([zeta ** k for k in range(12)])
    dist_map, dim = numpy_ball(steps, 4.0, basis, SCALE, 2_000_000,
                               log=lambda *a: None)
    log(f"ball: {len(dist_map)} points  [{time.time()-t0:.0f}s]")

    B = np.zeros(dim, dtype=np.int16)
    B[0] = 2 * SCALE
    tofield = lambda vb: tuple(Fraction(int(c), SCALE)
                               for c in np.frombuffer(vb, dtype=np.int16))

    # ---- T6 core, exports, pair query ----
    T6 = carve_T(dist_map, dim, B, 6)
    g6 = UDG(fld, [tofield(vb) for vb in T6])
    log(f"T6: {g6.summary()}")
    core8 = kcore_ids(g6, 8)
    sub = g6.subgraph(core8)
    pairs = dist_pairs(sub, 4)
    log(f"8-core: {sub.summary()}, {len(pairs)} d2-pairs")

    inst = ColourInstance(sub, 5)
    for u, v in pairs:
        for c in range(5):
            inst.add_clause([-inst.var(u, c), -inst.var(v, c)])
    with open("t6core_dial.cnf", "w") as fh:
        fh.write(f"p cnf {inst.pool.top} {len(inst.clauses)}\n")
        for cl in inst.clauses:
            fh.write(" ".join(map(str, cl)) + " 0\n")
    log(f"exported t6core_dial.cnf ({inst.pool.top} vars, "
        f"{len(inst.clauses)} clauses) — launch kissat on it")

    A_id = sub.index_of(tofield(bytes(np.zeros(dim, dtype=np.int16))))
    B_id = sub.index_of(tofield(bytes(B)))
    log(f"designated pair in core: A={A_id} B={B_id}")
    if A_id is not None and B_id is not None:
        verdict, dt = directed_query(sub, 5, A_id, B_id, budget_s=900)
        log(f"directed query (0,B) at k=5: {verdict.upper()}  [{dt:.0f}s]")

    for k in (9, 10):
        ids = kcore_ids(g6, k)
        log(f"T6 {k}-core: {len(ids)} vertices")

    # ---- T7 material ----
    t0 = time.time()
    T7 = carve_T(dist_map, dim, B, 7)
    log(f"|T7| = {len(T7)}  [{time.time()-t0:.0f}s]")
    if len(T7) <= 120_000:
        t0 = time.time()
        g7 = UDG(fld, [tofield(vb) for vb in T7])
        log(f"T7 graph: {g7.summary()}  [{time.time()-t0:.0f}s]")
        for k in (8, 9, 10):
            ids = kcore_ids(g7, k)
            log(f"T7 {k}-core: {len(ids)} vertices")
    else:
        log("T7 too large for this rung's graph build — noted for rung 2")

    log("rung 1b complete")


if __name__ == "__main__":
    main()
