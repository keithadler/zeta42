"""Campaign A+B: the distance-2 forcing hunt in the heptagonal world,
measured by MaxSAT (exact minimum mono-pair count), not sampling.

Facts this stands on (see FINDINGS.md):
  - chi({1,2}) >= 6 is known (Exoo--Ismailescu), so a unit-only gadget forcing
    a designated distance-2 pair non-mono at k=5 completes a chi(R^2) >= 6
    proof exactly as B-at-phi would.
  - The zeta_42 arc world (84 arcs over denominator 7) is the only arithmetic
    measured to sustain forcing density (Haugland's G1 template works there
    at k=4); distance-2 pairs are native to it.
  - The only aggregate strain ever observed at k=5 was at d=2 (the 510 took
    28s to dodge; everything else 0.0s).

The dial (B): minimum number of monochromatic d=2 pairs over ALL proper
5-colourings, via RC2 MaxSAT.  Exact, monotone under growth, bias-free.
"Aggregate forcer" == minimum > 0.

Phases:
  0  MaxSAT baseline on the 510 (expect minimum 0; calibrates cost)
  1  numpy BFS of the arc lattice to radius 4 (~1M points)
  2  T-sets between (0,0) and (2,0), density filters, cores
  3  dials on each candidate: directed pair query + MaxSAT minimum
"""

import time
from fractions import Fraction

import numpy as np

from pysat.formula import WCNF, IDPool
from pysat.examples.rc2 import RC2

MAXSAT_VERTEX_CAP = 6500


# ---------------------------------------------------------------- dial ----

def min_mono_pairs(graph, pairs, k=5, log=print):
    """Exact minimum number of monochromatic `pairs` over proper k-colourings.

    Instrument hierarchy: a plain SAT call answers "is the minimum 0?" far
    faster than core-guided MaxSAT proves it.  MaxSAT runs only on strain
    (SAT says no zero-mono colouring exists), which is when its exact value
    is worth paying for.
    """
    if len(graph) > MAXSAT_VERTEX_CAP:
        log(f"  [dial skipped: {len(graph)} vertices > cap {MAXSAT_VERTEX_CAP}]")
        return None
    from colour import ColourInstance
    inst = ColourInstance(graph, k)
    for u, v in pairs:
        for c in range(k):
            inst.add_clause([-inst.var(u, c), -inst.var(v, c)])
    t0 = time.time()
    sat, _ = inst.solve()
    if sat:
        log(f"  minimum mono pairs (k={k}): 0 (SAT dodge)  [{time.time()-t0:.0f}s]")
        return 0
    log(f"  STRAIN: no zero-mono colouring (UNSAT {time.time()-t0:.0f}s) "
        f"— computing exact minimum via MaxSAT")
    pool = IDPool()
    x = lambda v, c: pool.id(("x", v, c))
    m = lambda i: pool.id(("m", i))
    w = WCNF()
    for v in range(len(graph)):
        w.append([x(v, c) for c in range(k)])
    for i, j in graph.edges:
        for c in range(k):
            w.append([-x(i, c), -x(j, c)])
    for idx, (u, v) in enumerate(pairs):
        for c in range(k):
            w.append([-x(u, c), -x(v, c), m(idx)])
        w.append([-m(idx)], weight=1)
    t0 = time.time()
    with RC2(w) as rc2:
        model = rc2.compute()
        cost = rc2.cost if model is not None else None
    log(f"  MaxSAT minimum mono pairs (k={k}): {cost}  [{time.time()-t0:.0f}s]")
    return cost


# ------------------------------------------------------------- lattice ----

def numpy_ball(steps_int, radius_float, basis_complex, scale, max_pts, log=print):
    """BFS ball of integer step-vectors under a float radius cap, vectorised."""
    t0 = time.time()
    steps = np.array(steps_int, dtype=np.int16)
    dim = steps.shape[1]
    seen = {bytes(np.zeros(dim, dtype=np.int16))}
    frontier = np.zeros((1, dim), dtype=np.int16)
    dist = {bytes(frontier[0]): 0}
    r = 0
    while len(seen) < max_pts:
        r += 1
        cand = (frontier[:, None, :] + steps[None, :, :]).reshape(-1, dim)
        cand = np.unique(cand, axis=0)
        # radius filter, vectorised through the complex embedding
        z = cand.astype(np.float64) @ basis_complex / scale
        cand = cand[np.abs(z) <= radius_float]
        new = []
        for row in cand:
            key = bytes(row)
            if key not in seen:
                seen.add(key)
                dist[key] = r
                new.append(row)
        if not new:
            break
        frontier = np.array(new, dtype=np.int16)
        log(f"  radius {r}: +{len(new)} (total {len(seen)})  [{time.time()-t0:.0f}s]")
    return dist, dim


def main():
    log = lambda *a: print(*a, flush=True)

    # ---- phase 0: baseline dial on the 510 ----
    log("phase 0: MaxSAT baseline (510, d=2 pairs, k=5)")
    from parts510 import load
    from twodist import dist_pairs
    _, g510, _ = load("510")
    pairs510 = dist_pairs(g510, 4)
    log(f"  510: {len(pairs510)} distance-2 pairs")
    min_mono_pairs(g510, pairs510, 5, log=log)

    # ---- phase 1: heptagonal arc ball, radius 4 ----
    log("\nphase 1: arc-lattice ball, radius 4")
    from hept_g1 import arc_lattice, SCALE
    fld, arcs, steps = arc_lattice()
    import cmath
    zeta = cmath.exp(2j * cmath.pi / 42)
    basis = np.array([zeta ** k for k in range(12)])
    dist_map, dim = numpy_ball(steps, 4.0, basis, SCALE, 2_000_000, log=log)
    log(f"  ball complete: {len(dist_map)} points")

    # ---- phase 2: T-sets between (0,0) and (2,0) ----
    log("\nphase 2: T-sets for the designated pair (0,0)-(2,0)")
    B = np.zeros(dim, dtype=np.int16)
    B[0] = 2 * SCALE
    Bkey = bytes(B)
    dB = dist_map.get(Bkey)
    log(f"  d(0, B) = {dB}")

    # distance from B = distance of (v - B) from 0 by symmetry of the arc set
    def dist_from_B(vbytes):
        v = np.frombuffer(vbytes, dtype=np.int16)
        return dist_map.get(bytes((v - B).astype(np.int16)))

    for n in (6, 7):
        t0 = time.time()
        T = [vb for vb, d0 in dist_map.items()
             if d0 is not None and (dfb := dist_from_B(vb)) is not None
             and d0 + dfb <= n]
        log(f"  |T{n}| = {len(T)}  [{time.time()-t0:.0f}s]")
        if len(T) < 50 or len(T) > 60000:
            continue

        # build exact UDG on T, then cores
        from udg import UDG
        tofield = lambda vb: tuple(
            Fraction(int(c), SCALE)
            for c in np.frombuffer(vb, dtype=np.int16))
        pts = [tofield(vb) for vb in T]
        t0 = time.time()
        g = UDG(fld, pts)
        log(f"  T{n} graph: {g.summary()}  [{time.time()-t0:.0f}s]")

        adj = g.adjacency()
        for core_k in (8, 9, 10):
            keep = set(range(len(g)))
            changed = True
            while changed:
                changed = False
                for v in list(keep):
                    if sum(1 for u in adj[v] if u in keep) < core_k:
                        keep.discard(v)
                        changed = True
            log(f"  {core_k}-core: {len(keep)} vertices")
            if not keep:
                break
            sub = g.subgraph(sorted(keep))
            A_id = sub.index_of(tofield(bytes(np.zeros(dim, dtype=np.int16))))
            B_id = sub.index_of(tofield(Bkey))
            log(f"    pair present: A={A_id is not None} B={B_id is not None}")
            # phase 3 dials
            from twodist import dist_pairs as dp
            d2pairs = dp(sub, 4)
            log(f"    distance-2 pairs in core: {len(d2pairs)}")
            if d2pairs:
                min_mono_pairs(sub, d2pairs, 5, log=log)
            if A_id is not None and B_id is not None:
                from rigidity import directed_query
                verdict, dt = directed_query(sub, 5, A_id, B_id, budget_s=300)
                log(f"    designated pair directed (k=5): {verdict.upper()} [{dt:.0f}s]")

    log("\ncampaign rung complete")


if __name__ == "__main__":
    main()
