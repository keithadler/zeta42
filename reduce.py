"""Shrink a non-k-colourable graph using UNSAT cores.

Greedy vertex deletion is hopeless here: each "may I drop this vertex?" query
is an UNSAT solve, and on G that is ~7 minutes a vertex.  A pass over 1581
vertices would take about half a year.

The cheap route is to ask the solver what it actually used.  Give every vertex
an activation literal a_v, and make only its "needs a colour" clause
conditional:

    (not a_v) or x_v0 or x_v1 or ... or x_v(k-1)

Deactivating a_v then lets v take no colour at all, which satisfies all of its
edge clauses vacuously -- exactly as if v had been deleted.  Solve under the
assumptions {a_v : v}, and on UNSAT the solver reports a *core*: the subset of
those assumptions sufficient for the contradiction.  Every vertex outside the
core is provably droppable, all in a single solve.

Iterating to a fixpoint is the same idea Heule used to get de Grey's 1581 down
to 553.  No symmetry breaking here: fixing a clique's colours would constrain
vertices the core may have dropped, which could make the reduced graph look
non-colourable when it is not.  Every result is re-verified independently.
"""

import sys
import time

from pysat.formula import IDPool
from pysat.solvers import Cadical153

from colour import k_colourable


def core_pass(graph, k):
    """One core extraction.  Returns the surviving vertex ids, or None."""
    pool = IDPool()
    n = len(graph)

    def x(v, c):
        return pool.id(("x", v, c))

    def a(v):
        return pool.id(("a", v))

    clauses = []
    for v in range(n):
        clauses.append([-a(v)] + [x(v, c) for c in range(k)])
    for i, j in graph.edges:
        for c in range(k):
            clauses.append([-x(i, c), -x(j, c)])

    lit_to_v = {a(v): v for v in range(n)}
    with Cadical153(bootstrap_with=clauses) as s:
        if s.solve(assumptions=[a(v) for v in range(n)]):
            return None                       # k-colourable: nothing to reduce
        core = s.get_core() or []
    keep = sorted({lit_to_v[abs(l)] for l in core if abs(l) in lit_to_v})
    return keep


def core_reduce(graph, k, max_rounds=20, log=print):
    """Iterate core extraction to a fixpoint, verifying as we go."""
    g = graph
    history = [len(g)]
    for rnd in range(1, max_rounds + 1):
        t0 = time.time()
        keep = core_pass(g, k)
        dt = time.time() - t0
        if keep is None:
            log(f"  round {rnd}: graph is {k}-colourable — stopping")
            return g, history
        if len(keep) >= len(g):
            log(f"  round {rnd}: no reduction ({len(g)} vertices)  [{dt:.0f}s]")
            break
        sub = g.subgraph(keep)
        # never trust the core without an independent check
        ok, _ = k_colourable(sub, k)
        if ok:
            log(f"  round {rnd}: core gave a {k}-colourable graph — rejecting, stopping")
            break
        log(f"  round {rnd}: {len(g)} -> {len(sub)} vertices, "
            f"{len(sub.edges)} edges  [{dt:.0f}s]  verified non-{k}-colourable")
        g = sub
        history.append(len(g))
    return g, history


def main():
    from degrey_g import build_G
    k = 4
    log_every = sys.stdout

    print("building G ...", flush=True)
    G = build_G()
    print(f"G: {G.summary()}\n", flush=True)

    print("core reduction:", flush=True)
    t0 = time.time()
    small, history = core_reduce(G, k, log=lambda s: print(s, flush=True))
    print(f"\nchain: {' -> '.join(map(str, history))}")
    print(f"final: {small.summary()}  in {time.time()-t0:.0f}s")
    print(f"reference: smallest known is 509 (Parts, 2019); de Grey's was 1581")

    # persist the reduced graph with exact coordinates
    import json
    fld = small.fld
    out = {
        "vertices": len(small),
        "edges": len(small.edges),
        "chi_at_least": k + 1,
        "points_exact": [[list(map(str, p[0])), list(map(str, p[1]))]
                         for p in small.points],
        "points_float": [fld.pfloat(p) for p in small.points],
        "edge_list": small.edges,
    }
    path = f"reduced-{len(small)}v.json"
    with open(path, "w") as fh:
        json.dump(out, fh)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
