"""Independence number of Parts' 510-vertex graph.

Why: chi_f(plane) >= chi_f(G) >= n/alpha(G) for any finite unit-distance G.
The current record lower bound for the fractional chromatic number of the
plane is "strictly above 4" (June 2026, via a graph with independence ratio
below 1/4).  If alpha(510) <= 127, the 510's ratio also beats 1/4 and the
minimised 5-chromatic lineage is fractionally strong; if alpha >= 128 the two
notions of strength come apart.  Either answer settles the alignment question.

Method: greedy for a lower bound, then SAT with a cardinality constraint,
descending one step at a time from the first interesting bound.  Each query
"is there an independent set of size >= k" is a fresh instance with a
sequential-counter encoding; UNSAT at k pins alpha = k - 1 once a set of size
k - 1 is in hand.
"""

import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

from parts510 import load


def greedy_lb(g):
    """Greedy independent set (min-degree order) for a starting lower bound."""
    adj = [set(a) for a in g.adjacency()]
    alive = set(range(len(g)))
    picked = []
    while alive:
        v = min(alive, key=lambda u: len(adj[u] & alive))
        picked.append(v)
        alive -= adj[v] | {v}
    return picked


def exists_independent(g, k, budget_conflicts=None, log=print):
    """Is there an independent set of size >= k?  True/False/None(budget)."""
    n = len(g)
    pool = IDPool()
    lits = [pool.id(("a", v)) for v in range(n)]
    clauses = [[-lits[i], -lits[j]] for i, j in g.edges]
    card = CardEnc.atleast(lits=lits, bound=k, vpool=pool,
                           encoding=EncType.seqcounter)
    t0 = time.time()
    with Cadical153(bootstrap_with=clauses + card.clauses) as s:
        if budget_conflicts:
            s.conf_budget(budget_conflicts)
            r = s.solve_limited()
        else:
            r = s.solve()
        model = s.get_model() if r else None
    dt = time.time() - t0
    log(f"  independent set >= {k}: {r}  [{dt:.1f}s]")
    if r and model:
        pos = set(l for l in model if l > 0)
        return True, [v for v in range(n) if lits[v] in pos]
    return r, None


def main():
    _, g, _ = load("510")
    n = len(g)
    print(f"graph: {g.summary()}", flush=True)

    lb_set = greedy_lb(g)
    lb = len(lb_set)
    print(f"greedy independent set: {lb}  (ratio {lb/n:.4f})", flush=True)
    print(f"threshold of interest: alpha <= 127 <=> ratio < 1/4 "
          f"<=> chi_f >= {n}/127 = {n/127:.4f}", flush=True)

    # verify greedy set is truly independent, exactly
    es = set(g.edges)
    assert not any((min(a, b), max(a, b)) in es
                   for i, a in enumerate(lb_set) for b in lb_set[i + 1:])

    best_set = lb_set
    k = max(lb + 1, 128) if lb < 128 else lb + 1
    # first settle the interesting threshold, then walk to the exact value
    for probe in sorted(set([128, k])):
        r, found = exists_independent(g, probe, log=lambda m: print(m, flush=True))
        if r:
            best_set, lb = found, probe
    if lb >= 128:
        print(f"alpha >= {lb}: ratio >= {lb/n:.4f} >= 1/4 -- "
              f"fractionally weak; alignment FAILS", flush=True)
        # still walk up to the exact alpha for the record
        k = lb + 1
        while True:
            r, found = exists_independent(g, k, log=lambda m: print(m, flush=True))
            if not r:
                print(f"alpha({n}) = {k-1}", flush=True)
                break
            best_set, k = found, k + 1
    else:
        print(f"independent set of 128 ruled out: alpha <= 127, "
              f"chi_f(plane) >= {n/127:.4f} via this graph", flush=True)
        k = 127
        while True:
            r, found = exists_independent(g, k, log=lambda m: print(m, flush=True))
            if r:
                print(f"alpha({n}) = {k}  (ratio {k/n:.4f}, "
                      f"chi_f >= {n/k:.4f})", flush=True)
                break
            k -= 1
    print("done", flush=True)


if __name__ == "__main__":
    main()
