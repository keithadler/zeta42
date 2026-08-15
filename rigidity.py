"""Rigidity sweep: measure how much freedom a graph's k-colourings leave.

For every vertex pair (u, v), the question is whether some proper k-colouring
makes u and v the same colour.  Pairs that never can are *forced pairs* -- the
raw material of chromatic-number lower bounds (spindle a forced pair and the
chromatic number climbs).  Pairs that barely can are where a forcing gadget
might be grown.

Three phases, cheapest first, so almost all pairs are disposed of for nearly
nothing:

  A. sampling   generate diverse k-colourings; any pair seen monochromatic in
                any sample is answered ("can be mono") for free.
  B. directed   each surviving pair gets its own SAT query with
                colour(u) == colour(v) forced.  SAT kills it; UNSAT or
                timeout advances it.
  C. patterns   for the rare survivors, enumerate the essentially different
                ways the graph supports the pair being mono, projected onto
                the pair's neighbourhood, up to a cap.  The count is the
                stiffness measure: 0 = forced pair; small = hill-climbable.

Usage:
    ./.venv/bin/python rigidity.py --graph 510 --k 5 --samples 300
    ./.venv/bin/python rigidity.py --graph g1 --k 4 --pair-check   # seeded control
"""

import argparse
import json
import os
import time
from itertools import combinations

from pysat.formula import IDPool
from pysat.solvers import Cadical153

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------ base clauses --

def colour_clauses(graph, k, pool):
    def x(v, c):
        return pool.id(("x", v, c))

    cls = []
    for v in range(len(graph)):
        cls.append([x(v, c) for c in range(k)])
    for i, j in graph.edges:
        for c in range(k):
            cls.append([-x(i, c), -x(j, c)])
    return cls, x


def model_colouring(model, n, k, x):
    pos = set(l for l in model if l > 0)
    return [next(c for c in range(k) if x(v, c) in pos) for v in range(n)]


# ------------------------------------------------------------- phase A ------

def sample_colourings(graph, k, n_samples, log=print):
    """Diverse proper k-colourings via solve-and-block.  Returns colourings."""
    import random
    rng = random.Random(0)
    pool = IDPool()
    cls, x = colour_clauses(graph, k, pool)
    n = len(graph)
    out = []
    t0 = time.time()
    with Cadical153(bootstrap_with=cls) as s:
        for i in range(n_samples):
            # random assumptions push each sample into a different region of
            # colouring space; plain solve-and-block returns near-clones
            assume = [x(rng.randrange(n), rng.randrange(k)) for _ in range(8)]
            if not s.solve(assumptions=assume):
                if not s.solve():
                    log(f"  colouring space exhausted after {len(out)} samples")
                    break
            col = model_colouring(s.get_model(), n, k, x)
            out.append(col)
            s.add_clause([-x(v, col[v]) for v in range(n)])
    log(f"  {len(out)} samples in {time.time()-t0:.1f}s")
    return out


def mono_seen(graph, colourings):
    """Set of pairs observed monochromatic in at least one colouring."""
    seen = set()
    n = len(graph)
    for col in colourings:
        byc = {}
        for v, c in enumerate(col):
            byc.setdefault(c, []).append(v)
        for verts in byc.values():
            for a, b in combinations(verts, 2):
                seen.add((a, b))
    return seen


# ------------------------------------------------------------- phase B ------

def directed_query(graph, k, u, v, budget_s=60):
    """Can u and v share a colour in some proper k-colouring?

    Returns ('sat'|'unsat'|'timeout', seconds).
    """
    pool = IDPool()
    cls, x = colour_clauses(graph, k, pool)
    for c in range(k):
        cls.append([-x(u, c), x(v, c)])
        cls.append([x(u, c), -x(v, c)])
    t0 = time.time()
    with Cadical153(bootstrap_with=cls) as s:
        # pysat CaDiCaL has no wall-clock budget; conflict budget approximates
        s.conf_budget(int(budget_s * 40_000))       # ~40k conflicts/s observed
        r = s.solve_limited()
    dt = time.time() - t0
    if r is True:
        return "sat", dt
    if r is False:
        return "unsat", dt
    return "timeout", dt


# ------------------------------------------------------------- phase C ------

def pattern_count(graph, k, u, v, cap=100, budget_s=600, log=print):
    """How many essentially different local patterns support u,v being mono?

    A pattern is the colouring restricted to N(u) ∪ N(v) ∪ {u, v}, up to
    nothing -- colour permutations are broken by fixing colour(u) = 0.
    Enumerate by solve-and-block on the projection.  Returns (count, complete):
    complete=False means the cap or budget was hit.
    """
    adj = graph.adjacency()
    region = sorted(set(adj[u]) | set(adj[v]) | {u, v})
    pool = IDPool()
    cls, x = colour_clauses(graph, k, pool)
    for c in range(k):
        cls.append([-x(u, c), x(v, c)])
        cls.append([x(u, c), -x(v, c)])
    cls.append([x(u, 0)])                       # break colour permutation
    t0 = time.time()
    count = 0
    with Cadical153(bootstrap_with=cls) as s:
        while count < cap and time.time() - t0 < budget_s:
            if not s.solve():
                return count, True              # exhausted: exact count
            model = set(l for l in s.get_model() if l > 0)
            # block this projected pattern
            s.add_clause([-x(w, c) for w in region for c in range(k)
                          if x(w, c) in model])
            count += 1
    return count, False


# --------------------------------------------------------------- sweep ------

def sweep(graph, k, n_samples=300, budget_s=60, log=print):
    n = len(graph)
    all_pairs = n * (n - 1) // 2
    log(f"sweep: {n} vertices, k={k}, {all_pairs} pairs")

    log("phase A: sampling")
    cols = sample_colourings(graph, k, n_samples, log=log)
    seen = mono_seen(graph, cols)
    # adjacent pairs can never be mono; exclude them from 'survivors'
    adjacent = set(graph.edges)
    survivors = [(a, b) for a, b in combinations(range(n), 2)
                 if (a, b) not in seen and (a, b) not in adjacent]
    log(f"  seen-mono: {len(seen)}  adjacent: {len(adjacent)}  "
        f"survivors: {len(survivors)}")

    log("phase B: directed queries")
    results = []
    t0 = time.time()
    for i, (u, v) in enumerate(survivors):
        verdict, dt = directed_query(graph, k, u, v, budget_s)
        if verdict != "sat":
            results.append({"pair": [u, v], "verdict": verdict, "seconds": dt})
            log(f"  [{i+1}/{len(survivors)}] ({u},{v}): {verdict.upper()} {dt:.1f}s")
        elif dt > 5:
            results.append({"pair": [u, v], "verdict": "slow-sat", "seconds": dt})
    log(f"  done in {time.time()-t0:.0f}s; {len(results)} interesting pair(s)")

    log("phase C: pattern counts for non-sat survivors")
    for r in results:
        if r["verdict"] in ("unsat", "timeout"):
            u, v = r["pair"]
            cnt, complete = pattern_count(graph, k, u, v, log=log)
            r["patterns"] = cnt
            r["patterns_complete"] = complete
            log(f"  ({u},{v}): {cnt}{'' if complete else '+'} patterns")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="510",
                    help="510/517/529/553 (CNP-SAT), 'g1', or 'degreyG'")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--samples", type=int, default=300)
    ap.add_argument("--budget", type=float, default=60)
    ap.add_argument("--pair-check", action="store_true",
                    help="only run the seeded control pair (g1: A,B)")
    args = ap.parse_args()

    if args.graph == "g1":
        from fractions import Fraction
        from hept_g1 import load_G1, SCALE
        from udg import UDG
        fld, core, B, _ = load_G1(log=print)
        tofield = lambda v: tuple(Fraction(x, SCALE) for x in v)
        g = UDG(fld, [tofield(v) for v in sorted(core)])
        control = (g.index_of(tofield((0,) * 12)), g.index_of(tofield(B)))
    elif args.graph == "degreyG":
        from degrey_g import build_G
        g, control = build_G(), None
    else:
        from parts510 import load
        _, g, _ = load(args.graph)
        control = None
    print(f"graph {args.graph}: {g.summary()}, k={args.k}")

    if args.pair_check:
        assert control and None not in control, "no control pair for this graph"
        u, v = control
        print(f"seeded control: pair ({u},{v}) expecting zero patterns at k={args.k}")
        cnt, complete = pattern_count(g, args.k, u, v)
        print(f"patterns: {cnt}, complete: {complete}")
        print("PASS" if (cnt, complete) == (0, True) else "FAIL")
        return

    results = sweep(g, args.k, args.samples, args.budget)
    path = os.path.join(HERE, f"rigidity-{args.graph}-k{args.k}.json")
    with open(path, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
