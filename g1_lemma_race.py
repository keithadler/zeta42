"""Second attempt at G1's forced-pair lemma, racing the original encoding.

Differences from the first run (still going on another core):
  - symmetry breaking ON: fixing a greedy clique's colours is sound here
    because colour permutations still act freely on colourings satisfying
    colour(A) == colour(B), so it prunes without changing satisfiability.
  - output is flushed line by line, so progress is visible from the log.
  - G1 comes from the cache (and this run writes the cache if absent).

Expected result: UNSAT, i.e. no 4-colouring of G1 gives A and B the same
colour.  Whichever encoding finishes first, the other becomes a cross-check.
"""

import sys
import time
from fractions import Fraction

from hept_g1 import load_G1, SCALE
from udg import UDG
from colour import ColourInstance, k_colourable


def log(msg):
    print(msg, flush=True)


def main():
    t0 = time.time()
    fld, core, B, _ = load_G1(log=log)
    pts_int = sorted(core)

    def tofield(v):
        return tuple(Fraction(x, SCALE) for x in v)

    g = UDG(fld, [tofield(v) for v in pts_int])
    log(f"G1: {g.summary()}  [{time.time()-t0:.0f}s]")

    A_id = g.index_of(tofield((0,) * 12))
    B_id = g.index_of(tofield(B))
    assert A_id is not None and B_id is not None
    log(f"A -> vertex {A_id}, B -> vertex {B_id}")

    t0 = time.time()
    ok, _ = k_colourable(g, 4)
    log(f"unconstrained 4-colourable: {ok}  [{time.time()-t0:.1f}s]  (expected True)")

    inst = ColourInstance(g, 4, break_symmetry=True)
    for c in range(4):
        inst.add_clause([-inst.var(A_id, c), inst.var(B_id, c)])
        inst.add_clause([inst.var(A_id, c), -inst.var(B_id, c)])
    log("solving forced colour(A) == colour(B) ...")
    t0 = time.time()
    sat, _ = inst.solve()
    log(f"forced-equal instance: {'SAT' if sat else 'UNSAT'}  [{time.time()-t0:.1f}s]")
    log(f"lemma {'HOLDS' if not sat else 'FAILS'}: forced pair "
        f"{'confirmed' if not sat else 'REFUTED -- investigate!'}")


if __name__ == "__main__":
    main()
