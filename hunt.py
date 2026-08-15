"""Search harness: generate candidate unit-distance graphs, screen for chi >= 6.

The asymmetry that makes this worth running: SAT is *fast* when a colouring
exists.  Almost every candidate is 5-colourable and gets rejected in
milliseconds, so a sweep costs little.  The expensive case -- UNSAT -- is
exactly the case we are hoping for, and it only happens on a hit.

A recipe is a seed point set plus a list of rotations, closed under repeated
application to a bounded depth.  This generalises de Grey's construction
(triangular lattice + two chord rotations) to arbitrary rotation families.

Known limitation, measured: plain rotation-closure does not reach chi = 5.
De Grey's own rotation pair closed to 1303 vertices is still 4-colourable.
His 5-chromatic graph needs the graph M substituted into each of L's 52
H-copies, which is not implemented here.
"""

import argparse
import itertools
import json
import os
import time

from field import get_field, chord_rotation
from udg import UDG, tri_lattice, lattice_pt
from colour import k_colourable, verify_colouring

HITS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hits")


# --------------------------------------------------------------- recipes ----

class Recipe:
    """Seed points closed under a set of rotations, to a bounded depth."""

    def __init__(self, name, primes, seed_radius2, rotations, depth=2,
                 max_points=4000):
        self.name = name
        self.primes = tuple(sorted(primes))
        self.seed_radius2 = seed_radius2
        self.rotations = rotations       # list of (radius2, centre_lattice)
        self.depth = depth
        self.max_points = max_points

    def field(self):
        return get_field(self.primes)

    def build(self):
        """Returns a UDG, or None if it blew past max_points."""
        fld = self.field()
        pts = list(tri_lattice(fld, self.seed_radius2))
        ops = []
        for r2, (ca, cb) in self.rotations:
            try:
                cos_t, sin_t = chord_rotation(fld, r2)
            except ValueError:
                return None          # rotation not expressible in this field
            ops.append((cos_t, sin_t, lattice_pt(fld, ca, cb)))

        seen = set(pts)
        frontier = pts
        for _ in range(self.depth):
            new = []
            for cos_t, sin_t, centre in ops:
                for p in frontier:
                    q = fld.rotate(p, cos_t, sin_t, centre)
                    if q not in seen:
                        seen.add(q)
                        new.append(q)
                        if len(seen) > self.max_points:
                            return None
            if not new:
                break
            pts.extend(new)
            frontier = new
        return UDG(fld, pts)

    def describe(self):
        rots = ", ".join(f"r²={r2}@{c}" for r2, c in self.rotations)
        return (f"{self.name}: field Q{self.primes}, seed r²≤{self.seed_radius2}, "
                f"rotations [{rots}], depth {self.depth}")

    def to_dict(self):
        return {"name": self.name, "primes": list(self.primes),
                "seed_radius2": self.seed_radius2,
                "rotations": [[r, list(c)] for r, c in self.rotations],
                "depth": self.depth}


def field_for_rotations(radius2_list, extra=(3,)):
    """Which primes must the field contain for these chord rotations?

    A rotation moving radius² = r² by chord 1 has sin t = √(4r² - 1)/(2r²) up
    to square factors, so the square-free part of 4r² - 1 must be present.
    """
    primes = set(extra)
    for r2 in radius2_list:
        rem, d = 4 * r2 - 1, 2
        while d * d <= rem:
            while rem % (d * d) == 0:
                rem //= d * d
            d += 1
        if rem > 1:
            primes.add(rem)
    return tuple(sorted(primes))


# ---------------------------------------------------------------- screen ----

def screen(graph):
    """Fast reject ladder.  Returns (chi_lower_bound, note)."""
    t0 = time.time()
    for k in (4, 5, 6):
        ok, col = k_colourable(graph, k)
        if ok:
            assert verify_colouring(graph, col)[0], "solver returned a bad colouring"
            label = {4: "4-colourable", 5: "5-chromatic!", 6: "6-CHROMATIC"}[k]
            return k, f"{label} ({time.time()-t0:.2f}s)"
    return 7, f"not 6-colourable ({time.time()-t0:.2f}s)"


def minimise(graph, k, rounds=3):
    """Greedily drop vertices while the graph stays non-k-colourable.

    This is how de Grey's 20425 became 1581 and Heule's became 553 -- worth
    having ready, because a raw hit is always far bigger than it needs to be.
    """
    keep = list(range(len(graph)))
    for _ in range(rounds):
        progress = False
        for v in list(keep):
            trial = [u for u in keep if u != v]
            ok, _ = k_colourable(graph.subgraph(trial), k)
            if not ok:                     # still non-k-colourable: drop v
                keep = trial
                progress = True
        if not progress:
            break
    return graph.subgraph(keep)


def record_hit(graph, recipe, chi):
    os.makedirs(HITS_DIR, exist_ok=True)
    path = os.path.join(HITS_DIR, f"{recipe.name}-chi{chi}-{len(graph)}v.json")
    fld = graph.fld
    data = {
        "recipe": recipe.to_dict(),
        "chi_at_least": chi,
        "vertices": len(graph),
        "edges": len(graph.edges),
        "points_exact": [[list(map(str, p[0])), list(map(str, p[1]))]
                         for p in graph.points],
        "points_float": [fld.pfloat(p) for p in graph.points],
        "edge_list": graph.edges,
    }
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1)
    return path


# ----------------------------------------------------------------- sweep ----

def sweep(recipes, target=6, verbose=True):
    hits = []
    for r in recipes:
        g = r.build()
        if g is None:
            if verbose:
                print(f"  skip  {r.name}: too large or field mismatch")
            continue
        if len(g) < 7 or not g.edges:
            continue
        chi, note = screen(g)
        if verbose:
            mark = "HIT " if chi >= target else "    "
            print(f"  {mark}{r.name}: {g.summary()} -> chi >= {chi}  [{note}]", flush=True)
        if chi >= target:
            print(f"        recorded: {record_hit(g, r, chi)}", flush=True)
            hits.append((r, g, chi))
    return hits


def default_recipes(limit=None):
    """A sweep over chord rotations about lattice centres."""
    radii = [3, 4, 7, 12, 13, 16, 21, 28]
    centres = [(0, 0), (1, 0), (2, 0), (1, 1), (2, -1), (-2, 0), (3, 0)]
    out = []
    for seed_r2 in (4, 7, 12):
        for r2a, r2b in itertools.combinations(radii, 2):
            for ca in centres:
                for cb in centres:
                    if ca == cb and r2a == r2b:
                        continue
                    primes = field_for_rotations([r2a, r2b])
                    if len(primes) > 4:          # keep the field tractable
                        continue
                    name = f"s{seed_r2}-r{r2a}@{ca[0]},{ca[1]}-r{r2b}@{cb[0]},{cb[1]}"
                    out.append(Recipe(name, primes, seed_r2,
                                      [(r2a, ca), (r2b, cb)], depth=2))
    return out[:limit] if limit else out


def main():
    ap = argparse.ArgumentParser(description="hunt for high-chromatic unit-distance graphs")
    ap.add_argument("--limit", type=int, default=40, help="how many recipes (0 = all)")
    ap.add_argument("--target", type=int, default=6, help="record graphs with chi >= this")
    ap.add_argument("--max-points", type=int, default=4000)
    args = ap.parse_args()

    recipes = default_recipes(args.limit or None)
    for r in recipes:
        r.max_points = args.max_points
    print(f"sweeping {len(recipes)} recipes, recording chi >= {args.target}\n")
    t0 = time.time()
    hits = sweep(recipes, target=args.target)
    print(f"\ndone in {time.time()-t0:.1f}s; {len(hits)} hit(s)")


if __name__ == "__main__":
    main()
