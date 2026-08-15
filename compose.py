"""Composition experiment: does overlapping 5-chromatic cores create stiffness?

The rigidity screen showed every known 5-chromatic graph is uniformly floppy
at k=5 -- but every measured graph was a *minimised* object, and forcing has
always come from redundancy (de Grey's M = overlapping spindles, L = stacked
hexagons).  Nobody has packed overlapping 5-chromatic cores, because until
recently there was hardly one core to pack.

Experiment: take Parts' 510, union it with copies of itself rotated about one
of its own vertices by field-internal chord rotations.  The centre vertex is
shared, and every vertex at squared radius r2 from the centre lands at unit
distance from its origin (that is what a chord-1 rotation is), so the copies
are coupled by a sheaf of new unit edges through the shared region.  Then
re-run the frequency dial and compare against the baseline 510.

Field constraint: the 510 lives in Q(sqrt3, sqrt5, sqrt11), so the available
chord rotations are r2 = 3 (sin needs sqrt11), r2 = 4 (sqrt15) and r2 = 25
(sqrt99 = 3*sqrt11).  Rotations needing sqrt7 etc. would leave the field.

Dial (numpy): for each sampled colouring, mono |= indicator_c outer itself.
Never-mono count and min-frequency fall out of the accumulated matrix.  This
scales to ~5k vertices where the pure-Python pair loop cannot.
"""

import sys
import time
from fractions import Fraction

import numpy as np

from field import chord_rotation
from parts510 import load
from rigidity import sample_colourings
from udg import UDG

FIELD_ROTS = (3, 4, 25)      # r2 with sin(t) in Q(sqrt3, sqrt5, sqrt11)


def coupling_scores(fld, g, r2):
    """How many vertices lie at exactly squared distance r2 from each vertex.

    Each one becomes a fresh unit edge when the copy rotates about that
    vertex, so this counts the coupling strength of every candidate centre.
    """
    target = fld.rat(r2)
    scores = []
    for i, c in enumerate(g.points):
        n = sum(1 for p in g.points if fld.dist2(p, c) == target)
        scores.append((n, i))
    return sorted(scores, reverse=True)


def compose(fld, g, centre_id, r2, turns=1):
    """g union its images under 1..turns applications of the rotation."""
    cos_t, sin_t = chord_rotation(fld, r2)
    centre = g.points[centre_id]
    pts = list(g.points)
    cur = g.points
    for _ in range(turns):
        cur = [fld.rotate(p, cos_t, sin_t, centre) for p in cur]
        pts.extend(cur)
    return UDG(fld, pts)       # UDG dedups exactly


def dial(g, k=5, n_samples=200, log=print):
    """Frequency dial, numpy edition.  Returns (never_mono, min_freq, median)."""
    n = len(g)
    cols = sample_colourings(g, k, n_samples, log=lambda s: None)
    freq = np.zeros((n, n), dtype=np.uint16)
    for col in cols:
        a = np.asarray(col)
        for c in range(k):
            ind = (a == c)
            freq += np.outer(ind, ind).astype(np.uint16)
    iu = np.triu_indices(n, 1)
    f = freq[iu]
    adj = np.zeros((n, n), dtype=bool)
    for i, j in g.edges:
        adj[i, j] = True
    nonadj = ~adj[iu]
    f = f[nonadj]
    never = int((f == 0).sum())
    log(f"  samples={len(cols)}  never-mono={never}/{f.size}  "
        f"min-freq={int(f.min())}/{len(cols)}  median={int(np.median(f))}")
    return never, int(f.min()), int(np.median(f))


def main():
    fld, g, _ = load("510")
    print(f"seed: {g.summary()}", flush=True)

    print("\ncoupling scores (vertices at r2 from candidate centres):", flush=True)
    best = {}
    for r2 in FIELD_ROTS:
        scores = coupling_scores(fld, g, r2)
        best[r2] = scores[0]
        print(f"  r2={r2}: best centre {scores[0][1]} couples {scores[0][0]} "
              f"vertices; top5 {[s[0] for s in scores[:5]]}", flush=True)

    print("\nbaseline dial on the 510:", flush=True)
    dial(g, log=lambda s: print(s, flush=True))

    r2, (score, centre) = max(best.items(), key=lambda kv: kv[1][0])
    print(f"\ncomposing about vertex {centre} with r2={r2} "
          f"({score} coupled vertices per turn):", flush=True)
    for turns in (1, 2, 4):
        t0 = time.time()
        big = compose(fld, g, centre, r2, turns)
        print(f"  {turns+1} copies: {big.summary()}  "
              f"[built {time.time()-t0:.0f}s]", flush=True)
        dial(big, log=lambda s: print(s, flush=True))


if __name__ == "__main__":
    main()
