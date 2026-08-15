"""Stochastic pre-filter: a Markov chain on colourings.

The state is an assignment of colours to vertices; a step picks a vertex
involved in a conflict and recolours it, usually to whichever colour minimises
its conflicts, occasionally at random.  That is a random walk on colouring
space (TabuCol-style local search), and it finds proper colourings of large
graphs far faster than a complete solver does.

The catch, which is the whole reason this is a *filter* and not a method:
it is one-sided.  Finding a colouring proves chi <= k.  Failing to find one
proves nothing at all -- the walk may simply have missed it.  So a hunt looks
like this:

    local search finds a k-colouring   -> reject the candidate, cheaply, done
    local search fails                 -> escalate to SAT for a real answer

Since almost every candidate IS colourable, almost every candidate is disposed
of by the cheap half, and the expensive complete solver only ever runs on the
small set of things that might be interesting.
"""

import random


def local_search_colour(graph, k, max_steps=200_000, noise=0.03, seed=None):
    """Try to find a proper k-colouring by random walk.

    Returns (colouring, steps) on success, or (None, steps) if the walk ran out
    of steps -- which is NOT a proof that no colouring exists.
    """
    rng = random.Random(seed)
    n = len(graph)
    adj = graph.adjacency()
    if n == 0:
        return [], 0

    colour = [rng.randrange(k) for _ in range(n)]

    # conflict[v][c] = how many neighbours of v currently have colour c
    conflict = [[0] * k for _ in range(n)]
    for v in range(n):
        for u in adj[v]:
            conflict[v][colour[u]] += 1

    def bad_vertices():
        return [v for v in range(n) if conflict[v][colour[v]] > 0]

    bad = bad_vertices()
    for step in range(max_steps):
        if not bad:
            return colour, step
        v = bad[rng.randrange(len(bad))]
        if conflict[v][colour[v]] == 0:
            bad = bad_vertices()
            continue

        if rng.random() < noise:
            new = rng.randrange(k)
        else:
            best = min(conflict[v])
            new = rng.choice([c for c in range(k) if conflict[v][c] == best])

        old = colour[v]
        if new == old:
            new = (old + 1 + rng.randrange(k - 1)) % k if k > 1 else old
        colour[v] = new
        for u in adj[v]:
            conflict[u][old] -= 1
            conflict[u][new] += 1

        # refresh the worklist periodically rather than maintaining it exactly
        if step % 64 == 0:
            bad = bad_vertices()
        else:
            touched = [v] + adj[v]
            for u in touched:
                if conflict[u][colour[u]] > 0:
                    if u not in bad:
                        bad.append(u)
            bad = [u for u in bad if conflict[u][colour[u]] > 0]

    return None, max_steps


def screen_fast(graph, k, tries=3, max_steps=200_000):
    """Is `graph` k-colourable?  Returns True, or None for 'don't know'."""
    for t in range(tries):
        col, _ = local_search_colour(graph, k, max_steps=max_steps, seed=t)
        if col is not None:
            return True
    return None
