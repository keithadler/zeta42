"""SAT-based colouring queries on unit-distance graphs.

Everything routes through CaDiCaL via pysat.  The encoding is the standard
one-hot per vertex with clique symmetry breaking, which is what makes 4- and
5-colourability decidable in seconds on graphs of a few thousand vertices.
"""

from pysat.formula import IDPool
from pysat.solvers import Cadical153


class ColourInstance:
    """A k-colouring instance over a UDG, with optional extra constraints."""

    def __init__(self, graph, k, break_symmetry=True):
        self.graph = graph
        self.k = k
        self.pool = IDPool()
        self.clauses = []
        n = len(graph)

        # one-hot: every vertex gets at least one colour
        for v in range(n):
            self.clauses.append([self.var(v, c) for c in range(k)])

        # proper colouring: adjacent vertices differ
        for i, j in graph.edges:
            for c in range(k):
                self.clauses.append([-self.var(i, c), -self.var(j, c)])

        # a vertex needs no at-most-one clause for k-colourability, but fixing
        # a clique's colours prunes the k! colour permutations
        if break_symmetry:
            for pos, v in enumerate(self._greedy_clique()[:k]):
                self.clauses.append([self.var(v, pos)])

    def var(self, v, c):
        return self.pool.id(("x", v, c))

    def _greedy_clique(self):
        adj = [set(a) for a in self.graph.adjacency()]
        order = sorted(range(len(self.graph)), key=lambda v: -len(adj[v]))
        clique = []
        for v in order:
            if all(v in adj[u] for u in clique):
                clique.append(v)
        return clique

    # ------------------------------------------------------------ extras ----
    def forbid_monochromatic(self, triple):
        """Rule out all vertices of `triple` sharing a colour."""
        for c in range(self.k):
            self.clauses.append([-self.var(v, c) for v in triple])

    def add_clause(self, clause):
        self.clauses.append(clause)

    # ------------------------------------------------------------- solve ----
    def solve(self):
        """Returns (sat, colouring) -- colouring is a list of ints or None."""
        with Cadical153(bootstrap_with=self.clauses) as s:
            if not s.solve():
                return False, None
            model = set(l for l in s.get_model() if l > 0)
            colouring = []
            for v in range(len(self.graph)):
                col = next((c for c in range(self.k) if self.var(v, c) in model), None)
                colouring.append(col)
            return True, colouring


def k_colourable(graph, k):
    return ColourInstance(graph, k).solve()


def chromatic_number(graph, lo=1, hi=8):
    """Smallest k with a proper k-colouring (searching upward from lo)."""
    for k in range(lo, hi + 1):
        ok, _ = k_colourable(graph, k)
        if ok:
            return k
    return None


def verify_colouring(graph, colouring):
    """Independent check that a returned colouring really is proper."""
    for i, j in graph.edges:
        if colouring[i] == colouring[j]:
            return False, (i, j)
    return True, None
