# cnp — exact-arithmetic tooling for the chromatic number of the plane

A small research codebase for computational work on the Hadwiger–Nelson
problem (what is the chromatic number of the plane? Known: 5 ≤ χ(ℝ²) ≤ 7).

**What this is:** working infrastructure — exact arithmetic in the number
fields where unit-distance constructions live, SAT-based colouring queries,
and reproductions of the key published objects — plus a documented,
instrument-validated search that **found nothing**. The negative results and
the methods are the content. No bound is moved here.

**What this is not:** a contribution to the mathematics of the problem. See
[`FINDINGS.md`](FINDINGS.md) for the honest ledger and
[`note.md`](note.md) for the write-up, including one open question
(why the denominator-7 heptagonal arc system is coincidence-dense while the
denominator-11 pentagonal one is thin) that we believe is genuinely open.

## The one sharp statement

By results already in the literature (Parts 2020; Exoo–Ismailescu), the lower
bound χ(ℝ²) ≥ 6 would follow from a single finite object **B**: a
unit-distance graph with a designated pair at distance φ (or 2, or √3…, any
distance d with a published two-distance χ({1,d}) ≥ 6 result) that receives
distinct colours in every proper 5-colouring. Substitute B into each d-edge
of the published two-distance witness and the union is a 6-chromatic
unit-distance graph. This repo documents where B **is not**: every known
5-chromatic graph, their compositions, two arithmetic worlds' natural
construction families, and the purpose-built dense cores we could reach on a
laptop.

## Layout

Exact arithmetic (no floating-point edge decisions anywhere — floats only
prescreen candidate pairs; every edge is confirmed in the field):

| file | contents |
|---|---|
| `field.py` | real multiquadratic fields ℚ(√a, √b, …); chord rotations |
| `cyclo.py` | cyclotomic fields ℚ(ζₙ) via the power basis |
| `cycloext.py` | real quadratic extensions ℚ(ζₙ, √m) |
| `ext_imag.py` | imaginary quadratic extensions ℚ(ζₙ, √−D) |
| `udg.py` | unit-distance graphs over any of the above; exact edge detection |
| `colour.py` | SAT colouring queries (CaDiCaL via PySAT); verification |

Reproductions of published objects (all counts match the papers; see each
module's docstring for the construction and source):

| file | object |
|---|---|
| `degrey.py` | de Grey's J/K/L assembly + his mono-triple lemma (verified) |
| `degrey_g.py` | de Grey's 1581-vertex graph G; χ(G) = 5 verified here |
| `heptagon.py` | Haugland's 21-vertex ζ₄₂ seed and its 84 arcs |
| `hept_g1.py` | Haugland's G₁ (740v); forced-pair lemma verified, 2 engines |
| `parts510.py` | Parts' 510-vertex graph from Heule's CNP-SAT data, re-derived |
| `pentagon.py` | de Grey's pentagonal G₁₂₆ (two-distance world) |
| `g3_verify.py` | Haugland's G₃ (2131v) in ℚ(ζ₄₂, √11) |

Instruments and experiments (each module's docstring states what it
measured and what the result was):

| file | experiment |
|---|---|
| `rigidity.py` | pair-rigidity sweeps: sampling, directed queries, pattern counts |
| `twodist.py` | aggregate mono-distance forcing tests |
| `compose.py` | composition/coupling experiments on the 510 |
| `census.py` | cross-field unit-coincidence census |
| `alpha510.py` | independence number probe (fractional bound diagnostic) |
| `pent_b.py`, `pent_b2.py` | Haugland-template transplant to ℚ(ζ₅) (cores empty) |
| `zeta5_margin/` | ζ(5): sharper constants for Fauzan's proof — certified comparison measure (μ 260→191) plus wider inner prime range (μ →148) |
| `zeta20.py` | ramified-and-split test in ℚ(ζ₂₀) with native φ-pairs: dense 4-chromatic cores, no forcing |
| `nonagon.py` | ninefold analogue in ℚ(ζ₉): no seed exists; arc census; template cores empty or 3-chromatic |
| `spindle15.py` | Moser-twisted ζ₁₅ world (both φ-pairs and triangles) |
| `campaign_a.py`, `rung1b.py`, `t7prep.py`, `t7pair_prep.py` | the distance-2 campaign in the heptagonal world |
| `cnc.py`, `cnc2.py` | cube-and-conquer experiments (negative: not profitable here) |
| `reduce.py` | UNSAT-core graph shrinking |
| `hunt.py`, `fast.py` | early experiments (rotation closures; local search) |
| `draw.py` | SVG rendering from exact coordinates |
| `test_baseline.py` | run this first: every known-answer check in one script |

## Reproducing

```
python3 -m venv .venv
./.venv/bin/pip install python-sat numpy sympy mpmath
./.venv/bin/python test_baseline.py          # all published counts + lemmas, ~1s
./.venv/bin/python test_baseline.py --slow   # + chi(G)=5 UNSAT, ~7 min
```

kissat (used for the long UNSAT runs) installs via `brew install kissat` /
your package manager. External data: clone
[Heule's CNP-SAT](https://github.com/marijnheule/CNP-SAT) into
`external/CNP-SAT` for the 510/517/529/553 coordinates, and
[Heule's CnC](https://github.com/marijnheule/CnC) into `external/CnC` for
march_cu/iglucose (only needed by `cnc*.py`).

Large generated artifacts (CNF exports, caches, solver logs) are not
committed; every one regenerates from the scripts above.

## The deep end: ℚ(ζ₄₂)

Of every arithmetic we searched, exactly one sustains the density that
forcing arguments need: the heptagonal world of ℚ(ζ₄₂), whose 84 unit arcs
share denominator 7 and generate a coincidence-rich lattice (this is the
arithmetic underlying Haugland 2026). Everything else we tried — ℚ(ζ₅)'s
90-arc system, ℚ(ζ₁₅), quasicrystal slices — runs thin (measured; see
FINDINGS.md items 11–13). Exact facts from the far end of what a laptop
could build there:

- The arc lattice holds **10,223,809** exact points within radius 4
  (`campaign_a.py`); every coordinate a 12-tuple of rationals over
  denominator 7, every edge decided in the field.
- The path-set T₇ between the designated distance-2 pair (0,0)–(2,0)
  yields a **10-core of 18,959 vertices / 174,597 edges** with the pair
  surviving inside it — the densest level-5 arena we know how to build
  (`t7prep.py`). For contrast, ℚ(ζ₅)'s analogous cores are all empty.
- The T₆ 8-core (2,887v) **can** be 5-coloured avoiding all 1,456 of its
  distance-2 pairs — verdict SAT, so no aggregate forcer at that depth.
- **Two frontier instances ship with this repo, unresolved:**
  `t7pair_dial.cnf` (94,795 vars — *can the designated pair share a colour
  in any proper 5-colouring of the 10-core?* UNSAT here **is** the object B)
  and `t7core_dial.cnf` (*can the 10-core dodge all 37,032 distance-2
  pairs?*). Both survived ~9 unresolved hours of kissat on an M3 before we
  stopped; solver time is not evidence of anything — the instances are
  simply **open**, and they are the most concrete handle on B we can offer.
  Regenerate or scale them with `t7pair_prep.py` / `t7prep.py`.

## Instrument lessons (paid for in full, see FINDINGS.md)

1. Sampling frequency conflates "forbidden" with "unvisited"; only directed
   SAT queries and pattern enumeration are trustworthy at scale.
2. A fast aggregate UNSAT demands a core-reduction check: ours collapsed to
   the unit pentagon (K₅ in two-distance terms) — a triviality, caught
   before it was claimed.
3. Single-proof parallelism (cube-and-conquer, two designs) did not pay on
   these instances; parallelise across independent instances instead.
4. Solver runtime is a steering signal, not a result.

## Licence

MIT. If any of this is useful in an actual attack on the problem, that is
the best possible outcome — take it.
