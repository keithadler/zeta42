# Findings ledger

Measured results, in session order. Everything here was computed on this
machine with exact arithmetic; "measured" means reproducible from the named
script, not inferred.

## Reproductions (instrument validation)

1. **De Grey's assembly counts** — J 31/13, K 61/26, L 121/52 copies of H,
   matching arXiv:1804.02385 exactly. His key lemma (no 4-colouring of L keeps
   every H free of a mono √3-triple) verified UNSAT in 0.09s. `degrey.py`
2. **De Grey's G (1581v, 7877e)**: built from the paper's 7-step recipe, all
   intermediate counts match (39 → 397 → 791 → 1581). Not 4-colourable
   (421s UNSAT), 5-colourable (0.0s) — χ(G) = 5 verified end-to-end. `degrey_g.py`
3. **Haugland's 21-vertex seed** in Q(ζ₄₂): exact placement found by phase
   search; 84 distinct arc directions match his paper. `heptagon.py`, `cyclo.py`
4. **Haugland's G₁ (740v, 3985e)**: T₅ = 1042, T₆ = 12856, 7-core = 740/3985,
   all matching the paper. `hept_g1.py`
5. **Parts' 510 (via Heule's CNP-SAT repo)**: parsed exactly into
   Q(√3,√5,√11); our edge detection reproduces all 2504 published edges with
   zero discrepancies. `parts510.py`

5b. **Haugland's G₁ forced-pair lemma VERIFIED**: no 4-colouring of G₁ gives
   0 and i√3 the same colour — UNSAT in 27,607s (7h40m) with clique symmetry
   breaking (the unbroken encoding was killed unfinished at 10h22m; kissat
   cross-check pending). The heptagon pipeline is validated end-to-end.
   Calibration anchor: level-4 forced-pair proofs cost ~8 CPU-hours at 740v.
   `g1_lemma_race.py`
5c. **Haugland's G₂ (1066v, 6264e)** matches the paper; both rotation centres
   are roots of unity (ζ₃, ζ₆), so G₂ is native to Q(ζ₄₂). G₂'s published
   property (not yet re-proven here): (−1,0) and (1,0) are forced MONO in
   every 4-colouring.
5d. **Haugland's G₃ (2131v, 12530e)** — exact match to the paper, built in
   Q(ζ₄₂, √5), degree 24, in 6s (`cycloext.py` + unchanged `udg.py`).
   Spindle geometry verified exactly: the radius-2 chord rotation
   (cos 7/8, sin √15/8) about (−1,0) sends (1,0) to (3/4, √15/4), unit
   distance from (1,0). Direct non-4-colourability check queued (overnight
   class). The quadratic-over-cyclotomic wrapper is thereby production-tested:
   mixed-field exploration has no remaining entry cost.

## Negative results (each killed a strategy)

6. **Rotation-closures cap at χ = 4.** De Grey's own rotation pair closed to
   1303 vertices is still 4-colourable. Naive "lattice + rotations" families
   never reach χ = 5, let alone 6. `hunt.py`
7. **Local search loses to CDCL here.** A tuned Markov-chain colourer failed
   to 4-colour (in 300k steps) a graph CaDiCaL solves in 11ms. Structured
   instances favour complete solvers; also local search is one-sided (can
   never prove UNSAT). `fast.py`
8. **Uniform level-5 floppiness of every known 5-chromatic graph.** At k=5,
   zero never-mono pairs among 1.65M pairs across 510/517/529/553/G-1581
   (200 diverse samples each). Stiffest pair in the 510: mono in 75/1000
   samples, >100 supporting patterns (cap). No near-forced pairs exist in
   the multiquadratic shelf. `rigidity.py`
   *Positive control PASSED: 150 diverse proper 4-colourings of G₁, the
   proven-forced pair never monochromatic. The instrument is validated and
   these nulls are citable.*
9. **The 510 is fractionally weak.** α(510) ≥ 145 (greedy found 144 in
   seconds), ratio ≥ 0.284 vs the < 0.25 needed for fractional-record
   relevance. Discrete minimality ≠ fractional strength; the n/α dial is
   disqualified as a growth fitness function for this family. `alpha510.py`
10. **Composition does not stiffen (two mechanisms).**
    - Rotation-coupled fan: best pivot couples only 30 vertices/turn;
      2/3/5 copies leave the dial frozen (median 40/200). `compose.py`
    - Translation-welded chain: 146/510 shared vertices per copy (29%);
      dial still frozen across 2/3/5 copies.
    Overlapping loose 5-chromatic cores does not compound constraint at k=5.

11. **The B-hunt: every laptop-reachable route measured dead.** The program
    reduced (via Parts arXiv:2010.12656 + de Grey's G126) to one open object —
    B: a unit-only graph forcing a designated phi-pair non-mono in every
    5-colouring; substitute B into G31's 28 phi-edges and chi >= 6 follows.
    Routes tried and killed:
    - sqrt3-shelf, specific and aggregate forcing: floppy (results 8, and
      `twodist.py`: all shelf graphs dodge every sqrt3-pair instantly).
    - Pentagonal path-template (Haugland G1 lifted): 7/8-cores all EMPTY on
      both the 34-arc height-1 and the full 90-arc D=11 alphabets — the
      pentagon path-lattice cannot sustain template density. `pent_b.py`
    - Minkowski ladder n=7..13 (up to 2380v): min phi-pair frequency falls
      monotonically (77 -> 0 per 1000) BUT all stiff pairs have 50+ supporting
      patterns and every "never-mono" pair is SAT in 0.0s when asked directly.
      The gradient was sampling bias, not forcing.
    Positive discovery en route: Q(zeta_5) has exactly 90 unit arcs over
    denominator 11 (norm-equation census) — denser commensurable fuel than
    the 84-over-7 that powers the level-4 world. Unused by any construction
    so far.
11b. **The 90-arc W-analogue is thin.** Pairwise sums of the 90 D=11 arcs
    give 4141 points but only 8370 unit edges (mean degree ~4); the graph is
    4-colourable and nowhere near aggregate phi-forcing. Commensurability did
    not convert into unit-distance density in Q(zeta_5) the way it does over
    denominator 7 at level 4. The last unexplored lead of the phase.
12. **Instrument lesson: the frequency dial is unreliable at scale.** Sample
    mono-frequency conflates "forbidden" with "unvisited"; as graphs grow the
    sampler's coverage collapses and rare-but-easy events read as stiff.
    Pattern counts and directed queries are the only trustworthy measures.

13. **Single-proof parallelism does not pay on this instance class.** Two
    cube-and-conquer designs benchmarked on the G1 lemma CNF (57k vars;
    sequential yardsticks 7h40m pysat / ~2h15m kissat):
    - march_cu unbounded lookahead stalls (>30 min, no cubes); depth-12 cap
      fixes cubing to 2s / 4096 cubes.
    - Per-cube kissat (6 workers): ~1.2x wall-clock, ~5x WORSE per-core
      (startup + preprocessing per cube). `cnc.py`
    - Incremental iglucose stripes (6 persistent workers): killed at 85 min
      with zero stripes complete — at best parity, likely worse; iglucose's
      2013-era search loses more than incrementality saves. `cnc2.py`
    Conclusion: parallelise at the JOB level (independent instances across
    cores — screens, dials, parameter sweeps — where 6-8x is automatic), and
    run single hard proofs on one fast core under kissat.

13b. **The ninefold world Q(zeta_9) has no Haugland solution.** `nonagon.py`
    - No seed: a unit equilateral triangle with a vertex on each of three
      {9/k} circumcircles does not exist even under free rotation (closest
      miss 0.068; the n=7 control is exact; none for n=11, 13 either).
    - Best exact phasing in Q(zeta_18) (degree 6 -- zeta_6 comes free):
      {9/1},{9/2},{9/4}, 27v/81e, 6-regular, chi 3, and all 18 edge
      directions are roots of unity -- a lattice patch, no new arcs.
    - Arc census (height <= 3, method reproduces Q(zeta_5)'s 90 over 11):
      denominator 3: 0, 7: 36, 13: 36, 19: 468 (+18 roots).  3 is the only
      ramified prime and it ramifies in Q(sqrt-3) too, so it yields nothing;
      Q(zeta_42)'s 7 is ramified in zeta_7 AND split in zeta_3.
    - G1 template on 54-, 90- and 486-arc alphabets, pairs 0--i*sqrt3 and
      0--2, paths <= 6 (<= 5 for 486): every 7-core EMPTY.  Largest
      survivors (5-core 776v and 6-core 621v on 90 arcs; 6-core 369v on 54)
      are 3-chromatic, pair free at k=4 and 5.  The 486-arc system: no T5 vertex has 7 neighbours
      in T4.  Arc count is not density -- a data point for note.md section 4.

13c. **Ramified-and-split test: Q(zeta_20) PASSES on density, not on forcing.**
    `zeta20.py`.  Hypothesis from 13b: density needs a prime ramified in the
    n-part and split in the CM direction (7 in Q(zeta_42)); 5 in Q(zeta_20)
    is one (e=4 in Q(zeta_5), split in Q(i)).  Pass/fail fixed in advance.
    - Arcs over 5: exactly 160 (= 20 roots x a in +-1..+-4), 40 per depth;
      the doubly ramified 2 gives 0.  Counts match the ramification exactly.
    - G1 template, pair 0--phi (phi = |1+zeta_5| native): depth-1 arcs alone
      (60) give T5 = 1558 and a 7-core of 2132v/14986e; all 180 arcs give a
      7-core of 4580v and an 8-core of 2980v.  All cores CHI = 4 -- the first
      4-chromatic template cores outside the heptagonal world (Q(zeta_9),
      Q(zeta_5): empty or 3-chromatic).  Pair 0--2: 7-core 213v, chi 4.
    - No forcing: every pair is free at k = 4 and k = 5 (SAT, instantly),
      where heptagonal G1 is 5-chromatic and forced at k = 4.
    Q(zeta_20) has no unit triangles (no zeta_3); that is the prime suspect
    for the chi = 4 ceiling.  Next test: Q(zeta_60).

## Interpretation (running)

- Level-4 forcing primitives are tight (spindle: 7 vertices force 4 colours).
  The smallest level-5 primitive known is ~510 vertices — constraint density
  two orders thinner, and measured not to percolate across welds. A level-5
  forcing gadget, if it exists, plausibly needs a new *tight* primitive, not
  an assembly of known cores.
- The costs are asymmetric everywhere: finding colourings is ms, proving
  their absence is hours (G's UNSAT 421s at 1581v; G₁'s forced-pair lemma
  >10h across three encodings as of this writing).

14. **G₃ non-4-colourability: verification stopped by choice, not fault.**
    Our own kissat run on the exported instance (8524 vars, 52k clauses) was
    stopped at 7h46m without a verdict — a resource decision; the duration was
    within expected super-linear scaling from the lemma (2.7× the clauses).
    χ(G₃) = 5 stands on Haugland's published verification plus our structural
    match (2131/12530 exact) and 5-colourability (verified, 6s). The chain is
    reproduced end-to-end except this one UNSAT, which any future session can
    re-run from `g3_4col.cnf`.

## Open jobs at time of writing

- G₁ forced-pair lemma: original encoding, symmetry-broken encoding, kissat
  on exported DIMACS (`g1_lemma.cnf`) — all still running.
- Positive control: sampling G₁ at k=4, known forced pair must never go mono.
- α(510) exact value: walking up from 145.
