# ζ(5): larger proven margins for Fauzan's proof

Two independent improvements: a better comparison measure (real side, certified), and a
wider inner prime range (arithmetic side, exact integrals, proof modification). Combined:
decay 79.07 → 138.9 (M = 100000) and irrationality measure 260 → 148.

## Real side: a better comparison measure

Fauzan's proof (17 Sep 2026; working edition and Lean project at
[long-mathematics/zeta5-irrationality](https://github.com/long-mathematics/zeta5-irrationality))
ends with

    limsup K⁻² log Q_{K,M}(ζ(5))  ≤  A_M + U,

where `A_M` is the arithmetic (p-adic) constant and `U` the real one. `U` is certified with a
hand-built comparison measure: 16 nested arcsine pieces, `U = −1.366996`. This directory
replaces that measure with a better one. **Nothing else in the proof changes.**

## Result

| | paper | `measure64.json` | `measure256.json` |
|---|---|---|---|
| certified `U` (upper bound) | −1.366996 | **−1.383658** | **−1.384938** |
| cells checked | 684 | 31,813 | 113,890 |
| `Q_{40n,200}(ζ5) < exp(−c n²)` | c = 27.8 | c = 54 | **c = 56.5** |
| `Q_{40n,100000}(ζ5) < exp(−c n²)` | c = 79.07 | c = 105 | **c = 107.7** |
| irrationality measure: `\|ζ(5) − a/b\| > b^−μ` | μ = 260 | μ = 198 | **μ = 191** |

The margin `−(A_M + U)` goes from 0.0174 to 0.0354 per K² at M = 200 (about ×2) and from
0.0494 to 0.0674 at M = 100000 (about +36%). Run `python constants.py` for the exact rational
witnesses. For μ = 191 they are ε = 1077/16000, c = 55/4, ρc − 191 = −137/400 and
ε − λ/c = 7/176000.

**Ceiling.** The best that *any* comparison measure can give is the weighted equilibrium value,
`U* ≈ −1.38544` (`equilibrium.py`, converged over 200 to 1600 cells). The paper leaves 0.0184
per K² on the table. `measure256.json` recovers 0.0179 of it (97%), so the real side is now
essentially exhausted. The next gains have to come from `A_M`.

## What is checked, and how

`certify.py` makes the argument of the working edition's `verify_A` (Appendix A) rigorous for
any nested arcsine measure:

* rational data, exact mass 37/40, nested supports in (0, 2), every length > 1/225;
* `2U^ρ(t) − V(t) < M₀` on [0, 2]. Nesting makes `U^ρ` monotone between consecutive
  endpoints, and `V` has a single minimum at `q ∈ [q₋, q₊]`, so a cell passes if
  `2·max(U(l), U(r)) − inf V < M₀`. Failing cells are bisected;
* the tail condition of the confinement lemma at t = 2;
* the closed-form energy `I(ρ)`, `C_*`, and `U = λM₀ − I(ρ) + C_*`.

The arithmetic is python-flint `arb`: rigorous ball arithmetic, outward-rounded, at 200 bits.
**Control:** on the paper's own data (`paper16.json`, MIT, from long-mathematics) it reproduces
the published certificate exactly: 684 cells, worst bound −6.645002689, and
U = −1.366995564511.

The lemmas that depend on ρ are the comparison-measure proposition, the regularisation lemma
and the confinement lemma. The regularisation lemma's constant `60√ε` uses only
`Σ c_j = λ < 1` and `b_j − a_j > 1/225`, so it holds for any number of components. All three
are covered by the checks above.

## Caveats

* This is a Python certificate, not a Lean proof. To adopt it, the Lean project's interval
  checker would have to replay 31,813 or 113,890 cells instead of 684. That is more work but
  the same kind of work. `measure64.json` captures 90% of the gain at about a quarter of the cells.
* μ = 191 comes from the paper's own Appendix C argument with only ε changed. The relative-norm
  constant ρ = 6933/500 is already optimal for that argument: optimising `a, b` improves κ only
  from 13.8659 to 13.8635.
* The arithmetic constant `A_M` is taken from the working edition unchanged. It was not
  re-derived here.

## Arithmetic side: extending the inner range to K/p ≥ 5/2 (`arithmetic/`)

**Where the slack is.** For K = 40, 80, 120, `arithmetic/slack.py` compares, prime by prime,
the true valuation of the content of `F_K = S_K Δ_K` with what the proof credits. Two ranges
behave differently:

* **Outer primes with p/K ≥ 0.44:** the proof's outer bound (class weights plus the rank
  correction `min(r, z)`) is exact to within 1 at every prime (`gamma_out.py`). No slack.
* **1/3 < p/K < 0.433:** the true valuation beats the credited one by up to about 0.7·K per
  prime. The per-prime values agree at K = 80 and K = 120, so this is asymptotic, not
  finite-size. At p = 41 and 43 (K = 120) the truth even exceeds the outer weights *before*
  the rank loss, so the outer basis itself is the wrong tool just above K/3.

**Change.** Use the paper's *inner* bound, meaning its allocation basis, weights and zero class
(Proposition 4.1 of the working edition), for `2K/5 ≤ p ≤ K/3`, and the outer bound only for
`p > 2K/5`. In the inner proof, `x = K/p ≥ 3` is used only for the positivity estimate
`T − b_a > 2λx − 81/20 ≥ 3/2`. The argument needs only `> 0`, and `2λ·(5/2) − 81/20 = 23/40`.
The degree checks, the zero-source dominance and the extra-row case table do not use x ≥ 3.
In this range `m_N = 0` and `p² > 5K`, and the outer hypotheses still hold for `p > 2K/5`.
The limiting-function lemma and the prime sum go through unchanged with `x ∈ [5/2, M]`.

**Exact gain.** The arithmetic constant changes by

    ΔA = ∫_{5/2}^{3} R(x) x⁻³ dx − ∫_{1/3}^{2/5} (outer integrand) dy
       = 10433/48000 − 1421/6000 = −187/9600.

This is exact sympy integration on pieces where R is checked to be polynomial (`exactA.py`).
**Control:** the same code reproduces the paper's own exact integrals, ∫₃²⁰ R x⁻³ = 3224…/7535…
and I_out = 127751/96000 (`control.py`). The optimal split, 170/69, gives only −3203/163200.

**Empirical check.** `test_in.py` evaluates the exact finite bound γ_p^in, with every
hypothesis of the proposition enforced, against the true valuation at K = 40, 80, 120. There
are **0 violations**, including all 5 primes in the new window. A further run at K = 160
(`rows160.json`) adds p = 59, 61 in the window, again with no violation. The bound is exactly tight at
two of them, and there it beats the outer bound by 17–62 valuation units per prime.

**Combined result** (certified real side plus this change): margin 0.0868 per K² at
M = 100000; `Q < exp(−87.5 n²)` at M = 200 and `exp(−138.9 n²)` at M = 100000; and
`|ζ(5) − a/b| > b^−148` (ε = 1389/16000, c = 533/50).

**Status: a proof modification, not yet a written proof.** The integrals are exact and the
real side is certified. The claim that Proposition 4.1 extends to 5/2 ≤ K/p < 3 is argued above
step by step and checked numerically at small K. Small K cannot reach the admissible regime
K ≥ 200M², so it still needs to be written out in full and checked by the formalisation.

## Banded constructions

`banding/` certifies the real half for the two best banded profiles from branch
`claude/amazing-franklin-qz29gp` (R ≤ −4.71745 and −4.50609). It also shows how lossy their
still-unproven arithmetic half may be and still beat the record. See `banding/README.md`.

## Reproduce

```
pip install python-flint numpy scipy cvxpy
python certify.py paper16.json       # control, <1 s
python certify.py measure64.json     # ~25 s
python certify.py measure256.json    # ~5 min
python constants.py                  # theorem constants, exact
python equilibrium.py                # the ceiling U*
python search.py 256 out.json        # rebuild a measure (~10 min), then certify it
cd arithmetic
python test_in.py                    # extended inner bound vs truth, K <= 120 (~2 min)
python exactA.py                     # dA = -187/9600, exact
python control.py                    # reproduces the paper's exact integrals
```
