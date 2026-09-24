# ζ(7): does the ζ(5) Hankel-determinant method extend?

**Short answer: not directly. In every configuration tried the numbers point the wrong way, by a wide margin.**

Context: A. Fauzan's preprint *ζ(5) is irrational* (17 Sep 2026), formalised in Lean at
[mo271/zeta5](https://github.com/mo271/zeta5) (`#print axioms` → `propext, Classical.choice, Quot.sound`).

## The method in one paragraph

Take the positive weight `w(y) = y^s F^{(s-1)}(y) / ((s-1)!/2)` with `F = 1/(e^{2πy}-1)` and the
functional `μ_X` on rational functions of `t = y²`:

* `μ(t^k) = (-1)^k B_{2k+2} (2k+s)! / ((2k+2)! (s-1)!)`
* `μ(1/(t+j²)) = j^{s-1}(X - H_j^{(s)}) - 1/(s-1) + 1/(2j)`, where `X` stands for ζ(s)

The pole formula (Hermite's formula) was derived numerically here for s = 7 to 40 digits
(`hermite7.py`: `r_s(j) = -1/(s-1) + 1/(2j)`, which reproduces the paper's `-1/4 + 1/(2j)` at
s = 5). The polynomial moments were checked by quadrature. `Δ_K(X)` is the Hankel determinant of
`μ_X(W(t) t^{a+b} / D_K(t))` with `a, b < h = K - N`, where `D_m = ∏_{j≤m}(t+j²)` and the paper takes
`W = D_N^6`, `K = 40n`, `N = 3n`. `Δ_K(ζ(s)) > 0` is a Gram determinant, and `deg Δ_K = h` is linear in
`n`. So if the primitive integer polynomial `P_K` satisfies `log P_K(ζ(s)) ≤ -c n²`, then ζ(s) is
irrational. The height of `P_K` does not enter.

## Results (`logP/K²`: negative is good, and irrationality needs it to stay negative as K → ∞)

`hankel_s.py 5 40 3 6 1` reproduces the paper's ζ(5) value `log P = -265.13` exactly.

| s | numerator W | K | log Δ(ζ)/K² | −log content/K² | **log P(ζ)/K²** |
|---|---|---|---|---|---|
| 5 | 1 (N=0) | 40 | −2.842 | +2.832 | **−0.010** |
| 5 | D_3^6 (paper) | 40 | −1.148 | +0.982 | **−0.166** |
| 7 | 1 (N=0) | 40 | −2.741 | +3.711 | **+0.970** |
| 7 | D_3^6 (paper's shape) | 40 | | | **+0.767** (+0.846 at K=80) |
| 7 | best single level, D_5^5 | 40 | −0.411 | +1.165 | **+0.754** |
| 7 | best single level, D_10^5 | 80 | +0.021 | +0.805 | **+0.826** |
| 7 | best two-level, D_3^4·D_8^2 | 40 | −0.359 | +1.095 | **+0.737** |

Full grids: `python3 sweep.py 7 40 4,5,6,7 2,3,4,5,6`, `python3 twolevel.py`.

## Why it fails

The real side is almost the same for s = 5 and s = 7 (−2.84 vs −2.74 per K² with no numerator). The
arithmetic side costs about 0.88 K² more for s = 7, because the denominators of `H_j^{(7)}` carry
`d_j^7` instead of `d_j^5`. The numerator `W` only trades one side against the other: raising `N·r`
improves the content about as fast as it spoils `Δ(ζ)`. So the best case stays about 0.74 K² in the
red, and it gets worse with K. ζ(5) itself has only a thin margin (−0.01 K² with no numerator), so
ζ(7) needs a new arithmetic input rather than a parameter tweak. Candidates are extra p-adic savings
in the style of Zudilin or Rivoal, well-poised symmetry, or a different weight.

These are small-K numerics (K ≤ 80), not asymptotic proofs in either direction.

## Round 2: searching the whole family

`profile.py` generalises the construction to any exponent profile: the measure
`w(y) · ∏_j (t+j²)^{e_j}` with `e_j ∈ {-1, 0, 1, …}`, where `e_j = -1` marks a pole. Numerator zeros
can sit below or above the poles or be interleaved with them, and `h` can be any number of rows up
to the number of poles. A change of basis only rescales `Δ`, so the primitive `P` is unchanged, and
this is the complete knob set for this weight. The code reproduces the paper's −265.13 again (the
paper's `D_3^6/D_40` is net exponent 5 on `j ≤ 3`). Scores below are `log P(ζ7)/h²`, which does not
depend on the size of the construction.

| Try | Best `log P/h²` |
|---|---|
| Hill-climb over 8-block profiles, exponents −1..6, varying h (`search.py`, 807 shapes, 9 restarts, `search_s7_log.txt`) | **+0.69** at 15 poles (a finite-size effect, see below). Every local minimum is positive. |
| Fewer rows than poles (`knobs.py`) | +0.73·K² at best (vs +0.75) |
| Numerator zeros above the poles, as in Ball–Rivoal (`knobs.py`) | worse: +1.09 to +8.4 |
| Numerator zeros at half-integers, `(4t+i²)` with `i` odd (`halfint.py`) | worse: +1.18 to +2.01 |
| No polynomial part (to avoid the Bernoulli cost) | impossible: the X coefficient `∑ res_j j⁶` vanishes when `R = O(t⁻⁴)`, so B is singular |

**Finite size.** For the plain construction (poles 1..K, h = K) the score rises steadily with K:

| K | 10 | 15 | 20 | 30 | 40 | 60 | 80 |
|---|---|---|---|---|---|---|---|
| ζ(5) | −0.34 | −0.14 | −0.08 | −0.06 | −0.01 | +0.03 | +0.07 |
| ζ(7) | +0.52 | +0.69 | +0.82 | +0.89 | +0.97 | +1.02 | +1.07 |

So small constructions always look better, and the search's best shapes are just small. Plain ζ(5)
also turns positive, which is why the paper needs `D_N^6`.

**Scaling in s** (plain construction, K = 40): ζ(3) **−0.95**, ζ(5) **−0.01**, ζ(7) **+0.97**,
ζ(9) **+1.86**. Each step of 2 in s costs about 0.93·K². The real side barely moves (−2.96, −2.84,
−2.74, −2.65); the arithmetic side grows. `diag.py` (a diagnostic, not a valid construction) puts
the extra cost in the denominators of `H_j^{(s)}`: removing them saves 1.21·K² for s = 5 and 2.03·K²
for s = 7. The 0.82 difference is essentially the whole ζ(7) gap.

**Bottom line.** Within this family the method proves ζ(3) easily and ζ(5) with a thin margin (the
paper's result), and it misses ζ(7) by roughly one full step. Closing the gap needs a new source of
p-adic savings on `∑ c_j j⁶ H_j^{(7)}` of about 0.8·h², which is more than any numerator profile tested
here provides. It also cannot come from dropping the polynomial part.

## Round 3: leads outside the original family

**Lead 1: more rows than poles** (`manyrows.py`). Gram positivity holds for any number of rows h,
while `deg Δ ≤ m = #poles`. So h > m is a legal, untested regime, and for fixed m the criterion only
asks for `P(ζ) → 0`. `Δ` is computed by exact interpolation of `det(A + xB)`.
* One pole (an Apéry-type linear form with positivity for free): `log P` grows with h even for ζ(3).
  Dead.
* h slightly above m **helps ζ(5)**. Plain ζ(5) (poles 1..m, `log P/m²`): m = 20 goes from −0.083
  (h = 20) to −0.211 (h = 24); m = 40 goes from −0.010 (h = 40) to −0.097 (h = 48). That gain is about
  the same size as the paper's `D_N^6` and could combine with it. It is a side lead for a larger ζ(5)
  margin.
* It does **not** help ζ(7): plain goes from +0.820 to +0.826 at best, and with `D_5^r` it gets worse.

**Lead 2: |Vandermonde|⁴ instead of |Vandermonde|²** (`beta4.py`). By de Bruijn's formula the h-fold
integral with `∏(t_i − t_k)^4` is `h!·Pf[(b−a) μ_X(t^{a+b−1} W/D)]`. It is still a positive polynomial
in X, and `log P_Pf = ½ log P_det` by Gauss's lemma. It is much worse: at degree 14, ζ(5) scores +1.15
and ζ(7) +2.80. The higher-index Bernoulli moments cost more than the stronger repulsion gains. Dead.

**Open lead: a tail-sum reformulation.** With Hermite's formula,
`μ(R) = X·∑_j res_j j⁶ − ∑_{k} S(k)/k⁷ + (rational)`, where `S(k) = ∑_{j≥k} res_j j⁶`. The ζ(7) gap
is the cost of the `1/k⁷`. Closing it would take rational functions whose tail sums `S(k)` carry about
0.8·h² of extra p-adic divisibility, the analogue of the well-poised symmetry in Ball–Rivoal and
Zudilin transported to this Hankel setting. That is a derivation to attempt, not a parameter sweep.

## Round 4: the full positive family for this weight

**The ζ(7) part of an entry is a Beukers-type integral.** By Hermite's formula,
`∑_j ρ_j j⁶ (ζ(7) − H_j^{(7)}) = (1/720) ∫₀¹ (−log z)⁶ P(z)/(1−z) dz` with `P(z) = ∑_j ρ_j j⁶ z^j`.
So the arithmetic depends only on the residues of `R`.

**Every admissible measure is known.** `R = N/D_K` is a positive measure exactly when `N ≥ 0` on
`[0,∞)`, which happens exactly when `N(y²) = |Π(iy)|²` for a real polynomial `Π`
(Fejér–Riesz / Markov–Lukács). The residues are then `c_j ∝ Π(j)Π(−j)/D'(−j²)`.
* `Π` with integer roots gives the product family `∏(t+i²)^{e_i}` (rounds 1–3).
* Half-integer roots give `(4t+i²)` (round 3, worse).
* Non-product `Π` (`fejer.py`: Legendre-, Apéry- and central-binomial-type sums, degree 3–7, K = 20)
  is **always much worse** than the integer-root product of the same degree. At degree 3 the scores
  are +0.97 to +1.42 against the control's +0.79 for ζ(7), and +0.03 to +0.50 against −0.15 for ζ(5).
* Residue shapes chosen directly (`designs.py`: squared binomials, Apéry-like, reflection-symmetric)
  mostly give `R` that changes sign, so they are not positive measures. That is checked exactly by
  square-free factorisation plus root counting in `residues.py`. The two that are positive are worse.

**Conclusion for this weight.** The integer-root product numerators searched in rounds 1–3 are the
arithmetic optimum of the whole positive family, and they miss ζ(7) by about 0.8·h². To go further
you need a **different weight or a different positive structure**: for example a two-dimensional
Gram/Andréief integral whose kernel produces ζ(7) with different arithmetic, or a weight whose pole
moments bring in `H_j^{(7)}` with built-in cancellation. The `w = y⁷F⁽⁶⁾` family itself is exhausted.

## Round 5: the remaining gaps in the family

* **Root 0 of Π**: a numerator factor `t^m`, i.e. the weight `y^{7+2m}F^{(6)}` (K = 20). For ζ(7),
  m = 1 changes nothing (+0.699 vs +0.698 with `D_2^4`), and larger m is steadily worse
  (m = 10: +1.86). For ζ(5), m = 1 is marginally better (−0.276 vs −0.268).
* **Other structures reduce to the measure.** A Gram matrix of rational functions `p_a/Q` against μ
  is the polynomial Gram against `μ/Q²`. Adding positive point masses only makes the Gram determinant
  larger. A positive integer combination of constructions is at least as large as each term at ζ(7),
  and a signed combination loses positivity. So for one-dimensional Gram structures, the measure is
  the only thing that matters, and the measure family is exhausted.
* **What remains open** is a genuinely two-dimensional positive measure (a non-product density on
  `(0,∞)²`) whose moments are linear in ζ(7) alone. The natural kernels give multiple zeta values
  instead, and no clean candidate is known to us. That would be new mathematics, not a search.

## Round 6: a different weight (poles at half-integers)

Hermite's formula holds for real `a > 0`. At `a = j/2` with `j` odd it gives
`j⁶((127/64)ζ(7) − 2 O_j) − 1/6 + 1/j`, where `O_j` is the sum of `1/m⁷` over odd `m ≤ j`. This was
checked to 40 digits, and it is still linear in ζ(7) alone. Only the denominator 2 works: with 3 or
more, other L-values appear. After rescaling, this is the weight `y⁷F⁽⁶⁾(y/2)`, which decays at half
the rate, so it changes the analytic side for the first time in this search (`halfpoles.py`).

| m poles, h = m | integer poles (control) | half-integer spacing | odd half-integers only |
|---|---|---|---|
| 20 | +0.820 | +1.037 | +3.488 |
| 40 | +0.970 | +1.292 | +4.039 |

A slower-decaying weight is worse. A faster one (only even poles) raises the denominators to
`d_{2K}⁷` for K poles, so it is worse by construction. The paper's integer spacing is the optimum.

## Side result: ζ(5) with more rows (`zeta5_margin.py`)

The paper's construction (K = 40n, N = 3n, `D_N⁶`) with h = α(K − N), at n = 1:
α = 1 gives −265.13 (the paper's value, reproduced), α = 1.1 gives **−286.77**, α = 1.2 gives −266.76,
α = 1.3 gives −217.64. An extra factor `t` adds a little at α = 1 to 1.1 (−280.20 and −291.98). This
is one size only; it still needs checking at n = 2 and 3.

Requires `pip install python-flint mpmath sympy`.
