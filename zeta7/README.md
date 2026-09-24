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

Requires `pip install python-flint mpmath`.
