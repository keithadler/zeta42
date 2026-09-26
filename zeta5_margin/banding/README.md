# Banded constructions: the real half, certified

Banding comes from branch `claude/amazing-franklin-qz29gp` (`zeta7/zeta5_bands.py`,
`zeta5_deep.py`). It gives each band of small j its own numerator exponent and uses more rows
than poles. It beats Fauzan's construction on the actual values by about 30% (n ≤ 3), but it
has no proof. This directory proves the **real (analytic) half** for the two best banded
profiles.

## Result

For a profile with exponent density ε(u) (u = j/K) and rows h = λK, the field is
`V(t) = 2π√t − ∫₀¹ ε(u) log(t+u²) du`, and

    log Δ_K(ζ5) ≤ κ K² log K + R K² + o(K²),   κ = 2λ² + 2λ∫ε.

R is certified by a nested rational arcsine measure (`certify_general.py`, arb ball arithmetic).

| construction | λ | optimum R | **certified R** | measure |
|---|---|---|---|---|
| Fauzan (5 on j ≤ 3n) | 0.925 | −4.03855 | −4.03797 | `../measure256.json` |
| bands (5,4,2,1,0), t, rows ×1.15 | 1.00625 | −4.50859 | **−4.50609** | `b54210.json` (84,885 cells) |
| bands (4,4,3,2,0,0), rows ×1.20 | 1.02 | −4.71975 | **−4.71745** | `b443200.json` (85,552 cells) |

**Control.** On Fauzan's own profile and measure the general checker returns
−4.020031555 = U − C_*, the paper's certified value.

## What this does and does not show

R on its own does **not** say whether banding is better. Banding also changes κ (the
K² log K term), and κ is balanced by the arithmetic half, so only the total
`log P/K² = R + B` is comparable across constructions. B is the arithmetic constant, with the
log K terms cancelled.

What it does show is the bar the unproven arithmetic half has to clear.

| | actual log P/K² | arithmetic loss allowed to beat today's record (0.0868) | margin if its arithmetic proof is as tight as Fauzan's is now |
|---|---|---|---|
| Fauzan | −0.1253 (n = 3) | (0.0379 used) | 0.0868 |
| bands (4,4,3,2,0,0) ×1.20 | −0.1629 (n = 3) | **0.0738** | **≈ 0.123** |
| bands (5,4,2,1,0) t ×1.15 | −0.1686 (n = 2) | **0.0793** | **≈ 0.128** |

So an arithmetic proof for banding can be about **twice as lossy** as today's proof for Fauzan's
construction and still beat the record. If it is equally tight, the margin rises about 45% over
today's record. The actual values are small-n (K ≤ 120) and drift with K, so the right-hand
columns are estimates. The certified R values are exact bounds.

## Arithmetic half: generalised bounds and the resulting margin (`arith/`)

The proof's three prime ranges generalise to any band profile and any number of rows:

| primes | bound | check against the true valuations |
|---|---|---|
| p > K | scalar only (v_p(Δ) = 0) | v_p(Δ) = 0 at every such prime, K = 40 and 80 |
| K/2 < p ≤ K | `gout_general.py`: per-class outer count, numerator credit (e_a+1)/2, extra rows at weight 0, rank correction `r = 2h + 2 − 2p + deg W − #poles` | reproduces the paper's γ^out exactly; **exact** on banded (14 primes, zero slack) |
| p ≤ K/2 | `gin_general.py`: inner basis, each row weighted by the minimum of its source half-bounds | 0 violations, paper and banded; banded within 0–5 of the truth |

Asymptotic constant (`limits_fast.py`): sample the limiting functions with a large prime and
integrate. **Control:** on Fauzan's profile it reproduces R(x), the outer function, and
A = A_* − 187/9600 to within 0.0085 (the calibration below).

| construction | A (raw … calibrated) | certified U = R + C* | **provable margin** | vs today's record 0.0868 |
|---|---|---|---|---|
| bands (4,4,3,2,0,0) ×1.20 | 1.8105 … 1.8190 | −1.92891 | **0.110 … 0.118** | +27 … +36% |
| bands (5,4,2,1,0) ×1.15 | 1.7740 … 1.7825 | −1.87599 | **0.094 … 0.102** | +8 … +18% |

True arithmetic constant at finite K for (4,4,3,2,0,0): 1.5505 (K=40), 1.6848 (80), 1.6805 (120).

**Status.** The real half is certified. The arithmetic half is the paper's argument generalised
and tested against exact valuations. Its constant is a numerical integral of the limiting
functions, not yet an exact rational, and the bounds are not yet written as proofs. The
irrationality measure μ for banding also needs its own relative-norm constant ϱ (deg Q = #poles,
different numerator), which has not been derived.

Useful structural fact: at large K/p the inner function has mean slope F̄ = λ(4 − e₁), with e₁
the first band's exponent (Fauzan: e₁ = 5, F̄ = −λ). This term is a large part of A.

## Still missing

Written proofs of the three generalised bounds, exact integration of the banded limiting
functions (as done for the split in `../arithmetic/exactA.py`), and ϱ for the banded
construction.

## Reproduce

```
python equilibrium_general.py                           # optimum R for each profile
python certify_general.py b443200.json                  # ~10 min
python certify_general.py b54210.json
python search_general.py "bands (4,4,3,2,0,0) x1.20" 128 out.json   # rebuild a measure
```
