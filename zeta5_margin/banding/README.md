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

## Still missing

The arithmetic half for banded numerators: Fauzan's prime-by-prime valuation bounds (inner
allocation, outer classes and rank correction) redone with per-band exponents and h > #poles.
That is the research step. `../arithmetic/slack.py` and `gamma_in.py` are the tools to start
from.

## Reproduce

```
python equilibrium_general.py                           # optimum R for each profile
python certify_general.py b443200.json                  # ~10 min
python certify_general.py b54210.json
python search_general.py "bands (4,4,3,2,0,0) x1.20" 128 out.json   # rebuild a measure
```
