import sys
from limits_fast import A_const
ex = tuple(int(x) for x in sys.argv[1].split(',')); alpha = float(sys.argv[2]); lam = float(sys.argv[3])
A = A_const(ex, alpha, lam, 2.0)
print(f"{ex} x{alpha}: A_raw = {A:.5f}   calibrated (+0.0085 from the paper control) = {A + 0.0085:.5f}", flush=True)
