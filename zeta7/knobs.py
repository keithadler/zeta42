from profile import score, staircase
K = 40
def show(tag, e, h):
    lD, lc, lP = score(7, e, h)
    print(f"{tag:42s} h={h:2d} | logD/K^2 {lD/K**2:+.4f}  -logc/K^2 {-lc/K**2:+.4f}  logP/K^2 {lP/K**2:+.4f}", flush=True)
print("-- knob 1: fewer rows than poles")
for h in (35, 32, 28, 24, 20, 14):
    show("e=4 on j<=5, poles (5,40]", staircase(K, [(5, 4)], (5, 40)), h)
print("-- knob 2: numerator factors above the poles")
for M, r in [(45, 1), (50, 1), (60, 1), (80, 1), (50, 2), (60, 2), (80, 3)]:
    show(f"e=4 j<=5, poles (5,40], +{r} on (40,{M}]", staircase(K, [(5, 4)], (5, 40), [(M, r)]), 35)
