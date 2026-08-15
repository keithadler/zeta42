"""Positive control for the rigidity instrument.

G1 at k=4 has a PROVEN forced pair (two engines: pysat-CaDiCaL 7h40m and
kissat <=2h15m, both UNSAT).  The sampling instrument that produced the
level-5 floppiness nulls must, pointed at G1 at k=4, report that pair as
never monochromatic.  If it ever sees the pair mono, a sample was improper
(bug); if the pair's row is anything but zero, the nulls are not citable.

Every line is flushed; progress survives a crash or reboot.
"""

import time
from fractions import Fraction

from hept_g1 import load_G1, SCALE
from udg import UDG
from rigidity import sample_colourings
from colour import verify_colouring

N_SAMPLES = 150


def main():
    fld, core, B, _ = load_G1(log=lambda s: None)
    tofield = lambda v: tuple(Fraction(x, SCALE) for x in v)
    g = UDG(fld, [tofield(v) for v in sorted(core)])
    A_id = g.index_of(tofield((0,) * 12))
    B_id = g.index_of(tofield(B))
    print(f"control: G1 {g.summary()}, k=4, forced pair ({A_id},{B_id})",
          flush=True)

    t0 = time.time()
    mono_hits = 0
    done = 0
    # sample in small batches so progress is visible and crash-resistant
    for batch in range(N_SAMPLES // 10):
        cols = sample_colourings(g, 4, 10, log=lambda s: None)
        for col in cols:
            ok, bad = verify_colouring(g, col)
            assert ok, f"sampler returned an improper colouring: {bad}"
            if col[A_id] == col[B_id]:
                mono_hits += 1
        done += len(cols)
        print(f"  {done} samples, forced pair mono {mono_hits}x  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if mono_hits:
            print("CONTROL FAILED: instrument saw the proven-forced pair "
                  "monochromatic — every null result is suspect", flush=True)
            return
    print(f"CONTROL PASSED: {done} diverse proper 4-colourings, "
          f"forced pair never mono (as proven)", flush=True)


if __name__ == "__main__":
    main()
